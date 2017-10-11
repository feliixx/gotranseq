package main

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"os"
	"sync"

	"github.com/jessevdk/go-flags"

	nc "github.com/feliixx/gotranseq/NCBICode"
)

const (
	fastaID = '>'
	endLine = '\n'
	unknown = 'X'
	space   = ' '
	// 10MB ?
	bufferSize = 1000 * 1000 * 10
	nCode      = uint8(0)
	aCode      = uint8(1)
	cCode      = uint8(2)
	tCode      = uint8(3)
	gCode      = uint8(4)

	version  = "0.1"
	toolName = "gotranseq"
)

var (
	suffix = map[int][]byte{
		1: {'_', '1'},
		2: {'_', '2'},
		3: {'_', '3'},
		4: {'_', '4'},
		5: {'_', '5'},
		6: {'_', '6'},
	}
	letterCode = map[byte]uint8{
		'A': aCode,
		'C': cCode,
		'T': tCode,
		'G': gCode,
		'N': nCode,
	}
	spaceDelim = []byte{space}
)

type FastaSequence struct {
	ID       []byte
	Comment  []byte
	Sequence []uint8
}

func printErrorAndExit(err error) {
	fmt.Printf("error: %v\n", err)
	os.Exit(1)
}

func createMapCode(code int) (map[uint32]byte, error) {

	resultMap := map[uint32]byte{}
	twoLetterMap := map[string][]byte{}
	codonCode := [4]uint8{uint8(0), uint8(0), uint8(0), uint8(0)}
	_, ok := nc.Table[code]
	if !ok {
		return nil, fmt.Errorf("no table for code %v", code)
	}

	for k, v := range nc.Table[code] {
		// generate 3 letter code
		for i := 0; i < 3; i++ {
			codonCode[i] = letterCode[k[i]]
		}
		uint32Code := uint32(codonCode[0]) | uint32(codonCode[1])<<8 | uint32(codonCode[2])<<16
		resultMap[uint32Code] = v
		// generate 2 letter code
		m, ok := twoLetterMap[k[0:2]]
		// two letter codon is not present
		if !ok {
			twoLetterMap[k[0:2]] = []byte{v}
		} else {
			m = append(m, v)
			twoLetterMap[k[0:2]] = m
		}
	}
	for l, m := range twoLetterMap {
		uniqueAA := true
		for i := 0; i < len(m); i++ {
			if m[i] != m[0] {
				uniqueAA = false
			}
		}
		if uniqueAA {
			first := letterCode[l[0]]
			second := letterCode[l[1]]

			uint32Code := uint32(first) | uint32(second)<<8
			resultMap[uint32Code] = m[0]
		}
	}
	return resultMap, nil
}

type Options struct {
	Required `group:"required"`
	Optional `group:"optional"`
	General  `group:"general"`
}

type Required struct {
	Sequence string `short:"s" long:"sequence" value-name:"<filename>" description:"Nucleotide sequence(s) filename"`
	Outseq   string `short:"o" long:"outseq" value-name:"<filename>" description:"Protein sequence fileName"`
}

type Optional struct {
	Frame string `short:"f" long:"frame" value-name:"<code>" description:"frame"`
	Table int    `short:"t" long:"table" value-name:"<code>" description:"ncbi code to use" default:"0"`
}
type General struct {
	Help    bool `short:"h" long:"help" description:"show this help message"`
	Version bool `short:"v" long:"version" description:"print the tool version and exit"`
}

func main() {

	var options Options
	p := flags.NewParser(&options, flags.Default&^flags.HelpFlag)
	_, err := p.Parse()
	if err != nil {
		fmt.Printf("wrong args: %v, try %s --help for more informations\n", err, toolName)
		os.Exit(1)
	}
	if options.Help {
		fmt.Printf("%s version %s\n\n", toolName, version)
		p.WriteHelp(os.Stdout)
		os.Exit(0)
	}
	if options.Version {
		fmt.Printf("%s version version %s\n", toolName, version)
		os.Exit(0)
	}
	mapCode, err := createMapCode(0)
	if err != nil {
		printErrorAndExit(err)
	}

	if options.Sequence == "" {
		printErrorAndExit(fmt.Errorf("missing required parameter -s | -sequence. try %s --help for details", toolName))
	}
	f, err := os.Open(options.Sequence)
	if err != nil {
		printErrorAndExit(err)
	}
	defer f.Close()
	if options.Outseq == "" {
		printErrorAndExit(fmt.Errorf("missing required parameter -o | -outseq. try %s --help for details", toolName))
	}
	out, err := os.Create(options.Outseq)
	if err != nil {
		printErrorAndExit(err)
	}
	defer out.Close()

	fnaSequences := make(chan FastaSequence, 10)
	errs := make(chan error, 1)

	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(2)
	for nWorker := 0; nWorker < 2; nWorker++ {
		go func() {

			defer wg.Done()
			var translated bytes.Buffer
			var size int
			var codonCode uint32
			idx := make([]int, 3)

			for sequence := range fnaSequences {
				select {
				case <-ctx.Done():
					return
				default:
				}
				// forward mode
				for frame := 0; frame < 3; frame++ {
					translated.Write(sequence.ID)
					translated.Write(suffix[frame+1])

					if sequence.Comment != nil {
						translated.WriteByte(space)
						translated.Write(sequence.Comment)
					}
					translated.WriteByte(endLine)
					size = len(sequence.Sequence)
					for i := frame + 2; i < size; i += 3 {
						// all 3 letter code exists, so no need to check presence in map
						codonCode = uint32(sequence.Sequence[i-2]) | uint32(sequence.Sequence[i-1])<<8 | uint32(sequence.Sequence[i])<<16
						b, ok := mapCode[codonCode]
						if !ok {
							translated.WriteByte(unknown)
						} else {
							translated.WriteByte(b)
						}
					}
					if (size-frame)%3 == 2 {
						codonCode = uint32(sequence.Sequence[size-2]) | uint32(sequence.Sequence[size-1])<<8
						b, ok := mapCode[codonCode]
						if !ok {
							translated.WriteByte(unknown)
						} else {
							translated.WriteByte(b)
						}
					} else if (size-frame)%3 == 1 {
						translated.WriteByte(unknown)
					}
					translated.WriteByte(endLine)
				}
				// get the complementary sequence
				for i, n := range sequence.Sequence {
					switch n {
					case aCode:
						sequence.Sequence[i] = tCode
					case tCode:
						sequence.Sequence[i] = aCode
					case cCode:
						sequence.Sequence[i] = gCode
					case gCode:
						sequence.Sequence[i] = cCode
					default:
						//case N -> leave it
					}
				}
				// reverse it
				for i, j := 0, len(sequence.Sequence)-1; i < j; i, j = i+1, j-1 {
					sequence.Sequence[i], sequence.Sequence[j] = sequence.Sequence[j], sequence.Sequence[i]
				}

				switch len(sequence.Sequence) % 3 {
				case 0:
					idx[0] = 0
					idx[1] = 2
					idx[2] = 1
				case 1:
					idx[0] = 1
					idx[1] = 0
					idx[2] = 2
				case 2:
					idx[0] = 2
					idx[1] = 1
					idx[2] = 0
				}

				// complement mode
				for j, frame := range idx {
					translated.Write(sequence.ID)
					translated.Write(suffix[j+4])

					if sequence.Comment != nil {
						translated.WriteByte(space)
						translated.Write(sequence.Comment)
					}
					translated.WriteByte(endLine)

					size = len(sequence.Sequence)
					for i := frame + 2; i < size; i += 3 {
						// all 3 letter code exists, so no need to check presence in map
						codonCode = uint32(sequence.Sequence[i-2]) | uint32(sequence.Sequence[i-1])<<8 | uint32(sequence.Sequence[i])<<16
						b, ok := mapCode[codonCode]
						if !ok {
							translated.WriteByte(unknown)
						} else {
							translated.WriteByte(b)
						}
					}
					if (size-frame)%3 == 2 {
						codonCode = uint32(sequence.Sequence[size-2]) | uint32(sequence.Sequence[size-1])<<8
						b, ok := mapCode[codonCode]
						if !ok {
							translated.WriteByte(unknown)
						} else {
							translated.WriteByte(b)
						}
					} else if (size-frame)%3 == 1 {
						translated.WriteByte(unknown)
					}
					translated.WriteByte(endLine)
				}
				if translated.Len() > bufferSize {
					_, err = out.Write(translated.Bytes())
					if err != nil {
						select {
						case errs <- fmt.Errorf("fail to write to output file: %v", err):
						default:
						}
						cancel()
						return
					}
					translated.Reset()
				}
			}
			// some sequence left in the buffer
			if translated.Len() > 0 {
				_, err = out.Write(translated.Bytes())
				if err != nil {
					select {
					case errs <- fmt.Errorf("fail to write to output file: %v", err):
					default:
					}
					cancel()
					return
				}
			}
		}()
	}

	scanner := bufio.NewScanner(f)
	fcr := &FastaSequence{
		ID:       make([]byte, 0),
		Sequence: make([]uint8, 0),
	}

	var bufferedSequence bytes.Buffer
	var readError error
Loop:
	for scanner.Scan() {
		// read the line
		line := scanner.Bytes()
		// skip blank lines
		if len(line) == 0 {
			continue
		}
		// if the line starts with '>'; it's the ID of the sequence
		if line[0] == fastaID {
			if len(fcr.ID) > 0 {
				select {
				case <-ctx.Done(): // if an error occurred in one of the 'inserting' goroutines, close the channel
					break Loop
				default:
				}
				fastaSequence := FastaSequence{
					ID:       fcr.ID,
					Comment:  fcr.Comment,
					Sequence: make([]uint8, bufferedSequence.Len()),
				}
				for i, b := range bufferedSequence.Bytes() {
					uintCode, ok := letterCode[b]
					if !ok {
						readError = fmt.Errorf("invalid char in sequence %v: %v", string(fcr.ID), string(b))
						break Loop
					}
					fastaSequence.Sequence[i] = uintCode
				}
				fnaSequences <- fastaSequence
				bufferedSequence.Reset()
			}

			seqID := bytes.SplitN(line, spaceDelim, 2)
			fcr.ID = seqID[0]
			if len(seqID) > 1 {
				fcr.Comment = seqID[1]
			} else {
				fcr.Comment = nil
			}
		} else {
			bufferedSequence.Write(line)
		}
	}

	if readError != nil {
		printErrorAndExit(readError)
	}
	// don't forget tu push last sequence
	fastaSequence := FastaSequence{
		ID:       fcr.ID,
		Comment:  fcr.Comment,
		Sequence: make([]uint8, bufferedSequence.Len()),
	}
	for i, b := range bufferedSequence.Bytes() {
		uintCode, ok := letterCode[b]
		if !ok {
			readError = fmt.Errorf("invalid char in sequence %v: %v", string(fcr.ID), string(b))
			break
		}
		fastaSequence.Sequence[i] = uintCode
	}
	if readError != nil {
		printErrorAndExit(readError)
	}
	fnaSequences <- fastaSequence

	// close fasta sequence channel
	close(fnaSequences)
	//wait for goroutines to finish
	wg.Wait()
	if ctx.Err() != nil {
		printErrorAndExit(<-errs)
	} else {
		fmt.Println("Done!")
	}

}
