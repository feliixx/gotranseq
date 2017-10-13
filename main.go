package main

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"io"
	"os"
	"runtime"
	"sync"

	"github.com/jessevdk/go-flags"

	nc "github.com/feliixx/gotranseq/NCBICode"
)

const (
	// some constant for parsing fasta
	fastaID = '>'
	endLine = '\n'
	unknown = 'X'
	space   = ' '
	// some AA constant
	stopByte      = '*'
	cleanStopByte = 'X'
	// size of the buffer for writing to file
	maxBufferSize = 1000 * 1000 * 10
	// max line size for sequence
	maxLineSize = 60
	// uint8 code for supported nucleotides
	nCode = uint8(0)
	aCode = uint8(1)
	cCode = uint8(2)
	tCode = uint8(3)
	gCode = uint8(4)
	// general constant
	version  = "0.1"
	toolName = "gotranseq"
)

var (
	// suffix to append to sequenceID to keep track of the frame in
	// the output file
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

// FastaSequence stores a nucleic sequence and its meta-info
//
// fasta format is:
//
// >sequenceID some comments on sequence
// ACAGGCAGAGACACGACAGACGACGACACAGGAGCAGACAGCAGCAGACGACCACATATT
// TTTGCGGTCACATGACGACTTCGGCAGCGA
//
// see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
// section 1 for details
type FastaSequence struct {
	ID       []byte
	Comment  []byte
	Sequence []uint8
}

type fastaChannelFeeder struct {
	IDBuffer       *bytes.Buffer
	CommentBuffer  *bytes.Buffer
	SequenceBuffer *bytes.Buffer
	FastaChan      chan FastaSequence
}

type ioHandler struct {
	in  io.Reader
	out io.Writer
}

func (f *fastaChannelFeeder) sendFasta() error {
	// create a fastaSequence, comments is not required
	fastaSequence := FastaSequence{
		ID:       make([]byte, f.IDBuffer.Len()),
		Sequence: make([]uint8, f.SequenceBuffer.Len()),
	}
	// copy content of buffers to the new object
	copy(fastaSequence.ID, f.IDBuffer.Bytes())
	f.IDBuffer.Reset()

	if f.CommentBuffer.Len() != 0 {
		fastaSequence.Comment = make([]byte, f.CommentBuffer.Len())
		copy(fastaSequence.Comment, f.CommentBuffer.Bytes())
		f.CommentBuffer.Reset()
	}
	// convert the sequence of bytes to an array of uint8 codes,
	// so a codon (3 nucleotides | 3 bytes ) can be represented
	// as an uint32
	for i, b := range f.SequenceBuffer.Bytes() {
		switch b {
		case 'A':
			fastaSequence.Sequence[i] = aCode
		case 'C':
			fastaSequence.Sequence[i] = cCode
		case 'G':
			fastaSequence.Sequence[i] = gCode
		case 'T':
			fastaSequence.Sequence[i] = tCode
		case 'N':
			fastaSequence.Sequence[i] = nCode
		default:
			return fmt.Errorf("invalid char in sequence %v: %v", string(fastaSequence.ID), string(b))
		}
	}
	// push the sequence to a buffered channel
	f.FastaChan <- fastaSequence

	f.SequenceBuffer.Reset()
	return nil
}

func printErrorAndExit(err error) {
	fmt.Printf("error: %v\n", err)
	os.Exit(1)
}

// create the code map according to the selected table code
func createMapCode(code int, clean bool) (map[uint32]byte, error) {

	resultMap := map[uint32]byte{}
	twoLetterMap := map[string][]byte{}
	codonCode := [4]uint8{uint8(0), uint8(0), uint8(0), uint8(0)}

	// load the standard code
	m := nc.Standard
	// if we use a different code, load the difference map
	// and update the values
	if code != 0 {
		for k, v := range nc.TableDiff[code] {
			m[k] = v
		}
	}

	for k, v := range m {
		// generate 3 letter code
		for i := 0; i < 3; i++ {
			codonCode[i] = letterCode[k[i]]
		}
		// each codon is represented by an uint32:
		// each possible nucleotide is represented by an uint8 (255 possibility)
		// the three first bytes are the the code for each nucleotide
		// last byte is uint8(0)
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
	// if clean is specified, we want to replace all '*' by 'X' in the output
	// sequence, so replace all occurrences of '*' directly in the ref map
	if clean {
		for k, v := range resultMap {
			if v == stopByte {
				resultMap[k] = cleanStopByte
			}
		}
	}
	return resultMap, nil
}

// read a fata file, translate each sequence to the corresponding prot sequence in the specified frame
func (i *ioHandler) readSequenceAndTranslate(options Options) error {
	// get the codemap for codon <-> AA translation
	mapCode, err := createMapCode(options.Table, options.Clean)
	if err != nil {
		return err
	}
	// mask for required frames
	framesToGenerate := make([]int, 6)
	reverse := false

	switch options.Frame {
	case "1":
		framesToGenerate[0] = 1
	case "2":
		framesToGenerate[1] = 1
	case "3":
		framesToGenerate[2] = 1
	case "F":
		for i := 0; i < 3; i++ {
			framesToGenerate[i] = 1
		}
	case "-1":
		framesToGenerate[3] = 1
		reverse = true
	case "-2":
		framesToGenerate[4] = 1
		reverse = true
	case "-3":
		framesToGenerate[5] = 1
		reverse = true
	case "R":
		for i := 3; i < 6; i++ {
			framesToGenerate[i] = 1
		}
		reverse = true
	case "6":
		for i := range framesToGenerate {
			framesToGenerate[i] = 1
		}
		reverse = true
	default:
		return fmt.Errorf("wrong value for -f | --frame parameter: %s", options.Frame)
	}
	// a channel of fasta sequences that can be used from
	// multiple goroutines to parrallize the job
	fnaSequences := make(chan FastaSequence, 10)
	// a channel of error to get error from goroutine
	errs := make(chan error, 1)
	// use context to smoothly close all goroutines if
	// an error occurs
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(options.NumWorker)

	for nWorker := 0; nWorker < options.NumWorker; nWorker++ {

		go func() {

			defer wg.Done()
			// buffer to reduce calls to out.Write()
			var translated bytes.Buffer
			// length of the nucleic sequence
			var size int
			// nb of bytes since last '\n' char
			var currentLength int
			// code of the current codon
			var codonCode uint32
			// frame matrix in reverse mode because it depends on sequence
			// length, cf convention
			idx := make([]int, 3)
			// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
			// from right end of the sequence)
			var bytesToTrim int
			// how many times we run Translate loop on the sequence
			var nbRunLoop int

			for sequence := range fnaSequences {
				// if an error occured somewhere, return so
				// wg.Done() is called
				select {
				case <-ctx.Done():
					return
				default:
				}

				size = len(sequence.Sequence)
				nbRunLoop = 0

				idx[0] = 0
				idx[1] = 1
				idx[2] = 2

			Translate:
				// forward mode
				for i, frame := range idx {
					// only generate requested frames
					if framesToGenerate[frame] == 0 {
						continue
					}
					// sequence id should look like
					// >sequenceID_<frame> comments
					translated.Write(sequence.ID)

					if nbRunLoop >= 1 {
						// Translate has been run one already, so we are
						// in reverse mode
						translated.Write(suffix[i+4])
					} else {
						translated.Write(suffix[frame+1])
					}

					if sequence.Comment != nil {
						translated.WriteByte(space)
						translated.Write(sequence.Comment)
					}
					translated.WriteByte(endLine)

					currentLength = 0

					for i := frame + 2; i < size; i += 3 {
						// format sequence: should be 60 char long max
						if currentLength == maxLineSize {
							currentLength = 0
							translated.WriteByte(endLine)
							bytesToTrim++
						}
						// create an uint32 from the codon, to retrieve the corresponding
						// AA from the map
						codonCode = uint32(sequence.Sequence[i-2]) | uint32(sequence.Sequence[i-1])<<8 | uint32(sequence.Sequence[i])<<16
						b, ok := mapCode[codonCode]
						if !ok {
							// this may occur if the codon contains one or more
							// unknown nucleotide ('N')
							translated.WriteByte(unknown)
							bytesToTrim++
						} else {
							translated.WriteByte(b)
							if b == stopByte {
								bytesToTrim++
							} else {
								bytesToTrim = 0
							}
						}
						currentLength++
					}
					// the last codon is only 2 nucleotid long, try to guess
					// the corresponding AA
					if (size-frame)%3 == 2 {
						if currentLength == maxLineSize {
							translated.WriteByte(endLine)
							bytesToTrim++
						}
						codonCode = uint32(sequence.Sequence[size-2]) | uint32(sequence.Sequence[size-1])<<8
						b, ok := mapCode[codonCode]
						if !ok {
							translated.WriteByte(unknown)
							bytesToTrim++
						} else {
							translated.WriteByte(b)
							if b == stopByte {
								bytesToTrim++
							} else {
								bytesToTrim = 0
							}
						}
						// the last codon is only 1 nucleotid long, no way to guess
						// the corresponding AA
					} else if (size-frame)%3 == 1 {
						if currentLength == maxLineSize {
							translated.WriteByte(endLine)
							bytesToTrim++
						}
						translated.WriteByte(unknown)
						bytesToTrim++
					}
					if options.Trim && bytesToTrim > 0 {
						// remove the last bytesToTrim bytes of the buffer
						// as they are 'X', '*' or '\n'
						translated.Truncate(translated.Len() - bytesToTrim)
					}
					bytesToTrim = 0
					translated.WriteByte(endLine)
				}
				nbRunLoop++

				// if in reverse mode, reverse-complement the sequence
				if reverse && nbRunLoop < 2 {
					// get the complementary sequence.
					// Basically, switch
					//   A <-> T
					//   C <-> G
					// N is not modified
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
					// reverse the sequence
					for i, j := 0, len(sequence.Sequence)-1; i < j; i, j = i+1, j-1 {
						sequence.Sequence[i], sequence.Sequence[j] = sequence.Sequence[j], sequence.Sequence[i]
					}

					if !options.Alternative {
						// Staden convention: Frame -1 is the reverse-complement of the sequence
						// having the same codon phase as frame 1. Frame -2 is the same phase as
						// frame 2. Frame -3 is the same phase as frame 3
						//
						// use the matrix to keep track of the forward frame as it depends on the
						// length of the sequence
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
					}
					// run the same loop, but with the reverse-complemented sequence
					goto Translate
				}
				// if the buffer holds more than 10MB of data,
				// write it to output file and reset the buffer
				if translated.Len() > maxBufferSize {
					_, err := i.out.Write(translated.Bytes())
					if err != nil {
						// if this failed, push the error to the error channel so we can return
						// it to the user
						select {
						case errs <- fmt.Errorf("fail to write to output file: %v", err):
						default:
						}
						// close the context to tell other running goroutines
						// to stop
						cancel()
						// call wg.Done()
						return
					}
					translated.Reset()
				}
			}
			// some sequences left in the buffer
			if translated.Len() > 0 {
				_, err := i.out.Write(translated.Bytes())
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
	scanner := bufio.NewScanner(i.in)

	var readError error

	feeder := &fastaChannelFeeder{
		IDBuffer:       bytes.NewBuffer(make([]byte, 0)),
		CommentBuffer:  bytes.NewBuffer(make([]byte, 0)),
		SequenceBuffer: bytes.NewBuffer(make([]byte, 0)),
		FastaChan:      fnaSequences,
	}

	// tag the loop so we can break it from anywhere
Loop:
	// read the file line by line
	for scanner.Scan() {
		line := scanner.Bytes()
		// Fasta files with blanks lines are incorrect
		if len(line) == 0 {
			readError = fmt.Errorf("invalid fasta file: empty lines are forbidden")
			break Loop
		}
		// if the line starts with '>'; it's the ID of the sequence
		if line[0] == fastaID {
			// for the first sequence, the buffer for the ID is empty
			if feeder.IDBuffer.Len() > 0 {
				// if an error occurred in one of the 'inserting' goroutines,
				// break the loop
				select {
				case <-ctx.Done():
					break Loop
				default:
				}
				readError = feeder.sendFasta()
				if readError != nil {
					break Loop
				}
			}
			// parse the ID of the sequence. ID is formatted like this:
			// >sequenceID comments
			l := bytes.SplitN(line, spaceDelim, 2)
			feeder.IDBuffer.Write(l[0])
			// if there is two arrays returned, the sequence has comment
			if len(l) > 1 {
				feeder.CommentBuffer.Write(l[1])
			}
		} else {
			// if the line doesn't start with '>', then it's a part of the
			// nucleotide sequence, so write it to the buffer
			feeder.SequenceBuffer.Write(line)
		}
	}

	// don't forget to push last sequence
	readError = feeder.sendFasta()

	// if an error occured during the parsing of the fasta file,
	// return the error to trigger cancel()
	// so we can smoothly terminate all goroutines
	if readError != nil {
		return readError
	}

	// close fasta sequence channel
	close(fnaSequences)
	// wait for goroutines to finish
	wg.Wait()
	// if cancel() has been called from one of the goroutines,
	// then there must be an error in the error channel, so
	// return it
	if ctx.Err() != nil {
		return <-errs
	}
	return nil
}

// Options struct to store command line args
type Options struct {
	Required `group:"required"`
	Optional `group:"optional"`
	General  `group:"general"`
}

// Required struct to store required command line args
type Required struct {
	Sequence string `short:"s" long:"sequence" value-name:"<filename>" description:"Nucleotide sequence(s) filename"`
	Outseq   string `short:"o" long:"outseq" value-name:"<filename>" description:"Protein sequence fileName"`
}

// Optional struct to store required command line args
type Optional struct {
	Frame       string `short:"f" long:"frame" value-name:"<code>" description:"Frame to translate. Possible values:\n  [1, 2, 3, F, -1, -2, -3, R, 6]\n F: forward three frames\n R: reverse three frames\n 6: all 6 frames\n" default:"1"`
	Table       int    `short:"t" long:"table" value-name:"<code>" description:"NCBI code to use, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1 for details. Available codes: \n 0: Standard code\n 2: The Vertebrate Mitochondrial Code\n 3: The Yeast Mitochondrial Code\n 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n 5: The Invertebrate Mitochondrial Code\n 6: The Ciliate, Dasycladacean and Hexamita Nuclear Code\n 9: The Echinoderm and Flatworm Mitochondrial Code\n 10: The Euplotid Nuclear Code\n 11: The Bacterial, Archaeal and Plant Plastid Code\n 12: The Alternative Yeast Nuclear Code\n 13: The Ascidian Mitochondrial Code\n 14: The Alternative Flatworm Mitochondrial Code\n16: Chlorophycean Mitochondrial Code\n 21: Trematode Mitochondrial Code\n22: Scenedesmus obliquus Mitochondrial Code\n 23: Thraustochytrium Mitochondrial Code\n 24: Pterobranchia Mitochondrial Code\n 25: Candidate Division SR1 and Gracilibacteria Code\n 26: Pachysolen tannophilus Nuclear Code\n 29: Mesodinium Nuclear\n 30: Peritrich Nuclear\n" default:"0"`
	Clean       bool   `short:"c" long:"clean" description:"Replace stop codon '*' by 'X'"`
	Alternative bool   `short:"a" long:"alternative" description:"Define frame '-1' as using the set of codons starting with the last codon of the sequence"`
	Trim        bool   `short:"T" long:"trim" description:"Removes all 'X' and '*' characters from the right end of the translation. The trimming process starts at the end and continues until the next character is not a 'X' or a '*'"`
	NumWorker   int    `short:"n" long:"numcpu" value-name:"<n>" description:"Number of threads to use, default is number of CPU"`
}

// General struct to store required command line args
type General struct {
	Help    bool `short:"h" long:"help" description:"Show this help message"`
	Version bool `short:"v" long:"version" description:"Print the tool version and exit"`
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
	if options.Sequence == "" {
		printErrorAndExit(fmt.Errorf("missing required parameter -s | -sequence, try %s --help for details", toolName))
	}
	if options.Outseq == "" {
		printErrorAndExit(fmt.Errorf("missing required parameter -o | -outseq, try %s --help for details", toolName))
	}
	if options.Table != 0 {
		_, ok := nc.TableDiff[options.Table]
		if !ok {
			printErrorAndExit(fmt.Errorf("invalid table code: %v, try %s --help for details", options.Table, toolName))
		}
	}
	if options.NumWorker == 0 {
		options.NumWorker = runtime.NumCPU()
	}
	in, err := os.Open(options.Sequence)
	if err != nil {
		printErrorAndExit(fmt.Errorf("Could not read from input file %v: %v", options.Sequence, err))
	}
	defer in.Close()
	out, err := os.Create(options.Outseq)
	if err != nil {
		printErrorAndExit(fmt.Errorf("Could not write to output file file %v: %v", options.Outseq, err))
	}
	defer out.Close()

	ioHandler := ioHandler{
		in:  in,
		out: out,
	}

	err = ioHandler.readSequenceAndTranslate(options)
	if err != nil {
		printErrorAndExit(err)
	}
}
