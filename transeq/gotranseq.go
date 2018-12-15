package transeq

import (
	"bufio"
	"bytes"
	"context"
	"encoding/binary"
	"fmt"
	"io"
	"sync"

	"github.com/feliixx/gotranseq/ncbicode"
)

// Options struct to store command line args
type Options struct {
	Required `group:"required"`
	Optional `group:"optional"`
	General  `group:"general"`
}

// Required struct to store required command line args
type Required struct {
	Sequence string `short:"s" long:"sequence" value-name:"<filename>" description:"Nucleotide sequence(s) filename"`
	Outseq   string `short:"o" long:"outseq" value-name:"<filename>" description:"Protein sequence filename"`
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

var letterCode = map[byte]uint8{
	'A': aCode,
	'C': cCode,
	'T': tCode,
	'G': gCode,
	'N': nCode,
	'U': uCode,
}

const (
	stopByte      = '*'
	cleanStopByte = 'X'

	// uint8 code for supported nucleotides
	nCode = uint8(0)
	aCode = uint8(1)
	cCode = uint8(2)
	tCode = uint8(3)
	uCode = uint8(3)
	gCode = uint8(4)
)

// create the code map according to the selected table code
func createMapCode(code int, clean bool) (map[uint32]byte, error) {

	resultMap := map[uint32]byte{}
	twoLetterMap := map[string][]byte{}

	tmpCode := make([]uint8, 4)

	codeMap, err := ncbicode.LoadTableCode(code)
	if err != nil {
		return nil, err
	}

	for codon, aaCode := range codeMap {
		// generate 3 letter code
		for i := 0; i < 3; i++ {
			tmpCode[i] = letterCode[codon[i]]
		}
		// each codon is represented by an unique uint32:
		// each possible nucleotide is represented by an uint8 (255 possibility)
		// the three first bytes are the the code for each nucleotide
		// last byte is uint8(0)
		// example:
		// codon 'ACG' ==> uint8(1) | uint8(2) | uint8(4) | uint8(0)
		uint32Code := uint32(tmpCode[0]) | uint32(tmpCode[1])<<8 | uint32(tmpCode[2])<<16
		resultMap[uint32Code] = aaCode

		// generate 2 letter code
		codes, ok := twoLetterMap[codon[0:2]]
		if !ok {
			twoLetterMap[codon[0:2]] = []byte{aaCode}
		} else {
			twoLetterMap[codon[0:2]] = append(codes, aaCode)
		}
	}
	for twoLetterCodon, codes := range twoLetterMap {
		uniqueAA := true
		for i := 0; i < len(codes); i++ {

			if codes[i] != codes[0] {
				uniqueAA = false
			}
		}
		if uniqueAA {
			first := letterCode[twoLetterCodon[0]]
			second := letterCode[twoLetterCodon[1]]

			uint32Code := uint32(first) | uint32(second)<<8
			resultMap[uint32Code] = codes[0]
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

func computeFrames(options Options) (frames []int, reverse bool, err error) {

	frames = make([]int, 6)
	reverse = false

	switch options.Frame {
	case "1":
		frames[0] = 1
	case "2":
		frames[1] = 1
	case "3":
		frames[2] = 1
	case "F":
		for i := 0; i < 3; i++ {
			frames[i] = 1
		}
	case "-1":
		frames[3] = 1
		reverse = true
	case "-2":
		frames[4] = 1
		reverse = true
	case "-3":
		frames[5] = 1
		reverse = true
	case "R":
		for i := 3; i < 6; i++ {
			frames[i] = 1
		}
		reverse = true
	case "6":
		for i := range frames {
			frames[i] = 1
		}
		reverse = true
	default:
		err = fmt.Errorf("wrong value for -f | --frame parameter: %s", options.Frame)
	}
	return frames, reverse, err
}

const (
	endLine = '\n'
	unknown = 'X'

	// size of the buffer for writing to file
	maxBufferSize = 1024 * 1024 * 30
	// max line size for sequence
	maxLineSize = 60
	// Length of the array to store code/bytes
	// uses gCode because it's the biggest uint8 of all codes
	arrayCodeSize = (uint32(gCode) | uint32(gCode)<<8 | uint32(gCode)<<16) + 1
	// suffixes ta add to sequence id for each frame
	suffixes = "123456"
)

// Translate read a fata file, translate each sequence to the corresponding prot sequence in the specified frame
func Translate(inputSequence io.Reader, out io.Writer, options Options) error {

	mapCode, err := createMapCode(options.Table, options.Clean)
	if err != nil {
		return err
	}

	arrayCode := make([]byte, arrayCodeSize)
	for k, v := range mapCode {
		arrayCode[k] = v
	}

	framesToGenerate, reverse, err := computeFrames(options)
	if err != nil {
		return err
	}

	fnaSequences := make(chan encodedSequence, 10)
	errs := make(chan error, 1)

	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(options.NumWorker)

	for nWorker := 0; nWorker < options.NumWorker; nWorker++ {

		go func() {

			defer wg.Done()

			startPosition := make([]int, 3)
			translated := bytes.NewBuffer(nil)

			for sequence := range fnaSequences {

				select {
				case <-ctx.Done():
					return
				default:
				}

				frameIndex := 0
				startPosition[0] = 0
				startPosition[1] = 1
				startPosition[2] = 2

				idSize := int(binary.LittleEndian.Uint32(sequence[0:4]))
				nuclSeqLength := len(sequence) - idSize

			Translate:
				for _, startPos := range startPosition {

					if framesToGenerate[frameIndex] == 0 {
						frameIndex++
						continue
					}

					// sequence id should look like
					// >sequenceID_<frame> comment
					idEnd := bytes.IndexByte(sequence[4:idSize], ' ')
					if idEnd != -1 {
						translated.Write(sequence[4 : 4+idEnd])
						translated.WriteByte('_')
						translated.WriteByte(suffixes[frameIndex])
						translated.Write(sequence[4+idEnd : idSize])
					} else {
						translated.Write(sequence[4:idSize])
						translated.WriteByte('_')
						translated.WriteByte(suffixes[frameIndex])
					}
					translated.WriteByte(endLine)

					// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
					// from right end of the sequence)
					bytesToTrim := 0
					currentLength := 0

					// read the sequence 3 letters at a time, starting at a specific position
					// corresponding to the frame
					for pos := startPos + 2 + idSize; pos < len(sequence); pos += 3 {

						if currentLength == maxLineSize {
							currentLength = 0
							translated.WriteByte(endLine)
							bytesToTrim++
						}
						// create an uint32 from the codon, to retrieve the corresponding
						// AA from the map
						codonCode := uint32(sequence[pos-2]) | uint32(sequence[pos-1])<<8 | uint32(sequence[pos])<<16

						b := arrayCode[codonCode]
						if b == byte(0) {
							translated.WriteByte(unknown)
							bytesToTrim++
						} else {
							translated.WriteByte(b)
							if b == stopByte || b == cleanStopByte {
								bytesToTrim++
							} else {
								bytesToTrim = 0
							}
						}
						currentLength++
					}

					// the last codon is only 2 nucleotid long, try to guess
					// the corresponding AA
					if (nuclSeqLength-startPos)%3 == 2 {

						if currentLength == maxLineSize {
							translated.WriteByte(endLine)
							currentLength = 0
							bytesToTrim++
						}
						codonCode := uint32(sequence[len(sequence)-2]) | uint32(sequence[len(sequence)-1])<<8

						b := arrayCode[codonCode]
						if b == byte(0) {
							translated.WriteByte(unknown)
							bytesToTrim++
						} else {
							translated.WriteByte(b)
							if b == stopByte || b == cleanStopByte {
								bytesToTrim++
							} else {
								bytesToTrim = 0
							}
						}
						currentLength++
					}

					// the last codon is only 1 nucleotid long, no way to guess
					// the corresponding AA
					if (nuclSeqLength-startPos)%3 == 1 {
						if currentLength == maxLineSize {
							currentLength = 0
							translated.WriteByte(endLine)
							bytesToTrim++
						}
						translated.WriteByte(unknown)
						currentLength++
						bytesToTrim++
					}

					if options.Trim && bytesToTrim > 0 {
						// remove the last bytesToTrim bytes of the buffer
						// as they are 'X', '*' or '\n'
						translated.Truncate(translated.Len() - bytesToTrim)
						currentLength -= bytesToTrim
					}

					if currentLength != 0 {
						translated.WriteByte(endLine)
					}
					frameIndex++
				}

				if reverse && frameIndex < 6 {

					// get the complementary sequence.
					// Basically, switch
					//   A <-> T
					//   C <-> G
					// N is not modified
					for i, n := range sequence[idSize:] {

						switch n {
						case aCode:
							sequence[i+idSize] = tCode
						case tCode:
							// handle both tCode and uCode
							sequence[i+idSize] = aCode
						case cCode:
							sequence[i+idSize] = gCode
						case gCode:
							sequence[i+idSize] = cCode
						default:
							//case N -> leave it
						}
					}
					// reverse the sequence
					for i, j := idSize, len(sequence)-1; i < j; i, j = i+1, j-1 {
						sequence[i], sequence[j] = sequence[j], sequence[i]
					}

					if !options.Alternative {
						// Staden convention: Frame -1 is the reverse-complement of the sequence
						// having the same codon phase as frame 1. Frame -2 is the same phase as
						// frame 2. Frame -3 is the same phase as frame 3
						//
						// use the matrix to keep track of the forward frame as it depends on the
						// length of the sequence
						switch nuclSeqLength % 3 {
						case 0:
							startPosition[0] = 0
							startPosition[1] = 2
							startPosition[2] = 1
						case 1:
							startPosition[0] = 1
							startPosition[1] = 0
							startPosition[2] = 2
						case 2:
							startPosition[0] = 2
							startPosition[1] = 1
							startPosition[2] = 0
						}
					}
					// run the same loop, but with the reverse-complemented sequence
					goto Translate
				}

				if translated.Len() > maxBufferSize {
					_, err := out.Write(translated.Bytes())
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
				pool.Put(sequence)
			}

			if translated.Len() > 0 {
				_, err := out.Write(translated.Bytes())
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

	readSequenceFromFasta(ctx, inputSequence, fnaSequences)

	wg.Wait()

	select {
	case err, ok := <-errs:
		if ok {
			return err
		}
	default:
	}
	return nil
}

func readSequenceFromFasta(ctx context.Context, inputSequence io.Reader, fnaSequences chan encodedSequence) {

	feeder := &fastaChannelFeeder{
		idBuffer:       bytes.NewBuffer(make([]byte, 0)),
		commentBuffer:  bytes.NewBuffer(make([]byte, 0)),
		sequenceBuffer: bytes.NewBuffer(make([]byte, 0)),
		fastaChan:      fnaSequences,
	}

	scanner := bufio.NewScanner(inputSequence)

	// fasta format is:
	//
	// >sequenceID some comments on sequence
	// ACAGGCAGAGACACGACAGACGACGACACAGGAGCAGACAGCAGCAGACGACCACATATT
	// TTTGCGGTCACATGACGACTTCGGCAGCGA
	//
	// see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
	// section 1 for details
Loop:
	for scanner.Scan() {

		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {

			if feeder.idBuffer.Len() > 0 {
				select {
				case <-ctx.Done():
					break Loop
				default:
				}
				feeder.sendFasta()
			}
			feeder.reset()

			// parse the ID of the sequence. ID is formatted like this:
			// >sequenceID comments
			seqID := bytes.SplitN(line, []byte{' '}, 2)
			feeder.idBuffer.Write(seqID[0])

			if len(seqID) > 1 {
				feeder.commentBuffer.WriteByte(' ')
				feeder.commentBuffer.Write(seqID[1])
			}
		} else {
			// if the line doesn't start with '>', then it's a part of the
			// nucleotide sequence, so write it to the buffer
			feeder.sequenceBuffer.Write(line)
		}
	}

	// don't forget to push last sequence
	select {
	case <-ctx.Done():
	default:
		feeder.sendFasta()
	}

	close(fnaSequences)
}

// a type to hold an encoded fasta sequence
//
//	s[0:4] stores the size of the sequence id + the size of the comment as an uint32 (little endian)
//  s[4:idSize] stores the sequence id, and the comment id there is one
//  s[idSize:] stores the nucl sequence
type encodedSequence []byte

var pool = sync.Pool{
	New: func() interface{} {
		return make(encodedSequence, 512)
	},
}

func getSizedSlice(idSize, requiredSize int) encodedSequence {
	s := pool.Get().(encodedSequence)
	binary.LittleEndian.PutUint32(s[0:4], uint32(idSize))

	for len(s) < requiredSize {
		s = append(s, byte(0))
	}
	return s[0:requiredSize]
}

func (f *fastaChannelFeeder) sendFasta() {

	idSize := 4 + f.idBuffer.Len() + f.commentBuffer.Len()
	requiredSize := idSize + f.sequenceBuffer.Len()

	s := getSizedSlice(idSize, requiredSize)

	if f.commentBuffer.Len() > 0 {
		copy(s[idSize-f.commentBuffer.Len():idSize], f.commentBuffer.Bytes())
	}

	copy(s[4:4+f.idBuffer.Len()], f.idBuffer.Bytes())

	// convert the sequence of bytes to an array of uint8 codes,
	// so a codon (3 nucleotides | 3 bytes ) can be represented
	// as an uint32
	for i, b := range f.sequenceBuffer.Bytes() {

		switch b {
		case 'A':
			s[i+idSize] = aCode
		case 'C':
			s[i+idSize] = cCode
		case 'G':
			s[i+idSize] = gCode
		case 'T', 'U':
			s[i+idSize] = tCode
		case 'N':
			s[i+idSize] = nCode
		default:
			fmt.Printf("WARNING: invalid char in sequence %s: %s, ignoring", s[4:4+idSize], string(b))
		}
	}
	f.fastaChan <- s
}

type fastaChannelFeeder struct {
	idBuffer       *bytes.Buffer
	commentBuffer  *bytes.Buffer
	sequenceBuffer *bytes.Buffer
	fastaChan      chan encodedSequence
}

func (f *fastaChannelFeeder) reset() {
	f.idBuffer.Reset()
	f.sequenceBuffer.Reset()
	f.commentBuffer.Reset()
}
