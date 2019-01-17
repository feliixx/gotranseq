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

// Options struct to store required command line args
type Options struct {
	Frame       string `short:"f" long:"frame" value-name:"<code>" description:"Frame to translate. Possible values:\n  [1, 2, 3, F, -1, -2, -3, R, 6]\n F: forward three frames\n R: reverse three frames\n 6: all 6 frames\n" default:"1"`
	Table       int    `short:"t" long:"table" value-name:"<code>" description:"NCBI code to use, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1 for details. Available codes: \n 0: Standard code\n 2: The Vertebrate Mitochondrial Code\n 3: The Yeast Mitochondrial Code\n 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n 5: The Invertebrate Mitochondrial Code\n 6: The Ciliate, Dasycladacean and Hexamita Nuclear Code\n 9: The Echinoderm and Flatworm Mitochondrial Code\n 10: The Euplotid Nuclear Code\n 11: The Bacterial, Archaeal and Plant Plastid Code\n 12: The Alternative Yeast Nuclear Code\n 13: The Ascidian Mitochondrial Code\n 14: The Alternative Flatworm Mitochondrial Code\n16: Chlorophycean Mitochondrial Code\n 21: Trematode Mitochondrial Code\n22: Scenedesmus obliquus Mitochondrial Code\n 23: Thraustochytrium Mitochondrial Code\n 24: Pterobranchia Mitochondrial Code\n 25: Candidate Division SR1 and Gracilibacteria Code\n 26: Pachysolen tannophilus Nuclear Code\n 29: Mesodinium Nuclear\n 30: Peritrich Nuclear\n" default:"0"`
	Clean       bool   `short:"c" long:"clean" description:"Replace stop codon '*' by 'X'"`
	Alternative bool   `short:"a" long:"alternative" description:"Define frame '-1' as using the set of codons starting with the last codon of the sequence"`
	Trim        bool   `short:"T" long:"trim" description:"Removes all 'X' and '*' characters from the right end of the translation. The trimming process starts at the end and continues until the next character is not a 'X' or a '*'"`
	NumWorker   int    `short:"n" long:"numcpu" value-name:"<n>" description:"Number of threads to use, default is number of CPU"`
}

var letterCode = map[byte]uint8{
	'A': aCode,
	'C': cCode,
	'T': tCode,
	'G': gCode,
	'N': nCode,
	'U': uCode,
}

// create the code map according to the selected table code
func createArrayCode(tableCode int, clean bool) ([arrayCodeSize]byte, error) {

	var twoLetterMap = map[string][]byte{}
	var codes [arrayCodeSize]byte
	for i := range codes {
		codes[i] = unknown
	}

	codeMap, err := ncbicode.LoadTableCode(tableCode)
	if err != nil {
		return codes, err
	}

	for codon, aaCode := range codeMap {

		if !(clean && aaCode == stop) {
			// codon is always a 3 char string, for example 'ACG'
			// each  nucleotide of the codon is represented by an uint8
			n1, n2, n3 := letterCode[codon[0]], letterCode[codon[1]], letterCode[codon[2]]
			index := uint32(n1) | uint32(n2)<<8 | uint32(n3)<<16
			codes[index] = aaCode
		}

		// in some case, all codon for an AA will start with the same
		// two nucleotid
		// for example:
		// GTC -> 'V'
		// GTG -> 'V'
		aaCodeArray, ok := twoLetterMap[codon[:2]]
		if !ok {
			twoLetterMap[codon[:2]] = []byte{aaCode}
		} else {
			twoLetterMap[codon[:2]] = append(aaCodeArray, aaCode)
		}
	}

	for twoLetterCodon, aaCodeArray := range twoLetterMap {

		aaCode := aaCodeArray[0]
		uniqueAA := true
		for _, c := range aaCodeArray {
			if c != aaCode {
				uniqueAA = false
				break
			}
		}
		if uniqueAA {

			if clean && aaCode == stop {
				continue
			}
			n1, n2 := letterCode[twoLetterCodon[0]], letterCode[twoLetterCodon[1]]
			index := uint32(n1) | uint32(n2)<<8
			codes[index] = aaCode
		}
	}
	return codes, nil
}

var frameMap = map[string][6]int{
	"1":  {1, 0, 0, 0, 0, 0},
	"2":  {0, 1, 0, 0, 0, 0},
	"3":  {0, 0, 1, 0, 0, 0},
	"F":  {1, 1, 1, 0, 0, 0},
	"-1": {0, 0, 0, 1, 0, 0},
	"-2": {0, 0, 0, 0, 1, 0},
	"-3": {0, 0, 0, 0, 0, 1},
	"R":  {0, 0, 0, 1, 1, 1},
	"6":  {1, 1, 1, 1, 1, 1},
}

func computeFrames(frameName string) (frames [6]int, reverse bool, err error) {

	frames, ok := frameMap[frameName]
	if !ok {
		return frames, reverse, fmt.Errorf("wrong value for -f | --frame parameter: %s", frameName)
	}

	switch frameName {
	case "-1", "-2", "-3", "R", "6":
		reverse = true
	}
	return frames, reverse, nil
}

const (
	mb = 1 << (10 * 2)
	// size of the buffer for writing to file
	maxBufferSize = 5 * mb
)

// Translate read a fata file, translate each sequence to the corresponding prot sequence in the specified frame
func Translate(inputSequence io.Reader, out io.Writer, options Options) error {

	framesToGenerate, reverse, err := computeFrames(options.Frame)
	if err != nil {
		return err
	}

	codes, err := createArrayCode(options.Table, options.Clean)
	if err != nil {
		return err
	}

	fnaSequences := make(chan encodedSequence, 100)
	errs := make(chan error, 1)

	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(options.NumWorker)

	for nWorker := 0; nWorker < options.NumWorker; nWorker++ {

		go func() {

			defer wg.Done()

			var startPos [3]int
			w := newWriter()

			for sequence := range fnaSequences {

				select {
				case <-ctx.Done():
					return
				default:
				}

				frameIndex := 0
				startPos[0], startPos[1], startPos[2] = 0, 1, 2

				idSize := int(binary.LittleEndian.Uint32(sequence[0:4]))
				nuclSeqLength := len(sequence) - idSize

			Translate:
				for _, startPos := range startPos {

					if framesToGenerate[frameIndex] == 0 {
						frameIndex++
						continue
					}

					w.writeID(sequence[4:idSize], frameIndex)
					w.newLine()
					w.toTrim = 0

					// read the sequence 3 letters at a time, starting at a specific position
					// corresponding to the frame
					for pos := startPos + 2 + idSize; pos < len(sequence); pos += 3 {
						index := uint32(sequence[pos-2]) | uint32(sequence[pos-1])<<8 | uint32(sequence[pos])<<16
						w.writeAA(codes[index])
					}

					// the last codon is only 2 nucleotid long, try to guess
					// the corresponding AA
					if (nuclSeqLength-startPos)%3 == 2 {
						index := uint32(sequence[len(sequence)-2]) | uint32(sequence[len(sequence)-1])<<8
						w.writeAA(codes[index])
					}

					// the last codon is only 1 nucleotid long, no way to guess
					// the corresponding AA
					if (nuclSeqLength-startPos)%3 == 1 {
						w.writeAA(unknown)
					}

					if options.Trim && w.toTrim > 0 {
						w.Trim()
					}

					if w.currentLineLen != 0 {
						w.newLine()
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
							startPos[0], startPos[1], startPos[2] = 0, 2, 1
						case 1:
							startPos[0], startPos[1], startPos[2] = 1, 0, 2
						case 2:
							startPos[0], startPos[1], startPos[2] = 2, 1, 0
						}
					}
					// run the same loop, but with the reverse-complemented sequence
					goto Translate
				}

				if w.buf.Len() > maxBufferSize {
					w.flush(out, cancel, errs)
				}
				pool.Put(sequence)
			}
			w.flush(out, cancel, errs)
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

// fasta format is:
//
// >sequenceID some comments on sequence
// ACAGGCAGAGACACGACAGACGACGACACAGGAGCAGACAGCAGCAGACGACCACATATT
// TTTGCGGTCACATGACGACTTCGGCAGCGA
//
// see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
// section 1 for details
func readSequenceFromFasta(ctx context.Context, inputSequence io.Reader, fnaSequences chan encodedSequence) {

	idSize := 0
	buf := bytes.NewBuffer(make([]byte, 0, 512))
	scanner := bufio.NewScanner(inputSequence)
Loop:
	for scanner.Scan() {

		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {

			if buf.Len() > 0 {
				select {
				case <-ctx.Done():
					break Loop
				default:
				}
				fnaSequences <- encodeSequence(buf, idSize)
			}
			idSize = len(line)
		}
		// if the line doesn't start with '>', then it's a part of the
		// nucleotide sequence, so write it to the buffer
		buf.Write(line)
	}
	fnaSequences <- encodeSequence(buf, idSize)

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

func getSizedSlice(size int) encodedSequence {

	s := pool.Get().(encodedSequence)
	for len(s) < size {
		s = append(s, byte(0))
	}
	return s[:size]
}

func encodeSequence(buf *bytes.Buffer, idSize int) encodedSequence {

	s := getSizedSlice(4 + buf.Len())

	binary.LittleEndian.PutUint32(s[0:4], uint32(4+idSize))
	copy(s[4:], buf.Bytes())
	buf.Reset()

	// convert the sequence of bytes to an array of uint8 codes,
	// so a codon (3 nucleotides | 3 bytes ) can be represented
	// as an uint32
	for i, b := range s[4+idSize:] {

		switch b {
		case 'A':
			s[4+idSize+i] = aCode
		case 'C':
			s[4+idSize+i] = cCode
		case 'G':
			s[4+idSize+i] = gCode
		case 'T', 'U':
			s[4+idSize+i] = tCode
		case 'N':
			s[4+idSize+i] = nCode
		default:
			fmt.Printf("WARNING: invalid char in sequence %s: %s, ignoring", s[4+idSize:], string(s[i+4+idSize]))
		}
	}
	return s
}
