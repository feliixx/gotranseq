package transeq

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"io"
	"sync"

	"github.com/feliixx/gotranseq/ncbicode"
)

const (
	// nCode has to be 0 in order to compute two-letters code
	nCode uint8 = iota
	aCode
	cCode
	tCode
	gCode
	uCode = tCode

	// Length of the array to store codon <-> AA correspondance
	// uses gCode because it's the biggest uint8 of all codes
	arrayCodeSize = (uint32(gCode) | uint32(gCode)<<8 | uint32(gCode)<<16) + 1

	maxSeqLength = 100 * mb
)

var letterCode = map[byte]uint8{
	'A': aCode,
	'C': cCode,
	'T': tCode,
	'G': gCode,
	'N': nCode,
	'U': uCode,
}

func createCodeArray(tableCode int, clean bool) ([arrayCodeSize]byte, error) {

	var codes [arrayCodeSize]byte
	for i := range codes {
		codes[i] = unknown
	}

	twoLetterMap := map[string][]byte{}
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
		// two nucleotide, for example:
		// GTC -> 'V'
		// GTG -> 'V'
		aaCodeArray, ok := twoLetterMap[codon[:2]]
		if !ok {
			twoLetterMap[codon[:2]] = []byte{aaCode}
		} else {
			if aaCode != aaCodeArray[0] {
				twoLetterMap[codon[:2]] = append(aaCodeArray, aaCode)
			}
		}
	}

	for twoLetterCodon, aaCodeArray := range twoLetterMap {

		aaCode := aaCodeArray[0]
		if len(aaCodeArray) == 1 && !(clean && aaCode == stop) {

			n1, n2 := letterCode[twoLetterCodon[0]], letterCode[twoLetterCodon[1]]
			index := uint32(n1) | uint32(n2)<<8
			codes[index] = aaCode
		}
	}
	return codes, nil
}

func computeFrames(frameName string) (frames [6]int, reverse bool, err error) {

	var frameMap = map[string]struct {
		frames  [6]int
		reverse bool
	}{
		"1":  {[6]int{1, 0, 0, 0, 0, 0}, false},
		"2":  {[6]int{0, 1, 0, 0, 0, 0}, false},
		"3":  {[6]int{0, 0, 1, 0, 0, 0}, false},
		"F":  {[6]int{1, 1, 1, 0, 0, 0}, false},
		"-1": {[6]int{0, 0, 0, 1, 0, 0}, true},
		"-2": {[6]int{0, 0, 0, 0, 1, 0}, true},
		"-3": {[6]int{0, 0, 0, 0, 0, 1}, true},
		"R":  {[6]int{0, 0, 0, 1, 1, 1}, true},
		"6":  {[6]int{1, 1, 1, 1, 1, 1}, true},
	}

	f, ok := frameMap[frameName]
	if !ok {
		return frames, false, fmt.Errorf("wrong value for -f | --frame parameter: %s", frameName)
	}
	return f.frames, f.reverse, nil
}

// Translate read a fasta file and translate each sequence to the corresponding prot sequence
// with the specified options
func Translate(inputSequence io.Reader, out io.Writer, options Options) error {

	framesToGenerate, reverse, err := computeFrames(options.Frame)
	if err != nil {
		return err
	}

	codes, err := createCodeArray(options.Table, options.Clean)
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

			w := newWriter(codes, framesToGenerate, reverse, options.Alternative, options.Trim)

			for sequence := range fnaSequences {

				select {
				case <-ctx.Done():
					return
				default:
				}

				w.translate(sequence)

				if len(w.buf) > maxBufferSize {
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
func readSequenceFromFasta(ctx context.Context, inputSequence io.Reader, fnaSequences chan<- encodedSequence) {

	scanner := bufio.NewScanner(inputSequence)
	scanner.Buffer(make([]byte, 0, 4096), maxSeqLength)

	buf := bytes.NewBuffer(make([]byte, 0, 4096))
	headerSize := 0

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
				fnaSequences <- newEncodedSequence(buf, headerSize)
			}
			buf.Reset()
			headerSize = len(line)
		}
		buf.Write(line)
	}

	fnaSequences <- newEncodedSequence(buf, headerSize)

	close(fnaSequences)
}
