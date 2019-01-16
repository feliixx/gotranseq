package transeq

import (
	"bytes"
	"context"
	"fmt"
	"io"
)

const (
	// nCode has to be 0 in order to compute two-letters code
	nCode uint8 = iota
	aCode
	cCode
	tCode
	gCode
	uCode = tCode

	// Length of the array to store code/bytes
	// uses gCode because it's the biggest uint8 of all codes
	arrayCodeSize = (uint32(gCode) | uint32(gCode)<<8 | uint32(gCode)<<16) + 1
)

var letterCode = map[byte]uint8{
	'A': aCode,
	'C': cCode,
	'T': tCode,
	'G': gCode,
	'N': nCode,
	'U': uCode,
}

type writer struct {
	buf            *bytes.Buffer
	currentLineLen int
	// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
	// from right end of the sequence)
	toTrim int
	codes  [arrayCodeSize]byte
}

func newWriter(codeMap map[string]byte, clean bool) *writer {
	return &writer{
		buf:   bytes.NewBuffer(make([]byte, 0, maxBufferSize)),
		codes: createArrayCode(codeMap, clean),
	}
}

// suffixes ta add to sequence id for each frame
const suffixes = "123456"

// sequence id should look like
// >sequenceID_<frame> comment
func (w *writer) writeID(seqHeader []byte, frameIndex int) {
	end := bytes.IndexByte(seqHeader, ' ')
	if end != -1 {
		w.buf.Write(seqHeader[:end])
		w.buf.WriteByte('_')
		w.buf.WriteByte(suffixes[frameIndex])
		w.buf.Write(seqHeader[end:])
	} else {
		w.buf.Write(seqHeader)
		w.buf.WriteByte('_')
		w.buf.WriteByte(suffixes[frameIndex])
	}
}

const (
	// max line size for sequence
	maxLineSize = 60

	stop    = '*'
	unknown = 'X'
)

func (w *writer) writeAA(codonCode uint32) {

	if w.currentLineLen == maxLineSize {
		w.newLine()
	}
	aaCode := w.codes[codonCode]
	w.buf.WriteByte(aaCode)
	w.currentLineLen++

	if aaCode == stop || aaCode == unknown {
		w.toTrim++
	} else {
		w.toTrim = 0
	}
}

func (w *writer) newLine() {
	w.buf.WriteByte('\n')
	w.currentLineLen = 0
	w.toTrim++
}

// remove the last toTrim bytes of the buffer
// as they are 'X', '*' or '\n'
func (w *writer) Trim() {
	w.buf.Truncate(w.buf.Len() - w.toTrim)
	w.currentLineLen -= w.toTrim
}

func (w *writer) flush(out io.Writer, cancel context.CancelFunc, errs chan error) {

	_, err := out.Write(w.buf.Bytes())
	if err != nil {
		select {
		case errs <- fmt.Errorf("fail to write to output file: %v", err):
		default:
		}
		cancel()
	}
	w.buf.Reset()
}

// create the code map according to the selected table code
func createArrayCode(codeMap map[string]byte, clean bool) [arrayCodeSize]byte {

	var codes [arrayCodeSize]byte
	for i := range codes {
		codes[i] = unknown
	}

	twoLetterMap := map[string][]byte{}

	for codon, aaCode := range codeMap {

		if !(clean && aaCode == stop) {

			// codon is always a 3 char string, for example 'ACG'
			// each  nucleotide of the codon is represented by an uint8
			n1, n2, n3 := letterCode[codon[0]], letterCode[codon[1]], letterCode[codon[2]]

			// convert the codon to an unique uint32:
			uint32Code := uint32(n1) | uint32(n2)<<8 | uint32(n3)<<16
			codes[uint32Code] = aaCode
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

			uint32Code := uint32(n1) | uint32(n2)<<8
			codes[uint32Code] = aaCode
		}
	}
	return codes
}
