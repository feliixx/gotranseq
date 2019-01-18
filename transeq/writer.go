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

type writer struct {
	buf            *bytes.Buffer
	currentLineLen int
	// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
	// from right end of the sequence)
	toTrim int
}

func newWriter() *writer {
	return &writer{
		buf: bytes.NewBuffer(make([]byte, 0, maxBufferSize)),
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
	maxLineSize = 60

	stop    = '*'
	unknown = 'X'
)

func (w *writer) writeAA(aa byte) {

	if w.currentLineLen == maxLineSize {
		w.newLine()
	}

	w.buf.WriteByte(aa)
	w.currentLineLen++

	if aa == stop || aa == unknown {
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
