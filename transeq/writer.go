package transeq

import (
	"bytes"
	"context"
	"fmt"
	"io"
)

const (
	// suffixes ta add to sequence id for each frame
	suffixes = "123456"
	// max line size for sequence
	maxLineSize = 60
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
		buf: bytes.NewBuffer(nil),
	}
}

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

func (w *writer) addByte(b byte) {

	if w.currentLineLen == maxLineSize {
		w.newLine()
	}

	w.buf.WriteByte(b)
	w.currentLineLen++
	if b == stop || b == unknown {
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
