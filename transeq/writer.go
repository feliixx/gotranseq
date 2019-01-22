package transeq

import (
	"bytes"
	"context"
	"fmt"
	"io"
)

const (
	mb = 1 << (10 * 2)
	// size of the buffer for writing to file
	maxBufferSize = 5 * mb
	// suffixes ta add to sequence id for each frame
	suffixes    = "123456"
	maxLineSize = 60
	stop        = '*'
	unknown     = 'X'
)

type writer struct {
	buf            []byte
	currentLineLen int
	// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
	// from right end of the sequence)
	toTrim int
}

func newWriter() *writer {
	return &writer{
		buf: make([]byte, 0, maxBufferSize),
	}
}

// sequence id should look like
// >sequenceID_<frame> comment
func (w *writer) writeHeader(seqHeader []byte, frameIndex int) {
	end := bytes.IndexByte(seqHeader, ' ')
	if end != -1 {
		w.buf = append(w.buf, seqHeader[:end]...)
		w.buf = append(w.buf, '_', suffixes[frameIndex])
		w.buf = append(w.buf, seqHeader[end:]...)
	} else {
		w.buf = append(w.buf, seqHeader...)
		w.buf = append(w.buf, '_', suffixes[frameIndex])
	}
}

func (w *writer) writeAA(aa byte) {

	if w.currentLineLen == maxLineSize {
		w.newLine()
	}
	w.buf = append(w.buf, aa)
	w.currentLineLen++

	if aa == stop || aa == unknown {
		w.toTrim++
	} else {
		w.toTrim = 0
	}
}

func (w *writer) newLine() {
	w.buf = append(w.buf, '\n')
	w.currentLineLen = 0
	w.toTrim++
}

func (w *writer) Trim() {
	w.buf = w.buf[:len(w.buf)-w.toTrim]
	w.currentLineLen -= w.toTrim
}

func (w *writer) flush(out io.Writer, cancel context.CancelFunc, errs chan error) {
	_, err := out.Write(w.buf)
	if err != nil {
		select {
		case errs <- fmt.Errorf("fail to write to output file: %v", err):
		default:
		}
		cancel()
	}
	w.buf = w.buf[0:0]
}
