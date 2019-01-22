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
	codes            [arrayCodeSize]byte
	buf              []byte
	currentLineLen   int
	startPos         [3]int
	frameIndex       int
	framesToGenerate [6]int
	reverse          bool
	alternative      bool
	trim             bool
	// if in trim mode, nb of bytes to trim (nb of successive 'X', '*' and '\n'
	// from right end of the sequence)
	toTrim int
}

func newWriter(codes [arrayCodeSize]byte, framesToGenerate [6]int, reverse, alternative, trim bool) *writer {
	return &writer{
		codes:            codes,
		buf:              make([]byte, 0, maxBufferSize),
		startPos:         [3]int{0, 1, 2},
		framesToGenerate: framesToGenerate,
		reverse:          reverse,
		alternative:      alternative,
		trim:             trim,
	}
}

func (w *writer) reset() {
	w.frameIndex = 0
	if w.reverse && !w.alternative {
		w.startPos[0], w.startPos[1], w.startPos[2] = 0, 1, 2
	}
}

func (w *writer) translate(sequence encodedSequence) {

	w.reset()

	w.translate3Frames(sequence)

	if w.reverse {

		if !w.alternative {
			// Staden convention: Frame -1 is the reverse-complement of the sequence
			// having the same codon phase as frame 1. Frame -2 is the same phase as
			// frame 2. Frame -3 is the same phase as frame 3
			//
			// use the matrix to keep track of the forward frame as it depends on the
			// length of the sequence
			switch sequence.nuclSeqSize() % 3 {
			case 0:
				w.startPos[0], w.startPos[1], w.startPos[2] = 0, 2, 1
			case 1:
				w.startPos[0], w.startPos[1], w.startPos[2] = 1, 0, 2
			case 2:
				w.startPos[0], w.startPos[1], w.startPos[2] = 2, 1, 0
			}
		}
		sequence.reverseComplement()
		w.translate3Frames(sequence)
	}
}

func (w *writer) translate3Frames(sequence encodedSequence) {

	for _, startPos := range w.startPos {

		if w.framesToGenerate[w.frameIndex] == 0 {
			w.frameIndex++
			continue
		}
		w.writeHeader(sequence.header())
		w.newLine()

		// read the sequence 3 letters at a time, starting at a specific position
		// corresponding to the frame
		for pos := sequence.headerSize() + startPos; pos < len(sequence)-2; pos += 3 {
			index := uint32(sequence[pos]) | uint32(sequence[pos+1])<<8 | uint32(sequence[pos+2])<<16
			w.writeAA(w.codes[index])
		}

		switch (sequence.nuclSeqSize() - startPos) % 3 {
		case 2:
			// the last codon is only 2 nucleotid long, try to guess
			// the corresponding AA
			index := uint32(sequence[len(sequence)-2]) | uint32(sequence[len(sequence)-1])<<8
			w.writeAA(w.codes[index])
		case 1:
			// the last codon is only 1 nucleotid long, no way to guess
			// the corresponding AA
			w.writeAA(unknown)
		}

		w.trimAndReturn()
		w.frameIndex++
	}
}

// sequence id should look like
// >sequenceID_<frame> comment
func (w *writer) writeHeader(seqHeader []byte) {
	end := bytes.IndexByte(seqHeader, ' ')
	if end != -1 {
		w.buf = append(w.buf, seqHeader[:end]...)
		w.buf = append(w.buf, '_', suffixes[w.frameIndex])
		w.buf = append(w.buf, seqHeader[end:]...)
	} else {
		w.buf = append(w.buf, seqHeader...)
		w.buf = append(w.buf, '_', suffixes[w.frameIndex])
	}
}

func (w *writer) writeAA(aa byte) {

	if w.currentLineLen == maxLineSize {
		w.newLine()
	}
	w.buf = append(w.buf, aa)
	w.currentLineLen++

	if w.trim {
		if aa == stop || aa == unknown {
			w.toTrim++
		} else {
			w.toTrim = 0
		}
	}
}

func (w *writer) newLine() {
	w.buf = append(w.buf, '\n')
	w.currentLineLen = 0

	if w.trim {
		w.toTrim++
	}
}

func (w *writer) trimAndReturn() {
	if w.toTrim > 0 {
		w.buf = w.buf[:len(w.buf)-w.toTrim]
		w.currentLineLen -= w.toTrim
	}

	if w.currentLineLen != 0 {
		w.newLine()
	}
	w.toTrim = 0
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
