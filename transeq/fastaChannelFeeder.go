package transeq

import (
	"bytes"
	"fmt"
)

type fastaChannelFeeder struct {
	idBuffer       *bytes.Buffer
	commentBuffer  *bytes.Buffer
	sequenceBuffer *bytes.Buffer
	fastaChan      chan encodedSequence
}

func newFastaChannelFeeder(fnaSequences chan encodedSequence) *fastaChannelFeeder {
	return &fastaChannelFeeder{
		idBuffer:       bytes.NewBuffer(nil),
		commentBuffer:  bytes.NewBuffer(nil),
		sequenceBuffer: bytes.NewBuffer(nil),
		fastaChan:      fnaSequences,
	}
}

func (f *fastaChannelFeeder) reset() {
	f.idBuffer.Reset()
	f.sequenceBuffer.Reset()
	f.commentBuffer.Reset()
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
