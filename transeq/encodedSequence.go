package transeq

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"sync"
)

// a type to hold an encoded fasta sequence
//
// s[0:4] stores the size of the sequence header (sequence id + comment) as an uint32 (little endian)
// s[4:headerSize] stores the sequence header
// s[headerSize:] stores the nucl sequence
type encodedSequence []byte

func newEncodedSequence(buf *bytes.Buffer, headerSize int) encodedSequence {

	s := getSizedSlice(4 + buf.Len())
	// reserve 4 bytes to store the header size as an uint32
	headerSize += 4
	binary.LittleEndian.PutUint32(s[0:4], uint32(headerSize))
	copy(s[4:], buf.Bytes())

	for i, n := range s[headerSize:] {
		switch n {
		case 'A':
			s[headerSize+i] = aCode
		case 'C':
			s[headerSize+i] = cCode
		case 'G':
			s[headerSize+i] = gCode
		case 'T', 'U':
			s[headerSize+i] = tCode
		case 'N':
			s[headerSize+i] = nCode
		default:
			fmt.Printf("WARNING: invalid char in sequence %s: %s, ignoring", s[headerSize:], string(s[headerSize+i]))
		}
	}
	return s
}

var pool = sync.Pool{
	New: func() interface{} {
		return make(encodedSequence, 512)
	},
}

func getSizedSlice(size int) encodedSequence {
	s := pool.Get().(encodedSequence)
	if cap(s) < size {
		s = make([]byte, size)
	}
	return s[:size]
}

func (s encodedSequence) header() []byte {
	return s[4:s.headerSize()]
}

func (s encodedSequence) headerSize() int {
	return int(binary.LittleEndian.Uint32(s[0:4]))
}

func (s encodedSequence) nuclSeqSize() int {
	return len(s) - s.headerSize()
}

func (s encodedSequence) reverseComplement() {

	headerSize := s.headerSize()
	// get the complementary sequence.
	// Basically, switch
	//   A <-> T
	//   C <-> G
	for i, n := range s[headerSize:] {
		switch n {
		case aCode:
			s[headerSize+i] = tCode
		case tCode:
			// handle both tCode and uCode
			s[headerSize+i] = aCode
		case cCode:
			s[headerSize+i] = gCode
		case gCode:
			s[headerSize+i] = cCode
		default:
			//case N -> leave it
		}
	}
	// reverse the sequence
	for i, j := headerSize, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
