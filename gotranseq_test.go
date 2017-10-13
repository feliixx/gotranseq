package main

import (
	"bytes"
	"os"
	"testing"

	"github.com/stretchr/testify/assert"
)

var (
	tbytes       = [11]byte{'A', 'C', 'T', 'I', 'G', 'T', 'A', 'T', 'A', 'C', 'K'}
	resultBuffer = bytes.NewBuffer(make([]byte, 0))
	iohandler    = ioHandler{out: resultBuffer}
)

func TestMain(m *testing.M) {
	i, err := os.Open("test.fna")
	if err != nil {
		os.Exit(1)
	}
	defer i.Close()
	iohandler.in = i
	os.Exit(m.Run())
}

func TestStandard6frames(t *testing.T) {
	expected := `>Contig538_YY1-1_E08_(3)_1
NFADGFRVEGAGGLVKQHRLGFHRQGTRNRHPLLLAAGEHRRX
>Contig538_YY1-1_E08_(3)_2
TSPTVSGSRALVGSSNNIAWGFIARARAIATRCCWPPESIDG
>Contig538_YY1-1_E08_(3)_3
LRRRFPGRGRWWARQTTSPGVSSPGHAQSPPAAAGRRRASTG
>Contig538_YY1-1_E08_(3)_4
PSMLSGGQQQRVAIARALAMKPQAMLFDEPTSALDPETVGEV
>Contig538_YY1-1_E08_(3)_5
PVDALRRPAAAGGDCACPGDETPGDVV*RAHQRPRPGNRRRSX
>Contig538_YY1-1_E08_(3)_6
RRCSPAASSSGWRLRVPWR*NPRRCCLTSPPAPSTRKPSAKX
`

	options := Options{
		Optional: Optional{
			Frame:       "6",
			Table:       0,
			Clean:       false,
			Alternative: false,
			Trim:        false,
			NumWorker:   1,
		},
	}

	err := iohandler.readSequenceAndTranslate(options)
	assert.Nil(t, err)
	assert.Equal(t, expected, string(resultBuffer.Bytes()), "should be equal")
}

func BenchmarkMap(b *testing.B) {
	index := 0
	fail := 0
	success := 0
	for n := 0; n < b.N; n++ {
		_, ok := letterCode[tbytes[index]]
		if !ok {
			fail++
		} else {
			success++
		}
		index++
		if index == len(tbytes) {
			index = 0
		}
	}
}

func BenchmarkSwitch(b *testing.B) {
	index := 0
	fail := 0
	success := 0
	for n := 0; n < b.N; n++ {
		switch tbytes[index] {
		case 'A':
			success++
		case 'C':
			success += 2
		case 'G':
			success += 3
		case 'T':
			success += 4
		case 'N':
			success--
		default:
			fail++
		}
		index++
		if index == len(tbytes) {
			index = 0
		}
	}
}
