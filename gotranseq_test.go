package main

import (
	"testing"
)

var (
	tbytes = [11]byte{'A', 'C', 'T', 'I', 'G', 'T', 'A', 'T', 'A', 'C', 'K'}
)

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
