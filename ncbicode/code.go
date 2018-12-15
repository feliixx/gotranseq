// Package ncbicode stores codon <-> AA
// translation.
//
// Relevant documentation:
//
//    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1
//
package ncbicode

import (
	"fmt"
	"strings"
)

var (
	standard = map[string]byte{
		"TTT": 'F',
		"TCT": 'S',
		"TAT": 'Y',
		"TGT": 'C',
		"TTC": 'F',
		"TCC": 'S',
		"TAC": 'Y',
		"TGC": 'C',
		"TTA": 'L',
		"TCA": 'S',
		"TAA": '*',
		"TGA": '*',
		"TTG": 'L',
		"TCG": 'S',
		"TAG": '*',
		"TGG": 'W',
		"CTT": 'L',
		"CCT": 'P',
		"CAT": 'H',
		"CGT": 'R',
		"CTC": 'L',
		"CCC": 'P',
		"CAC": 'H',
		"CGC": 'R',
		"CTA": 'L',
		"CCA": 'P',
		"CAA": 'Q',
		"CGA": 'R',
		"CTG": 'L',
		"CCG": 'P',
		"CAG": 'Q',
		"CGG": 'R',
		"ATT": 'I',
		"ACT": 'T',
		"AAT": 'N',
		"AGT": 'S',
		"ATC": 'I',
		"ACC": 'T',
		"AAC": 'N',
		"AGC": 'S',
		"ATA": 'I',
		"ACA": 'T',
		"AAA": 'K',
		"AGA": 'R',
		"ATG": 'M',
		"ACG": 'T',
		"AAG": 'K',
		"AGG": 'R',
		"GTT": 'V',
		"GCT": 'A',
		"GAT": 'D',
		"GGT": 'G',
		"GTC": 'V',
		"GCC": 'A',
		"GAC": 'D',
		"GGC": 'G',
		"GTA": 'V',
		"GCA": 'A',
		"GAA": 'E',
		"GGA": 'G',
		"GTG": 'V',
		"GCG": 'A',
		"GAG": 'E',
		"GGG": 'G',
	}

	vertebrateMitochondrialDiff = map[string]byte{
		"AGA": '*',
		"AGG": '*',
		"AUA": 'M',
		"UGA": 'W',
	}
	// TODO absent code ?
	yeastMitochondrialDiff = map[string]byte{
		"AUA": 'M',
		"CUU": 'T',
		"CUC": 'T',
		"CUA": 'T',
		"CUG": 'T',
		"UGA": 'W',
	}
	moldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasmaDiff = map[string]byte{
		"UGA": 'W',
	}
	invertebrateMitochondrialDiff = map[string]byte{
		"AGA": 'S',
		"AGG": 'S',
		"AUA": 'M',
		"UGA": 'W',
	}
	ciliateDasycladaceanHexamitaDiff = map[string]byte{
		"UAA": 'Q',
		"UAG": 'Q',
	}
	echinodermFlatwormMitochondrialDiff = map[string]byte{
		"AAA": 'N',
		"AGA": 'S',
		"AGG": 'S',
		"UGA": 'W',
	}

	euplotidDiff = map[string]byte{
		"UGA": 'C',
	}

	bacterialArchaealPlantPlastidDiff = map[string]byte{
		"TTT": 'F',
		"TCT": 'S',
		"TAT": 'Y',
		"TGT": 'C',
		"TTC": 'F',
		"TCC": 'S',
		"TAC": 'Y',
		"TGC": 'C',
		"TTA": 'L',
		"TCA": 'S',
		"TAA": '*',
		"TGA": '*',
		"TTG": 'L',
		"TCG": 'S',
		"TAG": '*',
		"TGG": 'W',
		"CTT": 'L',
		"CCT": 'P',
		"CAT": 'H',
		"CGT": 'R',
		"CTC": 'L',
		"CCC": 'P',
		"CAC": 'H',
		"CGC": 'R',
		"CTA": 'L',
		"CCA": 'P',
		"CAA": 'Q',
		"CGA": 'R',
		"CTG": 'L',
		"CCG": 'P',
		"CAG": 'Q',
		"CGG": 'R',
		"ATT": 'I',
		"ACT": 'T',
		"AAT": 'N',
		"AGT": 'S',
		"ATC": 'I',
		"ACC": 'T',
		"AAC": 'N',
		"AGC": 'S',
		"ATA": 'I',
		"ACA": 'T',
		"AAA": 'K',
		"AGA": 'R',
		"ATG": 'M',
		"ACG": 'T',
		"AAG": 'K',
		"AGG": 'R',
		"GTT": 'V',
		"GCT": 'A',
		"GAT": 'D',
		"GGT": 'G',
		"GTC": 'V',
		"GCC": 'A',
		"GAC": 'D',
		"GGC": 'G',
		"GTA": 'V',
		"GCA": 'A',
		"GAA": 'E',
		"GGA": 'G',
		"GTG": 'V',
		"GCG": 'A',
		"GAG": 'E',
		"GGG": 'G',
	}

	alternativeYeastDiff = map[string]byte{
		"CUG": 'S',
	}

	ascidianMitochondrialDiff = map[string]byte{
		"AGA": 'G',
		"AGG": 'G',
		"AUA": 'M',
		"UGA": 'W',
	}

	alternativeFlatwormMitochondrialDiff = map[string]byte{
		"AAA": 'N',
		"AGA": 'S',
		"AGG": 'S',
		"UAA": 'Y',
		"UGA": 'W',
	}
	chlorophyceanMitochondrialDiff = map[string]byte{
		"TAG": 'L',
	}

	trematodeMitochondrialDiff = map[string]byte{
		"TGA": 'W',
		"ATA": 'M',
		"AGA": 'S',
		"AGG": 'S',
		"AAA": 'N',
	}

	scenedesmusObliquusMitochondrialDiff = map[string]byte{
		"TCA": '*',
		"TAG": 'L',
	}

	thraustochytriumMitochondrialDiff = map[string]byte{
		"TTT": 'F',
		"TCT": 'S',
		"TAT": 'Y',
		"TGT": 'C',
		"TTC": 'F',
		"TCC": 'S',
		"TAC": 'Y',
		"TGC": 'C',
		"TTA": '*',
		"TCA": 'S',
		"TAA": '*',
		"TGA": '*',
		"TTG": 'L',
		"TCG": 'S',
		"TAG": '*',
		"TGG": 'W',
		"CTT": 'L',
		"CCT": 'P',
		"CAT": 'H',
		"CGT": 'R',
		"CTC": 'L',
		"CCC": 'P',
		"CAC": 'H',
		"CGC": 'R',
		"CTA": 'L',
		"CCA": 'P',
		"CAA": 'Q',
		"CGA": 'R',
		"CTG": 'L',
		"CCG": 'P',
		"CAG": 'Q',
		"CGG": 'R',
		"ATT": 'I',
		"ACT": 'T',
		"AAT": 'N',
		"AGT": 'S',
		"ATC": 'I',
		"ACC": 'T',
		"AAC": 'N',
		"AGC": 'S',
		"ATA": 'I',
		"ACA": 'T',
		"AAA": 'K',
		"AGA": 'R',
		"ATG": 'M',
		"ACG": 'T',
		"AAG": 'K',
		"AGG": 'R',
		"GTT": 'V',
		"GCT": 'A',
		"GAT": 'D',
		"GGT": 'G',
		"GTC": 'V',
		"GCC": 'A',
		"GAC": 'D',
		"GGC": 'G',
		"GTA": 'V',
		"GCA": 'A',
		"GAA": 'E',
		"GGA": 'G',
		"GTG": 'V',
		"GCG": 'A',
		"GAG": 'E',
		"GGG": 'G',
	}

	pterobranchiaMitochondrialDiff = map[string]byte{
		"AGA": 'S',
		"AGG": 'K',
		"UGA": 'W',
	}

	candidateDivisionSR1GracilibacteriaDiff = map[string]byte{
		"UGA": 'G',
	}

	pachysolenTannophilusDiff = map[string]byte{
		"CUG": 'A',
	}

	mesodiniumDiff = map[string]byte{
		"UAA": 'Y',
		"UAG": 'Y',
	}

	peritrichDiff = map[string]byte{
		"UAA": 'E',
		"UAG": 'E',
	}

	diffs = map[int]map[string]byte{
		VertebrateMitochondrial: vertebrateMitochondrialDiff,
		YeastMitochondrial:      yeastMitochondrialDiff,
		MoldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasma: moldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasmaDiff,
		InvertebrateMitochondrial:                                   invertebrateMitochondrialDiff,
		CiliateDasycladaceanHexamita:                                ciliateDasycladaceanHexamitaDiff,
		EchinodermFlatwormMitochondrial:                             echinodermFlatwormMitochondrialDiff,
		Euplotid:                                                    euplotidDiff,
		BacterialArchaealPlantPlastid:                               bacterialArchaealPlantPlastidDiff,
		AlternativeYeast:                                            alternativeYeastDiff,
		AscidianMitochondrial:                                       ascidianMitochondrialDiff,
		AlternativeFlatwormMitochondrial:                            alternativeFlatwormMitochondrialDiff,
		ChlorophyceanMitochondrial:                                  chlorophyceanMitochondrialDiff,
		TrematodeMitochondrial:                                      trematodeMitochondrialDiff,
		ScenedesmusObliquusMitochondrial:                            scenedesmusObliquusMitochondrialDiff,
		ThraustochytriumMitochondrial:                               thraustochytriumMitochondrialDiff,
		PterobranchiaMitochondrial:                                  pterobranchiaMitochondrialDiff,
		CandidateDivisionSR1Gracilibacteria:                         candidateDivisionSR1GracilibacteriaDiff,
		PachysolenTannophilus:                                       pachysolenTannophilusDiff,
		Mesodinium:                                                  mesodiniumDiff,
		Peritrich:                                                   peritrichDiff,
	}
)

const (
	Standard                                                    = 0
	VertebrateMitochondrial                                     = 2
	YeastMitochondrial                                          = 3
	MoldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasma = 4
	InvertebrateMitochondrial                                   = 5
	CiliateDasycladaceanHexamita                                = 6
	EchinodermFlatwormMitochondrial                             = 9
	Euplotid                                                    = 10
	BacterialArchaealPlantPlastid                               = 11
	AlternativeYeast                                            = 12
	AscidianMitochondrial                                       = 13
	AlternativeFlatwormMitochondrial                            = 14
	ChlorophyceanMitochondrial                                  = 16
	TrematodeMitochondrial                                      = 21
	ScenedesmusObliquusMitochondrial                            = 22
	ThraustochytriumMitochondrial                               = 23
	PterobranchiaMitochondrial                                  = 24
	CandidateDivisionSR1Gracilibacteria                         = 25
	PachysolenTannophilus                                       = 26
	Mesodinium                                                  = 29
	Peritrich                                                   = 30
)

// LoadTableCode returns a map of condon <-> AA
func LoadTableCode(code int) (map[string]byte, error) {

	tableCodon := map[string]byte{}
	for codon, aaCode := range standard {
		tableCodon[codon] = aaCode
	}

	if code != 0 {

		tableDiff, ok := diffs[code]
		if !ok {
			return nil, fmt.Errorf("invalid table code: %v", code)
		}

		for codon, aaCode := range tableDiff {

			codon = strings.Replace(codon, "U", "T", -1)
			tableCodon[codon] = aaCode
		}
	}
	return tableCodon, nil
}
