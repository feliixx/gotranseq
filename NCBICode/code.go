// Package NCBICode stores codon <-> AA
// translation.
//
// Relevant documentation:
//
//    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1
//
package NCBICode

var (
	// Standard the standard code
	Standard = map[string]byte{
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
	// ************************************
	// diff from the standard code
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

	// TableDiff store the map of diff from standard code
	TableDiff = map[int]map[string]byte{
		2:  vertebrateMitochondrialDiff,
		3:  yeastMitochondrialDiff,
		4:  moldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasmaDiff,
		5:  invertebrateMitochondrialDiff,
		6:  ciliateDasycladaceanHexamitaDiff,
		9:  echinodermFlatwormMitochondrialDiff,
		10: euplotidDiff,
		11: bacterialArchaealPlantPlastidDiff,
		12: alternativeYeastDiff,
		13: ascidianMitochondrialDiff,
		14: alternativeFlatwormMitochondrialDiff,
		16: chlorophyceanMitochondrialDiff,
		21: trematodeMitochondrialDiff,
		22: scenedesmusObliquusMitochondrialDiff,
		23: thraustochytriumMitochondrialDiff,
		24: pterobranchiaMitochondrialDiff,
		25: candidateDivisionSR1GracilibacteriaDiff,
		26: pachysolenTannophilusDiff,
		29: mesodiniumDiff,
		30: peritrichDiff,
	}
)
