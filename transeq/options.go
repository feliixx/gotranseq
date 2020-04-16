package transeq

// Options struct to store required command line args
type Options struct {
	Frame       string `short:"f" long:"frame" value-name:"<code>" description:"Frame to translate. Possible values:\n  [1, 2, 3, F, -1, -2, -3, R, 6]\n F: forward three frames\n R: reverse three frames\n 6: all 6 frames\n" default:"1"`
	Table       int    `short:"t" long:"table" value-name:"<code>" description:"NCBI code to use, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1 for details. Available codes: \n 0: Standard code\n 2: The Vertebrate Mitochondrial Code\n 3: The Yeast Mitochondrial Code\n 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n 5: The Invertebrate Mitochondrial Code\n 6: The Ciliate, Dasycladacean and Hexamita Nuclear Code\n 9: The Echinoderm and Flatworm Mitochondrial Code\n 10: The Euplotid Nuclear Code\n 11: The Bacterial, Archaeal and Plant Plastid Code\n 12: The Alternative Yeast Nuclear Code\n 13: The Ascidian Mitochondrial Code\n 14: The Alternative Flatworm Mitochondrial Code\n16: Chlorophycean Mitochondrial Code\n 21: Trematode Mitochondrial Code\n22: Scenedesmus obliquus Mitochondrial Code\n 23: Thraustochytrium Mitochondrial Code\n 24: Pterobranchia Mitochondrial Code\n 25: Candidate Division SR1 and Gracilibacteria Code\n 26: Pachysolen tannophilus Nuclear Code\n 29: Mesodinium Nuclear\n 30: Peritrich Nuclear\n" default:"0"`
	Clean       bool   `short:"c" long:"clean" description:"Replace stop codon '*' by 'X'"`
	Alternative bool   `short:"a" long:"alternative" description:"Define frame '-1' as using the set of codons starting with the last codon of the sequence"`
	Trim        bool   `short:"T" long:"trim" description:"Removes all 'X' and '*' characters from the right end of the translation. The trimming process starts at the end and continues until the next character is not a 'X' or a '*'"`
	NumWorker   int    `short:"n" long:"numcpu" value-name:"<n>" description:"Number of worker to use (default: number of CPU)"`
}
