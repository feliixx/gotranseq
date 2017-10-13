# Gotranseq

Translate nucleic acid sequences to their corresponding peptide sequences. 
Like EMBOSS transeq, but written in go 

## Purpose 

EMBOSS transeq is a great tool, but can be quite painfull for some use cases, 
because it silently truncate the sequence ID if it contains chars like **':'**, 
or rename the sequence ID if it contains chars like **'|'**

This tool is an attempt to solve this porblem. It's also way faster than EMBOSS
transeq because it can be parrallelized: 

benchmark on ubuntu 16.04, machine with 2 CPU Intel(R) Core(TM)2 Duo CPU 3.00GHz
with a 189MB fasta file: 

```
#EMBOSS transeq
time transeq -sequence file.fna -outseq out.faa -frame 6  
41,82s user 0,76s system 85% cpu 49,696 total

#gotranseq
time ./gotranseq --sequence file.fna --outseq out.faa --frame 6
18,30s user 1,03s system 123% cpu 15,628 total
```

## Installation

download the latest binary from the [release page](https://github.com/feliixx/gotranseq/releases) 

or install with **go get** (require to have go installed on your machine)

```
go get -u "github.com/feliixx/gotranseq"
```

## Usage 

use `gotranseq --help` to print the help: 

```
gotranseq version 0.1

Usage:
  gotranseq

required:
  -s, --sequence=<filename>    Nucleotide sequence(s) filename
  -o, --outseq=<filename>      Protein sequence filename

optional:
  -f, --frame=<code>           Frame to translate. Possible values:
                               [1, 2, 3, F, -1, -2, -3, R, 6]
                               F: forward three frames
                               R: reverse three frames
                               6: all 6 frames
                               (default: 1)
  -t, --table=<code>           NCBI code to use, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1 for
                               details. Available codes:
                               0: Standard code
                               2: The Vertebrate Mitochondrial Code
                               3: The Yeast Mitochondrial Code
                               4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
                               5: The Invertebrate Mitochondrial Code
                               6: The Ciliate, Dasycladacean and Hexamita Nuclear Code
                               9: The Echinoderm and Flatworm Mitochondrial Code
                               10: The Euplotid Nuclear Code
                               11: The Bacterial, Archaeal and Plant Plastid Code
                               12: The Alternative Yeast Nuclear Code
                               13: The Ascidian Mitochondrial Code
                               14: The Alternative Flatworm Mitochondrial Code
                               16: Chlorophycean Mitochondrial Code
                               21: Trematode Mitochondrial Code
                               22: Scenedesmus obliquus Mitochondrial Code
                               23: Thraustochytrium Mitochondrial Code
                               24: Pterobranchia Mitochondrial Code
                               25: Candidate Division SR1 and Gracilibacteria Code
                               26: Pachysolen tannophilus Nuclear Code
                               29: Mesodinium Nuclear
                               30: Peritrich Nuclear
                               (default: 0)
  -c, --clean                  Replace stop codon '*' by 'X'
  -a, --alternative            Define frame '-1' as using the set of codons starting with the last codon of the sequence
  -T, --trim                   Removes all 'X' and '*' characters from the right end of the translation. The trimming process starts at the
                               end and continues until the next character is not a 'X' or a '*'
  -n, --numcpu=<n>             Number of threads to use, default is number of CPU

general:
  -h, --help                   Show this help message
  -v, --version                Print the tool version and exit
```