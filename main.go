package main

import (
	"fmt"
	"os"
	"runtime"

	"github.com/feliixx/gotranseq/transeq"
	"github.com/jessevdk/go-flags"
)

const toolName = "gotranseq"

// Version of gotranseq. Should be linked via ld_flags when compiling
// use this to set version to last known tag:
//
//  go build -ldflags "-X main.Version=$(git describe --tags $(git rev-list --tags --max-count=1))"
//
// use debug.ReadBuildInfo() when https://github.com/golang/go/issues/37475 is implemented
var Version string

// GlobalOptions struct to store command line args
type GlobalOptions struct {
	Required        `group:"required"`
	transeq.Options `group:"optional"`
	General         `group:"general"`
}

// Required struct to store required command line args
type Required struct {
	Sequence string `short:"s" long:"sequence" value-name:"<filename>" description:"Nucleotide sequence(s) filename"`
	Outseq   string `short:"o" long:"outseq" value-name:"<filename>" description:"Protein sequence filename"`
}

// General struct to store required command line args
type General struct {
	Help    bool `short:"h" long:"help" description:"Show this help message"`
	Version bool `short:"v" long:"version" description:"Print the tool version and exit"`
}

func run(options GlobalOptions) error {

	if options.Sequence == "" {
		return fmt.Errorf("missing required parameter -s | -sequence, try %s --help for details", toolName)
	}
	if options.Outseq == "" {
		return fmt.Errorf("missing required parameter -o | -outseq, try %s --help for details", toolName)
	}

	if options.NumWorker == 0 {
		options.NumWorker = runtime.NumCPU()
	}

	in, err := os.Open(options.Sequence)
	if err != nil {
		return err
	}
	defer in.Close()

	out, err := os.Create(options.Outseq)
	if err != nil {
		return err
	}
	defer out.Close()

	return transeq.Translate(in, out, options.Options)
}

func main() {

	var options GlobalOptions
	p := flags.NewParser(&options, flags.Default&^flags.HelpFlag)
	p.Usage = "--sequence file.fna --outseq out.faa"
	_, err := p.Parse()
	if err != nil {
		fmt.Printf("wrong arguments: %v, try %s --help for more informations\n", err, toolName)
		os.Exit(1)
	}
	if options.Help {
		fmt.Printf("%s %s\n\n", toolName, Version)
		p.WriteHelp(os.Stdout)
		os.Exit(0)
	}
	if options.Version {
		fmt.Printf("%s %s\n", toolName, Version)
		os.Exit(0)
	}

	err = run(options)
	if err != nil {
		fmt.Printf("fail to translate file:\n%v", err)
	}
}
