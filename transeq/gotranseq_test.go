package transeq_test

import (
	"bytes"
	"encoding/json"
	"io/ioutil"
	"strings"
	"testing"

	"github.com/feliixx/gotranseq/transeq"
	"github.com/jessevdk/go-flags"
)

type testParameters struct {
	Options  string `json:"options"`
	Expected string `json:"expected"`
}

func TestAllOptions(t *testing.T) {

	data, err := ioutil.ReadFile("testdata/data.json")
	if err != nil {
		t.Error(err)
	}

	var param []testParameters
	err = json.Unmarshal(data, &param)
	if err != nil {
		t.Error(err)
	}

	resultBuffer := bytes.NewBuffer(make([]byte, 0))
	iohandler := transeq.IOHandler{
		Out: resultBuffer,
	}

	fasta, err := ioutil.ReadFile("testdata/test.fna")
	if err != nil {
		t.Error(err)
	}

	for _, p := range param {

		t.Run(p.Options, func(t *testing.T) {

			fastaReader := bytes.NewReader(fasta)
			iohandler.In = fastaReader

			opts, err := getOptionsAndName(p.Options)
			opts.NumWorker = 1
			if err != nil {
				t.Error(err)
			}

			err = iohandler.ReadSequenceAndTranslate(opts)
			if err != nil {
				t.Error(err)
			}

			if want, got := p.Expected, resultBuffer.String(); want != got {
				t.Errorf("expected %s\nbut got\n%s\n", want, got)
			}

			resultBuffer.Reset()
		})
	}
}

func getOptionsAndName(opts string) (options transeq.Options, err error) {

	// convert emboss transeq flag (single '-' prefix) to gotranseq flags
	flagOpts := strings.Split(opts, " ")
	for i, flag := range flagOpts {
		if strings.HasPrefix(flag, "-") {
			flagOpts[i] = strings.Replace(flag, "-", "--", 1)
		}
	}
	_, err = flags.ParseArgs(&options, flagOpts)

	return options, err
}
