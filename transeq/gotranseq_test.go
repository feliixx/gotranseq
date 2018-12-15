package transeq_test

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"strconv"
	"strings"
	"testing"

	"github.com/feliixx/gotranseq/transeq"
	"github.com/jessevdk/go-flags"
)

func TestAllOptions(t *testing.T) {

	data, err := ioutil.ReadFile("testdata/data.json")
	if err != nil {
		t.Error(err)
	}

	var tests []struct {
		Options  string `json:"options"`
		Expected string `json:"expected"`
	}

	err = json.Unmarshal(data, &tests)
	if err != nil {
		t.Error(err)
	}

	inputBytes, err := ioutil.ReadFile("testdata/test.fna")
	if err != nil {
		t.Error(err)
	}

	out := bytes.NewBuffer(make([]byte, 0, 2*1024))

	for _, test := range tests {

		t.Run(test.Options, func(t *testing.T) {

			opts, err := getOptionsAndName(test.Options)
			opts.NumWorker = 1
			if err != nil {
				t.Error(err)
			}

			in := bytes.NewReader(inputBytes)
			err = transeq.Translate(in, out, opts)
			if err != nil {
				t.Error(err)
			}

			if want, got := test.Expected, out.String(); want != got {
				t.Fail()
				//	t.Errorf("different number of lines\nexpected %s\nbut got\n%s\n", want, got)

				//	compareByline(t, want, got)
			}

			out.Reset()
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

func compareByline(t *testing.T, want, got string) {

	wantLines := strings.Split(want, "\n")
	gotLines := strings.Split(got, "\n")

	if len(wantLines) != len(gotLines) {
		t.Errorf("different number of lines\nexpected %s\nbut got\n%s\n", want, got)
	}

	diffs := map[int]string{}

	for i := range wantLines {
		if wantLines[i] != gotLines[i] {
			diffs[i] = fmt.Sprintf("\nwant:\n%s\ngot:\n%s\n", strconv.Quote(wantLines[i]), strconv.Quote(gotLines[i]))
		}
	}

	t.Errorf("found differences in lines\n%v\n", diffs)
}
