package transeq_test

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"strings"
	"testing"

	"github.com/feliixx/gotranseq/transeq"
	"github.com/jessevdk/go-flags"
	"github.com/stretchr/testify/assert"
)

type testParameters struct {
	Options  string `json:"options"`
	Expected string `json:"expected"`
}

func getOptionsAndName(opts string) (string, transeq.Options, error) {
	list := strings.Split(opts, " ")
	for i := range list {
		list[i] = strings.Join([]string{string("-"), list[i]}, "")
	}
	var options transeq.Options
	_, err := flags.ParseArgs(&options, list)
	if err != nil {
		return "", options, err
	}

	return strings.Replace(strings.Join(list, "_"), "--", "", -1), options, nil
}

func TestAllOptions(t *testing.T) {
	// read the config file to run test
	data, err := ioutil.ReadFile("testdata/data.json")
	assert.Nil(t, err)
	var param []testParameters
	err = json.Unmarshal(data, &param)
	assert.Nil(t, err)

	resultBuffer := bytes.NewBuffer(make([]byte, 0))
	iohandler := transeq.IOHandler{
		Out: resultBuffer,
	}
	// fasta sequences to translate
	fasta, err := ioutil.ReadFile("testdata/test.fna")
	assert.Nil(t, err)

	for _, p := range param {

		fastaReader := bytes.NewReader(fasta)
		iohandler.In = fastaReader

		name, opts, err := getOptionsAndName(p.Options)
		opts.NumWorker = 1
		assert.Nil(t, err)
		fmt.Printf("running test %s\n", name)
		err = iohandler.ReadSequenceAndTranslate(opts)
		assert.Nil(t, err)
		assert.Equal(t, p.Expected, string(resultBuffer.Bytes()), name)
		resultBuffer.Reset()
	}
}