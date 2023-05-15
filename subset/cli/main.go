package main

import (
	"bufio"
	"bytes"
	"encoding/base64"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"net/http"
	"os"
	"strings"

	"github.com/schollz/progressbar/v3"
)

const usage string = `hfsubset - Hydrofabric Subsetter

Usage:
  hfsubset [OPTIONS] identifiers...
  hfsubset (-h | --help)

Examples:
  hfsubset -l divides,nexus        \
           -o ./divides_nexus.gpkg \
           -r "pre-release"        \
           -t hl                   \
           "Gages-06752260"

  hfsubset -o ./poudre.gpkg -t hl "Gages-06752260"

Options:
`

type SubsetRequest struct {
	id      []string
	id_type *string
	layers  *string
	version *string
	output  *string
}

type SubsetResponse struct {
	data []byte
}

func (opts *SubsetRequest) Layers() []string {
	split := strings.Split(*opts.layers, ",")
	for i, v := range split {
		split[i] = strings.TrimSpace(v)
	}
	return split
}

func (opts *SubsetRequest) IDs() []string {
	for i, v := range opts.id {
		opts.id[i] = strings.TrimSpace(v)
	}
	return opts.id
}

func (opts *SubsetRequest) MarshalJSON() ([]byte, error) {
	var key string
	jsonmap := make(map[string]any)

	switch *opts.id_type {
	case "hf":
		key = "id"
		break
	case "hl":
		key = "hl_id"
		break
	case "comid":
		key = "comid"
		break
	default:
		panic("type " + *opts.id_type + " not supported; only one of: hf, hl, comid")
	}

	jsonmap["layers"] = opts.Layers()
	jsonmap[key] = opts.IDs()
	jsonmap["version"] = "pre-release"
	return json.Marshal(jsonmap)
}

func makeRequest(lambda_endpoint string, opts *SubsetRequest, bar *progressbar.ProgressBar) *SubsetResponse {
	var uri string = lambda_endpoint + "/2015-03-31/functions/function/invocations"
	payload, err := opts.MarshalJSON()
	if err != nil {
		panic(err)
	}

	reader := bytes.NewReader(payload)

	bar.Describe("[1/4] waiting for response")
	req, err := http.Post(uri, "application/json", reader)
	if err != nil {
		panic(err)
	}
	defer req.Body.Close()

	bar.Describe("[2/4] reading hydrofabric subset")
	resp := new(SubsetResponse)
	b := new(bytes.Buffer)
	buffer := bufio.NewWriter(b)
	_, err = io.Copy(buffer, req.Body)
	if err != nil {
		panic(err)
	}

	r := b.Bytes()
	if r[0] == '"' && r[len(r)-1] == '"' {
		r = r[1 : len(r)-1]
	}

	bar.Describe("[3/4] decoding base64")
	rr := bytes.NewReader(r)
	gpkg := base64.NewDecoder(base64.StdEncoding, rr)
	resp.data, err = io.ReadAll(gpkg)
	if err != nil {
		panic(err)
	}

	return resp
}

func writeToFile(request *SubsetRequest, response *SubsetResponse, bar *progressbar.ProgressBar) int {
	f, err := os.Create(*request.output)
	if err != nil {
		panic(err)
	}

	bar.Describe(fmt.Sprintf("[4/4] writing to %s", *request.output))
	w := bufio.NewWriter(f)
	mw := io.MultiWriter(w, bar)
	n, err := mw.Write(response.data)
	if err != nil {
		panic(err)
	}

	return n
}

func main() {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr, usage)
		flag.PrintDefaults()
	}

	opts := new(SubsetRequest)
	opts.id_type = flag.String("t", "hf", "One of: \"hf\", \"hl\", or \"comid\"")
	opts.layers = flag.String("l", "all", "Comma-delimited list of layers to subset.\nEither \"all\" or one or more of:\n    \"divides\", \"nexus\", \"flowpaths\",\n    \"network\", \"hydrolocations\"")
	opts.version = flag.String("r", "pre-release", "Hydrofabric version")
	opts.output = flag.String("o", "hydrofabric.gpkg", "Output file name")
	quiet := flag.Bool("quiet", false, "Disable progress bar")
	flag.Parse()

	opts.id = flag.Args()
	if len(opts.id) == 0 {
		flag.Usage()
		return
	}

	if *opts.layers == "all" {
		*opts.layers = "divides,nexus,flowpaths,network,hydrolocations"
	}

	bar := progressbar.NewOptions(3,
		progressbar.OptionSetWidth(15),
		progressbar.OptionSetDescription("[0/4] sending http request"),
		progressbar.OptionShowBytes(false),
		progressbar.OptionSetVisibility(*quiet),
	)

	var endpoint string
	if v, ok := os.LookupEnv("HFSUBSET_ENDPOINT"); ok {
		endpoint = v
	} else {
		// TODO: Change to AWS endpoint
		endpoint = "http://localhost:9000"
	}

	resp := makeRequest(endpoint, opts, bar)
	response_size := len(resp.data)
	bytes_written := writeToFile(opts, resp, bar)
	bar.Finish()
	println() // so progress bar doesn't show up

	if bytes_written != response_size {
		panic(fmt.Sprintf("wrote %d bytes out of %d bytes to %s", bytes_written, response_size, *opts.output))
	}
}
