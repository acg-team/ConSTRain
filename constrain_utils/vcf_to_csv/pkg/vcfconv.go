package vcfconv

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/brentp/vcfgo"
	"github.com/maverbiest/vcfconv/pkg/internal/consts"
)

func RunFiles(vcfPath string, csvPath string) {
	vcfReadCloser, csvFile := initFiles(vcfPath, csvPath)
	defer vcfReadCloser.Close()
	defer csvFile.Close()

	fmt.Println("Writing variants from", vcfPath, "to", csvPath)

	vcfReader, err := vcfgo.NewReader(vcfReadCloser, false)
	if err != nil {
		log.Fatal(err)
	}
	csvWriter := csv.NewWriter(csvFile)
	defer csvWriter.Flush()

	if err := writeCsv(vcfReader, csvWriter); err != nil {
		log.Fatal(err)
	}
}

func RunDirs(vcfDir string, csvDir string) {
	vcfPaths := vcfsFromDir(vcfDir)
	if len(vcfPaths) == 0 {
		return
	}
	var wg sync.WaitGroup
	wg.Add(len(vcfPaths))

	for i, csvPath := range makeOutputPaths(vcfPaths, csvDir) {
		go func() {
			defer wg.Done()

			vcfPath := &vcfPaths[i]
			vcfReadCloser, csvFile := initFiles(*vcfPath, csvPath)
			defer vcfReadCloser.Close()
			defer csvFile.Close()

			fmt.Println("Writing variants from", *vcfPath, "to", csvPath)

			vcfReader, err := vcfgo.NewReader(vcfReadCloser, false)
			if err != nil {
				log.Fatal(err)
			}
			defer vcfReader.Close()

			csvWriter := csv.NewWriter(csvFile)
			defer csvWriter.Flush()

			if err := writeCsv(vcfReader, csvWriter); err != nil {
				log.Fatal(err)
			}

			fmt.Println("Finished ", *vcfPath)
		}()
	}

	wg.Wait()
}

type vcfGzipCloser struct {
	gzipReader io.Closer
	file       io.Closer
}

func (c *vcfGzipCloser) Close() error {
	if err := c.gzipReader.Close(); err != nil {
		return err
	}
	return c.file.Close()
}

func ReaderWithCloser(filePath string) (io.ReadCloser, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}

	ext := filepath.Ext(filePath)
	if ext == ".gz" {
		gzReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close() // Close the file if gzip.NewReader fails
			return nil, err
		}

		return &struct {
			io.Reader
			io.Closer
		}{
			Reader: bufio.NewReader(gzReader),
			Closer: &vcfGzipCloser{
				gzipReader: gzReader,
				file:       file,
			},
		}, nil
	}

	// For regular files
	return struct {
		io.Reader
		io.Closer
	}{
		Reader: bufio.NewReader(file),
		Closer: file,
	}, nil
}

func initFiles(vcfPath string, csvPath string) (io.ReadCloser, *os.File) {
	vcfReader, err := ReaderWithCloser(vcfPath)
	if err != nil {
		log.Fatal(err)
	}

	csvFile, err := os.Create(csvPath)
	if err != nil {
		log.Fatal(err)
	}

	return vcfReader, csvFile
}

func vcfsFromDir(vcfDir string) []string {
	entries, err := os.ReadDir(vcfDir)
	if err != nil {
		log.Fatal(err)
	}

	vcfPaths := make([]string, 0)
	for _, entry := range entries {
		if entry.Type().IsRegular() && (strings.HasSuffix(entry.Name(), ".vcf") || strings.HasSuffix(entry.Name(), ".vcf.gz")) {
			vcfPaths = append(vcfPaths, fmt.Sprintf("%s%c%s", vcfDir, os.PathSeparator, entry.Name()))
		}
	}

	return vcfPaths
}

func makeOutputPaths(vcfPaths []string, csvDir string) []string {
	fileInfo, err := os.Stat(csvDir)
	if err != nil {
		log.Fatal(err)
	} else if !fileInfo.IsDir() {
		log.Fatalf("cannot use path '%s' as --outdir, not a directory", csvDir)
	}

	csvPaths := make([]string, len(vcfPaths))
	for i, path := range vcfPaths {
		basename := filepath.Base(path)
		basename = strings.TrimSuffix(basename, ".vcf.gz")
		basename = strings.TrimSuffix(basename, ".vcf")
		csvPaths[i] = fmt.Sprintf("%s%c%s%s", csvDir, os.PathSeparator, basename, ".csv")
	}

	return csvPaths
}

func writeCsv(vcfReader *vcfgo.Reader, csvWriter *csv.Writer) error {
	if len(vcfReader.Header.SampleNames) != 1 {
		return errors.New("only VCF files with one sample are supported")
	}

	csvWriter.Write(consts.GetHeader())
	for {
		variant := vcfReader.Read()
		if variant == nil {
			break
		}
		outLine := make([]string, 6)

		strId := fmt.Sprintf("%s_%d", variant.Chromosome, variant.Pos-1)
		outLine[consts.GetColIdx("str_id")] = strId

		sample := variant.Samples[0]
		cn := getIntFormatField(variant, sample, "CN")
		dp := getIntFormatField(variant, sample, "DP")
		dpNorm := float64(dp) / float64(cn)
		outLine[consts.GetColIdx("copy_number")] = strconv.Itoa(cn)
		outLine[consts.GetColIdx("depth")] = strconv.Itoa(dp)
		outLine[consts.GetColIdx("depth_norm")] = strconv.FormatFloat(dpNorm, 'f', -1, 64)

		freqs := getStringFormatField(variant, sample, "FREQS")
		if freqs != "" {
			freqs = parseFreqString(freqs)
		}
		outLine[consts.GetColIdx("frequencies")] = freqs

		ft := getStringFormatField(variant, sample, "FT")
		if ft == "PASS" {
			gt := getStringFormatField(variant, sample, "REPLEN")
			gt = parseGtString(gt)
			outLine[consts.GetColIdx("genotype")] = gt
		}

		csvWriter.Write(outLine)
	}

	return csvWriter.Error()
}

func getIntFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) int {
	res, err := variant.GetGenotypeField(sample, tag, -1)
	if err != nil {
		log.Fatal(err)
	}
	return res.(int)
}

func getStringFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) string {
	res, err := variant.GetGenotypeField(sample, tag, "")
	if err != nil {
		return ""
	}
	return res.(string)
}

func parseFreqString(freqs string) string {
	res := ""
	for _, item := range strings.Split(freqs, "|") {
		split := strings.Split(item, ",")
		res = fmt.Sprintf("%s%s: %s,", res, split[0], split[1])
	}
	return fmt.Sprintf("{%s}", res)
}

func parseGtString(gt string) string {
	return fmt.Sprintf("[%s]", gt)
}
