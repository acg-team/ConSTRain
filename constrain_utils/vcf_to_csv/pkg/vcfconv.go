package vcfconv

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/brentp/vcfgo"
	"github.com/maverbiest/vcfconv/pkg/internal/consts"
)

func RunFile(vcfPath string, csvPath string) {
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

func RunDir(vcfDir string, csvDir string, recursive bool) {
	vcfPaths := vcfsFromDir(vcfDir, recursive)
	if len(vcfPaths) == 0 {
		log.Fatalf("no VCF files found under --directory '%s'", vcfDir)
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

func vcfsFromDir(vcfDir string, recursive bool) []string {
	vcfPaths := make([]string, 0)
	if recursive {
		err := filepath.Walk(vcfDir, func(path string, info os.FileInfo, err error) error {
			if err != nil {
				log.Fatal(err)
			}
			if pathIsVcfFile(path, info.Mode()) {
				vcfPaths = append(vcfPaths, path)
			}

			return nil
		})
		if err != nil {
			log.Fatal(err)
		}
	} else {
		entries, err := os.ReadDir(vcfDir)
		if err != nil {
			log.Fatal(err)
		}
		for _, entry := range entries {
			if pathIsVcfFile(entry.Name(), entry.Type()) {
				vcfPaths = append(vcfPaths, fmt.Sprintf("%s%c%s", vcfDir, os.PathSeparator, entry.Name()))
			}
		}
	}

	checkDuplicatePaths(vcfPaths)

	return vcfPaths
}

func pathIsVcfFile(path string, mode fs.FileMode) bool {
	if mode.IsRegular() && (strings.HasSuffix(path, ".vcf") || strings.HasSuffix(path, ".vcf.gz")) {
		return true
	}
	return false
}

func checkDuplicatePaths(vcfPaths []string) {
	seen := make(map[string]string)
	for _, path := range vcfPaths {
		basename := vcfPathBasename(path)
		other, ok := seen[basename]
		if ok {
			log.Fatalf("VCF files '%s' and '%s' would both make CSV file %s.csv", path, other, basename)
		}
		seen[basename] = path
	}
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
		basename := vcfPathBasename(path)
		csvPaths[i] = filepath.Join(csvDir, fmt.Sprintf("%s.csv", basename))
	}

	return csvPaths
}

func vcfPathBasename(vcfPath string) string {
	basename := filepath.Base(vcfPath)
	basename = strings.TrimSuffix(basename, ".vcf.gz")
	basename = strings.TrimSuffix(basename, ".vcf")

	return basename
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
