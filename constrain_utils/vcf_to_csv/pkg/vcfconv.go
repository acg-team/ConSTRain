package vcfconv

import (
	"encoding/csv"
	"errors"
	"fmt"
	"io/fs"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/brentp/vcfgo"
	"github.com/maverbiest/vcfconv/pkg/internal/consts"
	"github.com/maverbiest/vcfconv/pkg/internal/job"
)

func RunFile(vcfPath string, csvPath string) {
	job, err := job.NewCsvConversion(csvPath, vcfPath)
	if err != nil {
		log.Fatalf("error setting up files: %s", err)
	}
	defer job.Cleanup()

	fmt.Println("Writing variants from", vcfPath, "to", csvPath)

	csvWriter := csv.NewWriter(job.CsvFile)
	defer csvWriter.Flush()

	if err := writeCsv(job.VcfReader, csvWriter); err != nil {
		log.Fatal(err)
	}
}

func RunDir(vcfDir string, csvDir string, nWorkers int, recursive bool) {
	vcfPaths := vcfsFromDir(vcfDir, recursive)
	if len(vcfPaths) == 0 {
		log.Fatalf("no VCF files found under --directory '%s'", vcfDir)
	}

	var wg sync.WaitGroup
	wg.Add(len(vcfPaths))

	jobs := make(chan *job.CsvConversion, len(vcfPaths))
	for i := 1; i <= nWorkers; i++ {
		go worker(jobs, &wg)
	}

	for i, csvPath := range makeOutputPaths(vcfPaths, csvDir) {
		vcfPath := &vcfPaths[i]
		job, err := job.NewCsvConversion(csvPath, *vcfPath)
		if err != nil {
			log.Fatalf("error setting up files: %s", err)
		}

		jobs <- job
	}
	close(jobs)

	wg.Wait()
}

func worker(jobs <-chan *job.CsvConversion, wg *sync.WaitGroup) {
	for j := range jobs {
		fmt.Println("Creating output file", j.CsvPath)
		csvWriter := csv.NewWriter(j.CsvFile)
		if err := writeCsv(j.VcfReader, csvWriter); err != nil {
			log.Fatalf("error writing file %s: %s", j.CsvPath, err)
		}

		j.Cleanup()
		csvWriter.Flush()
		wg.Done()
	}
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
		cn, err := getIntFormatField(variant, sample, "CN")
		if err != nil {
			continue
		}
		dp, err := getIntFormatField(variant, sample, "DP")
		if err != nil {
			continue
		}
		dpNorm := float64(dp) / float64(cn)
		outLine[consts.GetColIdx("copy_number")] = strconv.Itoa(cn)
		outLine[consts.GetColIdx("depth")] = strconv.Itoa(dp)
		outLine[consts.GetColIdx("depth_norm")] = strconv.FormatFloat(dpNorm, 'f', -1, 64)

		freqs, err := getStringFormatField(variant, sample, "FREQS")
		if err != nil {
			continue
		}
		freqs = parseFreqString(freqs)
		outLine[consts.GetColIdx("frequencies")] = freqs

		ft, err := getStringFormatField(variant, sample, "FT")
		if err != nil {
			return err
		}
		if ft == "PASS" {
			gt, err := getStringFormatField(variant, sample, "REPLEN")
			if err != nil {
				return err
			}
			gt = parseGtString(gt)
			outLine[consts.GetColIdx("genotype")] = gt
		}

		csvWriter.Write(outLine)
	}

	return csvWriter.Error()
}

func getIntFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) (int, error) {
	res, err := variant.GetGenotypeField(sample, tag, -1)
	if err != nil {
		return 0, err
		// log.Fatal(err)
	}
	return res.(int), nil
}

func getStringFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) (string, error) {
	res, err := variant.GetGenotypeField(sample, tag, "")
	if err != nil {
		return "", err
	}
	return res.(string), nil
}

func parseFreqString(freqs string) string {
	if freqs == "" {
		return ""
	}
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
