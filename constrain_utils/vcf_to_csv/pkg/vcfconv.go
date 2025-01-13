package vcfconv

import (
	"encoding/csv"
	"fmt"
	"log"
	"sync"

	csvLocal "github.com/maverbiest/vcfconv/pkg/internal/csv"
	"github.com/maverbiest/vcfconv/pkg/internal/fileutil"
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

	if err := csvLocal.WriteCsv(job.VcfReader, csvWriter); err != nil {
		log.Fatal(err)
	}
}

func RunDir(vcfDir string, csvDir string, nWorkers int, recursive bool) {
	vcfPaths := fileutil.VcfsFromDir(vcfDir, recursive)
	if len(vcfPaths) == 0 {
		log.Fatalf("no VCF files found under --directory '%s'", vcfDir)
	}

	var wg sync.WaitGroup
	wg.Add(len(vcfPaths))

	jobs := make(chan *job.CsvConversion, len(vcfPaths))
	for i := 1; i <= nWorkers; i++ {
		go worker(jobs, &wg)
	}

	for i, csvPath := range fileutil.MakeOutputPaths(vcfPaths, csvDir) {
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
		if err := csvLocal.WriteCsv(j.VcfReader, csvWriter); err != nil {
			log.Fatalf("error writing file %s: %s", j.CsvPath, err)
		}

		csvWriter.Flush()
		j.Cleanup()
		wg.Done()
	}
}
