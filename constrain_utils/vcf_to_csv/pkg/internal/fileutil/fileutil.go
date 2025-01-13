package fileutil

import (
	"fmt"
	"io/fs"
	"log"
	"os"
	"path/filepath"
	"strings"
)

// Get paths to all VCF files encountered in a directory.
// Optionally: recursively walk subdirectories of the specified
// directory
func VcfsFromDir(vcfDir string, recursive bool) []string {
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

func MakeOutputPaths(vcfPaths []string, csvDir string) []string {
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
