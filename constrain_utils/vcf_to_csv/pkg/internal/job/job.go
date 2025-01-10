package job

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"

	"github.com/brentp/vcfgo"
)

type CsvConversion struct {
	CsvPath        string
	CsvFile        *os.File
	VcfPath        string
	VcfInnerReader io.ReadCloser
	VcfReader      *vcfgo.Reader
}

func NewCsvConversion(csvPath, vcfPath string) (*CsvConversion, error) {
	csvFile, err := os.Create(csvPath)
	if err != nil {
		return nil, err
	}
	innerReader, err := NewVcfReadCloser(vcfPath)
	if err != nil {
		csvFile.Close()
		return nil, err
	}
	vcfReader, err := vcfgo.NewReader(innerReader, false)
	if err != nil {
		csvFile.Close()
		innerReader.Close()
		return nil, err
	}

	return &CsvConversion{
		CsvPath:        csvPath,
		CsvFile:        csvFile,
		VcfPath:        vcfPath,
		VcfInnerReader: innerReader,
		VcfReader:      vcfReader,
	}, nil
}

func (c *CsvConversion) Cleanup() error {
	if err := c.CsvFile.Close(); err != nil {
		return err
	}
	if err := c.VcfReader.Close(); err != nil {
		return err
	}
	return c.VcfInnerReader.Close()
}

func NewVcfReadCloser(vcfPath string) (io.ReadCloser, error) {
	file, err := os.Open(vcfPath)
	if err != nil {
		return nil, err
	}

	ext := filepath.Ext(vcfPath)
	// we just assume gzip based on extension
	// could check magic bytes instead
	if ext == ".gz" {
		gzipReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close() // Close the file if gzip.NewReader fails
			return nil, err
		}

		return &struct {
			io.Reader
			io.Closer
		}{
			Reader: bufio.NewReader(gzipReader),
			Closer: &gzipCloser{
				gzipReader: gzipReader,
				file:       file,
			},
		}, nil
	}

	// For regular files
	return &struct {
		io.Reader
		io.Closer
	}{
		Reader: bufio.NewReader(file),
		Closer: file,
	}, nil
}

type gzipCloser struct {
	gzipReader io.Closer
	file       io.Closer
}

func (g *gzipCloser) Close() error {
	if err := g.gzipReader.Close(); err != nil {
		return err
	}
	return g.file.Close()
}
