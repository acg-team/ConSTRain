package csv

import (
	"encoding/csv"
	"errors"
	"fmt"
	"strconv"

	"github.com/brentp/vcfgo"
	"github.com/maverbiest/vcfconv/pkg/internal/consts"
	"github.com/maverbiest/vcfconv/pkg/internal/vcf"
)

func WriteCsv(vcfReader *vcfgo.Reader, csvWriter *csv.Writer) error {
	if len(vcfReader.Header.SampleNames) != 1 {
		return errors.New("only VCF files with one sample are supported")
	}

	csvWriter.Write(consts.GetHeader())
	for {
		variant := vcfReader.Read()
		if variant == nil {
			break
		}

		sample := variant.Samples[0]
		ft, err := vcf.GetStringFormatField(variant, sample, "FT")
		if err != nil {
			return err
		}
		if consts.IsSkipTag(ft) {
			continue
		}

		outLine := make([]string, 6)
		strId := fmt.Sprintf("%s_%d", variant.Chromosome, variant.Pos-1)
		outLine[consts.GetColIdx("str_id")] = strId

		cn, err := vcf.GetIntFormatField(variant, sample, "CN")
		if err != nil {
			continue
		}
		dp, err := vcf.GetIntFormatField(variant, sample, "DP")
		if err != nil {
			continue
		}
		dpNorm := float64(dp) / float64(cn)
		outLine[consts.GetColIdx("copy_number")] = strconv.Itoa(cn)
		outLine[consts.GetColIdx("depth")] = strconv.Itoa(dp)
		outLine[consts.GetColIdx("depth_norm")] = strconv.FormatFloat(dpNorm, 'f', -1, 64)

		freqs, err := vcf.GetStringFormatField(variant, sample, "FREQS")
		if err == nil {
			freqs = vcf.ParseFreqString(freqs)
			outLine[consts.GetColIdx("frequencies")] = freqs
		} else {
			outLine[consts.GetColIdx("frequencies")] = ""

		}

		if ft == "PASS" {
			gt, err := vcf.GetStringFormatField(variant, sample, "REPLEN")
			if err != nil {
				return err
			}
			gt = vcf.ParseGtString(gt)
			outLine[consts.GetColIdx("genotype")] = gt
		} else {
			outLine[consts.GetColIdx("genotype")] = ""
		}

		csvWriter.Write(outLine)
	}

	return csvWriter.Error()
}
