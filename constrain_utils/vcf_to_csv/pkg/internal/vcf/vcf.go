package vcf

import (
	"fmt"
	"strings"

	"github.com/brentp/vcfgo"
)

func GetIntFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) (int, error) {
	res, err := variant.GetGenotypeField(sample, tag, -1)
	if err != nil {
		return 0, err
	}
	return res.(int), nil
}

func GetStringFormatField(variant *vcfgo.Variant, sample *vcfgo.SampleGenotype, tag string) (string, error) {
	res, err := variant.GetGenotypeField(sample, tag, "")
	if err != nil {
		return "", err
	}
	return res.(string), nil
}

func ParseFreqString(freqs string) string {
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

func ParseGtString(gt string) string {
	return fmt.Sprintf("[%s]", gt)
}
