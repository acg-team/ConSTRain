package consts

import (
	"log"
	"slices"
)

var header = []string{"str_id", "copy_number", "frequencies", "genotype", "depth", "depth_norm"}

var skipTags = []string{"UNDEF", "DPZERO", "CNZERO", "CNMISSING"}

func GetHeader() []string {
	return header
}

func GetSkipTags() []string {
	return skipTags
}

func IsSkipTag(tag string) bool {
	for _, s := range skipTags {
		if s == tag {
			return true
		}
	}
	return false
}

func GetColIdx(variable string) int {
	idx := slices.IndexFunc(header, func(s string) bool {
		return s == variable
	})

	if idx == -1 {
		log.Fatalf("variable '%s' does not exist on writer", variable)
	}

	return idx
}
