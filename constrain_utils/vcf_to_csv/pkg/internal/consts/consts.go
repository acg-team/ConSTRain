package consts

import (
	"log"
	"slices"
)

var header = []string{"str_id", "copy_number", "frequencies", "genotype", "depth", "depth_norm"}

func GetHeader() []string {
	return header
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
