package cmd

import (
	"log"
	"runtime"

	vcfconv "github.com/maverbiest/vcfconv/pkg"
	"github.com/spf13/cobra"
)

var (
	rootCmd = &cobra.Command{
		Use:   "vcf-to-csv",
		Short: "Create CSV files from ConSTRain VCF output",
		Run: func(cmd *cobra.Command, args []string) {
			vcfPath, _ := cmd.Flags().GetString("vcf")
			csvPath, _ := cmd.Flags().GetString("output")
			vcfDir, _ := cmd.Flags().GetString("directory")
			csvDir, _ := cmd.Flags().GetString("outdir")

			nThreads, _ := cmd.Flags().GetInt64("threads")
			if nThreads > 0 {
				runtime.GOMAXPROCS(int(nThreads))
			} else if nThreads == 0 {
				log.Fatal("--threads must be greater than 0 (or -1 to use all available threads)")
			}

			file := vcfPath != "" && csvPath != ""
			dir := vcfDir != "" && csvDir != ""
			if file && !dir {
				vcfconv.RunFiles(vcfPath, csvPath)
			} else if !file && dir {
				vcfconv.RunDirs(vcfDir, csvDir)
			} else {
				panic("either `--vcf` and `--output` or `--directory` and `--outdir` need to be set")
			}
		},
	}
)

func Execute() error {
	return rootCmd.Execute()
}

func init() {
	rootCmd.PersistentFlags().StringP("vcf", "v", "", "path to VCF file")
	rootCmd.PersistentFlags().StringP("output", "o", "", "output handle to use for CSV file")

	rootCmd.PersistentFlags().StringP("directory", "d", "", "directory containing VCF files. A CSV will be created for each")
	rootCmd.PersistentFlags().StringP("outdir", "u", "", "directory where output cnvs will be generated")

	rootCmd.PersistentFlags().Int64P("threads", "t", -1, "maximum number of threads to use. Set to -1 to use all available threads (default: -1)")
}
