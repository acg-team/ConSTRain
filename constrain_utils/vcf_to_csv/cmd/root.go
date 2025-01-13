package cmd

import (
	"log"
	"runtime"

	vcfconv "github.com/maverbiest/vcfconv/pkg"
	"github.com/spf13/cobra"
)

const VERSION = "0.4.3"

var (
	fileCmd = &cobra.Command{
		Use:   "file",
		Short: "Convert a single ConSTRain VCF output file to a CSV file",
		Run: func(cmd *cobra.Command, args []string) {
			vcfPath, _ := cmd.Flags().GetString("vcf")
			csvPath, _ := cmd.Flags().GetString("output")

			vcfconv.RunFile(vcfPath, csvPath)
		},
	}
	dirCmd = &cobra.Command{
		Use:   "dir",
		Short: "Convert all ConSTRain VCF output files in a directory to CSV files",
		Run: func(cmd *cobra.Command, args []string) {
			vcfDir, _ := cmd.Flags().GetString("directory")
			csvDir, _ := cmd.Flags().GetString("outdir")
			recursive, _ := cmd.Flags().GetBool("recursive")

			nThreads, _ := cmd.Flags().GetInt("threads")
			nCPU := runtime.NumCPU()
			nWorkers := 0
			if nThreads > 0 {
				nWorkers = min(nThreads, nCPU)
			} else if nThreads == -1 {
				nWorkers = nCPU
			} else if nThreads == 0 {
				log.Fatal("--threads must be greater than 0 (or -1 to use all available CPUs)")
			}

			vcfconv.RunDir(vcfDir, csvDir, nWorkers, recursive)
		},
	}
	rootCmd = &cobra.Command{
		Use:     "vcf-to-csv",
		Short:   "Create CSV files from ConSTRain VCF output",
		Version: VERSION,
	}
)

func Execute() error {
	rootCmd.AddCommand(fileCmd)
	rootCmd.AddCommand(dirCmd)

	return rootCmd.Execute()
}

func init() {
	fileCmd.Flags().StringP("vcf", "v", "", "path to VCF file")
	fileCmd.MarkFlagRequired("vcf")
	fileCmd.Flags().StringP("output", "o", "", "output handle to use for CSV file")
	fileCmd.MarkFlagRequired("output")

	dirCmd.Flags().StringP("directory", "d", "", "directory to search for VCF files")
	dirCmd.MarkFlagRequired("directory")

	dirCmd.Flags().BoolP("recursive", "r", false, "recursively search for VCF files in subdirectories of --directory as well (default false)")

	dirCmd.Flags().StringP("outdir", "o", "", "directory where output cnvs will be generated")
	dirCmd.MarkFlagRequired("outdir")

	rootCmd.PersistentFlags().IntP("threads", "t", -1, "maximum number of threads to use. Set to -1 to use all available threads")
}
