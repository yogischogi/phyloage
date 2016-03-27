// Package phyloage calculates time estimates based on a phylogenetic tree.
package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"strings"
	"time"

	"github.com/yogischogi/phyloage/phylotree"
	"github.com/yogischogi/phylofriend/genetic"
	"github.com/yogischogi/phylofriend/genfiles"
)

func main() {
	// Command line flags.
	var (
		treein     = flag.String("treein", "", "Input filename for phylogenetic tree (.txt).")
		treeout    = flag.String("treeout", "", "Output filename for phylogenetic tree in TXT format.")
		cal        = flag.Float64("cal", 1, "Calibration factor for TMRCA calculation.")
		personsin  = flag.String("personsin", "", "Input filename (.txt or .csv) or directory.")
		mrin       = flag.String("mrin", "", "Filename for the import of mutation rates.")
		gentime    = flag.Float64("gentime", 1, "Generation time in years.")
		inspect    = flag.String("inspect", "", "Comma separated list of SNP names to search for.")
		statistics = flag.Bool("statistics", false, "Prints marker statistics.")
		method     = flag.String("method", "parsimony", "Method to calculate modal haplotypes: phylofriend or parsimony.")
		stage      = flag.Int("stage", 3, "Processing stage for parsimony algorithm: 1, 2 or 3.")
		trace      = flag.String("trace", "", "Comma separated list of STR names to print out trace information.")
	)
	flag.Parse()

	var (
		persons       []*genetic.Person
		mutationRates genetic.YstrMarkers
		stat          *genetic.MarkerStatistics
		err           error
	)

	// Load phylogenetic tree from file.
	if *treein == "" {
		fmt.Printf("No filename for input tree specified.\n")
		os.Exit(1)
	}
	tree, err := phylotree.NewFromFile(*treein)
	if err != nil {
		fmt.Printf("Error reading tree from file, %v.\n", err)
		os.Exit(1)
	}

	// Read mutation rates from file.
	if *mrin != "" {
		mutationRates, err = genfiles.ReadMutationRates(*mrin)
		if err != nil {
			fmt.Printf("Error reading mutation rates %v.\n", err)
			os.Exit(1)
		}
	} else {
		// Use default values.
		mutationRates = genetic.DefaultMutationRates()
	}

	// Load genetic sample results.
	if *personsin != "" {
		filenames := strings.Split(*personsin, ",")
		for _, filename := range filenames {
			var pers []*genetic.Person
			fileInfo, err := os.Stat(filename)
			switch {
			case err != nil:
				fmt.Printf("Error, something is wrong with personsin, %v.\n", err)
				os.Exit(1)
			case fileInfo.IsDir():
				pers, err = genfiles.ReadPersonsFromDir(filename)
			case strings.HasSuffix(strings.ToLower(filename), ".csv"):
				pers, err = genfiles.ReadPersonsFromCSV(filename, 0)
			default:
				pers, err = genfiles.ReadPersonsFromTXT(filename)
			}
			if err != nil {
				fmt.Printf("Error loading persons data %v.\n", err)
				os.Exit(1)
			}
			persons = append(persons, pers...)
		}
		tree.InsertPersons(persons)

		// Calculate marker statistics.
		if *statistics == true || *method == "parsimony" {
			stat = genetic.NewStatistics(persons)
		}

		// Print marker statistics.
		if *statistics == true {
			fmt.Print(stat.String())

			// XXX Temporary code.
			// WriteToFile(stat)
		}

		// Calculate modal haplotypes.
		switch *method {
		case "phylofriend":
			tree.CalculateModalHaplotypes()
		case "parsimony":
			tree.CalculateModalHaplotypesParsimony(stat, *stage)
		default:
			fmt.Printf("Error, unknown method %q to calculate modal haplotypes.\n", *method)
			os.Exit(1)
		}

		tree.CalculateDistances(mutationRates, genetic.Distance)
	}

	// Calculate the age of this clade and all subclades.
	// If the STR-Count is provided in the original tree input
	// file the calculation can be performed even without sample
	// data.
	tree.CalculateAge(*gentime, *cal)

	// Save resulting tree to file or print it out.
	if *treeout != "" {
		date := time.Now().Format("2006 Jan 2")
		var buffer bytes.Buffer
		buffer.WriteString("// This tree was created by the phyloage program: https://github.com/yogischogi/phyloage\n")
		buffer.WriteString("// Command used:\n// ")
		for _, arg := range os.Args {
			buffer.WriteString(arg)
			buffer.WriteString(" ")
		}
		buffer.WriteString("\n")
		buffer.WriteString("// " + date + "\n\n")
		buffer.WriteString(tree.String())
		err := ioutil.WriteFile(*treeout, buffer.Bytes(), os.ModePerm)
		if err != nil {
			fmt.Printf("Error writing tree to file, %v.\n", err)
			os.Exit(1)
		}
	} else {
		fmt.Printf("%v\n", tree)
	}

	// Print tree with values of specified STRs.
	if *trace != "" {
		snps := strings.Split(*trace, ",")
		fmt.Printf("%s", tree.Trace(snps))
	}

	// Search for SNPs and print out information about the matching subclades.
	if *inspect != "" {
		searchTerms := strings.Split(*inspect, ",")
		fmt.Printf("%s", tree.Inspect(searchTerms))
	}
}

// XXX Temporary method to determine stable marker set.
func WriteToFile(statistics *genetic.MarkerStatistics) {
	filename := "mutrates.txt"
	minFreq := 1.0
	nValuesMin := 2
	nValuesMax := 3
	stats := statistics.Select(minFreq, nValuesMin, nValuesMax)
	err := ioutil.WriteFile(filename, []byte(stats.MutationRates()), os.ModePerm)
	if err != nil {
		fmt.Printf("Error writing mutation rates to file, %v.\n", err)
		os.Exit(1)
	}
}
