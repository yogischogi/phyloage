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
		treein    = flag.String("treein", "", "Input filename for phylogenetic tree (.txt).")
		treeout   = flag.String("treeout", "", "Output filename for phylogenetic tree in TXT format.")
		cal       = flag.Float64("cal", 1, "Calibration factor for TMRCA calculation.")
		personsin = flag.String("personsin", "", "Input filename (.txt or .csv) or directory.")
		mrin      = flag.String("mrin", "", "Filename for the import of mutation rates.")
		gentime   = flag.Float64("gentime", 1, "Generation time in years.")
		idcol     = flag.Int("idcol", 1, "Column number for IDs in CSV input file.")
	)
	flag.Parse()

	var (
		persons       []*genetic.Person
		mutationRates genetic.YstrMarkers
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
		fileInfo, err := os.Stat(*personsin)
		switch {
		case err != nil:
			fmt.Printf("Error, something is wrong with personsin, %v.\n", err)
			os.Exit(1)
		case fileInfo.IsDir():
			persons, err = genfiles.ReadPersonsFromDir(*personsin)
		case strings.HasSuffix(strings.ToLower(*personsin), ".csv"):
			persons, err = genfiles.ReadPersonsFromCSV(*personsin, *idcol-1)
		default:
			persons, err = genfiles.ReadPersonsFromTXT(*personsin)
		}
		if err != nil {
			fmt.Printf("Error loading persons data %v.\n", err)
			os.Exit(1)
		}
		tree.InsertPersons(persons)
		tree.CalculateModalHaplotypes()
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
}
