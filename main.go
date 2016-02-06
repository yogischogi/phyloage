// Package phyloage calculates time estimates based on a phylogenetic tree.
package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"time"

	"github.com/yogischogi/phyloage/phylotree"
)

func main() {
	var (
		treein  = flag.String("treein", "", "Input filename for phylogenetic tree (.txt).")
		treeout = flag.String("treeout", "", "Output filename for phylogenetic tree in TXT format.")
		cal     = flag.Float64("cal", 1, "Calibration factor for TMRCA calculation.")
	)
	flag.Parse()

	// Load tree from file and calculate ages.
	if *treein == "" {
		fmt.Printf("No filename for input tree specified.\n")
		os.Exit(1)
	}
	tree, err := phylotree.NewFromFile(*treein)
	if err != nil {
		fmt.Printf("Error reading tree from file, %v.\n", err)
		os.Exit(1)
	}
	tree.CalculateAge(*cal)

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
