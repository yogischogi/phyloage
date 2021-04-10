# Phyloage

Phyloage is an experimental program to calculate time
estimates (TMRCA, Time To Most Recent Common Ancestor)
from SNP based phylogenetic trees.

It uses STR mutation counts for each node of the phylogenetic
tree and thus tries to circumvent the traditional shortcomings
of STR counting, which is often heavily influenced by younger
family branches.

The program can also be used to generate time estimates from
the number of SNP mutations on each branch of a phylogenetic
tree. The algorithm

* uses weighted averages.
* performs a top down recalculation of the input tree so that
  subclades can never be older than it's parent.
* calculates 95% confidence intervals analytically.


## Documentation

* [User Guide](https://github.com/yogischogi/phyloage/blob/master/doc/phyloage.pdf?raw=true)
* [Source Code](http://godoc.org/github.com/yogischogi/phyloage)

