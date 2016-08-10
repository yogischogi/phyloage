package phylotree

import (
	"math"
)

// calculateModalHaplotypesMaxParsimony calculates modal haplotypes
// for this clade and all of it's sublcades that satisfy the maximum
// parsimony criterion. Because this is often not possible for all
// values, those values are set to Uncertain.
func (c *Clade) calculateModalHaplotypesMaxParsimony(isInfiniteAlleles bool) {
	// Calculate maximum parsimony for each marker.
	for i, _ := range c.Person.YstrMarkers {
		c.calculateMaxParsimony(i, isInfiniteAlleles)
	}
}

// calculateMaxParsimony calculates the most parsimonious value
// of a marker for this clade and all of it's subclades.
// If the method does not yield a clear result for a specific
// marker value, that value is set to Uncertain.
func (c *Clade) calculateMaxParsimony(marker int, isInfiniteAlleles bool) {
	// Calculate modal value using only downstream samples
	// and subclades.
	var values []float64
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			value := c.Samples[i].Person.YstrMarkers[marker]
			if value != 0 {
				values = append(values, value)
			}
		}
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].calculateMaxParsimony(marker, isInfiniteAlleles)
		value := c.Subclades[i].Person.YstrMarkers[marker]
		if value != 0 {
			values = append(values, value)
		}
	}
	modal := maxParsimony(values, isInfiniteAlleles)
	c.Person.YstrMarkers[marker] = modal

	// If we got a clear result, recalculate Uncertain values
	// in all subclades, but this time also include the parent
	// node in the calculation.
	if modal != Uncertain && modal != 0 {
		for i, _ := range c.Subclades {
			if c.Subclades[i].Person.YstrMarkers[marker] == Uncertain {
				c.Subclades[i].recalculateMaxParsimony(marker, isInfiniteAlleles, c)
			}
		}
	}
}

// recalculateMaxParsimony does nearly the same as calculateMaxParsimony,
// but this time the marker value of the parent clade is included.
// This can yield to a clear result, if the the child value can not
// be calculated from it's own child values, but the parent value is
// clear because of parallel subclades.
func (c *Clade) recalculateMaxParsimony(marker int, isInfiniteAlleles bool, parent *Clade) {
	var values []float64
	values = append(values, parent.Person.YstrMarkers[marker])
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			value := c.Samples[i].Person.YstrMarkers[marker]
			if value != 0 {
				values = append(values, value)
			}
		}
	}
	for i, _ := range c.Subclades {
		value := c.Subclades[i].Person.YstrMarkers[marker]
		if value != 0 {
			values = append(values, value)
		}
	}
	modal := maxParsimony(values, isInfiniteAlleles)
	c.Person.YstrMarkers[marker] = modal

	// If we got a clear result, recalculate Uncertain values
	// in all subclades.
	if modal != Uncertain && modal != 0 {
		for i, _ := range c.Subclades {
			if c.Subclades[i].Person.YstrMarkers[marker] == Uncertain {
				c.Subclades[i].recalculateMaxParsimony(marker, isInfiniteAlleles, c)
			}
		}
	}
}

// maxParsimony returns the value from values that satisfies
// the maximum parsimony criterion.
// To calculate the distance, the stepwise mutation model is used.
// If no unique result can be found, the result is Uncertain.
//
// I have compared multiple variations of this function using
// YFull tree 4.03 and data from the M343 xU106 xP312 project.
// At the time of the comparison (April 2016) the version of the
// method, using the stepwise mutation model, yielded the best
// results. I have checked TMRCA and formed estimates for a
// selection off different clades.
func maxParsimony(values []float64, isInfiniteAlleles bool) float64 {
	// stepwiseDist is the distance between two mutational values
	// using the stepwise mutation model.
	var stepwiseDist = func(a, b float64) float64 {
		return math.Abs(a - b)
	}

	// infiniteDist is the distance between two mutational values
	// using the infinite alleles model.
	var infiniteDist = func(a, b float64) float64 {
		if a == b {
			return 0
		} else {
			return 1
		}
	}

	// singleDist is the distance between two mutational values
	var singleDist func(a, b float64) float64
	if isInfiniteAlleles == true {
		singleDist = infiniteDist
	} else {
		singleDist = stepwiseDist
	}

	// totalDist is the number of mutations neccessary to reach
	// all values from x.
	var totalDist = func(x float64, values []float64) float64 {
		dist := 0.0
		for _, v := range values {
			dist += singleDist(x, v)
		}
		return dist
	}

	// Use only positive values for calculation.
	vals := make([]float64, 0, len(values))
	for _, v := range values {
		if v > 0 {
			vals = append(vals, v)
		}
	}
	// Test all values if one of them satisfies the minimum
	// distance criterion (maximum parsimony).
	var result float64 = 0
	minDist := math.Inf(1)
	for _, x := range vals {
		distance := totalDist(x, vals)
		if distance < minDist {
			minDist = distance
			result = x
		} else if distance == minDist && x != result {
			// Two different values satisfy the minimum condition.
			result = Uncertain
		}
	}
	return result
}
