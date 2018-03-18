package phylotree

// valueSigma holds a value, it's standard deviation and
// a weight for weighted averages.
type valueSigma struct {
	value float64
	// sigma2 is the square of the standard deviation.
	sigma2 float64
	// weight is the weight for weighted averages.
	weight float64
}

// avgCalculator calculates a weighted average and it's standard deviation.
type avgCalculator struct {
	entries []valueSigma
	size    float64
}

// add adds a value and it's squared standard deviation to the calculator.
func (a *avgCalculator) add(value, sigma2 float64) {
	a.entries = append(a.entries, valueSigma{value: value, sigma2: sigma2})
	a.size++
}

// avg returns the weighted average and the square of the standard deviation
// of all entries in the calculator.
func (a *avgCalculator) avg() (average, sigma2 float64) {
	if a.size == 0 {
		return 0, 0
	}
	weightsTotal := 0.0
	for _, e := range a.entries {
		weightsTotal += e.value / e.sigma2
	}
	// Calculate the weights for each entry.
	for i, _ := range a.entries {
		a.entries[i].weight = a.entries[i].value / (weightsTotal * a.entries[i].sigma2)
	}
	// Weighted average.
	for _, e := range a.entries {
		average += e.value * e.weight
	}
	// Standard deviation to the sqare.
	for _, e := range a.entries {
		sigma2 += e.sigma2 * e.weight * e.weight
	}
	return average, sigma2
}
