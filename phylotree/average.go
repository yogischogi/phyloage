package phylotree

import (
	"math"
)

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

// confidenceIntervalsNormal returns the lower and the upper
// bound for a 95% confidence interval and a number of mutations.
// average is the average number of mutations and sigma2 it's variance.
func (a *avgCalculator) confidenceIntervals(average, sigma2 float64) (lower, upper float64) {
	s := math.Sqrt(sigma2)
	// Noise to signal
	ns := average / s
	// Number of mutations for normalized function.
	m := ns * ns
	low, high := confidenceIntervalsNormal(m)
	c := average / m
	return c * low, c * high
}

// confidenceIntervalsNormal returns the lower and the upper
// bound for a 95% confidence interval and a number of mutations m.
// The functions assumes a normalized Poisson function.
// For m values < 15 the Poisson distribution is used.
// For higher values the result uses a Gaussian distribution.
func confidenceIntervalsNormal(m float64) (lower, upper float64) {
	// Lower and upper bounds for Poisson distribution.
	poisLows := [...]float64{0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8}
	poisHighs := [...]float64{5, 7, 8, 10, 11, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24}
	if m < 15 {
		// Poisson range
		i := int(m)
		return poisLows[i], poisHighs[i]
	} else {
		// Gauss case.
		lower = m - 2.0*math.Sqrt(m+1.0) + 2.0
		upper = m + 2.0*math.Sqrt(m+1.0) + 2.0
		return
	}
}
