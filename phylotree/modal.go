package phylotree

import (
	"math"

	"github.com/yogischogi/phylofriend/genetic"
)

// CalculateModalHaplotypesParsimony calculates all modal haplotypes
// for this clade, using a method of maximum parsimony.
//
// Processing stages:
//
//  0: Return without doing anything.
//
//  1: Calculate average haplotypes using real numbers.
//
//  2: Mark results that do not have a nearest neighbor among
//     real mutation values as Uncertain.
//     This stage is only for visualization and debugging.
//
//  3: Replace uncertain values in the top node by using the
//     nearest and smallest real mutation neighbor.
//     Recalculate the tree top down to find values for previously
//     uncertain values.
func (c *Clade) CalculateModalHaplotypesParsimony(statistics *genetic.MarkerStatistics, processingStage int) {
	if processingStage < 1 {
		return
	}
	if processingStage >= 1 {
		// Calculate average haplotypes using real numbers.
		c.calculateAverageHaplotypes()
	}
	if processingStage == 2 {
		// Mark results that do not have a nearest neighbor among
		// real mutation values as Uncertain.
		// This stage is only for visualization and debugging.
		c.constrainHaplotypes(statistics, false)
	}
	if processingStage >= 3 {
		// Force a haplotype without uncertain values for the top node.
		constrainHaplotype(c.Person, statistics, true)

		// Recalculate values for uncertain values
		// using child and parent haplotypes.
		for i, _ := range c.Subclades {
			c.Subclades[i].recalculateModalHaplotypes(c, statistics)
		}
	}
}

// calculateAverageHaplotypes calculates the average haplotype for
// this clade and all subclades.
func (c *Clade) calculateAverageHaplotypes() {
	// Create a list of haplotypes from samples and subclades.
	persons := make([]*genetic.Person, 0)
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			persons = append(persons, c.Samples[i].Person)
		}
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].calculateAverageHaplotypes()
		if c.Subclades[i].Person != nil {
			persons = append(persons, c.Subclades[i].Person)
		}
	}
	// Calculate modal haplotype for the list.
	modal := averageHaplotype(persons)
	modal.ID = c.SNPs[0]
	modal.Name = c.SNPs[0]
	modal.Label = c.SNPs[0]
	c.Person = modal
}

// averageHaplotype calculates the average haplotype for a group of persons.
// This method works with real numbers and uses the average
// value as the modal value, while values of real mutatations
// can only be whole numbers.
//
// Criteria for the modal value of a mutation:
//
// 	1 value: return value, because the modal of a
//		single haplotype is the haplotype itself.
//
//  >2 values: return the average of all values > 0.
func averageHaplotype(persons []*genetic.Person) *genetic.Person {
	modal := new(genetic.Person)
	switch len(persons) {
	case 0:
		// Return set of empty values.
		return modal
	case 1:
		// Return the person itself.
		modal = persons[0]
		return modal
	default:
		// Calculate modal value for each marker.
		for marker := 0; marker < len(persons[0].YstrMarkers); marker++ {
			count := 0.0
			sum := 0.0
			for _, person := range persons {
				value := person.YstrMarkers[marker]
				if value > 0 {
					sum += value
					count++
				}
			}
			if count > 0 {
				modal.YstrMarkers[marker] = sum / count
			}
		}
		return modal
	}
}

// constrainHaplotypes maps calculated marker values to real
// world marker values using the marker statistics.
// If forceResult == true, all values are mapped to some real marker.
// If forceReslut == false, markers that do not have a real world
// closest neighbor are left Uncertain.
func (c *Clade) constrainHaplotypes(statistics *genetic.MarkerStatistics, forceResult bool) {
	// Create a list of modal haplotypes.
	persons := make([]*genetic.Person, 0)
	if c.Person != nil {
		persons = append(persons, c.Person)
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].constrainHaplotypes(statistics, forceResult)
		if c.Subclades[i].Person != nil {
			persons = append(persons, c.Subclades[i].Person)
		}
	}
	// Contrain haplotypes to real world values.
	for i, _ := range persons {
		constrainHaplotype(persons[i], statistics, forceResult)
	}
}

// constrainHaplotype mappes the marker values of person to the
// closest real world marker values from the marker statistics.
func constrainHaplotype(person *genetic.Person, statistics *genetic.MarkerStatistics, forceResult bool) {
	for i, _ := range person.YstrMarkers {
		closest, isUnique := closestKey(person.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
		switch {
		case isUnique || forceResult:
			person.YstrMarkers[i] = closest
		default:
			person.YstrMarkers[i] = Uncertain
		}
	}
}

// closest Key returns the key of the mutations map that is closest
// to the target value. If two keys are equally close, isUnique = false.
// The keys of the mutations map must hold the mutational values.
// if target <= 0, the target value itself is returned, because
// a valid mutational value must always be positive.
func closestKey(target float64, mutations map[float64]int) (closest float64, isUnique bool) {
	if target <= 0 {
		return target, true
	}
	var smallestMutation float64
	var greatestMutation float64
	lowDist := math.Inf(1)
	highDist := math.Inf(1)
	for mutation, _ := range mutations {
		dist := mutation - target
		switch {
		case dist < 0:
			dist = math.Abs(dist)
			if dist < lowDist {
				lowDist = dist
				smallestMutation = mutation
			}
		case dist > 0:
			dist = math.Abs(dist)
			if dist < highDist {
				highDist = dist
				greatestMutation = mutation
			}
		case dist == 0:
			lowDist = 0
			highDist = 0
			smallestMutation = mutation
			greatestMutation = mutation
		}
	}
	switch {
	case smallestMutation == greatestMutation:
		return smallestMutation, true
	case lowDist == highDist:
		// Prefer to return the smallest mutation that is
		// close to the target.
		return smallestMutation, false
	default:
		if lowDist < highDist {
			return smallestMutation, true
		} else {
			return greatestMutation, true
		}
	}
}

// recalculateModalHaplotypes recalculates modal haplotypes that
// contain uncertain values. It takes the average of the child
// haplotypes and the parent haplotype. The result is mapped to
// the closest set of real marker values.
func (c *Clade) recalculateModalHaplotypes(parent *Clade, statistics *genetic.MarkerStatistics) {
	// Create a list of haplotypes for calculation.
	persons := make([]*genetic.Person, 0)
	if parent != nil && parent.Person != nil {
		persons = append(persons, parent.Person)
	}
	for i, _ := range c.Subclades {
		if c.Subclades[i].Person != nil {
			persons = append(persons, c.Subclades[i].Person)
		}
	}
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			persons = append(persons, c.Samples[i].Person)
		}
	}
	recalc := averageHaplotype(persons)
	replaceUncertainValues(c.Person, recalc, statistics)

	for i, _ := range c.Subclades {
		c.Subclades[i].recalculateModalHaplotypes(c, statistics)
	}
}

// replaceUncertainValues replaces all uncertain marker values in target
// with values from source. The result is mapped to the closest real
// world marker values using the marker statistics.
func replaceUncertainValues(target, source *genetic.Person, statistics *genetic.MarkerStatistics) {
	for i, _ := range target.YstrMarkers {
		_, isUnique := closestKey(target.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
		if isUnique == false {
			newValue, _ := closestKey(source.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
			target.YstrMarkers[i] = newValue
		}
	}
}
