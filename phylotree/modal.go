package phylotree

import (
	"math"

	"github.com/yogischogi/phylofriend/genetic"
)

// resultOption specifies how float values should be mapped
// to marker values.
type mappingOption int

const (
	// mapAll maps all values to real world marker values.
	mapAll = iota
	// mapCertain maps only values to real word marker values,
	// if the original value has a clearly defined nearest neighbor
	// among the real world marker values.
	mapCertain
	// markUncertain maps values with a clearly defined nearest neighbor
	// to real world marker values and marks unclear values as -1.
	markUncertain
)

// CalculateModalHaplotypesParsimony calculates all modal haplotypes
// for this clade, using a method of maximum parsimony.
//
// Processing stages:
//
//  0: Return without doing anything.
//
//  1: Calculate haplotypes that satisfy the maximum
//		parsimony criterion.
//
//  2: Calculate average haplotypes using real numbers.
//
//  3: Mark results that do not have a nearest neighbor among
//     real mutation values as Uncertain.
//     This stage is only for visualization and debugging.
//
//  4: Replace uncertain values in the top node by using the
//     nearest and smallest real mutation neighbor.
//     Recalculate the tree top down to find values for previously
//     uncertain values.
func (c *Clade) CalculateModalHaplotypesParsimony(statistics *genetic.MarkerStatistics, processingStage int) {
	if processingStage < 1 {
		return
	}
	if processingStage >= 1 {
		// Make sure all node Persons are != nil.
		c.populateWithDummies()

		// Calculate haplotypes that satisfy the maximum
		// parsimony criterion.
		c.calculateModalHaplotypesMaxParsimony()
	}
	if processingStage >= 2 {
		// Calculate average haplotypes using real numbers.
		c.calculateHaplotypes(averageHaplotype)
	}
	if processingStage == 3 {
		// Mark results that do not have a nearest neighbor among
		// real mutation values as Uncertain.
		// This stage is only for visualization and debugging.
		c.constrainHaplotypes(statistics, markUncertain)
	}
	if processingStage >= 4 {
		// Map all markers with certain nearest neighbors to real world marker values.
		c.constrainHaplotypes(statistics, mapCertain)

		// Force a haplotype without uncertain values for the top node.
		constrainHaplotype(c.Person, statistics, mapAll)

		// Recalculate values for uncertain values
		// using child and parent haplotypes.
		for i, _ := range c.Subclades {
			c.Subclades[i].recalculateModalHaplotypes(c, statistics)
		}
	}
}

// populateWithDummies adds a person with 0 values to this
// clades and all of it's subclades. Persons are untouched.
func (c *Clade) populateWithDummies() {
	c.Person = &genetic.Person{
		ID:    c.SNPs[0],
		Name:  c.SNPs[0],
		Label: c.SNPs[0]}
	for i, _ := range c.Subclades {
		c.Subclades[i].populateWithDummies()
	}
}

// calculateHaplotypes calculates the haplotype for this clade
// and all subclades, using the function haplotype to perform
// the calculation.
// Only Uncertain values are replaced in the nodes of the Clade.
func (c *Clade) calculateHaplotypes(haplotype func(persons []*genetic.Person) *genetic.Person) {
	// Create a list of haplotypes from samples and subclades.
	persons := make([]*genetic.Person, 0)
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			persons = append(persons, c.Samples[i].Person)
		}
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].calculateHaplotypes(haplotype)
		persons = append(persons, c.Subclades[i].Person)
	}
	// Calculate result and replace Uncertain values in this Clade.
	modal := haplotype(persons)
	replaceUncertains(c.Person, modal)
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
func (c *Clade) constrainHaplotypes(statistics *genetic.MarkerStatistics, mapping mappingOption) {
	// Create a list of modal haplotypes.
	persons := make([]*genetic.Person, 0)
	if c.Person != nil {
		persons = append(persons, c.Person)
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].constrainHaplotypes(statistics, mapping)
		if c.Subclades[i].Person != nil {
			persons = append(persons, c.Subclades[i].Person)
		}
	}
	// Contrain haplotypes to real world values.
	for i, _ := range persons {
		constrainHaplotype(persons[i], statistics, mapping)
	}
}

// constrainHaplotype mappes the marker values of person to the
// closest real world marker values from the marker statistics.
func constrainHaplotype(person *genetic.Person, statistics *genetic.MarkerStatistics, mapping mappingOption) {
	for i, _ := range person.YstrMarkers {
		closest, isUnique := closestKey(person.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
		switch {
		case isUnique || mapping == mapAll:
			person.YstrMarkers[i] = closest
		case !isUnique && mapping == markUncertain:
			person.YstrMarkers[i] = Uncertain
		default:
			// Do nothing.
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
	replaceUncertainsWithMapping(c.Person, recalc, statistics)

	for i, _ := range c.Subclades {
		c.Subclades[i].recalculateModalHaplotypes(c, statistics)
	}
}

// replaceUncertains replaces all uncertain marker values in target
// with values from source.
func replaceUncertains(target, source *genetic.Person) {
	for i, _ := range target.YstrMarkers {
		if target.YstrMarkers[i] == Uncertain {
			target.YstrMarkers[i] = source.YstrMarkers[i]
		}
	}
}

// replaceUncertainsWithMapping replaces all uncertain marker values in target
// with values from source. The result is mapped to the closest real
// world marker values using the marker statistics.
func replaceUncertainsWithMapping(target, source *genetic.Person, statistics *genetic.MarkerStatistics) {
	for i, _ := range target.YstrMarkers {
		_, isUnique := closestKey(target.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
		if isUnique == false {
			newValue, _ := closestKey(source.YstrMarkers[i], statistics.Markers[i].ValuesOccurrences)
			target.YstrMarkers[i] = newValue
		}
	}
}
