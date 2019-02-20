// Package phylotree implements a phylogenetic tree based on SNP and Y-STR mutations.
package phylotree

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
	"unicode"

	"github.com/yogischogi/phylofriend/genetic"
)

// Uncertain is used for uncertain or unknown values.
const Uncertain = -1

// Element is a node or leaf of the tree.
type Element struct {
	SNPs []string
	// STRCount is the number of unique STR mutations for this element.
	STRCount float64
	// Person may be a real person from sample data
	// or a virtual ancestor (modal haplotype).
	// This may be nil.
	Person *genetic.Person
}

func newElement() Element {
	snps := make([]string, 0)
	return Element{SNPs: snps, STRCount: Uncertain}
}

func (e *Element) AddSNP(name string) {
	e.SNPs = append(e.SNPs, name)
}

// Contains checks if one of this element's SNPs
// equals searchTerm.
func (e *Element) Contains(searchTerm string) bool {
	result := false
	for _, snp := range e.SNPs {
		if strings.ToLower(snp) == strings.ToLower(searchTerm) {
			result = true
			break
		}
	}
	return result
}

// Details returns a detailed string representation of this element.
func (e *Element) Details() string {
	var buffer bytes.Buffer
	buffer.WriteString(e.String())
	buffer.WriteString("\r\n")
	if e.Person != nil {
		buffer.WriteString(e.Person.YstrMarkers.String())
	}
	return buffer.String()
}

func (e *Element) String() string {
	var buffer bytes.Buffer
	hasWritten := false
	n := len(e.SNPs)
	// Write SNPs
	if n > 0 {
		for i := 0; i < n-1; i++ {
			buffer.WriteString(e.SNPs[i])
			buffer.WriteString(", ")
		}
		buffer.WriteString(e.SNPs[n-1])
		hasWritten = true
	}
	// Write STR-Count.
	if e.STRCount >= 0 {
		if hasWritten {
			buffer.WriteString(fmt.Sprintf(", STR-Count: %.0f", e.STRCount))
		} else {
			buffer.WriteString(fmt.Sprintf("STR-Count: %.0f", e.STRCount))
		}
	}
	return buffer.String()
}

// strDetails returns a textual representation (names and values)
// of the Y-STR markers specified by indices.
func (e *Element) strDetails(indices []int) string {
	var buffer bytes.Buffer
	for _, i := range indices {
		if e.Person != nil {
			name := genetic.YstrMarkerTable[i].InternalName
			value := e.Person.YstrMarkers[i]
			buffer.WriteString(fmt.Sprintf(" %s: %g,", name, value))
		}
	}
	return buffer.String()
}

// Sample represents the genetic sample of an individual.
// This is a leaf of the phylogenetic tree.
type Sample struct {
	Element
	// ID od this sample, usually the kit number.
	ID string
}

func newSample() Sample {
	return Sample{Element: newElement()}
}

// newSample creates a new Sample from a textual representation.
// Format: id:SampleID, SNP1, SNP2, STR-Count: 11
// Only the "id:" field is mandatory.
func newSampleFromText(text string) (Sample, error) {
	result := newSample()
	tokens := strings.Split(text, ",")
	for _, token := range tokens {
		token = strings.TrimSpace(token)
		switch {
		case strings.HasPrefix(token, "id:"):
			result.ID = strings.TrimSpace(token[3:])
		case strings.HasPrefix(token, "STR-Count:"):
			strCount := strings.TrimSpace(token[10:])
			count, err := strconv.ParseFloat(strCount, 64)
			if err != nil {
				msg := fmt.Sprintf("could not convert STR-Count to float: %s", strCount)
				return result, errors.New(msg)
			}
			result.STRCount = count
		default:
			result.AddSNP(token)
		}
	}
	return result, nil
}

// Contains checks if one of this sample's SNPs or the ID
// equals searchTerm.
func (s *Sample) Contains(searchTerm string) bool {
	return s.Element.Contains(searchTerm) || strings.ToLower(searchTerm) == strings.ToLower(s.ID)
}

func (s *Sample) String() string {
	if s.Element.String() != "" {
		return fmt.Sprintf("id:%s, %s", s.ID, s.Element.String())
	} else {
		return fmt.Sprintf("id:%s", s.ID)
	}
}

func (s *Sample) Details() string {
	return fmt.Sprintf("id:%s, %s", s.ID, s.Element.Details())
}

// Clade models a haplogroup or subclade.
// A clade represents a node of the phylogenetic tree.
type Clade struct {
	Element
	Samples   []Sample
	Subclades []Clade
	// AgeSTR shows when this Clade has formed ybp
	// according to a calculation using Y-STR mutations.
	AgeSTR float64
	// STRCountDownstream is the average number of downstream STR mutations.
	STRCountDownstream float64
	// Sigma2 is the squared standard deviation of STRCountDownstream.
	Sigma2 float64
	// TMRCA_STR is the time to the most recent common Ancestor
	// for all downstream samples.
	TMRCA_STR float64
	// TMRCAlower and TMRCAupper are the lower and upper bounds
	// of the 95% confidence interval.
	TMRCAlower float64
	TMRCAupper float64
}

// newClade creates a new Clade from a textual representation.
// Format: SNP1, SNP2, STR-Count: 11
// "STR-Count:" is optional.
func newClade(text string) (Clade, error) {
	result := Clade{
		Element:            newElement(),
		AgeSTR:             Uncertain,
		STRCountDownstream: Uncertain,
		TMRCA_STR:          Uncertain}
	tokens := strings.Split(text, ",")
	for _, token := range tokens {
		token = strings.TrimSpace(token)
		switch {
		case strings.HasPrefix(token, "STR-Count:"):
			strCount := strings.TrimSpace(token[10:])
			count, err := strconv.ParseFloat(strCount, 64)
			if err != nil {
				msg := fmt.Sprintf("could not convert STR-Count to float: %s", strCount)
				return result, errors.New(msg)
			}
			result.STRCount = count
		case strings.HasPrefix(token, "TMRCA:"):
			// Ignore because this TMRCA has to be newly calculated.
		default:
			result.AddSNP(token)
		}
	}
	return result, nil
}

// NewFromFile parses a text file to create a tree.
// The return value Clade is the root node of the tree.
func NewFromFile(filename string) (*Clade, error) {
	lines := make([]lineInfo, 0)

	// Open file
	infile, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer infile.Close()

	// Read lines
	lineNo := 0
	scanner := bufio.NewScanner(infile)
	for scanner.Scan() {
		lineNo++
		text := stripComments(scanner.Text())
		if text != "" {
			indent := countSpaces(text)
			lines = append(lines, lineInfo{lineNo: lineNo, indent: indent, text: text})
		}
	}
	if scanner.Err() != nil {
		return nil, scanner.Err()
	}
	if len(lines) == 0 {
		return nil, errors.New("empty file, nothing to do")
	}

	// Build tree by parsing lines.
	root, err := newClade(lines[0].text)
	if err != nil {
		return nil, errors.New(fmt.Sprintf("invalid root element, %s", err))
	}
	err = parseTree(&root, lines[0].indent, lines[1:])
	if err != nil {
		return nil, errors.New(fmt.Sprintf("parsing tree, %s", err))
	}
	return &root, nil
}

func (c *Clade) AddSample(sample Sample) {
	if c.Samples == nil {
		c.Samples = make([]Sample, 0)
	}
	c.Samples = append(c.Samples, sample)
}

func (c *Clade) AddSubclade(clade Clade) {
	if c.Subclades == nil {
		c.Subclades = make([]Clade, 0)
	}
	c.Subclades = append(c.Subclades, clade)
}

// Persons returns a list of all persons who belong to this clade.
// This includes the calculated modal haplotypes.
func (c *Clade) Persons() []*genetic.Person {
	persons := make([]*genetic.Person, 0, 50)
	if c.Person != nil {
		persons = append(persons, c.Person)
	}
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			persons = append(persons, c.Samples[i].Person)
		}
	}
	for i, _ := range c.Subclades {
		persons = append(persons, c.Subclades[i].Persons()...)
	}
	return persons
}

// InsertPersons traverses the tree and adds the appropriate
// person to a leaf if the ID of the sample and the person's ID
// are identical.
func (c *Clade) InsertPersons(persons []*genetic.Person) {
	// Create hash map of persons' IDs.
	personsMap := make(map[string]*genetic.Person)
	for i, _ := range persons {
		personsMap[persons[i].ID] = persons[i]
	}
	// Search samples matching person IDs.
	for i, _ := range c.Samples {
		if person, exists := personsMap[c.Samples[i].ID]; exists {
			c.Samples[i].Person = person
		}
	}
	// Search subclades for matching person IDs.
	for i, _ := range c.Subclades {
		c.Subclades[i].InsertPersons(persons)
	}
}

// CalculateModalHaplotypes calculates the modal haplotype for
// this clade and all subclades.
func (c *Clade) CalculateModalHaplotypes() {
	// Create a list of haplotypes from samples and subclades.
	persons := make([]*genetic.Person, 0)
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			persons = append(persons, c.Samples[i].Person)
		}
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].CalculateModalHaplotypes()
		if c.Subclades[i].Person != nil {
			persons = append(persons, c.Subclades[i].Person)
		}
	}
	// Calculate modal haplotype for the list.
	modal := genetic.ModalHaplotype(persons)
	modal.ID = c.SNPs[0]
	modal.Name = c.SNPs[0]
	modal.Label = c.SNPs[0]
	c.Person = modal
}

// CalculateDistances calculated the genetic distances between
// the modal haplotype of this clade and it's downstream members.
func (c *Clade) CalculateDistances(mutationRates genetic.YstrMarkers, distance genetic.DistanceFunc) {
	if c.Person == nil {
		return
	}
	// Calculate genetic distances between modal haplotype
	// and samples or subclades.
	c.STRCount = 0
	for i, _ := range c.Samples {
		if c.Samples[i].Person != nil {
			ystr1 := c.Samples[i].Person.YstrMarkers
			ystr2 := c.Person.YstrMarkers
			c.Samples[i].STRCount = distance(ystr1, ystr2, mutationRates)
		}
	}
	for i, _ := range c.Subclades {
		c.Subclades[i].CalculateDistances(mutationRates, distance)
		if c.Subclades[i].Person != nil {
			ystr1 := c.Subclades[i].Person.YstrMarkers
			ystr2 := c.Person.YstrMarkers
			c.Subclades[i].STRCount = distance(ystr1, ystr2, mutationRates)
		}
	}
}

// CalculateAge calculates the age and TMRCA for this Clade.
// It fills the following variables insise Clade:
// TMRCA_STR, AgeSTR, STRCountDownstream.
// gentime is the generation time in years.
// calibration is a calibration factor that is multiplied
// to the result.
// offset is added to all calculated ages to account for the ages
// of living persons. YFull currently uses an offset of 60 years.
func (c *Clade) CalculateAge(gentime, calibration, offset float64) {
	var avgCalc avgCalculator
	// Count STR mutations for samples.
	// average value
	avgSamples := 0.0
	// sigma squared
	sigma2Samples := 0.0
	nSamples := float64(len(c.Samples))
	if nSamples > 0 {
		for i, _ := range c.Samples {
			if c.Samples[i].STRCount > 0 {
				avgSamples += c.Samples[i].STRCount
			}
		}
		avgSamples /= nSamples
		sigma2Samples = avgSamples / nSamples
		if sigma2Samples > 0 {
			avgCalc.add(avgSamples, sigma2Samples)
		}
	}
	// Count STR mutations for subclades.
	for i, _ := range c.Subclades {
		c.Subclades[i].CalculateAge(gentime, calibration, offset)
		subcladeSTRs := c.Subclades[i].STRCount + c.Subclades[i].STRCountDownstream
		subcladeSigma2 := c.Subclades[i].STRCount + c.Subclades[i].Sigma2
		if subcladeSigma2 > 0 {
			avgCalc.add(subcladeSTRs, subcladeSigma2)
		}
	}
	// Calculate average number of mutations.
	if avgCalc.size > 0 {
		c.STRCountDownstream, c.Sigma2 = avgCalc.avg()
		c.TMRCA_STR = c.STRCountDownstream*gentime*calibration + offset
		c.AgeSTR = (c.STRCount+c.STRCountDownstream)*gentime*calibration + offset
		lower, upper := avgCalc.confidenceIntervals(c.STRCountDownstream, c.Sigma2)
		c.TMRCAlower = lower*gentime*calibration + offset
		c.TMRCAupper = upper*gentime*calibration + offset
	}
}

func (c *Clade) String() string {
	var buffer bytes.Buffer
	c.prettyPrint(&buffer, 0)
	return buffer.String()
}

// prettyPrint prints a formatted version of the clade c
// into buffer. indent is the indentation for the root node.
func (c *Clade) prettyPrint(buffer *bytes.Buffer, indent int) {
	// Write this Element.
	for i := 0; i < indent; i++ {
		buffer.WriteString("\t")
	}
	buffer.WriteString(c.Element.String())

	// Write time estimates.
	if c.STRCountDownstream >= 0 {
		buffer.WriteString(
			fmt.Sprintf(", STRs Downstream: %.0f, formed: %.0f, TMRCA: %.0f, CI:[%.0f, %.0f]",
				c.STRCountDownstream, c.AgeSTR, c.TMRCA_STR, c.TMRCAlower, c.TMRCAupper))
	}
	buffer.WriteString("\r\n")

	// Write Samples.
	for _, sample := range c.Samples {
		for i := 0; i < indent+1; i++ {
			buffer.WriteString("\t")
		}
		buffer.WriteString(sample.String())
		buffer.WriteString("\r\n")
	}
	// Write Subclades.
	for _, clade := range c.Subclades {
		clade.prettyPrint(buffer, indent+1)
	}
}

// Inspect looks at this clade and all subclades.
// If any of the tree nodes' SNPs match one of the search
// terms, a string representation of the element is added
// to the result.
func (c *Clade) Inspect(searchTerms []string) string {
	// Create hash map containing results.
	results := make(map[string]string)
	for _, term := range searchTerms {
		results[term] = ""
	}

	results = c.searchFor(results)

	// Return representation of the findings.
	var buffer bytes.Buffer
	for _, term := range searchTerms {
		buffer.WriteString(results[term])
	}
	return buffer.String()
}

// searchFor searches for SNPs in this clade and it's subclades.
// The SNPs are defined as keys in the results map.
// The values of the results map are string representations of
// the matching clades.
func (c *Clade) searchFor(results map[string]string) map[string]string {
	// Search for searchTerms.
	for key, _ := range results {
		if c.Contains(key) {
			results[key] = c.Details()
		}
	}
	for i, _ := range c.Samples {
		for key, _ := range results {
			if c.Samples[i].Contains(key) {
				results[key] = c.Samples[i].Details()
			}
		}
	}
	for i, _ := range c.Subclades {
		for key, _ := range results {
			if c.Subclades[i].Contains(key) {
				results[key] = c.Subclades[i].Element.Details()
			}
		}
		results = c.Subclades[i].searchFor(results)
	}
	return results
}

// Trace returns a nicely formatted tree containing information
// (names and values) about the Y-STR markers specified by STRs.
func (c *Clade) Trace(STRs []string) string {
	// Look for STR indices.
	var indices []int
	for _, str := range STRs {
		str := strings.ToLower(str)
		for _, marker := range genetic.YstrMarkerTable {
			if str == strings.ToLower(marker.InternalName) ||
				str == strings.ToLower(marker.FTDNAName) ||
				str == strings.ToLower(marker.YFullName) {
				indices = append(indices, marker.Index)
				break
			}
		}
	}

	// Build tree with STR values.
	var buffer bytes.Buffer
	c.tracePrint(&buffer, 0, indices)
	return buffer.String()
}

// tracePrint creates the formatted tree for Trace.
func (c *Clade) tracePrint(buffer *bytes.Buffer, indent int, STRindices []int) {
	// Write this Element.
	for i := 0; i < indent; i++ {
		buffer.WriteString("\t")
	}
	buffer.WriteString(c.Element.String())
	buffer.WriteString(",")
	buffer.WriteString(c.Element.strDetails(STRindices))
	buffer.WriteString("\r\n")

	// Write Samples.
	for _, sample := range c.Samples {
		for i := 0; i < indent+1; i++ {
			buffer.WriteString("\t")
		}
		buffer.WriteString(sample.String())
		buffer.WriteString(",")
		buffer.WriteString(sample.Element.strDetails(STRindices))
		buffer.WriteString("\r\n")
	}
	// Write Subclades.
	for _, clade := range c.Subclades {
		clade.tracePrint(buffer, indent+1, STRindices)
	}
}

// Subclade returns the subclade that contains searchTerm.
func (c *Clade) Subclade(cladeName string) *Clade {
	var result *Clade
	if c.contains(cladeName) {
		result = c
	} else {
		for i, _ := range c.Subclades {
			if c.Subclades[i].contains(cladeName) {
				result = &c.Subclades[i]
				break
			} else {
				result = c.Subclades[i].Subclade(cladeName)
				if result != nil {
					break
				}
			}
		}
	}
	return result
}

// contains checks if the text representation of this clade
// contains cladeName.
// The method only returns true if the found cladeName is
// terminated by a non digit. This should make sure that the
// found cladeName is not part of a longer SNP name.
func (c *Clade) contains(cladeName string) bool {
	for _, snp := range c.SNPs {
		if strings.ToLower(snp) == strings.ToLower(cladeName) {
			return true
		}
	}
	return false
}

// lineInfo is a helper struct for parsing a tree in text format.
type lineInfo struct {
	lineNo int
	indent int
	text   string
}

// parseTree parses a tree in text format with white space indentations.
// The function works recursively and adds all new subclades and samples
// to the parent clade. indent is the indentation of the parent clade
// in the text file.
func parseTree(parent *Clade, indent int, lines []lineInfo) error {
	childIndent := -1
	for i, _ := range lines {
		switch {
		case lines[i].indent <= indent:
			// Return if indentation shows beginning of next block.
			return nil
		case childIndent == -1 && lines[i].indent > indent:
			// Determine the indentation of the child block.
			childIndent = lines[i].indent
			fallthrough
		case lines[i].indent == childIndent:
			// Parse child elements.
			if strings.Contains(lines[i].text, "id:") {
				// Child is Sample element.
				sample, err := newSampleFromText(lines[i].text)
				if err != nil {
					msg := fmt.Sprintf("line: %d, %s", lines[i].lineNo, err)
					return errors.New(msg)
				}
				parent.AddSample(sample)
			} else {
				// Child is Clade element.
				clade, err := newClade(lines[i].text)
				if err != nil {
					msg := fmt.Sprintf("line: %d, %s", lines[i].lineNo, err)
					return errors.New(msg)
				}
				parseTree(&clade, lines[i].indent, lines[i+1:])
				parent.AddSubclade(clade)
			}
		}
	}
	return nil
}

// stripComments removes comments from a line of text.
// A comment starts with // like in many programming languages.
// If a line contains only a comment or only whitespace characters
// an empty string is returned.
func stripComments(line string) string {
	result := line
	// Remove comment from line.
	idxComment := strings.Index(line, "//")
	if idxComment >= 0 {
		result = line[0:idxComment]
	}
	// Check if line contains only whitespaces.
	isEmpty := true
	for _, c := range result {
		if unicode.IsSpace(c) == false {
			isEmpty = false
			break
		}
	}
	if isEmpty {
		return ""
	} else {
		return result
	}
}

// countSpaces counts the white spaces at the beginning
// of a line.
func countSpaces(line string) int {
	result := 0
	for _, c := range line {
		if unicode.IsSpace(c) {
			result++
		} else {
			break
		}
	}
	return result
}
