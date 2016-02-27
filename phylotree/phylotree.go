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

// Element is a node or a leaf of the tree.
type Element struct {
	SNPs     []string
	STRCount float64
	// Person may be a real person from sample data
	// or a virtual ancestor (modal haplotype).
	Person *genetic.Person
}

// Sample represents the genetic sample of an individual.
// This is a leaf of the phylogenetic tree.
type Sample struct {
	Element
	// ID od this sample, usually the kit number.
	ID string
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
	// TMRCA_STR is the time to the most recent common Ancestor
	// for all downstream samples.
	TMRCA_STR float64
}

func (e *Element) AddSNP(name string) {
	if e.SNPs == nil {
		e.SNPs = make([]string, 0)
	}
	e.SNPs = append(e.SNPs, name)
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
	if e.STRCount > 0 {
		if hasWritten {
			buffer.WriteString(fmt.Sprintf(", STR-Count: %.0f", e.STRCount))
		} else {
			buffer.WriteString(fmt.Sprintf("STR-Count: %.0f", e.STRCount))
		}
	}
	return buffer.String()
}

func (s *Sample) String() string {
	if s.Element.String() != "" {
		return fmt.Sprintf("id:%s, %s", s.ID, s.Element.String())
	} else {
		return fmt.Sprintf("id:%s", s.ID)
	}
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

func (c *Clade) CalculateDistances(mutationRates genetic.YstrMarkers, distance genetic.DistanceFunc) {
	if c.Person == nil {
		return
	}
	// Calculate genetic distances between modal haplotype
	// and samples or subclades.
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
func (c *Clade) CalculateAge(gentime, calibration float64) {
	var cal = gentime * calibration
	var nChilds float64 = 0
	var count float64 = 0
	// Count STR mutations for samples.
	for i, _ := range c.Samples {
		if c.Samples[i].STRCount > 0 {
			count += c.Samples[i].STRCount
			nChilds++
		}
	}
	// Count STR mutations for subclades.
	for i, _ := range c.Subclades {
		c.Subclades[i].CalculateAge(gentime, cal)
		subcladeSTRs := c.Subclades[i].STRCount + c.Subclades[i].STRCountDownstream
		if subcladeSTRs > 0 {
			count += subcladeSTRs
			nChilds++
		}
	}
	// Calculate average number of mutations.
	if nChilds > 0 {
		c.STRCountDownstream = count / nChilds
		c.TMRCA_STR = c.STRCountDownstream * cal
	}
	c.AgeSTR = (c.STRCount + c.STRCountDownstream) * cal
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
	if c.STRCountDownstream > 0 {
		buffer.WriteString(
			fmt.Sprintf(", STRs Downstream: %.f, formed: %.f, TMRCA: %.f",
				c.STRCountDownstream, c.AgeSTR, c.TMRCA_STR))
	}
	buffer.WriteString("\n")

	// Write Samples.
	for _, sample := range c.Samples {
		for i := 0; i < indent+1; i++ {
			buffer.WriteString("\t")
		}
		buffer.WriteString(sample.String())
		buffer.WriteString("\n")
	}
	// Write Subclades.
	for _, clade := range c.Subclades {
		clade.prettyPrint(buffer, indent+1)
	}
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
	for i, _ := range lines {
		childIndent := lines[0].indent
		switch {
		case lines[i].indent < childIndent:
			// Return if indentation shows beginning of next block.
			return nil
		case lines[i].indent == childIndent:
			// Parse child elements.
			if strings.Contains(lines[i].text, "id:") {
				// Child is Sample element.
				sample, err := newSample(lines[i].text)
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

// newClade creates a new Clade from a textual representation.
// Format: SNP1, SNP2, STR-Count: 11
// "STR-Count:" is optional.
func newClade(text string) (Clade, error) {
	var result Clade
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
		default:
			result.AddSNP(token)
		}
	}
	return result, nil
}

// newSample creates a new Sample from a textual representation.
// Format: id:SampleID, SNP1, SNP2, STR-Count: 11
// Only the "id:" field is mandatory.
func newSample(text string) (Sample, error) {
	var result Sample
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
