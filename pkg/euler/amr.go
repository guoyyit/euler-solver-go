package euler

// RefineCriterion decides where finer grids are needed.
type RefineCriterion interface {
	// Flag returns true for cells (indexed on lv) that need refinement.
	Flag(lv *Level) []bool
}

// AMRConfig controls refinement, coarsening and interpolation behaviour.
type AMRConfig struct {
	Ratio            int
	MaxLevels        int
	NGhost           int
	MergeBufferCells int
	RefineThreshold  float64
	CoarsenThreshold float64
	CriterionEps     float64
	Criterion        RefineCriterion // optional override; nil => built-in jump criterion
}

// DefaultAMRConfig returns conservative, reproducible AMR settings.
func DefaultAMRConfig() AMRConfig {
	return AMRConfig{
		Ratio:            2,
		MaxLevels:        1,
		NGhost:           NGhost,
		MergeBufferCells: 1,
		RefineThreshold:  0.05,
		CoarsenThreshold: 0.01,
		CriterionEps:     1e-12,
	}
}

func (c AMRConfig) withDefaults() AMRConfig {
	d := DefaultAMRConfig()
	if c.Ratio != 0 {
		d.Ratio = c.Ratio
	}
	if c.MaxLevels != 0 {
		d.MaxLevels = c.MaxLevels
	}
	if c.NGhost != 0 {
		d.NGhost = c.NGhost
	}
	if c.MergeBufferCells != 0 {
		d.MergeBufferCells = c.MergeBufferCells
	}
	if c.RefineThreshold != 0 {
		d.RefineThreshold = c.RefineThreshold
	}
	if c.CoarsenThreshold != 0 {
		d.CoarsenThreshold = c.CoarsenThreshold
	}
	if c.CriterionEps != 0 {
		d.CriterionEps = c.CriterionEps
	}
	if c.Criterion != nil {
		d.Criterion = c.Criterion
	}
	return d
}

// uniqueCellsAMR is implemented in amr_step.go / amr_couple.go.
