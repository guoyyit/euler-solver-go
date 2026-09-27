package euler

// RefineCriterion selects the feature the refinement estimator senses.
type RefineCriterion string

const (
	RefineDensityGradient  RefineCriterion = "density-gradient"
	RefinePressureGradient RefineCriterion = "pressure-gradient"
)

// RefineConfig configures adaptive mesh refinement and coarsening. MinLevel
// and MaxLevel are RELATIVE to the initial tessellation (0 = level-0 cells);
// MaxLevel=3 allows up to 8x refinement.
type RefineConfig struct {
	Enabled   bool
	MinLevel  int
	MaxLevel  int
	RefineAt  float64 // score strictly above which a cell is split
	CoarsenAt float64 // score at or below which a calm sibling pair merges
	Criterion RefineCriterion
	Epsilon   float64 // normalization floor (default 1e-12)
}

func (c RefineConfig) eps() float64 {
	if c.Epsilon > 0 {
		return c.Epsilon
	}
	return 1e-12
}
func (c RefineConfig) minPhys() int { return RootLevel + c.MinLevel }
func (c RefineConfig) maxPhys() int { return RootLevel + c.MaxLevel }

// RefineScore is a scale-free jump indicator for a cell given its immediate
// primitive neighbors:
//
//	(2|fc-fl| + 2|fc-fr|) / (|fl|+2|fc|+|fr| + eps(1+|fc|))
//
// It is O(dx |grad f|/(|f|+eps)) for smooth data and O(1) at a shock, contact
// or rarefaction, so a fixed RefineAt threshold selects only genuinely
// discontinuous regions and the indicator converges under mesh refinement.
func RefineScore(g Gas, cfg RefineConfig, left, center, right Primitive) float64 {
	sel := func(p Primitive) float64 {
		if cfg.Criterion == RefinePressureGradient {
			return p.P
		}
		return p.Rho
	}
	fl, fc, fr := sel(left), sel(center), sel(right)
	num := 2*absF(fc-fl) + 2*absF(fc-fr)
	den := absF(fl) + 2*absF(fc) + absF(fr) + cfg.eps()*(1+absF(fc))
	return num / den
}

// span is a prospective leaf: either one original cell (a>=0,b<0) or a merged
// sibling pair (a,b are the two leaf indices). Lv is its target physical level.
type span struct {
	a, b   int
	lv     int
	id     uint64
	x0, x1 float64
}

// Regrid performs one deterministic, conservative adaptation pass:
//
//  1. compute refinement scores (periodic-aware / one-sided at reflecting);
//  2. mark high-score cells to split and calm sibling pairs to merge;
//  3. build prospective spans, then iterate a 2:1 size-balance closure that
//     only RAISES span levels (a merged span raised to its child level simply
//     un-merges; a cell span raised adds a split), until no two neighboring
//     spans differ by more than one level;
//  4. materialize left-to-right.
//
// Boundary cells are locked to MinLevel: for periodic boundaries this keeps
// the two sides of the domain cut geometrically identical (an exact periodic
// interface), and for reflecting walls it keeps a clean symmetric stencil.
//
// Splits conserve integrals exactly; merges use exact volume averaging, so the
// operation conserves mass, momentum and energy to machine precision.
func Regrid(g Gas, m *Grid, periodic bool, cfg RefineConfig) error {
	if !cfg.Enabled {
		return nil
	}
	if cfg.MaxLevel < cfg.MinLevel {
		cfg.MaxLevel = cfg.MinLevel
	}
	m.sortCells()
	n := len(m.Cells)
	if n < 2 {
		return nil
	}
	prim := make([]Primitive, n)
	for i, c := range m.Cells {
		p, err := g.ConservativeToPrimitive(c.U)
		if err != nil {
			return err
		}
		prim[i] = p
	}
	nb := func(i, d int) Primitive {
		j := i + d
		if periodic {
			return prim[(j+n)%n]
		}
		if j < 0 || j >= n {
			return prim[i] // one-sided at reflecting/free edges
		}
		return prim[j]
	}
	scores := make([]float64, n)
	for i := range m.Cells {
		scores[i] = RefineScore(g, cfg, nb(i, -1), prim[i], nb(i, +1))
	}
	mergeParent := map[uint64]bool{}
	splitWant := map[uint64]bool{}
	for i, c := range m.Cells {
		edge := i == 0 || i == n-1
		maxAllowed := cfg.maxPhys()
		if edge {
			maxAllowed = cfg.minPhys() // boundary lock
		}
		if c.Level >= cfg.minPhys() && c.Level < maxAllowed && scores[i] > cfg.RefineAt {
			splitWant[c.ID] = true
		}
		if !edge && c.Level > cfg.minPhys() && i+1 < n {
			bc := m.Cells[i+1]
			if bc.Level == c.Level && ParentID(bc.ID) == ParentID(c.ID) &&
				scores[i] <= cfg.CoarsenAt && scores[i+1] <= cfg.CoarsenAt {
				mergeParent[ParentID(c.ID)] = true
			}
		}
	}

	// Build prospective spans.
	spans := make([]span, 0, n)
	for i := 0; i < n; {
		c := m.Cells[i]
		if i+1 < n {
			bc := m.Cells[i+1]
			if bc.Level == c.Level && ParentID(bc.ID) == ParentID(c.ID) &&
				mergeParent[ParentID(c.ID)] {
				spans = append(spans, span{a: i, b: i + 1, lv: c.Level - 1,
					id: ParentID(c.ID), x0: c.X0, x1: bc.X1})
				i += 2
				continue
			}
		}
		lv := c.Level
		if splitWant[c.ID] {
			lv++
		}
		spans = append(spans, span{a: i, b: -1, lv: lv, id: c.ID, x0: c.X0, x1: c.X1})
		i++
	}

	// 2:1 closure (raises levels only).
	ns := len(spans)
	getlv := func(k int) int {
		if periodic {
			return spans[(k+ns)%ns].lv
		}
		if k < 0 {
			return spans[0].lv
		}
		if k >= ns {
			return spans[ns-1].lv
		}
		return spans[k].lv
	}
	for {
		changed := false
		for k := range spans {
			lv := spans[k].lv
			if getlv(k-1)-lv >= 2 || getlv(k+1)-lv >= 2 {
				spans[k].lv++
				changed = true
			}
		}
		if !changed {
			break
		}
	}

	// Materialize.
	out := make([]*Cell, 0, n*2)
	for _, sp := range spans {
		c := m.Cells[sp.a]
		if sp.b >= 0 {
			if sp.lv >= c.Level {
				// Merge was cancelled by the closure: emit the two children.
				// A further split (sp.lv > c.Level) is handled below per child.
				bc := m.Cells[sp.b]
				if sp.lv > c.Level {
					if err := appendSplit(g, m, &out, sp.a); err != nil {
						return err
					}
					if err := appendSplit(g, m, &out, sp.b); err != nil {
						return err
					}
				} else {
					out = append(out, c, bc)
				}
				continue
			}
			// genuine merge
			pc, err := mergeCells(g, c, m.Cells[sp.b])
			if err != nil {
				return err
			}
			out = append(out, pc)
			continue
		}
		if sp.lv > c.Level {
			if err := appendSplit(g, m, &out, sp.a); err != nil {
				return err
			}
			continue
		}
		out = append(out, c)
	}
	m.Cells = out
	m.sortCells()
	return nil
}

// appendSplit materializes the split of leaf index idx (one level), using a
// gradient-based conservative split when real in-grid neighbors exist.
func appendSplit(g Gas, m *Grid, out *[]*Cell, idx int) error {
	cells := m.Cells
	var l, r *Cell
	var err error
	if idx > 0 && idx < len(cells)-1 {
		l, r, err = splitCellGrad(g, cells[idx], cells[idx-1], cells[idx+1])
	} else {
		l, r, err = splitCell(g, cells[idx])
	}
	if err != nil {
		return err
	}
	*out = append(*out, l, r)
	return nil
}
