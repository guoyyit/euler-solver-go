package euler

import "fmt"

// NGhost is the number of ghost layers on each side. MUSCL reconstruction
// needs one neighbour on each side; two layers keep regridding simple.
const NGhost = 2

// Interval is a contiguous run of cells [Lo,Hi) in a packed array.
type Interval struct {
	Lo int
	Hi int
}

// Empty reports whether the interval contains no cells.
func (iv Interval) Empty() bool { return iv.Hi <= iv.Lo }

// Level holds one refinement level. dX is the uniform mesh spacing. U/Q are
// packed arrays of length nSpan+2*NGhost spanning IV (which may include
// coarse-injected gap cells). Active lists the fine-global intervals that are
// genuinely refined; everything else inside IV is a masked gap that must not
// enter totals, restriction or reflux.
type Level struct {
	DX     float64
	IV     Interval   // allocated (span) range on this level
	Active []Interval // genuinely refined intervals (fine global scale)
	U      []Conserved
	Q      []Primitive
	Time   float64
}

// NewLevel allocates a level covering global indices [lo,hi) with spacing dx.
func NewLevel(lo, hi int, dx float64) *Level {
	n := hi - lo
	return &Level{
		DX: dx,
		IV: Interval{Lo: lo, Hi: hi},
		U:  make([]Conserved, n+2*NGhost),
		Q:  make([]Primitive, n+2*NGhost),
	}
}

// NCells returns the number of interior cells.
func (lv *Level) NCells() int { return lv.IV.Hi - lv.IV.Lo }

// at maps a global cell index i inside IV to its packed-array position.
func (lv *Level) at(i int) int { return NGhost + i - lv.IV.Lo }

// Mesh describes the refinement hierarchy; level l has spacing
// X1-X0 / (NBase * Ratio^l).
type Mesh struct {
	Levels []*Level
	Ratio  int
	X0, X1 float64
	NBase  int
}

// BaseDX returns the level-0 mesh spacing.
func (m *Mesh) BaseDX() float64 { return (m.X1 - m.X0) / float64(m.NBase) }

// DX returns the spacing on level l.
func (m *Mesh) DX(l int) float64 { return m.BaseDX() / powInt(m.Ratio, l) }

// CellCenter returns the coordinate of global cell i at spacing dx.
func (m *Mesh) CellCenter(i int, dx float64) float64 {
	return m.X0 + (float64(i)+0.5)*dx
}

func powInt(base, exp int) float64 {
	r := 1
	for k := 0; k < exp; k++ {
		r *= base
	}
	return float64(r)
}

// syncQ recomputes every primitive state from the conserved states and
// returns an error at the first non-physical cell.
func (lv *Level) syncQ(eos *EOS) error {
	for p := range lv.U {
		q, err := eos.ConsToPrim(lv.U[p])
		if err != nil {
			return fmt.Errorf("level dx=%v packed cell %d: %w", lv.DX, p, err)
		}
		lv.Q[p] = q
	}
	return nil
}
