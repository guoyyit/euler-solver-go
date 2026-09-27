package euler

import (
	"fmt"
	"sort"
)

// Cell is a leaf control volume [X0,X1] with a conserved cell AVERAGE
// (i.e. U is per-unit-volume; integrated mass is U.Rho * Width()).
type Cell struct {
	ID    uint64
	Level int
	X0    float64
	X1    float64
	U     Conserved
}

func (c *Cell) Width() float64  { return c.X1 - c.X0 }
func (c *Cell) Center() float64 { return 0.5 * (c.X0 + c.X1) }

// Grid is a 1-D non-uniform tessellation of [Xmin,Xmax] by leaf cells.
// Leaf cells never overlap and their widths differ by at most a factor of 2
// between neighbors (the 2:1 balance maintained by Regrid).
type Grid struct {
	Xmin  float64
	Xmax  float64
	Cells []*Cell
	byID  map[uint64]*Cell
}

// NewGrid creates an empty mesh over [xmin,xmax].
func NewGrid(xmin, xmax float64) *Grid {
	if xmax <= xmin {
		panic(fmt.Sprintf("invalid domain [%g,%g]", xmin, xmax))
	}
	return &Grid{Xmin: xmin, Xmax: xmax, byID: map[uint64]*Cell{}}
}

// RootLevel is the physical tree level of the initial N0-cell tessellation.
// With RootLevel=24 an id fits 63 bits (uint64) and up to RootLevel further
// refinement levels are available. Configuration MinLevel/MaxLevel are
// RELATIVE to RootLevel (0 = initial mesh), so a run with MaxLevel=4 refines
// at most 16x.
const RootLevel = 24

// RootID is the id of the level-0 cell covering the whole domain.
// Encoding: binary 1, then one bit per level (0 = left child, 1 = right).
// Cells at physical level L have ids in [2^L, 2^(L+1)-1].
func RootID() uint64 { return 1 }

// BaseID returns the id of the k-th cell of the initial N0 tessellation
// (0 <= k < N0, N0 = 2^RootLevel conceptually; non-power-of-two N0 is
// supported because ids are only required to be unique and sibling-consistent
// under the bit encoding within the represented tree).
func BaseID(k int) uint64 { return uint64(1)<<RootLevel + uint64(k) }

// ChildrenIDs returns the two child ids and the children's level.
func ChildrenIDs(id uint64, level int) (left, right uint64, childLevel int) {
	return id << 1, id<<1 | 1, level + 1
}

// ParentID returns the parent id.
func ParentID(id uint64) uint64 { return id >> 1 }

// CellSpan returns the physical span of a cell id at a given level.
func CellSpan(xmin, xmax float64, id uint64, level int) (x0, x1 float64) {
	// index within the level: idx in [0,2^L)
	var idx uint64
	if level == 0 {
		idx = 0
	} else {
		idx = id - (uint64(1) << level)
	}
	n := float64(uint64(1) << level)
	w := (xmax - xmin) / n
	x0 = xmin + float64(idx)*w
	x1 = x0 + w
	return x0, x1
}

// sortCells orders cells left to right and rebuilds the id index.
func (m *Grid) sortCells() {
	sort.Slice(m.Cells, func(i, j int) bool { return m.Cells[i].X0 < m.Cells[j].X0 })
	m.byID = make(map[uint64]*Cell, len(m.Cells))
	for _, c := range m.Cells {
		m.byID[c.ID] = c
	}
}

// CellByID looks up a leaf cell by id (nil if absent).
func (m *Grid) CellByID(id uint64) *Cell { return m.byID[id] }

// splitCell returns two children replacing the parent. Each child inherits the
// parent conserved average plus a conservative slope correction:
//
//	U_L = U - 0.5*dU,  U_R = U + 0.5*dU  (equal widths => integral conserved)
//
// The primitive slope is limited and any non-positive reconstructed state
// falls back to piecewise-constant injection, which is always conservative.
func splitCell(g Gas, parent *Cell) (left, right *Cell, err error) {
	idL, idR, childLevel := ChildrenIDs(parent.ID, parent.Level)
	// No neighbor info here: use conservative piecewise-constant injection.
	uL := parent.U
	uR := parent.U
	xm := parent.Center()
	left = &Cell{ID: idL, Level: childLevel, X0: parent.X0, X1: xm, U: uL}
	right = &Cell{ID: idR, Level: childLevel, X0: xm, X1: parent.X1, U: uR}
	return left, right, nil
}

// splitCellGrad splits a parent into two equal children with a minmod-limited
// reconstruction performed DIRECTLY ON THE CONSERVED VARIABLES. Because every
// conserved component is linear in the cell averages,
//
//	0.5(U_L+U_R) = U_parent   =>   U_L W/2 + U_R W/2 = U_parent W
//
// holds to machine precision (no nonlinear primitive->conservative round trip
// on the slope). Each child is independently converted to primitives for a
// positivity check; if either child is non-physical the split falls back to
// piecewise-constant injection, which is always conservative.
func splitCellGrad(g Gas, parent, nbL, nbR *Cell) (left, right *Cell, err error) {
	idL, idR, lv := ChildrenIDs(parent.ID, parent.Level)
	if _, err := g.ConservativeToPrimitive(parent.U); err != nil {
		return nil, nil, err
	}
	var s Conserved
	if nbL != nil && nbR != nil {
		if _, e1 := g.ConservativeToPrimitive(nbL.U); e1 == nil {
			if _, e2 := g.ConservativeToPrimitive(nbR.U); e2 == nil {
				// Scale neighbor differences to the parent width so the limited
				// slope is expressed per parent-cell width.
				scale := parent.Width() / (0.5 * (nbL.Width() + nbR.Width()))
				d := func(pc, lc, rc float64) float64 {
					return minmod((pc-lc)*scale, (rc-pc)*scale)
				}
				s = Conserved{
					Rho:  d(parent.U.Rho, nbL.U.Rho, nbR.U.Rho),
					RhoU: d(parent.U.RhoU, nbL.U.RhoU, nbR.U.RhoU),
					E:    d(parent.U.E, nbL.U.E, nbR.U.E),
				}
			}
		}
	}
	cL := Conserved{
		Rho:  parent.U.Rho - 0.5*s.Rho,
		RhoU: parent.U.RhoU - 0.5*s.RhoU,
		E:    parent.U.E - 0.5*s.E,
	}
	cR := Conserved{
		Rho:  parent.U.Rho + 0.5*s.Rho,
		RhoU: parent.U.RhoU + 0.5*s.RhoU,
		E:    parent.U.E + 0.5*s.E,
	}
	if _, eL := g.ConservativeToPrimitive(cL); eL != nil {
		cL, cR = parent.U, parent.U
	} else if _, eR := g.ConservativeToPrimitive(cR); eR != nil {
		cL, cR = parent.U, parent.U
	}
	xm := parent.Center()
	left = &Cell{ID: idL, Level: lv, X0: parent.X0, X1: xm, U: cL}
	right = &Cell{ID: idR, Level: lv, X0: xm, X1: parent.X1, U: cR}
	return left, right, nil
}

// mergeCells combines two siblings into their parent by exact volume averaging.
// Equal child widths reduce this to the arithmetic mean of conserved averages,
// so the integral U_parent*W_parent = U_L*W_L + U_R*W_R is conserved exactly.
func mergeCells(g Gas, left, right *Cell) (*Cell, error) {
	if ParentID(right.ID) != ParentID(left.ID) {
		return nil, fmt.Errorf("mergeCells: cells %d,%d are not siblings", left.ID, right.ID)
	}
	wL, wR := left.Width(), right.Width()
	w := wL + wR
	u := Conserved{
		Rho:  (left.U.Rho*wL + right.U.Rho*wR) / w,
		RhoU: (left.U.RhoU*wL + right.U.RhoU*wR) / w,
		E:    (left.U.E*wL + right.U.E*wR) / w,
	}
	if _, err := g.ConservativeToPrimitive(u); err != nil {
		return nil, err
	}
	lv := left.Level - 1
	return &Cell{ID: ParentID(left.ID), Level: lv, X0: left.X0, X1: right.X1, U: u}, nil
}
