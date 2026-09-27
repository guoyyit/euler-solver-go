package euler

import "fmt"

// prolongCell splits coarse conserved cell C into r fine conserved children by
// piecewise-constant injection of the conserved density. Because each child
// width is dx/r and all r children share the parent value,
// sum_children U_child * (dx/r) = U_parent * dx exactly: mass, momentum and
// energy are conserved bit-for-bit (multiplication then integer repeated
// addition by the caller reconstructs the parent total).
func prolongCell(c Conserved, r int) []Conserved {
	out := make([]Conserved, r)
	for k := range out {
		out[k] = c
	}
	return out
}

// restrictChildren averages r fine conserved siblings into one coarse cell.
// The averaging weights sum to 1 exactly, so
// U_coarse * dx = sum_k U_fine_k * dx/r.
func restrictChildren(kids []Conserved) Conserved {
	var c Conserved
	w := 1.0 / float64(len(kids))
	for _, u := range kids {
		c.R += u.R
		c.M += u.M
		c.E += u.E
	}
	c.R *= w
	c.M *= w
	c.E *= w
	return c
}

// injectLevel allocates/populates fine level l+1 over intervals ivs (given on
// the fine global index scale), drawing data from coarse parent cells. New
// cells are filled conservatively; cells already present keep their data.
func injectLevel(m *Mesh, coarse, fine *Level, ivs []Interval, eos *EOS, bc BoundaryKind) error {
	r := m.Ratio
	for _, iv := range ivs {
		if (iv.Lo%r) != 0 || (iv.Hi%r) != 0 {
			return fmt.Errorf("fine interval %v not ratio-aligned", iv)
		}
		for fi := iv.Lo; fi < iv.Hi; fi++ {
			ci := fi / r
			if ci < coarse.IV.Lo || ci >= coarse.IV.Hi {
				continue // fine buffer region without a local parent (filled later by BC/coarse elsewhere)
			}
			kids := prolongCell(coarse.U[coarse.at(ci)], r)
			fine.U[fine.at(fi)] = kids[fi%r]
		}
	}
	FillGhostCons(fine.U, NGhost, NGhost+fine.NCells(), NGhost, bc)
	if err := fine.syncQ(eos); err != nil {
		return err
	}
	return nil
}
