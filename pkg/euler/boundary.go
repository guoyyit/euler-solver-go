package euler

// BoundaryKind identifies how a physical boundary is treated.
type BoundaryKind int

const (
	// BoundaryPeriodic identifies opposite physical edges, so waves leave one
	// side and re-enter the other. Total conserved quantities are unchanged.
	BoundaryPeriodic BoundaryKind = iota
	// BoundaryReflecting is a solid, impermeable wall implemented by ghost
	// mirroring: rho and p are reflected identically while u changes sign,
	// giving zero face-normal velocity and zero boundary mass/energy flux.
	BoundaryReflecting
)

// FillGhost overwrites the nG ghost layers outside interior [lo,hi) of the
// packed primitive array w. The operation is an exact reflection/copy of the
// interior states, hence introduces no extrapolation and no oscillation at
// the boundary.
func FillGhost(w []Primitive, lo, hi, nG int, kind BoundaryKind) {
	switch kind {
	case BoundaryPeriodic:
		for k := 1; k <= nG; k++ {
			w[lo-k] = w[hi-k]     // left ghost <- wrapped right interior
			w[hi+k-1] = w[lo+k-1] // right ghost <- wrapped left interior
		}
	case BoundaryReflecting:
		for k := 1; k <= nG; k++ {
			l := w[lo+k-1]
			w[lo-k] = Primitive{Rho: l.Rho, U: -l.U, P: l.P}
			r := w[hi-k]
			w[hi+k-1] = Primitive{Rho: r.Rho, U: -r.U, P: r.P}
		}
	}
}

// FillGhostCons is the conserved-variable counterpart of FillGhost. At a
// reflecting wall momentum changes sign while mass and energy are mirrored.
func FillGhostCons(u []Conserved, lo, hi, nG int, kind BoundaryKind) {
	switch kind {
	case BoundaryPeriodic:
		for k := 1; k <= nG; k++ {
			u[lo-k] = u[hi-k]
			u[hi+k-1] = u[lo+k-1]
		}
	case BoundaryReflecting:
		for k := 1; k <= nG; k++ {
			l := u[lo+k-1]
			u[lo-k] = Conserved{R: l.R, M: -l.M, E: l.E}
			r := u[hi-k]
			u[hi+k-1] = Conserved{R: r.R, M: -r.M, E: r.E}
		}
	}
}
