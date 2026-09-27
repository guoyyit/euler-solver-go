package euler

// BoundaryType selects a boundary treatment.
type BoundaryType string

const (
	BoundaryPeriodic   BoundaryType = "periodic"
	BoundaryReflecting BoundaryType = "reflecting" // solid wall, u=0
)

// GhostState returns the ghost primitive state mirrored outside a boundary.
//
//	Periodic:   the ghost copies the interior state (the caller maps it to the
//	            cell on the opposite side of the domain).
//	Reflecting: (rho,u,p)_ghost = (rho,-u,p), so the resulting interface flux
//	            has zero mass flux and zero momentum work: the wall carries no
//	            artificial source, and the reconstruction introduces no
//	            non-physical oscillation (it is an exact symmetry).
//
// side is -1 for the left boundary and +1 for the right boundary; it is
// unused for the symmetric reflection but kept for clarity.
func GhostState(b BoundaryType, side int, q Primitive) Primitive {
	switch b {
	case BoundaryReflecting:
		return Primitive{Rho: q.Rho, U: -q.U, P: q.P}
	default:
		return Primitive{Rho: q.Rho, U: q.U, P: q.P}
	}
}
