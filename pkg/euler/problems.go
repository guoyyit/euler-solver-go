package euler

// ProblemID names a built-in initial/boundary configuration.
type ProblemID int

const (
	// ProblemSod is the classic Sod shock tube (1,1,1)|(0.125,0,0.1).
	ProblemSod ProblemID = iota
	// ProblemLax is the Lax shock tube with a right-moving contact and shock.
	ProblemLax
	// ProblemDoubleRarefaction drives two strong rarefactions apart; it is a
	// stringent positivity test (near-vacuum in the centre).
	ProblemDoubleRarefaction
)

// Problem describes one standard Riemann configuration on [0,1).
type Problem struct {
	ID       ProblemID
	Gamma    float64
	Boundary BoundaryKind
	XDisc    float64 // discontinuity location
	QL, QR   Primitive
	Init     func(x float64) Primitive
}

func riemannProblem(id ProblemID, gamma float64, bc BoundaryKind, xd float64, ql, qr Primitive) Problem {
	return Problem{
		ID:       id,
		Gamma:    gamma,
		Boundary: bc,
		XDisc:    xd,
		QL:       ql,
		QR:       qr,
		Init: func(x float64) Primitive {
			if x < xd {
				return ql
			}
			return qr
		},
	}
}

// GetProblem returns the named standard problem.
func GetProblem(id ProblemID) Problem {
	switch id {
	case ProblemLax:
		return riemannProblem(ProblemLax, 1.4, BoundaryPeriodic, 0.5,
			Primitive{Rho: 0.445, U: 0.698, P: 3.528},
			Primitive{Rho: 0.5, U: 0, P: 0.571})
	case ProblemDoubleRarefaction:
		return riemannProblem(ProblemDoubleRarefaction, 1.4, BoundaryReflecting, 0.5,
			Primitive{Rho: 1, U: -2, P: 0.4},
			Primitive{Rho: 1, U: 2, P: 0.4})
	default:
		return riemannProblem(ProblemSod, 1.4, BoundaryPeriodic, 0.5,
			Primitive{Rho: 1, U: 0, P: 1},
			Primitive{Rho: 0.125, U: 0, P: 0.1})
	}
}
