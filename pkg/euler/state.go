package euler

import "errors"

// Physical errors that can be produced while converting or updating states.
var (
	// ErrNegativeDensity is returned when a conserved state implies rho <= 0.
	ErrNegativeDensity = errors.New("euler: non-positive density")
	// ErrNegativePressure is returned when a conserved state implies p <= 0.
	ErrNegativePressure = errors.New("euler: non-positive pressure")
)

// Primitive describes the primitive variables of a 1-D ideal gas cell:
// density Rho, velocity U, pressure P.
type Primitive struct {
	Rho float64
	U   float64
	P   float64
}

// Conserved describes the conserved variables: mass R, momentum M,
// total energy per unit volume E.
type Conserved struct {
	R float64
	M float64
	E float64
}

// EOS holds the material constants of a calorically perfect gas.
type EOS struct {
	Gamma float64 // ratio of specific heats, > 1
}

// NewEOS returns an EOS with the given ratio of specific heats.
func NewEOS(gamma float64) *EOS { return &EOS{Gamma: gamma} }

// PrimToCons converts primitive to conserved variables.
func (e *EOS) PrimToCons(q Primitive) Conserved {
	return Conserved{
		R: q.Rho,
		M: q.Rho * q.U,
		E: q.P/(e.Gamma-1.0) + 0.5*q.Rho*q.U*q.U,
	}
}

// ConsToPrim converts conserved to primitive variables. It returns
// ErrNegativeDensity/ErrNegativePressure for non-physical states.
func (e *EOS) ConsToPrim(u Conserved) (Primitive, error) {
	if !(u.R > 0.0) || isNaN(u.R) {
		return Primitive{}, ErrNegativeDensity
	}
	vel := u.M / u.R
	kin := 0.5 * u.M * u.M / u.R
	press := (e.Gamma - 1.0) * (u.E - kin)
	if !(press > 0.0) || isNaN(press) {
		return Primitive{}, ErrNegativePressure
	}
	return Primitive{Rho: u.R, U: vel, P: press}, nil
}

// PrimFlux returns the Euler flux evaluated from a primitive state.
func (e *EOS) PrimFlux(q Primitive) Conserved {
	ru := q.Rho * q.U
	return Conserved{
		R: ru,
		M: ru*q.U + q.P,
		E: q.U * (e.PrimToCons(q).E + q.P),
	}
}

// Flux returns the Euler flux F(U) for a conserved state. The state must be
// physical; callers validate through ConsToPrim before invoking it.
func (e *EOS) Flux(u Conserved) Conserved {
	q, err := e.ConsToPrim(u)
	if err != nil {
		return Conserved{}
	}
	return e.PrimFlux(q)
}

// SoundSpeed returns the sound speed sqrt(gamma p / rho).
func (e *EOS) SoundSpeed(q Primitive) float64 {
	return sqrt(e.Gamma * q.P / q.Rho)
}

func isNaN(x float64) bool { return x != x }

func sqrt(x float64) float64 {
	// Local sqrt wrapper keeps all numerical kernels in one package.
	return sqrtMath(x)
}
