package euler

import (
	"errors"
	"fmt"
	"math"
)

// Conserved variables: density, momentum, total energy (per unit volume).
// The state vector is U = [rho, rho*u, E]^T with E = rho*(e + 0.5*u^2).
type Conserved struct {
	Rho  float64
	RhoU float64
	E    float64
}

// Primitive variables: density, velocity, pressure.
type Primitive struct {
	Rho float64
	U   float64
	P   float64
}

// Gas is the ideal-gas equation of state p = (gamma-1)*rho*e.
type Gas struct {
	Gamma float64 // ratio of specific heats
}

// ErrNonPhysical is returned whenever a state has non-positive density/pressure.
var ErrNonPhysical = errors.New("non-physical state: density or pressure is not positive")

// NewGas validates the ratio of specific heats.
func NewGas(gamma float64) (Gas, error) {
	if !(gamma > 1) {
		return Gas{}, fmt.Errorf("gamma must be > 1, got %g", gamma)
	}
	return Gas{Gamma: gamma}, nil
}

// PrimitiveToConservative maps (rho,u,p) -> (rho,rho*u,E).
// A non-positive density or pressure is rejected.
func (g Gas) PrimitiveToConservative(p Primitive) (Conserved, error) {
	if !(p.Rho > 0) || !(p.P > 0) {
		return Conserved{}, fmt.Errorf("%w: rho=%g p=%g", ErrNonPhysical, p.Rho, p.P)
	}
	return Conserved{
		Rho:  p.Rho,
		RhoU: p.Rho * p.U,
		E:    p.P/(g.Gamma-1) + 0.5*p.Rho*p.U*p.U,
	}, nil
}

// ConservativeToPrimitive maps (rho,rho*u,E) -> (rho,u,p) and rejects states
// whose density or derived pressure is non-positive.
func (g Gas) ConservativeToPrimitive(c Conserved) (Primitive, error) {
	if !(c.Rho > 0) {
		return Primitive{}, fmt.Errorf("%w: rho=%g", ErrNonPhysical, c.Rho)
	}
	u := c.RhoU / c.Rho
	kin := 0.5 * c.Rho * u * u
	p := (g.Gamma - 1) * (c.E - kin)
	if !(p > 0) {
		return Primitive{}, fmt.Errorf("%w: p=%g (rho=%g u=%g E=%g)", ErrNonPhysical, p, c.Rho, u, c.E)
	}
	return Primitive{Rho: c.Rho, U: u, P: p}, nil
}

// SoundSpeed returns c = sqrt(gamma*p/rho).
func (g Gas) SoundSpeed(p Primitive) float64 {
	return math.Sqrt(g.Gamma * p.P / p.Rho)
}

// TotalEnthalpy returns H = (E + p)/rho.
func (g Gas) TotalEnthalpy(c Conserved, p Primitive) float64 {
	return (c.E + p.P) / c.Rho
}
