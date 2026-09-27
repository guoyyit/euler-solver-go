package euler

import (
	"errors"
	"testing"
)

// A hand-built grid carrying a non-physical conserved state must be rejected
// when the solver first converts it to primitives (during dt/RHS evaluation),
// rather than producing garbage.
func TestAdvanceRejectsNegativeState(t *testing.T) {
	cfg := testConfig(64, BoundaryPeriodic)
	sim, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	// Corrupt one cell: kinetic energy exceeds total energy -> negative pressure.
	c := sim.Grid.Cells[10]
	c.U.Rho = 1.0
	c.U.RhoU = 10.0
	c.U.E = 1.0
	if err := sim.Advance(); !errors.Is(err, ErrNonPhysical) {
		t.Fatalf("expected ErrNonPhysical, got %v", err)
	}
}
