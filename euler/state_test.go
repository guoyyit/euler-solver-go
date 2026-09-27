package euler

import (
	"errors"
	"math"
	"testing"
)

func TestConservativePrimitiveRoundTrip(t *testing.T) {
	g, _ := NewGas(1.4)
	cases := []Primitive{
		{1, 0, 1},
		{0.125, -3.7, 0.1},
		{2.5, 12.34, 40.7},
		{1e-6, 1e3, 1e-6},
	}
	for _, p := range cases {
		c, err := g.PrimitiveToConservative(p)
		if err != nil {
			t.Fatalf("P->C: %v", err)
		}
		q, err := g.ConservativeToPrimitive(c)
		if err != nil {
			t.Fatalf("C->P: %v", err)
		}
		for _, got := range [3]float64{q.Rho - p.Rho, q.U - p.U, q.P - p.P} {
			if math.Abs(got) > 1e-12*math.Max(1, math.Abs(got))+1e-14 {
				t.Fatalf("round trip mismatch %+v -> %+v", p, q)
			}
		}
	}
}

func TestRejectNegativeState(t *testing.T) {
	g, _ := NewGas(1.4)
	badP := []Primitive{{Rho: -1, P: 1}, {Rho: 0, P: 1}, {Rho: 1, P: -2}, {Rho: 1, P: 0}}
	for _, p := range badP {
		if _, err := g.PrimitiveToConservative(p); !errors.Is(err, ErrNonPhysical) {
			t.Fatalf("expected ErrNonPhysical for %+v, got %v", p, err)
		}
	}
	badC := []Conserved{{Rho: -1, E: 1}, {Rho: 0, E: 1}, {Rho: 1, RhoU: 100, E: 1}}
	for _, c := range badC {
		if _, err := g.ConservativeToPrimitive(c); !errors.Is(err, ErrNonPhysical) {
			t.Fatalf("expected ErrNonPhysical for %+v, got %v", c, err)
		}
	}
}

func TestSoundSpeed(t *testing.T) {
	g, _ := NewGas(1.4)
	if got := g.SoundSpeed(Primitive{Rho: 1, P: 1}); math.Abs(got-math.Sqrt(1.4)) > 1e-15 {
		t.Fatalf("soundspeed = %v", got)
	}
}
