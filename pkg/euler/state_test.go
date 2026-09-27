package euler

import (
	"math"
	"testing"
)

func TestPrimConsRoundtrip(t *testing.T) {
	eos := NewEOS(1.4)
	cases := []Primitive{
		{1, 0, 1},
		{1.7, -2.3, 3.5},
		{0.125, 5, 0.4},
	}
	for _, q := range cases {
		u := eos.PrimToCons(q)
		q2, err := eos.ConsToPrim(u)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if math.Abs(q2.Rho-q.Rho) > 1e-14 || math.Abs(q2.U-q.U) > 1e-14 || math.Abs(q2.P-q.P) > 1e-14 {
			t.Fatalf("roundtrip mismatch: %v -> %v -> %v", q, u, q2)
		}
	}
}

func TestConsToPrimRejectsBadStates(t *testing.T) {
	eos := NewEOS(1.4)
	if _, err := eos.ConsToPrim(Conserved{R: -1, M: 0, E: 2}); err != ErrNegativeDensity {
		t.Fatalf("want ErrNegativeDensity, got %v", err)
	}
	// rho=1, u=0 but total energy negative -> negative pressure.
	if _, err := eos.ConsToPrim(Conserved{R: 1, M: 0, E: -1}); err != ErrNegativePressure {
		t.Fatalf("want ErrNegativePressure, got %v", err)
	}
}

func TestSoundSpeed(t *testing.T) {
	eos := NewEOS(1.4)
	if got := eos.SoundSpeed(Primitive{Rho: 1, P: 1}); got != math.Sqrt(1.4) {
		t.Fatalf("sound speed = %v", got)
	}
}
