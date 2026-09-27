package euler

import "testing"

func TestFillGhostPeriodic(t *testing.T) {
	lo, hi, nG := 2, 6, 2
	w := make([]Primitive, hi+2*nG-lo)
	for i := lo; i < hi; i++ {
		w[i] = Primitive{Rho: float64(i), U: 0, P: 1}
	}
	FillGhost(w, lo, hi, nG, BoundaryPeriodic)
	if w[0] != w[4] || w[1] != w[5] || w[6] != w[2] || w[7] != w[3] {
		t.Fatalf("periodic ghosts wrong: %v", w)
	}
}

func TestFillGhostReflecting(t *testing.T) {
	lo, hi, nG := 2, 5, 2
	w := make([]Primitive, hi+2*nG-lo)
	w[2] = Primitive{Rho: 1, U: 3, P: 1}
	w[3] = Primitive{Rho: 1, U: 2, P: 1}
	w[4] = Primitive{Rho: 1, U: 1, P: 1}
	FillGhost(w, lo, hi, nG, BoundaryReflecting)
	if w[1].U != -3 || w[0].U != -2 || w[1].Rho != 1 || w[1].P != 1 {
		t.Fatalf("left reflection wrong: %+v", w[:2])
	}
	if w[5].U != -1 || w[6].U != -2 || w[5].Rho != 1 {
		t.Fatalf("right reflection wrong: %+v", w[5:])
	}
}

// At a reflecting wall the face between the first interior cell and its
// mirror ghost is stationary-symmetric, so the mass flux must vanish.
func TestReflectingWallZeroMassFlux(t *testing.T) {
	eos := NewEOS(1.4)
	q := Primitive{Rho: 2, U: 1.5, P: 3}
	ghost := Primitive{Rho: q.Rho, U: -q.U, P: q.P}
	for _, sc := range []FluxScheme{FluxRusanov, FluxHLL} {
		f := NewFlux(sc).Flux(eos, ghost, q)
		if abs(f.R) > 1e-14 {
			t.Fatalf("%v: mass flux through wall = %v", sc, f.R)
		}
	}
}
