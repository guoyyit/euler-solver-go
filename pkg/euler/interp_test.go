package euler

import "testing"

func TestProlongRestrictConservative(t *testing.T) {
	parent := Conserved{R: 1.3, M: -0.7, E: 2.5}
	for _, r := range []int{2, 4} {
		kids := prolongCell(parent, r)
		// Conservative integral check: sum kids * dx/r == parent * dx.
		var sum Conserved
		for _, k := range kids {
			sum.R += k.R / float64(r)
			sum.M += k.M / float64(r)
			sum.E += k.E / float64(r)
		}
		if sum != parent {
			t.Fatalf("r=%d prolong integral %v != parent %v", r, sum, parent)
		}
		back := restrictChildren(kids)
		if back != parent {
			t.Fatalf("r=%d restrict %v != parent %v", r, back, parent)
		}
	}
}

func TestUnionAndPartition(t *testing.T) {
	u := unionIntervals([]Interval{{1, 3}, {3, 5}, {8, 9}, {2, 4}})
	if len(u) != 2 || u[0] != (Interval{1, 5}) || u[1] != (Interval{8, 9}) {
		t.Fatalf("union wrong: %v", u)
	}
	e := expandIntervals([]Interval{{4, 6}}, 1, 10)
	if len(e) != 1 || e[0] != (Interval{3, 7}) {
		t.Fatalf("expand wrong: %v", e)
	}
}
