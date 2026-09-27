package euler

import (
	"encoding/binary"
	"fmt"
	"io"
)

// snapshotVersion is bumped whenever the on-disk layout changes.
const snapshotVersion = 1

// magic identifies an euler1d checkpoint.
var snapshotMagic = [8]byte{'E', 'U', 'L', 'E', 'R', '1', 'D', 0}

// WriteSnapshot writes s to w in a versioned, little-endian,
// architecture-independent format. The layout contains no maps, pointers or
// time-dependent data, so identical solver state yields identical bytes and an
// exact Resume.
func WriteSnapshot(w io.Writer, s *Solver) error {
	enc := newBinWriter(w)
	enc.bytes(snapshotMagic[:])
	enc.u64(snapshotVersion)

	// Config.
	enc.f64(s.cfg.X0)
	enc.f64(s.cfg.X1)
	enc.i64(int64(s.cfg.NBase))
	enc.i64(int64(s.cfg.Boundary))
	enc.i64(int64(s.cfg.Flux))
	enc.i64(int64(s.cfg.Limiter))
	enc.f64(s.cfg.CFL)
	enc.f64(s.cfg.Gamma)
	enc.i64(int64(s.cfg.AMR.Ratio))
	enc.i64(int64(s.cfg.AMR.MaxLevels))
	enc.i64(int64(s.cfg.AMR.NGhost))
	enc.i64(int64(s.cfg.AMR.MergeBufferCells))
	enc.f64(s.cfg.AMR.RefineThreshold)
	enc.f64(s.cfg.AMR.CoarsenThreshold)
	enc.f64(s.cfg.AMR.CriterionEps)

	// Runtime.
	enc.i64(s.step)
	enc.f64(s.t)
	enc.i64(int64(len(s.mesh.Levels)))
	for _, lv := range s.mesh.Levels {
		enc.i64(int64(lv.IV.Lo))
		enc.i64(int64(lv.IV.Hi))
		enc.f64(lv.DX)
		enc.f64(lv.Time)
		enc.i64(int64(len(lv.Active)))
		for _, iv := range lv.Active {
			enc.i64(int64(iv.Lo))
			enc.i64(int64(iv.Hi))
		}
		n := len(lv.U)
		enc.i64(int64(n))
		for i := 0; i < n; i++ {
			enc.f64(lv.U[i].R)
			enc.f64(lv.U[i].M)
			enc.f64(lv.U[i].E)
		}
	}
	return enc.err
}

// ReadSnapshot rebuilds a solver exactly as it was when written.
func ReadSnapshot(r io.Reader) (*Solver, error) {
	dec := newBinReader(r)
	magic := make([]byte, 8)
	dec.bytes(magic)
	if dec.err == nil && string(magic) != string(snapshotMagic[:]) {
		return nil, fmt.Errorf("snapshot: bad magic")
	}
	v := dec.u64()
	if dec.err != nil {
		return nil, dec.err
	}
	if v != snapshotVersion {
		return nil, fmt.Errorf("snapshot: unsupported version %d", v)
	}

	cfg := Config{}
	cfg.X0 = dec.f64()
	cfg.X1 = dec.f64()
	cfg.NBase = int(dec.i64())
	cfg.Boundary = BoundaryKind(dec.i64())
	cfg.Flux = FluxScheme(dec.i64())
	cfg.Limiter = Limiter(dec.i64())
	cfg.CFL = dec.f64()
	cfg.Gamma = dec.f64()
	cfg.AMR = DefaultAMRConfig()
	cfg.AMR.Ratio = int(dec.i64())
	cfg.AMR.MaxLevels = int(dec.i64())
	cfg.AMR.NGhost = int(dec.i64())
	cfg.AMR.MergeBufferCells = int(dec.i64())
	cfg.AMR.RefineThreshold = dec.f64()
	cfg.AMR.CoarsenThreshold = dec.f64()
	cfg.AMR.CriterionEps = dec.f64()
	if dec.err != nil {
		return nil, dec.err
	}

	s := &Solver{cfg: cfg, eos: NewEOS(cfg.Gamma)}
	s.step = dec.i64()
	s.t = dec.f64()
	nLev := int(dec.i64())
	m := &Mesh{Ratio: cfg.AMR.Ratio, X0: cfg.X0, X1: cfg.X1, NBase: cfg.NBase}
	for l := 0; l < nLev; l++ {
		lo := int(dec.i64())
		hi := int(dec.i64())
		dx := dec.f64()
		lvTime := dec.f64()
		lv := &Level{DX: dx, IV: Interval{Lo: lo, Hi: hi}, Time: lvTime}
		nAct := int(dec.i64())
		for a := 0; a < nAct; a++ {
			lv.Active = append(lv.Active, Interval{Lo: int(dec.i64()), Hi: int(dec.i64())})
		}
		n := int(dec.i64())
		lv.U = make([]Conserved, n)
		lv.Q = make([]Primitive, n)
		for i := 0; i < n; i++ {
			lv.U[i] = Conserved{R: dec.f64(), M: dec.f64(), E: dec.f64()}
		}
		if dec.err != nil {
			return nil, dec.err
		}
		if err := lv.syncQ(s.eos); err != nil {
			return nil, err
		}
		m.Levels = append(m.Levels, lv)
	}
	s.mesh = m
	return s, nil
}

// Snapshot serialises the complete solver state.
func (s *Solver) Snapshot() ([]byte, error) {
	var buf growBuffer
	if err := WriteSnapshot(&buf, s); err != nil {
		return nil, err
	}
	return buf.b, nil
}

// Restore rebuilds a solver from Snapshot bytes.
func Restore(data []byte) (*Solver, error) {
	return ReadSnapshot(&sliceReader{b: data})
}

type binWriter struct {
	w   io.Writer
	err error
}

func newBinWriter(w io.Writer) *binWriter { return &binWriter{w: w} }

func (b *binWriter) raw(p []byte) {
	if b.err != nil {
		return
	}
	_, b.err = b.w.Write(p)
}

func (b *binWriter) u64(v uint64) {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], v)
	b.raw(buf[:])
}

func (b *binWriter) i64(v int64) { b.u64(uint64(v)) }
func (b *binWriter) f64(v float64) {
	b.u64(mathFloat64bits(v))
}
func (b *binWriter) bytes(p []byte) { b.raw(p) }

type binReader struct {
	r   io.Reader
	err error
}

func newBinReader(r io.Reader) *binReader { return &binReader{r: r} }

func (b *binReader) fill(n int) []byte {
	buf := make([]byte, n)
	if b.err != nil {
		return buf
	}
	_, b.err = io.ReadFull(b.r, buf)
	return buf
}

func (b *binReader) u64() uint64 {
	return binary.LittleEndian.Uint64(b.fill(8))
}
func (b *binReader) i64() int64     { return int64(b.u64()) }
func (b *binReader) f64() float64   { return mathFloat64frombits(b.u64()) }
func (b *binReader) bytes(p []byte) { _, b.err = io.ReadFull(b.r, p) }

type growBuffer struct{ b []byte }

func (g *growBuffer) Write(p []byte) (int, error) {
	g.b = append(g.b, p...)
	return len(p), nil
}

type sliceReader struct{ b []byte }

func (s *sliceReader) Read(p []byte) (int, error) {
	if len(s.b) == 0 {
		return 0, io.EOF
	}
	n := copy(p, s.b)
	s.b = s.b[n:]
	return n, nil
}
