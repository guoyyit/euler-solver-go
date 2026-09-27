package euler

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// SnapshotFormatVersion is bumped when the on-disk layout changes.
const SnapshotFormatVersion = 1

// Snapshot is a full, deterministic checkpoint of a simulation.
// All floating point fields are serialized verbatim (Go's encoding/json writes
// the shortest round-trip decimal for float64), so a save/load cycle preserves
// every bit and a resumed run is identical to an uninterrupted one.
type Snapshot struct {
	FormatVersion int        `json:"format_version"`
	Config        Config     `json:"config"`
	Time          float64    `json:"time"`
	Step          int64      `json:"step"`
	Cells         []SnapCell `json:"cells"`
}

// SnapCell is the serialized form of a leaf cell.
type SnapCell struct {
	ID    uint64  `json:"id"`
	Level int     `json:"level"`
	X0    float64 `json:"x0"`
	X1    float64 `json:"x1"`
	Rho   float64 `json:"rho"`
	RhoU  float64 `json:"rhou"`
	E     float64 `json:"e"`
}

// Snapshot serializes the current simulation state.
func (s *Simulation) Snapshot() *Snapshot {
	cells := make([]SnapCell, len(s.Grid.Cells))
	for i, c := range s.Grid.Cells {
		cells[i] = SnapCell{
			ID: c.ID, Level: c.Level, X0: c.X0, X1: c.X1,
			Rho: c.U.Rho, RhoU: c.U.RhoU, E: c.U.E,
		}
	}
	return &Snapshot{
		FormatVersion: SnapshotFormatVersion,
		Config:        s.Cfg,
		Time:          s.Time,
		Step:          s.Step,
		Cells:         cells,
	}
}

// encodeSnapshot returns deterministic JSON with stable field order.
func encodeSnapshot(s *Snapshot) ([]byte, error) {
	return json.MarshalIndent(s, "", "  ")
}

// SaveSnapshot writes a checkpoint atomically: it writes to a temp file in the
// same directory and renames it into place, so a crash can never leave a
// torn/partial checkpoint behind.
func (s *Simulation) SaveSnapshot(path string) error {
	return SaveSnapshot(path, s.Snapshot())
}

// SaveSnapshot writes an arbitrary snapshot atomically.
func SaveSnapshot(path string, snap *Snapshot) error {
	data, err := encodeSnapshot(snap)
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0o644); err != nil {
		return err
	}
	if err := os.Rename(tmp, path); err != nil {
		return err
	}
	return nil
}

// LoadSnapshot reads and validates a checkpoint.
func LoadSnapshot(path string) (*Snapshot, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var snap Snapshot
	if err := json.Unmarshal(data, &snap); err != nil {
		return nil, err
	}
	if snap.FormatVersion != SnapshotFormatVersion {
		return nil, fmt.Errorf("unsupported snapshot format version %d (want %d)",
			snap.FormatVersion, SnapshotFormatVersion)
	}
	if len(snap.Cells) == 0 {
		return nil, fmt.Errorf("snapshot contains no cells")
	}
	return &snap, nil
}
