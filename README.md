# Illiad — Distance, Angle and Dihedral Calculator (MDAnalysis)

Calculates atomic distances, group distances, angles and dihedral angles along a molecular dynamics trajectory, from a multi-frame PDB file.

## Context
A Python rewrite of an equivalent tool originally written in C++ — the `.dat` output format and the group accumulation logic reproduce the `output_dist_stat` and `output_grp_stat` routines of the earlier version, so existing notebooks and plots keep reading the same files.

## What it does
- Sweeps every frame of the trajectory and computes, for each residue other than the reference one, the distance matrix between all reference atoms and all atoms of the other residue
- Accumulates the mean, standard deviation and frequency of each atom pair across the frames
- Groups chemically equivalent atoms (e.g. the hydrogens of a methyl) and computes the same statistics per group pair, when a groups file is provided
- Computes angles (atom triplets) and dihedral angles (quadruplets), in degrees, frame by frame
- Computes the frequency of occurrence of each angle and dihedral value (histogram)
- Allows restricting the analysis to a time window and to a maximum distance

## Requirements
- Python 3
- NumPy
- MDAnalysis

```bash
pip install numpy MDAnalysis
```

## Usage
```bash
python src/illiad.py <file.pdb> -r <reference_residue> [options]
```

### Arguments
- `pdb` — input PDB file (required)
- `-r, --reference` — name of the reference residue (required)
- `-b, --begin` — start time in ps (default: 0)
- `-e, --end` — end time in ps (default: -1, no limit)
- `-d, --max-dist` — maximum distance in Å; distances greater than or equal to the limit are discarded (default: -1, no limit)
- `-g, --groups` — groups definition file
- `-o, --output` — output file prefix (default: PDB filename without extension)
- `--bin-width` — bin width in degrees for the angle and dihedral frequency histograms (default: 10)
- `--dt` — time between frames in ps; frame *i* is at time *i*·dt and `-b`/`-e` use that time. Without it, the trajectory time is used, which for multi-model PDBs is MDAnalysis' 1 ps default (so `-b`/`-e` are effectively frame indices)
- `--db` — SQLite file to also store the results in (created if missing); each run is added to it, so several runs can be compared (see [Database](#database))
- `--all` — declared, but has no effect yet (see Limitations)

### Example
```bash
python src/illiad.py teste_dados/traj10frames.pdb -r SEH -g teste_dados/grupos.txt -o traj
```

## Groups file
A text file split into sections. Every section whose name is **not** `angle`, `dihedral` or `energy` is read as a residue; inside it, each line starts with the **group label** followed by the names of the atoms belonging to that group. Labels are assigned an integer id in order of appearance.

```
[ SEH ]
1    H0G H0H
3    H0F
a    H0M H0Q H0N H0O H0P

[ CMB ]
1    H3Q H3X H44 H30
2    H3J H3R H3Y H45

[ angle ]
CMB O0R  SEH C08  SEH H0P

[ dihedral ]
SEH H0O  SEH C07  SEH C08  SEH H0P
```

In the `[angle]` and `[dihedral]` sections, each atom can be given in two ways:
- by its atom index in the PDB (e.g. `42`)
- by the residue name + atom name pair (e.g. `SEH C08`)

Lines with an atom count other than 3 (angle) or 4 (dihedral) are skipped with a warning in the terminal, as are selections that match no atom.

## Output
- **`<prefix>_H-dist_<RESIDUE>.dat`** — one file per residue: atom–atom distances, with the columns `Id1 Label1 Id2 Label2 Dist Sigma Frequency`
- **`<prefix>_H-dist_<RESIDUE>_GROUP.dat`** — same columns, aggregated per group pair (only written when a groups file is given)
- **`<prefix>_angles.dat`** — mean and standard deviation of each angle, followed by a second block with each angle's value frame by frame
- **`<prefix>_dihedrals.dat`** — the same, for the dihedrals

- **`<prefix>_angles_hist.dat`** — frequency of occurrence of each angle: bins of `--bin-width` degrees from 0 to 180°, with the columns `Bin_start Bin_end Bin_center` followed by `Angle_<i>_count` and `Angle_<i>_freq(%)` for each angle
- **`<prefix>_dihedrals_hist.dat`** — the same, for the dihedrals, with bins from −180 to 180°

If the bin width does not divide the range evenly, the last bin is narrower.

The angle and dihedral files carry **two headers**: first `Angle_ID / Average / StdDev` (one line per angle), then `Frame Angle_0 Angle_1 ...` (one line per frame). When reading these files the two blocks must be handled separately.

## Database
With `--db <file>.db` the results are also stored in a SQLite database (code in `src/illiad_db.py`). The `.dat` files are still written. Running with the same file again adds a new run instead of overwriting:

```bash
python src/illiad.py teste_dados/teste.pdb -r SEH -g teste_dados/grupos.txt -d 3 -o teste_d3 --db illiad.db
```

Tables, all linked to `runs` by `run_id`:
- **`runs`** — one row per run: input files, `-b`/`-e`, `-d`, `--dt`, `--bin-width`, number of frames processed and date
- **`atom_distances`** / **`group_distances`** — the same columns as the `_H-dist_*.dat` files, plus `resname`
- **`geometry_stats`** — mean and standard deviation of each angle (`kind = 'angle'`) and dihedral (`kind = 'dihedral'`), with the atoms that define it
- **`geometry_series`** — value of each angle/dihedral frame by frame, with the frame time in ps
- **`histograms`** — frequency of occurrence, the same data as the `_hist.dat` files

Example queries:

```sql
-- runs stored so far
SELECT id, pdb, max_dist, n_frames, created_at FROM runs;

-- most frequent group contacts in run 2
SELECT label1, label2, dist, frequency FROM group_distances
WHERE run_id = 2 ORDER BY frequency DESC;

-- frames in which angle 0 is between 100 and 140 degrees
SELECT frame, time_ps, value FROM geometry_series
WHERE run_id = 2 AND kind = 'angle' AND geom_id = 0 AND value BETWEEN 100 AND 140;
```

From Python or the notebook: `pd.read_sql("SELECT ...", sqlite3.connect("illiad.db"))`. To browse the tables without writing SQL, [DB Browser for SQLite](https://sqlitebrowser.org/) opens the file directly.

## Test data
In `teste_dados/`:
- `traj10frames.pdb` — a 10-frame slice, for quick iteration
- `teste.pdb` — the full trajectory (~177 MB)
- `f57360_menor_ep.pdb` — lowest potential energy frame
- `grupos.txt` — example groups file (residues SEH and CMB, with one angle and two dihedrals)
- `illiad_md00.ipynb` — analysis notebook: reads the group `.dat` into a DataFrame, reorders the x-axis labels and plots the distances with a normalized color scale

## Limitations / next steps
- The `--all` flag is declared in argparse but not used yet: per-frame distances are not saved, only the means
- The `[energy]` section is recognized by the groups file parser, but no energy calculation has been implemented
- No automated tests
- Reads PDB only — other trajectory formats supported by MDAnalysis have not been tested
