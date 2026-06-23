# ILLIAD

A tool for calculating distances, angles, and dihedral angles in molecular dynamics trajectories using MDAnalysis.

## Features

- **Atomic distances**: Calculates distances between atoms belonging to different residues
- **Group statistics**: Groups atoms together and calculates average distances between groups
- **Angles**: Calculates angles between three atoms
- **Dihedral angles**: Calculates dihedral angles between four atoms
- **Temporal and spatial filters**: Select specific time windows and maximum distance cutoffs

## Requirements

- Python 3.x
- NumPy
- MDAnalysis

```bash
pip install numpy MDAnalysis
```

## Usage

```bash
python testeCodeConverter1.py <pdb_file> -r <reference_residue> [options]
```

### Required Arguments

- `pdb`: Input PDB file
- `-r, --reference`: Name of the reference residue

### Optional Arguments

- `-b, --begin`: Start time in ps (default: 0)
- `-e, --end`: End time in ps (default: -1, no limit)
- `-d, --max-dist`: Maximum distance in Ångström (default: -1, no limit)
- `-g, --groups`: Group definition file
- `-o, --output`: Output file prefix (default: PDB file name)
- `--all`: Save distances for all frames

### Group Definition File

A configuration file that defines atom groups, custom angles, and dihedrals. Each `[Section]` corresponds either to a residue name (defining atom groups for that residue) or to one of the reserved keywords `angle` / `dihedral`:

```
[GroupName1]
atom1 atom2 atom3

[GroupName2]
atom4 atom5

[angle]
atom1 atom2 atom3

[dihedral]
atom1 atom2 atom3 atom4
```

- Lines under a residue-name section assign each listed atom name to that group (the first token is the group label, the remaining tokens are atom names belonging to it).
- Lines under `[angle]` must contain exactly 3 atoms; lines under `[dihedral]` must contain exactly 4 atoms.

#### Atom Specification

Atoms (in `[angle]` and `[dihedral]` sections) can be specified in two ways:

- **By atom index**: a plain integer, matched via `id <n>` (e.g. `42`)
- **By residue + atom name**: two tokens, `<resname> <atomname>` (e.g. `SEH C08`)

These two styles can be mixed freely within the same line.

## Outputs

1. **`*_H-dist_<residue>.dat`**: Atom-to-atom distances, per residue
2. **`*_H-dist_<residue>_GROUP.dat`**: Average distances between defined groups
3. **`*_angles.dat`**: Angle averages plus full frame-by-frame values
4. **`*_dihedrals.dat`**: Dihedral angle averages plus full frame-by-frame values

## Example

```bash
python testeCodeConverter1.py trajectory.pdb -r SOL -g groups.cfg -o output -b 0 -e 100
```

## Output File Structure

- **`*_H-dist_<residue>.dat`** / **`*_H-dist_<residue>_GROUP.dat`**: each row contains the IDs and labels of the two atoms/groups, the average distance, the standard deviation, and the number of occurrences (frequency) within the cutoff.
- **`*_angles.dat`** / **`*_dihedrals.dat`**: starts with a summary table (ID, average, standard deviation) for each defined angle/dihedral, followed by a frame-by-frame table with one column per angle/dihedral.
