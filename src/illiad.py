import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_dihedrals
from MDAnalysis.analysis import distances as mda_distances
import argparse

def _is_atom_index(token):
    """Returns True if token is a pure integer (atom index)."""
    try:
        int(token)
        return True
    except ValueError:
        return False

def _parse_atom_tokens(tokens, u):
    """
    Parses tokens as atoms, either by index (42) or resname+atomname (SEH C08).
    Returns a list of AtomGroups, one per atom found.
    """
    atoms = []
    i = 0
    while i < len(tokens):
        if _is_atom_index(tokens[i]):
            ag = u.select_atoms(f"id {tokens[i]}")
            if len(ag) == 0:
                print(f"  WARNING: no atom found with id {tokens[i]}")
            else:
                atoms.append(ag)
            i += 1
        else:
            if i + 1 >= len(tokens):
                print(f"  WARNING: expected atomname after resname '{tokens[i]}'")
                break
            resname, atomname = tokens[i], tokens[i + 1]
            ag = u.select_atoms(f"resname {resname} and name {atomname}")
            if len(ag) == 0:
                print(f"  WARNING: no atom found for resname={resname} name={atomname}")
            else:
                atoms.append(ag)
            i += 2
    return atoms


def parse_groups_file(filename, u):
    """
    Parses the groups file and returns:
      configs['groups']  : resname -> {atomname -> group_label}
      configs['types']   : resname -> {group_label -> int_id}   (order of appearance)
      configs['angles']  : list of AtomGroups (3 atoms each)
      configs['dihedrals']: list of AtomGroups (4 atoms each)
    """
    configs = {
        'angles': [],
        'dihedrals': [],
        'groups': {},   # resname -> {atomname -> group_label}
        'types':  {},   # resname -> {group_label -> int_id}
    }

    current_section = None
    rname = None

    if not filename:
        return configs

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('['):
                section = line.replace('[', '').replace(']', '').strip()
                if section in ("angle", "dihedral", "energy"):
                    current_section = section
                else:
                    current_section = "group"
                    rname = section
                    configs['groups'][rname] = {}
                    configs['types'][rname]  = {}
                continue

            tokens = line.split()
            if not tokens:
                continue

            if current_section == "angle":
                atoms = _parse_atom_tokens(tokens, u)
                if len(atoms) == 3:
                    # Store atom indices from the original universe (positions update each frame)
                    indices = [ag.indices[0] for ag in atoms]
                    configs['angles'].append(indices)
                else:
                    print(f"  WARNING: angle line has {len(atoms)} atoms (expected 3): '{line}'")

            elif current_section == "dihedral":
                atoms = _parse_atom_tokens(tokens, u)
                if len(atoms) == 4:
                    indices = [ag.indices[0] for ag in atoms]
                    configs['dihedrals'].append(indices)
                else:
                    print(f"  WARNING: dihedral line has {len(atoms)} atoms (expected 4): '{line}'")

            elif current_section == "group":
                gname    = tokens[0]
                atomnames = tokens[1:]
                # assign integer id in order of appearance (like types[rname][gname]=++gnumber)
                if gname not in configs['types'][rname]:
                    configs['types'][rname][gname] = len(configs['types'][rname]) + 1
                for aname in atomnames:
                    configs['groups'][rname][aname] = gname

    return configs


# ─────────────────────────────────────────────
# Frequency of occurrence (histogram)
# ─────────────────────────────────────────────
def write_histogram(path, results, label, lo, hi, width):
    """
    Writes the frequency of occurrence of each angle/dihedral in bins of `width` degrees
    covering [lo, hi]. One Count and one Freq(%) column per angle/dihedral.
    """
    edges = np.arange(lo, hi, width)
    edges = np.append(edges, hi)   # last bin is narrower when width does not divide the range

    with open(path, 'w') as f:
        header = f"{'Bin_start':>10} {'Bin_end':>10} {'Bin_center':>10}"
        for i in range(len(results)):
            header += f" {f'{label}_{i}_count':>16} {f'{label}_{i}_freq(%)':>18}"
        f.write(header + "\n")

        columns = []
        for data in results:
            counts, _ = np.histogram(data, bins=edges)
            n = counts.sum()
            freqs = 100.0 * counts / n if n > 0 else np.zeros(len(counts))
            columns.append((counts, freqs))

        for b in range(len(edges) - 1):
            start, end = edges[b], edges[b + 1]
            row = f"{start:10.2f} {end:10.2f} {(start + end) / 2:10.2f}"
            for counts, freqs in columns:
                row += f" {counts[b]:16d} {freqs[b]:18.2f}"
            f.write(row + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate distances, angles, and dihedrals using MDAnalysis.")
    parser.add_argument("pdb",                               help="Input PDB file")
    parser.add_argument("-r", "--reference",   required=True, help="Reference residue name")
    parser.add_argument("-b", "--begin",       type=float, default=0,  help="Start time (ps)")
    parser.add_argument("-e", "--end",         type=float, default=-1, help="End time (ps)")
    parser.add_argument("-d", "--max-dist",    type=float, default=-1, help="Maximum distance (Å)")
    parser.add_argument("-g", "--groups",                    help="Groups definition file")
    parser.add_argument("-o", "--output",                    help="Output prefix")
    parser.add_argument("--all", action="store_true",        help="Save all frame distances")
    parser.add_argument("--bin-width", type=float, default=10,
                        help="Bin width (degrees) for the angle/dihedral frequency histograms")
    parser.add_argument("--dt", type=float,
                        help="Time between frames (ps). Frame i is at time i*dt. "
                             "Without it, the trajectory time is used (for multi-model "
                             "PDBs that is MDAnalysis' 1 ps default, i.e. the frame index)")
    args = parser.parse_args()

    if args.bin_width <= 0:
        parser.error("--bin-width must be greater than 0")
    if args.dt is not None and args.dt <= 0:
        parser.error("--dt must be greater than 0")

    # ── Load universe ──────────────────────────────────────────────────────────
    u      = mda.Universe(args.pdb)
    prefix = args.output if args.output else args.pdb.rsplit('.', 1)[0]

    # ── Parse groups ───────────────────────────────────────────────────────────
    configs = parse_groups_file(args.groups, u)
    groups  = configs['groups']   # resname -> {atomname -> group_label}
    types   = configs['types']    # resname -> {group_label -> int_id}

    ref_name     = args.reference
    other_resnames = [r for r in np.unique(u.atoms.resnames) if r != ref_name]

    # ── Selections and pair keys (topology is fixed, so built once) ───────────
    # Each (ref atom, other atom) pair maps to an integer key; per-frame
    # distances are then accumulated with np.bincount (count, sum, sum of squares).
    ref_atoms_all = u.select_atoms(f"resname {ref_name}")
    ref_names_u, ref_codes = np.unique(ref_atoms_all.names, return_inverse=True)

    # group label of each atom as an integer code (-1 = atom not in any group)
    ref_grp_of   = groups.get(ref_name, {})
    ref_grps_u   = sorted(set(ref_grp_of.values()))
    ref_grp_code = np.array([ref_grps_u.index(ref_grp_of[n]) if n in ref_grp_of else -1
                             for n in ref_atoms_all.names], dtype=int)

    sel = {}   # other_res -> selection, key matrices and accumulators
    for other_res in other_resnames:
        ag = u.select_atoms(f"resname {other_res}")
        oth_names_u, oth_codes = np.unique(ag.names, return_inverse=True)
        n_keys = len(ref_names_u) * len(oth_names_u)

        oth_grp_of   = groups.get(other_res, {})
        oth_grps_u   = sorted(set(oth_grp_of.values()))
        oth_grp_code = np.array([oth_grps_u.index(oth_grp_of[n]) if n in oth_grp_of else -1
                                 for n in ag.names], dtype=int)
        n_grp_keys = len(ref_grps_u) * len(oth_grps_u)
        grp_key = np.where((ref_grp_code[:, None] >= 0) & (oth_grp_code[None, :] >= 0),
                           ref_grp_code[:, None] * len(oth_grps_u) + oth_grp_code[None, :],
                           -1)

        sel[other_res] = {
            'atoms':      ag,
            'names':      oth_names_u,
            'pair_key':   ref_codes[:, None] * len(oth_names_u) + oth_codes[None, :],
            'n':          np.zeros(n_keys, dtype=np.int64),
            'sum':        np.zeros(n_keys),
            'sumsq':      np.zeros(n_keys),
            'grps':       oth_grps_u,
            'grp_key':    grp_key,
            'grp_n':      np.zeros(n_grp_keys, dtype=np.int64),
            'grp_sum':    np.zeros(n_grp_keys),
            'grp_sumsq':  np.zeros(n_grp_keys),
        }

    angle_results   = [[] for _ in range(len(configs['angles']))]
    dihedral_results = [[] for _ in range(len(configs['dihedrals']))]

    if args.dt is not None:
        print(f"Frame time: frame index * {args.dt} ps (--dt)")
    else:
        print(f"Frame time: trajectory time (dt = {u.trajectory.dt} ps); "
              f"use --dt if -b/-e should be in real time")

    print(f"Processing {len(u.trajectory)} frames...")
    n_frames = 0

    for ts in u.trajectory:
        time = ts.frame * args.dt if args.dt is not None else ts.time
        if time < args.begin:
            continue
        if args.end >= 0 and time > args.end:
            break
        n_frames += 1

        # ── Per-atom and group distances ───────────────────────────────────────
        for other_res, s in sel.items():
            dist_mat = mda_distances.distance_array(
                ref_atoms_all.positions, s['atoms'].positions)

            if args.max_dist > 0:
                mask = dist_mat < args.max_dist
                d, keys, gkeys = dist_mat[mask], s['pair_key'][mask], s['grp_key'][mask]
            else:
                d, keys, gkeys = dist_mat.ravel(), s['pair_key'].ravel(), s['grp_key'].ravel()

            n_keys = len(s['n'])
            s['n']     += np.bincount(keys, minlength=n_keys)
            s['sum']   += np.bincount(keys, weights=d,     minlength=n_keys)
            s['sumsq'] += np.bincount(keys, weights=d * d, minlength=n_keys)

            # ── Group accumulation (mirrors C++ logic) ─────────────────────────
            in_grp = gkeys >= 0
            if in_grp.any():
                gk, gd = gkeys[in_grp], d[in_grp]
                n_grp_keys = len(s['grp_n'])
                s['grp_n']     += np.bincount(gk, minlength=n_grp_keys)
                s['grp_sum']   += np.bincount(gk, weights=gd,      minlength=n_grp_keys)
                s['grp_sumsq'] += np.bincount(gk, weights=gd * gd, minlength=n_grp_keys)

        # ── Angles ────────────────────────────────────────────────────────────
        for i, indices in enumerate(configs['angles']):
            pos = u.atoms[indices].positions  # positions update each frame
            ang = calc_angles(pos[0], pos[1], pos[2])
            angle_results[i].append(np.degrees(ang))

        # ── Dihedrals ─────────────────────────────────────────────────────────
        for i, indices in enumerate(configs['dihedrals']):
            pos = u.atoms[indices].positions  # positions update each frame
            dih = calc_dihedrals(pos[0], pos[1], pos[2], pos[3])
            dihedral_results[i].append(np.degrees(dih))
        

    print(f"Processed {n_frames} frames.")

    # ── Output 1: per-residue atom-atom distance statistics ───────────────────
    # One file per other_resname  (mirrors output_dist_stat in C++)
    res_files = {}
    for other_res in other_resnames:
        path = f"{prefix}_H-dist_{other_res}.dat"
        res_files[other_res] = open(path, 'w')
        res_files[other_res].write(
            f" Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n")
        print(f"Writing atom-atom distances to: {path}")

    # collect unique ref/other atom names per residue for stable integer ids
    ref_atom_ids  = {a.name: i for i, a in enumerate(ref_atoms_all)}
    other_atom_ids = {}
    for other_res in other_resnames:
        other_atom_ids[other_res] = {
            a.name: i for i, a in enumerate(u.select_atoms(f"resname {other_res}"))
        }

    def _mean_std(n, total, total_sq):
        """Mean and population standard deviation from count, sum and sum of squares."""
        avg = total / n
        return avg, np.sqrt(max(total_sq / n - avg * avg, 0.0))

    # keys are ordered ref-name-major, and both name lists come sorted from
    # np.unique, so iterating the keys in order gives rows sorted by (ref, other)
    for other_res, s in sel.items():
        n_oth = len(s['names'])
        for key in np.nonzero(s['n'])[0]:
            ref_aname = ref_names_u[key // n_oth]
            oth_aname = s['names'][key % n_oth]
            n = s['n'][key]
            avg, sd = _mean_std(n, s['sum'][key], s['sumsq'][key])
            id1  = ref_atom_ids.get(ref_aname, 0)
            id2  = other_atom_ids[other_res].get(oth_aname, 0)
            res_files[other_res].write(
                f"{id1:4d} {ref_aname:>8s} {id2:4d} {oth_aname:>8s} "
                f"{avg:8.3f} {sd:9.3f} {n:10d}\n")

    for f in res_files.values():
        f.close()

    # ── Output 2: per-residue GROUP distance statistics ───────────────────────
    # One file per other_resname  (mirrors output_grp_stat in C++)
    if groups and any(s['grp_n'].any() for s in sel.values()):
        for other_res, s in sel.items():
            path = f"{prefix}_H-dist_{other_res}_GROUP.dat"
            print(f"Writing group distances to: {path}")
            with open(path, 'w') as f:
                f.write(
                    f" Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n")

                n_oth = len(s['grps'])
                for key in np.nonzero(s['grp_n'])[0]:
                    ref_grp = ref_grps_u[key // n_oth]
                    oth_grp = s['grps'][key % n_oth]
                    n = s['grp_n'][key]
                    avg, sd = _mean_std(n, s['grp_sum'][key], s['grp_sumsq'][key])
                    # integer ids from types dict (order of appearance in groups file)
                    ni    = types.get(ref_name,  {}).get(ref_grp, 0)
                    nj    = types.get(other_res, {}).get(oth_grp, 0)
                    f.write(
                        f"{ni:4d} {ref_grp:>8s} {nj:4d} {oth_grp:>8s} "
                        f"{avg:8.3f} {sd:9.3f} {n:10d}\n")

   # ── Output 3: angle statistics ────────────────────────────────────────────
    if configs['angles']:
        # 1. Arquivo original de Médias
        path_avg = f"{prefix}_angles.dat"
        with open(path_avg, 'w') as f:
            f.write("Angle_ID    Average    StdDev\n")
            for i, data in enumerate(angle_results):
                if data:
                    f.write(f"{i:<11} {np.mean(data):.4f}    {np.std(data):.4f}\n")
            header = "Frame   " + "   ".join([f"Angle_{i}" for i in range(len(configs['angles']))])
            f.write(header + "\n")
            
            # Escreve os dados frame por frame
            for frame_idx in range(n_frames):
                row = f"{frame_idx + 1:<8}"
                for i in range(len(configs['angles'])):
                    if frame_idx < len(angle_results[i]):
                        row += f" {angle_results[i][frame_idx]:.4f}"
                f.write(row + "\n")
        print(f"Angles averages saved to: {path_avg}")

        path_hist = f"{prefix}_angles_hist.dat"
        write_histogram(path_hist, angle_results, "Angle", 0, 180, args.bin_width)
        print(f"Angles frequency of occurrence saved to: {path_hist}")


    # ── Output 4: dihedral statistics ─────────────────────────────────────────
    if configs['dihedrals']:
        # 1. Arquivo original de Médias
        path_avg = f"{prefix}_dihedrals.dat"
        with open(path_avg, 'w') as f:
            f.write("Dihedral_ID    Average    StdDev\n")
            for i, data in enumerate(dihedral_results):
                if data:
                    arr = np.array(data)
                    f.write(f"{i:<14} {np.mean(arr):.4f}    {np.std(arr):.4f}\n")
       
            header = "Frame   " + "   ".join([f"Dihedral_{i}" for i in range(len(configs['dihedrals']))])
            f.write(header + "\n")
            
            for frame_idx in range(n_frames):
                row = f"{frame_idx + 1:<8}"
                for i in range(len(configs['dihedrals'])):
                    if frame_idx < len(dihedral_results[i]):
                        row += f" {dihedral_results[i][frame_idx]:.4f}"
                f.write(row + "\n")
        print(f"Dihedrals frame-by-frame saved to: {path_avg}")

        path_hist = f"{prefix}_dihedrals_hist.dat"
        write_histogram(path_hist, dihedral_results, "Dihedral", -180, 180, args.bin_width)
        print(f"Dihedrals frequency of occurrence saved to: {path_hist}")
    print("Done.")


if __name__ == "__main__":
    main()