import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_dihedrals
from MDAnalysis.analysis import distances as mda_distances
import argparse
from collections import defaultdict

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

# ─────────────────────────────────────────────
# Parse groups file
# ─────────────────────────────────────────────
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
# Main
# ─────────────────────────────────────────────
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
    args = parser.parse_args()

    # ── Load universe ──────────────────────────────────────────────────────────
    u      = mda.Universe(args.pdb)
    prefix = args.output if args.output else args.pdb.rsplit('.', 1)[0]

    # ── Parse groups ───────────────────────────────────────────────────────────
    configs = parse_groups_file(args.groups, u)
    groups  = configs['groups']   # resname -> {atomname -> group_label}
    types   = configs['types']    # resname -> {group_label -> int_id}

    ref_name     = args.reference
    other_resnames = [r for r in np.unique(u.atoms.resnames) if r != ref_name]

    # ── Storage ────────────────────────────────────────────────────────────────
    # distances[res_idx][ref_atom_idx][other_atom_idx] -> list of floats
    # We'll store per-atom distances as (ref_atom_name, other_resname, other_atom_name) -> [dist]
    atom_dist_data   = defaultdict(list)   # (ref_atom_name, other_res, other_atom_name) -> [dist]

    # group_stats: ref_grp -> other_res -> other_grp -> {'sum': , 'n': , 'data': []}
    group_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(
        lambda: {'sum': 0.0, 'n': 0, 'data': []}
    )))

    angle_results   = [[] for _ in range(len(configs['angles']))]
    dihedral_results = [[] for _ in range(len(configs['dihedrals']))]

    ref_atoms_all = u.select_atoms(f"resname {ref_name}")

    print(f"Processing {len(u.trajectory)} frames...")
    n_frames = 0

    for ts in u.trajectory:
        if ts.time < args.begin:
            continue
        if args.end >= 0 and ts.time > args.end:
            break
        n_frames += 1

        # ── Per-atom distances ─────────────────────────────────────────────────
        for other_res in other_resnames:
            other_atoms_all = u.select_atoms(f"resname {other_res}")

            dist_mat = mda_distances.distance_array(
                ref_atoms_all.positions, other_atoms_all.positions)

            for ri, ref_at in enumerate(ref_atoms_all):
                for oi, oth_at in enumerate(other_atoms_all):
                    d = dist_mat[ri, oi]

                    if args.max_dist > 0 and d >= args.max_dist:
                        continue

                    key = (ref_at.name, other_res, oth_at.name)
                    atom_dist_data[key].append(d)

                    # ── Group accumulation (mirrors C++ logic) ─────────────────
                    if groups:
                        ref_grp = groups.get(ref_name, {}).get(ref_at.name)
                        oth_grp = groups.get(other_res, {}).get(oth_at.name)
                        if ref_grp and oth_grp:
                            entry = group_stats[ref_grp][other_res][oth_grp]
                            entry['sum']  += d
                            entry['n']    += 1
                            entry['data'].append(d)

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

    for (ref_aname, other_res, oth_aname), data in sorted(atom_dist_data.items()):
        arr  = np.array(data)
        n    = len(arr)
        avg  = arr.mean()
        sd   = arr.std()
        id1  = ref_atom_ids.get(ref_aname, 0)
        id2  = other_atom_ids[other_res].get(oth_aname, 0)
        res_files[other_res].write(
            f"{id1:4d} {ref_aname:>8s} {id2:4d} {oth_aname:>8s} "
            f"{avg:8.3f} {sd:9.3f} {n:10d}\n")

    for f in res_files.values():
        f.close()

    # ── Output 2: per-residue GROUP distance statistics ───────────────────────
    # One file per other_resname  (mirrors output_grp_stat in C++)
    if groups and group_stats:
        for other_res in other_resnames:
            path = f"{prefix}_H-dist_{other_res}_GROUP.dat"
            print(f"Writing group distances to: {path}")
            with open(path, 'w') as f:
                f.write(
                    f" Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n")

                for ref_grp in sorted(group_stats.keys()):
                    if other_res not in group_stats[ref_grp]:
                        continue
                    for oth_grp in sorted(group_stats[ref_grp][other_res].keys()):
                        entry = group_stats[ref_grp][other_res][oth_grp]
                        n     = entry['n']
                        if n == 0:
                            continue
                        avg   = entry['sum'] / n
                        arr   = np.array(entry['data'])
                        sd    = arr.std()
                        # integer ids from types dict (order of appearance in groups file)
                        ni    = types.get(ref_name,  {}).get(ref_grp, 0)
                        nj    = types.get(other_res, {}).get(oth_grp, 0)
                        f.write(
                            f"{ni:4d} {ref_grp:>8s} {nj:4d} {oth_grp:>8s} "
                            f"{avg:8.3f} {sd:9.3f} {n:10d}\n")

    # ── Output 3: angle statistics ────────────────────────────────────────────
    if configs['angles']:
        path = f"{prefix}_angles_AVG.dat"
        with open(path, 'w') as f:
            f.write("Angle_ID    Average    StdDev\n")
            for i, data in enumerate(angle_results):
                if data:
                    f.write(f"{i:<11} {np.mean(data):.4f}    {np.std(data):.4f}\n")
        print(f"Angles saved to: {path}")

    # ── Output 4: dihedral statistics ─────────────────────────────────────────
    if configs['dihedrals']:
        path = f"{prefix}_dihedrals_AVG.dat"
        with open(path, 'w') as f:
            f.write("Dihedral_ID    Average    StdDev\n")
            for i, data in enumerate(dihedral_results):
                if data:
            
                    arr = np.array(data)
                    f.write(f"{i:<14} {np.mean(arr):.4f}    {np.std(arr):.4f}\n")
        print(f"Dihedrals saved to: {path}")

    print("Done.")


if __name__ == "__main__":
    main()