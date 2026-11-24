"""
Leitura da trajetória de saída do DICE (outname.xyz) para o ILLIAD.

O .xyz do DICE só traz símbolo do elemento + x y z, sem resname nem nome de
átomo. A topologia (nomes usados no grupos.txt e no -r) vem de uma de duas fontes:
  - o próprio .txt (ljname) do DICE + nmol + resnames (build_topology_from_txt),
    com nomes de átomo gerados ou copiados de um PDB de cada molécula; ou
  - um PDB/GRO do sistema inteiro com os átomos na MESMA ORDEM do .txt.
As coordenadas sempre vêm do .xyz.

Formato de cada configuração no outname.xyz:
           N
 Configuration number :    12000L =   49.8950  49.8950  49.8950
  C       13.23250      -14.49750       -1.86750
  ...
"""

# Comentários do código gerados com auxílio de IA (Claude, Anthropic) para facilitar a leitura por terceiros.


import re
import warnings
import numpy as np
import MDAnalysis as mda

_CONF_RE = re.compile(r"Configuration number\s*:\s*(\d+)")
_BOX_RE = re.compile(r"L\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)")


def scan_dice_xyz(path):
    
    steps, boxes, symbols = [], [], None
    with open(path) as f:
        frame = 0
        while True:
            line = f.readline()
            if not line.strip():
                break
            natoms = int(line)
            comment = f.readline()
            m = _CONF_RE.search(comment)
            steps.append(int(m.group(1)) if m else frame)
            b = _BOX_RE.search(comment)
            boxes.append([float(x) for x in b.groups()] if b else [np.nan] * 3)

            syms = [f.readline().split()[0] for _ in range(natoms)]
            if symbols is None:
                symbols = syms
            elif syms != symbols:
                bad = next(i for i, (a, b) in enumerate(zip(syms, symbols)) if a != b) \
                    if len(syms) == len(symbols) else None
                raise ValueError(
                    f"Configuração {steps[-1]} (frame {frame}) tem uma sequência de "
                    f"elementos diferente da primeira"
                    + (f" (1ª diferença no átomo {bad + 1})" if bad is not None else
                       f" ({len(syms)} átomos em vez de {len(symbols)})")
                    + ". O DICE pode ter trocado moléculas de tipos diferentes de lugar; "
                      "a topologia fixa do PDB não vale para esta trajetória.")
            frame += 1
    return np.array(steps), np.array(boxes), symbols


def _check_symbols(u, symbols):
    """Confere, átomo a átomo, se o símbolo do .xyz bate com o nome do átomo na topologia."""
    mismatches = []
    for a, s in zip(u.atoms, symbols):
        if s.upper() == "XX":          # átomo dummy do DICE (Z = 104)
            continue
        if not a.name.upper().startswith(s.upper()):
            mismatches.append(f"{a.id}: {a.resname} {a.name} x '{s}'")
    if mismatches:
        raise ValueError(
            f"{len(mismatches)} átomo(s) da topologia não batem com o elemento do .xyz "
            f"(a ordem dos átomos no PDB tem que ser a mesma do .txt do DICE). "
            f"Primeiros: {'; '.join(mismatches[:5])}")


_SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
            "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
            "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
            "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr").split()


def _symbol(z):
    if z == 104:                       # átomo dummy do DICE
        return "Xx"
    return _SYMBOLS[z] if 0 < z < len(_SYMBOLS) else "X"


def read_dice_txt(path):
    """
    Lê o .txt (ljname) do DICE e retorna, por tipo de molécula, a lista de
    números atômicos na ordem do arquivo. Linhas extras do formato de
    perturbação termodinâmica (depois dos m tipos) são ignoradas.
    """
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    n_types = int(lines[1].split()[0])
    types, i = [], 2
    for _ in range(n_types):
        na = int(lines[i].split()[0])      # aceita "11 (Furfural ...)"
        types.append([int(l.split()[1]) for l in lines[i + 1:i + 1 + na]])
        i += 1 + na
    return types


def _names_from_pdb(path, resname, expected_z):
    """
    Nomes dos átomos de uma molécula a partir de um PDB (mesma ordem do .txt).
    O PDB pode ter só essa molécula ou o sistema todo (ex.: um frame da MD);
    nesse caso usa o primeiro resíduo com esse resname.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(path)
    ag = u.atoms
    if len(ag) != len(expected_z):
        sel = u.select_atoms(f"resname {resname}")
        if len(sel) == 0:
            raise ValueError(f"{path} não tem resíduo {resname} e tem {len(ag)} átomos "
                             f"(a molécula no .txt tem {len(expected_z)}).")
        ag = sel.residues[0].atoms
    if len(ag) != len(expected_z):
        raise ValueError(f"{resname} em {path} tem {len(ag)} átomos, mas a molécula no "
                         f".txt tem {len(expected_z)}.")
    names = list(ag.names)
    bad = [f"{n} x Z={z}" for n, z in zip(names, expected_z)
           if z != 104 and not n.upper().startswith(_symbol(z).upper())]
    if bad:
        raise ValueError(f"Os nomes de {resname} em {path} não seguem a ordem do .txt "
                         f"(elemento não bate). Primeiros: {'; '.join(bad[:5])}")
    return names


def build_topology_from_txt(txt, nmol, resnames, names_pdb=None):
    """
    Monta a topologia (Universe sem coordenadas úteis) a partir do .txt do DICE.
      nmol      : número de moléculas de cada tipo, na ordem do .txt (como no input do DICE)
      resnames  : nome de resíduo de cada tipo (ex.: CMB SEH WAT)
      names_pdb : {resname: pdb de uma molécula} para reaproveitar os nomes de
                  átomo de um PDB existente (ex.: o mesmo usado na MD). Tipos sem
                  PDB recebem nomes gerados: elemento + posição na molécula (C1, H2...).
    """
    types = read_dice_txt(txt)
    if not (len(types) == len(nmol) == len(resnames)):
        raise ValueError(f"O .txt tem {len(types)} tipos de molécula, mas foram dados "
                         f"{len(nmol)} valores de nmol e {len(resnames)} resnames.")
    names_pdb = names_pdb or {}

    names, elements, res_of_atom, res_names = [], [], [], []
    for zs, n, rname in zip(types, nmol, resnames):
        if rname in names_pdb:
            mol_names = _names_from_pdb(names_pdb[rname], rname, zs)
        else:
            mol_names = [f"{_symbol(z)}{k + 1}" for k, z in enumerate(zs)]
        for _ in range(n):
            res_names.append(rname)
            names += mol_names
            elements += [_symbol(z) for z in zs]
            res_of_atom += [len(res_names) - 1] * len(zs)

    u = mda.Universe.empty(len(names), n_residues=len(res_names),
                           atom_resindex=np.array(res_of_atom), trajectory=True)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("types", elements)
    u.add_TopologyAttr("elements", elements)
    u.add_TopologyAttr("resnames", res_names)
    u.add_TopologyAttr("resids", np.arange(1, len(res_names) + 1))
    u.add_TopologyAttr("ids", np.arange(1, len(names) + 1))
    return u


def load_dice_universe(topology, xyz):
    """
    Cria o Universe com nomes/resíduos da topologia e coordenadas do .xyz do DICE.
    `topology` é um PDB/GRO (mesma ordem de átomos do .txt) ou um Universe já
    montado por build_topology_from_txt.
    Retorna (universe, steps, boxes); steps[i] é o número de configuração do frame i.
    """
    steps, boxes, symbols = scan_dice_xyz(xyz)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if isinstance(topology, mda.Universe):
            u = topology
            if len(u.atoms) != len(symbols):
                raise ValueError(f"O .txt × nmol dá {len(u.atoms)} átomos e o .xyz tem "
                                 f"{len(symbols)}. Confira os valores de nmol.")
            u.load_new(xyz, format="XYZ")
        else:
            u = mda.Universe(topology, xyz, format="XYZ")
    if len(symbols) != len(u.atoms):
        raise ValueError(f"A topologia tem {len(u.atoms)} átomos e o .xyz tem {len(symbols)}.")
    _check_symbols(u, symbols)
    if len(u.trajectory) != len(steps):
        raise ValueError(f"Lidos {len(u.trajectory)} frames pelo MDAnalysis, "
                         f"mas {len(steps)} configurações no arquivo.")
    return u, steps, boxes
