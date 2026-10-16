"""
Stores the results of an illiad.py run in a SQLite database.

Each call to save_run adds one row to `runs` (the run parameters) and the
results of that run to the other tables, all linked by run_id. Running
illiad.py several times with the same --db accumulates runs in one file,
so they can be compared with SQL.
"""
import sqlite3
from datetime import datetime

import numpy as np

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id          INTEGER PRIMARY KEY,
    created_at  TEXT NOT NULL,
    pdb         TEXT NOT NULL,
    reference   TEXT NOT NULL,
    groups_file TEXT,
    output      TEXT,
    begin_ps    REAL,
    end_ps      REAL,       -- -1 = no limit
    max_dist    REAL,       -- -1 = no cutoff
    dt          REAL,       -- NULL = trajectory time was used
    bin_width   REAL,
    n_frames    INTEGER
);

CREATE TABLE IF NOT EXISTS atom_distances (
    run_id    INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    resname   TEXT,
    id1       INTEGER,
    label1    TEXT,
    id2       INTEGER,
    label2    TEXT,
    dist      REAL,
    sigma     REAL,
    frequency INTEGER
);

CREATE TABLE IF NOT EXISTS group_distances (
    run_id    INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    resname   TEXT,
    id1       INTEGER,
    label1    TEXT,
    id2       INTEGER,
    label2    TEXT,
    dist      REAL,
    sigma     REAL,
    frequency INTEGER
);

CREATE TABLE IF NOT EXISTS geometry_stats (
    run_id   INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    kind     TEXT CHECK (kind IN ('angle', 'dihedral')),
    geom_id  INTEGER,
    atoms    TEXT,          -- e.g. "CMB O0R (12) - SEH C08 (197) - SEH H0P (215)"
    average  REAL,
    stddev   REAL
);

CREATE TABLE IF NOT EXISTS geometry_series (
    run_id   INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    kind     TEXT CHECK (kind IN ('angle', 'dihedral')),
    geom_id  INTEGER,
    frame    INTEGER,       -- 1-based, counting only processed frames (as in the .dat)
    time_ps  REAL,
    value    REAL
);

CREATE TABLE IF NOT EXISTS histograms (
    run_id     INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    kind       TEXT CHECK (kind IN ('angle', 'dihedral')),
    geom_id    INTEGER,
    bin_start  REAL,
    bin_end    REAL,
    bin_center REAL,
    count      INTEGER,
    freq       REAL         -- percent of the processed frames
);

CREATE INDEX IF NOT EXISTS idx_atom_distances_run  ON atom_distances(run_id);
CREATE INDEX IF NOT EXISTS idx_group_distances_run ON group_distances(run_id);
CREATE INDEX IF NOT EXISTS idx_geometry_stats_run  ON geometry_stats(run_id);
CREATE INDEX IF NOT EXISTS idx_geometry_series_run ON geometry_series(run_id, kind, geom_id);
CREATE INDEX IF NOT EXISTS idx_histograms_run      ON histograms(run_id, kind, geom_id);
"""


def _distance_rows(run_id, rows):
    # numpy integers are not accepted by sqlite3, so every value is converted
    return [(run_id, res, int(id1), str(l1), int(id2), str(l2), float(d), float(sd), int(n))
            for res, id1, l1, id2, l2, d, sd, n in rows]


def save_run(path, args, n_frames, frame_times, atom_rows, group_rows, geometries):
    """
    Saves one run and returns its id.

    atom_rows / group_rows : (resname, id1, label1, id2, label2, dist, sigma, frequency)
    geometries : one dict per kind ('angle' / 'dihedral') with
        'atoms'  : description of the atoms of each angle/dihedral
        'series' : per-angle list of values, one per processed frame
        'edges', 'hist' : histogram bin edges and per-angle (counts, freqs)
    """
    con = sqlite3.connect(path)
    try:
        con.execute("PRAGMA foreign_keys = ON")
        con.executescript(SCHEMA)
        with con:   # one transaction: a run is saved completely or not at all
            cur = con.execute(
                "INSERT INTO runs (created_at, pdb, reference, groups_file, output, begin_ps,"
                " end_ps, max_dist, dt, bin_width, n_frames)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (datetime.now().isoformat(timespec='seconds'), args.pdb, args.reference,
                 args.groups, args.output, args.begin, args.end, args.max_dist, args.dt,
                 args.bin_width, n_frames))
            run_id = cur.lastrowid

            con.executemany("INSERT INTO atom_distances VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                            _distance_rows(run_id, atom_rows))
            con.executemany("INSERT INTO group_distances VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                            _distance_rows(run_id, group_rows))

            for g in geometries:
                kind = g['kind']
                for i, data in enumerate(g['series']):
                    if data:
                        con.execute("INSERT INTO geometry_stats VALUES (?, ?, ?, ?, ?, ?)",
                                    (run_id, kind, i, g['atoms'][i],
                                     float(np.mean(data)), float(np.std(data))))
                    con.executemany(
                        "INSERT INTO geometry_series VALUES (?, ?, ?, ?, ?, ?)",
                        [(run_id, kind, i, f + 1, float(t), float(v))
                         for f, (t, v) in enumerate(zip(frame_times, data))])

                edges = g['edges']
                for i, (counts, freqs) in enumerate(g['hist']):
                    con.executemany(
                        "INSERT INTO histograms VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        [(run_id, kind, i, float(edges[b]), float(edges[b + 1]),
                          float((edges[b] + edges[b + 1]) / 2), int(counts[b]), float(freqs[b]))
                         for b in range(len(edges) - 1)])
    finally:
        con.close()
    return run_id


# ─────────────────────────────────────────────────────────────────────────────
# Guia das tabelas
# ─────────────────────────────────────────────────────────────────────────────
# O banco tem 6 tabelas. `runs` é a central e descreve cada execução; as outras
# 5 guardam os resultados e apontam para a execução pela coluna run_id.
#
#                       runs (1 linha por execução)
#                         │ id
#      ┌────────────┬─────┴──────┬───────────────┬────────────────┐
#  atom_distances group_distances geometry_stats geometry_series histograms
#
# 1. runs: as execuções
#    Cada vez que o illiad.py roda com --db entra uma linha nova. Ela registra
#    COMO os resultados foram obtidos, para saber o que está sendo comparado.
#      id                        número da execução (é o run_id das outras tabelas)
#      created_at                data e hora da execução
#      pdb, reference,           arquivo de entrada, resíduo de referência (-r),
#      groups_file, output       arquivo de grupos (-g) e prefixo (-o)
#      begin_ps, end_ps          intervalo de tempo (-b, -e); -1 = até o fim
#      max_dist                  corte de distância (-d); -1 = sem corte
#      dt                        --dt; NULL = tempo da trajetória
#      bin_width                 largura dos bins dos histogramas
#      n_frames                  quantos frames foram processados
#
# 2. atom_distances: distâncias átomo a átomo   (= *_H-dist_<RES>.dat)
#    Cada linha é um par de átomos (um da referência, outro do outro resíduo),
#    resumido sobre todos os frames.
#      resname                   resíduo do segundo átomo
#      id1, label1               índice e nome do átomo da referência
#      id2, label2               índice e nome do átomo do outro resíduo
#      dist, sigma               média e desvio padrão da distância (Å)
#      frequency                 número de frames em que o par ficou abaixo do corte
#    Ex.: H0Q – H3B, dist 2.63, frequency 3752 → a menos de 3 Å em 37,5% dos
#    10001 frames. Com -d, dist e sigma consideram só os frames em contato, e
#    pares que nunca chegaram ao corte não aparecem.
#
# 3. group_distances: distâncias entre grupos   (= *_H-dist_<RES>_GROUP.dat)
#    Mesmas colunas, mas cada linha é um par de grupos do arquivo de grupos:
#    as distâncias de todos os átomos de um grupo contra todos os do outro.
#    Aqui frequency NÃO é número de frames: é o número de pares de átomos em
#    contato somado sobre os frames, e pode passar do total de frames.
#    id1/id2 são os números dos grupos na ordem em que aparecem no arquivo de
#    grupos (ex.: o grupo `a` é o 7º do SEH, por isso id1 = 7).
#
# 4. geometry_stats: resumo dos ângulos e diedros   (1º bloco de *_angles.dat
#    e *_dihedrals.dat). Uma linha por ângulo/diedro.
#      kind                      'angle' ou 'dihedral'
#      geom_id                   número do ângulo/diedro (0, 1…), na ordem do arquivo
#      atoms                     átomos que o definem, com o id de cada um
#      average, stddev           média e desvio padrão (graus)
#    A coluna atoms só existe no banco: evita voltar ao arquivo de grupos para
#    lembrar o que é o "diedro 1".
#
# 5. geometry_series: valor frame a frame   (2º bloco de *_angles.dat e
#    *_dihedrals.dat). Uma linha por ângulo × frame; é a maior tabela.
#      kind, geom_id             qual ângulo/diedro
#      frame                     número do frame processado (começa em 1)
#      time_ps                   tempo do frame (ps)
#      value                     valor naquele frame (graus)
#    Serve para perguntas de "quando": em que momentos o ângulo esteve numa
#    faixa, como evolui no tempo, se há transições.
#
# 6. histograms: frequência de ocorrência   (= *_hist.dat)
#    Uma linha por ângulo × intervalo (bin).
#      kind, geom_id             qual ângulo/diedro
#      bin_start, bin_end,       limites e centro do intervalo (graus)
#      bin_center
#      count                     quantos frames caíram no intervalo
#      freq                      count em % dos frames processados
#    A distribuição fica pronta para plotar; também dá para recalculá-la a
#    partir de geometry_series com outra largura de bin.
#
# Qual tabela usar:
#   Com que parâmetros rodei isso?               → runs
#   Quais átomos fazem contato, e quanto?        → atom_distances
#   Quais regiões das moléculas interagem?       → group_distances
#   Qual o valor típico de um ângulo/diedro?     → geometry_stats
#   Quando o ângulo esteve em tal faixa?         → geometry_series
#   Como os valores se distribuem?               → histograms