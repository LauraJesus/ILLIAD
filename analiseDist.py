import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import distances
import os
//NAO TERMINADA
def analisar_distancias(caminho_pdb, lista_pares_atomicos):
    nomeArquivo = "traj_w_nj_m_c_rtt.pdb"
    nomePasta = "teste_dados"

    caminho_pdb = os.path.join(nomePasta, nomeArquivo)
    if not os.path.exists(caminho_pdb):
        caminho_pdb = os.path.join("..", nomePasta, nomeArquivo)
    
    print(f"----Análise de Distâncias para o arquivo: {caminho_pdb}----")

    try:
        u = mda.Universe(caminho_pdb)
        print(f"PDB '{nomeArquivo}' carregado.")
    except FileNotFoundError:
        print(f"Erro: Arquivo PDB não encontrado em {caminho_pdb}")
    except Exception as e:
        print(f"Erro ao carregar o PDB: {e}")
        return
    
    grupo_ref = u.select_atoms(f"resname {res_referencia}")

    grupo_outros = u.select_atoms(f"not resname {res_referencia}")


