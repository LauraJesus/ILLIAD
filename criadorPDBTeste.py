import MDAnalysis as mda
import os

nomeArquivo = "traj_w_nj_m_c_rtt.pdb"
nomePasta = "teste_dados"
caminho_pdb = os.path.join(nomePasta, nomeArquivo)

if not os.path.exists(caminho_pdb):
    caminho_pdb = os.path.join("..", nomePasta, nomeArquivo)

try:
    u = mda.Universe(caminho_pdb)
    print(f"PDB '{nomeArquivo}' carregado.")
    print(f"Número de átomos: {len(u.atoms)}")
    print(f"Número de resíduos: {len(u.residues)}")
except FileNotFoundError:
    print(f"Erro: Arquivo PDB não encontrado em {caminho_pdb}")
except Exception as e:
    print(f"Erro ao carregar o PDB: {e}")