import os
import MDAnalysis as mda

def ler_arquivo_grupos(caminho_arquivo):
  
    dados = {
        'angles': [],
        'dihedrals': [],
        'groups': {}
    }
    
    modo_atual = None
    grupo_atual_nome = None

    print(f"Lendo arquivo: {caminho_arquivo}")

    try:
        with open(caminho_arquivo, 'r') as f:
            for linha in f:
                linha = linha.strip()
                
               
                if not linha: 
                    continue
                
               
                if linha.startswith('['):
                    section = linha.replace('[', '').replace(']', '').strip()
                    
                    if section == "angle":
                        modo_atual = "angle"
                    elif section == "dihedral":
                        modo_atual = "dihedral"
                    elif section == "group":
                        modo_atual = "group"
                        grupo_atual_nome = None 
                    else:
                        modo_atual = "unknown"
                    continue

              
                if modo_atual == "angle":
                    
                    atomos = linha.split()
                    if len(atomos) >= 3:
                        dados['angles'].append(atomos[:3])
                
                elif modo_atual == "dihedral":
                    
                    atomos = linha.split()
                    if len(atomos) >= 4:
                        dados['dihedrals'].append(atomos[:4])
                
                elif modo_atual == "group":
                    
                    if grupo_atual_nome is None:
                        grupo_atual_nome = linha
                        dados['groups'][grupo_atual_nome] = []
                    else:
                       
                        atomos = linha.split()
                        dados['groups'][grupo_atual_nome].extend(atomos)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em {caminho_arquivo}")
        return None
    except Exception as e:
        print(f"Erro ao ler o arquivo: {e}")
        return None
    return dados


if __name__ == "__main__":
    caminho_teste = os.path.join("teste_dados", "grupos.txt")
 
    if not os.path.exists(caminho_teste):
        caminho_teste = os.path.join("..", "teste_dados", "grupos.txt")

    resultado = ler_arquivo_grupos(caminho_teste)
    
    if resultado:
        print(f"Angulos encontrados: {len(resultado['angles'])}")
        print(f"Exemplo de angulo: {resultado['angles'][1]}")
        print(f"Grupos encontrados: {list(resultado['groups'].keys())}")
        print(f"Atomos no grupo 'Agua': {resultado['groups']['Agua']}")