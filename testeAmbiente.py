import MDAnalysis as mda
import numpy as np

print("Ambiente configurado com sucesso!")
print(f"Versão do MDAnalysis: {mda.__version__}")
print(f"Versão do NumPy: {np.__version__}")

v1 = np.array([1.0, 2.0, 3.0])
v2 = np.array([4.0, 5.0, 6.0])
distancia = np.linalg.norm(v1 - v2)
print(f"Distância Euclidiana teste: {distancia:.2f}")