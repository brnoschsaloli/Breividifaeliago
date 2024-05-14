"""
TransCal, APS7
Feita por Breno Schneider, David Conselvan, Rafael Paves, Thiago Victoriano (Grupo 11)
"""

import math
import numpy as np

def calular_matriz_senos(c, s): # recebe um angulo em radianos e devolve a matriz de rigidez
    return [
        [c**2, c*s, -c**2, -c*s],
        [c*s, s**2, -c*s, -s**2],
        [-c**2, -c*s, c**2, c*s],
        [-c*s, -s**2, c*s, s**2]
    ]

E = 200 * 10**9
A = 6 * 10**(-5)

n_pontos = 16
n_elementos = 29

graus_restritos = [0*2, 0*2 + 1, 15*2 + 1]

pontos = [
    (0, 0),
    (0.05, 0),
    (0.05, 0.05),
    (0.1, 0),
    (0.1, 0.05),
    (0.15, 0),
    (0.15, 0.05),
    (0.2, 0),
    (0.2, 0.05),
    (0.25, 0),
    (0.25, 0.05),
    (0.3, 0),
    (0.3, 0.05),
    (0.35, 0),
    (0.35, 0.05),
    (0.4, 0)
]

elementos = [
    {"pontos": [0, 1]}, # elementos horizontais baixo
    {"pontos": [1, 3]},
    {"pontos": [3, 5]},
    {"pontos": [5, 7]},
    {"pontos": [7, 9]},
    {"pontos": [9, 11]},
    {"pontos": [11, 13]},
    {"pontos": [13, 15]},
    {"pontos": [2, 4]}, # elementos horizontais cima
    {"pontos": [4, 6]},
    {"pontos": [6, 8]},
    {"pontos": [8, 10]},
    {"pontos": [10, 12]},
    {"pontos": [12, 14]},
    {"pontos": [1, 2]}, # elementos verticais
    {"pontos": [3, 4]},
    {"pontos": [5, 6]},
    {"pontos": [7, 8]},
    {"pontos": [9, 10]},
    {"pontos": [11, 12]},
    {"pontos": [13, 14]},
    {"pontos": [0, 2]}, # elementos diagonais
    {"pontos": [2, 3]},
    {"pontos": [4, 5]},
    {"pontos": [6, 7]},
    {"pontos": [7, 10]},
    {"pontos": [9, 12]},
    {"pontos": [11, 14]},
    {"pontos": [14, 15]},
]

K = np.zeros((n_pontos*2, n_pontos*2))
k_g = np.zeros((n_pontos*2 - len(graus_restritos), n_pontos*2 - len(graus_restritos)))

for elemento in elementos:
    p1 = pontos[elemento["pontos"][0]]
    p2 = pontos[elemento["pontos"][1]]
    elemento["l"] = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
    elemento["s"] = (p2[1] - p1[1])/elemento["l"]
    elemento["c"] = (p2[0] - p1[0])/elemento["l"]
    elemento["k"] = calular_matriz_senos(elemento["c"], elemento["s"])
    elemento["graus"] = {0: elemento["pontos"][0]*2, 1: elemento["pontos"][0]*2 + 1, 2: elemento["pontos"][1]*2, 3: elemento["pontos"][1]*2 + 1}
    for i in range(4):
        for j in range(4):
            K[elemento["graus"][i]][elemento["graus"][j]] += (E*A/elemento["l"]) * elemento["k"][i][j]
            if (elemento["graus"][i] not in graus_restritos and elemento["graus"][j] not in graus_restritos):
                k_g[elemento["graus"][i] - 2][elemento["graus"][j] - 2] += (E*A/elemento["l"]) * elemento["k"][i][j]

P = [
    0, # 1
    0,
    0, # 2
    0,
    0, # 3
    0,
    0, # 4
    0,
    0, # 5
    0,
    0, # 6
    -100,
    0, # 7
    0,
    0, # 8
    -100,
    0, # 9
    0,
    0, # 10
    -100,
    0, # 11
    0,
    0, # 12
    0,
    0, # 13
    0,
    0, # 14
    0,
    0, # 15
]

u = np.linalg.solve(k_g, P)

U = np.zeros((n_pontos*2, ))
indice = 0
for i in range(n_pontos*2):
    if i not in graus_restritos:
        U[i] = u[i - indice]
    else:
        indice += 1

for i in range(len(U)):
    if i % 2 == 0:
        print(f"Ponto {pontos[i//2]}, eixo x --> u{i} = {U[i]:.8f} m")
    else:
        print(f"Ponto {pontos[i//2]}, eixo y --> u{i} = {U[i]:.8f} m")

print("-------------------------")
for elemento in elementos:
    matriz_deformacao = [-elemento["c"], -elemento["s"], elemento["c"], elemento["s"]]
    tensao = 0
    for i in range(4):
        tensao += matriz_deformacao[i] * elemento["graus"][i]
    tensao = (E/elemento["l"]) * tensao
    print(f"Elemento com nós em {pontos[elemento['pontos'][0]]} e {pontos[elemento['pontos'][1]]} --> tensão = {tensao} Pa")
print("-------------------------")

R1x = 0 # Nó mais da esquerda
R1y = 0 # Nó mais da esquerda
R2y = 0 # Nó mais da direita
for i in range(n_pontos*2):
    R1x += K[0][i] * U[i]
    R1y += K[1][i] * U[i]
    R2y += K[31][i] * U[i]

print(f"R1x = {R1x} N")
print(f"R1y = {R1y} N")
print(f"R2y = {R2y} N")

def gauss_seidel(A, b, x0=None, max_iterations=1000, tolerance=1e-10):
    """
    Resolve o sistema linear Ax = b usando o método de Gauss-Seidel.

    :param A: Matriz de coeficientes.
    :param b: Vetor de termos independentes.
    :param x0: Palpite inicial para a solução.
    :param max_iterations: Número máximo de iterações.
    :param tolerance: Tolerância para a convergência.
    :return: Vetor de solução x.
    """
    n = len(b)
    x = x0 if x0 is not None else np.zeros_like(b)
    for _ in range(max_iterations):
        x_new = np.copy(x)
        for i in range(n):
            s1 = np.dot(A[i, :i], x_new[:i])
            s2 = np.dot(A[i, i+1:], x[i+1:])
            x_new[i] = (b[i] - s1 - s2) / A[i, i]
        if np.linalg.norm(x_new - x, ord=np.inf) < tolerance:
            return x_new
        x = x_new
    raise ValueError("A solução não convergiu após o máximo de iterações.")


def verificar_falhas(tensoes, propriedades, elementos, compressive_loads):
    falhas_tensao = []
    falhas_flambagem = []
    for i, tensao in enumerate(tensoes):
        area = propriedades['area'][i]
        comprimento = ...  # Calcule o comprimento do elemento i
        carga_critica = np.pi**2 * propriedades['E'][i] * (area / (comprimento**2))
        
        # Verificar falha por tensão
        if abs(tensao) > propriedades['tensao_tracao'][i]:
            falhas_tensao.append((i, tensao, "Tração"))
        if tensao < -propriedades['tensao_compressao'][i]:
            falhas_tensao.append((i, tensao, "Compressão"))
        
        # Verificar falha por flambagem (apenas para cargas compressivas)
        if i in compressive_loads and -compressive_loads[i] > carga_critica:
            falhas_flambagem.append((i, compressive_loads[i], carga_critica))
            
    return falhas_tensao, falhas_flambagem
