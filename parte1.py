'''
DADOS DE ENTRADA:

    - Coordenada dos nos;
    - Conectividade dos elementos (Incidencia);
    - Propriedades dos elementos;
    - Numeracao dos graus de liberdade (gdl);
    - Grau de liberdade com restricao;
    - Grau de liberdade com cargas aplicadas e intensidade no sentido da carga.

PRE-PROCESSAMENTO:
    - Montagem das matrizes dos elementos;
    - Superposicao das matrizes - Rigidez global da estrutura;
    - Montagem do vetor de carga global da estrutura;
    - Aplicacao das condicoes de contorno.

SOLUCAO:
    - Solucao do sistema de equacoes.

POS-PROCESSAMENTO:
    - Deslocamento nos nos;
    - Tensao em cada elemento;
    - Reacoes de apoio nos nos com restricao.

'''

import math
import numpy as np

def calular_matriz_senos(c, s): # recebe um angulo em radianos e devolve a matriz de rigidez
    return [
        [c**2, c*s, -c**2, -c*s],
        [c*s, s**2, -c*s, -s**2],
        [-c**2, -c*s, c**2, c*s],
        [-c*s, -s**2, c*s, s**2]
    ]

# print(calular_matriz_senos(np.pi/2)[0][0], np.cos(np.pi/2), np.cos(np.pi/2)**2)

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

# print(K[31])

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

# print(u)
print(U)

# print(k_g[0])
# print(K[31])

for elemento in elementos:
    matriz_deformacao = [-elemento["c"], -elemento["s"], elemento["c"], elemento["s"]]
    tensao = 0
    for i in range(4):
        tensao += matriz_deformacao[i] * elemento["graus"][i]
    tensao = (E/elemento["l"]) * tensao

R1x = 0
R1y = 0
R2y = 0
for i in range(n_pontos*2):
    R1x += K[0][i] * U[i]
    R1y += K[1][i] * U[i]
    R2y += K[31][i] * U[i]

print(f"R1x = {R1x} N")
print(f"R1y = {R1y} N")
print(f"R2y = {R2y} N")
