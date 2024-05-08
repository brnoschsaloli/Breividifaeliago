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

from numpy.linalg import solve
import numpy as np

def dados_entrada():
    # Exemplo de dados de entrada. Estes devem ser adaptados para seu caso específico.
    nos = np.array([
        [0, 0],       # no 1
        [0.05, 0],    # no 2  
        [0.1, 0],     # no 3
        [0.15, 0],    # no 4
        [0.2, 0],     # no 5
        [0.25, 0],    # no 6
        [0.3, 0],     # no 7
        [0.35, 0],    # no 8
        [0.4, 0],     # no 9
        [0.05, 0.05], # no 10
        [0.1, 0.05],  # no 11
        [0.15, 0.05], # no 12
        [0.2, 0.05],  # no 13
        [0.25, 0.05], # no 14
        [0.3, 0.05],  # no 15
        [0.35,0.05]   # no 16
    ])
    
    incidencia = np.array([  #conectividade dos elementos (nos)
        [1, 2],         #incidencia
        [1, 10],
        [2, 3],
        [2, 10],
        [3, 10],
        [3, 11],
        [3, 4],
        [4, 11],
        [4, 12],
        [4, 5],
        [5, 12],
        [5, 13],
        [5, 14],
        [5, 6],
        [6, 14],
        [6, 15],
        [6, 7],
        [7, 15],
        [7, 16],
        [7, 8],
        [8, 16],
        [8, 9],
        [9, 16],
        [10, 11],
        [11, 12],
        [12, 13],
        [13, 14],
        [14, 15],
        [15, 16]
    ])
    
    propriedades = {
        'area': [6e-2] * 5,  # Área da seção transversal em m² para todos os incidencia
        'E': [200e9] * 5,    # Módulo de elasticidade em Pascal
        'tensao_tracao': [400e6] * 5,
        'tensao_compressao': [300e6] * 5
    }
    
    gdl = {
        'restricoes': [0, 1, 17],  # Graus de liberdade restritos
        'cargas': {
            12: [0, -100],  
            13: [0, -100],
            14: [0, -100]  
        }
    }
    
    return nos, incidencia, propriedades, gdl

def montar_matriz_global(nos, elementos, propriedades):
    n_gdl = nos.shape[0] * 2  # Cada nó tem 2 graus de liberdade (x, y)
    K_global = np.zeros((n_gdl, n_gdl))

    for elemento in elementos:
        n1, n2 = elemento - 1  # Convertendo para índice base-0
        x1, y1 = nos[n1]
        x2, y2 = nos[n2]
        L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        cos = (x2 - x1) / L
        sin = (y2 - y1) / L
        k = propriedades['E'][0] * propriedades['area'][0] / L
        
        # Matriz de rigidez para um elemento em coordenadas locais
        k_matrix = k * np.array([
            [cos**2, cos*sin, -cos**2, -cos*sin],
            [cos*sin, sin**2, -cos*sin, -sin**2],
            [-cos**2, -cos*sin, cos**2, cos*sin],
            [-cos*sin, -sin**2, cos*sin, sin**2]
        ])
        
        # Adicionar à matriz global
        indices = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for i in range(4):
            for j in range(4):
                K_global[indices[i], indices[j]] += k_matrix[i, j]
    
    return K_global


def aplicar_condicoes_contorno(K, gdl):
    # Aplicar condições de contorno
    restricoes = [0, 1]  # Restrições em x e y para o nó 1
    restricoes.append(17)  # Restrição em y para o nó 9 (índice 17 se considerarmos 0-indexed para 2 gdl por nó)

    for i in restricoes:
        K[i, :] = 0
        K[:, i] = 0
        K[i, i] = 1
    
    return K


def montar_vetor_forca(gdl, n_gdl):
    F = np.zeros(n_gdl)
    for indice, carga in gdl['cargas'].items():
        F[indice] = carga[0]
        F[indice + 1] = carga[1]
    
    return F

def main():
    nos, incidencia, propriedades, gdl = dados_entrada()
    n_gdl = nos.shape[0] * 2
    K_global = montar_matriz_global(nos, incidencia, propriedades)
    K_global = aplicar_condicoes_contorno(K_global, gdl)
    F = montar_vetor_forca(gdl, n_gdl)
    
    # Solução do sistema
    deslocamentos = np.linalg.solve(K_global, F)
    print("Deslocamentos nos nós:", deslocamentos)

if __name__ == '__main__':
    main()
