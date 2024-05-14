import numpy as np

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
