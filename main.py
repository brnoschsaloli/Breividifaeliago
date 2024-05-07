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