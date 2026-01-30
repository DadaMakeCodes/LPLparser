import numpy as np

ringvar = ['x', 'y', 'z']

#Ordinamento di Grado (DegLex)
Mdeglex = [[1.0, 1.0, 0.0, 0.0],
     [1.0, 0.0, 1.0, 0.0],
     [1.0, 0.0, 0.0, 0.1],] 

#Ordinamento Lessicografico (Lex) con ordine invertito
Mantilex = [[0.0, 0.0, 1.0],
     [0.0, 1.0, 0.0],
     [1.0, 0.0, 0.0],]

M = Mdeglex 


def lexGraph() -> list[list[float]]:
    n = len(ringvar)
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

