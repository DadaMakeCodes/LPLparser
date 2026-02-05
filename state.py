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

M = []


def lexGraph() -> list[list[float]]:
    n = len(ringvar)
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def Matrix_Prod(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    return (np.array(A) @ np.array(B)).tolist()

def Matrix_degLex() -> list[list[float]]:
    n = len(ringvar)
    M = []
    for i in range(n):
        row = [1.0] + [1.0 if i == j else 0.0 for j in range(n)]
        M.append(row)
    return M

def Matrix_AntiLex() -> list[list[float]]:
    n = len(ringvar)
    M = []
    for i in range(n):
        row = [1.0 if j == (n - 1 - i) else 0.0 for j in range(n)]
        M.append(row)
    return M

def Matrix_AntiDegLex() -> list[list[float]]:
    return Matrix_Prod( Matrix_AntiLex(), Matrix_degLex())
