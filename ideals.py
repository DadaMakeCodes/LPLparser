import Ordinamento
from polynomial import Polinomio
from typing import List
from monomial import Monomio
import state

def generateIdeal() -> List[Polinomio]:
    I = []
    m = 1
    while m != '\\':
        m = str(input())
        if m != '\\':
            I.append(Ordinamento.normal_pol_form(Polinomio.parse(m)))
    return I

def Isum(I: List[Polinomio], J: List[Polinomio]) -> List[Polinomio]:
    return Ordinamento.reduced_grobner_basis(I + J)

def Imul(I: List[Polinomio], J: List[Polinomio]) -> List[Polinomio]:
    prods = []
    for f in I:
        for g in J:
            prods.append(Ordinamento.normal_pol_form(f * g))
    return prods

def isEqual(I: List[Polinomio], J: List[Polinomio]) -> bool:
    GbI = Ordinamento.reduced_grobner_basis(I)
    GbJ = Ordinamento.reduced_grobner_basis(J)
    return GbI == GbJ

def is_pol_in_ideal(f: Polinomio, I: List[Polinomio]) -> bool:
    reduced_f = Ordinamento.reduce_by_set(f, I)
    return reduced_f.terms == []

def isIinJ(I: List[Polinomio], J: List[Polinomio]) -> bool:
    for f in I:
        if not is_pol_in_ideal(f, J):
            return False
    return True

def isEqual2(I: List[Polinomio], J: List[Polinomio]) -> bool:
    return isIinJ(I, J) and isIinJ(J, I)

print('Inserire variabili:')
state.ringvar = str(input()).split(',')
#print('Inserire numero variabili:')
#state.ringvar = Ordinamento.numerated_variables(str(input()), int(input()))
print('Variabili:', state.ringvar)
#print('M:', state.M)

I = generateIdeal()
J = Ordinamento.reduced_grobner_basis(I)
print(f'Grobner basis of I: ({", ".join(str(p) for p in J)})')
print('Inserire polinomio da testare per appartenenza all\'ideale I:')
f = Polinomio.parse(str(input()))
if is_pol_in_ideal(f, J):
    print(f'Il polinomio {f} appartiene all\'ideale I.')
else:
    print(f'Il polinomio {f} non appartiene all\'ideale I.')    
# print(normal_pol_form(p1))
