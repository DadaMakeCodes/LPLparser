from __future__ import annotations
from dataclasses import dataclass
from typing import FrozenSet, List, Tuple, Iterable, Optional
from fractions import Fraction
import re
import numpy as np
import copy

from polynomial import Polinomio
from monomial import Monomio
import Parser
import state
import sys


def isGreater(a: Monomio, b: Monomio) -> bool:
    # print('Comparing monomials...')
    if state.M == []:
        # print('Using lex order')
        if a.grado() == b.grado():
            return a.coeff > b.coeff
        else: return a.grado() > b.grado()
    else:
        # print('Using matrix order')
        v = np.array(a.grado()) @ state.M
        # print('v:', v)
        w = np.array(b.grado()) @ state.M
        # print('w:', w)
        if tuple(v) == tuple(w):
            return a.coeff > b.coeff
        else: return tuple(v) > tuple(w)


def sort_monomials(monomials: List[Monomio]) -> List[Monomio]:
    """Ordina i monomi usando quicksort con isGreater come criterio."""
    if len(monomials) <= 1:
        return monomials.copy()
    
    def quicksort(arr: List[Monomio]) -> List[Monomio]:
        if len(arr) <= 1:
            return arr
        
        pivot_idx = len(arr) // 2
        pivot = arr[pivot_idx]
        left = []
        right = []
        equal = [pivot]
        
        for i, mono in enumerate(arr):
            if i == pivot_idx:
                continue
            if isGreater(mono, pivot):
                left.append(mono)
            elif isGreater(pivot, mono):
                right.append(mono)
            else:
                equal.append(mono)
        
        return quicksort(left) + equal + quicksort(right)
    
    return quicksort(monomials)

def normal_pol_form(p: Polinomio) -> Polinomio:
    # Get all monomials from the polynomial
    monomials = list(p.terms)
    sorted_monomials = sort_monomials(monomials)
    return Polinomio(sorted_monomials)

def lead_term(p: Polinomio) -> Monomio:
    if not p.terms:
        return Monomio.parse('0')
    normal_p = normal_pol_form(p)
    # Return the first monomial (which is the lead term after sorting)
    if normal_p.terms:
        return normal_p.terms[0]
    return Monomio.parse('0')

def lead_coeff(p: Polinomio) -> Optional[Fraction]:
    lt = lead_term(p)
    if lt is None:
        return None
    return lt.coeff
def lead_monomial(p: Polinomio) -> Optional[Monomio]:
    lt = lead_term(p)
    if lt is None:
        return None
    return Monomio(Fraction(1), dict(lt.powers))


def is_divisible_mons(a: Monomio, b: Monomio) -> bool:
    if a.is_zero() or b.is_zero():
        return False
    if not (b.powers.keys() - a.powers.keys()):
        boolean = True
        for v in b.powers.keys():
            boolean = (boolean and (b.powers[v] <=a.powers[v]))
        return boolean
    return False


def is_reducible_by(a: Polinomio, b: Polinomio) -> bool:
    lt_b = lead_term(b)
    if lt_b is None:
        return False
    # Check if any monomial in a is divisible by the lead term of b
    for m in a.terms:
        div = is_divisible_mons(m, lt_b)
        if div:
            return True
    return False

def division_mon(a: Monomio, b: Monomio) -> Monomio:
    if is_divisible_mons(a, b) != True:
        return Monomio.parse('0')
    
    # Divide i coefficienti
    if b.coeff == 0:
        return Monomio.parse('0')
    new_coeff = a.coeff / b.coeff  # divisione in virgola mobile
    
    # Divide le potenze
    new_powers = {}
    for k in a.powers.keys():
        b_exp = b.powers.get(k, 0)
        new_exp = a.powers[k] - b_exp
        if new_exp > 0:
            new_powers[k] = new_exp
    
    return Monomio(new_coeff, new_powers)


def order_pols(pols: List[Polinomio]) -> List[Polinomio]:
    """
    Ordina una lista di polinomi usando quick sort.
    Polinomio a > polinomio b se lead_term(a) > lead_term(b) secondo isGreater.
    """
    if len(pols) <= 1:
        return pols.copy()
    
    def quicksort(arr: List[Polinomio]) -> List[Polinomio]:
        if len(arr) <= 1:
            return arr
        pivot_idx = len(arr) // 2
        pivot = arr[pivot_idx]
        pivot_lead = lead_term(pivot)
        left = []
        right = []
        equal = [pivot]
        
        for i, pol in enumerate(arr):
            if i == pivot_idx:
                continue
            pol_lead = lead_term(pol)
            if isGreater(pol_lead, pivot_lead):
                left.append(pol)
            elif isGreater(pivot_lead, pol_lead):
                right.append(pol)
            else:
                equal.append(pol)
        return quicksort(left) + equal + quicksort(right)
    
    return quicksort(pols)


def reduce_pol_by_pol(f: Polinomio, g: Polinomio) -> Polinomio:
    """
    Riduce il polinomio f usando il polinomio g.
    Riduce i termini di f che sono divisibili per il lead term di g.
    """
    lt_g = lead_term(g)
    if lt_g is None or lt_g.is_zero():
        return f
    
    h = normal_pol_form(f)  # resto da ridurre
    r = Polinomio.zero()     # risultato accumulato (termini non riducibili)
    
    while h.terms:
        h = normal_pol_form(h)
        lt_h = lead_term(h)
        
        if is_divisible_mons(lt_h, lt_g):
            # Il lead term di h è divisibile per il lead term di g
            quotient = division_mon(lt_h, lt_g)
            h = normal_pol_form(h - Polinomio([quotient]) * g)
        else:
            # Il lead term di h non è riducibile, spostalo in r
            r = normal_pol_form(r + Polinomio([lt_h]))
            h = normal_pol_form(h - Polinomio([lt_h]))
    
    return r


def reduce_by_set(f: Polinomio, G: List[Polinomio]) -> Polinomio:
    reduced_f = copy.deepcopy(f)
    changed = True
    while changed:
        changed = False
        for g in G:
            new_reduced = reduce_pol_by_pol(reduced_f, g)
            if new_reduced != reduced_f:
                reduced_f = new_reduced
                changed = True
                break
    return reduced_f

def gcd_mon(a: Monomio, b: Monomio) -> Monomio:
    """
    Calcola il massimo comun divisore (GCD) di due monomi.
    Il GCD ha esponente minimo per ogni variabile.
    """
    gcd_powers = {}
    common_vars = set(a.powers.keys()).intersection(set(b.powers.keys()))
    for var in common_vars:
        gcd_powers[var] = min(a.powers[var], b.powers[var])
    return Monomio(Fraction(1), gcd_powers)


def are_coprime(a: Monomio, b: Monomio) -> bool:
    """
    Verifica se due monomi sono coprimi (GCD = 1).
    Due monomi sono coprimi se non hanno variabili in comune.
    """
    gcd = gcd_mon(a, b)
    return not gcd.powers  # True se gcd.powers è vuoto (= costante 1)


def mcm_mon(a: Monomio, b: Monomio) -> Monomio:
    """
    Calcola il minimo comune multiplo (LCM) di due monomi.
    Il LCM ha esponente massimo per ogni variabile.
    """
    c_powers = {}
    all_vars = set(a.powers.keys()).union(set(b.powers.keys()))
    for var in all_vars:
        power_a = a.powers.get(var, 0)
        power_b = b.powers.get(var, 0)
        c_powers[var] = max(power_a, power_b)
    return Monomio(Fraction(1), c_powers)


def apply_buchberger_criteria(pairs: List[Tuple[int, int]], G: List[Polinomio]) -> List[Tuple[int, int]]:
    """
    Applica i criteri di Buchberger per eliminare coppie ridondanti.
    
    Criterio 1 (Prodotto): Se gcd(LM(f), LM(g)) = 1, scarta (f,g).
    Criterio 2 (Catena): Se esiste h tale che LM(h) | lcm(LM(f), LM(g))
                         e le coppie (f,h) e (g,h) non sono in pairs, scarta (f,g).
    """
    filtered_pairs = []
    
    for i, j in pairs:
        lt_i = lead_term(G[i])
        lt_j = lead_term(G[j])
        
        # Criterio 1: prodotto (coprimi)
        if are_coprime(lt_i, lt_j):
            continue  # Scarta questa coppia
        
        # Criterio 2: catena
        lcm_ij = mcm_mon(lt_i, lt_j)
        discard = False
        
        for k in range(len(G)):
            if k == i or k == j:
                continue
            
            lt_k = lead_term(G[k])
            
            # Verifica se LM(k) divide lcm(LM(i), LM(j))
            if is_divisible_mons(lcm_ij, lt_k):
                # Verifica se le coppie (i,k) e (j,k) non sono in pairs
                pair_ik_in = (min(i, k), max(i, k)) in pairs or (max(i, k), min(i, k)) in pairs
                pair_jk_in = (min(j, k), max(j, k)) in pairs or (max(j, k), min(j, k)) in pairs
                
                if not pair_ik_in and not pair_jk_in:
                    # Entrambe le coppie sono state già processate
                    discard = True
                    break
        
        if not discard:
            filtered_pairs.append((i, j))
    
    return filtered_pairs


def sPol(f: Polinomio, g: Polinomio) -> Polinomio:
    m = mcm_mon(lead_term(f), lead_term(g))
    p = Polinomio([division_mon(m, lead_term(f))])*f - Polinomio([division_mon(m, lead_term(g))])*g
    return p


def total_degree(m: Monomio) -> int:
    """Calcola il grado totale di un monomio (somma di tutti gli esponenti)."""
    return sum(m.powers.values())


def poly_degree(p: Polinomio) -> int:
    """Calcola il grado totale di un polinomio (massimo tra i gradi dei monomi)."""
    if not p.terms:
        return -1  # Polinomio zero
    return max(total_degree(m) for m in p.terms)


def sugar_degree_of_spol(i: int, j: int, G: List[Polinomio], sugar: List[int]) -> int:
    """
    Calcola il sugar degree dell'S-polinomio S(G[i], G[j]).
    
    Sugar degree di S(f,g) = max(deg(lcm(LM(f), LM(g))) - deg(LM(f)) + sugar(f),
                                  deg(lcm(LM(f), LM(g))) - deg(LM(g)) + sugar(g))
    """
    lt_i = lead_term(G[i])
    lt_j = lead_term(G[j])
    lcm = mcm_mon(lt_i, lt_j)
    
    deg_lcm = total_degree(lcm)
    deg_lt_i = total_degree(lt_i)
    deg_lt_j = total_degree(lt_j)
    
    sugar_spol = max(deg_lcm - deg_lt_i + sugar[i], 
                     deg_lcm - deg_lt_j + sugar[j])
    
    return sugar_spol


def grobner_basis(F: List[Polinomio]) -> List[Polinomio]:
    """
    Calcola la base di Gröbner usando l'algoritmo di Buchberger con:
    - Criteri di Buchberger per filtrare coppie inutili
    - Strategia sugar per selezione ottimale delle coppie
    """
    G = copy.deepcopy(F)
    
    # Inizializza sugar degrees: per i generatori originali è il loro grado totale
    sugar = [poly_degree(f) for f in G]
    
    # Crea coppie iniziali con il loro sugar degree
    # Formato: (sugar_degree, i, j)
    pairs = []
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            sugar_spol = sugar_degree_of_spol(i, j, G, sugar)
            pairs.append((sugar_spol, i, j))
    
    while pairs:
        # Ordina le coppie per sugar degree (strategia sugar)
        pairs.sort(key=lambda x: x[0])
        
        # Applica i criteri di Buchberger per filtrare coppie inutili
        # Estrai solo gli indici per apply_buchberger_criteria
        index_pairs = [(i, j) for _, i, j in pairs]
        filtered_index_pairs = apply_buchberger_criteria(index_pairs, G)
        
        # Ricostruisci pairs con sugar degrees
        filtered_pairs = []
        for i, j in filtered_index_pairs:
            # Trova il sugar degree corrispondente
            sugar_spol = sugar_degree_of_spol(i, j, G, sugar)
            filtered_pairs.append((sugar_spol, i, j))
        
        pairs = filtered_pairs
        
        if not pairs:
            break
        
        # Prendi la coppia con sugar degree minimo
        sugar_spol, i, j = pairs.pop(0)
        
        S = sPol(G[i], G[j])
        R = reduce_by_set(S, G)
        
        if not R.terms:
            continue
        
        # Aggiungi il nuovo polinomio
        G.append(R)
        new_idx = len(G) - 1
        
        # Il sugar degree del nuovo polinomio è il sugar dell'S-polinomio
        sugar.append(sugar_spol)
        
        # Aggiungi nuove coppie con il nuovo polinomio
        for k in range(new_idx):
            new_sugar_spol = sugar_degree_of_spol(k, new_idx, G, sugar)
            pairs.append((new_sugar_spol, k, new_idx))
    
    return order_pols(G)

def normalized_form(F: List[Polinomio]) -> List[Polinomio]:
    """Normalizza i polinomi in F in modo che il coefficiente del termine principale sia 1"""
    normalized_F = []
    for p in F:
        lt = lead_term(p)
        if lt is None or lt.coeff == 0:
            normalized_F.append(p)
            continue
        coeff_inv = Fraction(1) / lt.coeff
        normalized_terms = [Monomio(m.coeff * coeff_inv, m.powers) for m in p.terms]
        normalized_F.append(Polinomio(normalized_terms))
    return normalized_F

def remove_zero_polys(F: List[Polinomio]) -> List[Polinomio]:
    """Rimuove i polinomi nulli dalla lista F"""
    return [p for p in F if p.terms]

def reduced_grobner_basis(F: List[Polinomio]) -> List[Polinomio]:
    """Calcola la base di Gröbner ridotta"""
    # Passo 1: Calcola la base di Gröbner
    G = grobner_basis(F)
    # Passo 2: Riduci ogni polinomio nella base rispetto agli altri
    copyed_G = []
    while copyed_G != G:
        copyed_G = copy.deepcopy(G)
        for i in range(len(G)):
            others = G[:i] + G[i+1:]
            G[i] = reduce_by_set(G[i], others)
    return normalized_form(remove_zero_polys(order_pols(G)))

def polinomial_sum(F: List[Polinomio]) -> Polinomio:
    result = Polinomio.parse('0')
    for p in F:
        result = normal_pol_form(result + p)
    return result

def polinomial_product(F: List[Polinomio]) -> Polinomio:
    result = Polinomio.parse('1')
    for p in F:
        result = normal_pol_form(result * p)
    return result
"""
n = Monomio.parse(str(input()))
m = Monomio.parse(str(input()))
print(is_divisible_mons(n, m))
#print(division_mon(n, m))
print(mcm_mon(n, m))
"""

print('Inserire variabili:')
state.ringvar = str(input()).split(',')

"""
f = lead_term(Polinomio.parse(input()))
g = lead_term(Polinomio.parse(input()))
sPol(f,g)
"""





