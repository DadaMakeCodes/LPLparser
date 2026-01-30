from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, List
from fractions import Fraction

from lark import Lark, Transformer, v_args
import state
import numpy as np

# ------------------ Grammar (FIXED PRECEDENCE) ------------------
# Precedenza desiderata:
#   ^  (più forte)
#   unario +/-
#   *
#
# Quindi: -x^2 => -(x^2)  (NON (-x)^2)
MONOM_DSL_GRAMMAR = r"""
?start: expr

?expr: product

?product: unary ( "*" unary )*                 -> mul

?unary: "+" unary                              -> uplus
      | "-" unary                              -> uminus
      | power

?power: atom ( "^" UNSIGNED_INT )?             -> pow

?atom: UNSIGNED_INT                            -> int_lit
     | VAR                                     -> var
     | "(" expr ")"                            -> group

VAR: /[a-zA-Z][a-zA-Z0-9_]*/
UNSIGNED_INT: /0|[1-9]\d*/

%import common.WS
%ignore WS
"""


# ------------------ Transformer ------------------
@v_args(inline=True)
class MonomTransformer(Transformer):
    @staticmethod
    def _norm(m):
        c, p = int(m[0]), dict(m[1])
        if c == 0:
            return [0, {}]
        p = {v: int(e) for v, e in p.items() if int(e) != 0}
        return [c, p]

    @staticmethod
    def _mul(a, b):
        a = MonomTransformer._norm(a)
        b = MonomTransformer._norm(b)
        c1, p1 = a
        c2, p2 = b
        if c1 == 0 or c2 == 0:
            return [0, {}]
        out_c = c1 * c2
        out_p = dict(p1)
        for v, e in p2.items():
            out_p[v] = out_p.get(v, 0) + e
        return MonomTransformer._norm([out_c, out_p])

    @staticmethod
    def _pow(m, n):
        m = MonomTransformer._norm(m)
        n = int(n)
        if n < 0:
            raise ValueError("Esponenti negativi non supportati nei monomi.")
        c, p = m

        if n == 0:
            # 0^0 esplicito => errore
            if c == 0 and not p:
                raise ValueError("0^0 non è definito nel DSL.")
            return [1, {}]

        if c == 0:
            return [0, {}]

        out_c = c ** n
        out_p = {v: e * n for v, e in p.items()}
        return MonomTransformer._norm([out_c, out_p])

    # ---- terminals ----
    def UNSIGNED_INT(self, tok):
        return int(tok)

    def int_lit(self, n):
        # n arriva già come int via UNSIGNED_INT()
        return [int(n), {}]

    def var(self, tok):
        return [1, {str(tok): 1}]

    def group(self, expr):
        return self._norm(expr)

    # ---- unary ----
    def uplus(self, m):
        return self._norm(m)

    def uminus(self, m):
        m = self._norm(m)
        return [-m[0], m[1]]

    # ---- ops ----
    def mul(self, first, *rest):
        out = first
        for m in rest:
            out = self._mul(out, m)
        return self._norm(out)

    def pow(self, base, exp=None):
        if exp is None:
            return self._norm(base)
        return self._pow(base, exp)


def isGreater(a: str, b: str) -> bool:
    return isMonGreater(Monomio(coeff=Fraction(1), powers={a: 1}), Monomio(coeff=Fraction(1), powers={b: 1}))
    
    

def isMonGreater(a: Monomio, b: Monomio) -> bool:
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
    """Ordina i monomi usando quicksort con isMonGreater come criterio."""
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
            if isMonGreater(mono, pivot):
                left.append(mono)
            elif isMonGreater(pivot, mono):
                right.append(mono)
            else:
                equal.append(mono)
        
        return quicksort(left) + equal + quicksort(right)
    
    return quicksort(monomials)



# ------------------ Monomio class ------------------
@dataclass(frozen=True, slots=True)
class Monomio:
    coeff: Fraction
    powers: Dict[str, int]

    _parser = Lark(
        MONOM_DSL_GRAMMAR,
        parser="lalr",
        transformer=MonomTransformer(),
    )

    def __post_init__(self):
        c = Fraction(self.coeff)
        p = dict(self.powers or {})

        if c == 0:
            object.__setattr__(self, "coeff", Fraction(0))
            object.__setattr__(self, "powers", {})
            return

        cleaned = {}
        for v, e in p.items():
            e = int(e)
            if e > 0:
                cleaned[v] = e
            elif e < 0:
                raise ValueError("Esponenti negativi non supportati nei monomi.")
        object.__setattr__(self, "coeff", c)
        object.__setattr__(self, "powers", cleaned)

    @classmethod
    def parse(cls, s: str) -> "Monomio":
        coeff, powers = cls._parser.parse(s)
        powers_dict = dict(powers)
        
        # Verifica che tutte le variabili siano in state.ringvar
        invalid_vars = [v for v in powers_dict.keys() if v not in state.ringvar]
        if invalid_vars:
            raise ValueError(
                f"Variabili non ammesse nel monomio: {invalid_vars}. "
                f"Variabili ammesse: {state.ringvar}"
            )
        
        return cls(Fraction(coeff), powers_dict)

    def to_list(self) -> list:
        return [self.coeff, dict(self.powers)]

    def ordered_items(self) -> Tuple[Tuple[str, int], ...]:
        """Ordina le coppie (variabile, esponente) usando quicksort con isGreater."""
        items = list(self.powers.items())
        
        if len(items) <= 1:
            return tuple(items)
        
        def quicksort(arr: List[Tuple[str, int]]) -> List[Tuple[str, int]]:
            if len(arr) <= 1:
                return arr
            
            pivot_idx = len(arr) // 2
            pivot = arr[pivot_idx]
            left = []
            right = []
            equal = [pivot]
            
            for i, item in enumerate(arr):
                if i == pivot_idx:
                    continue
                if isGreater(item[0], pivot[0]):
                    left.append(item)
                elif isGreater(pivot[0], item[0]):
                    right.append(item)
                else:
                    equal.append(item)
            
            return quicksort(left) + equal + quicksort(right)
        
        return tuple(quicksort(items))

    def __str__(self) -> str:
        if not self.powers:
            return str(self.coeff)

        parts = []
        for v, e in self.ordered_items():
            parts.append(v if e == 1 else f"{v}^{e}")
        core = "*".join(parts)

        if self.coeff == 1:
            return core
        if self.coeff == -1:
            return f"-{core}"
        return f"{self.coeff}*{core}"

    def __repr__(self) -> str:
        return f"Monomio(coeff={self.coeff}, powers={dict(self.ordered_items())})"

    def __neg__(self) -> "Monomio":
        return Monomio(-self.coeff, self.powers)

    def __mul__(self, other: "Monomio") -> "Monomio":
        if not isinstance(other, Monomio):
            return NotImplemented
        if self.coeff == 0 or other.coeff == 0:
            return Monomio(Fraction(0), {})
        out_c = self.coeff * other.coeff
        out_p = dict(self.powers)
        for v, e in other.powers.items():
            out_p[v] = out_p.get(v, 0) + e
        return Monomio(out_c, out_p)

    def __pow__(self, n: int) -> "Monomio":
        n = int(n)
        if n < 0:
            raise ValueError("Esponenti negativi non supportati nei monomi.")
        if n == 0:
            if self.coeff == 0 and not self.powers:
                raise ValueError("0^0 non è definito nel DSL.")
            return Monomio(Fraction(1), {})
        if self.coeff == 0:
            return Monomio(Fraction(0), {})
        out_c = self.coeff ** n
        out_p = {v: e * n for v, e in self.powers.items()}
        return Monomio(out_c, out_p)

    def is_zero(self) -> bool:
        return self.coeff == 0 and not self.powers
    
    def signature(self) -> Tuple[Tuple[str, int], ...]:
        return self.ordered_items()
    
    def grado(self) -> list:
        # Restituisce vettore completo di esponenti per tutte le variabili in state.ringvar
        return [self.powers.get(v, 0) for v in state.ringvar]
