from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, List
from fractions import Fraction

from lark import Lark, Transformer, v_args

from monomial import Monomio
import monomial
import state


POLY_DSL_GRAMMAR = r"""
?start: expr
?expr: sum

?sum: product ((PLUS|MINUS) product)*        -> sum_chain

?product: unary ( STAR unary )*              -> mul_chain

?unary: PLUS unary                           -> uplus
      | MINUS unary                          -> uminus
      | power

?power: atom ( CIRCUMFLEX UNSIGNED_INT )?    -> pow

?atom: UNSIGNED_INT                          -> int_lit
     | VAR                                   -> var
     | "(" expr ")"                          -> group

PLUS: "+"
MINUS: "-"
STAR: "*"
CIRCUMFLEX: "^"

VAR: /[a-zA-Z][a-zA-Z0-9_]*/
UNSIGNED_INT: /0|[1-9]\d*/

%import common.WS
%ignore WS
"""


def normal_pol_form(p: Polinomio) -> Polinomio:
    # Get all monomials from the polynomial
    monomials = list(p.terms)
    sorted_monomials = monomial.sort_monomials(monomials)
    return Polinomio(sorted_monomials)

@dataclass(frozen=True, slots=True)
class Polinomio:
    terms: List[Monomio]
    _parser = None

    def __post_init__(self):
        # Combine like terms (monomials with the same signature)
        sig_to_coeff: Dict[Tuple[Tuple[str, int], ...], Fraction] = {}
        for m in (self.terms or []):
            sig = m.signature()
            sig_to_coeff[sig] = sig_to_coeff.get(sig, Fraction(0)) + m.coeff
        
        # Create list of non-zero monomials
        cleaned = [Monomio(c, dict(sig)) for sig, c in sig_to_coeff.items() if c != 0]
        object.__setattr__(self, "terms", cleaned)

    # ---- factory ----
    @classmethod
    def zero(cls) -> "Polinomio":
        return cls([])

    @classmethod
    def const(cls, n) -> "Polinomio":
        return cls.from_monom(Monomio(Fraction(n), {}))

    @classmethod
    def var(cls, name: str) -> "Polinomio":
        return cls.from_monom(Monomio(Fraction(1), {name: 1}))

    @classmethod
    def from_monom(cls, m: Monomio) -> "Polinomio":
        if m.coeff == 0:
            return cls.zero()
        return cls([m])

    @classmethod
    def parse(cls, s: str) -> "Polinomio":
        # Nota: durante lo sviluppo è più sicuro NON cachare.
        # Se vuoi cache, ok, ma assicurati di ricreare il parser dopo modifiche alla grammatica.
        parser = Lark(
            POLY_DSL_GRAMMAR,
            parser="lalr",
            transformer=PolinomTransformer(cls),
        )
        return parser.parse(s)

    # ---- algebra ----
    def __add__(self, other: "Polinomio") -> "Polinomio":
        # Combine all monomials from both polynomials
        all_terms = list(self.terms) + list(other.terms)
        return Polinomio(all_terms)

    def __sub__(self, other: "Polinomio") -> "Polinomio":
        # Add negated terms from other
        all_terms = list(self.terms) + [-m for m in other.terms]
        return Polinomio(all_terms)

    def __neg__(self) -> "Polinomio":
        return Polinomio([-m for m in self.terms])

    def __mul__(self, other: "Polinomio") -> "Polinomio":
        if not self.terms or not other.terms:
            return Polinomio.zero()

        # Multiply each monomial from self with each monomial from other
        products = []
        for m_a in self.terms:
            for m_b in other.terms:
                products.append(m_a * m_b)
        return Polinomio(products)

    def __pow__(self, n: int) -> "Polinomio":
        n = int(n)
        if n < 0:
            raise ValueError("Esponenti negativi non supportati nei polinomi.")
        if n == 0:
            if not self.terms:
                raise ValueError("0^0 non è definito nel DSL.")
            return Polinomio.const(1)
        if not self.terms:
            return Polinomio.zero()

        # exp by squaring
        result = Polinomio.const(1)
        base = self
        e = n
        while e > 0:
            if e & 1:
                result = result * base
            base = base * base
            e >>= 1
        return result

    # ---- rendering ----
    @staticmethod
    def _sort_key(sig: Tuple[Tuple[str, int], ...]):
        total_deg = sum(e for _, e in sig)
        return (-total_deg, sig)

    def __str__(self) -> str:
        if not self.terms:
            return "0"

        # Sort monomials by degree (descending) and then by signature
        #items = sorted(self.terms, key=lambda m: self._sort_key(m.signature()))
        items = normal_pol_form(self).terms
        out = []

        for i, m in enumerate(items):
            abs_m = Monomio(abs(m.coeff), m.powers)
            term_str = str(abs_m)

            if i == 0:
                out.append(term_str if m.coeff >= 0 else "-" + term_str)
            else:
                out.append(("+" if m.coeff >= 0 else "-") + term_str)

        return "".join(out)


@v_args(inline=True)
class PolinomTransformer(Transformer):
    def __init__(self, poly_cls: type[Polinomio]):
        super().__init__()
        self.P = poly_cls

    def UNSIGNED_INT(self, tok):
        return int(tok)

    def int_lit(self, n):
        return self.P.const(int(n))

    def var(self, tok):
        return self.P.var(str(tok))

    def group(self, expr):
        return expr

    def uplus(self, x):
        return x

    def uminus(self, x):
        return -x

    def mul_chain(self, first, *rest):
        out = first
        for x in rest:
            if getattr(x, "type", None) == "STAR":
                continue
            out = out * x
        return out

    def sum_chain(self, first, *tail):
        # tail: (PLUS|MINUS, product, PLUS|MINUS, product, ...)
        if len(tail) % 2 != 0:
            raise ValueError(f"sum_chain: coda non in coppie op/term: {tail}")

        out = first
        i = 0
        while i < len(tail):
            op = tail[i].type   # "PLUS" o "MINUS"
            term = tail[i + 1]
            out = out + term if op == "PLUS" else out - term
            i += 2
        return out

    def pow(self, base, *rest):
        if not rest:
            return base
        exp = rest[-1]
        return base ** int(exp)

