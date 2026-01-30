from __future__ import annotations
from dataclasses import dataclass
from typing import FrozenSet, List, Tuple, Iterable, Optional
import re

from polynomial import Polinomio
import state


def _split_top_level_commas(s: str) -> List[str]:
    parts: List[str] = []
    buf: List[str] = []
    depth = 0

    for ch in s:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise ValueError("Parentesi non bilanciate.")
            buf.append(ch)
        elif ch == "," and depth == 0:
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
        else:
            buf.append(ch)

    if depth != 0:
        raise ValueError("Parentesi non bilanciate.")

    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)
    return parts


def _poly_used_vars(p: Polinomio) -> FrozenSet[str]:
    used = set()
    for m in p.terms:
        for var in m.powers.keys():
            used.add(var)
    return frozenset(used)


@dataclass(frozen=True, slots=True)
class PolynomialRing:
    domain: str
    variables: FrozenSet[str]

    @classmethod
    def parse(cls, s: str) -> "PolynomialRing":
        s = s.strip()

        m = re.fullmatch(r"""\\?([A-Za-z]+)\s*\[\s*(.*?)\s*\]\s*""", s)
        if not m:
            raise ValueError(f"Formato anello non valido: {s!r} (atteso tipo \\\\R[x,y])")

        dom = m.group(1)
        varlist = m.group(2).strip()
        if not varlist:
            raise ValueError("Lista variabili vuota nell'anello.")

        vars_raw = [v.strip() for v in varlist.split(",")]
        if any(not v for v in vars_raw):
            raise ValueError("Lista variabili contiene elementi vuoti.")

        bad = [v for v in vars_raw if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", v)]
        if bad:
            raise ValueError(f"Variabili non valide: {bad}")

        return cls(domain=dom, variables=frozenset(vars_raw))


@dataclass(frozen=True, slots=True)
class Ambiente:
    ring: PolynomialRing

    @classmethod
    def from_ring_string(cls, ring_s: str) -> "Ambiente":
        return cls(PolynomialRing.parse(ring_s))

    def parse_poly(self, s: str) -> Polinomio:
        p = Polinomio.parse(s)
        used = _poly_used_vars(p)
        extra = sorted(used - self.ring.variables)
        if extra != []:
            raise ValueError(
                f"Il polinomio usa variabili non ammesse {extra}. "
                f"Variabili ammesse: {sorted(self.ring.variables)}"
            )
        return p

    def parse_ideal(self, s: str) -> "Ideale":
        return Ideale.parse(s, env=self)


@dataclass(frozen=True, slots=True)
class Ideale:
    env: Ambiente
    generators: Tuple[Polinomio, ...]

    @classmethod
    def parse(cls, s: str, env: Ambiente) -> "Ideale":
        raw = s.strip()
        if not (raw.startswith("<") and raw.endswith(">")):
            raise ValueError(f"Ideale non valido: {s!r} (atteso formato <p1,...,pn>)")

        inside = raw[1:-1].strip()
        if not inside:
            raise ValueError("Ideale senza generatori: < > non è ammesso.")
        
        gen_strs = _split_top_level_commas(inside)
        gens = tuple(env.parse_poly(gs) for gs in gen_strs)
        return cls(env=env, generators=gens)

    def __str__(self) -> str:
        return "<" + ", ".join(str(g) for g in self.generators) + ">"


