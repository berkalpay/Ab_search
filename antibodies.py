import regex as re
from dataclasses import dataclass, field
from typing import NamedTuple, List, Tuple, Callable, Iterable


AA_ALPHABET = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}  # fmt: skip


def filter_inplace(iterable: List, function: Callable) -> int:
    pre_len = len(iterable)
    filtered = list(filter(function, iterable))
    post_len = len(filtered)
    iterable[:] = filtered
    return pre_len - post_len


def kabat_to_cdr(kabat: int) -> int:
    return kabat - 95 + 3


@dataclass(frozen=True)
class Gene:
    name: str
    allele: str = field(default=None)

    def __str__(self):
        if self.allele is not None:
            return f"{self.name}*{self.allele}"
        else:
            return self.name

    def is_a(self, other) -> bool:
        return self.name == other.name and (
            self.allele == other.allele or other.allele is None
        )

    def is_one_of(self, others) -> bool:
        return any(self.is_a(other) for other in others)


class Ab(NamedTuple):
    chain: str
    productive: bool
    v: Gene
    d: Gene
    j: Gene
    cdr: str
    isotype: str


@dataclass(frozen=True)
class Class:
    name: str
    vs: List[Gene] = field(default=None)
    ds: List[Gene] = field(default=None)
    js: List[Gene] = field(default=None)
    cdr_length: Tuple[int, int] = field(default=None)
    cdr_signature: Tuple[int, str] = field(default=None)

    def __str__(self):
        return self.name

    @property
    def cdr_signature_re(self):
        if self.cdr_signature is None:
            return None

        starting_kabat, signature = self.cdr_signature

        # Define what to look for before and after the CDR signature
        if starting_kabat is None:
            pre_signature = "(.*)"
        else:
            n_aa_before = kabat_to_cdr(starting_kabat) - 1
            pre_signature = ".{" + str(n_aa_before) + "}"

        return re.compile(f"({pre_signature})({signature})")

    def gene_match(self, ab: Ab, segment: str) -> bool:
        class_genes = getattr(self, f"{segment}s")
        ab_gene = getattr(ab, segment)
        return class_genes is None or (
            ab_gene is not None and ab_gene.is_one_of(class_genes)
        )

    def matches(self, antibodies: Iterable) -> List[Ab]:
        matches = filter(
            lambda ab: all(self.gene_match(ab, segment) for segment in ["v", "d", "j"]),
            antibodies,
        )
        if self.cdr_length is not None:
            matches = filter(
                lambda ab: self.cdr_length[0] <= len(ab.cdr) <= self.cdr_length[1],
                matches,
            )
        if self.cdr_signature_re is not None:
            matches = filter(
                lambda ab: self.cdr_signature_re.match(ab.cdr) is not None, matches
            )
        return list(matches)
