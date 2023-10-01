import os
from collections import Counter
from typing import NamedTuple, List, Tuple
import pandas as pd

from antibodies import Ab, Class


class SummaryEntry(NamedTuple):
    class_name: str
    subject: str
    replicate: int
    matches: int
    total: int
    clonotype_matches: int
    total_clonotypes: int


class ClonotypeEntry(NamedTuple):
    class_name: str
    subject: str
    replicate: int
    count: int
    clonotype: Ab


def extract_match_data(
    antibodies: List[Ab], classes: List[Class], subject: str, replicate: int
) -> Tuple[List[SummaryEntry], List[ClonotypeEntry]]:
    summary_entries = []
    clonotype_entries = []
    ab_counter = Counter(antibodies)
    for cl in classes:
        matches = cl.matches(ab_counter.keys())
        summary_entries.append(
            SummaryEntry(cl.name, subject, replicate,
                         sum(ab_counter[ab] for ab in matches), len(antibodies),
                         len(matches), len(ab_counter))  # fmt: skip
        )
        for ab in matches:
            clonotype_entries.append(
                ClonotypeEntry(cl.name, subject, replicate, ab_counter[ab], ab)
            )
    return summary_entries, clonotype_entries


def save_summary(entries: List[SummaryEntry]) -> None:
    append_df = pd.DataFrame(entries).rename(columns={"class_name": "class"})
    summary_df = pd.read_csv("results/summary.tsv", sep="\t", dtype={"subject": str})
    summary_df = pd.concat([summary_df, append_df], ignore_index=True)
    summary_df.to_csv("results/summary.tsv", sep="\t", index=False)


def save_clonotypes(entries: List[ClonotypeEntry]) -> None:
    with open("results/matches.tsv", "a") as f:
        for entry in entries:
            f.write("\t".join(map(str, (*entry[:-1], *entry.clonotype))) + "\n")


def initialize_results_file(path: str, columns: Tuple) -> None:
    pd.DataFrame(None, columns=columns).to_csv(path, sep="\t", index=False)


def initialize_results_files() -> None:
    os.makedirs("results", exist_ok=True)
    initialize_results_file("results/summary.tsv", ("class", *SummaryEntry._fields[1:]))
    initialize_results_file(
        "results/matches.tsv", ("class", *ClonotypeEntry._fields[1:-1], *Ab._fields)
    )


def identify_unsearched_classes(
    subject: str, replicate: int, classes: List[Class]
) -> List[Class]:
    summary_df = pd.read_csv("results/summary.tsv", sep="\t", dtype={"subject": str})
    searched_class_names = set(
        summary_df[
            (summary_df["subject"] == subject) & (summary_df["replicate"] == replicate)
        ]["class"]
    )
    return [cl for cl in classes if cl.name not in searched_class_names]
