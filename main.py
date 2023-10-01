import os
import itertools
from typing import List

from ab_classes import classes
from antibodies import Gene, Ab, AA_ALPHABET, filter_inplace
from file_output import (
    initialize_results_files,
    identify_unsearched_classes,
    extract_match_data,
    save_summary,
    save_clonotypes,
)


def filter_briney_abs(antibodies: List[Ab]) -> List[Ab]:
    # Filter by antibody characteristics
    filter_inplace(
        antibodies,
        lambda ab: ab.chain == "heavy" and ab.productive and ab.isotype == "IgM",
    )

    # Filter antibodies with non-amino acid symbols in CDR
    n_x = filter_inplace(antibodies, lambda ab: "X" not in ab.cdr)
    n_stop = filter_inplace(antibodies, lambda ab: "*" not in ab.cdr)
    assert all(all(aa in AA_ALPHABET for aa in ab.cdr) for ab in antibodies)
    n_short = filter_inplace(antibodies, lambda ab: len(ab.cdr) >= 6)
    # Counts don't consider intersections of exclusion criteria
    print(f"CDRs: {n_x} unresolved | {n_stop} w/ stop codons | {n_short} short")

    check_gene_names(antibodies)

    return antibodies


def read_briney_abs(annotation_fn: str) -> List[Ab]:
    def _read_gene(s: str):
        return None if s == "-" else Gene(*s.split("*"))

    antibodies = []
    with open(annotation_fn) as f:
        lines = [line.split(",") for line in f.read().splitlines()[1:]]
        for line in lines:
            ab = Ab(
                chain=line[0],
                productive=line[1] == "yes",
                v=_read_gene(line[2]),
                d=_read_gene(line[3]),
                j=_read_gene(line[4]),
                cdr=line[5],
                isotype=line[6],
            )
            antibodies.append(ab)

    return antibodies


def check_gene_names(antibodies: List[Ab]) -> None:
    """Check that each Ab has IgH genes and alleles as expected."""
    for ab in antibodies:
        assert ab.v and ab.v.name.startswith("IGHV") and len(ab.v.allele) == 2
        if ab.d is not None:
            assert ab.d.name.startswith("IGHD") and len(ab.d.allele) == 2
        assert ab.j and ab.j.name.startswith("IGHJ") and len(ab.j.allele) == 2


if __name__ == "__main__":
    if not os.path.exists("results/summary.tsv"):
        initialize_results_files()

    # Iterate over each subject and replicate
    subjects = ("316188", "326650", "326651", "326713", "326737",
                "326780", "326797", "326907", "327059", "D103")  # fmt: skip
    for subject, replicate in itertools.product(subjects, range(1, 19)):
        # Locate annotated results for this subject and replicate and skip if it doesn't exist
        annotation_fn = f"data/{subject}/{replicate}.csv"
        if not os.path.exists(annotation_fn):
            print(f"SKIPPING subject {subject} replicate {replicate}: no data")
            continue

        # If class already searched for this subject and replicate, don't search that class
        unsearched_classes = identify_unsearched_classes(subject, replicate, classes)
        if not unsearched_classes:
            continue

        print(f"SEARCHING subject {subject} replicate {replicate}")
        antibodies = read_briney_abs(annotation_fn)
        antibodies = filter_briney_abs(antibodies)

        summary_entries, clonotype_entries = extract_match_data(
            antibodies, unsearched_classes, subject, replicate
        )
        save_summary(summary_entries)
        save_clonotypes(clonotype_entries)
