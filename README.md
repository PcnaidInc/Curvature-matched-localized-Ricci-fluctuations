# Black-hole interior quantum-fluctuation research

**Status:** author-review manuscript and reproducibility material. Not independently peer reviewed.
Prepared repository files; this README does not assert that a remote repository has been published.

## Manuscript

* `reports/BlackHoles_Revised_Manuscript.pdf`: identified article plus supplementary methods.
* `manuscript/`: editable LaTeX, BibTeX, compiled bibliography, and vector figures.
* Revision date: 28 September 2026. The four declaration headings and explicit figure/table callouts
  have been added. Scientific equations, table values, and scope are unchanged.

The author is Abdul Badran (Independent Researcher). Conceptual questions and project direction are
his contributions. Mathematics, code, numerical execution, analysis, and drafting used substantial
AI assistance, as disclosed in the manuscript. AI systems are not authors.

## Reproducibility studies

| Folder | Scope |
|---|---|
| `studies/localized-ricci` | Dataset supporting the primary manuscript: one minimally coupled scalar, prepared state, finite Schwarzschild interior, sampled Ricci projections. |
| `studies/curvature-coupling` | Separate fixed-state coupling control; not silently combined with the original dataset. |
| `studies/optical-reference` | Subsequent causal optical benchmark; not a calculated new black-hole quantum variance. |

Each folder contains code, configurations, results, accepted and rejected checks, and a file manifest.
Do not delete failed pilots or treat them as accepted results. Each `reproduce.py` separates hashing,
small checks, saved-result assembly, and full computation where applicable. Dependencies are listed
per study. This publication-preparation step did not rerun a production scientific calculation.

## Scope

These results do not solve a singularity endpoint, establish universal quantum dominance, model a
rotating quantum ring, or calculate singularity hit/bypass probabilities. The original auxiliary state
preparation is not an astrophysical collapse model. Numerical checks are not total physical error bars.

## Data location and reuse

The verified repository URL and exact commit/tag must be added to the manuscript after publication.
No DOI is assigned by creating a GitHub repository. A later archival DOI may reference an exact release.
No blanket software/data reuse licence is granted by this preparation; see `RIGHTS_AND_REUSE.md`.

The anonymous CQG backup, account screenshots, contact information, manuscript-submission forms,
internal Jira records, and private author audit are intentionally **not** included here.
