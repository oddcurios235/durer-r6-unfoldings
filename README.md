# Two Original-Edge Unfoldings for a Family of Parallelogram-Base Prismatoids at Every Positive Height

**Y. T. Miah** · Mathematical preprint dated 5 November 2026

[Read the manuscript (PDF)](main.pdf) · [LaTeX source](main.tex)

## About

This preprint constructs two explicitly defined, nonoverlapping original-edge unfoldings for a family of convex prismatoids at **every positive height**. The lower base is a parallelogram, and the entire upper polygon is strictly convex, with two anchored sides parallel to the lower edge directions and a corner chain monotone in the corresponding oblique coordinates. The anchor length is at least both lower side lengths. Both specified cut trees unfold the complete polyhedron while retaining its quadrilateral lateral facets without diagonal cuts.

The argument uses Pinciu's theorem for a face and its direct edge-neighbours, together with new separation arguments for the remaining faces and lower base. In the broader prismatoid literature, Radons established an edge-unfolding result for nested prismatoids, while Bian, Demaine and Madhukara treated tall and certain rectangular-base prismatoids. The present result concerns the stated, restricted family; the paper also gives an unbounded rational family whose two bases have disjoint horizontal projections.

## Contents

- [`main.pdf`](main.pdf) and [`main.tex`](main.tex): the eight-page manuscript and its complete, self-contained LaTeX source, including the bibliography.
- [`code/`](code/): six standard-library Python files for finite exact checks, including the geometric checker used by both development suites.
- [`certificates/`](certificates/): four JSON files containing the reference cases, verification results and structural/algebraic certificates.

## Reproducible exact checks

Use Python **3.8 or later**, without third-party packages. From the repository root run:

```sh
python code/run_all_checks.py
```

Run with assertions enabled; **do not use `python -O`**, since assertions implement verification conditions. The combined command runs the two complete-development suites and the algebraic and structural checks. They compare recomputed results with the supplied JSON certificates. The main suites cover 676 and 748 developments (1,424 in total, from 712 rational inputs) and check 31,780 and 37,396 complete face pairs, respectively. The code checks supporting facets, cut-tree structure, planar developments and convex-polygon separation using exact rational arithmetic and outward interval bounds.

The two companion checks cover algebraic identities, rational geometric placements and a structural witness. The geometric checker, `code/prismatoid_checker.py`, implements hull, cut-tree, development and face-separation tests used by the supplied suites. These finite computations complement the analytic proof of the theorem for all positive heights.

The component checks may also be run separately:

```sh
python code/r6_main_checks.py --check
python code/r6_additional_checks.py --check
python code/algebra_checks.py --check
python code/structural_checks.py --check
```

## Build the manuscript

With a LaTeX installation containing pdfLaTeX and the packages named in `main.tex`, run from the repository root:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The bibliography is contained in `main.tex`; BibTeX is not needed. The supplied `main.pdf` is the designated PDF for this source package.

## Citation

Y. T. Miah, *Two Original-Edge Unfoldings for a Family of Parallelogram-Base Prismatoids at Every Positive Height* (manuscript dated 5 November 2026). Machine-readable citation metadata is provided in [`CITATION.cff`](CITATION.cff).

Zenodo preprint: [10.5281/zenodo.23993991](https://doi.org/10.5281/zenodo.23993991).

## Licences

The original manuscript, repository documentation and four JSON certificates are available under **CC BY 4.0**; the six Python files in `code/` are available under **BSD 3-Clause**. The file-by-file scope and terms are in [`LICENSE.md`](LICENSE.md).
