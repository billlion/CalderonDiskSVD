# Multiplicity-normalized linearized Calderon data on the disk

Working LaTeX project intended for submission to **Inverse Problems**.

## Main files

- `main.tex` - manuscript source.
- `references.bib` - BibTeX database; every entry is still marked for human audit.
- `main.pdf` - compiled working draft (16 pages in the bundled local fallback style).
- `main.bbl` - generated bibliography, included for portability.
- `figures/` - the two figures used in the paper, plus one additional exploratory figure.
- `iopart-num-titles.bst` - a modified bibtex style like the iopart but giving names of papers
- `code/calderon_hyperbolic_model.m` - exploratory MATLAB finite-section calculation.
- `code/check_global_commutator.wl` - Mathematica check of the global kernel identity.
- `TODO.md` - pre-submission mathematical, numerical, bibliographic and publication checks.

## IOP class file

The current IOP Publishing class is `iopjournal.cls`. The source begins with

```tex
\IfFileExists{iopjournal.cls}
  {\documentclass{iopjournal}}
  {\documentclass{iopjournal-draft}}
```

The bundled `iopjournal-draft.cls` is **not** an official IOP file. It is only a small local fallback that allows the working draft to compile when the official class is unavailable. Before submission, download the current IOP LaTeX template from

```text
https://publishingsupport.iopscience.iop.org/wp-content/uploads/2025/07/ioplatextemplate.zip
```

and place `iopjournal.cls` and the accompanying bibliography-style/support files in this directory (or install them in the local TeX tree). Then rebuild and inspect the official-layout PDF.

IOP does not require authors to use its class at initial submission, but the project is structured to use it as soon as the class is present.

## Building

A standard build is

```bash
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

or simply

```bash
make
```

The included `main.bbl` also makes it possible to inspect the manuscript if BibTeX is temporarily unavailable.

## Status of the draft

This is a research working draft, not a submission-ready manuscript. Red `TODO` boxes are deliberately visible. They identify, among other things:

- the exact Fourier/sign/normalization audit from the original linearized boundary map;
- treatment of the zero boundary mode and the conductivity-to-Schrodinger convention;
- the distinction between multiplicity normalization and physical/noise weighting;
- independent verification of the global PDE commutator;
- the exact primary source and normalization for the Berezin/Mehler-Fock multiplier;
- a modest but properly documented numerical convergence/residual study;
- bibliography, authorship, CRediT, repository and AI-disclosure checks.

The numerical section is intentionally brief. The current mathematical paper is meant to establish the measurement reduction, normal kernel, norm, non-compactness, global commutator, hyperbolic/Berezin identification and continuous spectral multiplier. A more extensive numerical-analysis contribution can be developed separately.

## AI-assisted development

The source contains a provisional disclosure of substantial ChatGPT assistance. It should be revised only after the human verification record, final model/version information, author list and contribution statement have been settled. The named human authors remain responsible for every theorem, calculation, citation, code result and originality claim.
