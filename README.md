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


The numerical section is intentionally brief. The current mathematical paper is meant to establish the measurement reduction, normal kernel, norm, non-compactness, global commutator, hyperbolic/Berezin identification and continuous spectral multiplier. A more extensive numerical-analysis contribution can be developed separately.

