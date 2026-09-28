#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
BIBTEX="${BIBTEX:-bibtex}"
pdflatex -interaction=nonstopmode -halt-on-error article.tex
"$BIBTEX" article
pdflatex -interaction=nonstopmode -halt-on-error article.tex
pdflatex -interaction=nonstopmode -halt-on-error article.tex
pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
