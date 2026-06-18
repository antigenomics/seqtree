# Compile the E-value derivation. Requires a TeX distribution (latexmk + pdflatex + bibtex).
evalue.pdf: evalue.tex refs.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error evalue.tex

.PHONY: clean
clean:
	latexmk -C evalue.tex
