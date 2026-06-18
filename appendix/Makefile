# Compile the E-value derivation. Requires a TeX distribution (latexmk + pdflatex + bibtex)
# and rsvg-convert (librsvg) to convert the referenced SVG plots to PDF for \includegraphics.
FIGS = epitope_detection.pdf

evalue.pdf: evalue.tex refs.bib $(FIGS)
	latexmk -pdf -interaction=nonstopmode -halt-on-error evalue.tex

# Benchmark SVGs (bench/figures) -> PDF for inclusion as figures.
epitope_detection.pdf: ../bench/figures/epitope_detection.svg
	rsvg-convert -f pdf -o $@ $<

.PHONY: clean
clean:
	latexmk -C evalue.tex
	rm -f $(FIGS)
