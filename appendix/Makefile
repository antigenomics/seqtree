# Compile the E-value derivation. Requires a TeX distribution (latexmk + lualatex + bibtex)
# and rsvg-convert (librsvg) to convert the referenced SVG plots to PDF for \includegraphics.
# lualatex (not pdflatex) is required for the OldStandard OpenType text+math fonts; the
# compiler is pinned in .latexmkrc ($pdf_mode = 4) so a bare `latexmk` also picks it up.
FIGS = epitope_detection.pdf evalue_matrix.pdf mhc1_rocpr.pdf mhc2_rocpr.pdf

evalue.pdf: evalue.tex refs.bib $(FIGS)
	latexmk -lualatex -interaction=nonstopmode -halt-on-error evalue.tex

# Benchmark SVGs (bench/figures) -> PDF for inclusion as figures.
%.pdf: ../bench/figures/%.svg
	rsvg-convert -f pdf -o $@ $<

.PHONY: clean
clean:
	latexmk -C evalue.tex
	rm -f $(FIGS)
