# Build the manuscript and check fonts.  Figures come from ../figures (run the fig*.py scripts there first).
main.pdf: main.tex refs.bib ../figures/fig1_model.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	@pdffonts main.pdf | awk 'NR>2 && ($$2 ~ /Type/ && $$3 == "3" || $$0 ~ / no /){bad=1; print} END{if(bad) print "FONT PROBLEM"; else print "fonts ok"}'
	@grep -c "undefined" main.log || true
clean:
	latexmk -C main.tex
