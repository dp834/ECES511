LATEX_SRC  = Homework-4.tex
PDF_OUTPUT = Homework-4


pdf: $(LATEX_SRC)
	pdflatex -jobname $(PDF_OUTPUT) $(LATEX_SRC)

clean:
	@rm -rf $(PDF_OUTPUT).pdf *.aux *.log &2> /dev/null || true
