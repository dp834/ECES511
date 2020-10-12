

MATLAB_FILES = Homework2_Special.m
MATLAB_MAIN  = Homework2_Special.m

IMG_OUTPUT = images/

LATEX_SRC  = Homework-special-1.tex
PDF_OUTPUT = Homework-special-1

run: $(MATLAB_FILES) $(MATLAB_MAIN)
	mkdir -p $(IMG_OUTPUT)
	matlab -nodisplay -nosplash -batch "run('$(MATLAB_MAIN)')"

pdf: $(LATEX_SRC)
	pdflatex -jobname $(PDF_OUTPUT) $(LATEX_SRC)

clean:
	@rm -rf $(IMG_OUTPUT)/* $(PDF_OUTPUT).pdf *.aux *.log &2> /dev/null || true
