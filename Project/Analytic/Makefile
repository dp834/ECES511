

MATLAB_FILES = analytic_sol.m
MATLAB_MAIN  = analytic_sol.m

IMG_OUTPUT = images/

LATEX_SRC  = analytic_sol.tex
PDF_OUTPUT = Analytic_Solution

run: $(MATLAB_FILES) $(MATLAB_MAIN)
	mkdir -p $(IMG_OUTPUT)
	matlab -nodisplay -nosplash -batch "run('$(MATLAB_MAIN)')"

pdf: $(LATEX_SRC)
	pdflatex -jobname $(PDF_OUTPUT) $(LATEX_SRC)

clean:
	@rm -rf $(IMG_OUTPUT)/* $(PDF_OUTPUT).pdf *.aux *.log *.nav *.toc *.snm *.out &2> /dev/null || true
