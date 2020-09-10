pdflatex ./figure_3.tex
pdflatex ./figure_4.tex
pdflatex ./figure_5.tex
pdflatex -output-directory='../' figures_tables.tex
SetFile -a V *.aux *.log *.synctex.gz *.bbl *.blg *.out *.toc *.acn *.glo *.lof *.lot *.maf *.mtc* *.ntn *.xdy *.flg
open ../figures_tables.pdf
