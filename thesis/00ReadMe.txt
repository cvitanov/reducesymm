siminos/thesis/00ReadMe.txt
$Author$ $Date$
----------------------------------------------------------------------------

Evangelos Siminos thesis
========================

./thesis/
	pdflatex thesis
	
../bibtex/siminos.bib
	one bibliography for all projects

../figs/
	one directory for thesis and blog figures
	.eps figures needed for compatibility with articles, ChaosBook

../figsSrc/
	one directory for source programs for all figures
	(or .txt files explain how to regenerate them if needed)

../inputs/
	one directory for all Siminos specific macros

--------------------------------------------------------------------------
TO FIX:

--------------------------------------------------------------------------
NOTES:

-- \figurespagetrue option generates hyperref errors
   if there are math symbols in the \caption[...]

Time stamp:
----------

To have svn time-stamp file "someFile.type", include the contents of
	thesis/chapter/svnHeader.txt
into the file, and then
	svn propset svn:keywords "Date Author" someFile.type

--------------------------------------------------------------------------
HISTORY:


						Predrag Jun 26 2008
	rescued flotsam from siminos/blog/*.tex for inclusion
	into thesis


						Predrag Oct 10 2007:
	installed  thesis configuration files (from Halcrow's setup)


--------------------------------------------------------------------------
FIXED:

-- made hyperlinks active
-- made figure references [..], not (..) as is GaTech default
   by option \usepackage[square]{natbib}





