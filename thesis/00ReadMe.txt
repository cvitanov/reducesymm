siminos/thesis/00ReadMe.txt
$Author: predrag $ $Date: 2007-06-10 02:45:16 -0400 (Sun, 10 Jun 2007) $
----------------------------------------------------------------------------

Evangelos Siminos thesis
========================

./thesis/
	pdflatex thesis
	yupee!

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

-- make hyperlinks active
-- \figurespagetrue option generates hyperref errors

--------------------------------------------------------------------------
NOTES:

Time stamp:
----------

To have svn time-stamp file "someFile.type", include the contents of
	thesis/chapter/svnHeader.txt
into the file, and then
	svn propset svn:keywords "Date Author" someFile.type

--------------------------------------------------------------------------
HISTORY:

	Predrag Oct 10 2007:
	installed  thesis configuration files (from Halcrow's setup)
