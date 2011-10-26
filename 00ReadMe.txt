siminos/00ReadMe.txt
$Author$ 
$Date$
--------------------------------------------------------

Evangelos Siminos thesis, publications and working notes
========================================================

thesis/
	pdflatex thesis
	
rpo_ks/
    Cvitanovi\'c, Davidchack and Siminos
    "State space geometry of a spatio-temporally chaotic
     Kuramoto-Sivashinsky flow"

CLE/00ReadMe.txt
    Siminos and Cvitanovi\'c
    "Continuous symmetry reduction
        and return maps for higher-dimensional flows"

bibtex/siminos.bib
	one bibliography for all projects

figs/
    one directory for thesis and blog figures .eps figures
    needed for compatibility with articles, ChaosBook.
    When creating a new version, please use the SAME name (so
    it propagates to publications etc. as well, without messing
    around with renaming it)

Fig/
    copied from ChaosBook.org (do not edit here)

figsSrc/
	one directory for source programs for all figures
	(or .txt files explain how to regenerate them if needed
     from the other repository, siminos????, with code)

inputs/
	one directory for all Siminos specific macros

thesis/

www/
    Siminos home pages (tres elegantes)

talks/
    all symmetry reduction talks
    Tufts10/    Predrag seminar, 2010-04-02

posters/

scripts/
    potentially useful scripts 

blog/
    all matters pertaining to symmetry reduction reading

lyapunov/
    all matters pertaining to 'covariant Lyapunov vectors'

froehlich/
    Stefan Froehlich blog, papers 

chao/
    Chao's blog

---------------------------------------------------------------
TO FIX:
- 2011-03-16: emaildict for Chao Shi reply-to address

---------------------------------------------------------------
NOTES:

-- maths classification for a paper about Lorenz system:
   MSC: Primary: 37C45, 37D40; Secondary: 37D45


Time stamp:
----------

To have svn time-stamp file "someFile.type", include the contents of
	thesis/chapter/svnHeader.txt
into the file, and then
	svn propset svn:keywords "Date Author" someFile.type

----------------------------------------------------------------
HISTORY:
						Predrag Jul 26 2008
    created this file

----------------------------------------------------------------
FIXED:
						Evangelos Dec 4 2010	
	emaildict for correct reply-to address

-- made ...
