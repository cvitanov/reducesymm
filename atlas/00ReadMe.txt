siminos/atlas/00ReadMe.txt
$Author: predrag $
$Date: 2011-08-07 09:53:16 -0400 (Sun, 07 Aug 2011) $
------------------------------------------------------------

Slice & chart letter
====================
	Reduction of continuous symmetries of chaotic flows
	by the method of slices

Predrag Cvitanovic

-------------------------------------------------------------
Predrag ver. 0.1, 						Dec 10 2011

-------------------------------------------------------------

- ../figs/
	all figs
- ../figsSrc/
	all fig source programs
- ../bibtex/siminos.bib

----- NOTES        ------------------------------------------

-------------------------------------------------------------

arXiv resubmission (here for future use)
------------------
arxiv article-id:       1101.????v2     2011-04-28
     for details, please see/edit ../arXiv-v2/00ReadMe.txt

JFM resubmission  (here for future use)
----------------
                                        2011-04-28
JFM referee responses, edits
----------------------------
Predrag:						      Apr 1 2009
    generate colored comments for referees,
        PCedit\{...\}, etc in defslice.tex,
    then go back to B&W for the final resubmission

JFM submission (here for future use)
--------------

arXiv submission (here for future use)
----------------
arxiv article-id:       1101.????
     for details, please see/edit ../arxiv-v1/00ReadMe.txt

===========================================================
Processing
----------

make sure that you are in repository siminos/atlas/ , then
for a B&W compile:
	\draftfalse \colorfigsfalse   in inputs/type.tex
for color version:
	\draftfalse %\colorfigsfalse
	then:
./update
    or
pdflatex main; bibtex main; pdflatex main; pdflatex main

Things to fix
-------------
    [when fixed, move the line to the FIXED section at the end]

= see 00edits.txt


Style files
-----------

Preparation of figures
----------------------

==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= Working titles
-- Frank Zappa and The Mothers of Invention
