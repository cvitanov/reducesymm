siminos/atlas/Chaos-v1/00ReadMe.txt
$Author$
$Date$
------------------------------------------------------------
submit before April 16, 2012

Slice & chart pedagogical introduction
====================
	Reduction of continuous symmetries of chaotic flows
	by the method of slices

Predrag Cvitanovic et al

-------------------------------------------------------------
Predrag ver. 0.2, 						Feb 22 2012

Focus issue for Chaos: "50 Years of Chaos"
To submit your paper to Chaos, please go to our web-based peer review and
submission system at http://chaos.peerx-press.org. The website is
currently open to accept papers for this focus issue, and you may submit
your papers at any point through the manuscript due date,

                April 16, 2012.

Please note the following items requested at submission:
    (1) a cover letter
        siminos/atlas/letters/submLett.doc
        title, authors, contact information, and mention this Focus issue
    (2) the article file in LaTeX
    (3) individual figure files in high-quality PDF formats

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

Chaos resubmission  (here for future use)
----------------
                                        2011-04-28
Chaos referee responses, edits
----------------------------
Predrag:						      Apr 1 2009
    generate colored comments for referees,
        PCedit\{...\}, etc in defslice.tex,
    then go back to B&W for the final resubmission

Chaos submission (here for future use)
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
