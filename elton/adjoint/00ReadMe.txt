pipes/mapit/00ReadMe.txt
$Author: predrag $
$Date: 2015-05-04 10:57:54 -0400 (Mon, 04 May 2015) $
------------------------------------------------------------

Pipe flows slicing article
==========================
	Reduction of continuous symmetries of chaotic flows
	by the method of slices

Ashley Willis
Predrag Cvitanovic
...

-------------------------------------------------------------
proofs received:                      Feb 25 2013
    [x] pipes/slice/jfm-v3/flm1300075_PRF.pdf
Ashley requested to be featured       Feb 11 2013
    [x] pipes/slice/reviews/FocusOnFluids.txt
JFM received JFM-12-S-0339.R2         Jan 29 2013
Ashley JFM submitted  ver. 5.0        Jan 29 2013
    pipes/slice/jfm-final/
    [x] set up referee 1 comments     Jan 27 2013
        in pipes/blog/sliceRev.tex
    [ ] complete and sign a transfer of copyright form
    [x] accommodate suggestions of referee 1 in the
        production version of the paper. Referee 1:
        pipes/slice/reviews/Referee1-2.pdf
    [ ] any figures to appear in colour online are legible and sensible
        when printed in black and white.
JFM gave us 2 weeks deadline          Jan 23 2013
                           deliver by Feb  6 2013
JFM accepted:                         Jan 05 2013
Ashley JFM   resubmitted  ver. 4.0    Oct 31 2012
    pipes/slice/jfm-v2/
Predrag JFM resub. start  ver. 3.0    Jul 30 2012
JFM referee responses, edits          Jun 29 2012
    pipes/slice/reviews/
    generate colored comments for referees,
        PCedit\{...\}, etc in defslice.tex,
    then go back to B&W for the final resubmission
Ashley JFM submitted ver. 2.1,          Apr  2 2012
    read jfm-v1/ToEditor.txt, fromEditor.txt
Predrag ACHKW11.tex  ver. 2.1,          Mar 19 2012
Ashley arXiv submit  ver. 2.0,          Mar 16 2012
Predrag ACHKW11.tex  ver. 1.1,          Mar  5 2012
    merged all sections, ready to finalize
    but cannot before we get through footnotes
Predrag ACHKW11.tex  ver. 1.0,          Mar  4 2012
Ashley, Predrag ver. 0.4, 				Aug  7 2011
    abstract
Predrag ver. 0.3, 						May 13 2011
	subdivided into sections, now: 	pdflatex main
Predrag ver. 0.2, 						May 10 2011
Predrag ver. 0.1, 						Apr 28 2011

-------------------------------------------------------------

- ../figs/
	all figs
- ../figsSrc/
	all fig source programs
- ../bibtex/pipes.bib
	the one and only pipes.bib bibliography file
	enter new entries at the top (read info there)

----- NOTES        ------------------------------------------
Ashes    2011-05-12: There appears to be a jfm started,
	with a 1-page pdf.  The tex looks more substantial
	but isn't building.
Humbledt 2011-05-12: pdflatex slice  gives me 17 pages.
	There are some erros, I'll fix them

Ashes 2011-05-12: Which files or sections of files
	should we be working on?
Humbledt 2011-05-12: Edit whatever you want. The rule number 1 -
	always
cd  ~[PATH]/pipes/; svn up
	before you start editing, then
cd  ~[PATH]/pipes/; svn ci -m"descripive message"
	often - the minimizies conflicts. Subversions is good
	at merging files. You do not highlight changes in the
	early stages - you can always see what was changed by
	using diff. You probably want to install
kdesvn package
	or
http://rabbitvcs.org/ (have not tried that one)
	the graphical interfaces is very handy.

Humbledt 2011-05-12:
	For starters - do some very minor editing, and then
cd  ~[PATH]/pipes/; svn ci -m"this is my first test commit"
	so we are sure that uploading works as well

-------------------------------------------------------------

arXiv resubmission (here for future use)
------------------
arxiv article-id:       1101.????v2     2011-04-28
     for details, please see/edit ../arXiv-v2/00ReadMe.txt

JFM resubmission  (here for future use)
----------------

JFM submission (here for future use)
--------------

arXiv submission (here for future use)
----------------
arxiv article-id:       1101.????
     for details, please see/edit ../arxiv-v1/00ReadMe.txt

===========================================================
Processing
----------

make sure that you are in repository pipes/slice/ , then
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

The LaTeX2e style file jfm.cls together with a guide to its use
old/jfm2egui.tex, sample pages old/jfm2esam.tex
and a guide to editorial style and notation old/jfm2enot.tex
are from
   ftp://ftp.cup.cam.ac.uk/pub/texarchive/journals/latex/jfm-cls/

A JFM bibtex style file can be found, together with
sample files bibsamp.tex and bibsamp.bib, in the directories
   /pub/texarchive/journals/latex/jfm-cls/bst

Preparation of figures
----------------------

Space is at a premium so figures should be as small as possible
while displaying clearly all the information required. Large areas
of white space should be avoided and an aspect ratio chosen that
makes full use of the width of the page rather than being
long and narrow (the width available is 13 cm).

Line widths must be no smaller than 0.3 pt. Lines any thinner than
this will not reproduce well.

Figures published in JFM are relettered in the Journal's typeface,
so lettering inserted by authors will be replaced. If figures are
reduced in size please ensure that the lettering is still
easily readable.
Figure files

The preferred format for figure files is .eps or .tiff
at resolution 1200 dpi for lines, 600 dpi for greyscale
and 300 dpi for colour (which preferably should also be
in CMYK - cyan magenta yellow black - format). However,
most standard image formats such as pct, ppm, png, psd, Word,
ppt, CorelDraw, ChemDraw, AutoCAD can also be used, but not
customized output of software not designed for publishing purposes
such as Matlab, nor PDF.

Figures to be printed in black and white must be submitted
as black and white files.

==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= Working titles
-- Frank Zappa and The Mothers of Invention
-- Predrag 2011-04-28: halcrow/n00bs/n00bs.tex as a starting template
