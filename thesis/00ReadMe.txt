siminos/thesis/00ReadMe.txt
$Author$ $Date$
-----------------------------------------------------------------

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

-----------------------------------------------------------------
TO FIX:

size of
    figs/ks_L22_long_orbit.eps (also in rpo_ks/figs)
    figs/ks_prox_eq_tw.eps
and is few others is disgusting - ~2MB

						Predrag  Mar  1 2009
Vladimir Ilyich Lenin:
``One Step Forward, Two Steps Back''
        dasbuch/book/Figs/lorenzPolar1.eps     70KB
        siminos/figs/lorenzPolar1.eps       3,320KB
            now much inferior resolution
        dasbuch/book/Figs/lorenz2Poinc.eps     73KB
        siminos/figs/lorenz2Poinc.eps       2,116KB
            a bit better resolution
        dasbuch/book/Figs/lorenz2Poinc2D.eps    2KB
        siminos/figs/lorenz2Poinc2D.eps        47KB
            better, but does it need 25 x larger size?
        dasbuch/book/Figs/lorenzPolarM...il1.eps  49KB
        siminos/figs/lorenzPolarManifDetail1.eps 798KB
            now much inferior resolution
        dasbuch/book/Figs/lorenzSaddle0.eps    53KB
        siminos/figs/lorenzSaddle0.eps      3,298KB
            now somewhat inferior resolution
            next one is thesis only, so far:
        siminos/figs/ksSO2inv145eqbTo0.eps    180KB
            inferior resolution
For a sense of scale: ChaosBook.pdf had about 800 pages
and 3MB before Rytis, Halcrow and Siminos figures which
have taken it over 10MB already.

						Evangelos  Mar  1 2009
The only useful observation is that gs displays the figures
correctly when invoked as
    gs -dEPSCrop lorenzSaddle0.eps
The resolution is screwed by epswrite (used by eps2eps), so we
cannot use it as output device for gs to get a correct eps. Do
you know any gs tricks?

						Predrag  Feb 27 2009
if you are using figures in
    siminos/thesis/figs/
please move them to the canonical location
    ../figs/
then
    svn rm siminos/thesis/figs/
(makes it easier to keep track of them)

						Predrag  Jun 27 2008

if we modify siminos-thesis.cls, we'll modify
	    \documentclass{gatech-thesis} to \documentclass{siminos-thesis} in
\input ../inputs/setupThesis    %% logical chores (nothing to edit).
Currently siminos-thesis.cls is the same as gatech-thesis.cls.


_________________________________________________________________
Name: svn:keywords
   + Author Date


-----------------------------------------------------------------
NOTES:

--  \figurespagetrue option generates hyperref errors
    if there are math symbols in the \caption[...]

--  When getting irrelevant output from svn diff, use
        svn diff -x -w
    (ignore all white change.)

Time stamp:
----------

To have svn time-stamp file "someFile.type", include the contents of
	thesis/chapter/svnHeader.txt
into the file, and then
	svn propset svn:keywords "Date Author" someFile.type

-----------------------------------------------------------------
HISTORY:

						Predrag Jun 26 2008
	rescued flotsam from siminos/blog/*.tex for inclusion
	into thesis

						Predrag Oct 10 2007:
	installed  thesis configuration files (from Halcrow's setup)


-----------------------------------------------------------------
FIXED:

-- made hyperlinks active
-- made figure references [..], not (..) as is GaTech default
   by option \usepackage[square]{natbib}
