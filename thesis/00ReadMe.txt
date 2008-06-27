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

						Predrag  Jun 27 2008
Sorry about this - I'm quite sure I did not touch text like
`The complex Lorenz equations ...', and do not understand why
it shows up as edited under svn diff. No need for you to recheck
the yesterday's bunch of edits in detail - nothing of conceptual interest
happened, just it is easier for me too proofread if paper printout
pages have file and date info on them, that was the main change.

if we modify siminos-thesis.cls, we'll modify
	    \documentclass{gatech-thesis} to \documentclass{siminos-thesis} in
\input ../inputs/setupThesis    %% logical chores (nothing to edit).
Currently siminos-thesis.cls is the same as gatech-thesis.cls.

						Vaggelis Jun 26 2008
svn diff examples.tex -r 903:904
Index: examples.tex
===================================================================
--- examples.tex        (revision 903)
+++ examples.tex        (revision 904)
@@ -1,4 +1,7 @@
-% Chapter Introduction, section Example dynamical systems used throughout the thesis. Label: s:exampleIntro
+\renewcommand{\inputfile}{\version\ - edited 2007-03-11 examples}
+% Chapter Introduction, section Example dynamical systems used throughout the thesis.
+% $Author$ $Date$
+% Label: s:exampleIntro

 In this section we briefly introduce some dynamical systems that will be used as simple examples
 to demonstrate various concepts in later chapters.
@@ -17,7 +20,7 @@

 \subsection{Complex Lorenz equations}

-The complex Lorenz equations were introduced by Gibbon and McGuinness, \refref{GibMcCLE82}, as a low dimensional model
+The complex Lorenz equations were introduced by Gibbon and McGuinness, \refref{GibMcCLE82}, as a low dimensional model
 of baroclinic instability in the atmosphere. As the name suggest the equations turned out to be a complex generalization
 of Lorenz equations:
 \beq
@@ -47,6 +50,5 @@
 \eeq
 where $z_1,z_2\in \mathbf{C}$ and $\mu_j$ and $e_{jk}$ real parameters. Details on the motivation
 of those equation will be given in a later chapter. %This system corresponds to the
-%first few terms in the center manifold reduction of a $O(2)$-symmetric partial differential
+%first few terms in the center manifold reduction of a $O(2)$-symmetric partial differential
 %equation near a codimension two bifurcation.
-

Property changes on: examples.tex
___________________________________________________________________
Name: svn:keywords
   + Author Date


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





