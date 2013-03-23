% siminos/figSrc/00ReadMe.txt
% $Author$ $Date$

Siminos, Froehlich, ... source programs for figures
===================================================

TO ALL:
    please save, clearly document here the source code for your figures

read also
    halcrow/figSrc/00ReadMe.txt
    dasbuch/book/FigSrc/00ReadMe.txt

----------------------------------------------------------------
Converting eps to pdf
                                       Predrag 2013-03-22
robjhyndman.com/hyndsight/converting-eps-to-pdf
Now that there is a way to sync a pdf file and tex file in both
directions, the only remaining reason to use dvi files is when the
graphics are in eps format. However, that problem has also been
solved for those using Mik­TeX 2.8 or TeXLive 2009. In Mik­TeX 2.8,
simply include the package epstopdf along with graphicx. (As noted
in the comments below, even this step is not necessary in TeXLive
2009.) Then when you use pdflatex, the eps files will be
automatically converted to pdf at compile time. (The
conversion only happens the first time you process the file, and
is skipped if there is already a pdf file with the same name.) For
example:
    \documentclass{article}
    \usepackage{graphicx,epstopdf}
    \begin{document}
    \includegraphics[width=\textwidth]{fig1}
    \end{document}
Then even though the only graphics file available is
fig1.eps, this will still be processed ok using pdflatex or
pdftexify. On the first pass, a new file called
fig1-eps-coverted-to.pdf is created and inserted at the
appropriate place.

to convert *.ps to *.eps:
www.online-convert.com  (there are also many scripts)

                                        Predrag 2012-04-20
To Daniel
	siminos/figSrc/Rossler No Arrows/
	the figs look nice, but not tightly cropped (acres of white)
[x] trajectorynoarrows.png
[ ] nearequilibriumnoarrows.png
   make the same size as other 3, to facilitate inkscaping
[x] bothsectionnoarrows.png
[x] farequilibriumnoarrows.png

[ ] remebr to remove siminos/figSrc/Rossler Arrows/

----------------------------------------------------------------
                                        Sarah Flynn 2012-04-12

www.sagemath.org
    generates cool graphics in python
www.texample.net/tikz/
    can generate animations that can be embedded into pdf files
www.texample.net/tikz/examples/tag/animations/

----------------------------------------------------------------
                                        Borrero 2012-04-09

siminos/figSrc/matlab/Wurst01.m
    An updated version of TWand01.m. Draws 3D surface to represent group
    orbit of 01 rpo. This figured is referenced as fig:CLf01group in the
    atlas blog/paper. Comments in the code should be pretty self
    explanatory. Has some ancillary code saved with it:
    ComplexLorenzEOM.m, gCLE.m, GroupOrbit.m, and mArrow3.m, so don't go
    deleting those.

----------------------------------------------------------------
                                        Borrero 2012-03-23

siminos/figSrc/matlab/TWand01.m
    generates a plot of one period of the 01 relative periodic orbit for
    the 5D complex Lorenz system along with the group orbit of the
    relative equilibrium TW_0 (or at least their projection in the x1,
    x2, z subspace). Comments in the code should make playing around with
    it fairly self explanatory.
    Once figure is generated in Matlab it is up to you to figure out how
    you want to export it.

----------------------------------------------------------------
                                        Predrag 2012-01-06
you can hand-write equations
    webdemo.visionobjects.com/equation.html?locale=default
get LaTeX back!

----------------------------------------------------------------
                                        Predrag 2012-02-06
a 2010 overview over LaTeX graphics:
    www.math.nus.edu.sg/aslaksen/cs/textrix.pdf
----------------------------------------------------------------
                                        Predrag 2012-01-04
MapleWorks.com/software-products/overview-software-products/android-note-taking-app/
    MaplePaint seems easy to use from Samsung Galaxy Android tablet,
    creates compact svg vector graphics, which inkscape converts
    to pdf

----------------------------------------------------------------
                                        Predrag 2010-12-21

*.pdf vector graphic files are much larger than
*.png raster graphics, so deferred these Mathematica generated
files to the journal submission:

froehlich/CNSNS-v1/dthetanearsing.pdf
		   Fullspace.pdf
		   singpass1.pdf
		   dthetasing.pdf
		   RedTrajNoPlane1.pdf

----------------------------------------------------------------
                                        Stefan  2010-12-20

'final' figures for the ../froehlich/slice article

../figs/Fullspace.png
RedTrajNoPlane1.png
RedTrajNoPlane2.png
RedTrajPlane1.png
RedTrajPlane2.png
RedTrajPlane3.png
dthetanearsing.png
dthetasing.png
singpass1.png

----------------------------------------------------------------
                                        Evangelos 2011-01-06

  repository vaggelis/sys_utils/epstopdf_dir
			      sys_utils/makefile

Added script that applies epstopdf on all eps files
in a specified directory. Usually results in huge savings in
figure size, without loss of quality. Use it for arxiv (with
pdflatex), keep eps files for APS submission.

----------------------------------------------------------------
                                        Evangelos 2011-03-15

    I've recently started using matplotlib for visualization:
matplotlib.sourceforge.net
    It is a matlab-like python library, that can be used either
    interactively or through scripts. Excellent! Perhaps it only
    lags behind matlab and mathematica in 3D capabilities.

    I've also noticed asymptote
asymptote.sourceforge.net
    produces high-quality figures in 2D and 3D. I think it is intended
    mostly for drawing, rather than for data visualization. I haven't
    tried it.

----------------------------------------------------------------
                                        Predrag 2010-04-02

talked to Zbigniew Nitecki,  www.tufts.edu/~znitecki who is
writing a honors calculus texbook. He does everything in
pstricks, says, use

tug.org/PSTricks/
en.wikipedia.org/wiki/PSTricks
www.ctan.org/tex-archive/help/Catalogue/entries/pst-3d.html

and for 3D with shadows:
melusine.eu.org/syracuse/pstricks/

also, get the German edition of
www.amazon.com/PSTricks-Graphics-PostScript-TeX-LaTeX/dp/1906860130/


					Vaggelis Mar 28 2009

To have latex correctly display figures with white space extending
beyond the bounding box use the option clip=true in the \includegraphics
command:

\includegraphics[width=w, clip=true]{file}

Works with the graphicx package, the same should be achieved
with use of \includegraphics* in graphics package.

					Predrag Feb  9 2009

    dasbuch/book/Fig/ lorenz2Poinc.eps
                      lorenzPolarManifDetail1
                      lorenzSaddle0
                      lorenzPolar1.eps
have never worked - as can be seen in
    siminos/figSrc/lorentzWoes.pdf
(and looking at them in the book.dvi is also instructive) they
take lots of white acreage: in Fig 3.7 (a) is partially
covered, and so is the page runner at the top; in Fig 4.7 (a),
(b) are covered, and right figure covers part of the left one;
in Fig 10.8 (a) and the page runner are covered. Can you fix
that. Also, copy/edit into your thesis all Lorentz and related
stuff from dasbuch/book/book.tex; you do not have much time to
finalize the thesis, this will help a lot. If you have
documented source programs for these, I do not see where. The
closest I come is:

Author: siminos 2008-01-10  Revision: 865
  testing/flows/lorenzPolar.nb
  testing/flows/lorenz.nb
Log: Plotting Lorenz figures.

Vaggelis fix:
------------
Use:

eps2eps -dEPSCrop -dEmbedAllFonts=true file.eps corrected_file.eps

The crucial option is -dEPSCrop, -dEmbedAllFonts=true can help if
there are problems with fonts, eg mathematica fonts. The output
has to be a different file. Not sure whether the fix is universal.


-----------------------------------------------------------------
                                        Vaggelis Oct 9 2009

Lorenz sources

source for
	siminos/figs/lorenzAttr.eps
    (lorenzAttrChaosbook.eps is Chaosbook version with EQ_i instead of E_i)
is in
	vaggelis/testing/flows/lorenzPoincFlow.nb


-----------------------------------------------------------------
					Vaggelis Oct 1 2009

source for
	siminos/rpo_ks/figs_pst/splitting.eps
used in
	figSrc/rpo_ks_figs/ks22E2-E3hetero.tex
is in
	[local slow backup]/home/vasimos/sandbox/PropIntEquilDet/
files:	
	integrNew.nb
	equilFamily.f90
	main.f90

                                        Predrag Sep 13 2009

sources for
        siminos/figs/lyapSpec.eps
        siminos/figs/lyapSpecRscld.eps
used in
        siminos/blog/blog.tex
are
        lyapSpec.gnu
        lyapSpecRscld.gnu


					Predrag Feb 10 2007

source for
   	rpo_ks/figs/ksBifDiag.eps
is in
	programs/production/ksBif/
source for
	rpo_ks/figs/ks22TurbConn_xfig.eps
is in	
	programs/production/KS22.0/rpo/
	use plotEnergyRPO.nb and convert to .fig with:
	pstoedit -f xfig -usebbfrominput -adt ks22TurbConn.eps ks22TurbConn.fig
	Edited in xfig and exported in .eps format.
source for
	rpo_ks/figs/*wKS22equil.eps
is in
	programs/production/KS22.0/equil/
	file: plotProfile.nb
	use pstoedit to convert to .fig and then export to .eps
source for
	rpo_ks/figs/L22-eqvaEigenvalues.eps
is in
	programs/production/KS22.0/equil
	1st version: 	mathematica notebook plotEigenvalues.nb
			Converted to .fig with pstoedit
	2nd version:	gnuplot script plotEigenvalues.gnuplot
source for
	rpo_ks/figs/equilSpatial.eps
is in	
	programs/production/KS22.0/equil
	file: plotSpatial.nb
	

Processing
==========


figure too big?
===============
					Evangelos	30 Aug 2007

Use ImageMagick to convert eps to png:

 # convert +antialias file.eps file.png

Option +antialias will reduce number of colors used. Usually this greatly reduces size
without significant loss in quality over true-color.

Then use program bmeps ( download from bmeps.sourceforge.net/, be sure to install
all required programs including the developer versions before compiling) to generate
reduced size eps:

 # bmeps -leps3 -tpng file.png new_file.eps

Options -leps3 will generate level-3 eps. The reduction in size is anywhere from 5 to 50
times (!). On the other hand one trades scalability of both graphics and fonts.

A moderate reduction in size (by a factor of 2) of some eps figures , particularly those
create by programs such as Mathematica, can be achieved if one uses pstoedit to convert
them to fig files:

 # pstoedit -f xfig -usebbfrominput -adt file.eps file.fig

The files are then loaded to xfig and exported in eps. This procedure also fixes a problem
with some Mathematica fonts not being embedded in the mapping file of the figures (-adt option
converts non-recognized fonts to polygons).

To create a png with a transparent background, this works pretty well:

 # (cat file.eps ; echo showpage ) | gs -q -dEPSCrop -dBATCH -dNOPAUSE -sDEVICE=pngalpha -sOutputFile=file.png -f -

To get a level-3 eps with transparent background from this png use:

 # bmeps -leps3,c.i.m=y,i.m.t.l=0 -tpng  -a file.png new_file.eps

The resulting files do not always display correctly with gs.

Question: How can I use bitmapped graphics while keeping scalable fonts?

The following answer is Mathematica specific regarding the
creation of files but one should be able to use it with other
programs. It cannot be used to separate text and graphics in
already existing images, one has to export those separately:

Export three figures from Mathematica:
fig.eps: The original figure
figB.png: The part of the figure that can be bitmapped
figC.eps: The part of the figure that has to remain in vector format

Plot all of them with the same options (use function Absolute options
to get the info from the original image) and do not specify size or
resolution when exporting.

Wrap figB.png into an eps file:

 # bmeps -leps3 -tpng figB.png figB.eps

Use the original figure to obtain the correct bounding box for the eps files:

 # gs -sDEVICE=bbox -dNOPAUSE -dBATCH fig.eps

Replace the bounding box in figB.eps and figC.eps with the output of the above.
Use pstricks to superimpose the two figures. A sample fig.tex file should look like:

	\documentclass{article}
	\usepackage{graphicx}
	\usepackage{pstricks}
	\pagestyle{empty}

	\begin{document}


	\rput(0,0){\includegraphics{plDetGraphB.eps}}
	\rput(0,0){\includegraphics{plDetTextC.eps}}


	\end{document}

Then

 # latex fig.tex
 # dvips -E fig.dvi -o figNew.eps

If dvips complains about missing fonts one can convert the
missing fonts in figC.eps to curves using:

eps2eps figC.eps figCfixed.eps

or

pstoedit -f ps -adt -usebbfrominput figC.eps figCfixed.eps

and after editing the bounding box in either case use
figCfixed.eps in pstricks. This will generally result in even
smaller file size so it might be a good idea to use eps2eps
anyway with figC.eps.

					Mason Porter 	20 Aug 2003

 One downloads jpeg2ps

 One converts all .ps files to .jpg to make them small (so they
 aren't as sharp, but they're _much_ smaller).

 One then converts the jpgs back to .ps with jpeg2ps blah.jpg > blah.ps

 (These are jpgs with .ps 'wrappers' to make latex think they're .ps
 files.  They take up about twice the space of a jpg file, so--in
 particular--the 11 meg .ps file is a 1 meg file on the arxiv.)

--------------

					Nicolas		27 Feb 2002
used gimp on linux to generate
    174895 Feb 22 18:13 Fig/standard1.2.eps
from
   1629686 Jul  3  2001 OldFig/standard1.2-070301.eps
decreased by factor 10!

--------------

					Predrag		27 May 2007

ps2png energyBalance_pst.eps energyBlncKS.png produces coarse image
gimp  energyBlncKS.png
	save as energyBlncKS.eps	  % poor, but good for a talk
	size reduced by a factor of 10

Better;
 1382283 energyBalance_pst.eps
gimp energyBalance_pst.eps - full antialias, resolution 100
	save as
   41603 energyBlncKS.png
gimp energyBlncKS.png - save as
  133561 energyBlncKS.eps
	size reduced by a factor of 10


NOTES
=====
epstopdf --nogs connEP.eps | gs -q -dBATCH -dNOPAUSE -sDEVICE=pngalpha -sOutputFile=test.png -f -

Fix these:
==========
