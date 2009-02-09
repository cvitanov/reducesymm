siminos/figSrc/00ReadMe.txt
% $Author$ $Date$

Siminos source programs for figures
===================================

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
Added:
  testing/flows/lorenzPolar.nb
Modified:
  testing/flows/lorenz.nb
Log: Plotting Lorenz figures.

-----------------------------------------------------------------
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

Then use program bmeps ( download from http://bmeps.sourceforge.net/, be sure to install
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
