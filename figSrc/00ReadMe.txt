siminos/figSrc/00ReadMe.txt
% $Author$ $Date$

Siminos source programs for figures
===================================
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
					Mason Porter 	20 Aug 2003 

 One downloads jpeg2ps

 One converts all .ps files to .jpg to make them small (so they aren't as
 sharp, but they're _much_ smaller).

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


Fix these:
==========

