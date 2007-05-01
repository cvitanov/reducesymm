siminos/figSrc
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
	program plotEigenvalues.nb
	Converted to .fig with pstoedit
	

Processing
----------



Fix these:
--------------------

