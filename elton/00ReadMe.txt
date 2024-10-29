https://github.com/cvitanov/reducesymm folder [elton]         

--------------------------------------------------------------------------
moved to GitHub                                                 2024-10-29
Do not edit CNS subversion repo   elton       any longer:
svn $Author: predrag $ $Date: 2018-12-07 13:04:16 -0500 (Fri, 07 Dec 2018) $
--------------------------------------------------------------------------

John R Elton            
 - fluids blog/blog.tex Fall 2007, Spring, Summer 2008, Fall 2024
 - elton/blog/EltonBlog.tex
 - GaTech Spring 2007 PHYS-4699 undergraduate project

Mohammad Farazmand
 - Fall 2014 - ... LCS + mixing project

Adam Fox
 - Fall 2013 tori determination project

--------------------------------------------------------------------------
TO DO
  [ ] change \bibliography{../bibtex/elton} to reducesym *.bib  2024-10-29

--------------------------------------------------------------------------

blog/
    all matters pertaining to
        mixing, invariant tori, adjoint descent for 2D Kolmogorov,
        Elton lands
	
bibtex/elton.bib
	one bibliography for all projects

inputs/
	one directory for all LaTeX macros

figs/
    one directory for all figures
    When creating a new version, please use the SAME name (so
    it propagates to publications etc. as well, without messing
    around with renaming it)

Fig/
    copied from ChaosBook.org (do not edit here)

figsSrc/
	one directory for source programs for all figures
	(or .txt files explain how to regenerate them, if needed
     from the other repositories, for example pipes/)

talks/
    all mixed up talks and posters

FoxCvi14/
	Adam & Predrag's KS torus paper

--------------------------------------------------------------------------

NOTES:

 						Predrag Oct 29 2013
added
    torusmake   torusmake2
torusorb.f input data sets:                             STATUS:

c      open(2,file='antorb4b.dat',status='old')         DO NOT HAVE
      open(3,file='test4.dat',status='old')             HAVE
c      open(3,file='torusks412a.dat',status='old')      HAVE
      open(4,file='torusks412b.dat',status='unknown')   EMPTY FILE
c      open(5,file='torussmF.dat',status='unknown')     HAVE
c      open(3,file='torustest1.dat',status='unknown')   HAVE
        call torusfast(omg,x,xx,jb,y,ct2,dh,t4,t1)      HAVE
        call antprop(L,d,nu,b,dydx,yscal,yerr,...       HAVE
        call antpropj(L,d,nu,a,dydx,yscal,yerr,...      HAVE
        call four1(n*m,isign,x2)                        HAVE

it also calls Numerical Recipes routines (can install if needed)
        numrec_f/recipes/dftint.f  , included in dftcor.f


 						Predrag Oct 16 2013
 elton/y-lan/papers/kuramot/
    are Yuheng Lan's fortran and data files for KS tori, copied from
26637 2005-03-05 06:10 /export/home/y-lan/papers/kuramot/
    torusks412a.dat is used by torusorb.f , the figures are in:
 y-lan/tori/torusks412a1.pdf, torusks412a2.pdf
    are torus data on  the Poincar\'e section
	the other data file is empty:
    0 2005-03-05 14:26 /export/home/y-lan/papers/kuramot/torusks412b.dat

 						Predrag Dec 19 2007
to refresh homepages (after updating your svn repository), log in zero
	cd public_html
	svn up
images/ 		contains webpage figs
figs/ 			contains projects' figs


                                                Predrag Jan 17 2007
Maintain only ONE bib file for all fluid dynamics related:
	elton.bib 	(initially a copy of Gibson's fluid.bib)

                                                Predrag Aug 22 2006
Time stamp:
-----------
To have svn time-stamp file "someFile.type", include the contents of
        svnHeader.txt
into the file, and then
        svn propset svn:keywords "Date Author" someFile.type

--------------------------------------------------------------------------
HISTORY:

                                                Predrag Jan 17 2007
installed MikTeX, WinEdt, ghostscript on John's notebook

                                                Predrag Jan 13 2007
created project LaTeX template
