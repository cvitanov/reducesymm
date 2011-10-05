siminos/chao/matlab/ruslan/00ReadMe.txt
$Author: ruslan $
$Date: Wed Oct  5 13:33:43 GMTDT 2011 $
---------------------------------------------------------------

Ruslan 2011-03-29
My KS related Matlab code for Chao's enjoyment (or pain)...
	chao/matlab/ruslan/

ksdupo.m is the main file.  It's a mess, since I use Matlab cells, but it contains the track record of what I've done or tried to do with PPOs and RPOs in KS.

ks22f90h25t100.mat file contains data on all the PPOs and RPOs I found with T < 100.  I also have a file ks22f90h25.mat with all the 60k orbits, but it's 50MB large, so I did not commit it to the repository.  If you want I can do it, although I'm sure you'll have enough on your hands with the first 480 orbits. :-)

If you load ks22f90h25t100.mat, you'll see L = 22 and two structure variables: ppo and rpo:

>> ppo
ppo = 
1x240 struct array with fields:
    a  
    T
    r
    e
    nhit
    a1
    T1
    r1
    e1

>> rpo
rpo = 
1x239 struct array with fields:
    a
    T
    s
    r
    e
    nhit
    a1
    T1
    s1
    r1
    e1

Here what the fields contain:

>> rpo(1)
ans = 
       a: [30x1 double]  - initial condition of RPO with 16 modes and h = 0.25
       T: 16.316         - period
       s: 2.8638         - shift (only for RPO)
       r: 2.6419e-015    - residual after convergence
       e: [30x1 double]  - eigenvalues
    nhit: 95             - number of times this orbit was detected by my f90 code that searched for PPOs and RPOs.
      a1: [62x1 double]  \
      T1: 16.315          \
      s1: 2.8634           - same as above, but with 32 modes and h = 0.1
      r1: 3.0708e-014     /
      e1: [62x1 double]  /

About some of the subroutines that are used in ksdupo:

ksfmetd2.m  - KS ODE solver (you can see how it's used in ksdupo.m)

ksfmms5.m - KS ODE solver + shift/reflection + multiple shooting

ksfmms3.m - same as above but simpler: for generating return map, stability matrix (df) etc.  At the bottom of ksdupo.m file you'll see how I used it to check that indeed the orbits that Evangelos asked me about converge correctly.  You can use ksfmetd2 to actually plot the orbits, as I do in the section of ksdupo.m called

%% Plot individual RPOs and PPOs from ks22f90 with the first 4 FMs in polar coordinates (19-Jan-2010)

Evangelos, I don't know why numbers in the *.dat files don't work, but they are the older versions, so I'd stop using them.


