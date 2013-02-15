How to use 'lyap-upo.cpp'

* If you run the program without argument, it works in the interactive mode. You can specify the system, the number of time steps, what to record, etc. interactively.

* You can give these specifications as the arguments of the program, exactly in the order you would be asked in the interactive mode. If you launch the program from a shell script, you should use this "silent" mode.

* The first four #define lines are options that you can switch on if you want:

OMPSW makes the program run in a parallel mode using openMP (of course you need to pass the openmp option to the compiler as well). OMPSW=1 is optimal when you compute many Lyapunov modes, while OMPSW=2 is good if you compute only a few Lyapunov modes for a very large system.

OMPRSW speeds up the parallel program even more (to be used always with OMPSW), at the cost of reproducibility of the results. In other words, you will have different numerical errors (round-off errors) each time you run the program, even if you use the same parameter values.

ARMASW must be on if you need C++ library "Armadillo" (linear algebra). It is required if you want to compute angle between vectors / subspaces.

If FFTWSW is on, the program uses FFTW for Fourier transforms. Otherwise, the built-in code taken from the numerical recipes is used.
