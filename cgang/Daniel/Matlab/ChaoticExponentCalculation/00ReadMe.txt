Files in this folder can be used to calculate the Lyapunov exponents for 
the Portner-Knobloch system and to optimize parameters for chaos (or whatever)
based on those exponents. This code was developed by Daniel Borrero between 8/20/2013 and 8/22/2013
and is based on the code by V. Govorukhin available at
http://www.mathworks.com/matlabcentral/fileexchange/4628-calculation-lyapunov-exponents-for-ode

lyapunov.m - Employs variational method of Wolf et al. Physica D Vol. 16, pp. 285-317, 1985
to calculate Lyapunov exponents

lyapunov_exp.mat - stores the Lyapunov exponents for the set of parameters saved in params.mat

params.mat - Whenever run_lyap.m decides that a set of parameters is more "chaotic" than what
was being used before, the new set is saved here

PKeqInv - Defines the extended ODE system for the PK system and sets up the variational equation
that lyapunov.m uses.

run_lyap.m - Iterates lyapunov.m with slight modifications to the parameters saved in params.mat
If Lyapunov exponents calculated using the modified parameters satisfy a user defined criterion of
"better" (e.g., leading exponent is larger with the modified set of parameters than with the original 
set), it saves the new set of parameters to params.mat and the Lyapunov exponents to lyapunov_exp.mat.