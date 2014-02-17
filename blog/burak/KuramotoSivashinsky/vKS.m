function vel = vKS(t, x)
    %"""
    %Velocity function in the full state space.
	
    %Arguments:
    %x : vector of mode values
		%v = [Re[a1] Im[a1] Re[a2] Im[a2] ... Re[aN-1], Im[aN-1]]
    %N : number of non zero modes to be used
    %"""
	
    L = 22;
    Ltilde = 22/pi;
    N = length(x)/2; %num of nonzero modes
    #From state space variables to complex modes:
    a = zeros(N,1);
    a(:) = x(1:2:2*N-1) .+ x(1:2:2*N)*j;

    #Velocity function:
    adot = zeros(N,1);
    #adot[0] = 0, range starts from k = 1
    for k=1:N
	
	qk = k/Ltilde;
	
	nonlinearterm = 0;
    	
	for m = -N:N
	    if m < 0 
		am = conj(a(-m));
	    elseif m == 0
		am = 0;
	    else
	    	am = a(m);
	    end
	    
	    if (k-m) < 0 && -N <= (k-m)
		akMinm = conj(a(m-k));
				
	    elseif (k-m) > 0 && (k-m) <= N
		akMinm = a(k-m);
	    else
		akMinm = 0;
	    end
	    
	    nonlinearterm += am*akMinm;
	end
	 
	adot(k) = (qk**2-qk**4)*a(k) - j*(qk/2)*nonlinearterm;    
    
    end
    
    vel = zeros(N*2,1);
    vel(1:2:2*N-1) = real(adot);
    vel(2:2:2*N) = imag(adot);
