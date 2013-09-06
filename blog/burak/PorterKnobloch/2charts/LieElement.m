function xp = LieElement(phi,x)
%Function applies the group transportation with the parameter theta and
%returns the transformed vector

xp = [cos(phi)  sin(phi)  0  		  0;
	 -sin(phi)  cos(phi)  0  		  0;
	  0  		0  		  cos(2*phi)  sin(2*phi);
	  0  		0 		 -sin(2*phi)  cos(2*phi)] * x;
