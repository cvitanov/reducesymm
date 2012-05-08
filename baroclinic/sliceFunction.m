function [F,Fp]=sliceFunction(a1,a1p,b1,b1p,rf,nx,ny) 
        aux1=0;
        aux2=0;
        aux=0;
     %Construct ap+Tap for rf_i
        for n=1:2
            if n==1
                u=a1;
                up=a1p;
            elseif n==2
                u=b1;
                up=b1p;
            end
            for k=0:nx-1
                k1=k+1;
                for l=0:ny-1
                    l1=l+1;
                    aux=k*conj(u(l1,k1))*up(l1,k1)*exp(-2*pi*1i*rf*k/nx);
                    aux1=aux1+aux;
%                     aux2=aux2+k*aux;
                end
            end
        end       
        F1=-(2*pi*1i/nx)*aux1;
%         Fp=real(-(4*pi^2/nx^2)*aux2);
        F=real(F1);