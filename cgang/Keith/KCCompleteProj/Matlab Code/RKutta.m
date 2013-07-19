function [ y2 ] = RKutta(funfcn, x, dt)
    k1 = funfcn(x);
    k2 = funfcn(x+k1*dt/2);
    k3 = funfcn(x+k2*dt/2);
    k4 = funfcn(x+k3*dt);
    y2 = x +1/6*(k1+2*k2+2*k3+k4)*dt;
end