%clear all
%clc

ti=45.21335;
tf=48.40836;

Jac = Jacobian(ti, tf);
FloquetMultipliers = eig(Jac)
