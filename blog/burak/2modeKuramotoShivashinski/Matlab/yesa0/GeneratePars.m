function void = GeneratePars(Ltilda, e2, a0)

q = @(i) i/Ltilda;

mu1 = q(1)^2 - q(1)^4
q1 = q(1)
mu2 = q(2)^2 - q(2)^4
q2= q(2)
e2 = e2
a0 = a0

pars(1) = mu1;
pars(2) = q1;
pars(3) = mu2;
pars(4) = q2;
pars(5) = e2;
pars(6) = a0;

save pars.mat pars;
