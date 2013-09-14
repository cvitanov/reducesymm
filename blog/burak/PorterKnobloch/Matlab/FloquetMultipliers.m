function Lambda = FloquetMultipliers(ti,tf)

Lambda = eig(Jacobian(ti, tf));
