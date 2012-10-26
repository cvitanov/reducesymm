(* ::Package:: *)

b = 0.6
d=0.001
eps = 0.000000000001
nmax=8
z = 0.404566786121645
lam = 1./(z*(1-z)*(1-b*z))
f[x] = lam*x*(1-b*x)*(1-x)
border[3] = z
border[1] =  1.
border[2] = 0.
prev=4
For[n=1,n<nmax,n++,
        ctr=3;
        For[l=3,l<prev,l++,
                dumb = x/.FindRoot[border[l]-f[x]==0,{x,0.}];
                If[dumb>=border[2],
                	ctr=ctr+1;
                	border[ctr] = dumb];
                dumber =  x/.FindRoot[border[l]-f[x]==0,{x,1.}];
                If[dumber<=border[1],
                	ctr=ctr+1;
                	border[ctr]= dumber]];
	Print[ctr];
	v = Array[border,ctr];
	prev=ctr+1;	
	pts = Sort[v];
	For[h=1,h<ctr+1,h++,
		a[h]=Extract[pts,h]];	
	For[i=1,i<ctr,i++,
                dt[i]=a[i+1]-a[i]];
        For[j=1,j<ctr,j++,
                For[i=1,i<ctr,i++,
                        elem[i,j]=1/(Sqrt[4*Pi*d]*dt[i])*NIntegrate[Exp[-(y-f[x])^2/(4*d)],{y,a[j],a[j+1]},{x,a[i],a[i+1]},Method-> GaussKronrod];
			If[elem[i,j]<eps, 
				elem[i,j]=0]]];
	For[i=1,i<ctr,i++,
		mass = elem[i,1];
		For[j=2,j<ctr,j++,
			mass = Max[elem[i,j],mass]];
		mxm[i] = mass];		
	For[i=1,i<ctr,i++,	
		For[j=1,j<ctr,j++,
				If[elem[i,j]==0,		
				rett[i,j]={GrayLevel[1.],Rectangle[{a[i],a[j]},{a[i+1],a[j+1]}]},
				rett[i,j]={GrayLevel[(1.-elem[i,j]/mxm[i])],Rectangle[{a[i],a[j]},{a[i+1],a[j+1]}]}]]];
	A = Array[elem,{ctr-1,ctr-1}];
	ev = Eigenvalues[A];
	dagev = Eigenvectors[Transpose[A]];
	evctrs = Eigenvectors[A];
	area=0.;
	For[i=1,i<ctr,i++,
		freq[i] = Abs[Extract[dagev,{1,i}]];
		area= area + Abs[freq[i]]];
	Print[area];
	liap=0.;
	For[i=1,i<ctr,i++,
		ist[i]= {Hue[0.7],Rectangle[{a[i],0.},{a[i+1],Abs[freq[i]]/(area*dt[i])}]};
		xmed=0.5*(a[i]+a[i+1]);
		liap = liap + Abs[freq[i]]*Log[Abs[lam*(1-2*b*xmed-2*xmed+3*b*xmed^2)]]];
	liap = liap/area;
	density = Array[ist,ctr-1];
	Print[Show[Graphics[density],AspectRatio->1,Axes->Automatic]];
	esc[n] =  Extract[ev,1];
	logesc[n] = -Log[esc[n]]; 
	conv[n] = Extract[ev,2];
	logconv[n] = -Log[Abs[conv[n]]];
	Print[logesc[n]];
	Print[liap];
	R = Array[rett,{ctr-1,ctr-1}];
	Show[Graphics[R,AspectRatio->1,Axes->True,AxesOrigin->{0.,0.}]]]
For[n=1,n<nmax-1,n++,
	gamma[n]=Log[Abs[esc[nmax-1]-esc[n]]];
	delta[n]=Log[Abs[conv[nmax-1]-conv[n]]]]
logs = Array[logesc,nmax-1]
rates = Array[gamma,nmax-2]
seconds = Array[logconv,nmax-1]
convsec = Array[delta,nmax-2] 
ListPlot[logs,PlotStyle->PointSize[0.03]]
ListPlot[rates,PlotStyle->PointSize[0.03]]
ListPlot[seconds,PlotStyle->PointSize[0.03]]
ListPlot[convsec,PlotStyle->PointSize[0.03]]














