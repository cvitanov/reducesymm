(* ::Package:: *)

b = 0.6
nmax=7
eps = 0.000000000001
z = 0.4045667835969583
lam = 1./(z*(1-b*z)*(1-z))
f[x] = lam*x*(1-b*x)*(1-x)
preim[1] = z
border[1] = z
border[2] =  lam*z*(1-b*z)*(1-z)
border[3] = lam*border[2]*(1-b*border[2])*(1-border[2])
ctr= 1
k=2
cum=1
For[n=1,n<nmax,n++,	
	cum=cum+k;
	k=k*2]
For[l=1,l<cum+1,l++,
        preim[l+ctr] = x/.FindRoot[preim[l]-f[x]==0,{x,0.}];
        preim[l+ctr+1]= x/.FindRoot[preim[l]-f[x]==0,{x,1.}];	
	ctr=ctr+1]
k=2
For[n=1,n<nmax,n++,
	For[p=1,p<k+1,p++,
		border[k+p+1]=preim[k+p-1]];
	v = Array[border,2*k+1];	
	pts = Sort[v];
	Print[pts];
	For[h=1,h<2*k+2,h++,
		a[h]=Extract[pts,h]];	
	For[i=1,i<2*k+1,i++,
                dt[i]=a[i+1]-a[i]];
        For[j=1,j<2*k+1,j++,
		uj = Min[x/.FindRoot[a[j]-f[x]==0,{x,0.}],x/.FindRoot[a[j+1]-f[x]==0,{x,0.}]];
                wj = Max[x/.FindRoot[a[j]-f[x]==0,{x,0.}],x/.FindRoot[a[j+1]-f[x]==0,{x,0.}]];
                xj = Min[x/.FindRoot[a[j]-f[x]==0,{x,1.}],x/.FindRoot[a[j+1]-f[x]==0,{x,1.}]];
                yj = Max[x/.FindRoot[a[j]-f[x]==0,{x,1.}],x/.FindRoot[a[j+1]-f[x]==0,{x,1.}]];
                For[i=1,i<2*k+1,i++,
                        If[wj>=a[i] && uj<=a[i+1],
                                dmn[i,j] = (Min[wj,a[i+1]]-Max[uj,a[i]])/dt[i],
                                dmn[i,j] = 0.];
                        If[yj>=a[i] && xj<=a[i+1],
                                dmn[i,j] =  (Min[yj,a[i+1]]-Max[xj,a[i]])/dt[i]]]];
        For[i=1,i<2*k+1,i++,
                For[j=1,j<2*k+1,j++,
                                rett[i,j]={GrayLevel[1.-dmn[i,j]],Rectangle[{a[i],a[j]},{a[i+1],a[j+1]}]}]];
	B = Array[dmn,{2*k,2*k}];
	H = Transpose[B];
	ev = Eigenvalues[B];
	evctrs = Eigenvectors[B];
	lev = Eigenvalues[H];
	lvcts = Eigenvectors[H];
	area=0.;
	For[i=1,i<2*k+1,i++,
		freq[i] = Extract[lvcts,{1,i}];
		area= area + Abs[freq[i]]];
	For[i=1,i<2*k+1,i++,
		ist[i]= {Hue[0.7],Rectangle[{a[i],0.},{a[i+1],Abs[freq[i]]/(dt[i]*area)}]}];
	density = Array[ist,2*k];
	Show[Graphics[density],AspectRatio->1,Axes->Automatic];
	esc[n] =  Extract[ev,1];
	logesc[n] = -Log[esc[n]]; 
	Print[logesc[n]];
	conv[n] = Extract[ev,2];
	logconv[n] = -Log[Abs[conv[n]]];
	R = Array[rett,{2*k,2*k}];
	Show[Graphics[R,AspectRatio->1,Axes->True,AxesOrigin->{0.,0.}]];
	k=k*2]
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














