(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];


v[tt_]:={m1 z1r[tt]+a1 z1r[tt] (z1r[tt]^2+z1i[tt]^2)+b1 z1r[tt](z2r[tt]^2+z2i[tt]^2)+c1( z1r[tt] z2r[tt]+ z1i[tt] z2i[tt]),
m1 z1i[tt]+a1 z1i[tt] (z1r[tt]^2+z1i[tt]^2)+b1 z1i[tt](z2r[tt]^2+z2i[tt]^2)+c1( z1r[tt] z2i[tt]- z1i[tt] z2r[tt]),
m2 z2r[tt]+ e2 z2i[tt]+a2 z2r[tt](z1r[tt]^2+z1i[tt]^2)+b2 z2r[tt](z2r[tt]^2+z2i[tt]^2) +c2( z1r[tt]^2-z1i[tt]^2),
m2 z2i[tt]- e2 z2r[tt]+a2 z2i[tt](z1r[tt]^2+z1i[tt]^2)+b2 z2i[tt](z2r[tt]^2+z2i[tt]^2) +2 c2 z1r[tt]z1i[tt]}


vars[tt_]:={z1r[tt],z1i[tt],z2r[tt],z2i[tt]};


varsde={z1r,z1i,z2r, z2i};


d=Length[vars[t]];


solParIC[m1_,m2_,c1_,c2_,a1_,b1_,a2_,b2_,e2_,ic_List,tmin_,tmax_]:=Module[{eqs,sol,d},
d=Length[vars[t]];
eqs=Table[D[vars[t][[i]],t]==v[t][[i]],{i,1,d}];
sol[t1_,t2_]:=NDSolve[Flatten[{eqs,Table[vars[0][[i]]==ic[[i]],{i,1,d}] }],varsde,{t,t1,t2}];
varsde/.sol[tmin,tmax]//Flatten
]


sol2points[sol_]:=Module[{solt},
solt=InterpolatingFunctionCoordinates[sol[[1]] ]//First;
Table[sol[[i]][solt],{i,1,Length[sol]}]//Transpose
]


soltimes[sol_]:=Module[{solt},
solt=InterpolatingFunctionCoordinates[sol[[1]] ]//First
]


rot[th_]:={{Cos[th], Sin[th],0,0},{-Sin[th], Cos[th],0,0},{0,0,Cos[2th], Sin[2th]},{0,0,-Sin[2th],Cos[2th]}};


rot[th_]:={{Cos[th], Sin[th],0,0},{-Sin[th], Cos[th],0,0},{0,0,Cos[2th], Sin[2th]},{0,0,-Sin[2th],Cos[2th]}};


rotDer[th_]=D[rot[th],th];


rotGen=rotDer[0];


(*Matrix of variations*)


A[t_]= Table[D[v[t][[i]],vars[t][[j]] ],{i,1,Length[v[t]]},{j,1,Length[vars[t]]}];


(*Check for symmetry using the generator*)


checkSym[rotGen_]:=rotGen.v[t]-A[t].rotGen.vars[t]


(*Slice condition*)


sliceCond[x_List,templ_List,th_]:=(x.rot[-th]).(rotGen.templ)


sliceCond[x_List,templ_List]:=x.(rotGen.templ)


(*Rename this. Border condition*)


posCond[x_List,templ_List]:=templ.rotGen.rotGen.x;


(*Rotate single point to a chart fixed by a template*)


rot2chart[x_List,templ_List]:=Module[{thSol,xRot,nrm},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
Print[thSol];
xRot=(rot[th].x)/.thSol;
nrm=Map[Norm[#-templ]&,xRot,1];
xRot[[ Position[nrm,Min[nrm]]//First ]]//First
]


rot2chart2min[x_List,templ_List]:=Module[{thSol,xRot,nrm},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
xRot=(rot[th].x)/.thSol;
nrm=Map[Norm[#-templ]&,xRot,1];
If[Length[xRot]==4,xRot=Drop[xRot,Position[nrm,Min[nrm]]//First]];
nrm=Map[Norm[#-templ]&,xRot,1];
xRot[[ Position[nrm,Min[nrm]]//First ]]//First
]


(*Rotate single point to a chart fixed by a single Fourier mode template Subscript[y, fm]=0, Subscript[x, fm]>0 (provide index fm)*)


rot2chartExpl[x_List,fm_]:=Module[{th},
th=ArcTan[x[[2*fm-1]],x[[2*fm]]]/fm;
xRot=rot[th].x
]


rot2chartMax[x_List,templ_List]:=Module[{thSol,xRot,nrm},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
xRot=(rot[th].x)/.thSol;
nrm=Map[Norm[#-templ]&,xRot,1];
xRot[[ Position[nrm,Max[nrm]]//First ]]//First
]


(*Find angle that rotates single point to chart fixed by template*)


thetaChart[x_List,templ_List]:=Module[{thSol,xRot,nrm,minPos,theta},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
theta=th/.thSol;
xRot=(rot[th].x)/.thSol;
nrm=Map[Norm[#-templ]&,xRot,1];
minPos=Position[nrm,Min[nrm]]//First;
theta[[minPos]]//First
]


thetaChartAll[x_List,templ_List]:=Module[{thSol,xRot,nrm,minPos,theta},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
theta=th/.thSol
]


(*Find angle that rotates single point to chart fixed by Fourier mode fm*)


thetaChartExpl[x_List,fm_]:=ArcTan[x[[2*fm-1]],x[[2*fm]]]/fm;


thetaChart2[x_List,templ_List]:=Module[{thSol,xRot,theta,side,sidePlus,posPos},
thSol=Solve[sliceCond[x,templ,th]==0,th,Reals]/.C[1]->0;
theta=th/.thSol;
xRot=(rot[th].x)/.thSol;
side=templ.rotGen.rotGen.xRot;
sidePlus=Map[#>0&,side,1];
(*nrm=Map[Norm[#-templ]&,xRot,1];*)
posPos=Position[sidePlus,True];
theta[[posPos]]
]


rot2chartAngl[x_List,angl_,templ_List]:=rot[angl].x;


rot2chartAnglAll[x_List,angl_List,templ_List]:=Table[rot[angl][[i]].x,{i,1,Length[angl]}];


rot2slice[x_List,templ_List]:=Module[{xRot},
xRot=rot2chart[x,templ[[1]] ];
If[Sign[sliceCond[xRot,templ[[2]]] ]==Sign[sliceCond[templ[[1]],templ[[2]]]],{xRot,1},{rot2chart[x,templ[[2]]],2} ](*Check ridge crossing*)
]


templRelOrient[templ1_,templ2_]:=Re[Conjugate[templ1].templ2/(Norm[Conjugate[templ1]] Norm[templ2])]


rot2chartTraj[traj_List,templ_List]:=Map[rot2chart[#,templ]&,traj,1]


rot2chartTraj2min[traj_List,templ_List]:=Map[rot2chart2min[#,templ]&,traj,1]


rot2chartTrajMax[traj_List,templ_List]:=Map[rot2chartMax[#,templ]&,traj,1]


rot2chartTrajExpl[traj_List,fm_]:=Map[rot2chartExpl[#,fm]&,traj,1]


rot2sliceTraj[traj_List,templ_List]:=Map[rot2slice[#,templ]&,traj,1]


thetaChartTraj[traj_List,templ_List]:=Map[thetaChart[#,templ]&,traj,1]


thetaChartTrajExpl[traj_List,fm_]:=Map[thetaChartExpl[#,fm]&,traj,1]


thetaChartTrajAll[traj_List,templ_List]:=Map[thetaChartAll[#,templ]&,traj,1]


rot2chartTrajAngl[traj_List,angl_List,templ_List]:=Module[{out},
out=traj;
Do[
out[[i]]=rot2chartAngl[traj[[i]],angl[[i]],templ]
,{i,1,Length[angl]}
];
out
]


rot2chartTrajAnglAll[traj_List,angl_List,templ_List]:=Module[{out},
out=traj;
Do[
out[[i]]=rot2chartAnglAll[traj[[i]],angl[[i]],templ]
,{i,1,Length[angl]}
];
out
]


rot2chartTrajCont[traj_List,templ_List]:=Module[{trajRot,thSol,xRot,nrm},
trajRot=traj;
trajRot[[1]]=rot2chart[traj[[1]],templ];
Do[
thSol=Solve[sliceCond[traj[[i]],templ,th]==0,th,Reals]/.C[1]->0;
xRot=(rot[th].traj[[i]])/.thSol;
nrm=Map[Norm[#-trajRot[[i-1]] ]&,xRot,1];
trajRot[[i]]=xRot[[ Position[nrm,Min[nrm]]//First ]]//First;
,{i,2,Length[traj]}
];
trajRot
]


rot2chartTrajContBord[traj_List,templ_List]:=Module[{trajRot,thetaRot,thSol,xRot,nrm,bordTab,bordF},
trajRot=traj;
thetaRot=Table[Null,{Length[traj]}];
trajRot[[1]]=rot2chart[traj[[1]],templ];
thetaRot[[1]]=thetaChart[traj[[1]],templ];
trajRot[[2]]=rot2chart[traj[[2]],templ];(*Does not check cont. Fix this*)
thetaRot[[2]]=thetaChart[traj[[2]],templ];
Do[
thSol=Solve[sliceCond[traj[[i]],templ,th]==0,th,Reals]/.C[1]->0;
xRot=(rot[th].traj[[i]])/.thSol;
bordTab=Table[posCond[trajRot[[j]],templ],{j,i-2,i-1}];
bordF=Interpolation[bordTab];
(*nrm=Map[Norm[posCond[#,templ]-posCond[trajRot[[i-1]],templ] ]&,xRot,1];*)
nrm=Map[Norm[posCond[#,templ]-bordF[3] ]&,xRot,1];
trajRot[[i]]=xRot[[ Position[nrm,Min[nrm]]//First ]]//First;thetaRot[[i]]=th/.thSol[[Position[nrm,Min[nrm]]//First]]//First;
,{i,3,Length[traj]}
];
{trajRot,thetaRot}
]


Wop[th_]:=th-2*Pi*Round[th/(2Pi)]


itoh[psi_List]:=Module[{phi},
phi=psi;
Do[
phi[[m]]=phi[[1]];
Do[
phi[[m]]=phi[[m]]+Wop[psi[[n]]-psi[[n-1]]]
,{n,2,m}
]
,{m,2,Length[psi]}
];
phi
]


(*Integrate equations on slice fixed by templ.*)


solParICslice[m1_,m2_,c1_,c2_,a1_,b1_,a2_,b2_,e2_,ic_List,templ_List,tmin_,tmax_]:=Module[{eqs,sol,d, icsl},
icsl=rot2chart[ic,templ];
Print[icsl];
Print[icsl.(rotGen.templ)];
d=Length[vars[t]];
vtp[tt_]:=v[tt].(rotGen.templ);
ttp [tt_]:=(rotGen.{z1r[tt],z1i[tt],z2r[tt],z2i[tt]}).(rotGen.templ);
vt[tt_]:= v[tt]-(rotGen.{z1r[tt],z1i[tt],z2r[tt],z2i[tt]})vtp[tt]/ttp[tt];
eqs=Table[D[vars[t][[i]],t]==vt[t][[i]],{i,1,d}];
sol[t1_,t2_]:=NDSolve[Flatten[{eqs,Table[vars[0][[i]]==icsl[[i]],{i,1,d}] }],varsde,{t,t1,t2} ,Method->{"FixedStep",Method->"ExplicitRungeKutta"}, StartingStepSize->0.001,MaxSteps-> 1000000];
varsde/.sol[tmin,tmax]//Flatten
]


(*Integrate equations on slice fixed by templ, using rescalled time (general template)*)


solParICsliceTau[m1_,m2_,c1_,c2_,a1_,b1_,a2_,b2_,e2_,ic_List,templ_List,tmin_,tmax_]:=Module[{eqs,vars,varsde,sol,v,d, icsl},
Print[ic];
icsl=rot2chart[ic,templ];
Print[icsl];
Print[icsl.(rotGen.templ)];
d=Length[vars[t]];
vtp[tt_]:=v[tt].(rotGen.templ);
ttp [tt_]:=(rotGen.{z1r[tt],z1i[tt],z2r[tt],z2i[tt]}).(rotGen.templ);
vt[tt_]:=ttp[tt]* v[tt]-Sign[ttp[tt]](rotGen.{z1r[tt],z1i[tt],z2r[tt],z2i[tt]})vtp[tt];
eqs=Table[D[vars[t][[i]],t]==vt[t][[i]],{i,1,d}];
sol[t1_,t2_]:=NDSolve[Flatten[{eqs,Table[vars[0][[i]]==icsl[[i]],{i,1,d}] }],varsde,{t,t1,t2} (*,Method\[Rule]{"FixedStep",Method->"ExplicitRungeKutta"}*), StartingStepSize->0.001,MaxSteps-> 100000];
varsde/.sol[tmin,tmax]//Flatten
]


param3D[traj_List,time_List,proj_List]:=Module[{},
c1t=Interpolation[{time,traj[[All,proj[[1]]]]}//Transpose];
c2t=Interpolation[{time,traj[[All,proj[[2]]]]}//Transpose];
c3t=Interpolation[{time,traj[[All,proj[[3]]]]}//Transpose];
{c1t,c2t,c3t}
]


ip3D[traj_List,time_List]:=Module[{b1=traj[[All,1]]//Re,c1=traj[[All,1]]//Im,b2=traj[[All,2]]//Re,c2=traj[[All,2]]//Im},
c1t=Interpolation[{time,b1^2+c1^2}//Transpose];
c2t=Interpolation[{time,b2(b1^2-c1^2)+2b1 c1 c2}//Transpose];
c3t=Interpolation[{time,-2b1 b2 c1 +(b1^2-c1^2)c2}//Transpose];
{c1t,c2t,c3t}
]


plot3D[traj_List,proj_,opts_]:=Show[Graphics3D[Line[traj[[All,proj]]] ] ,opts]


plot2D[traj_List,proj_,opts_]:=Show[Graphics[Line[traj[[All,proj]]] ] ,opts]
