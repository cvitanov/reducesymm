(* ::Package:: *)

        (*      for Cvitanovic, Hansen Vattay paper 	23 aug 95 
	from mathematica: <<Fig_CmplxContour.m generates the eps file
        plots contours in the complex plane for
        quadratic mapping x^2-c, as well as some fix points
		*)

(*    might want to reset this:  $DefaultFont = {"Helevetica",7};*)

Clear[F,X,Contou]

R = 2.5
c = 0.5 - 0.7 I

F/: F[x_Complex,c_Complex]:= x^2 - c

X/: X[R_Real,phi_Real]:= R Sin[phi] + R Cos[phi] I

Contou/: Contou[R_Real,phi_Real]:= R Sin[phi] + R (0.1 Cos[5 phi+.1]+0.2 Cos[2 phi+2.1]+0.9) Cos[phi] I

(* Contour/: Contour[]:= (         *)
      Ccontour1 = ParametricPlot[
         	{
	  	{Re[Contou[R,phi]], Im[Contou[R,phi]]}, 
          	{Re[F[Contou[R,phi],c]], Im[F[Contou[R,phi],c]]} 
         	},{phi,0, 2 Pi},
	  			AspectRatio -> Automatic,
	  			Frame -> True,
	  			DisplayFunction -> Identity,
	  			PlotStyle -> {AbsoluteThickness[2]}
           	  ];
      FixPs = ReplaceAll[x, NSolve[((x^2 - c)^2 -c)^2 -c -x ==0,x] ];
      PerPoints1 = ListPlot[ 
		  Table[{Re[ FixPs][[i]],Im[ FixPs][[i]]},{i,Length[ FixPs]}]
          		,
	  		PlotStyle -> { AbsolutePointSize[5]},
	  		DisplayFunction -> Identity
                   ];
      ph = -.5;
      Arr1 ={Re[Contou[R,ph]], Im[Contou[R,ph]]}; 
      Arr2 ={Re[F[Contou[R,ph],c]], Im[F[Contou[R,ph],c]]}; 
      PerPoints2 = Graphics[ 
	  			{AbsoluteThickness[2],
		       Table[Circle[{Re[ FixPs][[i]],Im[ FixPs][[i]]},0.2]
                       		,{i,Length[ FixPs]}] ,
	  		         AbsolutePointSize[4],
		       Table[Point[{Re[ FixPs][[i]],Im[ FixPs][[i]]}]
                       		,{i,Length[ FixPs]}],
		       Line[{Arr1 +{0.4,0.05},Arr1,Arr1 +{0.3,0.3}}], 
		       Line[{Arr2 +{-0.4,0.1},Arr2,Arr2 +{-0.25,0.35}}],
                                }
                   ];
      Labeling = Graphics[{Text["C",{1.3,2}],
			   Text["f(C)",{-1.7,-4.3}]
                          },
                 ];
      Show[Ccontour1,  PerPoints2, Labeling,
	  		DisplayFunction -> $DisplayFunction 
      ]
      Export["trash.eps", Show[Ccontour1, PerPoints2,  Labeling]
      ];



