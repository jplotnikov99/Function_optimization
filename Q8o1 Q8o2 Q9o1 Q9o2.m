(* ::Package:: *)

(* ::Section:: *)
(*Generate code for testing*)


(* ::Input::Initialization:: *)
singlegrid[n_,m_,k_,V_,s_, vw_,x_,T_]:=
Block[{gamw = 1/Sqrt[1-vw^2] , dfs={},Vx=1},
f0=1/(Exp[w]+s);
df0=D[f0,{w,k}];
pwt=Sqrt[w^2-x^2];
pzt=gamw*(y*pwt-w*vw);
ET=gamw*(w-vw*y*pwt);
If[V==1,Vx*=(Abs[pzt]/Sqrt[pzt^2+x^2])^2*1/Sqrt[1-pzt^2/ET^2]];
integr= -3/(Pi^2*gamw)*T^(n-m-k+1)*Vx*pwt*pzt^n/ET^(m-1)*df0//Simplify;
(*Print[integr];*)
dfs =NIntegrate[integr,{w,x,Infinity},{y,-1,1},PrecisionGoal->6,AccuracyGoal->8,MinRecursion->5,MaxRecursion->8,WorkingPrecision->8];
Return[dfs];
]
D0[x_,s_,vw_,T_]:=singlegrid[0,0,1,0,s,vw,x,T];
D1[x_,s_,vw_,T_]:=singlegrid[1,1,1,0,s,vw,x,T];
D2[x_,s_,vw_,T_]:=singlegrid[2,2,1,0,s,vw,x,T];
Q1[x_,s_,vw_,T_]:=singlegrid[0,1,2,0,s,vw,x,T];
Q2[x_,s_,vw_,T_]:=singlegrid[1,2,2,0,s,vw,x,T];
Qe1[x_,s_,vw_,T_]:=singlegrid[0,1,1,0,s,vw,x,T];
Qe2[x_,s_,vw_,T_]:=singlegrid[1,2,1,0,s,vw,x,T];
Q8o1[x_,s_,vw_,T_]:=0.5Re[singlegrid[-1,1,1,1,s,vw,x,T]];
Q8o2[x_,s_,vw_,T_]:=0.5Re[singlegrid[0,2,1,1,s,vw,x,T]];
Q9o1[x_,s_,vw_,T_]:=0.25(Re[singlegrid[-1,3,1,1,s,vw,x,T]-Re[singlegrid[-1,2,2,1,s,vw,x,T]/Sqrt[1-vw^2]]]);
Q9o2[x_,s_,vw_,T_]:=0.25(Re[singlegrid[0,4,1,1,s,vw,x,T]-singlegrid[0,3,2,1,s,vw,x,T]/Sqrt[1-vw^2]]);
N0[mass_,s_,T_]:=4Pi NIntegrate[p^2/(Exp[Sqrt[p^2+mass^2]/T]+s),{p,0,Infinity}] (*Auxiliary function*)
Rbar[mass_,s_,vw_,T_]:= Pi*(1-vw^2)/N0[mass,s,T] T Re[NIntegrate[Log[Abs[(Sqrt[w^2-(mass/T)^2]-vw w)/(Sqrt[w^2-(x)^2]+vw w)]]/( Exp[w]+s),{w,mass/T,Infinity},MaxRecursion->20]];
K0[x_,s_,vw_,T_]:=-singlegrid[0,0,0,0,s,vw,x,T];
K4[x_,s_,vw_,T_]:=singlegrid[2,2,1,0,s,0,x,T];


(* ::Input::Initialization:: *)
singlegridint[n_,m_,k_,V_,s_, vw_,x_,T_]:=
Block[{gamw = 1/Sqrt[1-vw^2] , dfs={},Vx=1},
f0=1/(Exp[w]+s);
df0=D[f0,{w,k}];
pwt=Sqrt[w^2-x^2];pzt=gamw*(y*pwt-w*vw);
ET=gamw*(w-vw*y*pwt);
If[V==1,Vx*=(Abs[pzt]/Sqrt[pzt^2+x^2])^2*1/Sqrt[1-pzt^2/ET^2]];
integr= -3/(Pi^2*gamw)*T^(n-m-k+1)*Vx*pwt*pzt^n/ET^(m-1)*df0//Simplify;
(*Print[integr];*)
dfs =integr//Simplify;
Return[dfs];
]
D0int[x_,s_,vw_,T_]:=singlegridint[0,0,1,0,s,vw,x,T];
D1int[x_,s_,vw_,T_]:=singlegridint[1,1,1,0,s,vw,x,T];
D2int[x_,s_,vw_,T_]:=singlegridint[2,2,1,0,s,vw,x,T];
Q1int[x_,s_,vw_,T_]:=singlegridint[0,1,2,0,s,vw,x,T];
Q2int[x_,s_,vw_,T_]:=singlegridint[1,2,2,0,s,vw,x,T];
Qe1int[x_,s_,vw_,T_]:=singlegridint[0,1,1,0,s,vw,x,T];
Qe2int[x_,s_,vw_,T_]:=singlegridint[1,2,1,0,s,vw,x,T];
Q8o1int[x_,s_,vw_,T_]:=1/2Re[singlegridint[-1,1,1,1,s,vw,x,T]];
Q8o2int[x_,s_,vw_,T_]:=1/2Re[singlegridint[0,2,1,1,s,vw,x,T]];
Q9o1int[x_,s_,vw_,T_]:=1/4(Re[singlegridint[-1,3,1,1,s,vw,x,T]-Re[singlegridint[-1,2,2,1,s,vw,x,T]/Sqrt[1-vw^2]]]);
Q9o2int[x_,s_,vw_,T_]:=1/4(Re[singlegridint[0,4,1,1,s,vw,x,T]-singlegridint[0,3,2,1,s,vw,x,T]/Sqrt[1-vw^2]]);


(* ::Section:: *)
(*Q8o1 for fermions*)


(* ::Subsection:: *)
(* x \[RightArrow] 0*)


(* ::Code:: *)
(*$Assumptions= x>0&&z>0&&1>y>-1&&1>vw>0&&T>0&&s\[Element]Integers*)
(*temp= (Q8o1int[x,1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}//Simplify)/.{Abs[x_]->x}//FullSimplify*)
(*Q8o1lowX[x_,vw_,T_] = (Re[Integrate[AsymptoticIntegrate[temp,{y,-1,1},{x,0,1}],{z,0,Infinity}]]//Simplify)//TrigToExp//ComplexExpand//Refine[#,{vw \[Element] Reals, 1>vw>0}]&//Simplify*)


Q8o1lowX[x_,vw_,T_] =(3 vw (-2+x) (2 \[Pi] Sqrt[1-vw^2]-2 vw Log[vw]+vw Log[1-Sqrt[1-vw^2]]+vw Log[1+Sqrt[1-vw^2]]))/(16 \[Pi]^2 T^2);


(* ::Code:: *)
(*StringReplace[ToString[Q8o1lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow","Sqrt"->"sqrt"}]*)


(* ::Code:: *)
(*plot=Plot[NIntegrate[temp/.{x->0.001,T->1},{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]*)


(* ::Code:: *)
(*Show[plot,Plot[Q8o1lowX[1,vw,1],{vw,0,1},PlotStyle->Orange]]*)


(* ::Code:: *)
(*vwall=0.4;*)
(*ErrorPlotData=Table[{10^i,Abs[Q8o1lowX[10^i,vwall,1]/Q8o1[10^i,1,vwall,1]-1]},{i,-4.5,0,0.25}]*)
(*ListLogLogPlot[ErrorPlotData]*)


(* ::Subsection:: *)
(*Sanity check*)


(* ::Code:: *)
(*Q8o1[0.001,1,0.2,1]/Q8o1lowX[0.001,0.2,1]-1*)
(*Q8o1[0.0001,1,0.5,1]/Q8o1lowX[0.0001,0.5,1]-1*)
(*Q8o1[0.00001,1,0.999,1]/Q8o1lowX[0.00001,0.999,1]-1*)


(* ::Subsection:: *)
(* x \[RightArrow] \[Infinity]*)


(* ::Code:: *)
(*$Assumptions=x>0*)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,2}];*)
(*temp2=Integrate[temp,{y,-1,1}]]*)
(**)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,3}];*)
(*temp3=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,4}];*)
(*temp4=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,5}];*)
(*temp5=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,6}];*)
(*temp6=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=temp/.Abs[x_]->x//Simplify;*)
(*temp=Asymptotic[temp,{x,Infinity,7}];*)
(*temp7=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*vwall=0.1;*)
(*dataQ8o1 = Table[{x,Q8o1[x,1,vwall,1]},{x,0.1,20,0.5}];*)


(* ::Code:: *)
(*templist={temp2,temp3,temp4,temp5,temp6,temp7};*)


(* ::Code:: *)
(*Table[Q8o1Expansion[x_,vw_,T_, sterm, t]= Integrate[(Series[templist[[t]],{s,0,sterm}]//Normal)/.{s->1},{z,0,Infinity}],{sterm,0,3},{t,templist//Length}];*)


(* ::Code:: *)
(*Show[ListPlot[dataQ8o1,PlotStyle->Green,ImageSize->800,PlotRange->{-0.002,0.001}],Plot[Evaluate[Table[Q8o1Expansion[x,vwall,1, sterm, t],{sterm,0,3},{t,{1,3}}]],{x,1,10},PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,4}],PlotRange->All]]*)


(* ::Code:: *)
(*c=0;*)
(*Flatten[Table[c+=1;{c, "s =", sterm, " | t =", t},{t,{1,2,3,4,5}},{sterm,{0,1,2,3}}],1]//TableForm*)


(* ::Code:: *)
(*plotlist=Evaluate[Flatten[Table[Table[{i[[1]],Abs[(i[[2]]/Q8o1Expansion[i[[1]],vwall,1, sterm, t])-1]},{i,dataQ8o1}],{t,templist//Length},{sterm,{0,1,2,3}}],1]];*)
(*Manipulate[ListLogPlot[plotlist, *)
(*  PlotStyle -> (Opacity@Boole@MemberQ[x, #] & /@ Range@Length@plotlist),ImageSize->1000,PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,4}],PlotRange->All], *)
(*    {{x, {13},"Approximation"}, Dynamic@Range@Length@plotlist, TogglerBar}]*)


Q8o1highX[x_,vw_,T_]=(3 E^-x vw (-1+vw^2) (6996285-7416744 x+8 (-3440640 vw^10-12288 vw^8 (-1201+40 x)-768 vw^6 (34137+16 x (-137+8 x))-32 vw^4 (-245319+8 x (2385+16 x (-45+8 x)))+32 x^2 (11241+8 x (-399+16 x (3+8 x)))+vw^2 (4813125-32 x (-20211+8 x (2073+16 x (-41+8 x)))))))/(262144 \[Pi]^(3/2) T^2 Sqrt[2-2 vw^2] x^(11/2));


(* ::Code:: *)
(*StringReplace[ToString[Q8o1highX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow","Sqrt"->"sqrt"}]*)


(* ::Code:: *)
(*Show[Plot[{Q8o1highX[x,vwall,1],Q8o1lowX[x,vwall,1]},{x,0.001,10}],ListPlot[dataQ8o1,PlotStyle->Red,ImageSize->800,PlotLegends->{"Numerical"},PlotMarkers->Style[\[FilledCircle],20]]]*)


(* ::Section:: *)
(*Q8o2 for fermions*)


(* ::Subsection:: *)
(*x \[RightArrow] 0*)


(* ::Code:: *)
(*$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}*)
(*temp= (Q8o2int[x,1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}//Simplify)/.{Abs[x_]->x}//FullSimplify*)


(* ::Code:: *)
(*plot=Plot[NIntegrate[temp/.{x->0.001,T->1},{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]*)


(* ::Code:: *)
(*Q8o2lowX[x_,vw_,T_] = (Re[Integrate[AsymptoticIntegrate[temp,{y,-1,1},{x,0,1}],{z,0,Infinity}]]//Simplify)//TrigToExp//ComplexExpand//Refine[#,{vw \[Element] Reals, 1>vw>0}]&//Simplify*)
(*Q8o2[0.001,1,0.2,1]/Q8o2lowX[0.001,0.2,1]-1*)
(*Q8o2[0.0001,1,0.5,1]/Q8o2lowX[0.0001,0.5,1]-1*)
(*Q8o2[0.00001,1,0.999,1]/Q8o2lowX[0.00001,0.999,1]-1*)


Q8o2lowX[x_,vw_,T_] =-((3 Sqrt[1-vw^2] (-2+x))/(8 \[Pi] T^2));


(* ::Code:: *)
(*Q8o2lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify*)
(*StringReplace[ToString[Q8o2lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow"}]*)


(* ::Code:: *)
(*vwall=0.4;*)
(*ErrorPlotData=Table[{10^i,Abs[Q8o2lowX[10^i,vwall,1]/Q8o2[10^i,1,vwall,1]-1]},{i,-12,0,1}]*)
(*ListLogLogPlot[ErrorPlotData]*)


(* ::Code:: *)
(*Show[plot,Plot[Q8o2lowX[1,vw,1],{vw,0,1},PlotStyle->Orange]]*)


(* ::Subsection:: *)
(*Sanity check*)


(* ::Code:: *)
(*Q8o2[0.001,1,0.2,1]/Q8o2lowX[0.001,0.2,1]-1*)
(*Q8o2[0.0001,1,0.5,1]/Q8o2lowX[0.0001,0.5,1]-1*)
(*Q8o2[0.00001,1,0.999,1]/Q8o2lowX[0.00001,0.999,1]-1*)


(* ::Subsection:: *)
(* x \[RightArrow] \[Infinity]*)


(* ::Code:: *)
(*$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}*)


(* ::Code:: *)
(*Timing[temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,3}];*)
(*temp3=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,4}];*)
(*temp4=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,5}];*)
(*temp5=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,6}];*)
(*temp6=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*templist={temp3,temp4,temp5,temp6};*)


(* ::Code:: *)
(*vwall=0.1;*)
(*dataQ8o2 = Table[{x,Q8o2[x,1,vwall,1]},{x,0.1,20,0.5}];*)


(* ::Code:: *)
(*Table[Q8o2Expansion[x_,vw_,T_, sterm, t]= Integrate[(Series[templist[[t]],{s,0,sterm}]//Normal)/.{s->1},{z,0,Infinity}],{sterm,0,3},{t,4}];*)


(* ::Code:: *)
(*Show[ListPlot[dataQ8o2,PlotStyle->Green,ImageSize->800,PlotRange->{0,0.01}],Plot[Evaluate[Table[Q8o2Expansion[x,vwall,1, sterm, t],{sterm,0,3},{t,{1,3}}]],{x,1,10},PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,4}],PlotRange->All]]*)


(* ::Code:: *)
(*Q8o2Expansion[x,vw,T, 1, 1]*)


(* ::Code:: *)
(*Show[ListPlot[dataQ8o2,PlotStyle->Green,ImageSize->800,PlotRange->{0,0.01}],Plot[Evaluate[Table[Q8o2Expansion[x,vwall,1, sterm, t],{sterm,1,3},{t,3,4}]],{x,0.01,10},PlotRange->All]]*)


(* ::Code:: *)
(*c=0;*)
(*Flatten[Table[c+=1;{c, "s =", sterm, " | t =", t},{t,{1,2,3,4}},{sterm,{0,1,2,3}}],1]//TableForm*)


(* ::Code:: *)
(*plotlist=Evaluate[Flatten[Table[Table[{i[[1]],Abs[(i[[2]]/Q8o2Expansion[i[[1]],vwall,1, sterm, t])-1]},{i,dataQ8o2}],{t,{1,2,3,4}},{sterm,{0,1,2,3}}],1]];*)
(*Manipulate[ListLogPlot[plotlist, *)
(*  PlotStyle -> (Opacity@Boole@MemberQ[x, #] & /@ Range@Length@plotlist),ImageSize->1000,PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,4}],PlotRange->All], *)
(*    {{x, {12,13},"Approximation"}, Dynamic@Range@Length@plotlist, TogglerBar}]*)


(* ::Code:: *)
(*Show[LogPlot[{Q8o2Expansion[x,vwall,1, 3, 3],Q8o2Expansion[x,vwall,1, 0, 4],Q8o2lowX[x,vwall,1]},{x,0.01,10},PlotLegends->{"3, 3","0, 4", "around x"}],ListLogPlot[dataQ8o2,PlotStyle->Red,ImageSize->800,PlotLegends->{"Numerical"},PlotMarkers->Style[\[FilledCircle],20],PlotRange->{0,0.01}]]*)


Q8o2highX[x_,vw_,T_]=(1/(43486543872 Sqrt[2] \[Pi]^(3/2) T^2 x^5 Sqrt[
 x - vw^2 x]))E^(-4 x) (-1 + 
   vw^2) (-497664 E^(
    3 x) (247726080 vw^12 + 491520 vw^10 (-2323 + 56 x) + 
      6144 vw^8 (342669 - 17232 x + 640 x^2) + 
      256 vw^6 (-7863933 + 612888 x - 46976 x^2 + 3072 x^3) + 
      8 vw^4 (115513755 - 11358048 x + 1304832 x^2 - 200704 x^3 + 
         32768 x^4) + 
      8 (598299 - 277344 x + 141568 x^2 - 69632 x^3 + 32768 x^4) + 
      vw^2 (-137909835 + 16940376 x - 4157184 x^2 + 1640448 x^3 - 
         688128 x^4 + 262144 x^5)) - 
   2048 Sqrt[3] E^
    x (82575360 vw^12 + 163840 vw^10 (-2323 + 168 x) + 
      6144 vw^8 (114223 - 17232 x + 1920 x^2) + 
      256 vw^6 (-2621311 + 612888 x - 140928 x^2 + 27648 x^3) + 
      8 vw^4 (38504585 - 11358048 x + 3914496 x^2 - 1806336 x^3 + 
         884736 x^4) + 
      8 (199433 - 277344 x + 424704 x^2 - 626688 x^3 + 884736 x^4) + 
      3 vw^2 (-15323315 + 5646792 x - 4157184 x^2 + 4921344 x^3 - 
         6193152 x^4 + 7077888 x^5)) + 
   7776 Sqrt[2] E^(
    2 x) (247726080 vw^12 + 491520 vw^10 (-2323 + 112 x) + 
      6144 vw^8 (342669 - 34464 x + 2560 x^2) + 
      256 vw^6 (-7863933 + 1225776 x - 187904 x^2 + 24576 x^3) + 
      8 vw^4 (115513755 - 22716096 x + 5219328 x^2 - 1605632 x^3 + 
         524288 x^4) + 
      8 (598299 - 554688 x + 566272 x^2 - 557056 x^3 + 524288 x^4) + 
      vw^2 (-137909835 + 33880752 x - 16628736 x^2 + 13123584 x^3 - 
         11010048 x^4 + 8388608 x^5)) + 
   243 (247726080 vw^12 + 491520 vw^10 (-2323 + 224 x) + 
      6144 vw^8 (342669 - 68928 x + 10240 x^2) + 
      256 vw^6 (-7863933 + 2451552 x - 751616 x^2 + 196608 x^3) + 
      8 vw^4 (115513755 - 45432192 x + 20877312 x^2 - 12845056 x^3 + 
         8388608 x^4) + 
      8 (598299 - 1109376 x + 2265088 x^2 - 4456448 x^3 + 
         8388608 x^4) + 
      vw^2 (-137909835 + 67761504 x - 66514944 x^2 + 104988672 x^3 - 
         176160768 x^4 + 268435456 x^5)));


(* ::Code:: *)
(*Q8o2Expansion[x,vwall,1, 3, 3]/.{E^x_->exp[x]}//Simplify*)
(*StringReplace[ToString[Q8o2highX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow"}]*)


(* ::Section:: *)
(*Q9o1 for fermions*)


(* ::Subsection:: *)
(* x \[RightArrow] 0 *)


(* ::Code:: *)
(*$Assumptions= x>0&&z>0&&1>y>-1&&1>vw>0*)


(* ::Code:: *)
(*temp= (x Q9o1int[x,1,vw,T]//.{Re[x_]->x,w->x+z x,Abs[x_]->x}//Simplify)/.{Abs[x_]->x}//FullSimplify/.{Abs[x_]->x}*)
(*plot=Plot[NIntegrate[temp/.{x->0.0001,T->1},{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]*)


Q9o1lowX[x_,vw_,T_]= -0.0508343947899479 vw/T^4


(* ::Code:: *)
(*StringReplace[ToString[Q9o1lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow"}]*)


(* ::Code:: *)
(*Show[plot,Plot[Q9o1lowX[0.0001,vw,1],{vw,0,1},PlotStyle->Orange]]*)


(* ::Subsection:: *)
(*x \[RightArrow] \[Infinity]*)


(* ::Code:: *)
(*temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,3}];*)
(*temp3=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,4}];*)
(*temp4=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,5}];*)
(*temp5=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,6}];*)
(*temp6=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,7}];*)
(*temp7=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.Abs[x_]->x;*)
(*temp=Asymptotic[temp,{x,Infinity,8}];*)
(*temp8=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*templist={temp3,temp4,temp5,temp6,temp7, temp8};*)


(* ::Code:: *)
(*vwall=0.1;*)
(*dataQ9o1 = Table[{x,Q9o1[x,1,vwall,1]},{x,0.1,20,0.5}];*)


(* ::Code:: *)
(*Table[Print["s", sterm," | t ", t];Q9o1Expansion[x_,vw_,T_, sterm, t]= Integrate[(Series[templist[[t]],{s,0,sterm}]//Normal)/.{s->1},{z,0,Infinity}],{sterm,0,2},{t,templist//Length}];*)


(* ::Code:: *)
(*Show[ListPlot[dataQ9o1,PlotStyle->Green,ImageSize->800,PlotRange->{-0.01,0.01}],Plot[Evaluate[Table[Q9o1Expansion[x,vwall,1, sterm, t],{sterm,0,1},{t,{1,3}}]],{x,1,10},PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,templist//Length}],PlotRange->All]]*)


(* ::Code:: *)
(*c=0;*)
(*Flatten[Table[c+=1;{c, "s =", sterm, " | t =", t},{t,templist//Length},{sterm,{0,1,2}}],1]//TableForm*)


(* ::Code:: *)
(*plotlist=Evaluate[Flatten[Table[Table[{i[[1]],Abs[(i[[2]]/Q9o1Expansion[i[[1]],vwall,1, sterm, t])-1]},{i,dataQ9o1}],{t,templist//Length},{sterm,{0,1,2}}],1]];*)
(*Manipulate[ListLogPlot[plotlist, *)
(*  PlotStyle -> (Opacity@Boole@MemberQ[x, #] & /@ Range@Length@plotlist),ImageSize->1000,PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,templist//Length}],PlotRange->All], *)
(*    {{x, {17},"Approximation"}, Dynamic@Range@Length@plotlist, TogglerBar}]*)


(* ::Code:: *)
(*Show[ListPlot[dataQ9o1,PlotStyle->Green,ImageSize->800,PlotRange->{-0.01,0.01}],Plot[Evaluate[{Q9o1Expansion[x,vwall,1, 2, 5],Q9o1Expansion[x,vwall,1, 1, 6]}],{x,0.1,10},PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,templist//Length}],PlotRange->All]]*)


Q9o1highX[x_,vw_,T_]=(1/(536870912 \[Pi]^(3/2) T^4 Sqrt[2 - 2 vw^2] x^(
 15/2)))3 E^(-2 x) vw (-1 + 
   vw^2) (64 E^
    x (-11861352435 - 38289801216 vw^12 + 2014884240 x - 
      345989760 x^2 + 50515968 x^3 - 1277952 x^4 - 4718592 x^5 + 
      4194304 x^6 - 12582912 vw^10 (-15399 + 248 x) - 
      32768 vw^8 (12248175 - 414192 x + 9088 x^2) - 
      32768 vw^6 (-13266975 + 709128 x - 33408 x^2 + 1024 x^3) - 
      128 vw^4 (1895514939 - 143811168 x + 11393280 x^2 - 
         823296 x^3 + 32768 x^4) + 
      512 vw^2 (140415915 - 16800480 x + 2310400 x^2 - 315392 x^3 + 
         32768 x^4)) + 
   Sqrt[2] (11861352435 + 38289801216 vw^12 - 4029768480 x + 
      1383959040 x^2 - 404127744 x^3 + 20447232 x^4 + 150994944 x^5 - 
      268435456 x^6 + 12582912 vw^10 (-15399 + 496 x) + 
      32768 vw^8 (12248175 - 828384 x + 36352 x^2) + 
      32768 vw^6 (-13266975 + 1418256 x - 133632 x^2 + 8192 x^3) + 
      128 vw^4 (1895514939 - 287622336 x + 45573120 x^2 - 
         6586368 x^3 + 524288 x^4) - 
      512 vw^2 (140415915 - 33600960 x + 9241600 x^2 - 2523136 x^3 + 
         524288 x^4)));


(* ::Code:: *)
(*StringReplace[ToString[Q9o1highX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow","Sqrt"->"sqrt"}]*)


(* ::Subsection:: *)
(*Sanity check*)


(* ::Code:: *)
(*Q9o1[0.01,1,0.2,1]/Q9o1lowX[0.01,0.2,1]-1*)
(*Q9o1[0.01,1,0.5,1]/Q9o1lowX[0.01,0.5,1]-1*)
(*Q9o1[0.01,1,0.99,1]/Q9o1lowX[0.01,0.99,1]-1*)


(* ::Code:: *)
(*Show[Plot[{Q9o1lowX[x,vwall,1],Q9o1highX[x,vwall,1]},{x,0.01,10},PlotLegends->{"3, 3","0, 4", "around x"}],ListPlot[dataQ9o1,PlotStyle->Red,ImageSize->800,PlotLegends->{"Numerical"},PlotMarkers->Style[\[FilledCircle],20],PlotRange->{0,0.01}]]*)


(* ::Section:: *)
(*Q9o2 for fermions*)


(* ::Subsection:: *)
(*x \[RightArrow] 0*)


(* ::Code:: *)
(*$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}*)
(*temp= (x Q9o2int[x,1,vw,T]//.{Re[x_]->x,w->x+z x,Abs[x_]->x}//Simplify)/.{Abs[x_]->x}//FullSimplify*)
(*asym = Asymptotic[temp,{x,0,1}];*)
(*LimitAt0VW = Integrate[Asymptotic[asym,vw->0]//Simplify,{z,0,Infinity},{y,-1,1}]*)
(*Q9o2lowX[x_,vw_,T_]=LimitAt0VW Sqrt[1-vw^2]*)


(* ::Code:: *)
(*StringReplace[ToString[Q9o2lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow","Sqrt"->"sqrt"}]*)


(* ::Code:: *)
(*vwall=0.4;*)
(*ErrorPlotData=Table[{10^i,Abs[Q9o2lowX[10^i,vwall,1]/Q9o2[10^i,1,vwall,1]-1]},{i,-12,0,1}]*)
(*ListLogLogPlot[ErrorPlotData]*)


(* ::Code:: *)
(*plot=Plot[NIntegrate[temp/.{x->0.001,T->1},{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]*)


(* ::Code:: *)
(*Show[plot,Plot[Q9o2lowX[0.001,vw,1],{vw,0,1},PlotStyle->Directive[Dashed,Orange]]]*)


Q9o2lowX[x_,vw_,T_]=(3 (-2+\[Pi]) Sqrt[1-vw^2])/(16 \[Pi]^2 T^4 x);


(* ::Subsection:: *)
(*x \[RightArrow] \[Infinity]*)


(* ::Code:: *)
(*$Assumptions= x>0&&z>0&&1>y>-1&&1>vw>0&&T>0&&s\[Element]Integers*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,1}];*)
(*temp1=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,2}];*)
(*temp2=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,3}];*)
(*temp3=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,4}];*)
(*temp4=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,5}];*)
(*temp5=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,6}];*)
(*temp6=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,7}];*)
(*temp7=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*Timing[temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z};*)
(*temp=(temp/.Abs[x_]->x//Simplify)/.{Abs[x_]->x};*)
(*temp=Asymptotic[temp,{x,Infinity,8}];*)
(*temp8=Integrate[temp,{y,-1,1}]]*)


(* ::Code:: *)
(*templist={temp3,temp4,temp5,temp6,temp7, temp8};*)


(* ::Code:: *)
(*vwall=0.1;*)
(*dataQ9o2 = Table[{x,Q9o2[x,1,vwall,1]},{x,0.1,20,0.5}];*)


(* ::Code:: *)
(*Table[Print["s", sterm," | t ", t];Q9o2Expansion[x_,vw_,T_, sterm, t]= Integrate[(Series[templist[[t]],{s,0,sterm}]//Normal)/.{s->1},{z,0,Infinity}],{sterm,0,2},{t,templist//Length}];*)


(* ::Code:: *)
(*Show[ListPlot[dataQ9o2,PlotStyle->Green,ImageSize->800,PlotRange->{0,0.01}],Plot[Evaluate[Table[Q9o2Expansion[x,vwall,1, sterm, t],{sterm,0,1},{t,{1,3}}]],{x,1,10},PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,templist//Length}],PlotRange->All]]*)


(* ::Code:: *)
(*plotlist=Evaluate[Flatten[Table[Table[{i[[1]],Abs[(i[[2]]/Q9o2Expansion[i[[1]],vwall,1, sterm, t])-1]},{i,dataQ9o2}],{t,templist//Length},{sterm,{0,1,2}}],1]];*)
(*Manipulate[ListLogPlot[plotlist, *)
(*  PlotStyle -> (Opacity@Boole@MemberQ[x, #] & /@ Range@Length@plotlist),ImageSize->1000,PlotLegends->Table["s = " <> ToString[sterm] <> " | t " <>ToString[t],{sterm,0,3},{t,templist//Length}],PlotRange->All], *)
(*    {{x, {16},"Approximation"}, Dynamic@Range@Length@plotlist, TogglerBar}]*)


(* ::Code:: *)
(*c=0;*)
(*Flatten[Table[c+=1;{c, "s =", sterm, " | t =", t},{t,templist//Length},{sterm,{0,1,2}}],1]//TableForm*)


Q9o2highX[x_,vw_,T_]=-(1/(67108864 Sqrt[2] \[Pi]^(3/2) T^4 x^8 Sqrt[x-vw^2 x]))3 E^-x (-1+vw^2) (137514115399680 vw^16+188743680 vw^14 (-4320347+24024 x)+47185920 vw^12 (43382485-493164 x+7392 x^2)+491520 vw^10 (-5728208877+99970872 x-3111808 x^2+64512 x^3)+49152 vw^8 (46558313955-1110376650 x+54103856 x^2-2348800 x^3+71680 x^4)+384 vw^6 (-2866546545795+87757928760 x-5985312512 x^2+410904576 x^3-26640384 x^4+1310720 x^5)+32 vw^4 (9144081088965-342529804260 x+30764575776 x^2-2999276544 x^3+316563456 x^4-34471936 x^5+3145728 x^6)+8 (124058570085-6548070960 x+1120404864 x^2-231960576 x^3+56131584 x^4-15204352 x^5+4194304 x^6)+vw^2 (-35248384552995+1539976101480 x-181378098048 x^2+25188559872 x^3-4185489408 x^4+811859968 x^5-171966464 x^6+33554432 x^7));


(* ::Code:: *)
(*StringReplace[ToString[Q9o2lowX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow"}]*)
(*StringReplace[ToString[Q9o2highX[x,vw,T]/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,x)"->"exp(x)","Power"->"pow"}]*)


(* ::Code:: *)
(*Show[LogPlot[{Q9o2lowX[x,vwall,1],Q9o2highX[x,vwall,1]},{x,0.01,10},PlotLegends->{"3, 3","0, 4", "around x"}],ListLogPlot[dataQ9o2,PlotStyle->Red,ImageSize->800,PlotLegends->{"Numerical"},PlotMarkers->Style[\[FilledCircle],20],PlotRange->{0,0.01}]]*)


(* ::Subsection:: *)
(*Sanity check*)


(* ::Code:: *)
(*Q9o2[0.01,1,0.2,1]/Q9o2lowX[0.01,0.2,1]-1*)
(*Q9o2[0.0001,1,0.5,1]/Q9o2lowX[0.0001,0.5,1]-1*)
(*Q9o2[0.00001,1,0.999,1]/Q9o2lowX[0.00001,0.999,1]-1*)


(* ::Code:: *)
(*Q9o2[10,1,0.2,1]/Q9o2highX[10,0.2,1]-1*)
(*Q9o2[100,1,0.5,1]/Q9o2highX[100,0.5,1]-1*)
(*Q9o2[200,1,0.999,1]/Q9o2highX[200,0.999,1]-1*)


(* ::Section:: *)
(*To L AT EX (run the input cells first)*)


(* ::Code:: *)
(*ToLaTeX[exp_]:=StringReplace[ToString[exp/.{E^x_->exp[x]}//Simplify//CForm],{"Power(E,"->"exp(","Power"->"pow","Sqrt"->"sqrt","Log"->"log"}]*)


(* ::Code:: *)
(*Q8o1lowX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q8o1highX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q8o2lowX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q8o2highX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q9o1lowX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q9o1highX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q9o2lowX[x,vw,T]//ToLaTeX*)


(* ::Code:: *)
(*Q9o2highX[x,vw,T]//ToLaTeX*)


(* ::Section:: *)
(*Q9o1 x \[RightArrow] 0 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}
temp=x^3 Q9o1int[x,-1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
function[z_]=z^2 x;
temp2 = function'[z](temp/.{z->function[z]})//Simplify
asymp = Asymptotic[temp2,x->0]//Simplify
plot=Plot[NIntegrate[Asymptotic[asymp/.T->1,vw->0],{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]


(* ::Subsubsection:: *)
(*Limit at vw = 0*)


IntOverY = Integrate[Asymptotic[asymp,vw->0]//Simplify,{y,-1,1}]
Q9o1lowxlowvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


(* ::Subsubsection:: *)
(*Limit at vw = 1*)


IntOverY = Integrate[Asymptotic[asymp,vw->1]//Simplify,{y,-1,1}]
Q9o1lowxhighvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


guess = ((vw)(Q9o1lowxhighvw[0.00001,1,vw])(a + b vw + c vw ^2 + d vw^3 + f vw^4))
Q9o1lowx[x_,T_,vw_]= (guess/.FindFit[plot[[1,1,1,3,1,2,1]], guess,{a,b,c,d,e,f},vw])/(x^3 T^4)


Show[plot,Plot[{Q9o1lowxhighvw[0.00001,1,vw]},{vw,0,1},PlotStyle->Orange],Plot[0.00001^3 Q9o1lowx[0.00001,1,vw],{vw,0,1},PlotStyle->Green]]


Q9o1[0.001,-1,0.2,1]/Q9o1lowx[0.001,1,0.2]-1
Q9o1[0.0001,-1,0.5,1]/Q9o1lowx[0.0001,1,0.5]-1
Q9o1[0.00001,-1,0.99,1]/Q9o1lowx[0.00001,1,0.99]-1


(* ::Section:: *)
(*Q9o2 x \[RightArrow] 0 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}
temp=x^3 Q9o2int[x,-1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
function[z_]=z^2 x;
temp2 = function'[z](temp/.{z->function[z]})//Simplify
asymp = Asymptotic[temp2,x->0]//Simplify
plot=Plot[NIntegrate[Asymptotic[asymp/.T->1,vw->0],{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]


(* ::Subsubsection:: *)
(*Limit at vw = 0*)


IntOverY = Integrate[Asymptotic[asymp,vw->0]//Simplify,{y,-1,1}]
Q9o2lowxlowvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


(* ::Subsubsection:: *)
(*Limit at vw = 1*)


IntOverY = Integrate[Asymptotic[asymp,vw->1]//Simplify,{y,-1,1}]
Q9o2lowxhighvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


Q9o2lowxlowvw[0.00001,1,vw]


guess = Q9o2lowxlowvw[0.00001,1,vw] (a+b vw + c vw^2 + h vw^3)/(d+e vw + g vw^2) Q9o2lowxhighvw[0.00001,1,vw]
Q9o2lowx[x_,T_,vw_]= (guess/.FindFit[plot[[1,1,1,3,1,2,1]], guess,{a,b,c,d,e,g,h},vw])/(x^3 T^4)
Show[plot,Plot[{Q9o2lowxlowvw[0.00001,1,vw],Q9o2lowxhighvw[0.00001,1,vw]},{vw,0,1},PlotStyle->Orange,PlotRange->{0,0.04}],Plot[0.00001^3 Q9o2lowx[0.00001,1,vw],{vw,0,1},PlotStyle->Directive[Purple,Dashed]]]


Q9o2[0.001,-1,0.2,1]/Q9o2lowx[0.001,1,0.2]-1
Q9o2[0.0001,-1,0.5,1]/Q9o2lowx[0.0001,1,0.5]-1
Q9o2[0.00001,-1,0.99,1]/Q9o2lowx[0.00001,1,0.99]-1

