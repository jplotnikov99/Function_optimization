(* ::Package:: *)

(* ::Section:: *)
(*Generate code for testing*)


(* ::Input:: *)
(*singlegrid[n_,m_,k_,V_,s_, vw_,x_,T_]:=*)
(*Block[{gamw = 1/Sqrt[1-vw^2] , dfs={},Vx=1},*)
(*f0=1/(Exp[w]+s);*)
(*df0=D[f0,{w,k}];*)
(*pwt=Sqrt[w^2-x^2];*)
(*pzt=gamw*(y*pwt-w*vw);*)
(*ET=gamw*(w-vw*y*pwt);*)
(*If[V==1,Vx*=(Abs[pzt]/Sqrt[pzt^2+x^2])^2*1/Sqrt[1-pzt^2/ET^2]];*)
(*integr= -3/(Pi^2*gamw)*T^(n-m-k+1)*Vx*pwt*pzt^n/ET^(m-1)*df0//Simplify;*)
(*(*Print[integr];*)*)
(*dfs =NIntegrate[integr,{w,x,Infinity},{y,-1,1},PrecisionGoal->6,AccuracyGoal->8,MinRecursion->5,MaxRecursion->8,WorkingPrecision->8];*)
(*Return[dfs];*)
(*]*)
(*D0[x_,s_,vw_,T_]:=singlegrid[0,0,1,0,s,vw,x,T];*)
(*D1[x_,s_,vw_,T_]:=singlegrid[1,1,1,0,s,vw,x,T];*)
(*D2[x_,s_,vw_,T_]:=singlegrid[2,2,1,0,s,vw,x,T];*)
(*Q1[x_,s_,vw_,T_]:=singlegrid[0,1,2,0,s,vw,x,T];*)
(*Q2[x_,s_,vw_,T_]:=singlegrid[1,2,2,0,s,vw,x,T];*)
(*Qe1[x_,s_,vw_,T_]:=singlegrid[0,1,1,0,s,vw,x,T];*)
(*Qe2[x_,s_,vw_,T_]:=singlegrid[1,2,1,0,s,vw,x,T];*)
(*Q8o1[x_,s_,vw_,T_]:=0.5Re[singlegrid[-1,1,1,1,s,vw,x,T]];*)
(*Q8o2[x_,s_,vw_,T_]:=0.5Re[singlegrid[0,2,1,1,s,vw,x,T]];*)
(*Q9o1[x_,s_,vw_,T_]:=0.25(Re[singlegrid[-1,3,1,1,s,vw,x,T]-Re[singlegrid[-1,2,2,1,s,vw,x,T]/Sqrt[1-vw^2]]]);*)
(*Q9o2[x_,s_,vw_,T_]:=0.25(Re[singlegrid[0,4,1,1,s,vw,x,T]-singlegrid[0,3,2,1,s,vw,x,T]/Sqrt[1-vw^2]]);*)
(*N0[mass_,s_,T_]:=4Pi NIntegrate[p^2/(Exp[Sqrt[p^2+mass^2]/T]+s),{p,0,Infinity}] (*Auxiliary function*)*)
(*Rbar[mass_,s_,vw_,T_]:= Pi*(1-vw^2)/N0[mass,s,T] T Re[NIntegrate[Log[Abs[(Sqrt[w^2-(mass/T)^2]-vw w)/(Sqrt[w^2-(x)^2]+vw w)]]/( Exp[w]+s),{w,mass/T,Infinity},MaxRecursion->20]];*)
(*K0[x_,s_,vw_,T_]:=-singlegrid[0,0,0,0,s,vw,x,T];*)
(*K4[x_,s_,vw_,T_]:=singlegrid[2,2,1,0,s,0,x,T];*)


(* ::Input:: *)
(*singlegridint[n_,m_,k_,V_,s_, vw_,x_,T_]:=*)
(*Block[{gamw = 1/Sqrt[1-vw^2] , dfs={},Vx=1},*)
(*f0=1/(Exp[w]+s);*)
(*df0=D[f0,{w,k}];*)
(*pwt=Sqrt[w^2-x^2];pzt=gamw*(y*pwt-w*vw);*)
(*ET=gamw*(w-vw*y*pwt);*)
(*If[V==1,Vx*=(Abs[pzt]/Sqrt[pzt^2+x^2])^2*1/Sqrt[1-pzt^2/ET^2]];*)
(*integr= -3/(Pi^2*gamw)*T^(n-m-k+1)*Vx*pwt*pzt^n/ET^(m-1)*df0//Simplify;*)
(*(*Print[integr];*)*)
(*dfs =integr//Simplify;*)
(*Return[dfs];*)
(*]*)
(*D0int[x_,s_,vw_,T_]:=singlegridint[0,0,1,0,s,vw,x,T];*)
(*D1int[x_,s_,vw_,T_]:=singlegridint[1,1,1,0,s,vw,x,T];*)
(*D2int[x_,s_,vw_,T_]:=singlegridint[2,2,1,0,s,vw,x,T];*)
(*Q1int[x_,s_,vw_,T_]:=singlegridint[0,1,2,0,s,vw,x,T];*)
(*Q2int[x_,s_,vw_,T_]:=singlegridint[1,2,2,0,s,vw,x,T];*)
(*Qe1int[x_,s_,vw_,T_]:=singlegridint[0,1,1,0,s,vw,x,T];*)
(*Qe2int[x_,s_,vw_,T_]:=singlegridint[1,2,1,0,s,vw,x,T];*)
(*Q8o1int[x_,s_,vw_,T_]:=1/2Re[singlegridint[-1,1,1,1,s,vw,x,T]];*)
(*Q8o2int[x_,s_,vw_,T_]:=1/2Re[singlegridint[0,2,1,1,s,vw,x,T]];*)
(*Q9o1int[x_,s_,vw_,T_]:=1/4(Re[singlegridint[-1,3,1,1,s,vw,x,T]-Re[singlegridint[-1,2,2,1,s,vw,x,T]/Sqrt[1-vw^2]]]);*)
(*Q9o2int[x_,s_,vw_,T_]:=1/4(Re[singlegridint[0,4,1,1,s,vw,x,T]-singlegridint[0,3,2,1,s,vw,x,T]/Sqrt[1-vw^2]]);*)


(* ::Section:: *)
(*Q8o1 x = \[Infinity]*)


(* ::Input:: *)
(*temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z}*)
(*temp=temp/.Abs[x_]->x//Simplify*)
(*temp=Asymptotic[temp,x->Infinity]*)
(*temp=Integrate[temp,{y,-1,1}]*)
(*temp=temp/.{s->0}//FullSimplify*)
(*Q8o1highX[x_,vw_,T_]=Integrate[temp,{z,0,Infinity}]/.Sqrt[2-2 vw^2]->Sqrt[2] Sqrt[1- vw^2]*)


(* ::Section:: *)
(*Q8o2 x = \[Infinity]*)


(* ::Input:: *)
(*temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z}*)
(*temp=temp/.Abs[x_]->x//Simplify*)
(*temp=Asymptotic[temp,x->Infinity]*)
(*temp=Integrate[temp,{y,-1,1}]*)
(*temp=temp/.{s->0}//FullSimplify*)
(*Q8o2highX[x_,vw_,T_]=Integrate[temp,{z,0,Infinity}]/.Sqrt[2-2 vw^2]->Sqrt[2] Sqrt[1- vw^2]*)


(* ::Section:: *)
(*Q9o1 x = \[Infinity]*)


(* ::Input:: *)
(*temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z}*)
(*temp=temp/.Abs[x_]->x//Simplify*)
(*temp=Asymptotic[temp,{x,Infinity,2}]*)
(*temp=Integrate[temp,{y,-1,1}]//Simplify*)
(*temp=temp/.{s->0}//FullSimplify*)
(*Q9o1highX[x_,vw_,T_]=Integrate[temp,{z,0,Infinity}]/.{1/Sqrt[2-2 vw^2]->1/(Sqrt[2] Sqrt[1- vw^2])}*)


(* ::Section:: *)
(*Q9o2 x = \[Infinity]*)


(* ::Input:: *)
(*temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z}*)
(*temp=temp/.Abs[x_]->x//Simplify*)
(*temp=Asymptotic[temp,{x,Infinity,2}]*)
(*temp=Integrate[temp,{y,-1,1}]//Simplify*)
(*temp=temp/.{s->0}//FullSimplify*)
(*Q9o2highX[x_,vw_,T_]=Integrate[temp,{z,0,Infinity}]/.{1/Sqrt[2-2 vw^2]->1/(Sqrt[2] Sqrt[1- vw^2])}*)


(* ::Section:: *)
(*Sanity check for Q8o1/Q8o2/Q9o1/Q9o2*)


(* ::Input:: *)
(*Q8o1highX[10,0.5,1]/Q8o1[10,1,0.5,1]-1*)
(*Q8o2highX[10,0.5,1]/Q8o2[10,1,0.5,1]-1*)
(*Q9o1highX[10,0.5,1]/Q9o1[10,1,0.5,1]-1*)
(*Q9o2highX[10,0.5,1]/Q9o2[10,1,0.5,1]-1*)


(* ::Section:: *)
(*Q8o1 vw = 1 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}


temp=Q8o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
tempvw1 =Asymptotic[temp,vw->1]


tempvw1INTY = Integrate[tempvw1,{y,-1,1}]//Simplify


tempvw1INTYtemp = x^2 tempvw1INTY/.{z->z*x,s->-1}//Simplify
xtozeroAsym = Asymptotic[tempvw1INTYtemp,x->0]
Show[
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.01},{z,10^-10,1},PlotRange->All,PlotStyle->Red, GridLines -> {{{0.2592874577954832, Red}}, None}],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.001},{z,10^-10,1},PlotRange->All,PlotStyle->Orange],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.0001},{z,10^-10,1},PlotRange->All,PlotStyle->Green],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.00001},{z,10^-10,1},PlotRange->All,PlotStyle->Purple],
Plot[(xtozeroAsym)/.{vw->0.99,T->1},{z,10^-10,1},PlotRange->All,PlotStyle->Yellow]]


PartThatCanBeZero=(D[tempvw1INTYtemp,z]//FullSimplify)[[8]]


Asymptotic[PartThatCanBeZero,x->0]/.x->1


FindRoot[(Asymptotic[PartThatCanBeZero,x->0]/.x->1)==0,{z,0.2}]


Q8o1Vw1Lowx[x_,T_,vw_]=Integrate[xtozeroAsym,{z,0,Infinity}]/x


(* ::Subsection:: *)
(*Sanity check*)


Q8o1Vw1Lowx[0.001,1,0.99]/Q8o1[0.001,-1,0.99,1]-1
Q8o1Vw1Lowx[0.0001,1,0.99]/Q8o1[0.0001,-1,0.99,1]-1


(* ::Section:: *)
(*Q8o2 vw = 1 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}


temp=Q8o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
tempvw1 =Asymptotic[temp,vw->1]


tempvw1INTY = Integrate[tempvw1,{y,-1,1}]//Simplify


tempvw1INTYtemp = x^2 tempvw1INTY/.{z->z*x,s->-1}//Simplify
xtozeroAsym = Asymptotic[tempvw1INTYtemp,x->0]
Show[
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.01},{z,10^-10,1},PlotRange->All,PlotStyle->Red, GridLines -> {{{0.2592874577954832, Red}}, None}],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.001},{z,10^-10,1},PlotRange->All,PlotStyle->Orange],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.0001},{z,10^-10,1},PlotRange->All,PlotStyle->Green],
Plot[(tempvw1INTYtemp)/.{vw->0.99,T->1,x->0.00001},{z,10^-10,1},PlotRange->All,PlotStyle->Purple],
Plot[(xtozeroAsym)/.{vw->0.99,T->1},{z,10^-10,1},PlotRange->All,PlotStyle->Yellow]]


PartThatCanBeZero=(D[tempvw1INTYtemp,z]//FullSimplify)[[8]]


Asymptotic[PartThatCanBeZero,x->0]/.x->1


FindRoot[(Asymptotic[PartThatCanBeZero,x->0]/.x->1)==0,{z,0.2}]


Q8o2Vw1Lowx[x_,T_,vw_]=Integrate[xtozeroAsym,{z,0,Infinity}]/x


(* ::Subsection:: *)
(*Sanity check*)


Q8o2Vw1Lowx[0.001,1,0.999]/Q8o2[0.001,-1,0.999,1]-1
Q8o2Vw1Lowx[0.00001,1,0.999]/Q8o2[0.00001,-1,0.999,1]-1


(* ::Section:: *)
(*Q9o1 vw = 1 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}


temp=Q9o1int[x,s,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
tempvw1 =Asymptotic[temp,vw->1]


tempvw1INTY = x^3 Integrate[tempvw1,{y,-1,1}]/.s->-1//FullSimplify


symlog =
  {
   Function[x, Sign[x] (Log[Abs[x] + 1])],
   Function[y, Sign[y] (Exp[Abs[y]] - 1)]
   };



tempvw1INTY
parttosimplify = tempvw1INTY[[9]] * tempvw1INTY[[10]]
function[z_]=z^2 x;
tempvw1INTYtemp = function'[z](parttosimplify/.{z->function[z]})//FullSimplify
integrand = function'[z](tempvw1INTY/.{z->function[z]})//FullSimplify//TrigToExp//Simplify
Q9o1Vw1Lowx[x_,T_,vw_]=Integrate[Asymptotic[integrand,x->0],{z,0,Infinity}]/x^3


(* ::Subsection:: *)
(*Sanity check*)


Q9o1Vw1Lowx[0.001,1,.999]/Q9o1[0.001,-1,0.999,1]-1
Q9o1Vw1Lowx[0.0001,1,.999]/Q9o1[0.0001,-1,0.999,1]-1


(* ::Section:: *)
(*Q9o2 vw = 1 for fermions*)


$Assumptions= {x>0&&z>=0&&1>=y>=-1&&1>vw>0}


temp=Q9o2int[x,s,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
tempvw1 =Asymptotic[temp,vw->1]


tempvw1INTY = x^3 Integrate[tempvw1,{y,-1,1}]/.s->-1//FullSimplify


tempvw1INTY
parttosimplify = tempvw1INTY[[9]] * tempvw1INTY[[10]]
function[z_]=z^2 x;
tempvw1INTYtemp = function'[z](parttosimplify/.{z->function[z]})//FullSimplify
integrand = function'[z](tempvw1INTY/.{z->function[z]})//FullSimplify//TrigToExp//Simplify
Q9o2Vw1Lowx[x_,T_,vw_]=Integrate[Asymptotic[integrand,x->0],{z,0,Infinity}]/x^3


(* ::Subsection:: *)
(*Sanity check*)


Q9o2[0.001,-1,0.999,1]/Q9o2Vw1Lowx[0.001,1,0.999]-1
Q9o2[0.0001,-1,0.999,1]/Q9o2Vw1Lowx[0.0001,1,0.999]-1
Q9o2[0.00001,-1,0.999,1]/Q9o2Vw1Lowx[0.00001,1,0.999]-1


(* ::Section:: *)
(*Q8o1 x \[RightArrow] 0 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}
temp=x Q8o1int[x,-1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
function[z_]=z x;
temp2 = function'[z](temp/.{z->function[z]})//Simplify
asymp = Asymptotic[temp2,x->0]//Simplify
plot=Plot[NIntegrate[Asymptotic[asymp/.T->1,vw->0],{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]


guessfunction[vw_] := -vw*0.35 * Sqrt[1-vw^2] 
guess = (-vw*e * Sqrt[1-vw^2](a + b vw + c vw ^2 + d vw^3 + f vw^4))
Q8o1lowx[x_,T_,vw_]= (guess/.FindFit[plot[[1,1,1,3,1,2,1]], guess,{a,b,c,d,e,f},vw])/(x T^2)


(* ::Subsection:: *)
(*Sanity check*)


Q8o1[0.001,-1,0.2,1]/Q8o1lowx[0.001,1,0.2]-1
Q8o1[0.0001,-1,0.5,1]/Q8o1lowx[0.0001,1,0.5]-1
Q8o1[0.00001,-1,0.999,1]/Q8o1lowx[0.00001,1,0.999]-1


(* ::Section:: *)
(*Q8o2 x \[RightArrow] 0 for fermions*)


$Assumptions= {x>0&&z>0&&1>y>-1&&1>vw>0}
temp=x Q8o2int[x,-1,vw,T]//.{Re[x_]->x,w->x+z,Abs[x_]->x}
function[z_]=z x;
temp2 = function'[z](temp/.{z->function[z]})//Simplify
asymp = Asymptotic[temp2,x->0]//Simplify
plot=Plot[NIntegrate[Asymptotic[asymp/.T->1,vw->0],{y,-1,1},{z,0,Infinity}],{vw,0,1},PlotRange->All]


(* ::Subsubsection:: *)
(*Limit at vw = 0*)


IntOverY = Integrate[Asymptotic[asymp,vw->0]//Simplify,{y,-1,1}]
Q8o2lowxlowvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


(* ::Subsubsection:: *)
(*Limit at vw = 1*)


IntOverY = Integrate[Asymptotic[asymp,vw->1]//Simplify,{y,-1,1}]
Q8o2lowxhighvw[x_,T_,vw_]=Integrate[IntOverY,{z,0,Infinity}]


guess = ((Q8o2lowxhighvw[0.00001,1,vw])(Q8o2lowxhighvw[0.00001,1,vw]) (a+b vw + c vw^2 + h vw^3 + k vw^4)/(d+e vw + g vw^2))
Q8o2lowx[x_,T_,vw_]= (guess/.FindFit[plot[[1,1,1,3,1,2,1]], guess,{a,b,c,d,e,g,h,k},vw])/(x T^2)


Show[plot,Plot[{Q8o2lowxlowvw[0.00001,1,vw],Q8o2lowxhighvw[0.00001,1,vw]},{vw,0,1},PlotStyle->Orange],Plot[0.00001 Q8o2lowx[0.00001,1,vw],{vw,0,1},PlotStyle->Green]]


Q8o2[0.001,-1,0.2,1]/Q8o2lowx[0.001,1,0.2]-1
Q8o2[0.0001,-1,0.5,1]/Q8o2lowx[0.0001,1,0.5]-1
Q8o2[0.00001,-1,0.99,1]/Q8o2lowx[0.00001,1,0.99]-1


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

