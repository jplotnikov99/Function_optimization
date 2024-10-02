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
(*Q8o1highX[100,0.5,1]/Q8o1[100,1,0.5,1]-1*)
(*Q8o2highX[100,0.5,1]/Q8o2[100,1,0.5,1]-1*)
(*Q9o1highX[100,0.5,1]/Q9o1[100,1,0.5,1]-1*)
(*Q9o2highX[100,0.5,1]/Q9o2[100,1,0.5,1]-1*)



