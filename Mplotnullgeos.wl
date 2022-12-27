(* ::Package:: *)

plotgeo[b_,M_,plotstyle_,cf_:None]/;b>bcri[M]&&Abs[G[umin[M],b,M]]> 10^-16:=Block[{geof,tphi,uroot,fun,phi},
{geof,tphi,uroot}=geo[b,M];
fun[phi_]/;phi<tphi:=geof[phi];
fun[phi_]/;phi==tphi:=uroot;
fun[phi_]/;phi>tphi:=fun[2tphi-phi];
ParametricPlot[{1/fun[phi]Cos[phi],1/fun[phi]Sin[phi]},{phi,10^-16,2tphi},PlotStyle->plotstyle,ColorFunction->cf]
]
plotgeo[b_,M_,plotstyle_,cf_:None,rpm_:1]/;b<= bcri[M]&&Abs[G[umin[M],b,M]]>= 10^-7:=Block[{geof,tphi,uroot,phi,rp},
{geof,tphi,uroot}=geo[b,M];
rp=horizons[M][[-rpm]];
tphi=NIntegrate[1/Sqrt[G[u,b,M]],{u,0,1/rp}];
ParametricPlot[{1/geof[phi]Cos[phi],1/geof[phi]Sin[phi]},{phi,10^-16,tphi},PlotStyle->plotstyle,ColorFunction->cf]
]


plotgeo[b_,M_,plotstyle_,cf_:None,rpm_:1]/;b<= bcri[M]&&Abs[G[umin[M],b,M]]<10^-7&&G[umin[M],b,M]>10^-16:=
Block[{geof,tphi,tphip,uroot,bc,uc,Gc,v,sep,y,phi,rp,ep,u2,u1,err,n,rest,crest,u},
{geof,tphip,uroot}=geo[b,M];
rp=horizons[M][[-rpm]];
uc=umin[M];
bc=bcri[M];
ep=Sqrt[bc-b]/Sqrt[bc]^3;
u2=uc+ep;
u1=uc-ep;
tphi=NIntegrate[1/Re@Sqrt[G[u,b ,M]],{u,0,u1},MaxRecursion->12];
tphi=tphi+NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,u2,uroot},MaxRecursion->12];
Gc=CoefficientList[G[uc-v,bc,M],v][[3;;-1]].{1,v,v^2,v^3,v^4};
Gc=(v^2  Gc/.v->y sep)-1/bc^2+1/(bc-bc^3 sep^2)^2;
Gc=Series[Gc,{sep,0,6}]//Normal;
Gc=CoefficientList[Gc,sep][[3;;-1]];
Gc=Gc.(sep^#&/@Range[0,Length[Gc]-1 ]);
Gc=Series[1/Sqrt[Gc],{sep,0,3}]//Normal;
Gc=Gc/.sep->ep;
rest=NIntegrate[Gc,{y,-1,1}];
tphi=tphi+rest;
(*****************************************************************************************************)
ParametricPlot[{1/geof[phi]Cos[phi],1/geof[phi]Sin[phi]},{phi,10^-16,tphi},PlotStyle->plotstyle,ColorFunction->cf]
]

plotgeo[b_,M_,plotstyle_,cf_:None]/;Abs[G[umin[M],b,M]]<= 10^-16:=Block[{geof,tphi,uroot,phi,rp,fun,tphip},
{geof,tphi,uroot}=geo[b,M];
ParametricPlot[{1/uroot Cos[phi],1/uroot Sin[phi]},{phi,0,2Pi},PlotStyle->plotstyle,ColorFunction->cf]
]
(***************************************)
plotgeop[b_,M_,plotstyle_,cf_:None,rpm_:1]/;b<= bcri[M]&&Abs[G[umin[M],b,M]]< 10^-7&&G[umin[beta],b,beta]>10^-16:=Block[{geof,tphi,tphip,uc,bc,Gc,v,y,sep,ep,u2,u1,err,n,rest,crest,u,uroot,phi,rp,fun},
{geof,tphip,uroot}=geo[b,M];
rp=horizons[M][[-rpm]];
uc=umin[M];
bc=bcri[M];
ep=Sqrt[bc-b]/Sqrt[bc]^3;
u2=uc+ep;
u1=uc-ep;
tphi=NIntegrate[1/Re@Sqrt[G[u,b ,M]],{u,0,u1},MaxRecursion->12];
tphi=tphi+NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,u2,uroot},MaxRecursion->12];
Gc=CoefficientList[G[uc-v,bc,M],v][[3;;-1]].{1,v,v^2,v^3,v^4};
Gc=(v^2  Gc/.v->y sep)-1/bc^2+1/(bc-bc^3 sep^2)^2;
Gc=Series[Gc,{sep,0,6}]//Normal;
Gc=CoefficientList[Gc,sep][[3;;-1]];
Gc=Gc.(sep^#&/@Range[0,Length[Gc]-1 ]);
Gc=Series[1/Sqrt[Gc],{sep,0,3}]//Normal;
Gc=Gc/.sep->ep;
rest=NIntegrate[Gc,{y,-1,1}];
tphi=tphi+rest;
(*****************************************************************************************************)
fun[phi_]/;phi<tphip:=geof[phi];
fun[phi_]/;phi==tphip:=uroot;
fun[phi_]/;phi>tphip:=fun[2tphip-phi];
ParametricPlot[{1/fun[phi]Cos[phi],1/fun[phi]Sin[phi]},{phi,2tphip-tphi,2tphi},PlotStyle->plotstyle,ColorFunction->cf]
]

plotgeop[b_,M_,plotstyle_,cf_:None,rpm_:1]/;b<= bcri[M]&&Abs[G[umin[M],b,M]]>= 10^-7:=Block[{geof,tphi,uroot,phi,rp,fun,tphip},
{geof,tphi,uroot}=geo[b,M];
rp=horizons[M][[-rpm]];
tphip=NIntegrate[1/Sqrt[G[u,b,M]],{u,0,1/rp}];
fun[phi_]/;phi<tphi:=geof[phi];
fun[phi_]/;phi==tphi:=uroot;
fun[phi_]/;phi>tphi:=fun[2tphi-phi];
ParametricPlot[{1/fun[phi]Cos[phi],1/fun[phi]Sin[phi]},{phi,2tphi-tphip,2tphi},PlotStyle->plotstyle,ColorFunction->cf]
]
plotgeop[b_,M_,plotstyle_,cf_:None]/;Abs[G[umin[M],b,M]]<= 10^-16:=Block[{geof,tphi,uroot,phi,rp,fun,tphip},
{geof,tphi,uroot}=geo[b,M];
ParametricPlot[{1/uroot Cos[phi],1/uroot Sin[phi]},{phi,0,2Pi},PlotStyle->plotstyle,ColorFunction->cf]]
