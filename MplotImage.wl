(* ::Package:: *)

rofm[b_,M_]/;b>bcri[M]&&Abs[G[umin[M],b,M]]> 10^-16:=Block[{geof,tphi,uroot,fun,phi,tphim},
{geof,tphi,uroot}=geo[b,M];
fun[phi_]/;phi<tphi:=geof[phi];
fun[phi_]/;phi==tphi:=uroot;
fun[phi_]/;phi>tphi:=fun[2tphi-phi];
tphim=(2tphi/Pi-0.5)//Floor;
tphim=Range[0,tphim];
tphim=0.5 Pi+Pi tphim;
If[tphim=={},Return[{}]];
1/(fun/@tphim)
]
rofm[b_,M_]/;b<= bcri[M]&&G[umin[M],b,M]>= 10^-7:=Block[{geof,tphi,uroot,phi,rp,tphim},
{geof,tphi,uroot}=geo[b,M];
rp=horizons[M]//Last;
tphi=NIntegrate[1/Sqrt[G[u,b,M]],{u,0,1/rp}];
tphim=(tphi/Pi-0.5)//Floor;
tphim=Range[0,tphim];
tphim=0.5 Pi+Pi tphim;
If[tphim=={},Return[{}]];
1/(geof/@tphim)
]
rofm[b_,M_]/;b<= bcri[M]&&G[umin[M],b,M] <10^-7&&Abs[G[umin[M],b,M]]> 10^-16:=Block[
{geof,tphi,tphip,uroot,phi,rp,up,uc,bc,ep,Gc,v,y,sep,u2,u1,err,n,rest,crest,u,tphim},
{geof,tphip,uroot}=geo[b,M];
rp=horizons[M]//Last;
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
(*************************)
tphi=tphi+rest;
(******)
tphim=(tphi/Pi-0.5)//Floor;
tphim=Range[0,tphim];
tphim=0.5 Pi+Pi tphim;
If[tphim=={},Return[{}]];
1/(geof/@tphim)
]
rofm[b_,M_]/;Abs[G[umin[M],b,M]]<= 10^-16:={}(*because the light at photo sphere cannot be recieved by us*)


rofmp[b_,M_]/;b<= bcri[M]&&G[umin[M],b,M] >= 10^-7:=Block[{geof,tphi,tphip,uroot,phi,rp,tphim,tphim1},
{geof,tphip,uroot}=geo[b,M];
rp=horizons[M]//Last;
tphi=NIntegrate[1/Sqrt[G[u,b,M]],{u,0,1/rp},MaxRecursion->12];
(*the light come out the horizon from phi=2tphip-tphi*)
tphim1=1+((2tphip-tphi)/Pi-0.5)//Floor;
tphim=((2tphip)/Pi-0.5)//Floor;
tphim=Range[tphim1,tphim];
tphim=0.5 Pi+Pi tphim;
If[tphim=={},Return[{}]];
1/(geof[2tphip-#]&/@tphim)
]

rofmp[b_,M_]/;b<= bcri[M]&&G[umin[M],b,M] <10^-7&&Abs[G[umin[M],b,M]]> 10^-16:=Block[{geof,tphi,tphip,uroot,phi,rp,up,uc,bc,Gc,v,y,sep,ep,u2,u1,err,n,rest,crest,u,tphim,tphim1},
{geof,tphip,uroot}=geo[b,M];
rp=horizons[M]//Last;
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
(**********************************************************************)
tphi=tphi+rest;
(*****************************************************************************************************)
(*the light come out the horizon from phi=2tphip-tphi*)
tphim1=1+((2tphip-tphi)/Pi-0.5)//Floor;
tphim=((2tphip)/Pi-0.5)//Floor;
tphim=Range[tphim1,tphim];
tphim=0.5 Pi+Pi tphim;
If[tphim=={},Return[{}]];
1/(geof[2tphip-#]&/@tphim)
]
rofmp[b_,M_]/;Abs[G[umin[M],b,M]]<= 10^-16||b>bcri[M]:={}


Iob[Iemw_,Iemb_,b_,M_,"Wh"]:=Block[{rofmv,rofmpv,iobw,iobh},
rofmpv=rofmp[b,M];
iobw=((alpha M^-2/#^4-2/#+1)^2 Iemw[#])&/@rofmpv;
iobw=If[Length[iobw]<=3,Total[iobw],Total[iobw[[1;;3]]]]
]
Iob[Iemw_,Iemb_,b_,M_,"Bh"]:=Block[{rofmv,rofmpv,iobw,iobh},
rofmv=rofm[b,M];
iobh=((alpha M^-2/#^4-2/#+1)^2 Iemb[#])&/@rofmv;
If[Length[iobh]<=3,Total[iobh],Total[iobh[[1;;3]]]]
]
Iob[Iemw_,Iemb_,b_,M_,"Bh-Wh"]:=Block[{rofmv,rofmpv,iobw,iobh},
Iob[Iemw,Iemb,b,M,"Bh"]+Iob[Iemw,Iemb,b,M,"Wh"]
]
