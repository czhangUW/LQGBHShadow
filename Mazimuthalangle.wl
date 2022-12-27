(* ::Package:: *)

myEllipticPi[n_,th_,m_]:=If[th<=Pi/2,EllipticPi[n,th,m],2EllipticPi[n,Pi/2,m]-EllipticPi[n,Pi-th,m]]
mytheta[x_]:=If[x>=0,1,0]


phiout[M_,b_]/;b<= bcri[M]&&G[umin[M],b,M] >= 10^-7:=Block[{uroot,u,tphi,rp,tphip},
uroot=u/.NSolve[G[u,b,M]==0,u];
uroot=Cases[uroot,_?(Element[#,Reals]&&#>0&)]//Sort//First;
tphi=NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,0,uroot},MaxRecursion->12];
rp=horizons[M][[-1]];
tphip=NIntegrate[1/Sqrt[G[u,b,M]],{u,0,1/rp},MaxRecursion->12];
{2tphi-tphip,2tphi}
]
phiout[M_,b_]/;G[umin[M],b,M]>10^-16&&G[umin[M],b,M]<10^-7:=Block[
{uroot,u,uc,bc,ep,u2,u1,tphi,Gc,v,y,sep,err,n,rest,crest,rp,tphip},
uroot=turning[b,M];
uc=umin[M];
bc=bcri[M];
ep=Sqrt[bc-b]/(Sqrt[bc])^3;
u2=uc+ep;
u1=uc-ep;
tphi=NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,0,u1},MaxRecursion->12];
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
rp=horizons[M][[-1]];
tphip=NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,1/rp,uroot},MaxRecursion->12];
{tphi+tphip,2tphi}
]


gradient4logfit[M_]:=Block[{bc,uc,u,Gc,v,vroots,v2,v1,v3,v4,eq1,eq2,A1v,A2v,B1v,B2v,\[Lambda]v,\[Beta]v,k,n,A1,A2,B1,B2,\[Lambda],\[Beta],L,Cfp,CI,CIp,Ctot,um,up,vp,Cout,slopvout
},
bc=bcri[M];
(*uc=u/.NSolve[D[G[u,bc,betav],u]==0,u];
uc=Cases[uc,_?(Element[#,Reals]&&#>0&)]//Sort//First;*)
uc=umin[M];
Gc=CoefficientList[G[uc-v,bc,M],v][[3;;-1]].{1,v,v^2,v^3,v^4};
vroots=v/.Solve[Gc==0,v];
{v2,v1}=Cases[vroots,_?(Element[#,Reals]&)]//Sort;
{v3,v4}=Complement[vroots,{v1,v2}];
eq1=CoefficientList[(v1-v)(v-v2),v]==CoefficientList[A1v (v-\[Lambda]v)^2+B1v (v-\[Beta]v)^2,v]//Thread;
eq2=CoefficientList[v^2-2Re[v3]v+Re[v3]^2+Im[v3]^2,v]==CoefficientList[A2v (v-\[Lambda]v)^2+B2v (v-\[Beta]v)^2,v]//Thread;
{A1,A2,B1,B2,\[Lambda],\[Beta]}=({A1v,A2v,B1v,B2v,\[Lambda]v,\[Beta]v}/.NSolve[Join[eq1,eq2,{A1v<0}],{A1v,A2v,B1v,B2v,\[Lambda]v,\[Beta]v}])//First;
k=Sqrt[(A2 B1)/(A2 B1-A1 B2)];
n=(\[Beta]^2 B1)/(\[Lambda]^2 A1+\[Beta]^2 B1);
L[t_]:=Log[Abs[(Sqrt[A2 t^2+B2] Sqrt[A1 \[Lambda]^2+\[Beta]^2 B1]+Sqrt[A1 t^2+B1] Sqrt[A2 \[Lambda]^2+\[Beta]^2 B2])/(Sqrt[A2 t^2+B2] Sqrt[A1 \[Lambda]^2+\[Beta]^2 B1]-Sqrt[A1 t^2+B1] Sqrt[A2 \[Lambda]^2+\[Beta]^2 B2])]]/(Sqrt[A1 \[Lambda]^2+\[Beta]^2 B1] Sqrt[A2 \[Lambda]^2+\[Beta]^2 B2]);
Cfp=-k/(Sqrt[alpha]M^-1(\[Lambda]-\[Beta])\[Beta] Sqrt[A2 B1]) EllipticF[ArcCos[Sqrt[-A1/B1] (uc-\[Lambda])/(uc-\[Beta])],k^2];
Cfp=Cfp+(k n \[Lambda] A1)/(Sqrt[alpha]M^-1 \[Beta]^3 B1 Sqrt[A2 B1]) myEllipticPi[n,ArcCos[Sqrt[-A1/B1] (uc-\[Lambda])/(uc-\[Beta])],k^2]-1/(2Sqrt[alpha]M^-1) L[(uc-\[Lambda])/(uc-\[Beta])];
CI=Cfp+k/(Sqrt[alpha]M^-1(\[Lambda]-\[Beta])\[Beta] Sqrt[A2 B1]) EllipticF[ArcCos[Sqrt[-A1/B1] \[Lambda]/\[Beta]],k^2];
CIp=2myTh[-\[Lambda] \[Beta]]myEllipticPi[n,Pi/2,k^2];
CIp=CIp-Sign[\[Lambda] \[Beta]] myEllipticPi[k^2/n,ArcCos[Sqrt[-A1/B1]Abs[\[Lambda]/\[Beta]]],k^2];
CIp=CIp+Sign[\[Lambda] \[Beta]]EllipticF[ArcCos[Sqrt[-A1/B1]Abs[\[Lambda]/\[Beta]]],k^2];
CIp=(k n \[Lambda] A1)/(Sqrt[alpha]M^-1 \[Beta]^3 B1 Sqrt[A2 B1]) CIp;
CI=CI-CIp;
CI=CI+Log[Abs[(8 bc^3 alpha M^-2 (B1 \[Beta]^2+A1 \[Lambda]^2)^3 (B2 \[Beta]^2+A2 \[Lambda]^2)^3)/((A2 B1-A1 B2) (\[Beta]-\[Lambda])^2 (B1 B2 \[Beta]^4+2 A2 B1 \[Beta]^2 \[Lambda]^2+A1 A2 \[Lambda]^4))]]/(2Sqrt[alpha]M^-1 Sqrt[A1 \[Lambda]^2+B1 \[Beta]^2] Sqrt[A2 \[Lambda]^2+B2 \[Beta]^2]);
Ctot=2CI-Cfp;
{um,up}=horizons[M];
up=1/up;
vp=uc-up;
Cout=2CI-Cfp-k/(Sqrt[alpha]M^-1(\[Lambda]-\[Beta])\[Beta] Sqrt[A2 B1]) EllipticF[ArcCos[Sqrt[-A1/B1]Abs[(vp-\[Lambda])/(vp-\[Beta])]],k^2];
Cout=Cout+(k n \[Lambda] A1)/(Sqrt[alpha]M^-1 \[Beta]^3 B1 Sqrt[A2 B1]) myEllipticPi[n,ArcCos[Sqrt[-A1/B1] (vp-\[Lambda])/(vp-\[Beta])],k^2]-1/(2Sqrt[alpha]M^-1) L[(vp-\[Lambda])/(vp-\[Beta])];
Cout=-Cout+2Ctot;
Ctot=2Ctot;
slopvout=-(1/(Sqrt[alpha]M^-1 Sqrt[A1 \[Lambda]^2+B1 \[Beta]^2] Sqrt[A2 \[Lambda]^2+B2 \[Beta]^2]));
{{Cout,Ctot},slopvout}
]


criPhisofb[M_]:=Block[{plot,data,fit1,fit2,xrange,b1,b2,b1p},
plot=Plot[phiout[M,b ],{b, 10^-11,bcri[M]-10^-11},PlotRange->All];
data=Cases[plot,_Line,Infinity];
data=List@@@data//Flatten[#,1]&;
fit1=Interpolation[data[[1]]];
fit2=Interpolation[data[[2]]];
xrange=fit1[[1,1]];
b1=Table[bisection[fit1[#]-Pi/2-n Pi&,xrange[[1]],xrange[[2]],10^-8],{n,0,4}];
b1p=Table[bisection[fit2[#]-Pi/2-n Pi&,xrange[[1]],xrange[[2]],10^-8],{n,0,4}];
b2=bisection[fit2[#]-fit1[#]-Pi&,xrange[[1]],xrange[[2]],10^-8];
{b1p,b1,b2}
]


criPhisofb[M_,n1_,n2_]:=Block[{x,plot,data,fit1,fit2,xrange,b1,b2,b1p},
plot=Plot[phiout[M,b ],{b, 10^-11,bcri[M]-10^-5},PlotRange->All];
data=Cases[plot,_Line,Infinity];
data=List@@@data//Flatten[#,1]&;
fit1=Interpolation[data[[1]]];
fit2=Interpolation[data[[2]]];
xrange=fit1[[1,1]];
b1=Table[bisection[fit1[#]-Pi/2-n Pi&,xrange[[1]],xrange[[2]],10^-8],{n,n1,n2}];
b1p=Table[bisection[fit2[#]-Pi/2-n Pi&,xrange[[1]],xrange[[2]],10^-8],{n,n1,n2}];
b2=bisection[fit2[#]-fit1[#]-Pi&,xrange[[1]],xrange[[2]],10^-8];
{b1p,b1,b2}
]


Lambda[\[Lambda]_]:=-1-(6 (-9+Sqrt[81-3 \[Lambda]^2]))/\[Lambda]^2


phitotLargeM[lam_,M_]:=Block[{coe,ang,term1,term2,adb,kk,coe2},
coe=Sqrt[alpha]/2 (2/alpha)^(1/3);
term1=0.5063315995556266 M^(-1/3);
adb=2 Lambda[lam]^(2/3.)-Sqrt[3(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+Lambda[lam]^(1/3.)+2;
adb=-adb/((1+Lambda[lam]^(1/3.))^2);
kk=1/4 ((Sqrt[3.](1+Lambda[lam]^(2/3)))/Sqrt[(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+2);
coe2=(2^(2/3.) 3^(1/4.) (Lambda[lam]/alpha)^(1/6))/(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)^(1/4.);
term2=coe2 EllipticF[ArcCos[adb],kk]-2.M^(-1/3);
2coe(term1+term2)
]



phitotLargeM[lam_,Infinity]:=Block[{coe,ang,term1,term2,adb,kk,coe2},
coe=Sqrt[alpha]/2 (2/alpha)^(1/3);
term1=0;
adb=2 Lambda[lam]^(2/3.)-Sqrt[3(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+Lambda[lam]^(1/3.)+2;
adb=-adb/((1+Lambda[lam]^(1/3.))^2);
kk=1/4 ((Sqrt[3.](1+Lambda[lam]^(2/3)))/Sqrt[(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+2);
coe2=(2^(2/3.) 3^(1/4.) (Lambda[lam]/alpha)^(1/6))/(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)^(1/4.);
term2=coe2 EllipticF[ArcCos[adb],kk]-0;
2coe(term1+term2)
]


deltaphiLargeM[lam_,Infinity]:=Block[{coe,term1,adb,kk,adbH},
coe=Sqrt[alpha]/2 (2/alpha)^(1/3);
coe=coe (2^(2/3.) 3^(1/4.) (Lambda[lam]/alpha)^(1/6))/(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)^(1/4.);
adb=2 Lambda[lam]^(2/3.)-Sqrt[3(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+Lambda[lam]^(1/3.)+2;
adb=-adb/((1+Lambda[lam]^(1/3.))^2);
kk=1/4 ((Sqrt[3.](1+Lambda[lam]^(2/3)))/Sqrt[(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+2);
term1=EllipticF[ArcCos[adb],kk];
adbH=1-(2Sqrt[3] Sqrt[(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)])/(Lambda[lam]^(2/3)+Sqrt[3] Sqrt[(Lambda[lam]^(2/3.)+Lambda[lam]^(4/3)+1)]+2Lambda[lam]^(1/3)+1);
term1=term1-EllipticF[ArcCos[adbH],kk];
coe term1
]


phioutLargeM[lam_,M_]:=phitotLargeM[lam,M]-deltaphiLargeM[lam,Infinity]
