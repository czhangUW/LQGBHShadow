(* ::Package:: *)

perturbateroot[pol_,perb_,r_,r0_,ep_,degroot_,n_:7]:=Block[{plp,v,fp,coe,eq,coev},
plp=pol/.r->r0+v//CoefficientList[#,v]&;
plp=plp[[degroot+1;;-1]].(v^#&/@Range[0,Length[plp]-degroot-1 ]);
plp=v^degroot plp;
fp=plp+perb/.v->(coe[#]ep^#&/@Range[n]//Total);
fp=Series[fp,{ep,0,n-2}]//Normal//Expand;
eq=CoefficientList[fp,ep][[degroot+1;;-1]]==0//Thread;
coev=coe/@Range[1,Length[eq]];
coev=coev/.NSolve[eq,coev];
coev=SortBy[coev,First];
coev=#.(ep^#&/@Range[1,Length[eq]])&/@coev;
coev=coev+r0
]


geo[b_,M_]/;Abs[G[umin[M],b,M]]>= 10^-7:=Block[{uroot,u,tphi,sol,x},
uroot=turning[b,M];
tphi=NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,0,uroot},MaxRecursion->12];
sol=NDSolve[{D[u[x],x]==Re@Sqrt[G[u[x],b,M]],u[0]==0},u,{x,0,tphi}];
{u/.sol[[1]],tphi,uroot}]
(***************************************************)
geo[b_,M_]/;G[umin[M],b,M]>10^-16&&G[umin[M],b,M]<10^-7:=Block[
{uroot,uc,bc,u1,ep,u2,u,tphi,Gc,v,y,sep,rest,sol},
uroot=turning[b,M];
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
sol=NDSolve[{D[u[x],x]==Re@Sqrt[G[u[x],b,M]],u[0]==0},u,{x,0,tphi}];
{u/.sol[[1]],tphi,uroot}]
(**************************************************)
geo[b_,M_]/;(G[umin[M],b,M]<-10^-16.&&G[umin[M],b,M]>= -10^-7):=Block[
{uroot,bc,ep,u,tphi,uc,Gc,v,y,sep,urootall,ucp,urootallp,rest1,err,n,rest,sol,x},
uroot=turning[b,M];
bc=bcri[M];
ep=Sqrt[b-bc]/Sqrt[bc^3];
tphi=NIntegrate[1/Re@Sqrt[G[u,b,M]],{u,0,uroot-ep},MaxRecursion->12];
uc=umin[M];
Gc=G[v+uc,bc,M]//CoefficientList[#,v]&;
Gc=Gc[[3;;-1]].{1,v,v^2,v^3,v^4};
urootall=v/.NSolve[Gc==0,v];
urootall=uc+urootall;
ucp=perturbateroot[G[u,bc,M],-1/bc^2+1/(bc+sep^2)^2,u,uc,sep,2];
urootallp=perturbateroot[G[u,bc,M],-1/bc^2+1/(bc+bc^3sep^2)^2,u,#,sep,1]&/@urootall//Flatten;
Gc=(Times@@(ucp[[1]]-#-y sep&/@urootallp))(-y sep)(ucp[[1]]-ucp[[2]]-y sep)//Expand;
Gc=Gc Coefficient[G[u,bc,M],u^6];
Gc=CoefficientList[Gc,sep][[3;;-1]];
Gc=Gc.(sep^#&/@Range[0,Length[Gc]-1 ]);
Gc=Series[1/Sqrt[Gc],{sep,0,3}]//Normal;
Gc=Gc/.sep->ep;
rest=NIntegrate[Gc,{y,0,1}]//Re;
tphi=tphi+rest;
sol=NDSolve[{D[u[x],x]==Re@Sqrt[G[u[x],b,M]],u[0]==0},u,{x,0,tphi}];
{u/.sol[[1]],tphi,uroot}
]
geo[b_,beta_]/;Abs[G[umin[beta],b,beta]]<= 10^-16:={{},10000Pi,umin[beta]}
