(* ::Package:: *)

\[Gamma]=2375/10000;alpha=16Sqrt[3]\[Pi] \[Gamma]^3;
(***Here we choose M as unit***)
G[u_,b_,M_]:=-alpha M^-2 u^6+2u^3-u^2+b^-2
horizons[M_]:=Block[{f,r},
f=r/.NSolve[1-2/r+alpha M^-2 1/r^4==0,r];
f=Cases[f,_?(Element[#,Reals]&)];
f=Sort[f]
]


bcri[M_]:=Block[{DG,u,uv},
DG=D[-alpha M^-2 u^6+2u^3-u^2,u];
uv=u/.Solve[DG==0,u];
uv=Cases[uv,_?(Element[#,Reals]&&#>0&)]//Sort//First;
1/Sqrt[alpha M^-2 uv^6-2uv^3+uv^2]
]
umin[M_]:=Block[{u,uv},
uv=u/.NSolve[D[-alpha M^-2 u^6+2u^3-u^2,u]==0,u];
uv=Cases[uv,_?(Element[#,Reals]&&#>0&)]//Sort//First
]

turning[b_,M_]:=Block[{u,uroot},
uroot=u/.NSolve[G[u,b,M]==0,u];
Cases[uroot,_?(Element[#,Reals]&&#>0&)]//Sort//First]


bisection[f_,x0_,x1_,ep_]:=Block[{err,y0,y1,root},
err=1;
y0=x0;
y1=x1;
While[Abs[err]>ep||err>0,err=f[0.5(y0+y1)];root=0.5(y0+y1);If[err f[y0]>0,y0=root,y1=root]];
root
]
