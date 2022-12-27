(* ::Package:: *)

(* ::Input:: *)
(**)


G1[w_,bv_,betav_]:=Re@Sqrt[G[w,bv,betav]](1+bv Re@Sqrt[G[w,bv,betav]])/bv


uvpr[r_,beta_,alpha_:1.1663]:=Block[{M,rp,rm,b0,ap,b1,am,W0,Wr},
M=Sqrt[(4 (beta^4) )/(1-beta^2)^3 alpha];
rp=M/2 (1+1/beta)(1+Sqrt[2beta-1]);
rm=M/2 (1+1/beta)(1-Sqrt[2beta-1]);
b0=( M^2 (-2+beta) (-1+beta)^3 (1+beta))/(4 beta^2 (2+beta^2));
ap=(M (1+beta)^2 (-1+Sqrt[-1+2 beta]+beta (3+beta+2 Sqrt[-1+2beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
b1=(M (-1+beta)^2 (-1+2 beta))/(2 beta (2+beta^2));
am=(M (1+beta)^2 (-1-Sqrt[-1+2 beta]+beta (3+beta-2 Sqrt[-1+2 beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
W0=((1-beta)^2 (beta+1))/(2beta^2) M^2;
Wr=r^2+(1/beta-1) M r+((1-beta)^2 (beta+1))/(2beta^2) M^2;
-Exp[r/ap-(2b0/( M)+b1(1-beta))/(ap(1-beta)Sqrt[1+2beta]) (ArcTan[(1-beta+2 beta r/( M))/((1-beta)Sqrt[2beta+1])]-ArcTan[1/Sqrt[1+2beta]])](r/rp-1)(r/rm-1)^(-am/ap) (Wr/W0)^(b1/(2ap))
]
coorv[r_,vpp_,beta_]:={1/2 (vpp-ArcTan[uvpr[r,beta]/Tan[vpp]]),1/2 (vpp+ArcTan[uvpr[r,beta]/Tan[vpp]])}
cooru[r_,upp_,beta_]:={1/2 (ArcTan[uvpr[r,beta]/Tan[upp]]-upp),1/2 (upp+ArcTan[uvpr[r,beta]/Tan[upp]])}


tildeuvpr[r_,beta_,alpha_:1.1663]:=Block[{M,rp,rm,b0,ap,b1,am,W0,Wr},
M=Sqrt[(4 (beta^4) )/(1-beta^2)^3 alpha];
rp=M/2 (1+1/beta)(1+Sqrt[2beta-1]);
rm=M/2 (1+1/beta)(1-Sqrt[2beta-1]);
b0=( M^2 (-2+beta) (-1+beta)^3 (1+beta))/(4 beta^2 (2+beta^2));
ap=(M (1+beta)^2 (-1+Sqrt[-1+2 beta]+beta (3+beta+2 Sqrt[-1+2beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
b1=(M (-1+beta)^2 (-1+2 beta))/(2 beta (2+beta^2));
am=(M (1+beta)^2 (-1-Sqrt[-1+2 beta]+beta (3+beta-2 Sqrt[-1+2 beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
W0=((1-beta)^2 (beta+1))/(2beta^2) M^2;
Wr=r^2+(1/beta-1) M r+((1-beta)^2 (beta+1))/(2beta^2) M^2;
Exp[r/-am-(2beta b0/( M)+b1(1-beta))/(-am(1-beta)Sqrt[2beta+1]) (ArcTan[(1-beta+2beta r/( M))/((1-beta)Sqrt[2beta+1])]-ArcTan[1/Sqrt[1+2beta]])](Wr/W0)^(b1/(-2am)) (1-r/rp)^(-ap/am) (r/rm-1)
]


myTh[x_]/;x!=0:=HeavisideTheta[x]
myTh[0]:=1

coortildv[r_,vpp_,beta_,alpha_:1.1663]/;vpp!=Pi/2:=Block[{M,ap,am,sgntildeup,sgntildevp,tildeup,upp},
M=Sqrt[(4 (beta^4) )/(1-beta^2)^3 alpha];
ap=(M (1+beta)^2 (-1+Sqrt[-1+2 beta]+beta (3+beta+2 Sqrt[-1+2beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
am=(M (1+beta)^2 (-1-Sqrt[-1+2 beta]+beta (3+beta-2 Sqrt[-1+2 beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
sgntildevp=Sign[vpp-Pi/2];
sgntildeup=Sign[tildeuvpr[r,beta]]sgntildevp;
upp=\[Pi] myTh[sgntildeup]-sgntildeup ArcTan[Abs[tildeuvpr[r,beta]]^(-am/ap)/ Tan[sgntildevp (myTh[sgntildevp] Pi-vpp)]];
{1/2 (vpp-upp),1/2 (vpp+upp )}]

coortildu[r_,upp_,beta_,alpha_:1.1663]/;upp!=Pi/2:=Block[{M,am,ap,sgntildeup,sgntildevp,tildeup,vpp},
M=Sqrt[(4 (beta^4) )/(1-beta^2)^3 alpha];
ap=(M (1+beta)^2 (-1+Sqrt[-1+2 beta]+beta (3+beta+2 Sqrt[-1+2beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
am=(M (1+beta)^2 (-1-Sqrt[-1+2 beta]+beta (3+beta-2 Sqrt[-1+2 beta])))/(4 beta Sqrt[-1+2 beta] (2+beta^2));
sgntildeup=Sign[upp-Pi/2];
sgntildevp=Sign[tildeuvpr[r,beta]]sgntildeup;
vpp=\[Pi] myTh[sgntildevp]-sgntildevp ArcTan[Abs[tildeuvpr[r,beta]]^(-am/ap)/ Tan[sgntildeup (myTh[sgntildeup] Pi-upp)]];
{1/2 (vpp-upp),1/2 (vpp+upp )}]


coorbarv[r_,vpp_,beta_]:={1/2 (vpp-ArcTan[uvpr[r,beta]/Tan[vpp-Pi]]-Pi),1/2 (vpp+ArcTan[uvpr[r,beta]/Tan[vpp-Pi]]+Pi)}
coorbaru[r_,upp_,beta_]:={1/2 (Pi+ArcTan[uvpr[r,beta]/Tan[upp-Pi]]-upp),1/2 (upp+ArcTan[uvpr[r,beta]/Tan[upp-Pi]]+Pi)}


geoInPenrose[bv_,betav_,iniv_:0,alpha_:1.1663]/;bv<bcri[betav]:=Block[{wroot,w,tv,tvu,rm,rp,vpp,vp,up,sol,geowv,M,ap,am,Bouncetildevp,Bouncetildup,u0,tu,geowu},
M=Sqrt[(4 (betav^4) )/(1-betav^2)^3 alpha];
ap=(M (1+betav)^2 (-1+Sqrt[-1+2 betav]+betav (3+betav+2 Sqrt[-1+2betav])))/(4 betav Sqrt[-1+2 betav] (2+betav^2));
am=(M (1+betav)^2 (-1-Sqrt[-1+2 betav]+betav (3+betav-2 Sqrt[-1+2 betav])))/(4 betav Sqrt[-1+2 betav] (2+betav^2));
{rm,rp}=horizons[betav];
vpp=NIntegrate[G1[w,bv,betav]^-1,{w,0,1/rp}];
vp=iniv+vpp;(*the vaule of v at horizon rp*)
wroot=w/.NSolve[G[w,bv,betav]==0,w];
wroot=Cases[wroot,_?(Element[#,Reals]&&#>0&)]//Sort//First;
tvu=NIntegrate[G1[w,bv,betav]^-1 ,{w,0,wroot}];
tv=tvu+iniv;(*the value of v at bounce*)
sol=NDSolve[{D[w[v],v]==G1[w[v],bv,betav],w[iniv]==0},w,{v,iniv,tv}];
geowv=w/.sol[[1]];
Bouncetildevp=-Exp[-(1/(2am)) tv];
Bouncetildup=tildeuvpr[1/geowv[tv],betav]/Bouncetildevp;
u0=2 am Log[ Bouncetildup];
tu=u0+tvu;
up=tu-vpp(*value of u at rp*);
sol=NDSolve[{D[w[u],u]==-G1[w[u],bv,betav],w[tu]==0},w,{u,u0,tu}];
geowu=w/.sol[[1]];
{{geowv,{iniv,tv}},{geowu,{u0,tu}},vp,up,betav}
]


geoInPenrose[bv_,betav_,iniv_:0,alpha_:1.1663]/;bv>bcri[betav]:=Block[{wroot,w,tv,tvu,sol,geowv,M,ap,am,Bouncebarvp,Bouncebarup,u0,tu,geowu},
M=Sqrt[(4 (betav^4) )/(1-betav^2)^3 alpha];
ap=(M (1+betav)^2 (-1+Sqrt[-1+2 betav]+betav (3+betav+2 Sqrt[-1+2betav])))/(4 betav Sqrt[-1+2 betav] (2+betav^2));
am=(M (1+betav)^2 (-1-Sqrt[-1+2 betav]+betav (3+betav-2 Sqrt[-1+2 betav])))/(4 betav Sqrt[-1+2 betav] (2+betav^2));
wroot=w/.NSolve[G[w,bv,betav]==0,w];
wroot=Cases[wroot,_?(Element[#,Reals]&&#>0&)]//Sort//First;
tvu=NIntegrate[G1[w,bv,betav]^-1,{w,0,wroot}];
tv=tvu+iniv;(*the value of v at bounce*)
sol=NDSolve[{D[w[v],v]==G1[w[v],bv,betav],w[iniv]==0},w,{v,iniv,tv}];
geowv=w/.sol[[1]];
Bouncebarvp=Exp[1/(2ap) tv];
Bouncebarup=uvpr[1/geowv[tv],betav]/Bouncebarvp;
u0=-2 ap Log[- Bouncebarup];
tu=u0+tvu;
sol=NDSolve[{D[w[u],u]==-G1[w[u],bv,betav],w[tu]==0},w,{u,u0,tu}];
geowu=w/.sol[[1]];
{{geowv,{iniv,tv}},{geowu,{u0,tu}},{},{},betav}
]
