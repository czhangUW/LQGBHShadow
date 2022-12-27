(* ::Package:: *)

maTexTicks[plotrangep_,xscale_,yscale_,tickslen_,ifxticks_,ifyticks_:True,fontsize_:1]:=Block[
{xticks,yticks,posx,posy,xticksp,yticksp,xtickslen,ytickslen},
xticks=Charting`ScaledTicks[xscale]@@plotrangep[[1]]//Transpose;
yticks=Charting`ScaledTicks[yscale]@@plotrangep[[2]]//Transpose;
posx=Position[xticks[[2]],Except[_Spacer],{1},Heads->False];
posy=Position[yticks[[2]],Except[_Spacer],{1},Heads->False];
xticksp=If[ifxticks,MapAt[MaTeX[ToString[#],Magnification->fontsize]&,xticks[[2]],posx],ReplacePart[xticks[[2]],posx->Spacer[List[0,0]]]];
yticksp=If[ifyticks,MapAt[MaTeX[ToString[#],Magnification->fontsize]&,yticks[[2]],posy],ReplacePart[yticks[[2]],posx->Spacer[List[0,0]]]];
xtickslen=Map[ReplacePart[#,1->tickslen  #[[1]]]&,xticks[[3]]];
ytickslen=Map[ReplacePart[#,1->tickslen #[[1]]]&,yticks[[3]]];
xticks=MapAt[xscale[[2]],xticks,{1}];
yticks=MapAt[yscale[[2]],yticks,{1}];
xticks=ReplacePart[xticks,{2->xticksp,3->xtickslen}]//Transpose;
yticks=ReplacePart[yticks,{2->yticksp,3->ytickslen}]//Transpose;
{xticks,yticks}
]


change2maTexTicks[oxticks_,oyticks_,tickslen_,fontsize_:1]:=Block[
{xticks,yticks,posx,posy,xticksp,yticksp,xtickslen,ytickslen},
xticks=oxticks//Transpose;
yticks=oyticks//Transpose;
posx=Position[xticks[[2]],Except[_Spacer],{1},Heads->False];
posy=Position[yticks[[2]],Except[_Spacer],{1},Heads->False];
xticksp=MapAt[MaTeX[ToString[#],Magnification->fontsize]&,xticks[[2]],posx];
yticksp=MapAt[MaTeX[ToString[#],Magnification->fontsize]&,yticks[[2]],posy];
xtickslen=Map[ReplacePart[#,1->tickslen  #[[1]]]&,xticks[[3]]];
ytickslen=Map[ReplacePart[#,1->tickslen #[[1]]]&,yticks[[3]]];
xticks=ReplacePart[xticks,{2->xticksp,3->xtickslen}]//Transpose;
yticks=ReplacePart[yticks,{2->yticksp,3->ytickslen}]//Transpose;
{xticks,yticks}
]



matexrules={Superscript[10,x_]:>MaTeX["10^{"<>ToString[x]<>"}"],
x_?NumberQ/;(Log10[x]==Floor[Log10[x]]&&x!=10&&x!=1):>MaTeX[ "10^{"<>ToString[Log10[x]//Floor]<>"}" ],
10->MaTeX[10],
1->MaTeX[1],
x_?NumericQ/;Log10[x]!=Floor[Log10[x]]:>MaTeX[ToString[x]],
Row[{x_,Superscript[10,y_]},_]/;x!=1:>Spacer[{0,0}],
Row[{x_,Superscript[10,y_]},_]/;x==1:>MaTeX["10^{"<>ToString[y]<>"}"],
NumberForm[x_,__]/;Log10[x]==Floor[Log10[x]]:>MaTeX[ "10^{"<>ToString[Log10[x]//Floor]<>"}" ],
(y:NumberForm[x_,__])/;Log10[x]!=Floor[Log10[x]]:>Spacer[{0,0}]
};


myplot[func_,xlabel_,ylabel_,lt_,rt_,plotrange_,tickslen_,ifxticks_,ifframe_,style_:Automatic,legend_:{}]:=Block[
{plot,plotrangep,xticks,yticks},
plot=Plot[func[t],{t,lt,rt},PlotRange->plotrange,Frame->ifframe,PlotStyle->style,PlotLegends->legend];
plotrangep=AbsoluteOptions[plot,PlotRange][[1,2]];
{xticks,yticks}=maTexTicks[plotrangep,{Identity,Identity},{Identity,Identity},tickslen,ifxticks];
If[ifframe, plot/.{Rule[FrameTicks,_]:>FrameTicks->{{yticks,Automatic},{xticks,Automatic}},Rule[FrameLabel,_]:>Rule[FrameLabel,{xlabel,ylabel}]},
plot/.{(Ticks->_):>Ticks->{xticks,yticks},Rule[AxesLabel,_]:>Rule[AxesLabel,{xlabel,ylabel}]}
]
]


myplotScientificNotionY[func_,xlabel_,ylabel_,lt_,rt_,plotrange_,tickslen_,ifxticks:True,style_:Automatic]:=Block[
{plot,plotrangep,xticks,yticks,posx,posy,xticksp,commonpower,yticksp,xtickslen,ytickslen},
plot=Plot[func,{t,lt,rt},PlotRange->plotrange,Frame->True];
plotrangep=AbsoluteOptions[plot,PlotRange][[1,2]];
xticks=Charting`ScaledTicks[{Identity,Identity}]@@plotrangep[[1]]//Transpose;
yticks=Charting`ScaledTicks[{Identity,Identity}]@@plotrangep[[2]]//Transpose;
posx=Position[xticks[[2]],Except[_Spacer],{1},Heads->False];
posy=Position[yticks[[2]],Except[_Spacer],{1},Heads->False];
xticksp=MapAt[MaTeX[ToString[#]]&,xticks[[2]],posx];
commonpower=DeleteCases[Log10[yticks[[1]][[#]]&/@posy]//Floor//Flatten,Indeterminate]//Min;
yticksp=Thread[Rule[posy//Flatten,10^(-commonpower)(yticks[[1]][[#]]&/@posy)//Flatten]];
yticksp=ReplacePart[yticks[[2]],yticksp];
yticksp=MapAt[MaTeX,yticksp,posy];
xtickslen=Map[ReplacePart[#,1->tickslen  #[[1]]]&,xticks[[3]]];
ytickslen=Map[ReplacePart[#,1->tickslen #[[1]]]&,yticks[[3]]];
xticks=ReplacePart[xticks,{2->xticksp,3->xtickslen}]//Transpose;
yticks=ReplacePart[yticks,{2->yticksp,3->ytickslen}]//Transpose;
plot=Plot[func,{t,lt,rt},PlotRange->plotrange,Frame->True,Axes->False,FrameLabel->{{ylabel,None},{xlabel,None}},FrameTicks->{{yticks,Automatic},{xticks,Automatic}},PlotStyle->style];
Labeled[plot,MaTeX["\\quad\\times 10^{"<>ToString[commonpower]<>"}"],{{Top,Left}},Spacings->{ 0,-0.3}]
]
