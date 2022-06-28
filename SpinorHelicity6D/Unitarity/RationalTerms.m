(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: RationalTerms *)
(* :Context: RationalTerms` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`Unitarity`RationalTerms`",
  {"RationalSeries`","SpinorHelicity6D`","SpinorHelicity6D`Unitarity`","Utils`"}];


GetBoxRationalTerm::usage = "";
GetTriangleRationalTerm::usage = "";
GetBubbleRationalTerm::usage = "";
GetBubbleBubRationalTerm::usage = "";
GetBubbleBubRationalTermOld::usage = "";
GetBubbleTriRationalTerm::usage = "";
GetBubbleTriRationalTermOld::usage = "";

GetBoxCoefficient::usage = "";
GetTriangleCoefficient::usage = "";
GetBubbleBubCoefficient::usage = "";
GetBubbleTriCoefficient::usage = "";


Begin["`Private`"];


(* ::Section:: *)
(* Rational Terms *)


(* ::Subsection:: *)
(* Boxes *)


Options[GetBoxRationalTerm]={MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,LargeMuLimit->True,Assumptions->{},AnalyticVariables -> {},Gravity->False};
GetBoxRationalTerm[boxTerm_, opts : OptionsPattern[]] := Module[{res,muOrders},
  If[OptionValue[LargeMuLimit],
    (* If LargeMuLimit is true, we make the implicit redefinition that MuVar is Sqrt[MuVar], etc. to
      not have to deal with square roots *)
    If[OptionValue[Gravity],
      muOrders = {2,3,4};,
      muOrders = {2};
    ];
    res = Table[RationalSeriesCoefficient[
      RationalSeriesCoefficient[
        boxTerm,
        OptionValue[MuVar],-2*muOrder,AnalyticVariables->Join[{OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
      OptionValue[MuTildeVar],-2*muOrder,AnalyticVariables->Join[OptionValue[AnalyticVariables]]],
      {muOrder,muOrders}];
    (*I/2*SeriesCoefficient[boxTerm/.{OptionValue[MuVar]->1/OptionValue[MuVar],OptionValue[MuTildeVar]->1/OptionValue[MuTildeVar]},
    {OptionValue[MuVar],0,-4},{OptionValue[MuTildeVar],0,-4},Assumptions->OptionValue[Assumptions]]*)
  ];
  1/2*res
];


Options[GetBoxCoefficient]={Gravity->False};
GetBoxCoefficient[boxTerm_, opts : OptionsPattern[]] := Module[{},
  1/2*boxTerm
];


(* ::Subsection:: *)
(* Triangles *)


Options[GetTriangleRationalTerm]={MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,tVar->$LoopMomentumT,Assumptions->{},AnalyticVariables -> {},Gravity->False};

GetTriangleRationalTerm[triangleTerm_, OptionsPattern[]] := Module[{res,muOrders},
  If[OptionValue[Gravity],
    muOrders = {1,2,3};,
    muOrders = {1};
  ];
  res = Table[RationalSeriesCoefficient[
    RationalSeriesCoefficient[
      RationalSeriesCoefficient[
        triangleTerm,
        OptionValue[tVar],0,AnalyticVariables->Join[{OptionValue[MuVar],OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
      OptionValue[MuVar],-muOrder,AnalyticVariables->Join[{OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
    OptionValue[MuTildeVar],-muOrder,AnalyticVariables->Join[OptionValue[AnalyticVariables]]],
    {muOrder,muOrders}];
  1/2*res
  (*-1/2*SeriesCoefficient[(triangleTerm/.{OptionValue[tVar]->1/OptionValue[tVar],OptionValue[MuVar]->1/OptionValue[MuVar],OptionValue[MuTildeVar]->1/OptionValue[MuTildeVar]}),
 {OptionValue[tVar],0,0},{OptionValue[MuVar],0,-1},{OptionValue[MuTildeVar],0,-1},Assumptions->OptionValue[Assumptions]]*)
];


Options[GetTriangleCoefficient]={tVar->t,Assumptions->{},AnalyticVariables -> {},Gravity->False};

GetTriangleCoefficient[triangleTerm_, OptionsPattern[]] := Module[{},
  1/2*RationalSeriesCoefficient[
    triangleTerm,
    OptionValue[tVar],0,AnalyticVariables->OptionValue[AnalyticVariables]]
];


(* ::Subsection:: *)
(* Bubbles *)


(* ::Subsubsection:: *)
(* Bubble Bub*)


Options[GetBubbleBubRationalTerm]=
    {MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,tVar->$LoopMomentumT,yVar->$LoopMomentumY,Assumptions->{},AnalyticVariables->{},Gravity->False};

GetBubbleBubRationalTerm[bubbleBubTerm_,Yi_, opts : OptionsPattern[]]:=GetBubbleBubRationalTerm[bubbleBubTerm,Yi,0,opts];

GetBubbleBubRationalTerm[bubbleBubTerm_,Yi_,MuValue_?NumericQ, opts : OptionsPattern[]] := Module[{res,resComp,yOrders,muOrders},
  If[OptionValue[Gravity],
    yOrders = Range[0,-4,-1];
    muOrders = {1,2};,
    yOrders = Range[0,-2,-1];
    muOrders = {1};
  ];
  res = Table[RationalSeriesCoefficient[
    RationalSeriesCoefficient[
      RationalSeriesCoefficient[Yi .
          RationalSeriesCoefficient[
            bubbleBubTerm,
            OptionValue[yVar],yOrders,AnalyticVariables->Join[{OptionValue[tVar],OptionValue[MuVar],OptionValue[MuTildeVar]},OptionValue[AnalyticVariables],{}]],
        OptionValue[tVar],0,AnalyticVariables->Join[{OptionValue[MuVar],OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
      OptionValue[MuVar],-muOrder,AnalyticVariables->Join[{OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
    OptionValue[MuTildeVar],-muOrder,AnalyticVariables->Join[OptionValue[AnalyticVariables]]],
    {muOrder,muOrders}];
  res
];




Options[GetBubbleBubCoefficient]={tVar->$LoopMomentumT,yVar->$LoopMomentumY,Assumptions->{},AnalyticVariables->{},Gravity->False};

GetBubbleBubCoefficient[bubbleBubTerm_,Yi_, opts : OptionsPattern[]] := Module[{res,yOrders,muOrders},
  If[OptionValue[Gravity],
    yOrders = Range[0,-4,-1];,
    yOrders = Range[0,-2,-1];
  ];
  res = RationalSeriesCoefficient[Yi .
      RationalSeriesCoefficient[
        bubbleBubTerm,
        OptionValue[yVar],yOrders,AnalyticVariables->Join[{OptionValue[tVar]},OptionValue[AnalyticVariables]]],
    OptionValue[tVar],0,AnalyticVariables->OptionValue[AnalyticVariables]];
  res
];


(* ::Subsubsection:: *)
(*Bubble Triangle*)


Options[GetBubbleTriRationalTerm]={MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,tVar->t,Assumptions->{},AnalyticVariables->{},Gravity->False};

GetBubbleTriRationalTerm[expr_,ti_List, opts : OptionsPattern[]]:= Module[{res,tOrders,muOrders},
  If[OptionValue[Gravity],
    tOrders = Range[-1,-6,-1];
    muOrders = {1,2};,
    tOrders = Range[-1,-3,-1];
    muOrders = {1};
  ];
  res = Table[RationalSeriesCoefficient[
    RationalSeriesCoefficient[ti .
        RationalSeriesCoefficient[
          expr,
          OptionValue[tVar],tOrders,AnalyticVariables->Join[{OptionValue[MuVar],OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
      OptionValue[MuVar],-muOrder,AnalyticVariables->Join[{OptionValue[MuTildeVar]},OptionValue[AnalyticVariables]]],
    OptionValue[MuTildeVar],-muOrder,AnalyticVariables->Join[OptionValue[AnalyticVariables]]],
    {muOrder,muOrders}];
  1/2*res
];

GetBubbleTriRationalTerm[exprs_List,ti_List, opts : OptionsPattern[]] := GetBubbleTriRationalTerm[#,ti,opts]&/@exprs;


Options[GetBubbleTriCoefficient]={tVar->t,Assumptions->{},AnalyticVariables->{},Gravity->False};

GetBubbleTriCoefficient[expr_,ti_List, opts : OptionsPattern[]]:= Module[{res,tOrders},
  If[OptionValue[Gravity],
    tOrders = Range[-1,-6,-1];,
    tOrders = Range[-1,-3,-1];
  ];
  res = ti . RationalSeriesCoefficient[
    expr,
    OptionValue[tVar],tOrders,AnalyticVariables->OptionValue[AnalyticVariables]];
  1/2*res
];

GetBubbleTriCoefficient[exprs_List,ti_List, opts : OptionsPattern[]] := GetBubbleTriRationalTerm[#,ti,opts]&/@exprs;



End[]; (* `Private` *)

EndPackage[]