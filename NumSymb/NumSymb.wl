(*

    This file is part of SpinorHelicityAddons.

    SpinorHelicityAddons is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SpinorHelicityAddons is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SpinorHelicityAddons.  If not, see <http://www.gnu.org/licenses/>.

*)

(* ::Package:: *)

(* ::Title:: *)
(*NumSymb *)


BeginPackage["NumSymb`", "Utils`", "RationalSeries`"];

Unprotect["NumSymb`*"];
ClearAll["NumSymb`*"];
ClearAll["NumSymb`*"];

NumSymb::usage = "Function applied to unique sumbol that stands for a numerical value.";
NumSymbN::usage = "When applied to unique symbol of NumSymb, returns the numerical value associated to it.";
NumSymbEvaluated::usage = "Fully numerically evaluated NumSymb objects, where remaining variables are replaced with their hashes.";
NumSymbSimplified::usage = "Analytically simplified NumSymb.";
NumSymbLookUp::usage = "Look up table for NumSymb.";
NumSymbReplace::usage = "Replaces numeric coefficients with NumSymb objects.";
CoalesceNumSymb::usage = "Coalesces products (not yet sums) of NumSymb into a single one.";
NumSymbLookUp = <||>;

GetAllNumSymbs::usage = "";

SimplifyNumSymbs::usage = "";
SimplifyNumSymbExpression::usage = "";

$LeafCountLimit::usage = "";
$LeafCountLimit = 100;
sqrt::usage = "Sqrt placeholder during NumSymb simplifications.";


Begin["`Private`"];


(* ::Chapter:: *)
(*Main*)


(* ::Section:: *)
(*Creating NumSymbs*)


SetAttributes[NumSymb,{NumericFunction,Constant}];
SetAttributes[NumSymbN,{NumericFunction,Constant}];


genUnique[expr_]:=Module[{unique,possibleZero},
	If[KeyExistsQ[NumSymbLookUp,expr],
		(*Print["Found entry for ", expr];*)
		NumSymbLookUp[expr]
	,
		unique = Unique[];
		EvaluateNumSymb[unique,expr];
		
		(* 
			We are using assignments rather than AssociateTo for NumSymbLookUp, because the latter
			does not play well with being called in subkernels, c.f.
			mathematica.stackexchange.com/questions/140198
		*)
		
		If[PossibleZeroQ[NumSymbEvaluated[unique],Method->"ExactAlgebraics"],
			(*Print["Zero"];*)
			NumSymbEvaluated[unique]=0;
			NumSymbN[unique]=0; 
			NumSymb[unique]=0; 
			NumSymbLookUp[expr] = unique;
			(*AssociateTo[NumSymbLookUp,expr->unique];*)
		,
			(*Print["Nonzero"];*)
			NumSymbN[unique]=expr; 
			NumSymbLookUp[expr] = unique;
			(*AssociateTo[NumSymbLookUp,expr->unique];*)
		];
		unique
	]
];


EvaluateNumSymb[unique_,expr_]:=
	NumSymbEvaluated[unique] = RationalSeries`SampleEval[expr];


NumSymbReplace[expr_List,vars_]:= NumSymbReplace[#,vars]&/@expr;
NumSymbReplace[expr_] := NumSymbReplace[expr,{}];
NumSymbReplace[expr_,vars_] := expr //. {
		Verbatim[Power][a_/;And[!MatchQ[a,NumSymb[_]],AllTrue[vars,FreeQ[a,#]&]],b___] :> Power[NumSymb[genUnique[a]],b],
		Verbatim[Times][c___,a_/;And[!MatchQ[a,NumSymb[_]],AllTrue[vars,FreeQ[a,#]&]],b___] :> Times[c,NumSymb[genUnique[a]],b],
		Verbatim[Plus][c___,a_/;And[!MatchQ[a,NumSymb[_]],AllTrue[vars,FreeQ[a,#]&]],b___] :> Plus[c,NumSymb[genUnique[a]],b]
	};


(*NumSymb /: HoldPattern@Times[x:NumSymb[_],y:NumSymb[_]] := NumSymb[genUnique[x*y]];*)


coalescingRules = {
		(*HoldPattern[x:Rational[__]] \[RuleDelayed] NumSymb[genUnique[x]],*)
		HoldPattern@Plus[x:NumSymb[_],y:NumSymb[_]]:> NumSymb[genUnique[x+y]],
		HoldPattern@Plus[Times[z__,x:NumSymb[_]],Times[z__,y:NumSymb[_]]]  :> Times[z]*NumSymb[genUnique[x+y]],
		HoldPattern@Times[x:NumSymb[_],y:NumSymb[_]] :> NumSymb[genUnique[x*y]],
		HoldPattern@Times[x_?NumericQ,y:NumSymb[_]] /; !MatchQ[x,NumSymb[_]] :> NumSymb[genUnique[x*y]],
		HoldPattern@Times[x:NumSymb[_],Plus[y:NumSymb[_],b__]]:> Plus[NumSymb[genUnique[x*y]],Sequence@@(x*{b})],
		HoldPattern@Times[x:NumSymb[_],Plus[Times[y:NumSymb[_],c__],b__]]:> Plus[Times[NumSymb[genUnique[x*y]],c],Sequence@@(x*{b})],
		HoldPattern@Power[x:NumSymb[_],exp_] :> NumSymb[genUnique[Power[x,exp]]]
		};
	
CoalesceNumSymb[expr_]:= NestWhile[ReplaceRepeated[#,coalescingRules]&,expr,UnsameQ[##]&,2,10];


(* ::Section:: *)
(*NumSymb Simplification*)


GetAllNumSymbs[expr_] /; FreeQ[expr,NumSymb] := {}; 
GetAllNumSymbs[expr_]:=Module[{ns},
	ns = DeleteDuplicates@Cases[expr,NumSymb[_],{0,Infinity}];
	
	Table[{i,GetAllNumSymbs[i/.NumSymb->NumSymbN]},{i,ns}]
]


Options[SimplifyNestedNumSymbs]={Verbose->False}


SimplifyNestedNumSymbs[tree:{node_NumSymb,leafs_List},opts:OptionsPattern[]]:=Module[{},
	SimplifyNestedNumSymbs[#,PassRules[SimplifyNestedNumSymbs,SimplifyNestedNumSymbs,opts]]&/@leafs;
	SimplifyNumSymb[node,PassRules[SimplifyNestedNumSymbs,SimplifyNumSymb,opts]];
	
]
SimplifyNestedNumSymbs[{leaf_NumSymb},opts:OptionsPattern[]]:=Module[{},
	SimplifyNumSymb[leaf,PassRules[SimplifyNestedNumSymbs,SimplifyNumSymb,opts]];
]


Options[SimplifyNumSymb]={Verbose->False}


SimplifyNumSymb[NumSymb[id_],opts:OptionsPattern[]] /; MatchQ[NumSymbSimplified[id],_NumSymbSimplified] := Module[{ns = NumSymb[id],leafcount},
	
	If[OptionValue[Verbose],Print["Looking at \t", NumSymb[id],"..."]];
	
	If[(leafcount=LeafCount[ns//.NumSymb->NumSymbN]) > $LeafCountLimit,
		If[OptionValue[Verbose],Print["Splitting \t", NumSymb[id],", since LeafCount is ", leafcount]];
		ns = NumSymbReplace[ns/.NumSymb->NumSymbN];
		SimplifyNumSymbs[ns,PassRules[SimplifyNumSymb,SimplifyNumSymbs,opts]];
		If[OptionValue[Verbose],Print["Back to \t", NumSymb[id]]];,
		ns = (ns/.NumSymb->NumSymbN);
	];
	
	If[OptionValue[Verbose],Print["Simplifying \t", NumSymb[id],", now with LeafCount ", LeafCount[ns/.NumSymb->NumSymbSimplified]]];
	
	(*$debugNS = (ns/.NumSymb\[Rule]NumSymbN)/.NumSymb->NumSymbSimplified;*)
	
	(*NumSymbSimplified[id]=
	Simplify@Replace[
		Simplify[
			Together[
				NumSymb`Private`ResolveAbs@Refine[ns/.NumSymb\[Rule]NumSymbSimplified,Assumptions->(#>0&/@RationalSeries`$NumericVariables)]//.
					{Power[x_,-1/2]:>1/sqrt[x],Power[x_,1/2]:>sqrt[x],Sqrt->sqrt}
			]
		],
	{Plus[x__] /; !FreeQ[{x},_sqrt] \[RuleDelayed] Collect[Plus[x],Cases[Variables[{x}],_sqrt],Simplify]},Infinity]/.sqrt->Sqrt;*)
	
	(*NumSymbSimplified[id]= Simplify[Together[
		ResolveAbs@Refine[ns/.NumSymb->NumSymbSimplified,Assumptions->(#>0&/@RationalSeries`$NumericVariables)]//.
			{Power[x_,-1/2]:>1/sqrt[x],Power[x_,1/2]:>sqrt[x],Sqrt->sqrt}
	]]/.sqrt->Sqrt;*)
	
	NumSymbSimplified[id]= Simplify[
		ResolveAbs@Refine[ns/.NumSymb->NumSymbSimplified,Assumptions->(#>0&/@RationalSeries`$NumericVariables)]//.
			{Power[x_,-1/2]:>1/sqrt[x],Power[x_,1/2]:>sqrt[x],Sqrt->sqrt}
	]/.sqrt->Sqrt;
	
	NumSymbN[id]=NumSymbSimplified[id];
	
	];


ResolveAbs[expr_]:=expr/.Abs[f_]:>Which[RationalSeries`SampleEval[f]>0,f,RationalSeries`SampleEval[f]<0,-f,RationalSeries`SampleEval[f]==0,0];


Options[SimplifyNumSymbs]={Verbose->False}


SimplifyNumSymbs[expr_,opts:OptionsPattern[]]:=SimplifyNestedNumSymbs[#,PassRules[SimplifyNumSymbs,SimplifyNestedNumSymbs,opts]]&/@GetAllNumSymbs[expr];


Options[SimplifyNumSymbExpression]={Verbose->False}


SimplifyNumSymbExpression[expr_,opts:OptionsPattern[]]:=
Module[{intermediate},
	intermediate = CoalesceNumSymb@NumSymbReplace@expr;
	SimplifyNumSymbs[intermediate,PassRules[SimplifyNumSymbExpression,SimplifyNumSymbs,opts]];
	intermediate = intermediate/.NumSymb->NumSymbSimplified//Together//Simplify;
	Simplify[intermediate/.sqrt->Sqrt,Assumptions->(#>0&/@$NumericVariables)]
];


(* ::Chapter:: *)
(*Postamble*)


End[];
Protect["NumSymb`*"];
Unprotect[NumSymbLookUp,NumSymb,NumSymbN,NumSymbEvaluated,NumSymbSimplified,$LeafCountLimit];

EndPackage[];
