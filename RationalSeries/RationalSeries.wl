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
(*RationalSeries*)


(* ::Chapter:: *)
(*Preamble*)


(* :Title: RationalSeries.m -- a package to compute Series coefficients of rational functions efficiently *)

(* :Context: RationalSeries` *)

(* :Summary:

 *)

(* :Copyright: \[Copyright] 2020 by Sebastian Poegel *)

(* :Package Version: 0.1 *)

(* :Mathematica Version: 12.0 *)

(* :History:
   0.1 Copying functions from notebook
*)

(* :Keywords: *)

(* :Sources:
*)

(* :Warnings:

*)

(* :Limitations:

*)

(* :Discussion:

*)

(* :Requirements:
   Utils` for the definition of NumSymb,NumSymbN, NumSymbEvaluated. Theses may be merged here in the future.
   We also need PassRules from Utils`.
*)

(* :Todo:
*)


BeginPackage["RationalSeries`",{"NumSymb`","SpinorHelicity6D`","SpinorHelicity6D`Unitarity`",
	"SpinorHelicity6D`Amplitudes`","SpinorHelicity6D`TwoLoopAllPlusRationals`","Utils`"}];

Unprotect["RationalSeries`*"];
ClearAll["RationalSeries`*"];
ClearAll["RationalSeries`Private`*"];

SampleEval::usage = "Function responsible for consistently evaluating expressions numerically on a sample point.";
ZeroCoefficientQ::usage = "Test whether an expression is zero. It uses SampleEval to generate the numeric value.";

$NumericVariablesRules::usage = "Global variable to specify variables which will later take on numeric values";
$NumericVariablesRules = {};
$NumericVariables::usage = "";
$NumericVariables = {};

GetLowestOrder::usage = "";
LowestPowerSeries::usage = "";

RationalSeries::usage = "Computes the Series expansion of rational functions more or less efficiently.";
RationalSeriesCoefficient::usage = "Computes the Series expansion of rational functions more or less efficiently.";

$TotalTimingsRationalSeries::usage = "Log of the total runtime spent in the functions of RationalSeries.";
$TotalTimingsRationalSeries = <|"PossibleZeroQ"-><|"Total"->0,"Max"->0|>,"SampleEval"-><|"Total"->0,"Max"->0|>,
	"LowestPowerSeries"-><|"Total"->0,"Max"->0|>,"RationalSeries"-><|"Total"->0,"Max"->0|>|>;

$TotalCallsRationalSeries::usage = "Log of the total runtime spent in the functions.";
$TotalCallsRationalSeries = <|"CurrentCall"->Null,"PossibleZeroQ"->0,"SampleEval"->0,"RationalSeriesCoefficient"->0|>;

$ZeroCoefficientQData::usage  = "";
$SampleEvalData::usage  = "";
$ZeroCoefficientQData = <||>;
$SampleEvalData = <||>;

ResetTimingsCallsRationalSeries::usage = "";

$PossibleZeroQTimeOut = 10;

AnalyticVariables::usage = "";
NumericVariables::usage = "";

VarsSample::usage = "VarsSample[par] returns the numeric value that SampleEval is assuming for an analytic parameter par. VarsSample threads over lists.";


Begin["`Private`"];


(* ::Chapter:: *)
(*Main*)


(* ::Section:: *)
(*Add Timing*)


ResetTimingsCallsRationalSeries[] := (
	$TotalCallsRationalSeries = <|"CurrentCall"->Null,"PossibleZeroQ"->0,"SampleEval"->0,"RationalSeriesCoefficient"->0|>;
	$TotalTimingsRationalSeries = <|"PossibleZeroQ"-><|"Total"->0,"Max"->0|>,"SampleEval"-><|"Total"->0,"Max"->0|>,
		"LowestPowerSeries"-><|"Total"->0,"Max"->0|>,"RationalSeries"-><|"Total"->0,"Max"->0|>|>;
);


AddTimingRationalSeries[fctName_String]:=Function[timingList,
	$TotalTimingsRationalSeries[fctName,"Total"] += timingList[[1]];
	If[$TotalTimingsRationalSeries[fctName,"Max"] < timingList[[1]],
		$TotalTimingsRationalSeries[fctName,"Max"] = timingList[[1]];
	];
	
	timingList[[2]]
];


(* ::Section:: *)
(*Sample Evaluate*)

VarsSample[vars_List] := VarsSample/@vars;
VarsSample[var_] := (FromDigits[IntegerDigits[#][[;;7]]]&@Hash[var]);


SampleEval[expr_Plus,otherVars_ : {}]:=SampleEval[expr,otherVars]=SampleEval[#,otherVars]&/@expr;
SampleEval[expr_Times,otherVars_ : {}]:=SampleEval[expr,otherVars]=SampleEval[#,otherVars]&/@expr;
SampleEval[expr:Power[base_,exp_],otherVars_ : {}]:=SampleEval[expr,otherVars]=Power[SampleEval[base,otherVars],exp];


SampleEval[expr_,otherVars_ : {}] := SampleEval[expr,otherVars] = 
Module[
	{
		callTemp = $TotalCallsRationalSeries["CurrentCall"],
		result,
		exprPartialEval,
		vars,
		varsValues
	},
	$TotalCallsRationalSeries["CurrentCall"] = "SampleEval with ByteCount " <> ToString@ByteCount@expr;
	
	$TotalCallsRationalSeries["SampleEval"]++;

	vars = Join[otherVars,$NumericVariables];
	varsValues = VarsSample @ vars;

	result = AbsoluteTiming[
		expr//.Dispatch@Join[$NumericVariablesRules,{NumSymb`NumSymb[x_] :>
		NumSymb`NumSymbEvaluated[x]},{Global`tCoeff[_]:>1,Global`yCoeff[_]:>1},
		{AmplitudeCoefficient[_]:>1}, MapThread[Rule,{vars,varsValues}]
		]
	];

	
	If[MemberQ[$ContextPath, "SpinorHelicity6D`"]===True,
		result = ToNum[result];
	];
	(*If[MemberQ[Contexts[],"SpinorHelicity6D`TwoLoopAllPlusRationals`"],
		AppendTo[$SampleEvalData,Unique[]\[Rule]<|"Configuration"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[1]],
			"Loop1"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[2]],
			"Loop2"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[3]],
			"Timing"->result[[1]]|>];
	];*)
	
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	
	AddTimingRationalSeries["SampleEval"]@result
];


(* ::Section:: *)
(*Eliminate Zero Series Coefficients*)


PossibleZeroQWrapper[eval_]:=PossibleZeroQ[eval,Method->"ExactAlgebraics"];


ZeroCoefficientQ[expr_,otherVars_] := ZeroCoefficientQ[expr,otherVars] = 
Module[{eval,result,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "ZeroCoefficientQ";

	$TotalCallsRationalSeries["PossibleZeroQ"]++;
	
	$TotalCallsRationalSeries["CurrentCall"] = "ZeroCoefficientQ SampleEval";
	eval = SampleEval[expr,otherVars];
	If[!NumericQ[eval],
		Print["Warning: eval is not numeric, free vars: ", Variables@eval,
					"Configuration: ", SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration];
	];
	
	(*Global`testEval = eval;
	Global`testExpr = expr;
	Global`testVars = otherVars;
	DumpSave["ZeroCoefficientQEval.mx",{Global`testEval,Global`testExpr,Global`testVars}];*)
	
	$TotalCallsRationalSeries["CurrentCall"] =
     "ZeroCoefficientQ after eval, ByteCount: " <> ToString[ByteCount@eval] <>
         ", Hash: " <> ToString[Hash@eval];
	
	
	result=AbsoluteTiming[
			TimeConstrained[
				PossibleZeroQ[eval,Method->"ExactAlgebraics"],
				$PossibleZeroQTimeOut,
				PossibleZeroQ[RootReduce@eval,Method->"ExactAlgebraics"]
			]
			(*PossibleZeroQWrapper[eval]*)
		];
	
	(*If[MemberQ[Contexts[],"SpinorHelicity6D`TwoLoopAllPlusRationals`"],
		AppendTo[$ZeroCoefficientQData,Unique[]\[Rule]<|"Configuration"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[1]],
			"Loop1"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[2]],
			"Loop2"->SpinorHelicity6D`TwoLoopAllPlusRationals`$currentConfiguration[[3]],
			"Timing"->result[[1]]|>];
	];*)
	
	$TotalCallsRationalSeries["CurrentCall"] = "ZeroCoefficientQ after PossibleZeroQ, ByteCount: " <> ToString[ByteCount@eval];
	(*If[result[[1]] > 120, Put[result[[1]],eval,"LongZeroCoefficientQ.m"]];*)

	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	ClearAll[Global`testEval,Global`testExpr,Global`testVars];
	
	AddTimingRationalSeries["PossibleZeroQ"]@result
];


ZeroCoefficientQ[otherVars_List] := Function[expr,ZeroCoefficientQ[expr,otherVars]];
ZeroCoefficientQ[expr_] := ZeroCoefficientQ[expr,{}];


EliminateZeroSeriesCoefficients[series_SeriesData,otherVars_] /; series[[3]] =!= {} :=MapAt[EliminateZeroSeriesCoefficient[otherVars],series,{3,All}];
EliminateZeroSeriesCoefficients[expr_,_] := expr;


EliminateZeroSeriesCoefficient[expr_, otherVars_] /; ZeroCoefficientQ[expr,otherVars] := 0;
EliminateZeroSeriesCoefficient[expr_, otherVars_] := expr;

EliminateZeroSeriesCoefficient[otherVars_] := Function[expr,EliminateZeroSeriesCoefficient[expr,otherVars]];


(* ::Section:: *)
(*SeriesData Helpers*)


SeriesDataQ[expr_] := MatchQ[expr,_SeriesData];
SameSeriesVarQ[series__SeriesData] := SameQ@@((List@@@{series})[[All,{1,2,6}]]);
SameSeriesVarQ[___] := False;


addSeries[series_]/;MemberQ[series,0]:=addSeries[DeleteCases[series,0]];


addSeries[series_?(VectorQ[#,SeriesDataQ]&)] /; SameSeriesVarQ@@series := 
Module[{highestOrder = Min[series[[All,4]]],lowestOrder = Max[series[[All,5]]], order = Min[series[[All,5]]],
	coefficientLists,seriesVar,seriesPoint,seriesDenominator},
	
	{seriesVar,seriesPoint,seriesDenominator} = (List@@series[[1]])[[{1,2,6}]];
	coefficientLists = Take[
		Join[
			ConstantArray[0,#[[4]]-highestOrder],
			#[[3]],
			ConstantArray[0,(lowestOrder-#[[5]])+(#[[5]]-#[[4]]-Length@#[[3]])] 
		],
	order-highestOrder]&/@series;
SeriesData[seriesVar,seriesPoint,Plus@@coefficientLists,highestOrder,order,seriesDenominator]
];


powerSeries[expr_,expo_] := Power[expr,expo];

powerSeries[series_SeriesData,expo_Integer] := powerSeries[series,expo] =  
Module[{coefficientList,seriesVar,seriesPoint,seriesDenominator,orderHigh,orderLow,rules={},coeff},
	
	{seriesVar,seriesPoint,orderHigh,orderLow,seriesDenominator} = (List@@series)[[{1,2,4,5,6}]];
	
	coefficientList = Replace[series[[3]],{a_ /; a =!= 0 || a != 0 :> ((AppendTo[rules,coeff[#]->a]; coeff[#])&@Unique[])},{1}];

	Power[SeriesData[seriesVar,seriesPoint,coefficientList,orderHigh,orderLow,seriesDenominator],expo]/.rules
];


(* ::Section:: *)
(*LowestPowerSeries*)


Options[LowestPowerSeries]={Variables -> {}};


LowestPowerSeries[expr_,var_,extraOrders_, opts:OptionsPattern[]] /; FreeQ[expr,var] := LowestPowerSeries[expr,var,extraOrders, opts] = 
	If[ZeroCoefficientQ[expr,OptionValue[Variables]],
		0
	,
		(*SeriesData[var,Infinity,{NumSymb[genUnique[SampleEval[expr,OptionValue[Variables]]]]},0,1+extraOrders,1]*)
		SeriesData[var,Infinity,{SampleEval[expr,OptionValue[Variables]]},0,1+extraOrders,1]
	];
	
LowestPowerSeries[expr:var_,var_,extraOrders_, opts:OptionsPattern[]] := LowestPowerSeries[var,var,extraOrders, opts] = 
	SeriesData[var,Infinity,{1},-1,0+extraOrders,1];
	
LowestPowerSeries[expr:Power[var_,exp_],var_,extraOrders_, opts:OptionsPattern[]] := LowestPowerSeries[expr,var,extraOrders, opts] = 
	SeriesData[var,Infinity,{1},-exp,-exp+1+extraOrders,1];



LowestPowerSeries[expr_Times,var_,extraOrders_, opts:OptionsPattern[]] := LowestPowerSeries[expr,var,extraOrders, opts] = 
Module[{minSeries,res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Times: Series";
	
	minSeries = EliminateZeroSeriesCoefficients[
		Times@@(LowestPowerSeries[#,var,extraOrders,PassRules[LowestPowerSeries,LowestPowerSeries,opts]]&/@(List@@expr))
	,OptionValue[Variables]];

	$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Times: Check";
	res = CheckAndOutputLowestPowerSeries[expr,minSeries,var,extraOrders,PassRules[LowestPowerSeries,CheckAndOutputLowestPowerSeries,opts]];
	$TotalCallsRationalSeries["CurrentCall"]=callTemp;
	res
];


LowestPowerSeries[expr_Plus,var_,extraOrders_, opts:OptionsPattern[]] := LowestPowerSeries[expr,var,extraOrders, opts] = 
Module[{minSeries,res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Plus: Series";
	res = If[ZeroCoefficientQ[expr,Prepend[OptionValue[Variables],var]],
		0
	,
		$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Plus: Series";
		minSeries = 
			EliminateZeroSeriesCoefficients[
				addSeries[
					Replace[
						LowestPowerSeries[#,var,extraOrders,PassRules[LowestPowerSeries,LowestPowerSeries,opts]]&/@(List@@expr),
					0->Nothing,{1}]
				]
			,OptionValue[Variables]];
		$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Plus: Check";
		CheckAndOutputLowestPowerSeries[expr,minSeries,var,extraOrders,PassRules[LowestPowerSeries,CheckAndOutputLowestPowerSeries,opts]]
	];
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	res
];


LowestPowerSeries[expr:Power[arg_,exp_],var_,extraOrders_, opts:OptionsPattern[]] := LowestPowerSeries[expr,var,extraOrders, opts] = 
Module[{minExp,minSeries,res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Power: Series";
	minSeries = EliminateZeroSeriesCoefficients[
		powerSeries[#,exp]&@LowestPowerSeries[arg,var,extraOrders,PassRules[LowestPowerSeries,LowestPowerSeries,opts]]
	,OptionValue[Variables]];
	
	$TotalCallsRationalSeries["CurrentCall"] = "LowestPowerSeries Power: Check";
	res = CheckAndOutputLowestPowerSeries[expr,minSeries,var,extraOrders,PassRules[LowestPowerSeries,CheckAndOutputLowestPowerSeries,opts]];
	
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	res
]


Options[CheckAndOutputLowestPowerSeries] = Options[LowestPowerSeries];

CheckAndOutputLowestPowerSeries[expr_,minSeries_,var_,extraOrders_,opts:OptionsPattern[]]:=
If[(minSeries === 0) || (minSeries[[4]] < minSeries[[5]]),
	minSeries
,
	(*If[extraOrders>50, Print["Aborting, since we hit extraOrders = 50, so we are probably in a infinite loop"]; Abort[]];*)
	LowestPowerSeries[expr,var,extraOrders+1,PassRules[CheckAndOutputLowestPowerSeries,LowestPowerSeries,opts]]
];


(* ::Section:: *)
(*GetLowestOrder*)


Options[GetLowestOrder] = {Variables -> {}};
GetLowestOrder[expr_,var_,opts:OptionsPattern[]]:= GetLowestOrder[expr,var,opts] = 
Module[{res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	 $TotalCallsRationalSeries["CurrentCall"] = "GetLowestOrder";
	res = AddTimingRationalSeries["LowestPowerSeries"]@Timing[
		(If[#===0,Infinity,#[[4]]])&@LowestPowerSeries[expr,var,0,PassRules[GetLowestOrder,LowestPowerSeries,opts]]
	];
	
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	
	res
];


(* ::Section:: *)
(*RationalSeries*)


Options[RationalSeries]={AnalyticVariables->{},NumericVariables->{}};


RationalSeries[HoldPattern[SeriesData[x_,y_,coeffs_List,z_,w_,q_]],var_,order_, opts:OptionsPattern[]] :=
	SeriesData[x,y,
		RationalSeries[#,var,order,opts]&/@coeffs,
	z,w,q];


RationalSeries[expr_,var_,order_, opts:OptionsPattern[]] /; FreeQ[expr,var] := RationalSeries[expr,var,order, opts] = 
	If[ZeroCoefficientQ[expr,Join[{var},OptionValue[NumericVariables],OptionValue[AnalyticVariables]]],
		0
	,
		If[0 <= order,
			SeriesData[var,Infinity,{expr},0,order+1,1]
		,
			(* Mimicking behaviour of Series *)
			SeriesData[var,Infinity,{},0,0,1]
		]
	];
	
RationalSeries[expr:var_,var_,order_, opts:OptionsPattern[]] /; -1 <= order :=
	SeriesData[var,Infinity,{1},-1,order+1,1];
	
	(* Mimicking behaviour of Series *)
RationalSeries[expr:var_,var_,order_, opts:OptionsPattern[]] /; -1 > order :=
	SeriesData[var,Infinity,{},-1,-1,1];
	
RationalSeries[expr:Power[var_,exp_],var_,order_, opts:OptionsPattern[]] /; -exp <= order :=
	SeriesData[var,Infinity,{1},-exp,Max[-exp+1,order+1],1];

(* Mimicking behaviour of Series *)
RationalSeries[expr:Power[var_,exp_],var_,order_, opts:OptionsPattern[]] /; -exp > order :=
	SeriesData[var,Infinity,{},-exp,-exp,1];


RationalSeries[expr_Plus,var_,order_, opts:OptionsPattern[]] := RationalSeries[expr, var, order, opts] =
Module[{lowestPowers, deleteSummands, variables},
	variables = Join[OptionValue[NumericVariables], OptionValue[AnalyticVariables]];
	lowestPowers = GetLowestOrder[#, var, Variables -> variables] & /@ (List @@ expr);
	
	If[order < Min @@ lowestPowers || ZeroCoefficientQ[expr, Prepend[variables, var]],
		0
	,
		(* Removing summands whose order is lower than the requested one *)
		deleteSummands = Rule[#,0] & /@ Flatten@Position[lowestPowers, order < # &];

		EliminateZeroSeriesCoefficients[
			addSeries @ (RationalSeries[#, var, order, PassRules[RationalSeries,RationalSeries,opts]] & /@ (ReplacePart[List@@expr,deleteSummands])),variables
		] // EnsureLowestOrder[expr, var, Variables -> variables]
	]
];


RationalSeries[expr_Times,var_,order_, opts:OptionsPattern[]] := RationalSeries[expr,var,order, opts] = 
Module[{lowestPowers,powerSums},
	lowestPowers = GetLowestOrder[#,var,Variables->Join[OptionValue[NumericVariables],OptionValue[AnalyticVariables]]]&/@(List@@expr);
	(*Print["Times lowPow: ",lowestPowers];*)
	If[MemberQ[lowestPowers,Infinity],
		0
	,
		(* there is a minus in powerSums, since we need to match e.g. a term ~t^2 with something of ~1/t^2, so we have (-2)\[Rule](2) *)
		powerSums = Table[-Plus@@Delete[#,i],{i,1,Length@#}] & @ lowestPowers;
		(*Print["Times powSums: ",powerSums];*)
		EliminateZeroSeriesCoefficients[
			Times@@MapThread[RationalSeries[#1,var,order+#2,PassRules[RationalSeries,RationalSeries,opts]]&,{List@@expr,powerSums}]
		,Join[OptionValue[NumericVariables],OptionValue[AnalyticVariables]]]
	]
]


RationalSeries[expr:Power[arg_,exp_],var_,order_, opts:OptionsPattern[]] := RationalSeries[expr,var,order, opts] = 
Module[{lowestPower,res},
	lowestPower = GetLowestOrder[arg,var,Variables->Join[OptionValue[NumericVariables],OptionValue[AnalyticVariables]]];
	(*Print["Power: ",lowestPower];*)
	res=Which[lowestPower === Infinity && exp < 0,
		Print["Aborting: Zero argument raised to a negative exponent"];
		Indeterminate
	,order < lowestPower*exp,
		0
	,True,
		EliminateZeroSeriesCoefficients[
			powerSeries[RationalSeries[arg,var,order+lowestPower*(-exp+1),PassRules[RationalSeries,RationalSeries,opts]],exp]
		,Join[OptionValue[NumericVariables],OptionValue[AnalyticVariables]]]
	];
	res
];

(* This function is applied in RationalSeries after checking for zero coefficients. If the series is zero up to the
	 requested order, we still want to stay consistent with conventions, i.e. that in this case elements 4 and 5 give the
	 next non-zero order. This is exactly the result of GetLowestOrder, which
	 should have already been calculated and cached previously. *)
Options[EnsureLowestOrder] = Options[GetLowestOrder];

EnsureLowestOrder[expr_, var:Except[_Rule|_List], opts:OptionsPattern[]] := Function[series,EnsureLowestOrder[series, expr, var, opts]];

EnsureLowestOrder[series_SeriesData, expr_, var:Except[_Rule|_List], opts:OptionsPattern[]] /; series[[3]] === {} :=
	Module[{lowestOrder},
		lowestOrder = GetLowestOrder[expr,var,opts];
		ReplacePart[series, {4 -> lowestOrder, 5 -> lowestOrder}]
	];

EnsureLowestOrder[series_SeriesData, expr_, var:Except[_Rule|_List], opts:OptionsPattern[]] := series;


(* ::Section:: *)
(*RationalSeriesCoefficient*)


OrderInSeriesDataQ[series_SeriesData,order_Integer]:=And[!SameQ@@(Extract[series,{{4},{5}}]),IntervalMemberQ[Interval[Extract[series,{{4},{5}}]-{0,1}],order]];
ZeroInSeriesQ[series_SeriesData,order_Integer] := Extract[series,4] > order;
(*OrderInSeriesDataQ[series_,orders_List]:=And[!SameQ@@(Extract[series,{{4},{5}}]),
		Sequence@@(IntervalMemberQ[Interval[Extract[series,{{4},{5}}]-{0,1}],#]&/@orders)];*)


GetRationalSeriesCoefficient[series_,var_,powers_List] := GetRationalSeriesCoefficient[series,var,#]&/@powers;

GetRationalSeriesCoefficient[0,var_,power_Integer] := 0;
GetRationalSeriesCoefficient[series_SeriesData,var_,power_Integer] /; ZeroInSeriesQ[series,power] := 0;
GetRationalSeriesCoefficient[series_SeriesData,var_,power_Integer] /; OrderInSeriesDataQ[series,power]:=
Module[{hi,lo},
	{lo,hi} = Extract[series,{{4},{5}}];
	PadRight[Extract[series,{3}],hi-lo][[power-lo+1]]
];

GetRationalSeriesCoefficient[series_,var_,power_Integer] /; And[!MatchQ[series,_SeriesData],FreeQ[series,var]] := If[power === 0, series, 0];


Options[RationalSeriesCoefficient] = Options[RationalSeries];


RationalSeriesCoefficient[expr_List,args___] := RationalSeriesCoefficient[#,args]&/@expr;
RationalSeriesCoefficient[expr_Plus,args___] := RationalSeriesCoefficient[#,args]&/@expr;

RationalSeriesCoefficient[expr_,var_,power_, opts:OptionsPattern[]] := 
Module[{res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "RationalSeries";
	res = GetRationalSeriesCoefficient[AddTimingRationalSeries["RationalSeries"]@AbsoluteTiming@
			RationalSeries[expr,var,power, PassRules[RationalSeriesCoefficient,RationalSeries,opts]],
		var,power];
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	$TotalCallsRationalSeries["RationalSeriesCoefficient"]++;
	res
];
RationalSeriesCoefficient[expr_,var_,powers_List, opts:OptionsPattern[]] :=
Module[{res,callTemp = $TotalCallsRationalSeries["CurrentCall"]},
	$TotalCallsRationalSeries["CurrentCall"] = "RationalSeries";
	res = GetRationalSeriesCoefficient[AddTimingRationalSeries["RationalSeries"]@AbsoluteTiming@
			RationalSeries[expr,var,Max@powers, PassRules[RationalSeriesCoefficient,RationalSeries,opts]],
		var,powers];
	$TotalCallsRationalSeries["CurrentCall"] = callTemp;
	$TotalCallsRationalSeries["RationalSeriesCoefficient"]++;
	res
];


(* ::Chapter:: *)
(*Postamble*)


End[];
Protect["RationalSeries`*"];
Unprotect[$TotalCallsRationalSeries,$TotalTimingsRationalSeries,$NumericVariablesRules,$NumericVariables,SampleEval,RationalSeries,ZeroCoefficientQ,
		LowestPowerSeries,GetLowestOrder,$PossibleZeroQTimeOut,$SampleEvalData,$ZeroCoefficientQData];

EndPackage[];
