(* ::Package:: *)

(* ::Title:: *)
(*Analytic Tools*)


(* ::Chapter:: *)
(*Preamble*)



BeginPackage["SpinorHelicity6D`AnalyticTools`", "SpinorHelicity6D`","Utils`"];

Unprotect["SpinorHelicity6D`AnalyticTools`*"];
ClearAll["SpinorHelicity6D`AnalyticTools`*"];
ClearAll["SpinorHelicity6D`AnalyticTools`Private`*"];


(* ::Section:: *)
(*Exports*)


MultiCollect::usage = "";


ApplyAnalyticTransformationRules::usage = "Apply a series of sets of transformation rules on expression.
	Functions can be specified to be applied to the expression before a set, after a set, and after every application of a
	set of rules. The rules and functions have to be provided as a list of Associations in the format:
	<|\"Pre\"->PreFunction,\"Inter\"->InterFunction,\"Post\"->PostFunction,\"Rules\"->Rules|>";

(* ::Chapter:: *)
(*Main*)


Begin["`Private`"];    (* begin the private context (implementation part) *)


(* ::Section:: *)
(*MultiCollect*)


Options[MultiCollect]={};
MultiCollect[expr_,var:Except[_List,_],h:_:Identity,opts:OptionsPattern[]]:=Collect[expr,var,h];
MultiCollect[expr_,{var:Except[_List,_]},h:_:Identity,opts:OptionsPattern[]]:=MultiCollect[expr,var,h,opts];
MultiCollect[expr_,{var:Except[_List,_],varsRest__},h:_:Identity,opts:OptionsPattern[]]:=
	MultiCollect[#,{varsRest},h,opts]&/@MultiCollect[expr,var,h,opts];



(* ::Section:: *)
(*ApplyAnalyticTransformationRules*)


Options[ApplyAnalyticTransformationRules]={DefaultRules->{},MaxRecursion->10};


ApplyAnalyticTransformationRules[expr_,ruleAssociations_?(VectorQ[#,AssociationQ]&),opts:OptionsPattern[]]:=
		Fold[
			applyAnalyticTransformationRule[#1,#2,PassRules[ApplyAnalyticTransformationRules,applyAnalyticTransformationRule,opts]]&,
			expr,
			ruleAssociations
		];


(* ::Subsection:: *)
(*applyAnalyticTransformationRule*)


Options[applyAnalyticTransformationRule]=Options[ApplyAnalyticTransformationRules];


applyAnalyticTransformationRule[expr_,rulesAssoc_Association,opts:OptionsPattern[]]:=
    Module[{pre,inter,post,rules},
			{pre,inter,post} = Lookup[rulesAssoc,#,Identity]&/@{"Pre","Inter","Post"};
			rules = Lookup[rulesAssoc,"Rules",{}];

			post@applyRules[pre@expr,rules,inter,PassRules[applyAnalyticTransformationRule,applyRules,opts]]
		];


Options[applyRules]=Options[applyAnalyticTransformationRule];


applyRules[expr_,{},interFunction_,opts:OptionsPattern[]]:=expr//.OptionValue[DefaultRules];


applyRules[expr_,{rules__List},interFunction_,opts:OptionsPattern[]]:=
		Fold[
			applyRules[#1,#2,interFunction,PassRules[applyRules,applyRules,opts]]&,
			expr,
			{rules}
		];


applyRules[expr_,rules:{(_Rule|_RuleDelayed)..},interFunction_,opts:OptionsPattern[]]:=
		NestWhile[
			interFunction@*ReplaceRepeated[OptionValue[DefaultRules]]@*ReplaceRepeated[rules]@*ReplaceRepeated[OptionValue[DefaultRules]],
			expr,
			UnsameQ,
			2,
			OptionValue[MaxRecursion]
		];



(* ::Chapter:: *)
(*Postamble*)


End[];


Protect[ "SpinorHelicity6D`AnalyticTools`*" ];  
Unprotect[{}];

EndPackage[ ]; 
