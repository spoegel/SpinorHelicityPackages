(* ::Package:: *)

(* ::Title:: *)
(*Generate Finite Field Momenta*)


(* ::Chapter:: *)
(*Preamble*)


(* :Title: Amplitudes.m -- a package collecting various amplitudes *)

(* :Context: SpinorHelicity6D` *)

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
   SpinorHelicity6D.wl handling 6D Spinors, written by M. Huber and S. de Angelis
	SpinorFunctions.wl for definition of trace objects
*)

(* :Todo:
*)


BeginPackage["SpinorHelicity6D`MomentumTwistors`",
	"SpinorHelicity6D`SpinorFunctions`",
	"SpinorHelicity6D`",
	"SpinorHelicity6D`",
	"MomentumTwistorDefinitions`",
	"Utils`"];

Unprotect["SpinorHelicity6D`MomentumTwistors`*"];
ClearAll["SpinorHelicity6D`MomentumTwistors`*"];
ClearAll["SpinorHelicity6D`MomentumTwistors`Private`*"];


GenerateSpinorsFromTwistorMatrix::usage = "Sets up kinematics based on passed twistor matrix.";

DefineKinematicsFromSpinors::usage = "Sets up kinematic defintions given a set of 4D spinors.";

ModConditional::usage = "Evaluates normally when passing Modulus 0, otherwise evaluating on the finite field.";
RationalMod::usage = "Evaluates Power and Rational on the finite field.";
Mom4DVectorNFromMatrixMod::usage = "Finite field version of Mom4DVectorNFromMatrix of SpinorHelicity6D`";

MomentumTwistorParameterSolution::usage = "";

Begin["`Private`"];


(* ::Chapter:: *)
(*Main*)


(* ::Section::Closed:: *)
(*Finite Field Functions*)


ModConditional[modulus_Integer] := Function[expr,
	If[modulus === 0,
		expr
	,
		Mod[RationalMod[expr,modulus],modulus]
	]
];


RationalMod[expr_,modulus_Integer] := expr//.{Power[x_,y_]:>PowerMod[x,y,modulus],Rational[x_,y_]:>Mod[x*PowerMod[y,-1,modulus],modulus]};


(* note the change to +-+- signature *)
Mom4DVectorNFromMatrixMod[label_,prime_Integer][type_]/; prime =!= 0 :=Module[{comps=Evaluate[Mom4DN[label][type][$down]]},
   Mom4DVectorN[label][type][$up] = 1/2*{comps[[1,1]]+comps[[2,2]],-comps[[1,2]]-comps[[2,1]],-comps[[2,1]]+comps[[1,2]],comps[[2,2]]-comps[[1,1]]} // ModConditional[prime];
   Mom4DVectorN[label][type][$down] = {+1,-1,+1,-1}*Mom4DVectorN[label][type][$up] // ModConditional[prime];
];





(* ::Section::Closed:: *)
(*Setting up Kinematics from Lambda & LambdaTilde*)


Options[DefineKinematicsFromSpinors]={Modulus->0};

DefineKinematicsFromSpinors[label_,lambda_,lambdaTilde_,opts:OptionsPattern[]] := Module[{},
	SpinorUndotN[label][$lam][$down] = lambda;
	SpinorUndotN[label][$lam][$up] = Normal[LeviCivitaTensor[2]].lambda // ModConditional[OptionValue[Modulus]];
	SpinorUndotN[label][$mu][$down] = {0,0};
	SpinorUndotN[label][$mu][$up] = {0,0};
	
	SpinorDotN[label][$lam][$down] = lambdaTilde;
	SpinorDotN[label][$lam][$up] = Normal[LeviCivitaTensor[2]].lambdaTilde // ModConditional[OptionValue[Modulus]];
	SpinorDotN[label][$mu][$down] = {0,0};
	SpinorDotN[label][$mu][$up] = {0,0};
	
	(*Mom4DN[label][Null][$up] = KroneckerProduct[SpinorUndotN[label][$lam][$up],SpinorDotN[label][$lam][$up]] // ModConditional[OptionValue[Modulus]];
	Mom4DN[label][Null][$down] = KroneckerProduct[SpinorDotN[label][$lam][$down],SpinorUndotN[label][$lam][$down]] // ModConditional[OptionValue[Modulus]];
	Mom4DN[label][$flat][$up] = Mom4DN[label][Null][$up];
	Mom4DN[label][$flat][$down] = Mom4DN[label][Null][$down];*)
	
	Activate @ {
		Inactive[SetDelayed][SpinorUndotN[label][$lam][j_/;MemberQ[{1,2},j]][Null],Indexed[SpinorUndotN[label][$lam][$up],j]],
		Inactive[SetDelayed][SpinorUndotN[label][$lam][Null][j_/;MemberQ[{1,2},j]],Indexed[SpinorUndotN[label][$lam][$down],j]],
		Inactive[SetDelayed][SpinorDotN[label][$lam][Null][j_/;MemberQ[{1,2},j]],Indexed[SpinorDotN[label][$lam][$down],j]],
		Inactive[SetDelayed][SpinorDotN[label][$lam][j_/;MemberQ[{1,2},j]][Null],Indexed[SpinorDotN[label][$lam][$up],j]]
	};
	
	ExtramassN[label] = 0;
	ExtramasstildeN[label] = 0;
	
	If[IntegerQ[OptionValue[Modulus]] && OptionValue[Modulus]=!=0,	
		Mom4DVectorNFromMatrixMod[label,OptionValue[Modulus]][Null];
		Mom4DVectorNFromMatrixMod[label,OptionValue[Modulus]][$flat];
	,
		Mom4DVectorNFromMatrix[label][Null];
		Mom4DVectorNFromMatrix[label][$flat];
	];
];


(* ::Section::Closed:: *)
(*Main Twistor to Kinematics Function*)


Options[GenerateSpinorsFromTwistorMatrix]={Modulus->0};

GenerateSpinorsFromTwistorMatrix[labels_List, twistorMatrix_List, opts : OptionsPattern[]] /; Dimensions[twistorMatrix]=={4,Length@labels} := 
Module[{lambda,lambdaTilde},
(*	lambda = (Transpose@twistorMatrix)[[All, 1 ;; 2]];*)

(*	lambdaTilde = LambdaTildeFromTwistor[*)
(*		twistorMatrix,*)
(*		FilterRules[*)
(*			{*)
(*				opts,*)
(*				Options[GenerateSpinorsFromTwistorMatrix]*)
(*			},*)
(*			Options[LambdaTildeFromTwistor]*)
(*		]*)
(*	];*)

	{lambda,lambdaTilde} = LambdaLambdaTildeFromTwistor[
		twistorMatrix,
		PassRules[GenerateSpinorsFromTwistorMatrix,LambdaLambdaTildeFromTwistor,opts]
	];

	MapThread[
		DefineKinematicsFromSpinors[
			#1,
			#2,
			#3,
			FilterRules[
				{
					opts,
					Options[GenerateSpinorsFromTwistorMatrix]
				},
				Options[DefineKinematicsFromSpinors]]
		]&,
		{labels,lambda,lambdaTilde}
	];
	(*Print["Test: p_1 = ",momenta[[1]],". p_1^2 = ", Mod[Plus@@(momenta[[1]]*momenta[[1]]*{1,-1,1,-1}),prime]];
	Print["Test momentum conservation: Sum[p_i] = ", Mod[Plus@@momenta,prime]];*)
];



MomentumTwistorParameterSolution[n_]:=MomentumTwistorParameterSolutionTemplate[n]/.{
	$momTwistorInv -> S,
	$momTwistorSpaaProd -> SpinorAngleBracket,
	$momTwistorSpbbProd -> SpinorSquareBracket,
	$momTwistorSpaaChain[a_,b___,c_]:> Chain[$angle,a,{b},c,$angle],
	$momTwistorSpabChain[a_,b___,c_]:> Chain[$angle,a,{b},c,$square],
	$momTwistorSpbaChain[a_,b___,c_]:> Chain[$square,a,{b},c,$angle],
	$momTwistorSpbbChain[a_,b___,c_]:> Chain[$square,a,{b},c,$square]
};


(* ::Chapter:: *)
(*Postamble*)


End[];
Protect["SpinorHelicity6D`MomentumTwistors`"];

EndPackage[];
