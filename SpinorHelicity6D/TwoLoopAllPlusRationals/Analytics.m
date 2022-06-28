(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: Analytics *)
(* :Context: Analytics` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2022-03-02 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2022 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)


BeginPackage["SpinorHelicity6D`TwoLoopAllPlusRationals`Analytics`",
	     {
		     "SpinorHelicity6D`",
		     "SpinorHelicity6D`TwoLoopAllPlusRationals`",
		     "RationalSeries`",
		     "Utils`",
		     "NumSymb`",
		     "MomentumTwistorDefinitions`"
	     }
];

GenerateCutEvaluationFunction::usage = "GenerateCutEvaluationFunction[expr,parameters] returns a function that efficiently evaluates the expression on a kinematic point.";

Begin["`Private`"];

GenerateCutEvaluationFunction[expr_,parameters_List]:=
	Module[
		{numSymbs,numSymbsEval,pars,evalFunction,numSymbsEvalFunction,dummy},
		numSymbs = Join[{dummy},{expr//NumSymb`GetAllNumSymbs//Flatten//Reverse//DeleteDuplicates}]//Flatten;

		numSymbsEval = numSymbs /. NumSymb`NumSymb->NumSymb`NumSymbN /. MapThread[Rule,{parameters,Table[Indexed[pars,i],{i, Length@parameters}]}];
			
		numSymbsEvalFunction = GenerateOptimizedEvaluationFunction[numSymbsEval,parameters];
		evalFunction = GenerateOptimizedEvaluationFunction[expr,parameters];

		Function@@Hold[{pars},Block[{NumSymb},#1=RootReduce@#2[pars]; RootReduce@#3[pars]]]&[numSymbs,numSymbsEvalFunction,evalFunction]
	];

GenerateCutEvaluationFunction[expr_,n_Integer] := GenerateCutEvaluationFunction[expr, MomentumTwistorParameters[n]];

End[]; (* `Private` *)

EndPackage[]
