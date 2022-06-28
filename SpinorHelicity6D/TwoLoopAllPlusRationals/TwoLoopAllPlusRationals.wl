(* ::Package:: *)

(* ::Title:: *)
(*Two-Loop All-Plus Rational Terms Config *)


(* :Title: TwoLoopAllPlusRationals.m -- a package for computing the two-loop rational terms config-by-config *)

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


(* ::Chapter:: *)
(*Preamble*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`TwoLoopAllPlusRationals`",{
	"SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`",
	"SpinorHelicity6D`TwoLoopAllPlusRationals`Analytics`",
	"SpinorHelicity6D`TwoLoopAllPlusRationals`DrawConfigs`",
	"SpinorHelicity6D`TwoLoopAllPlusRationals`Configurations`",
	"SpinorHelicity6D`Unitarity`",
	"SpinorHelicity6D`SpinorFunctions`",
	"SpinorHelicity6D`Amplitudes`",
	"SpinorHelicity6D`",
	"RationalSeries`",
	"Utils`",
	"NumSymb`",
	"SpinorHelicity6D`BerendsGiele`",
	"SpinorHelicity6D`Unitarity`LoopMomentum`",
	"SpinorHelicity6D`Unitarity`OneLoopIntegrals`",
	"SpinorHelicity6D`Unitarity`ConfigurationTools`",
	"SpinorHelicity6D`Unitarity`RationalTerms`"
}];

Unprotect["SpinorHelicity6D`TwoLoopAllPlusRationals`*"];
ClearAll["SpinorHelicity6D`TwoLoopAllPlusRationals`*"];
ClearAll["SpinorHelicity6D`TwoLoopAllPlusRationals`Private`*"];

(* usage messages for the exported functions and the context itself *)
ComputeConfiguration::usage = ""; 

Recompute::usage = "";
LogDirectory::usage = "";
DataDirectory::usage = "";
ReuseResultQ::usage = "";

$rational::usage = "";
$rationalFull::usage = "";
$rationalFinal::usage = "";

$currentConfiguration::usage = "";
$currentConfiguration = Null;

CustomConfigs::usage = "";


$debugFirstLoopMom::usage = "";
$debugSecondLoopMom::usage = "";
$debugFirstTrees::usage = "";
$debugCenterTrees::usage = "";
$debugSecondTrees::usage = "";


$missingConfig::usage = "";

ShuffleConfigs::usage = "";
PostProcessing::usage = "";
NumSymbReplacement::usage = "";
ResultFormat::usage = "";
DebugLevel::usage = "";

$CutResult::usage = "Symbol used when exporting result in .mx format.";
WriteResult::usage = "Option for ComputeConfigurations. If set to \"MX\" or \"M\",
the result is exported as a .mx or .m file. For all other values, nothing is written to disk.";
ApplyIntegrals::usage = "Option for ComputeConfigurations. If set to False, only the integral
coefficient matrices are returned.";

ComputeConfigurations::usage = "";
ComputeConfiguration::missingAmp = "Amplitudes missing for configuration `1`";
ComputeConfiguration::EvalNoPars = "Warning: Output as EvaluationFunction request, but no analytic parameters passed. Maybe you forgot to set ";
EvaluationFunction::usage = "Postprocessing option for ComputeConfiguration and ComputeConfigurations to generate an numerical evaluation function as output.";

TwoLoopAllPlusRationalsTermConfigs::usage = "Deprecated, used ComputeConfigurations instead.";

$TwoLoopAllPlusRationalsStatus::usage = "Shows the current evaluation status.";
$TwoLoopAllPlusRationalsStatus = Null;

Begin["`Private`"];    (* begin the private context (implementation part) *)


(* ::Chapter:: *)
(*Parallel Configurations Main*)


(* ::Section:: *)
(*Loop Momenta*)


Options[FixLoopMomenta] = {
	ComputeLoopMomentumSpinors -> False,
	BubbleReference -> $bubbleReference,
	LoopMomentumPostProcessingFirstPass -> Identity,
	Gravity->False,
	ReduceSqrts->True
};


FixLoopMomenta[totalConfig_,opts:OptionsPattern[]] :=
	{FixFirstLoopMomentum[totalConfig,PassRules[FixLoopMomenta,FixFirstLoopMomentum,opts]],
	 FixSecondLoopMomentum[totalConfig,PassRules[FixLoopMomenta,FixSecondLoopMomentum,opts]]};


(* ::Subsection::Closed:: *)
(*Second Loop Momenta*)


Options[FixSecondLoopMomentum] = Options[FixLoopMomenta];


(* ::Subsubsection::Closed:: *)
(*Box*)


FixSecondLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,secondLoopBoxPattern] := Module[
{assocOut=<||>,l21=Unique["l21$"],l22=Unique["l22$"],l23=Unique["l23$"],l24=Unique["l24$"],
	Mu=Unique["Mu2$"],MuTilde=Unique["MuTilde2$"], config},
	
	config = ConstructSecondLoopConfig[totalConfig];
	
	FixBoxLoopMomentum6D[config,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l21,l22,l23,l24},LargeMuLimit->True,
			PassRules[FixSecondLoopMomentum,FixBoxLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Box";
	assocOut["Momenta"] = {l21,l22,l23,l24};
	assocOut["Variations"] = {$plus,$minus};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde|>;
		
	assocOut
];


(* ::Subsubsection::Closed:: *)
(*Triangle*)


FixSecondLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,secondLoopTrianglePattern] := Module[
{assocOut=<||>,l21=Unique["l21$"],l22=Unique["l22$"],l23=Unique["l23$"],
	Mu=Unique["Mu2$"],MuTilde=Unique["MuTilde2$"],t=Unique["t2$"], config},
	
	config = ConstructSecondLoopConfig[totalConfig];

	FixTriangleLoopMomentum6D[config,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l21,l22,l23},
			PassRules[FixSecondLoopMomentum,FixTriangleLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Triangle";
	assocOut["Momenta"] = {l21,l22,l23};
	assocOut["Variations"] = {$star,$unstar};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;
	
	assocOut
];


(* ::Subsubsection::Closed:: *)
(*Bubble*)


FixSecondLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,secondLoopBubblePattern] := Module[
{assocOut=<||>,l21=Unique["l21$"],l22=Unique["l22$"],
	Mu=Unique["Mu2$"],MuTilde=Unique["MuTilde2$"],t=Unique["t2$"],y=Unique["y2$"], config, yIntegrals},
	
	config = ConstructSecondLoopConfig[totalConfig];
	
	yIntegrals = FixBubbleLoopMomentum6D[config,OptionValue[BubbleReference],yVar->y,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l21,l22},
			PassRules[FixSecondLoopMomentum,FixBubbleLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Bubble";
	assocOut["Momenta"] = {l21,l22};
	assocOut["Variations"] = {Null};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t,"y"->y|>;
	assocOut["yIntegrals"] = ToNum@yIntegrals;
	
	assocOut
];


(* ::Subsubsection:: *)
(*BubbleTriangle*)


FixSecondLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,secondLoopBubbleTrianglePattern] := Module[
{assocOut=<||>,l21=Unique["l21$"],l22=Unique["l22$"],l23=Unique["l23$"],
	Mu=Unique["Mu2$"],MuTilde=Unique["MuTilde2$"],t=Unique["t2$"], baseConfig, subConfig, splittingMomenta, tIntegrals},
	
	baseConfig = ConstructSecondLoopBaseConfig[totalConfig];
	subConfig = ConstructSecondLoopConfig[totalConfig];
	
	splittingMomenta = GetSplittingMomenta[baseConfig,subConfig];
	
	tIntegrals = FixBubbleTriSubLoopMomentum6D[baseConfig,splittingMomenta,OptionValue[BubbleReference],tVar->t,
			MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l21,l22,l23},LargeTLimit->True,
			PassRules[FixSecondLoopMomentum,FixBubbleTriSubLoopMomentum6D,opts]];
	
	assocOut["Type"] = "BubbleTriangle";
	assocOut["Momenta"] = {l21,l22,l23};
	assocOut["Variations"] = {$plus,$minus};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;
	assocOut["tIntegrals"] = ToNum@tIntegrals;
	
	assocOut
];


(* ::Subsection:: *)
(*First Loop Momenta*)


Options[FixFirstLoopMomentum] = Options[FixLoopMomenta];


(* ::Subsubsection::Closed:: *)
(*Box*)


FixFirstLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,firstLoopBoxPattern] := Module[
{assocOut=<||>,l11=Unique["l11$"],l12=Unique["l12$"],l13=Unique["l13$"],l14=Unique["l14$"],
	Mu=Unique["Mu1$"],MuTilde=Unique["MuTilde1$"], config},
	
	config = ConstructFirstLoopConfig[totalConfig];
	
	FixBoxLoopMomentum6D[config,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l11,l12,l13,l14},LargeMuLimit->True,
			PassRules[FixFirstLoopMomentum,FixBoxLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Box";
	assocOut["Momenta"] = {l11,l12,l13,l14};
	assocOut["Variations"] = {$plus,$minus};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde|>;
		
	assocOut
];


(* ::Subsubsection::Closed:: *)
(*Triangle*)


FixFirstLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,firstLoopTrianglePattern] := Module[
{assocOut=<||>,l11=Unique["l11$"],l12=Unique["l12$"],l13=Unique["l13$"],
	Mu=Unique["Mu1$"],MuTilde=Unique["MuTilde1$"],t=Unique["t1$"], config},
	
	config = ConstructFirstLoopConfig[totalConfig];
	
	FixTriangleLoopMomentum6D[config,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l11,l12,l13},
			PassRules[FixFirstLoopMomentum,FixTriangleLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Triangle";
	assocOut["Momenta"] = {l11,l12,l13};
	assocOut["Variations"] = {$star,$unstar};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;
	
	assocOut
];


(* ::Subsubsection::Closed:: *)
(*Bubble*)


FixFirstLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,firstLoopBubblePattern] := Module[
{assocOut=<||>,l11=Unique["l11$"],l12=Unique["l12$"],
	Mu=Unique["Mu1$"],MuTilde=Unique["MuTilde1$"],t=Unique["t1$"],y=Unique["y1$"], config, yIntegrals},
	
	config = ConstructFirstLoopConfig[totalConfig];
	
	yIntegrals = FixBubbleLoopMomentum6D[config,OptionValue[BubbleReference],yVar->y,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l11,l12},
			PassRules[FixFirstLoopMomentum,FixBubbleLoopMomentum6D,opts]];
	
	assocOut["Type"] = "Bubble";
	assocOut["Momenta"] = {l11,l12};
	assocOut["Variations"] = {Null};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t,"y"->y|>;
	assocOut["yIntegrals"] = ToNum@yIntegrals;
	
	assocOut
];


(* ::Subsubsection::Closed:: *)
(*BubbleTriangle*)


FixFirstLoopMomentum[totalConfig_,opts:OptionsPattern[]] /; MatchQ[totalConfig,firstLoopBubbleTrianglePattern] := Module[
{assocOut=<||>,l11=Unique["l11$"],l12=Unique["l12$"],l13=Unique["l13$"],
	Mu=Unique["Mu1$"],MuTilde=Unique["MuTilde1$"],t=Unique["t1$"], baseConfig, subConfig, splittingMomenta, tIntegrals},

	baseConfig = ConstructFirstLoopBaseConfig[totalConfig];
	subConfig = ConstructFirstLoopConfig[totalConfig];
	
	splittingMomenta = GetSplittingMomenta[baseConfig,subConfig];
	
	tIntegrals = FixBubbleTriSubLoopMomentum6D[baseConfig,splittingMomenta,OptionValue[BubbleReference],tVar->t,
			MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l11,l12,l13},LargeTLimit->True,
			PassRules[FixFirstLoopMomentum,FixBubbleTriSubLoopMomentum6D,opts]];
	
	assocOut["Type"] = "BubbleTriangle";
	assocOut["Momenta"] = {l11,l12,l13};
	assocOut["Variations"] = {$plus,$minus};
	assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;
	assocOut["tIntegrals"] = ToNum@tIntegrals;
	
	assocOut
];


(* ::Section:: *)
(*Product of Trees*)


(* Return a list of the loop momenta in the order in which they appear in the loop. 
	In case of cut bubbles it may be necessary to rotate the momenta. *)
orderedFirstLoopMomenta[totalConfig_, firstLoopMom_]:=If[CutBottomFirstLoopQ[totalConfig],RotateRight@#,#]&@firstLoopMom["Momenta"];
orderedSecondLoopMomenta[totalConfig_, secondLoopMom_]:=If[CutTopSecondLoopQ[totalConfig],RotateRight@#,#]&@secondLoopMom["Momenta"];


Options[ConstructCenterTree] = {Gravity->False};
ConstructCenterTree[totalConfig_,firstLoopMom_,secondLoopMom_,opts:OptionsPattern[]]:=
Module[{plusTop,plusBottom, plusFirstInner, plusSecondInner, scalarLabelsLeft,scalarLabelsRight,TreeHead},

(*	TreeHead = If[OptionValue[Gravity]===True,Inactive@MTree,Inactive@ATree];*)
	TreeHead = Which[
		OptionValue[Gravity] === True,
		Inactive@MTree,
		ContainsAll[totalConfig["Center","Properties"],{"Gluon"}],
		Inactive@ATreeGluon,
		ContainsAll[totalConfig["Center","Properties"],{"Contact"}],
		Inactive@ATreeContact,
		ContainsAll[totalConfig["Center","Properties"],{"BerendsGiele"}],
		Inactive@ATreeBerendsGiele,
		True,
		Inactive@ATree
	];

	
	plusTop = ConstantArray["+",Length@totalConfig["Center"]["Top"]];
	plusBottom = ConstantArray["+",Length@totalConfig["Center"]["Bottom"]];
	plusFirstInner = ConstantArray["+",Length@totalConfig["Center"]["FirstInner"]];
	plusSecondInner = ConstantArray["+",Length@totalConfig["Center"]["SecondInner"]];
	
	scalarLabelsLeft = MapAt[pm,#,{All,2}]&@(Through[(orderedFirstLoopMomenta[totalConfig,firstLoopMom][[{1,-1}]])[#]]&/@firstLoopMom["Variations"]);
	scalarLabelsRight = MapAt[pm,#,{All,2}]&@(Through[(orderedSecondLoopMomenta[totalConfig,secondLoopMom][[{1,-1}]])[#]]&/@secondLoopMom["Variations"]);


	Which[
		BottomTwistedConfigQ[totalConfig],
			Outer[
				TreeHead
					["S",Sequence@@plusFirstInner,"s",Sequence@@plusTop,"S",Sequence@@plusSecondInner,"s",Sequence@@plusBottom]
					[#2[[2]],Sequence@@totalConfig["Center"]["FirstInner"],#1[[2]],Sequence@@totalConfig["Center"]["Top"],
					#2[[1]],Sequence@@totalConfig["Center"]["SecondInner"],#1[[1]],Sequence@@totalConfig["Center"]["Bottom"]]&,
				scalarLabelsLeft,
				scalarLabelsRight,
				1
			],
		TopTwistedConfigQ[totalConfig],
			Outer[
				TreeHead
					["S",Sequence@@plusFirstInner,"s",Sequence@@plusTop,"S",Sequence@@plusSecondInner,"s",Sequence@@plusBottom]
					[#1[[1]],Sequence@@totalConfig["Center"]["FirstInner"],#2[[1]],Sequence@@totalConfig["Center"]["Top"],
					#1[[2]],Sequence@@totalConfig["Center"]["SecondInner"],#2[[2]],Sequence@@totalConfig["Center"]["Bottom"]]&,
				scalarLabelsLeft,
				scalarLabelsRight,
				1
			],
		True,
			Outer[
				TreeHead
					["S",Sequence@@plusFirstInner,"S",Sequence@@plusTop,"s",Sequence@@plusSecondInner,"s",Sequence@@plusBottom]
					[#1[[1]],Sequence@@totalConfig["Center"]["FirstInner"],#1[[2]],Sequence@@totalConfig["Center"]["Top"],
					#2[[1]],Sequence@@totalConfig["Center"]["SecondInner"],#2[[2]],Sequence@@totalConfig["Center"]["Bottom"]]&,
				scalarLabelsLeft,
				scalarLabelsRight,
				1
			]
	]
];


Options[ConstructFirstLoopAmplitudes] = {Gravity->False};
ConstructFirstLoopAmplitudes[totalConfig_,firstLoopMom_Association,opts:OptionsPattern[]] /; 1+Length@totalConfig["FirstLoop"] === Length@firstLoopMom["Momenta"]:=
Module[{loopMomenta},
	loopMomenta = Through[(orderedFirstLoopMomenta[totalConfig,firstLoopMom])@#]&/@firstLoopMom[["Variations"]];
	ConstructLoopAmplitudes[totalConfig["FirstLoop"],#,PassRules[ConstructFirstLoopAmplitudes,ConstructLoopAmplitudes,opts]]&/@loopMomenta
];


Options[ConstructSecondLoopAmplitudes] = {Gravity->False};
ConstructSecondLoopAmplitudes[totalConfig_,secondLoopMom_Association,opts:OptionsPattern[]] /; 1+Length@totalConfig["SecondLoop"] === Length@secondLoopMom["Momenta"]:=
Module[{loopMomenta},
	loopMomenta = Through[(orderedSecondLoopMomenta[totalConfig,secondLoopMom])@#]&/@secondLoopMom[["Variations"]];
	ConstructLoopAmplitudes[totalConfig["SecondLoop"],#,PassRules[ConstructSecondLoopAmplitudes,ConstructLoopAmplitudes,opts]]&/@loopMomenta
];


(* ::Section:: *)
(*Rational Terms*)

evaluateAmplitudes[amps_] /; $forceBerendsGiele :=
	Module[{res},
	       $TwoLoopAllPlusRationalsStatus = "Evaluating trees";
	       res = Activate @ amps;
	       $TwoLoopAllPlusRationalsStatus = Null;
	       res
	];

evaluateAmplitudes[amps_] :=
	Module[{res},
	       $TwoLoopAllPlusRationalsStatus = "Evaluating trees";
	       res = ToNum @ Activate @ amps;
	       $TwoLoopAllPlusRationalsStatus = Null;
	       res
	];


Options[ComputeRationalTerm] = {DebugLevel->0,NumericVariables->{},DebugLogPath->"",Gravity->False};

ComputeRationalTerm[totalConfig_,firstTrees_,centerTree_,secondTrees_,firstLoopMom_,secondLoopMom_,opts : OptionsPattern[]] := 
Catch@Module[{secondLoopRational},
	
	secondLoopRational = Map[
		(PrintDebug["Second Loop...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->2]];
		
		ComputeSecondLoopRationalTerm[totalConfig, evaluateAmplitudes @ #,
			secondLoopMom,
			Values@firstLoopMom["Parameters"],
			PassRules[ComputeRationalTerm,ComputeSecondLoopRationalTerm,opts]])&,
		Inner[Times,centerTree,secondTrees,List],
	{2}];
		
	Transpose[
		Plus@@@Plus@@Map[
			(PrintDebug["First Loop...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->2]];
			
			ComputeFirstLoopRationalTerm[totalConfig, evaluateAmplitudes @ #,firstLoopMom,{},
				PassRules[ComputeRationalTerm,ComputeFirstLoopRationalTerm,opts]])&,
			Inner[Times,firstTrees,secondLoopRational,List],
		{3}]
	]
];


(* ::Subsection:: *)
(*First Loop Rational Term*)


Options[ComputeFirstLoopRationalTerm]={Gravity->False};


ComputeFirstLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,firstLoopBoxPattern] :=
	GetBoxRationalTerm[trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		LargeMuLimit->True,AnalyticVariables->analyticVariables,
		PassRules[ComputeFirstLoopRationalTerm,GetBoxRationalTerm,opts]];


ComputeFirstLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,firstLoopTrianglePattern] :=
	GetTriangleRationalTerm[trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],AnalyticVariables->analyticVariables,
		PassRules[ComputeFirstLoopRationalTerm,GetTriangleRationalTerm,opts]];


ComputeFirstLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,firstLoopBubblePattern] :=
	GetBubbleBubRationalTerm[trees,loopMom["yIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],yVar->loopMom["Parameters"]["y"],AnalyticVariables->analyticVariables,
		PassRules[ComputeFirstLoopRationalTerm,GetBubbleBubRationalTerm,opts]];


ComputeFirstLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,firstLoopBubbleTrianglePattern] :=
	GetBubbleTriRationalTerm[trees,loopMom["tIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],AnalyticVariables->analyticVariables,
		PassRules[ComputeFirstLoopRationalTerm,GetBubbleTriRationalTerm,opts]];


ComputeFirstLoopRationalTerm[___,opts:OptionsPattern[]]:=Throw[Null];


(* ::Subsection::Closed:: *)
(*Second Loop Rational Term*)


Options[ComputeSecondLoopRationalTerm]={Gravity->False};


ComputeSecondLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,secondLoopBoxPattern] :=
	GetBoxRationalTerm[trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		LargeMuLimit->True,AnalyticVariables->analyticVariables,
		PassRules[ComputeSecondLoopRationalTerm,GetBoxRationalTerm,opts]];


ComputeSecondLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,secondLoopTrianglePattern] :=
	GetTriangleRationalTerm[trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],AnalyticVariables->analyticVariables,
		PassRules[ComputeSecondLoopRationalTerm,GetTriangleRationalTerm,opts]];


ComputeSecondLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,secondLoopBubblePattern] :=
	GetBubbleBubRationalTerm[trees,loopMom["yIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],yVar->loopMom["Parameters"]["y"],AnalyticVariables->analyticVariables,
		PassRules[ComputeSecondLoopRationalTerm,GetBubbleBubRationalTerm,opts]];


ComputeSecondLoopRationalTerm[totalConfig_,trees_,loopMom_,analyticVariables_,opts:OptionsPattern[]] /; FreeQ[trees,ATree]&&MatchQ[totalConfig,secondLoopBubbleTrianglePattern] :=
	GetBubbleTriRationalTerm[trees,loopMom["tIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
		tVar->loopMom["Parameters"]["t"],AnalyticVariables->analyticVariables,
		PassRules[ComputeSecondLoopRationalTerm,GetBubbleTriRationalTerm,opts]];


ComputeSecondLoopRationalTerm[___,opts:OptionsPattern[]]:=Throw[Null];


(* ::Section::Closed:: *)
(*Config Hash*)


TotalConfigHash[totalConfig_]:=Hash[{totalConfig,Mom4DN[#][Null][$up]&/@GetAllMomenta@totalConfig},"SHA256"];


TotalConfigHashString[totalConfig_]:=ToString@TotalConfigHash[totalConfig];


(* ::Section::Closed:: *)
(*Writing to Config (Failed) File*)


WriteToFailedFile[message_String,failedFile_String]:=
With[{outOptions = Options[$Output]}, 
		Internal`WithLocalSettings[
			SetOptions[$Output,PageWidth->Infinity], 
			PutAppend[message,failedFile];,
			SetOptions[$Output,outOptions]
		];
	];


(* ::Section::Closed:: *)
(*Compute Configuration*)


Options[ComputeConfiguration] = {
	DebugLevel -> 0,
	NumericVariables -> {},
	PostProcessing -> False,
	DataDirectory -> FileNameJoin[{Directory[],"data"}],
	LogDirectory -> FileNameJoin[{Directory[],"logs"}],
	BubbleReference -> $bubbleReference,
	ComputeLoopMomentumSpinors -> True,
	Recompute -> False,
	LoopMomentumPostProcessingFirstPass -> Identity,
	NumSymbReplacement -> True,
	Gravity->False,
	Prefactor->1,
	AnalyticParameters->{},
	ReduceSqrts->True
};


ReuseResultQ[totalConfig_, opts:OptionsPattern[ComputeConfiguration]]:= 
	(! OptionValue[Recompute]) && FileExistsQ[
		FileNameJoin@{
			OptionValue[DataDirectory],
			TotalConfigHashString[totalConfig]
		}
	];


ComputeConfiguration[totalConfig_, opts : OptionsPattern[]] /; ReuseResultQ[totalConfig,opts] := 
Module[{out,cachedFunctions,cachedDownValues,dataFile, hash,failedFile},
	Print["Found file for ",totalConfig];
	(* saving default DownValues of caching functions *)
	cachedFunctions = {RationalSeries`LowestPowerSeries,RationalSeries`RationalSeries,RationalSeries`SampleEval,
					  RationalSeries`ZeroCoefficientQ,GetLowestOrder};
	cachedDownValues = DownValues/@cachedFunctions;
	
	
	dataFile=FileNameJoin@{OptionValue[DataDirectory],TotalConfigHashString[totalConfig]};
	Get[dataFile];
	
	$currentConfiguration = ExtractSymbolName@totalConfig;

	hash = TotalConfigHash[totalConfig];

	failedFile = FileNameJoin@{OptionValue[LogDirectory],"failed",ToString[hash]<> "_kernel"<>ToString[$KernelID] <> "_" <> $MachineName};

	WriteToFailedFile[ToString@totalConfig,failedFile];
	WriteToFailedFile["Found cached result...",failedFile];
	
	Which[
		!MatchQ[$rationalFinal,_Symbol],
		out = $rationalFinal;
	      ,
		!MatchQ[$rationalFull,_Symbol],
		out = $rationalFinal = PostprocessCutResult[ComputeConfiguration,$rationalFull,PassRules[ComputeConfiguration,PostprocessCutResult,opts]];
		
		DumpSave[dataFile,{$rational,$rationalFull,$rationalFinal,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];
	      ,
		!MatchQ[$rational,_Symbol],
		$rationalFull = OptionValue[Prefactor] *
				Switch[
					OptionValue[NumSymbReplacement],
					True,
					$rational//.NumSymb->NumSymbN,
					_,
					$rational
				];
		DumpSave[dataFile,{$rational,$rationalFull,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];
		
		out = $rationalFinal = PostprocessCutResult[ComputeConfiguration,$rationalFull,PassRules[ComputeConfiguration,PostprocessCutResult,opts]];
		
		DumpSave[dataFile,{$rational,$rationalFull,$rationalFinal,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];
	      ,
		_,
		Message[TwoLoopAllPlusRationalsTermConfigs::InvFile,totalConfig];
		out = $Failed[totalConfig];
	];

	DeleteFile[failedFile];
	$currentConfiguration = Null;
	
	(* resetting caching functions *)
	ClearAll[NumSymb,NumSymbN,NumSymbSimplified,NumSymbEvaluated,$rational,$rationalFull,$rationalFinal,$currentConfiguration];
	NumSymbLookUp = <||>;
	MapThread[(DownValues[#1] = #2;)&,{cachedFunctions,cachedDownValues}];
	
	out
];


ComputeConfiguration[totalConfig_,opts:OptionsPattern[]] :=
Module[{firstLoopMom, secondLoopMom, firstTrees, centerTree, secondTrees,
		out,cachedFunctions,cachedDownValues,dataFile,hash,failedFile},

	(* saving default DownValues of caching functions *)
	cachedFunctions = {RationalSeries`LowestPowerSeries,RationalSeries`RationalSeries,RationalSeries`SampleEval,
					  RationalSeries`ZeroCoefficientQ,GetLowestOrder};
	cachedDownValues = DownValues/@cachedFunctions;
	
	(* Clearing any information of previous runs *)
	ClearAll[NumSymb,NumSymbN,NumSymbSimplified,NumSymbEvaluated,$rational,$rationalFull,$rationalFinal,$currentConfiguration];
	NumSymbLookUp = <||>;
	
	hash = TotalConfigHash[totalConfig];
	
	dataFile = FileNameJoin@{OptionValue[DataDirectory],ToString@hash};
	failedFile = FileNameJoin@{OptionValue[LogDirectory],"failed",ToString[hash]<> "_kernel"<>ToString[$KernelID] <> "_" <> $MachineName};
	
	WriteToFailedFile[ToString@totalConfig,failedFile];
	WriteToFailedFile["Setting up loop-momenta and trees...",failedFile];

	PrintDebug["Fixing loop momenta...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	{firstLoopMom,secondLoopMom} = FixLoopMomenta[totalConfig,PassRules[ComputeConfiguration,FixLoopMomenta,opts]];
	
	$currentConfiguration = ExtractSymbolName@totalConfig;
	
	PrintDebug["Generating trees...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	
	(* Check whether all amplitudes are present, otherwise throw $missingConfig as result*)
	
	{firstTrees,centerTree,secondTrees} = 
		{
			ConstructFirstLoopAmplitudes[totalConfig,firstLoopMom,PassRules[ComputeConfiguration,ConstructFirstLoopAmplitudes,opts]],
			ConstructCenterTree[totalConfig,firstLoopMom,secondLoopMom,PassRules[ComputeConfiguration,ConstructCenterTree,opts]],
			ConstructSecondLoopAmplitudes[totalConfig,secondLoopMom,PassRules[ComputeConfiguration,ConstructSecondLoopAmplitudes,opts]]
		};

	If[OptionValue[DebugLevel] > 0,
		$debugFirstLoopMom = firstLoopMom;
		$debugSecondLoopMom = secondLoopMom;
		$debugFirstTrees = firstTrees;
		$debugCenterTrees = centerTree;
		$debugSecondTrees = secondTrees;
	];


	WriteToFailedFile["Computing rational terms...",failedFile];
	PrintDebug["Computing Rational Term...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];

	Check[
		$rational = ComputeRationalTerm[totalConfig,firstTrees,centerTree,secondTrees,firstLoopMom,secondLoopMom,
			PassRules[ComputeConfiguration,ComputeRationalTerm,opts]];,
		(
			DeleteFile[failedFile];
			Message[ComputeConfiguration::missingAmp, $currentConfiguration];
			Throw[$missingConfig[$currentConfiguration]];
		),
		ATree::noamp
	];

	PrintDebug["Finished ComputeRationalTerm, saving $rational...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	DumpSave[dataFile,{$rational,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];


	WriteToFailedFile["Replacing NumSymbs...",failedFile];
	PrintDebug["NumSymbs:",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	
	$rationalFull = OptionValue[Prefactor] * Switch[
		OptionValue[NumSymbReplacement],
		True,
			PrintDebug["Replacing",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->2]];
			$rational//.NumSymb->NumSymbN,
		_,
			PrintDebug["No replacing",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->2]];
			$rational
	];
	
	PrintDebug["Saving $rationalFull...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	DumpSave[dataFile,{$rational,$rationalFull,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];


	WriteToFailedFile["Post-processing...",failedFile];

	

	out = $rationalFinal = PostprocessCutResult[$rationalFull,OptionValue[PostProcessing],PassRules[ComputeConfiguration,PostprocessCutResult,opts]];


	PrintDebug["Saving $rationalFinal...",Append[PassRules[ComputeConfiguration,PrintDebug,opts],Indent->1]];
	DumpSave[dataFile,{$rational,$rationalFull,$rationalFinal,NumSymb,NumSymbN,NumSymbSimplified,NumSymbLookUp,$currentConfiguration}];

	WriteToFailedFile["Done. Cleaning up...",failedFile];

	If[OptionValue[DebugLevel]<=0,
		(* resetting caching functions *)
		ClearAll[NumSymb,NumSymbN,NumSymbSimplified,NumSymbEvaluated,$rational,$rationalFull,$rationalFinal,$currentConfiguration];
		NumSymbLookUp = <||>;
		MapThread[(DownValues[#1] = #2;)&,{cachedFunctions,cachedDownValues}];
	];

	DeleteFile[failedFile];
	$currentConfiguration = Null;

	out
];

(* ::Section::Closed:: *)
(* PostprocessCutResult *)

Options[PostprocessCutResult]=Options[ComputeConfiguration];

PostprocessCutResult[cutResult_,postProcessOption:True,opts:OptionsPattern[]]:=
	(
		PrintDebug["Postprocessing: RootReduce...",Append[PassRules[PostprocessCutResult,PrintDebug,opts],Indent->2]];
		RootReduce@$rationalFull
	);

PostprocessCutResult[cutResult_,postProcessOption:EvaluationFunction,opts:OptionsPattern[]]:=
	(
		PrintDebug["Postprocessing: EvaluationFunction...",Append[PassRules[PostprocessCutResult,PrintDebug,opts],Indent->2]];
		If[
			!MatchQ[OptionValue[AnalyticParameters],{__}],
			Message[ComputeConfiguration::EvalNoPars];
		];
		GenerateCutEvaluationFunction[$rationalFull,OptionValue[AnalyticParameters]]
	);

PostprocessCutResult[cutResult_,postProcessOption:_Function|_Symbol|_Composition,opts:OptionsPattern[]]:=
	(
		PrintDebug["Postprocessing: Function passed " <> ToString[postProcessOption] <> " via option PostProcessing...",Append[PassRules[PostprocessCutResult,PrintDebug,opts],Indent->2]];
		postProcessOption@cutResult
	);

PostprocessCutResult[cutResult_,postProcessOption:_,opts:OptionsPattern[]]:=
	(
		PrintDebug["No Postprocessing",Append[PassRules[PostprocessCutResult,PrintDebug,opts],Indent->2]];
		cutResult
	);


(* ::Section::Closed:: *)
(*Compute Rational Terms*)


Options[TwoLoopAllPlusRationalsTermConfigs] = Options[ComputeConfigurations];


TwoLoopAllPlusRationalsTermConfigs[momenta_List,opts:OptionsPattern[]]:=
	ComputeConfigurations[
		Join@@GenerateComputationTotalConfigs[momenta]/@CyclicPermutations[momenta],
	opts];


Options[ComputeConfigurations] = {
	DebugLevel->0,
	ComputeLoopMomentumSpinors->True,
	NumSymbReplacement -> True,
	PostProcessing -> False,
	LogDirectory -> FileNameJoin[{Directory[],"logs"}],
	DataDirectory -> FileNameJoin[{Directory[],"data"}],
	BubbleReference -> $bubbleReference,
	Recompute->False,
	LoopMomentumPostProcessingFirstPass -> Identity,
	ShuffleConfigs -> False,
	ResultFormat->Plus,
	Gravity->False,
	WriteResult->True,
	ApplyIntegrals->True,
	Prefactor->1,	
	AnalyticParameters->{},
	ReduceSqrts->True
};


ComputeConfigurations[totalConfigsInput_List,opts:OptionsPattern[]]:=
    Module[{totalConfigs,logDir,totalConfigsEval,result,momenta,logFile},
	
	SetOptions[$Output,FormatType->OutputForm];
	
	totalConfigs = Append[#,#[[2]]]&/@totalConfigsInput;
	
	If[OptionValue[ShuffleConfigs],
		totalConfigs = RandomSample[totalConfigs];
	];
	
	momenta = totalConfigs[[1,3]]//GetAllMomenta;	
	
	If[StringQ@OptionValue[LogDirectory],
		logDir = FileNameJoin@{OptionValue[LogDirectory],DateString[{"Year", "-", "Month", "-", "Day", "-", "Hour", "-", "Minute", "-", "Second" }] <> "_" <> 
				ToString[Length@momenta] <> "pt_kernel_" <> ToString[$KernelID] <> "_machine_" 
				<> ToString[$MachineName] <> "_momenta_" <> StringJoin@@StringReplace[ToString/@momenta," "->""]};

		CreateDirectory[logDir,CreateIntermediateDirectories -> True];
		CreateDirectory[FileNameJoin@{logDir,"failed"}];

		logFile = OpenWrite[FileNameJoin@{logDir,"mathematica_output.log"},FormatType -> OutputForm];
		AppendTo[$Output, logFile];
		AppendTo[$Messages, logFile];
	];

	If[!DirectoryQ@OptionValue[DataDirectory],
		CreateDirectory[OptionValue[DataDirectory],CreateIntermediateDirectories -> True];
	];
	
	totalConfigsEval = MapAt[
		ParallelSubmit[{logDir},
			Catch@ComputeConfiguration[#, LogDirectory -> logDir,
			PassRules[ComputeConfigurations,ComputeConfiguration,{opts}]]
		]&,
		totalConfigs,
		{All,3}
	];

	If[OptionValue[PostProcessing] === EvaluationFunction,
	   (* pure Coefficient matrix *)
	   result = ApplyMuIntegralEvaluationFunction[#,PassRules[ComputeConfigurations,ApplyMuIntegralEvaluationFunction,opts]]&/@WaitAll[totalConfigsEval],

	   
	   (* symmetry factor times integralsLeft.Coefficients.integralsRight *)
	   result = OptionValue[ResultFormat]@@ToNum[
		   (#[[1]]*ApplyMuIntegrals[#[[2]],#[[3]],PassRules[ComputeConfigurations,ApplyMuIntegrals,opts]]) & /@ WaitAll[totalConfigsEval]]
	];

	Switch[OptionValue[WriteResult],
		"MX",
		$CutResult = result;
		DumpSave[FileNameJoin[{logDir,DateString[{"Year", "-", "Month", "-", "Day", "-", "Hour", "-", "Minute", "-", "Second" }] <> "_" <> ToString[Length@momenta] <> "pt_result.mx"}],$CutResult];
		ClearAll[$CutResult];,
		"M",
		Put[result,FileNameJoin[{logDir,DateString[{"Year", "-", "Month", "-", "Day", "-", "Hour", "-", "Minute", "-", "Second" }] <> "_" <> ToString[Length@momenta] <> "pt_result.m"}]];,
		_,
		Null;
	];

	result
    ];


Options[ApplyMuIntegralEvaluationFunction]={Gravity->False,AnalyticParameters->{}};

ApplyMuIntegralEvaluationFunction[{sym_,config_,cutMatrix_},opts:OptionsPattern[]]:=
	Module[{matrixDummy,
		integralProd,
		pars,
		parameters = OptionValue[AnalyticParameters]
	       },
	       
	       integralProd=ToNum[ApplyMuIntegrals[config,matrixDummy,PassRules[ApplyMuIntegralEvaluationFunction,ApplyMuIntegrals,opts]]];

	       integralProd = integralProd/. MapThread[Rule,{parameters,Table[Indexed[pars,i],{i, Length@parameters}]}];
       
	       Function@@(Hold[{pars},Block[{matrixDummy},matrixDummy=#2[pars];#1*#3]]&[sym,cutMatrix,integralProd])
	];


(* ::Chapter:: *)
(*Postamble*)


End[];


Protect[ "SpinorHelicity6D`TwoLoopAllPlusRationals`*" ];    (* protect exported symbols *)
Unprotect[
	$FullConfigPattern,
	$rational,
	$rationalFull,
	$rationalFinal,
	$currentConfiguration,
	$debugFirstLoopMom,
	$debugSecondLoopMom,
	$debugFirstTrees,
	$debugCenterTrees,
	$debugSecondTrees,
	$CutResult,
	$TwoLoopAllPlusRationalsStatus
];

EndPackage[ ];  (* end the package context *)
