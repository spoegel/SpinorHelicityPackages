(* ::Package:: *)

(* ::Title:: *)
(*Helpers*)


(* :Title: Helpers.m -- all helper functions used in both the full permutation and config-by-config implementations *)

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
	TwoLoopAllPlusRationals.wl base
*)

(* :Todo:
*)


(* ::Chapter:: *)
(*Preamble*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`",{
	"SpinorHelicity6D`Unitarity`",
	"SpinorHelicity6D`SpinorFunctions`",
	"SpinorHelicity6D`Amplitudes`",
	"SpinorHelicity6D`TwoLoopAllPlusRationals`",
	"SpinorHelicity6D`",
	"RationalSeries`",
	"Utils`",
	"NumSymb`",
	"SpinorHelicity6D`Unitarity`LoopMomentum`",
	"SpinorHelicity6D`Unitarity`OneLoopIntegrals`",
	"SpinorHelicity6D`Unitarity`ConfigurationTools`",
	"SpinorHelicity6D`Unitarity`RationalTerms`"}];

Unprotect["SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`*"];
ClearAll["SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`*"];
ClearAll["SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`Private`*"];



(* ::Section::Closed:: *)
(*Debug*)


PrintDebug::usage = "";


DebugExprOut::usage = "";


createMemoryWatchdog::usage = "";


createRationalSeriesWatchdog::usage = "";


$LogFileTemplate = "";


$debugMode::usage = "";
$debugMode = False;


(* ::Section:: *)
(*Configurations*)


(* ::Subsection::Closed:: *)
(*Translate From Old to New TotalConfig Notation*)


TranslateOldTotalConfig::usage = "";


(* ::Subsection::Closed:: *)
(*Loop Patterns & TwistedQ*)


firstLoopBoxPattern::usage = "";
secondLoopBoxPattern::usage = "";


firstLoopTrianglePattern::usage = "";
secondLoopTrianglePattern::usage = "";


firstLoopBubblePattern::usage = "";
secondLoopBubblePattern::usage = "";


firstLoopBubbleTrianglePattern::usage = "";
secondLoopBubbleTrianglePattern::usage = "";


TwistedConfigQ::usage = "";


TopTwistedConfigQ::usage = "";
BottomTwistedConfigQ::usage = "";


(* ::Subsection::Closed:: *)
(*Marking Corners*)


GetFirstLoopMomenta::usage = "";


MarkBubbleCutCorner::usage = "";


MarkSecondLoopBubbleCut::usage = "";


MarkSecondLoopBubbleCut::usage = "";


UnmarkConfig::usage = "";


(* ::Subsection::Closed:: *)
(*Generating Configurations*)


GenConfigBare::usage = "";
GenConfigsBare::usage = "";
GenConfigPartitions::usage = "";
GenSubConfigs::usage = "";
GetSplittingMomenta::usage = "";


GenFirstLoopConfigs::usage = "";
GenSecondLoopConfigs::usage = "";
CompleteSecondLoopConfig::usage = "";
CutInLastCornerQ::usage = "";
CutInFirstCornerQ::usage = "";
CutCorners::usage = "";
GenTotalConfig::usage = "";

CreateTotalConfigTemplates::usage = "";
CreateTwistedCentralVertex::usage = "";
DistributeTemplateFreeToFreeInnerFreeOuter::usage = "";
DistributeSubleadingSingleTrace::usage = "";
FilterDoubleTriangleTemplates::usage = "";

(* ::Subsubsection::Closed:: *)
(*GenTotalConfig*)


GenAllTotalConfigs::usage = "";


CutSecondLoop::usage = "";


CutFirstLoop::usage = "";


genConfigRestore$cutBubbleLabels::usage = "";


CutBubbles::usage="";
CutBubble::usage="";


(* ::Subsection::Closed:: *)
(*Constructing Single Loop Configs for Loop Momenta*)


GetAllMomenta::usage = "";
GetFirstLoopMomenta::usage = "";
GetSecondLoopMomenta::usage = "";
GetCenterMomenta::usage = "";


(* ::Subsubsection::Closed:: *)
(*Second Loop*)


ConstructSecondLoopConfig::usage = "";


ConstructSecondLoopBaseConfig::usage = "";


(* ::Subsubsection:: *)
(*First Loop*)


ConstructFirstLoopConfig::usage = "";


ConstructFirstLoopBaseConfig::usage = "";





(* ::Section:: *)
(*Automatic Config Generation*)


GenerateTotalConfigs::usage = "";
GenerateComputationTotalConfigs::usage = "";

GenerateTotalConfigs::duplc = "Labels `1` contain duplicates.";
GenerateComputationTotalConfigs::duplc = "Labels `1` contain duplicates.";


(* ::Subsection:: *)
(*CutBubble Helpers*)


CutOuterVertex::usage = "";


CutTopFirstLoopQ::usage = "";
CutBottomFirstLoopQ::usage = "";


CutTopSecondLoopQ::usage = "";
CutBottomSecondLoopQ::usage = "";


fuseRestore$cutBubbleLabels::usage = "";


FuseSecondLoopCutBubble::usage = "";


FuseFirstLoopCutBubble::usage = "";


FuseCutBubbles::usage = "";


(* ::Subsection:: *)
(*Subleading Single Trace*)


GenerateSubLeadingSingleTraceConfigurations::usage = "";
GenerateComputeSubLeadingSingleTraceConfigurations::usage = "";


(* ::Section::Closed:: *)
(*PostProcessing*)


NestedReplaceNumSymbs::usage = "";



(* ::Section::Closed:: *)
(*Symmetry Factors*)


ComputeSymmetryFactor::usage = "";

$FullConfigPattern::usage = "Pattern to single out configurations to compute. Default: HoldPattern[_]";
$FullConfigPattern = HoldPattern[_];


(* ::Section::Closed:: *)
(*Integrals*)


ConfigIntegral::usage = "";
ApplyMuIntegrals::usage = "";

(* ::Chapter:: *)
(*Main*)


Begin["`Private`"];    (* begin the private context (implementation part) *)


(* ::Section::Closed:: *)
(*Debug*)


SetAttributes[PrintDebug,HoldFirst];
Options[PrintDebug] = {DebugLevel -> 0, Indent -> 0};


(* The Indent conditional is for backwards compatibility for StringRepeat *)
PrintDebug[args___,opts:OptionsPattern[]] /; $debugMode || OptionValue[PrintDebug,opts,DebugLevel] > 0 := Print[
	If[OptionValue[Indent] === 0, "", StringRepeat["\t",OptionValue[Indent]]],
	Sequence@@ReleaseHold/@{args}
];

PrintDebug[___,opts:OptionsPattern[]] := Null;


SetAttributes[DebugExprOut,HoldFirst];


DebugExprOut[expr_,fileName_String] /; $debugExpressionMode && $LogFileTemplate =!= "" := 
	Put[Evaluate[ToNum@expr],TemplateApply[$LogFileTemplate,<|"type"->fileName,"ending"->"ms"|>]];


createMemoryWatchdog[] :=
SessionSubmit@ScheduledTask[
	Module[{file=TemplateApply[$LogFileTemplate,<|"type"->"memory","ending"->"log"|>]},
		WriteString[file,
			"Kernel " <> ToString[$KernelID] <> " on " <> ToString[$MachineName] <> ":\n" <> 
			"\tMathematica Use: " <> ToString[MemoryInUse[]] <> "\n" <>
			"\tAvailable: " <> ToString[MemoryAvailable[]]<>
			"\n"];
		Close[file];
	];
,5];


createRationalSeriesWatchdog[] :=
SessionSubmit@ScheduledTask[
	Module[{file=TemplateApply[$LogFileTemplate,<|"type"->"rationalSeries","ending"->"log"|>]},
		WriteString[file,
			"Kernel " <> ToString[$KernelID] <> " on " <> ToString[$MachineName] <> ":\n" <> 
			"\tTimings Use: " <> ToString[$TotalTimingsRationalSeries] <> "\n" <>
			"\tCalls: " <> ToString[$TotalCallsRationalSeries] <>
			"\n"];
		Close[file];
	];
,5];


(* ::Section::Closed:: *)
(*Configurations*)


(* ::Subsection::Closed:: *)
(*Translate From Old to New TotalConfig Notation*)


translateCorner[corner_]:=<|"Outer"->corner/.{$cutBubble->Nothing},"Inner"->{},"Properties"->Cases[corner,$cutBubble]|>;


TranslateOldTotalConfig[totalConfigOld_]:=<|
	"FirstLoop"->translateCorner/@totalConfigOld[[1,1]],
	"SecondLoop"->translateCorner/@totalConfigOld[[2,1]],
	"Center"-><|"Top"->totalConfigOld[[1,2]]/.{$cutBubble->Nothing},"Bottom"->totalConfigOld[[2,2]]/.{$cutBubble->Nothing},"FirstInner"->{},"SecondInner"->{},
		"Properties"->DeleteDuplicates@Cases[Flatten[Join[totalConfigOld[[1,2]],totalConfigOld[[2,2]]]],$cutBubble]|>
|>;


(* ::Subsection::Closed:: *)
(*Loop Patterns & TwistedQ*)


firstLoopBoxPattern = <|___,"FirstLoop"->{_,_,_},___|> | {{{_,_,_},{___}},{{__},{___}}};
secondLoopBoxPattern = <|___,"SecondLoop"->{_,_,_},___|> | {{{__},{___}},{{_,_,_},{___}}};


firstLoopTrianglePattern = HoldPattern[<|___,"FirstLoop"->{a_,b_},___|> /; FreeQ[{a,b},$cutBubble]] | HoldPattern[{{{a_,b_},{___}},{{__},{___}}} /; FreeQ[{a,b},$cutBubble]];;
secondLoopTrianglePattern = HoldPattern[<|___,"SecondLoop"->{a_,b_},___|> /; FreeQ[{a,b},$cutBubble]] | HoldPattern[{{{__},{___}},{{a_,b_},{___}}} /; FreeQ[{a,b},$cutBubble]];


firstLoopBubblePattern = <|___,"FirstLoop"->{_},___|> | {{{_},{___}},{{__},{___}}};
secondLoopBubblePattern = <|___,"SecondLoop"->{_},___|> | {{{__},{___}},{{_},{___}}};


firstLoopBubbleTrianglePattern = HoldPattern[<|___,"FirstLoop"->{a_,b_},___|> /; !FreeQ[{a,b},$cutBubble]] | HoldPattern[{{{a_,b_},{___}},{{__},{___}}} /; !FreeQ[{a,b},$cutBubble]];
secondLoopBubbleTrianglePattern = HoldPattern[<|___,"SecondLoop"->{a_,b_},___|> /; !FreeQ[{a,b},$cutBubble]] | HoldPattern[{{{__},{___}},{{a_,b_},{___}}} /; !FreeQ[{a,b},$cutBubble]];


TwistedConfigQ[totalConfig_]:=TopTwistedConfigQ[totalConfig]||BottomTwistedConfigQ[totalConfig];


TopTwistedConfigQ[totalConfig_]:=totalConfig["Center","Twisted"]==="Top";
BottomTwistedConfigQ[totalConfig_]:=totalConfig["Center","Twisted"]==="Bottom";


(* ::Subsection:: *)
(*Combinatoric Functions*)




GenFirstLoopConfigs[nBlocks_,momenta_]:=Module[{configs},
	configs = Select[GenConfigsBare[nBlocks,momenta],MemberQ[#,{momenta[[1]],___}]&];
	configs = RotateLeft[#,Position[#,momenta[[1]]][[1,1]]-1]&/@configs;
	Select[configs,Length@Last[#]>1&]
];


GenSecondLoopConfigs[nCorners_,momenta_] := DeleteCases[GenConfigPartitions[nCorners+1,{"S",Sequence@@momenta,"Sb"}],{{__},{_},{__}}];


CompleteSecondLoopConfig[secondLoopConfig_List,pChannel_List] := Prepend[secondLoopConfig[[2;;-2]],Join[Last@secondLoopConfig,pChannel,First@secondLoopConfig]/.{"S"->Nothing,"Sb"->Nothing}];


CutInLastCornerQ[subconfig_,config_] := !MemberQ[subconfig,config[[-1]]];
CutInFirstCornerQ[subconfig_,config_] := !MemberQ[subconfig,config[[1]]];


CutCorners[subconfig_,config_] := Sort@Flatten[Position[subconfig,#]&/@Complement[subconfig,config]];


GenTotalConfig[configFirstLoop_,secondLoopConfig_] := {{configFirstLoop,Rest@First@secondLoopConfig},{secondLoopConfig[[2;;-2]],Most@Last@secondLoopConfig}};


(* ::Subsection::Closed:: *)
(*Constructing Single Loop Configs for Loop Momenta*)


GetCenterMomenta[totalConfig_]:=Join[#["Top"],#["SecondInner"],#["FirstInner"],#["Bottom"]]&@totalConfig["Center"];


GetFirstLoopMomenta[totalConfig_]:=Join[#["Outer"],#["Inner"]]&/@totalConfig["FirstLoop"];
GetSecondLoopMomenta[totalConfig_]:=Join[#["Outer"],#["Inner"]]&/@totalConfig["SecondLoop"];


GetAllMomenta[totalConfig_]:=Flatten@Through[{GetFirstLoopMomenta,GetCenterMomenta,GetSecondLoopMomenta}[totalConfig]];


(* ::Subsubsection::Closed:: *)
(*Second Loop*)


ConstructSecondLoopConfig[totalConfig_] := 
Module[{firstloop,secondLoopConfig},
	firstloop = Join[
		Flatten@GetFirstLoopMomenta@totalConfig,
		GetCenterMomenta@totalConfig
	];
	secondLoopConfig=Append[GetSecondLoopMomenta@totalConfig,firstloop];
	
	If[CutTopSecondLoopQ[totalConfig], secondLoopConfig = RotateLeft@secondLoopConfig;];
	
	secondLoopConfig	
];


ConstructSecondLoopBaseConfig[totalConfig_] /; MatchQ[totalConfig,secondLoopBubbleTrianglePattern] := 
	ConstructSecondLoopConfig@FuseSecondLoopCutBubble@totalConfig;


(* ::Subsubsection::Closed:: *)
(*First Loop*)


ConstructFirstLoopConfig[totalConfig_] := 
Module[{secondloop,firstLoopConfig},
	secondloop = Join[
		Flatten@GetSecondLoopMomenta@totalConfig,
		GetCenterMomenta@totalConfig
	];
	
	firstLoopConfig = Append[GetFirstLoopMomenta@totalConfig,secondloop];

	If[CutBottomFirstLoopQ[totalConfig], firstLoopConfig = RotateLeft@firstLoopConfig;];
	
	firstLoopConfig	
];


ConstructFirstLoopBaseConfig[totalConfig_] /; MatchQ[totalConfig,firstLoopBubbleTrianglePattern] := 
	ConstructFirstLoopConfig@FuseFirstLoopCutBubble@totalConfig;


(* ::Section:: *)
(*Automatic Config Generation*)


(* ::Subsection::Closed:: *)
(*Creating TotalConfigTemplates*)


VertexPartitions[vFirst_,vSecond_,n_] /; vFirst + vSecond > n := {};
VertexPartitions[vFirst_,vSecond_,n_] /; vFirst + vSecond == n := {Join[ConstantArray[1,vFirst],{0},ConstantArray[1,vSecond]]};
VertexPartitions[vFirst_,vSecond_,n_]:=Join[
	Flatten[Permutations/@IntegerPartitions[n,{vFirst+1+vSecond}],{1,2}],
	Insert[#,0,vFirst+1]&/@Flatten[Permutations/@IntegerPartitions[n,{vFirst+vSecond}],{1,2}]
];


LegDistributions[vFirst_,vSecond_,n_]:=
	RemoveSymmetricDuplicates[vFirst,vSecond]@RemoveScalessBubbles[vFirst,vSecond]@VertexPartitions[vFirst,vSecond,n];


RemoveScalessBubbles[1,1] := Cases[{_?(#>1&),_,_?(#>1&)}];
RemoveScalessBubbles[1,_] := Cases[{_?(#>1&),__}];
RemoveScalessBubbles[vFirst_,1] := Cases[{Repeated[_,{vFirst}],_,_?(#>1&)}];
RemoveScalessBubbles[_,_] := Identity;


CyclicSymmetricLegDistributionsQ[legDistribution1_,legDistribution2_,vBoth_]:=
	legDistribution1[[vBoth+1]]===legDistribution2[[vBoth+1]] &&
	legDistribution1[[;;vBoth]]===legDistribution2[[vBoth+2;;]] && 
	legDistribution2[[;;vBoth]]===legDistribution1[[vBoth+2;;]];
	
CyclicSymmetricLegDistributionsQ[vBoth_]:= Function[{legDistribution1,legDistribution2},CyclicSymmetricLegDistributionsQ[legDistribution1,legDistribution2,vBoth]];


RemoveSymmetricDuplicates[vBoth_,vBoth_]:=Function[legDistribution,DeleteDuplicates[legDistribution,CyclicSymmetricLegDistributionsQ[vBoth]]];
RemoveSymmetricDuplicates[vFirst_,vSecond_] := Identity;


CreateVertexTemplate[n_]:=
	ReplacePart[
		<|"Outer"->{},"Inner"->{},"Free"->Null,"VertexIndex"->Null,"Properties"->{}|>,
		{"Free"->n,"VertexIndex"->Unique[]}
	];


TotalConfigTemplateFromLegDistribution[legDistribution_,vFirst_,vSecond_]:=
	ReplacePart[
		<|
			"FirstLoop"->{},
			"SecondLoop"->{},
			"Center"-><|"Top"->{},"Bottom"->{},"FirstInner"->{},"SecondInner"->{},"Free"->Null,"VertexIndex"->Null,"Twisted"->Null,"Properties"->{}|>
		|>,
	{
		{"FirstLoop"}->Table[CreateVertexTemplate[legDistribution[[i]]],{i,1,vFirst}],
		{"SecondLoop"}->Table[CreateVertexTemplate[legDistribution[[i]]],{i,-vSecond,-1}],
		{"Center","Free"}->legDistribution[[vFirst+1]],{"Center","VertexIndex"}->Unique[]
	}
];


CreateTotalConfigTemplates[n_]:=
	Flatten@Table[
		TotalConfigTemplateFromLegDistribution[#,vFirst,vSecond]&/@LegDistributions[vFirst,vSecond,n],
	{vFirst,1,3},{vSecond,1,vFirst}];


(* ::Subsection:: *)
(*Filling TotalConfig Templates*)


(* ::Subsubsection::Closed:: *)
(*Creating Degenerate List of Vertices*)


LoopVertexListDegenerate[totalConfig_,loop:"FirstLoop"|"SecondLoop"]:=
	Flatten[
		Table[ConstantArray[{Key[loop],i},totalConfig[[loop,i,"Free"]]],{i,Length@totalConfig[loop]}],
	{1,2}];

CenterVertexListDegenerate[totalConfig_]:=ConstantArray[{Key["Center"]},totalConfig["Center","Free"]];


ConfigVertexListDegenerateOuterLoop[totalConfig_]:=
	Join[
		LoopVertexListDegenerate[totalConfig,"FirstLoop"],
		LoopVertexListDegenerate[totalConfig,"SecondLoop"],
		CenterVertexListDegenerate[totalConfig]
	];

ConfigVertexListDegenerateInnerLoop[loop:"FirstLoop"|"SecondLoop"] :=
    Function[totalConfig,ConfigVertexListDegenerateInnerLoop[totalConfig,loop]];
ConfigVertexListDegenerateInnerLoop[totalConfig_,loop:"FirstLoop"|"SecondLoop"] :=
	Join[
		LoopVertexListDegenerate[totalConfig,loop],
		CenterVertexListDegenerate[totalConfig]
	];


(* ::Subsubsection::Closed:: *)
(*Outer Assignment*)


GetOuterAssignments[totalConfigTemplate_,n_Integer]:=
	Tally/@DeleteDuplicates@Subsets[ConfigVertexListDegenerateOuterLoop@totalConfigTemplate,{n}];


AssignOuterLabels[totalConfigTemplate_,labels_]:=
	Flatten[
		AssignOuterTalliedLegs[totalConfigTemplate,labels,#]&/@GetOuterAssignments[totalConfigTemplate,Length@labels]
	];


AssignOuterTalliedLegs[totalConfigs_List,labelsOuter_List,talliedListOuter_List]:=AssignOuterTalliedLegs[#,labelsOuter,talliedListOuter]&/@totalConfigs;

AssignOuterTalliedLegs[totalConfig_Association,labelsOuter_List,{{position_,number_},rest___}] /; First@position === Key["Center"] := 
	AssignOuterTalliedLegs[#,Drop[labelsOuter,number],{rest}]&/@
	(ReplacePart[totalConfig,{
			{"Center","Top"}->#1,{"Center","Bottom"}->#2,{"Center","Free"}->totalConfig["Center","Free"]-number}]&@@@
		Table[{labelsOuter[[;;number]][[;;i]],labelsOuter[[;;number]][[i+1;;]]},{i,0,number}]);
		
AssignOuterTalliedLegs[totalConfig_Association,labelsOuter_List,{{position_,number_},rest___}] := 
	AssignOuterTalliedLegs[
		ReplacePart[totalConfig,{Append[position,"Outer"]->labelsOuter[[;;number]],Append[position,"Free"]->totalConfig[[Sequence@@position,"Free"]]-number}],
	Drop[labelsOuter,number],{rest}];
	
AssignOuterTalliedLegs[totalConfig_Association,{},{}]:=totalConfig;


(* ::Subsubsection::Closed:: *)
(*Inner Assignment*)


AssignInnerLoopLabels[totalConfigTemplate_,labels_,loop:"FirstLoop"|"SecondLoop"]:=
	Flatten[
		AssignInnerLoopTalliedLegs[totalConfigTemplate,labels,#,loop]&/@GetInnerLoopAssignments[totalConfigTemplate,Length@labels,loop]
	];


GetInnerLoopAssignments[totalConfigTemplate_,n_Integer,loop:"FirstLoop"|"SecondLoop"]:=
	Tally/@DeleteDuplicates@Subsets[ConfigVertexListDegenerateInnerLoop[loop]@totalConfigTemplate,{n}];


AssignInnerLoopTalliedLegs[totalConfigs_List,labelsInner_List,talliedListInner_List,loop:"FirstLoop"|"SecondLoop"]:=
	AssignInnerLoopTalliedLegs[#,labelsInner,talliedListInner,loop]&/@totalConfigs;

AssignInnerLoopTalliedLegs[totalConfig_Association,labelsInner_List,{{position_,number_},rest___},"FirstLoop"] /; First@position === Key["Center"] := 
	AssignInnerLoopTalliedLegs[
		ReplacePart[totalConfig,{Append[position,"FirstInner"]->labelsInner[[-number;;]],
			Append[position,"Free"]->totalConfig[[Sequence@@position,"Free"]]-number}],
	Drop[labelsInner,-number],{rest},"FirstLoop"];

AssignInnerLoopTalliedLegs[totalConfig_Association,labelsInner_List,{{position_,number_},rest___},"SecondLoop"] /; First@position === Key["Center"] := 
	AssignInnerLoopTalliedLegs[
		ReplacePart[totalConfig,{Append[position,"SecondInner"]->labelsInner[[-number;;]],Append[position,"Free"]->totalConfig[[Sequence@@position,"Free"]]-number}],
	Drop[labelsInner,-number],{rest},"SecondLoop"];

AssignInnerLoopTalliedLegs[totalConfig_Association,labelsInner_List,{{position_,number_},rest___},loop:"FirstLoop"|"SecondLoop"] := 
	AssignInnerLoopTalliedLegs[
		ReplacePart[totalConfig,{Append[position,"Inner"]->labelsInner[[-number;;]],Append[position,"Free"]->totalConfig[[Sequence@@position,"Free"]]-number}],
	Drop[labelsInner,-number],{rest},loop];
	
AssignInnerLoopTalliedLegs[totalConfig_Association,{},{},loop:"FirstLoop"|"SecondLoop"]:=totalConfig;


(* ::Subsubsection::Closed:: *)
(*Filling Template*)


ToggleTwist[totalConfig_]/; TopTwistedConfigQ[totalConfig]:=ReplacePart[{"Center","Twisted"}->"Bottom"]@totalConfig;
ToggleTwist[totalConfig_]/; BottomTwistedConfigQ[totalConfig]:=ReplacePart[{"Center","Twisted"}->"Top"]@totalConfig;
ToggleTwist[totalConfig_] := totalConfig;


FlipLoops[totalConfig_]:=FlipLoops[totalConfig]=ToggleTwist@ReplacePart[
	totalConfig,
	{
		{"FirstLoop"}->totalConfig["SecondLoop"],{"SecondLoop"}->totalConfig["FirstLoop"],
		{"Center","Top"}->totalConfig["Center","Bottom"],{"Center","Bottom"}->totalConfig["Center","Top"],
		{"Center","FirstInner"}->totalConfig["Center","SecondInner"],{"Center","SecondInner"}->totalConfig["Center","FirstInner"]
	}
];


TwistConfig[totalConfig_]:=
	ToggleTwist@(ReplacePart[
		#,
		{
			{"Center","Top"}->Reverse@#["Center","Top"],{"Center","Bottom"}->Reverse@#["Center","Bottom"],
			{"Center","FirstInner"}->Reverse@#["Center","SecondInner"],{"Center","SecondInner"}->Reverse@#["Center","FirstInner"]
		}
	]&@
	MapAt[
		ReplacePart[#,{{"Inner"}->Reverse@#["Outer"],{"Outer"}->Reverse@#["Inner"]}]&,
		{{"FirstLoop",All},{"SecondLoop",All}}
	]@totalConfig);


FlipInnerOuter[totalConfig_]/;TopTwistedConfigQ@totalConfig:=MapAt[ReplacePart[#,{{"Inner"}->#["Outer"],{"Outer"}->#["Inner"]}]&,totalConfig,{{"SecondLoop",All}}];
ReverseExteriorVertices[totalConfig_]/;TopTwistedConfigQ@totalConfig:=MapAt[Reverse,totalConfig,{{"SecondLoop"}}];

FlipInnerOuter[totalConfig_]/;BottomTwistedConfigQ@totalConfig:=MapAt[ReplacePart[#,{{"Inner"}->#["Outer"],{"Outer"}->#["Inner"]}]&,totalConfig,{{"FirstLoop",All}}];
ReverseExteriorVertices[totalConfig_]/;BottomTwistedConfigQ@totalConfig:=MapAt[Reverse,totalConfig,{{"FirstLoop"}}];


(* Rotates left the positions of Top, SecondInner, Bottom, and FirstInner *)
RotateLeftTwistedConfig[totalConfig_]:=RotateLeftTwistedConfig[totalConfig]=
	ReverseExteriorVertices@FlipInnerOuter@ReplacePart[
		#,
		{
			{"Center","Top"}->#["Center","SecondInner"],{"Center","SecondInner"}->#["Center","Bottom"],
			{"Center","Bottom"}->#["Center","FirstInner"],{"Center","FirstInner"}->#["Center","Top"],
			{"FirstLoop"}->#["SecondLoop"],{"SecondLoop"}->#["FirstLoop"]
		}
	]&@totalConfig;


(* Checks for equality under rotation of twisted central vertex *)
TwistRotatedDuplicatesQ[totalConfig1_,totalConfig2_] /; TwistedConfigQ[totalConfig1] && TwistedConfigQ[totalConfig2]:=
Or@@(
	SameQ[NormalizeAllLabels@totalConfig1,#]&/@
	NormalizeAllLabels/@Through[
		{
			RotateLeftTwistedConfig,FlipLoops@*RotateLeftTwistedConfig,
			RotateLeftTwistedConfig@*RotateLeftTwistedConfig,
			FlipLoops@*RotateLeftTwistedConfig@*RotateLeftTwistedConfig,
			RotateLeftTwistedConfig@*RotateLeftTwistedConfig@*RotateLeftTwistedConfig,
			FlipLoops@*RotateLeftTwistedConfig@*RotateLeftTwistedConfig@*RotateLeftTwistedConfig
		}[totalConfig2]
	]
);

TwistRotatedDuplicatesQ[___] := False;
RemoveTwistRotatedDuplicates[totalConfigs_List]:=
    DeleteDuplicates[totalConfigs,TwistRotatedDuplicatesQ];


CyclicDuplicateTotalConfigsQ[totalConfig1_,totalConfig2_]:= 
	NormalizeAllLabels@totalConfig1 === NormalizeAllLabels@FlipLoops@totalConfig2;
	
PlainDuplicateTotalConfigsQ[totalConfig1_,totalConfig2_]:= 
	NormalizeAllLabels@totalConfig1 === NormalizeAllLabels@totalConfig2;


TwistDuplicateTotalConfigsQ[totalConfig1_,totalConfig2_] /; TwistedConfigQ[totalConfig1] && TwistedConfigQ[totalConfig2] :=
	NormalizeAllLabels@totalConfig1 === TwistConfig@totalConfig2;
	
TwistDuplicateTotalConfigsQ[___] := False;


DuplicateTotalConfigsQ[totalConfig1_,totalConfig2_] := 
		PlainDuplicateTotalConfigsQ[totalConfig1,totalConfig2] ||
		CyclicDuplicateTotalConfigsQ[totalConfig1,totalConfig2] ||
		TwistDuplicateTotalConfigsQ[totalConfig1,totalConfig2];


RemoveDuplicateTotalConfigs[totalConfigs_List]:=DeleteDuplicates[totalConfigs,DuplicateTotalConfigsQ];


NormalizeAllLabels[totalConfig_] /; TwistedConfigQ[totalConfig] === False := NormalizeAllLabels[totalConfig] =
	NormalizeOuterOrdering[Sort@GetOuterMomenta[totalConfig]]@
	NormalizeInnerLoopOrdering[Sort@GetInnerLoopMomenta[totalConfig,"FirstLoop"],"FirstLoop"]@
	NormalizeInnerLoopOrdering[Sort@GetInnerLoopMomenta[totalConfig,"SecondLoop"],"SecondLoop"]@totalConfig;


NormalizeAllLabels[totalConfig_] /; TwistedConfigQ[totalConfig]:= NormalizeAllLabels[totalConfig] =
	NormalizeSingleTraceTwistedOrdering[Sort@GetSingleTraceTwistedMomenta[totalConfig]]@totalConfig;


GetInnerLoopMomenta[totalConfig_,"FirstLoop"]:=Flatten@
	Join[totalConfig[["FirstLoop",All,"Inner"]],totalConfig["Center","FirstInner"]];
GetInnerLoopMomenta[totalConfig_,"SecondLoop"]:=Flatten@
	Join[totalConfig[["SecondLoop",All,"Inner"]],totalConfig["Center","SecondInner"]];
	
GetInnerLoopMomenta[loop:"FirstLoop"|"SecondLoop"]:= Function[totalConfig,GetInnerLoopMomenta[totalConfig,loop]];


GetOuterMomenta[totalConfig_]:=
	Flatten@Join[
		totalConfig[["FirstLoop",All,"Outer"]],
		totalConfig["Center","Top"],
		totalConfig[["SecondLoop",All,"Outer"]],
		totalConfig["Center","Bottom"]
	];

GetSingleTraceTwistedMomenta[totalConfig_] /; TopTwistedConfigQ[totalConfig]:=
	Flatten@Join[
		totalConfig[["FirstLoop",All,"Outer"]],
		totalConfig["Center","SecondInner"],
		Flatten@Reverse@totalConfig[["SecondLoop",All,"Inner"]],
		totalConfig["Center","Top"],
		Flatten@Reverse@totalConfig[["FirstLoop",All,"Inner"]],
		totalConfig["Center","FirstInner"],
		totalConfig[["SecondLoop",All,"Outer"]],
		totalConfig["Center","Bottom"]
	];

GetSingleTraceTwistedMomenta[totalConfig_] /; BottomTwistedConfigQ[totalConfig]:=
	Flatten@Join[
		totalConfig[["FirstLoop",All,"Outer"]],
		totalConfig["Center","Top"],
		totalConfig[["SecondLoop",All,"Outer"]],
		totalConfig["Center","FirstInner"],
		Flatten@Reverse@totalConfig[["FirstLoop",All,"Inner"]],
		totalConfig["Center","Bottom"],
		Flatten@Reverse@totalConfig[["SecondLoop",All,"Inner"]],
		totalConfig["Center","SecondInner"]
	];


NormalizeOuterOrdering[totalConfig_,outerLabels_]:=totalConfig/.MapThread[Rule,{GetOuterMomenta[totalConfig],outerLabels}];
NormalizeOuterOrdering[outerLabels_]:=Function[totalConfig,NormalizeOuterOrdering[totalConfig,outerLabels]];

NormalizeInnerLoopOrdering[totalConfig_,innerLoopLabels_,loop:"FirstLoop"|"SecondLoop"]:=
	totalConfig/.MapThread[Rule,{GetInnerLoopMomenta[totalConfig,loop],innerLoopLabels}];
NormalizeInnerLoopOrdering[innerLoopLabels_,loop:"FirstLoop"|"SecondLoop"]:=
	Function[totalConfig,NormalizeInnerLoopOrdering[totalConfig,innerLoopLabels,loop]];
	
NormalizeSingleTraceTwistedOrdering[totalConfig_,labels_]:=totalConfig/.MapThread[Rule,{GetSingleTraceTwistedMomenta[totalConfig],labels}];
NormalizeSingleTraceTwistedOrdering[labels_]:=Function[totalConfig,NormalizeSingleTraceTwistedOrdering[totalConfig,labels]];


DropVertexIndexFreeKeys[totalConfigFilledTemplate_]:=
	MapAt[KeyDrop[{"VertexIndex","Free"}],totalConfigFilledTemplate,{{"FirstLoop",All},{"SecondLoop",All},{"Center"}}];


FillConfigurationTemplateInnerLoop[totalConfigTemplateSemiFilled_,_]:=totalConfigTemplateSemiFilled;

FillConfigurationTemplateInnerLoop[totalConfigTemplateSemiFilled_,"FirstLoop",innerLabels1_,innerLabels2___]:=
	FillConfigurationTemplateInnerLoop[#,"SecondLoop",innerLabels2]&/@
		AssignInnerLoopLabels[totalConfigTemplateSemiFilled,innerLabels1,"FirstLoop"];
		
FillConfigurationTemplateInnerLoop[totalConfigTemplateSemiFilled_,"SecondLoop",innerLabels1_,innerLabels2___]:=
	FillConfigurationTemplateInnerLoop[#,"FirstLoop",innerLabels2]&/@
		AssignInnerLoopLabels[totalConfigTemplateSemiFilled,innerLabels1,"SecondLoop"];


(* We only need to check for cyclic duplicates within the set of configs from a single template, since 
	we have already removed cyclic duplicate templates (cf. RemoveSymmetricDuplicates).
	Also, we are generating the two assignments of the inner traces at the same
	time, which explains the two FillConfigurationTemplateInnerLoop appearing *)
FillConfigurationTemplate[totalConfigTemplate_,outerLabels_List,innerLabels___]:=
	RemoveDuplicateTotalConfigs[
	NormalizeOuterOrdering[outerLabels]/@DropVertexIndexFreeKeys/@Flatten[
		Join[
			FillConfigurationTemplateInnerLoop[#,"FirstLoop",innerLabels],
			FillConfigurationTemplateInnerLoop[#,"SecondLoop",innerLabels]
		]&/@AssignOuterLabels[totalConfigTemplate,outerLabels]
	]
];


(* ::Subsection::Closed:: *)
(*Putting Templates and Filling Together*)


GenerateTotalConfigs[labelsOuter_,labelsInner___]/; DuplicateFreeQ@Join[labelsOuter,labelsInner]:=Module[{templates},
	templates = CreateTotalConfigTemplates[Length@Join[labelsOuter,labelsInner]];
	Flatten[FillConfigurationTemplate[#,labelsOuter,labelsInner]&/@templates]
];

GenerateTotalConfigs[labelsOuter_,labelsInner___] /; Message[GenerateTotalConfigs::duplc,{labelsOuter,labelsInner}] := Null;


SymmetricTotalConfig[totalConfig_]:=CyclicDuplicateTotalConfigsQ[FuseCutBubbles@totalConfig,FuseCutBubbles@totalConfig];


GenerateComputationTotalConfigs[labelsOuter_,labelsInner___] /; DuplicateFreeQ@Join[labelsOuter,labelsInner]:=
	{If[SymmetricTotalConfig@#,1/2,1],#}&/@
		CutBubbles@GenerateTotalConfigs[labelsOuter,labelsInner];

GenerateComputationTotalConfigs[labelsOuter_,labelsInner___] /; Message[GenerateComputationTotalConfigs::duplc,{labelsOuter,labelsInner}] := Null;


(* ::Subsection:: *)
(*Cutting of Bubbles*)


CutBubble[totalConfig_]:=
    Flatten[CutFirstLoopBubble/@CutSecondLoopBubble@totalConfig]/.{a___,$cutBubble,$cutBubble,b___}:>{a,$cutBubble,b};
CutBubbles[totalConfigs_List]:=Flatten[CutBubble/@totalConfigs];


CutFirstLoopBubble[totalConfig_] /; MatchQ[totalConfig,firstLoopBubblePattern] := Join@@Through[{List,CutOuterFirstLoop,CutFirstLoopTopCenter,CutFirstLoopBottomCenter}@totalConfig];
CutSecondLoopBubble[totalConfig_] /; MatchQ[totalConfig,secondLoopBubblePattern] := Join@@Through[{List,CutOuterSecondLoop,CutSecondLoopTopCenter,CutSecondLoopBottomCenter}@totalConfig];

CutFirstLoopBubble[totalConfig_] :={totalConfig};
CutSecondLoopBubble[totalConfig_] :={totalConfig};


(* ::Subsubsection:: *)
(*Cut Outer Loops*)


CutOuterFirstLoop[totalConfig_] /; MatchQ[totalConfig,firstLoopBubblePattern] := Module[{cutVertices},
	cutVertices = CutOuterVertex/@totalConfig["FirstLoop"];
	Flatten@Table[
		ReplacePart[totalConfig,{"FirstLoop"}->ReplacePart[totalConfig["FirstLoop"],{vertPos->Sequence@@vertexCut}]],
		{vertPos,Length@totalConfig["FirstLoop"]},
		{vertexCut,cutVertices[[vertPos]]}
	]
];


CutOuterSecondLoop[totalConfig_] /; MatchQ[totalConfig,secondLoopBubblePattern] := Module[{cutVertices},
	cutVertices = CutOuterVertex/@totalConfig["SecondLoop"];
	Flatten@Table[
		ReplacePart[totalConfig,{"SecondLoop"}->ReplacePart[totalConfig["SecondLoop"],{vertPos->Sequence@@vertexCut}]],
		{vertPos,Length@totalConfig["SecondLoop"]},
		{vertexCut,cutVertices[[vertPos]]}
	]
];


CutOuterVertex[vertex_Association] := Module[{cuts},
	cuts = Replace[
		Transpose/@Tuples[{Table[{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{___,{{},_},{_,{}},___}->Nothing,{___,{_,{}},{{},_},___}->Nothing},2
	]&[vertex["Outer"],vertex["Inner"]];

	{<|"Outer"->#[[2,1]],"Inner"->#[[1,2]],"Properties"->Append[vertex["Properties"],$cutBubble]|>,
		<|"Outer"->#[[1,1]],"Inner"->#[[2,2]],"Properties"->Append[vertex["Properties"],$cutBubble]|>}&/@cuts
];


(* ::Subsubsection::Closed:: *)
(*Second Loop Center*)


CutSecondLoopTopCenter[totalConfig_] := Module[{cutVertices},
	cutVertices = GetSecondLoopCenterTopCut@totalConfig["Center"];
	Table[
		MapAt[Prepend[Last@vertexCut],ReplacePart[totalConfig,{"Center"}->First@vertexCut],{"SecondLoop"}],
		{vertexCut,cutVertices}
	]
];


GetSecondLoopCenterTopCut[vertex_Association] /; vertex["Twisted"]=!="Top" :=Module[{cuts,twisted},
	twisted=If[KeyExistsQ[vertex,"Twisted"],vertex["Twisted"],Null];
	cuts = Replace[
		Transpose/@Tuples[{Table[{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[Reverse@{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Top"],vertex["SecondInner"]];
	{<|"Top"->#[[2,1]],"Bottom"->vertex["Bottom"],"FirstInner"->vertex["FirstInner"],
			"SecondInner"->#[[2,2]],"Twisted"->twisted,"Properties"->Append[vertex["Properties"],$cutBubble]|>,
		<|"Outer"->#[[1,1]],"Inner"->#[[1,2]],"Properties"->{$cutBubble}|>}&/@cuts
];


GetSecondLoopCenterTopCut[vertex_Association] /; vertex["Twisted"]==="Top" :=Module[{cuts},
	cuts = Replace[
		Transpose/@Tuples[{Table[Reverse@{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Top"],vertex["FirstInner"]];
	{<|"Top"->#[[2,1]],"Bottom"->vertex["Bottom"],"FirstInner"->#[[2,2]],
			"SecondInner"->vertex["SecondInner"],"Twisted"->vertex["Twisted"],"Properties"->Append[vertex["Properties"],$cutBubble]|>,
		<|"Outer"->#[[1,2]],"Inner"->#[[1,1]],"Properties"->{$cutBubble}|>}&/@cuts
];





CutSecondLoopBottomCenter[totalConfig_] := Module[{cutVertices},
	cutVertices = GetSecondLoopCenterBottomCut@totalConfig["Center"];
	Table[
		MapAt[Append[Last@vertexCut],ReplacePart[totalConfig,{"Center"}->First@vertexCut],{"SecondLoop"}],
		{vertexCut,cutVertices}
	]
];


GetSecondLoopCenterBottomCut[vertex_Association] /; vertex["Twisted"]=!="Bottom" :=Module[{cuts,twisted},
	twisted=If[KeyExistsQ[vertex,"Twisted"],vertex["Twisted"],Null];
	cuts = Replace[
		Transpose/@Tuples[{Table[Reverse@{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{___,{{},{}}}->Nothing,2
	]&[vertex["SecondInner"],vertex["Bottom"]];
	{<|"Top"->vertex["Top"],"Bottom"->#[[1,2]],"FirstInner"->vertex["FirstInner"],
		"SecondInner"->#[[1,1]],"Twisted"->twisted,"Properties"->Append[vertex["Properties"],$cutBubble]|>,
		<|"Outer"->#[[2,2]],"Inner"->#[[2,1]],"Properties"->{$cutBubble}|>}&/@cuts
];


GetSecondLoopCenterBottomCut[vertex_Association] /; vertex["Twisted"]==="Bottom" :=Module[{cuts},
	cuts = Replace[
		Transpose/@Tuples[{Table[{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[Reverse@{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{___,{{},{}}}->Nothing,2
	]&[vertex["FirstInner"],vertex["Bottom"]];
	{<|"Top"->vertex["Top"],"Bottom"->#[[1,2]],"FirstInner"->#[[1,1]],
		"SecondInner"->vertex["SecondInner"],"Twisted"->vertex["Twisted"],"Properties"->Append[vertex["Properties"],$cutBubble]|>,
		<|"Outer"->#[[2,1]],"Inner"->#[[2,2]],"Properties"->{$cutBubble}|>}&/@cuts
];


CutSecondLoopBottomCenter@<|
			"FirstLoop"->{<|"Outer"->{pFirst},"Inner"->{},"Properties"->{}|>},
			"SecondLoop"->{<|"Outer"->{pSecond},"Inner"->{},"Properties"->{}|>},
			"Center"-><|"Top"->{p1,p2},"Bottom"->{p5,p6},"FirstInner"->{p7,p8},"SecondInner"->{p3,p4},"Properties"->{}|>
		|>;


(* ::Subsubsection::Closed:: *)
(*First Loop Center*)


CutFirstLoopTopCenter[totalConfig_] := Module[{cutVertices},
	cutVertices = GetFirstLoopCenterTopCut@totalConfig["Center"];
	Table[
		MapAt[Append[First@vertexCut],ReplacePart[totalConfig,{"Center"}->Last@vertexCut],{"FirstLoop"}],
		{vertexCut,cutVertices}
	]
];


GetFirstLoopCenterTopCut[vertex_Association]  /; vertex["Twisted"]=!="Top" :=Module[{cuts,twisted},
	twisted=If[KeyExistsQ[vertex,"Twisted"],vertex["Twisted"],Null];
	cuts = Replace[
		Transpose/@Tuples[{Table[Reverse@{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Top"],vertex["FirstInner"]];
	{<|"Outer"->#[[1,1]],"Inner"->#[[1,2]],"Properties"->{$cutBubble}|>,
	<|"Top"->#[[2,1]],"Bottom"->vertex["Bottom"],"FirstInner"->#[[2,2]],
			"SecondInner"->vertex["SecondInner"],"Twisted"->twisted,"Properties"->Append[vertex["Properties"],$cutBubble]|>}&/@cuts
];


GetFirstLoopCenterTopCut[vertex_Association] /; vertex["Twisted"]==="Top" :=Module[{cuts},
	cuts = Replace[
		Transpose/@Tuples[{Table[{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[Reverse@{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Top"],vertex["SecondInner"]];
	{<|"Outer"->#[[1,2]],"Inner"->#[[1,1]],"Properties"->{$cutBubble}|>,
	<|"Top"->#[[2,1]],"Bottom"->vertex["Bottom"],"FirstInner"->vertex["FirstInner"],
			"SecondInner"->#[[2,2]],"Twisted"->vertex["Twisted"],"Properties"->Append[vertex["Properties"],$cutBubble]|>}&/@cuts
];


CutFirstLoopBottomCenter[totalConfig_] := Module[{cutVertices},
	cutVertices = GetFirstLoopCenterBottomCut@totalConfig["Center"];
	Table[
		MapAt[Prepend[First@vertexCut],ReplacePart[totalConfig,{"Center"}->Last@vertexCut],{"FirstLoop"}],
		{vertexCut,cutVertices}
	]
];


GetFirstLoopCenterBottomCut[vertex_Association] /; vertex["Twisted"]=!="Bottom" :=Module[{cuts,twisted},
	twisted=If[KeyExistsQ[vertex,"Twisted"],vertex["Twisted"],Null];
	cuts = Replace[
		Transpose/@Tuples[{Table[{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[Reverse@{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Bottom"],vertex["FirstInner"]];
	{<|"Outer"->#[[1,1]],"Inner"->#[[1,2]],"Properties"->{$cutBubble}|>,
		<|"Top"->vertex["Top"],"Bottom"->#[[2,1]],"FirstInner"->#[[2,2]],
			"SecondInner"->vertex["SecondInner"],"Twisted"->twisted,"Properties"->Append[vertex["Properties"],$cutBubble]|>}&/@cuts
];


GetFirstLoopCenterBottomCut[vertex_Association]/; vertex["Twisted"]==="Bottom" :=Module[{cuts},
	cuts = Replace[
		Transpose/@Tuples[{Table[Reverse@{#1[[i+1;;]],#1[[;;i]]},{i,0,Length@#1}],Table[{#2[[i+1;;]],#2[[;;i]]},{i,0,Length@#2}]}],
		{{{},{}},___}->Nothing,2
	]&[vertex["Bottom"],vertex["SecondInner"]];
	{<|"Outer"->#[[1,2]],"Inner"->#[[1,1]],"Properties"->{$cutBubble}|>,
		<|"Top"->vertex["Top"],"Bottom"->#[[2,1]],"FirstInner"->vertex["FirstInner"],
			"SecondInner"->#[[2,2]],"Twisted"->vertex["Twisted"],"Properties"->Append[vertex["Properties"],$cutBubble]|>}&/@cuts
];


(* ::Subsection::Closed:: *)
(*CutBubble Helpers*)


FuseCutBubbles[totalConfig_]:=FuseFirstLoopCutBubble@FuseSecondLoopCutBubble@totalConfig;


CutTopFirstLoopQ[totalConfig_]:= MemberQ[totalConfig["FirstLoop"][[-1]]["Properties"], $cutBubble] && Count[totalConfig["FirstLoop"],$cutBubble,Infinity] == 1;
CutBottomFirstLoopQ[totalConfig_]:= MemberQ[totalConfig["FirstLoop"][[1]]["Properties"], $cutBubble] && Count[totalConfig["FirstLoop"],$cutBubble,Infinity] == 1;


CutTopSecondLoopQ[totalConfig_]:= MemberQ[totalConfig["SecondLoop"][[1]]["Properties"], $cutBubble] && Count[totalConfig["SecondLoop"],$cutBubble,Infinity] == 1;
CutBottomSecondLoopQ[totalConfig_]:= MemberQ[totalConfig["SecondLoop"][[-1]]["Properties"], $cutBubble] && Count[totalConfig["SecondLoop"],$cutBubble,Infinity] == 1;


fuseRestore$cutBubbleLabels[totalConfig_] := 
Module[{temp = totalConfig},

	If[
		FreeQ[totalConfig["Center"]["Properties"],$cutBubble] &&
		(CutTopFirstLoopQ[totalConfig] || CutBottomFirstLoopQ[totalConfig] || CutBottomSecondLoopQ[totalConfig] || CutTopSecondLoopQ[totalConfig] ),
		AppendTo[temp["Center"]["Properties"],$cutBubble];
	];
	temp
];


FuseSecondLoopCutBubble[totalConfig_](* /; MatchQ[totalConfig,secondLoopBubbleTrianglePattern]*) := Module[{temp = totalConfig},
	temp = MapAt[
		Replace[#,
			{a___,<|"Outer"->out1_,"Inner"->in1_,"Properties"->prop1_|>,<|"Outer"->out2_,"Inner"->in2_,"Properties"->prop2_|>,d___} /; MemberQ[prop1,$cutBubble]&&MemberQ[prop2,$cutBubble] :>
				{a,<|"Outer"->Join[out1,out2],"Inner"->Join[in1,in2],"Properties"->Join[prop1,prop2]/.$cutBubble->Nothing|>,d}]&,
	temp,{"SecondLoop"}];
	(* fusing top *)
	If[CutTopSecondLoopQ[totalConfig],
		temp["SecondLoop"] = Rest@temp["SecondLoop"];
		If[TopTwistedConfigQ@totalConfig,
			temp["Center"]["Top"] = Join[totalConfig["SecondLoop"][[1]]["Inner"],totalConfig["Center"]["Top"]];
			temp["Center"]["FirstInner"] = Join[totalConfig["Center"]["FirstInner"],totalConfig["SecondLoop"][[1]]["Outer"]];
		,
			temp["Center"]["Top"] = Join[totalConfig["Center"]["Top"],totalConfig["SecondLoop"][[1]]["Outer"]];
			temp["Center"]["SecondInner"] = Join[totalConfig["SecondLoop"][[1]]["Inner"],totalConfig["Center"]["SecondInner"]];
		];
		temp["Center"]["Properties"] = temp["Center"]["Properties"]/.$cutBubble->Nothing;
	];
	
	(* fusing bottom *)
	If[CutBottomSecondLoopQ[totalConfig],
		temp["SecondLoop"] = Most@temp["SecondLoop"];
		If[BottomTwistedConfigQ@totalConfig,
			temp["Center"]["Bottom"] = Join[totalConfig["Center"]["Bottom"],totalConfig["SecondLoop"][[-1]]["Inner"]];
			temp["Center"]["FirstInner"] = Join[totalConfig["SecondLoop"][[-1]]["Outer"],totalConfig["Center"]["FirstInner"]];
		,
			temp["Center"]["Bottom"] = Join[totalConfig["SecondLoop"][[-1]]["Outer"],totalConfig["Center"]["Bottom"]];
			temp["Center"]["SecondInner"] = Join[totalConfig["Center"]["SecondInner"],totalConfig["SecondLoop"][[-1]]["Inner"]];
		];
		temp["Center"]["Properties"] = temp["Center"]["Properties"]/.$cutBubble->Nothing;
	];
	fuseRestore$cutBubbleLabels[temp]
];


FuseFirstLoopCutBubble[totalConfig_] (*/; MatchQ[totalConfig,firstLoopBubbleTrianglePattern]*) := Module[{temp = totalConfig},
	temp = MapAt[
		Replace[#,
			{a___,<|"Outer"->out1_,"Inner"->in1_,"Properties"->prop1_|>,<|"Outer"->out2_,"Inner"->in2_,"Properties"->prop2_|>,d___} /; MemberQ[prop1,$cutBubble]&&MemberQ[prop2,$cutBubble] :>
				{a,<|"Outer"->Join[out1,out2],"Inner"->Join[in1,in2],"Properties"->Join[prop1,prop2]/.$cutBubble->Nothing|>,d}]&,
	temp,{"FirstLoop"}];
	
	(* fusing top *)
	If[CutTopFirstLoopQ[totalConfig],
		temp["FirstLoop"] = Most@temp["FirstLoop"];
		If[TopTwistedConfigQ@totalConfig,
			temp["Center"]["Top"] = Join[totalConfig["Center"]["Top"],totalConfig["FirstLoop"][[-1]]["Inner"]];
			temp["Center"]["SecondInner"] = Join[totalConfig["FirstLoop"][[-1]]["Outer"],totalConfig["Center"]["SecondInner"]];
		,
			temp["Center"]["Top"] = Join[totalConfig["FirstLoop"][[-1]]["Outer"],totalConfig["Center"]["Top"]];
			temp["Center"]["FirstInner"] = Join[totalConfig["Center"]["FirstInner"],totalConfig["FirstLoop"][[-1]]["Inner"]];
		];
		temp["Center"]["Properties"] = temp["Center"]["Properties"]/.$cutBubble->Nothing;
	];
	
	(* fusing bottom *)
	If[CutBottomFirstLoopQ[totalConfig],
		temp["FirstLoop"] = Rest@temp["FirstLoop"];
		If[BottomTwistedConfigQ@totalConfig,
			temp["Center"]["Bottom"] = Join[totalConfig["FirstLoop"][[1]]["Inner"],totalConfig["Center"]["Bottom"]];
			temp["Center"]["SecondInner"] = Join[totalConfig["Center"]["SecondInner"],totalConfig["FirstLoop"][[1]]["Outer"]];
		,
			temp["Center"]["Bottom"] = Join[totalConfig["Center"]["Bottom"],totalConfig["FirstLoop"][[1]]["Outer"]];
			temp["Center"]["FirstInner"] = Join[totalConfig["FirstLoop"][[1]]["Inner"],totalConfig["Center"]["FirstInner"]];
		];
		temp["Center"]["Properties"] = temp["Center"]["Properties"]/.$cutBubble->Nothing;
	];
	
	fuseRestore$cutBubbleLabels[temp]
];


(* ::Subsection::Closed:: *)
(*Subleading Single Trace*)


CreateTwistedCentralVertices[totalConfigTemplates_]:=Join@@(CreateTwistedCentralVertex/@totalConfigTemplates);

CreateTwistedCentralVertex[totalConfigTemplate_Association]:=
	Through[{ReplacePart[{"Center","Twisted"}->"Top"],ReplacePart[{"Center","Twisted"}->"Bottom"]}[totalConfigTemplate]];


DistributeTemplateFreeToFreeInnerFreeOuter[totalConfigTemplates_List]:=Join@@(DistributeTemplateFreeToFreeInnerFreeOuter/@totalConfigTemplates);
DistributeTemplateFreeToFreeInnerFreeOuter[totalConfigTemplate_Association]:=
	DistributeLoopFreeToFreeInnerFreeOuter["First"]@
	DistributeLoopFreeToFreeInnerFreeOuter["Second"]@
	DistributeCenterFreeToFreeTopBottomLeftRight@totalConfigTemplate;


DistributeLoopFreeToFreeInnerFreeOuter[loop:"First"|"Second"]:=DistributeLoopFreeToFreeInnerFreeOuter[#,loop]&;
DistributeLoopFreeToFreeInnerFreeOuter[totalConfigTemplates_List,loop:"First"|"Second"]:=
	Join@@(DistributeLoopFreeToFreeInnerFreeOuter[loop]/@totalConfigTemplates);

DistributeLoopFreeToFreeInnerFreeOuter[totalConfigTemplate_Association,loop:"First"]:=
	ReplacePart[totalConfigTemplate,"FirstLoop"->#]&/@Tuples[DistributeVertexFreeToFreeInnerFreeOuter/@totalConfigTemplate["FirstLoop"]];
DistributeLoopFreeToFreeInnerFreeOuter[totalConfigTemplate_Association,loop:"Second"]:=
	ReplacePart[totalConfigTemplate,"SecondLoop"->#]&/@Tuples[DistributeVertexFreeToFreeInnerFreeOuter/@totalConfigTemplate["SecondLoop"]];


DistributeVertexFreeToFreeInnerFreeOuter[vertex_]:=
	Join[KeyDrop[vertex,"Free"],<|"FreeInner"->#,"FreeOuter"->vertex["Free"]-#|>]&/@ Range[0,vertex["Free"]];


DistributeCenterFreeToFreeTopBottomLeftRight[totalConfigTemplates_List]:=
	Join@@(DistributeCenterFreeToFreeTopBottomLeftRight/@totalConfigTemplates);
	
DistributeCenterFreeToFreeTopBottomLeftRight[totalConfigTemplate_Association]:=
	ReplacePart[totalConfigTemplate,"Center"->#]&/@CreateCenterVerticesFreeToFreeTopBottomLeftRight[totalConfigTemplate["Center"]];


CreateCenterVerticesFreeToFreeTopBottomLeftRight[center_Association]:=
	Join[KeyDrop[center,"Free"],<|"FreeTop"->#[[1]],"FreeBottom"->#[[2]],"FreeFirstInner"->#[[3]],"FreeSecondInner"->#[[4]]|>]&/@ 
		(Join@@(Permutations/@(PadRight[#,4]&/@IntegerPartitions[center["Free"],4])));





DistributeSubleadingSingleTrace[totalConfigTemplates_List,momenta_List]:= DistributeSubleadingSingleTrace[#,momenta]&/@totalConfigTemplates;

DistributeSubleadingSingleTrace[totalConfigTemplate_Association,momenta_List] /; totalConfigTemplate["Center","Twisted"]==="Bottom" :=
	DropVertexIndexFreeKeys@First[
		DistributeSingleTraceCenter["SecondInner"]@@
		DistributeSingleTraceLoop["SecondLoop","Inner"]@@
		DistributeSingleTraceCenter["Bottom"]@@
		DistributeSingleTraceLoop["FirstLoop","Inner"]@@
		DistributeSingleTraceCenter["FirstInner"]@@
		DistributeSingleTraceLoop["SecondLoop","Outer"]@@
		DistributeSingleTraceCenter["Top"]@@
		DistributeSingleTraceLoop["FirstLoop","Outer"]@@
			{totalConfigTemplate,momenta}
	];
	
DistributeSubleadingSingleTrace[totalConfigTemplate_Association,momenta_List] /; totalConfigTemplate["Center","Twisted"]==="Top" :=
	DropVertexIndexFreeKeys@First[
		DistributeSingleTraceCenter["Bottom"]@@
		DistributeSingleTraceLoop["SecondLoop","Outer"]@@
		DistributeSingleTraceCenter["FirstInner"]@@
		DistributeSingleTraceLoop["FirstLoop","Inner"]@@
		DistributeSingleTraceCenter["Top"]@@
		DistributeSingleTraceLoop["SecondLoop","Inner"]@@
		DistributeSingleTraceCenter["SecondInner"]@@
		DistributeSingleTraceLoop["FirstLoop","Outer"]@@
			{totalConfigTemplate,momenta}
	];


ClearAll[DistributeSingleTraceCenter];
DistributeSingleTraceCenter[pos:"Top"|"Bottom"|"SecondInner"|"FirstInner"]:=
Function[{totalConfigTemplate,momenta},
	Module[{freePos = "Free" <> pos, nFree},
		nFree=totalConfigTemplate["Center",freePos];
		{
			MapAt[KeyDrop[freePos],{"Center"}]@ReplacePart[
				totalConfigTemplate,
				{"Center",pos}->momenta[[;;nFree]]
			],
			momenta[[nFree+1;;]]
		}
	]
];


DistributeSingleTraceLoop[loop:"FirstLoop"|"SecondLoop","Outer"]:=Function[{totalConfigTemplate,momenta},
	Fold[AttachMomentaLoop[#1,#2,loop,"Outer"]&,{totalConfigTemplate,momenta},Range@Length@totalConfigTemplate[loop]]
];

DistributeSingleTraceLoop[loop:"FirstLoop"|"SecondLoop","Inner"]:=Function[{totalConfigTemplate,momenta},
	Fold[AttachMomentaLoop[#1,#2,loop,"Inner"]&,{totalConfigTemplate,momenta},Reverse@Range@Length@totalConfigTemplate[loop]]
];


AttachMomentaLoop[{template_Association,momenta_},pos_,loop:"FirstLoop"|"SecondLoop","Outer"]:=Module[{nFree=template[[loop,pos,"FreeOuter"]]},
	{
		MapAt[KeyDrop["FreeOuter"],{loop,pos}]@ReplacePart[
			template,
			{loop,pos,"Outer"}->momenta[[;;nFree]]
		],
		momenta[[nFree+1;;]]
	}
];

AttachMomentaLoop[{template_Association,momenta_},pos_,loop:"FirstLoop"|"SecondLoop","Inner"]:=Module[{nFree=template[[loop,pos,"FreeInner"]]},
	{
		MapAt[KeyDrop["FreeInner"],{loop,pos}]@ReplacePart[
			template,
			{loop,pos,"Inner"}->momenta[[;;nFree]]
		],
		momenta[[nFree+1;;]]
	}
];

FilterDoubleTriangleTemplates[True] :=
    Cases[secondLoopTrianglePattern]@*Cases[firstLoopTrianglePattern];
FilterDoubleTriangleTemplates[False] := Identity;


Options[GenerateSubLeadingSingleTraceConfigurations] = {TrianglesOnly -> False};

GenerateSubLeadingSingleTraceConfigurations[momenta_List,opts:OptionsPattern[]]:=
	RemoveTwistRotatedDuplicates @ Flatten[
		RemoveDuplicateTotalConfigs /@ (
			DistributeSubleadingSingleTrace[#,momenta] & /@
				DistributeTemplateFreeToFreeInnerFreeOuter /@
        	CreateTwistedCentralVertex /@
						FilterDoubleTriangleTemplates[OptionValue[TrianglesOnly]] @
							CreateTotalConfigTemplates[Length@momenta]
		)
	];



TwistedSymmetryFactor[totalConfig_]:=1/(1+Count[#,True])&@(And[!SameQ[totalConfig,#],SameQ[NormalizeAllLabels@totalConfig,NormalizeAllLabels@#]]&/@Through[
		{
			RotateLeftTwistedConfig,
			RotateLeftTwistedConfig@*RotateLeftTwistedConfig,
			RotateLeftTwistedConfig@*RotateLeftTwistedConfig@*RotateLeftTwistedConfig
		}[totalConfig]
	]);

Options[GenerateComputeSubLeadingSingleTraceConfigurations] = {TrianglesOnly -> False};
GenerateComputeSubLeadingSingleTraceConfigurations[momenta_List,opts:OptionsPattern[]]:=
	{TwistedSymmetryFactor[FuseCutBubbles@#],#}&/@
     (CutBubbles@GenerateSubLeadingSingleTraceConfigurations[
			 momenta,
			 PassRules[
				 GenerateComputeSubLeadingSingleTraceConfigurations,
				 GenerateSubLeadingSingleTraceConfigurations,opts]
		 		]
		 );


(* ::Section::Closed:: *)
(*PostProcessing*)


NestedReplaceNumSymbs = Function[expr,
	expr//.NumSymb->NumSymbN
];




(* ::Section::Closed:: *)
(*Symmetry Factors*)


ComputeSymmetryFactor[markedTotalConfig_]:=Module[{unmarkedTotalConfig,symmetryFactor},
	unmarkedTotalConfig = UnmarkTotalConfig@markedTotalConfig;
	
	symmetryFactor = Which[
		CompletelySymmetricTotalConfigQ[unmarkedTotalConfig],
			1/2
	,
		PreferredTotalConfigQ[unmarkedTotalConfig],
			1
	,
		True,
			0
	];
	
	If[!MatchQ[markedTotalConfig,$FullConfigPattern], symmetryFactor = 0];
	
	(*If[$FrontEnd === Null,
		Print["\t\t\tSymmetry Factor of ",markedTotalConfig ," is ", symmetryFactor];,
		Print["\t\t\tSymmetry Factor of ",DrawTotalConfig[markedTotalConfig]," is ", symmetryFactor];
	];*)
	
	(*If[Head[globalconfigList] =!= List, globalconfigList = {}];
	AppendTo[globalconfigList, {markedTotalConfig,symmetryFactor}];*)
	
	(*res *= CreateConstant[markedTotalConfig];*)

	symmetryFactor
];


(*ComputeSymmetryFactor[_]:=1*)


CreateConfigName[totalConfig_] := StringReplace[ToString[totalConfig],{"{"->"O","}"->"C",","->"K"," "->"","$cutBubble"->"cb"}];


CreateConstant[markedTotalConfig_] := Symbol["C"<>StringReplace[ToString[markedTotalConfig],{"{"->"o","}"->"c",","->"k"," "->""}]];


(*CreateConstant[{{{{1},{2},{3}},{4}},{{{5,$cutBubble},{6},{7}},{8}}}]*)


CompletelySymmetricTotalConfigQ[totalConfig_]:=SameQ@@MapAt[Length,totalConfig,{{All,1,All},{All,-1}}];


UnmarkTotalConfig[markedTotalConfig_]:=Module[{unmarkStepTotalConfig = markedTotalConfig},
	(* unmarking first loop *)
	unmarkStepTotalConfig = Which[
		Count[unmarkStepTotalConfig[[1,1]],$cutBubble,Infinity] == 2,
			MapAt[ReplaceAll[#,{a___,{b__,$cutBubble},{c__,$cutBubble},d___} :> {a,{b,c},d}]&,unmarkStepTotalConfig,1],
		MemberQ[markedTotalConfig[[1,1,1]],$cutBubble],
			MapAt[Join[#,UnmarkConfig[unmarkStepTotalConfig[[1,1,1]]]]&,Delete[unmarkStepTotalConfig,{1,1,1}],{-1,-1}],
		MemberQ[markedTotalConfig[[1,1,-1]],$cutBubble],
			MapAt[Join[UnmarkConfig[unmarkStepTotalConfig[[1,1,-1]]],#]&,Delete[unmarkStepTotalConfig,{1,1,-1}],{1,-1}],
		True,
			unmarkStepTotalConfig
	];
	
	(* unmarking second loop *)
	unmarkStepTotalConfig = Which[
		Count[unmarkStepTotalConfig[[2,1]],$cutBubble,Infinity] == 2,
			MapAt[ReplaceAll[#,{a___,{b__,$cutBubble},{c__,$cutBubble},d___} :> {a,{b,c},d}]&,unmarkStepTotalConfig,2],
		MemberQ[markedTotalConfig[[2,1,1]],$cutBubble],
			MapAt[Join[#,UnmarkConfig[unmarkStepTotalConfig[[2,1,1]]]]&,Delete[unmarkStepTotalConfig,{2,1,1}],{1,-1}],
		MemberQ[markedTotalConfig[[2,1,-1]],$cutBubble],
			MapAt[Join[UnmarkConfig[unmarkStepTotalConfig[[2,1,-1]]],#]&,Delete[unmarkStepTotalConfig,{2,1,-1}],{-1,-1}],
		True,
			unmarkStepTotalConfig
	];
	unmarkStepTotalConfig
];


PreferredTotalConfigQ[totalConfig_] :=
	Or[
		FirstCornerHeavyQ@@totalConfig,
		And[
			SymmetricBaseQ@totalConfig,
			GreaterEqual@@Length/@totalConfig[[All,2]]
		],
		And[
			!SymmetricBaseQ@totalConfig,
			CornerBalancedQ@@totalConfig,
			BaseFirstLoopLegHeavyQ@totalConfig
		]
	];


BaseFirstLoopLegHeavyQ[totalConfig_] := 
Module[{firstLoopConfigFlat, secondLoopConfigFlat, partialSumsFirstLoop, partialSumsSecondLoop},
	{firstLoopConfigFlat,secondLoopConfigFlat} = totalConfig[[All,1]];
	partialSumsFirstLoop = Length/@FoldList[Join,firstLoopConfigFlat];
	partialSumsSecondLoop = Length/@FoldList[Join,secondLoopConfigFlat];
	Greater@@FromDigits@*Reverse/@{partialSumsFirstLoop,partialSumsSecondLoop}
];


SymmetricBaseQ[totalConfig_] := CompletelySymmetricTotalConfigQ@ReplacePart[totalConfig,{{1,2}->{},{2,2}->{}}];


FirstCornerHeavyQ[firstLoopConfig_,secondLoopConfig_] := Length@firstLoopConfig[[1]] > Length@secondLoopConfig[[1]];
CornerBalancedQ[firstLoopConfig_,secondLoopConfig_] := Length@firstLoopConfig[[1]] == Length@secondLoopConfig[[1]];


(* ::Section::Closed:: *)
(*Integrals*)


Options[ConfigIntegral]={Gravity->False};


ConfigIntegral[totalConfig_,opts:OptionsPattern[]] := 
	{
		FirstLoopConfigIntegrals[totalConfig,PassRules[ConfigIntegral,FirstLoopConfigIntegrals,opts]],
		SecondLoopConfigIntegrals[totalConfig,PassRules[ConfigIntegral,FirstLoopConfigIntegrals,opts]]
	};


Options[ApplyMuIntegrals]={Gravity->False,ApplyIntegrals->True};


ApplyMuIntegrals[totalConfig_,coefficientMatrix_,opts:OptionsPattern[]]:=
	FirstLoopConfigIntegrals[totalConfig,PassRules[ApplyMuIntegrals,FirstLoopConfigIntegrals,opts]] .
	coefficientMatrix .
	SecondLoopConfigIntegrals[totalConfig,PassRules[ApplyMuIntegrals,FirstLoopConfigIntegrals,opts]];

ApplyMuIntegrals[_,missing_$missingConfig,opts:OptionsPattern[]]:=missing;


(* ::Subsection:: *)
(*FirstLoopConfigIntegral*)

Options[FirstLoopConfigIntegrals]=Options[ConfigIntegral];


FirstLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,firstLoopBoxPattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBoxMu4],Apply[IntegralBoxMu6],Apply[IntegralBoxMu8]}[ConstructFirstLoopConfig@totalConfig]],
	Through[{Apply[IntegralBoxMu4]}[ConstructFirstLoopConfig@totalConfig]]
];


FirstLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,firstLoopTrianglePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralTriangleMu2],Apply[IntegralTriangleMu4],Apply[IntegralTriangleMu6]}[ConstructFirstLoopConfig@totalConfig]],
	Through[{Apply[IntegralTriangleMu2]}[ConstructFirstLoopConfig@totalConfig]]
];


FirstLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,firstLoopBubblePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[ConstructFirstLoopConfig@totalConfig]],
	Through[{Apply[IntegralBubbleMu2]}[ConstructFirstLoopConfig@totalConfig]]
];


FirstLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,firstLoopBubbleTrianglePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[ConstructFirstLoopBaseConfig@totalConfig]],
	Through[{Apply[IntegralBubbleMu2]}[ConstructFirstLoopBaseConfig@totalConfig]]
];


(* ::Subsection:: *)
(*SecondLoopConfigIntegral*)


Options[SecondLoopConfigIntegrals]=Options[ConfigIntegral];


SecondLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,secondLoopBoxPattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBoxMu4],Apply[IntegralBoxMu6],Apply[IntegralBoxMu8]}[ConstructSecondLoopConfig@totalConfig]],
	Through[{Apply[IntegralBoxMu4]}[ConstructSecondLoopConfig@totalConfig]]
];


SecondLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,secondLoopTrianglePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralTriangleMu2],Apply[IntegralTriangleMu4],Apply[IntegralTriangleMu6]}[ConstructSecondLoopConfig@totalConfig]],
	Through[{Apply[IntegralTriangleMu2]}[ConstructSecondLoopConfig@totalConfig]]
];


SecondLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,secondLoopBubblePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[ConstructSecondLoopConfig@totalConfig]],
	Through[{Apply[IntegralBubbleMu2]}[ConstructSecondLoopConfig@totalConfig]]
];


SecondLoopConfigIntegrals[totalConfig_,opts:OptionsPattern[]]/;MatchQ[totalConfig,secondLoopBubbleTrianglePattern] := 
If[OptionValue[Gravity]===True,
	Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[ConstructSecondLoopBaseConfig@totalConfig]],
	Through[{Apply[IntegralBubbleMu2]}[ConstructSecondLoopBaseConfig@totalConfig]]
];


(* ::Chapter:: *)
(*Postamble*)


End[];


Protect[ "SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`*" ];    (* protect exported symbols *)
Unprotect[$LogFileTemplate];

EndPackage[ ];  (* end the package context *)
