(* ::Package:: *)

(* ::Title:: *)
(*Amplitudes.wl*)


(* ::Chapter:: *)
(*Preamble*)


(* :Title: Amplitudes.m -- a package collecting various amplitudes *)

(* :Context: SpinorHelicity6D` *)



BeginPackage["SpinorHelicity6D`Amplitudes`",
	{
		"SpinorHelicity6D`",
		"SpinorHelicity6D`Numerics`",
		"SpinorHelicity6D`SpinorFunctions`",
		"SpinorHelicity6D`BCFWAnalytics`",
		"Utils`",
		"SpinorHelicity6D`AnalyticTools`",
		"SpinorHelicity6D`BerendsGiele`"
	}
];

Unprotect["SpinorHelicity6D`Amplitudes`*"];
ClearAll["SpinorHelicity6D`Amplitudes`*"];
ClearAll["SpinorHelicity6D`Amplitudes`Private`*"];

(* usage messages for the exported functions and the context itself *)

ATree::usage = "Returns the tree amplitude with give (6D) helicities and the specified momenta, if such an amplitude exists in the collection.";
ATreeN::usage = "Returns amplitude in a form that is partially or completely evaluated, with some optimization.";
A1LAllPlusnPt::usage ="";
A1LOneMinus4Pt::usage = "";
A1LOneMinus5Pt::usage = "";
A1LOneMinus6Pt::usage = "";
ATree::noamp = "No amplitude available for ATree[`1`]";

ATreeGluon::usage = "";
ATreeContact::usage = "";
ATreeGluon::noamp = "No amplitude available for ATreeGluon[`1`]";
ATreeContact::noamp = "No amplitude available for ATreeContact[`1`]";

AOneLoop::usage = "";
AOneLoop::noamp = "No amplitude available for AOneLoop[`1`]";
AOneLoopGluon::usage = "";
AOneLoopGluon::noamp = "No amplitude available for AOneLoopGluon[`1`]";
AOneLoopContact::usage = "";
AOneLoopContact::noamp = "No amplitude available for AOneLoopContact[`1`]";

MTree::usage = "Returns the tree amplitude with give (6D) helicities and the specified momenta, if such an amplitude exists in the collection.";
MTree::noamp = "No amplitude available for MTree[`1`]";

MOneLoop::usage = "";
MOneLoop::noamp = "No amplitude available for MOneLoop[`1`]";

MTwoLoopRemainderFreiburg::usage = "";
MTwoLoopRationalFreiburg::usage = "";

A4GGGG6D::usage = "";
ATreeSPSScalarGluonShift::usage = "For BCFW computations";
ATreeSMSScalarGluonShift::usage = "For BCFW computations";

TriInt::usage ="";
A2LnPtDiv::usage ="";
A2LnPtDivPartial::usage ="";
F2ME::usage ="";
K2K4ST::usage ="";
cri::usage ="";
A2LnPtPoly::usage ="";
A2LnPtPolyPartial::usage ="";

A52LRatLitBadger::usage ="";
A52LRatLitBadgerPartial::usage ="";
A52LRatLitDunbar::usage ="";
A52LRatLitDunbarPartial::usage ="";
A62LRatLitBadger::usage ="";
A62LRatLitBadgerPartial::usage = "";
A62LRatLitDunbar::usage ="";
A62LRatLitDunbarPartial::usage = "";
A72LRatLitDunbar::usage ="";
A72LRatLitDunbarPartial::usage = "";

R241::usage = "Four-point two-Loop all-plus rational term with colour structure N_c^2 Tr[xxxx]";
R243::usage = "Four-point two-Loop all-plus rational term with colour structure N_c Tr[xx]Tr[xx]";
R241B::usage = "Four-point two-Loop all-plus rational term with colour structure N_c^0 Tr[xxxx]";
R251::usage = "Five-point two-Loop all-plus rational term with colour structure N_c^2 Tr[xxxxx]";
R253::usage = "Five-point two-Loop all-plus rational term with colour structure N_c Tr[xx]Tr[xxx]";
R251B::usage = "Five-point two-Loop all-plus rational term with colour structure N_c^0 Tr[xxxxx]";
R261::usage = "Six-point two-Loop all-plus rational term with colour structure N_c^2 Tr[xxxxxx]";
R263::usage = "Six-point two-Loop all-plus rational term with colour structure N_c Tr[xx]Tr[xxxx]";
R264::usage = "Six-point two-Loop all-plus rational term with colour structure N_c Tr[xxx]Tr[xxx]";
R2622::usage = "Six-point two-Loop all-plus rational term with colour structure Tr[xx]Tr[xx]Tr[xx]";
R261B::usage = "Six-point two-Loop all-plus rational term with colour structure N_c^0 Tr[xxxxxx]";
R271::usage = "Seven-point two-Loop all-plus rational term with colour structure N_c^2 Tr[xxxxxxx]";
R2n1B::usage = "";


VerifyKK2ScalarAmplitudes::usage = "";
VerifyKK4ScalarAmplitudes::usage = "";


$amplitudesNumericalOptimize::usage = "List of amplitudes that are optimized for numerical evaluation.";
TranslHel::usage = "Rules for translating between 6D helicity index and the usual convention."

RecomputeAmplitudes::usage = "Recomputes analytic forms of all defined amplitudes (for example if they cannot be loaded from a file.";
LoadAmplitudesFile::usage = "Loads file with amplitude definitions";
SaveAmplitudesToFile::usage = "Saves amplitudes to file";
Overwrite::usage = "Option for SaveAmplitudesToFile to overwrite any existing file found.";

$amplitudesLoaded::usage = "Variable that specifies, whether the amplitudes are available, either by loading a save file, or recomputation.";
$amplitudesFilePath::usage = "List of paths where to look for an amplitudes save file."; 

$amplitudesLoaded = False;

$AmplitudesReference::usage = "4D massless momentum that is used when a reference momentum is needed.";

KillMasses[{$AmplitudesReference}];
GenSpinors[
	{$AmplitudesReference},
	FourD->{$AmplitudesReference},
	RandomSpinors->True,
	Seed->2356,
	ParameterRange->50,
	DisplaySpinors->False
];

AmplitudeCoefficient::usage = "";

$simplifyBCFWAmplitudes = True;

$centerNGluonMax = 3;

$forceBerendsGiele::usage = "";
$forceBerendsGiele = False;


(* error messages for the exported objects *)



Begin["`Private`"];    (* begin the private context (implementation part) *)


(* ::Text:: *)
(*All Amplitudes are defined with a factor of (-i)*)


(* ::Chapter:: *)
(*Main*)


(* ::Section:: *)
(* $PackageDirectory *)

(* Saving directory containing Amplitudes package*)
$PackageDirectory = FileNameTake[$InputFileName,{1,-2}];

(* Default location to look for amplitudesDefinitions.mx is in the package directory*)
$amplitudesFilePath = {$PackageDirectory};

(* ::Section:: *)
(* GGT Support *)


(* Checking for GGT*)
$GGTAvailable = FindFile["GGT`"]=!=$Failed;

If[$GGTAvailable,
	Print["Found GGT, using for gluon amplitudes..."];
	<<GGT`;
];

GGTToSpinorHelicity6DRules[n_] := {
	GGT`GGTspaa[a_, b_] :> SpinorAngleBracket[a, b],
	GGT`GGTspbb[a_, b_] :> SpinorSquareBracket[a, b],
	GGT`GGTspaa[a_, b__, c_] :> Chain[$angle, a, {b}, c, $angle],
	GGT`GGTspbb[a_, b__, c_] :> Chain[$square, a, {b}, c, $square],
	GGT`GGTx[a_, b_] /; b > a :> PlusM @@ Range @@ {a, b - 1},
	GGT`GGTx[a_, b_] /; b < a :> PlusM @@ Mod[Range[a, n + b - 1], n, 1],
	GGT`GGTs[a_, b_] /; b > a :> S @@ Range @@ {a, b - 1},
	GGT`GGTs[a_, b_] /; b < a :> S @@ Mod[Range[a, n + b - 1], n, 1]
};


(* ::Section:: *)
(*Helper Function*)

TranslHel={{1,1}->"+",{2,2}->"-",{1,2}->"<",{2,1}->">"};


amplitudePartsReplacementRules[ampl_]:=Module[{bareObjects,rules,atoms,atomRules},
	bareObjects = Cases[ampl,_S6Mom|_S|_Chain|_SpinorAngleBracket|_SpinorSquareBracket,Infinity]//DeleteDuplicates;
	rules = {#,canonicalizeChain@#}&/@bareObjects;
	atoms = rules[[All,2]]//DeleteDuplicates;
	atomRules = Rule@@{#,Simplify@ToNum@#}&/@atoms;
	Rule@@@MapAt[ReplaceAll[atomRules],rules,{{All,2}}]//Dispatch
];

canonicalizeChain[chain_] := chain/.{Chain[x_,p1_,{pMiddle__},pn_,y_] /; Order[{p1,pMiddle,pn},Reverse@{p1,pMiddle,pn}] == -1 :> (-1)^(1+Length@{pMiddle})*Chain[y,pn,Reverse@{pMiddle},p1,x]};


ConjugateHelicities:={$angle->$square,$square->$angle,SpinorAngleBracket->SpinorSquareBracket,SpinorSquareBracket->SpinorAngleBracket}

simplifyBCFW[expr_] /; $simplifyBCFWAmplitudes := Simplify@expr;
simplifyBCFW[expr_] := expr;


(* ::Section:: *)
(*Translation Layer*)

ATree[argsA___][argsB___] /; $forceBerendsGiele :=
    ATreeBerendsGiele[argsA][argsB];

ATree[argsA___][argsB___] /;
    $forceBerendsGiele === "Adaptive" &&
		Count[{argsA},"S"] >= 2 &&
    Count[{argsA},"s"] >= 2 &&
		Count[{argsA},"+"] >= 2 :=
			ATreeBerendsGiele[argsA][argsB];

(* ::Subsection:: *)
(*Yang-Mills*)

(* ::Subsubsection::Closed:: *)
(*Pure Gluon*)

ATree[hels : ("+" | "-") ..][momenta__] /; $GGTAvailable && Length@{hels} === Length@{momenta} &&
     Length@{hels} > 3 :=
		ReplaceKinematics[
			GGT`GGTgluon[
				Length@{momenta},
				Flatten@Position[{hels}, "-", 1]
			] //. GGTToSpinorHelicity6DRules[Length@{momenta}],
			Range@Length@{momenta},
			{momenta}
		];

ATree[p1:"+"...,"-",p2:"+"...,"-",p3:"+"...][momenta__]/; Length@{p1,p2,p3}+2 === Length@{momenta} :=
		SpinorAngleBracket[{momenta}[[Length@{p1}+1]],{momenta}[[Length[{p1,
			p2}]+2]]]^4/PTF[momenta];

ATree[p1:"-"...,"+",p2:"-"...,"+",p3:"-"...][momenta__] /; Length@{p1,p2,p3}+2 === Length@{momenta} :=
		(-1)^(Length@{momenta})SpinorSquareBracket[{momenta}[[Length@{p1}+1]],{momenta}[[Length[{p1,
			p2}]+2]]]^4/AntiPTF[momenta];


ATree[hel1_,hel2_,hel3_][p1_,p2_,p3_] /; AllTrue[{hel1,hel2,hel3},MemberQ[{"+","-","<",">"},#]&] := A3GGG[hel1,hel2,hel3][p1,p2,p3];
ATree[hel1_,hel2_,hel3_,hel4_][p1_,p2_,p3_,p4_] /; AllTrue[{hel1,hel2,hel3,hel4},MemberQ[{"+","-","<",">"},#]&] := A4GGGG[hel1,hel2,hel3,hel4][p1,p2,p3,p4];


AOneLoop[a:"+"...,"-",b:"+"...][momenta__]/;Length@{a,b}===3 := A1LOneMinus4Pt@@RotateLeft[{momenta},Length@{a}];
AOneLoop[a:"+"...,"-",b:"+"...][momenta__]/;Length@{a,b}===4 := A1LOneMinus5Pt@@RotateLeft[{momenta},Length@{a}];
AOneLoop[a:"+"...,"-",b:"+"...][momenta__]/;Length@{a,b}===5 := A1LOneMinus6Pt@@RotateLeft[{momenta},Length@{a}];

(* ::Subsubsection:: *)
(*Two Scalar*)


ATree["S",hel2_,"S"][p1_,p2_,p3_] /; AllTrue[{hel2},MemberQ[{"+","-","<",">"},#]&] := A3SGS[hel2][p1,p2,p3];
ATree["S","S",hel3_][p1_,p2_,p3_] /; AllTrue[{hel3},MemberQ[{"+","-","<",">"},#]&] := ATree["S",hel3,"S"][p2,p3,p1];

ATree["S","+","+","S"][p1_,p2_,p3_,p4_] := A4SPPS[p1,p2,p3,p4];
ATree["S",hel2_,hel3_,"S"][p1_,p2_,p3_,p4_] /; AllTrue[{hel2,hel3},MemberQ[{"+","-","<",">"},#]&] := A4SGGS[hel2,hel3][p1,p2,p3,p4];
ATree["S","S",hel3_,hel4_][p1_,p2_,p3_,p4_] := ATree["S",hel3,hel4,"S"][p2,p3,p4,p1];
ATree["S","+","S","+"][p1_,p2_,p3_,p4_] := A4SPSP[p1,p2,p3,p4];
ATree["+","S","+","S"][p1_,p2_,p3_,p4_] := ATree["S","+","S","+"][p2,p3,p4,p1];

ATree["S","+","+","+","S"][p1_,p2_,p3_,p4_,p5_] := A5SPPPSBadger[p1,p2,p3,p4,p5];
ATree["S","S","+","+","+"][p1_,p2_,p3_,p4_,p5_] := ATree["S","+","+","+","S"][p2,p3,p4,p5,p1];
ATree["S","+","+","-","S"][p1_,p2_,p3_,p4_,p5_] := A5SPPMSBadger[p1,p2,p3,p4,p5];
ATree["S","-","+","+","S"][p1_,p2_,p3_,p4_,p5_] := -A5SPPMSBadger[p5,p4,p3,p2,p1];
ATree["S","+","-","+","S"][p1_,p2_,p3_,p4_,p5_] := A5SPMPSBadger[p1,p2,p3,p4,p5];
ATree["S","+","-","-","S"][p1_,p2_,p3_,p4_,p5_] := A5SPPMSBadger[p5,p4,p3,p2,p1]/.ConjugateHelicities;
ATree["S","-","-","+","S"][p1_,p2_,p3_,p4_,p5_] := ATree["S","+","-","-","S"][p5,p4,p3,p2,p1];
ATree["S",hel2_,hel3_,hel4_,"S"][p1_,p2_,p3_,p4_,p5_] /; AllTrue[{hel2,hel3,hel4},MemberQ[{"+","-"},#]&] := A5SGGGSBCFWGG[hel2,hel3,hel4][p1,p2,p3,p4,p5];
ATree["S",hel2_,hel3_,hel4_,"S"][p1_,p2_,p3_,p4_,p5_] /; MemberQ[{"<",">"},hel2]&&AllTrue[{hel3,hel4},MemberQ[{"+","-"},#]&] := A5SGGGSBCFWGG[hel2,hel3,hel4][p1,p2,p3,p4,p5];
ATree["S",hel2_,hel3_,hel4_,"S"][p1_,p2_,p3_,p4_,p5_] /; MemberQ[{"<",">"},hel3]&&AllTrue[{hel2,hel4},MemberQ[{"+","-"},#]&] := A5SGGGSBCFWGG[hel2,hel3,hel4][p1,p2,p3,p4,p5];
ATree["S",hel2_,hel3_,hel4_,"S"][p1_,p2_,p3_,p4_,p5_] /; MemberQ[{"<",">"},hel4]&&AllTrue[{hel2,hel3},MemberQ[{"+","-"},#]&] := -A5SGGGSBCFWGG[hel4,hel3,hel2][p5,p4,p3,p2,p1];

(*ATree["S","+","+","S","+"][p1_,p2_,p3_,p4_,p5_] := A5SPPSP[p1,p2,p3,p4,p5];*)
ATree["S","+","+","S","+"][p1_,p2_,p3_,p4_,p5_] :=(-ATree["S","+","+","+","S"][p1,p2,p3,p5,p4]-
	ATree["S","+","+","+","S"][p1,p2,p5,p3,p4]-ATree["S","+","+","+","S"][p1,p5,p2,p3,p4]);
ATree["S","+","S","+","+"][p1_,p2_,p3_,p4_,p5_] := ATree["S","+","+","S","+"][p3,p4,p5,p1,p2];

ATree["S","+","+","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPPS[p1,p2,p3,p4,p5,p6];
ATree["S","S","+","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPPS[p2,p3,p4,p5,p6,p1];
ATree["S","+","+","+","S","+"][p1_,p2_,p3_,p4_,p5_,p6_] := -A6SPPPPS[p1,p2,p3,p4,p6,p5]-A6SPPPPS[p1,p2,p3,p6,p4,p5]+
	-A6SPPPPS[p1,p2,p6,p3,p4,p5]-A6SPPPPS[p1,p6,p2,p3,p4,p5];
ATree["S","+","S","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := ATree["S","+","+","+","S","+"][p3,p4,p5,p6,p1,p2];
ATree["S","+","+","S","+","+"][p1_,p2_,p3_,p4_,p5_,p6_]:= A6SPPPPS[p1,p2,p3,p6,p5,p4]+A6SPPPPS[p1,p2,p6,p3,p5,p4]+A6SPPPPS[p1,p6,p2,p3,p5,p4]+
	+A6SPPPPS[p1,p2,p6,p5,p3,p4]+A6SPPPPS[p1,p6,p2,p5,p3,p4]+A6SPPPPS[p1,p6,p5,p2,p3,p4];

ATree["S","+","+","+","-","S"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPMS[p1,p2,p3,p4,p5,p6];
ATree["S","-","+","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPMS[p6,p5,p4,p3,p2,p1];
ATree["S","S","+","+","+","-"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPMS[p2,p3,p4,p5,p6,p1];
ATree["S","S","-","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPPMS[p1,p6,p5,p4,p3,p2];


ATree["S","+","+","-","+","S"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPMPS[p1,p2,p3,p4,p5,p6];
ATree["S","+","-","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPMPS[p6,p5,p4,p3,p2,p1];
ATree["S","S","+","+","-","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPMPS[p2,p3,p4,p5,p6,p1];
ATree["S","S","+","-","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPMPS[p1,p6,p5,p4,p3,p2];


ATree["S","+","+","+","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPPPS[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","+","+","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPPPS[p2,p3,p4,p5,p6,p7,p1];
ATree["S","+","+","+","+","S","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -Sum[ATree["S","+","+","+","+","+","S"]@@Insert[{p1,p2,p3,p4,p5,p6},p7,i],{i,2,6}];
ATree["S","+","S","+","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := ATree["S","+","+","+","+","S","+"][p3,p4,p5,p6,p7,p1,p2];

ATree["S","+","+","+","S","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPSPP[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","+","S","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPSPP[p4,p5,p6,p7,p1,p2,p3];


ATree["S",plus__,"S"][p1_,pGl__,pn_] /; {plus} === ConstantArray["+",Length@{pGl}] && Length@{pGl}>=5 := AnSnGS[plus][p1,pGl,pn];


AOneLoop["S","+","+","S"][p1_,p2_,p3_,p4_]:=A1LSSPP[p4,p1,p2,p3];
AOneLoop["S","S","+","+"][p1_,p2_,p3_,p4_]:=AOneLoop["S","+","+","S"][p2,p3,p4,p1];
AOneLoopGluon["S","S","+","+"][p1_,p2_,p3_,p4_]:=(S6Mom[p2,p3]-S6Mom[p1, p3])/(6SpinorAngleBracket[p3, p4]^2);
AOneLoopContact["S","S","+","+"][p1_,p2_,p3_,p4_]:=-1/2SpinorSquareBracket[p3,p4]/SpinorAngleBracket[p3,p4];

AOneLoop["S","+","+","+","S"][p1_,p2_,p3_,p4_,p5_]:=A1LSSPPP[p5,p1,p2,p3,p4];
AOneLoopContact["S","S","+","+","+"][p1_,p2_,p3_,p4_,p5_]:=-1/2*Chain[$square, p3, {p4, UnderBar[p1]}, p5, $square]/(S6Mom[p1, p5]*SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5]) -
		Chain[$square, p3, {UnderBar[p2], p4}, p5, $square]/(2*S6Mom[p2, p3]*SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5]) - SpinorSquareBracket[p3, p5]/(2*SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5])
AOneLoopContact["S","+","+","+","S"][p1_,p2_,p3_,p4_,p5_]:=AOneLoopContact["S","S","+","+", "+"][p5,p1,p2,p3,p4];

AOneLoop["S"]["S","+","+","+"][p1_, p2_,p3_,p4_,p5_]:=
		Chain[$square, p3, {UnderBar[p1], p4}, p5, $square]/(S6Mom[p1, p3]*SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5]) -
				Chain[$square, p5, {UnderBar[p1], p4}, p3, $square]/(S6Mom[p1, p5]*SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5]) +
				SpinorSquareBracket[p3, p5]/(SpinorAngleBracket[p3, p4]*SpinorAngleBracket[p4, p5]);


(* ::Subsubsection::Closed:: *)
(*Four Scalar*)


ATree[pols:PatternSequence["S",___,"s",___,"s",___,"S",tail___]][moms__] :=
    ATree[Sequence@@RotateRight[{pols},1+Length@{tail}]][Sequence@@RotateRight[{moms}, 1+Length@{tail}]];
ATreeGluon[pols:PatternSequence["S",___,"s",___,"s",___,"S",tail___]][moms__] :=
		ATreeGluon[Sequence@@RotateRight[{pols},1+Length@{tail}]][Sequence@@RotateRight[{moms}, 1+Length@{tail}]];
ATreeContact[pols:PatternSequence["S",___,"s",___,"s",___,"S",tail___]][moms__] :=
		ATreeContact[Sequence@@RotateRight[{pols},1+Length@{tail}]][Sequence@@RotateRight[{moms}, 1+Length@{tail}]];


(* ::Subsubsubsection::Closed:: *)
(* Four-Point	*)


ATree["S","S","s","s"][p1_,p2_,p3_,p4_] := S6Mom[p1,p3]/S6Mom[p1,p2];
ATree["S","s","s","S"][p1_,p2_,p3_,p4_] := S6Mom[p2,p4]/S6Mom[p2,p3];
ATree["S","s","S","s"][p1_,p2_,p3_,p4_] := 1;

ATreeGluon["S","S","s","s"][p1_,p2_,p3_,p4_] := -(1/2+S6Mom[p2,p3]/S6Mom[p1,p2]);
ATreeGluon["S","s","s","S"][p1_,p2_,p3_,p4_] := -(1/2+S6Mom[p3,p4]/S6Mom[p2,p3]);
ATreeGluon["S","s","S","s"][p1_,p2_,p3_,p4_] := 0;

ATreeContact["S","S","s","s"][p1_,p2_,p3_,p4_] := -(1/2);
ATreeContact["S","s","s","S"][p1_,p2_,p3_,p4_] := -(1/2);
ATreeContact["S","s","S","s"][p1_,p2_,p3_,p4_] := 1;




(* ::Subsubsubsection:: *)
(*Five-Point	*)

(* Full *)
ATree["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SSPss[p1,p2,p3,p4,p5];
ATree["S","S","s","s","+"][p1_,p2_,p3_,p4_,p5_] := A5SSPss[p3,p4,p5,p1,p2];

ATree["S","+","S","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SPSss[p1,p2,p3,p4,p5];
ATree["S","S","s","+","s"][p1_,p2_,p3_,p4_,p5_] := ATree["S","+","S","s","s"][p3,p4,p5,p1,p2];

ATree["S","+","s","S","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p1,p2,p3,p4,p5];
ATree["S","s","+","S","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p2,p3,p4,p5,p1];
ATree["S","s","S","+","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p3,p4,p5,p1,p2];
ATree["S","s","S","s","+"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p4,p5,p1,p2,p3];


(* Full *)
ATreeGluon["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SSPssGluon[p1,p2,p3,p4,p5];
ATreeGluon["S","S","s","s","+"][p1_,p2_,p3_,p4_,p5_] := A5SSPssGluon[p3,p4,p5,p1,p2];

ATreeGluon["S","+","S","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SPSssGluon[p1,p2,p3,p4,p5];
ATreeGluon["S","S","s","+","s"][p1_,p2_,p3_,p4_,p5_] := ATreeGluon["S","+","S","s","s"][p3,p4,p5,p1,p2];

ATreeGluon["S","+","s","S","s"][p1_,p2_,p3_,p4_,p5_] := 0;
ATreeGluon["S","s","+","S","s"][p1_,p2_,p3_,p4_,p5_] := 0;
ATreeGluon["S","s","S","+","s"][p1_,p2_,p3_,p4_,p5_] := 0;
ATreeGluon["S","s","S","s","+"][p1_,p2_,p3_,p4_,p5_] := 0;


(* Full *)
ATreeContact["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SSPssContact[p1,p2,p3,p4,p5];
ATreeContact["S","S","s","s","+"][p1_,p2_,p3_,p4_,p5_] := A5SSPssContact[p3,p4,p5,p1,p2];

ATreeContact["S","+","S","s","s"][p1_,p2_,p3_,p4_,p5_] := A5SPSssContact[p1,p2,p3,p4,p5];
ATreeContact["S","S","s","+","s"][p1_,p2_,p3_,p4_,p5_] := ATreeContact["S","+","S","s","s"][p3,p4,p5,p1,p2];

ATreeContact["S","+","s","S","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p1,p2,p3,p4,p5];
ATreeContact["S","s","+","S","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p2,p3,p4,p5,p1];
ATreeContact["S","s","S","+","s"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p3,p4,p5,p1,p2];
ATreeContact["S","s","S","s","+"][p1_,p2_,p3_,p4_,p5_] := A5SPsSs[p4,p5,p1,p2,p3];



(* ::Subsubsubsection::Closed:: *)
(*Six-Point	*)


ATree["S","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SSPPss[p1,p2,p3,p4,p5,p6];
ATree["S","S","s","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SSPPss[p3,p4,p5,p6,p1,p2];

ATree["S","S","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SSPssP[p1,p2,p3,p4,p5,p6];

ATree["S","+","S","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPSssP[p1,p2,p3,p4,p5,p6];
ATree["S","+","S","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPSssP[p3,p2,p1,p6,p5,p4];
ATree["S","S","s","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPSssP[p5,p4,p3,p2,p1,p6];
ATree["S","S","+","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPSssP[p4,p5,p6,p1,p2,p3];

ATree["S","S","s","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6sSSsPP[p6,p1,p2,p3,p4,p5];
ATree["S","+","+","S","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6sSSsPP[p4,p5,p6,p1,p2,p3];

ATree["S","+","S","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6sSPSsP[p6,p1,p2,p3,p4,p5];

ATree["S","+","+","s","S","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPsSs[p1,p2,p3,p4,p5,p6];
ATree["S","s","+","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPsSs[p2,p3,p4,p5,p6,p1];
ATree["S","s","S","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPsSs[p3,p4,p5,p6,p1,p2];
ATree["S","s","S","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPPsSs[p4,p5,p6,p1,p2,p3];

ATree["S","+","s","S","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSsP[p1,p2,p3,p4,p5,p6];
ATree["S","+","s","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSsP[p3,p4,p5,p6,p1,p2];
ATree["S","s","+","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSsP[p4,p5,p6,p1,p2,p3];
ATree["S","s","S","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSsP[p5,p6,p1,p2,p3,p4];

ATree["S","+","s","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSPs[p1,p2,p3,p4,p5,p6];
ATree["S","s","+","S","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SPsSPs[p2,p3,p4,p5,p6,p1];


(* ::Subsubsubsection::Closed:: *)
(*Seven-Point	*)


ATree["S","S","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPss[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","s","s","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPss[p3,p4,p5,p6,p7,p1,p2];

ATree["S","S","+","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPssPP[p5,p6,p7,p1,p2,p3,p4];
ATree["S","S","+","s","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPssPP[p1,p2,p3,p4,p5,p6,p7];


ATree["S","+","S","s","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSsPPs[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","+","S","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSsPPs[p5,p6,p7,p1,p2,p3,p4];

ATree["S","+","S","+","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPsPs[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","S","s","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPsPs[p4,p5,p6,p7,p1,p2,p3];

ATree["S","+","S","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPssP[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","+","s","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPssP[p4,p5,p6,p7,p1,p2,p3];

ATree["S","+","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPPss[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","+","+","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPSPPss[p7,p6,p5,p4,p3,p2,p1];
ATree["S","S","s","+","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPSPPss[p3,p4,p5,p6,p7,p1,p2];
ATree["S","+","S","s","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPSPPss[p3,p2,p1,p7,p6,p5,p4];

ATree["S","+","+","S","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPSPss[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","+","s","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPPSPss[p7,p6,p5,p4,p3,p2,p1];
ATree["S","S","s","+","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPSPss[p3,p4,p5,p6,p7,p1,p2];
ATree["S","+","+","S","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPPSPss[p4,p3,p2,p1,p7,p6,p5];

ATree["S","+","+","+","S","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPSss[p1,p2,p3,p4,p5,p6,p7];
ATree["S","S","s","+","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPSss[p3,p4,p5,p6,p7,p1,p2];


ATree["S","+","s","+","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsPSPs[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","s","S","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsPSPs[p4,p5,p6,p7,p1,p2,p3];
ATree["S","s","+","S","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsPSPs[p2,p3,p4,p5,p6,p7,p1];
ATree["S","+","s","+","S","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsPSPs[p6,p7,p1,p2,p3,p4,p5];


ATree["S","+","s","S","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSsPP[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","+","s","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSsPP[p4,p5,p6,p7,p1,p2,p3];
ATree["S","s","+","+","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSsPP[p5,p6,p7,p1,p2,p3,p4];
ATree["S","s","S","+","+","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSsPP[p6,p7,p1,p2,p3,p4,p5];

ATree["S","+","+","s","S","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPsSsPP[p1,p7,p6,p5,p4,p3,p2];
ATree["S","+","s","+","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPsSsPP[p3,p2,p1,p7,p6,p5,p4];
ATree["S","s","+","S","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPsSsPP[p4,p3,p2,p1,p7,p6,p5];
ATree["S","s","S","+","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := -A7SPsSsPP[p5,p4,p3,p2,p1,p7,p6];

ATree["S","+","s","S","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSPPs[p1,p2,p3,p4,p5,p6,p7];
ATree["S","+","+","s","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSPPs[p5,p6,p7,p1,p2,p3,p4];
ATree["S","s","+","S","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSPPs[p2,p3,p4,p5,p6,p7,p1];
ATree["S","s","+","+","S","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPsSPPs[p6,p7,p1,p2,p3,p4,p5];

ATree["S","+","+","+","s","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPsSs[p1,p2,p3,p4,p5,p6,p7];
ATree["S","s","+","+","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPsSs[p2,p3,p4,p5,p6,p7,p1];
ATree["S","s","S","+","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPsSs[p3,p4,p5,p6,p7,p1,p2];
ATree["S","s","S","s","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SPPPsSs[p4,p5,p6,p7,p1,p2,p3];

ATree["S","s","+","+","+","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
		-2*A8SSPPPPssContact[p1,p2,p3,p4,p5,p6,p7,p8];
ATree["S","+","+","+","+","s","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
		ATree["S","s","+","+","+","+","S","s"][p8,p1,p2,p3,p4,p5,p6,p7];
ATree["S","s","S","+","+","+","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
		ATree["S","s","+","+","+","+","S","s"][p2,p3,p4,p5,p6,p7,p8,p1];
ATree["S","s","S","s","+","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
		ATree["S","s","+","+","+","+","S","s"][p3,p4,p5,p6,p7,p8,p1,p2];


AOneLoop[plus__][pGl__] /; {plus} === ConstantArray["+",Length@{pGl}] && Length@{pGl}>=4 := A1LAllPlusnPt[pGl];


(* ::Subsubsection::Closed:: *)
(*Four Scalar Gluon-Contact Separated *)

ATreeGluon[pols:PatternSequence["S",___,"s",___,"s",___,"S",tail___]][moms__] :=
		ATreeGluon[Sequence@@RotateRight[{pols},1+Length@{tail}]][Sequence@@RotateRight[{moms}, 1+Length@{tail}]];


(* ::Subsubsubsection::Closed:: *)
(*Four-Point	*)

ATreeContact["S","S","s","s"][___]:=-1/2;
ATreeGluon["S","S","s","s"][p1_,p2_,p3_,p4_]:=-(1/2+S6Mom[p2,p3]/S6Mom[p1,p2])


(* ::Subsubsubsection::Closed:: *)
(*Five-Point	*)

(*ATreeContact["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_]:=*)
(*    -1/2*Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]/*)
(*        (S6Mom[p2,p3]*S6Mom[p3, p4]);*)

(*ATreeGluon["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_]:= -1/2*(S6Mom[p2, p4]**)
(*    (Chain[$square, p3, {UnderBar[p4], UnderBar[p5]}, p3, $square]*S6Mom[p2, p3] +*)
(*		 Chain[$square, p3, {UnderBar[p1], UnderBar[p2]}, p3, $square]*S6Mom[p3, p4]) -*)
(*		Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]**)
(*      (S6Mom[p1, p2]*S6Mom[p1, p4] -*)
(*				S6Mom[p1, p2]*S6Mom[p3, p4] -*)
(*        2*S6Mom[p2, p3]*S6Mom[p3, p4] -*)
(*				S6Mom[p2, p3]*S6Mom[p4, p5] +*)
(*				S6Mom[p2, p5]*S6Mom[p4, p5] +*)
(*        S6Mom[p1, p2]S6Mom[p4, p5])*)
(*			)/(S6Mom[p1, p2]*S6Mom[p2, p3]*S6Mom[p3, p4]*S6Mom[p4, p5]);*)


(*ATreeContact["S","+","S","s","s"][p1_,p2_,p3_,p4_,p5_]:=*)
(*		-1/2*Chain[$square, p2, {UnderBar[p1], UnderBar[p3]}, p2, $square]/*)
(*      (S6Mom[p1,p2]*S6Mom[p2, p3]);*)


(* ::Subsubsubsection::Closed:: *)
(*Six-Point	*)

ATreeContact["S","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] :=
    A6SSPPssContactScalarCurrents[p1,p2,p3,p4,p5,p6];
ATreeContact["S","S","s","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] :=
    A6SSPPssContactScalarCurrents[p3,p4,p5,p6,p1,p2];

ATreeContact["S","S","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_] :=
		A6SSPssPContactScalarCurrents[p1,p2,p3,p4,p5,p6];

ATreeGluon["S","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SSPPssGluon[p1,p2,p3,p4,p5,p6];
ATreeGluon["S","S","s","s","+","+"][p1_,p2_,p3_,p4_,p5_,p6_] := A6SSPPssGluon[p3,p4,p5,p6,p1,p2];


(* ::Subsubsubsection::Closed:: *)
(*Seven-Point	*)

ATreeContact["S","S","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPssContact[p1,p2,p3,p4,p5,p6,p7];
ATreeContact["S","S","s","s","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPssContact[p3,p4,p5,p6,p7,p1,p2];

ATreeGluon["S","S","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPssGluon[p1,p2,p3,p4,p5,p6,p7];
ATreeGluon["S","S","s","s","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] := A7SSPPPssGluon[p3,p4,p5,p6,p7,p1,p2];


(* ::Subsubsubsection::Closed:: *)
(*Eight-Point	*)

ATreeContact["S","S","+","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
    A8SSPPPPssContact[p1,p2,p3,p4,p5,p6,p7,p8];
ATreeContact["S","S","s","s","+","+","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_] :=
    A8SSPPPPssContact[p3,p4,p5,p6,p7,p8,p1,p2];

(* ::Subsection::Closed:: *)
(*Gravity*)


MTree["S","+","S"][p1_,p2_,p3_]:=MSppS[p1,p2,p3];
MTree["S","-","S"][p1_,p2_,p3_]:=MSmmS[p1,p2,p3];
MTree["S","S","-"][p1_,p2_,p3_]:=MTree["S","-","S"][p1,p3,p2];
MTree["S","S","+"][p1_,p2_,p3_]:=MTree["S","+","S"][p1,p3,p2];
MTree[hel1_,hel2_,hel3_][p1_,p2_,p3_] /; AllTrue[{hel1,hel2,hel3},MemberQ[{"+","-"},#]&] := -ATree[hel1,hel2,hel3][p1,p2,p3]^2;

MTree["S","+","+","S"][p1_,p2_,p3_,p4_]:=MSppPPS[p1,p2,p3,p4];
MTree["S","+","S","+"][p1_,p2_,p3_,p4_]:=MTree["S","+","+","S"][p1,p2,p4,p3];
MTree["-","+","+","-"][p1_,p2_,p3_,p4_]:=MmmPPppMM[p1,p2,p3,p4];

MTree["S","S","s","s"][p1_,p2_,p3_,p4_]:=MSSss[p1,p2,p3,p4];


MOneLoop["+","+","+","+"][k1_,k2_,k3_,k4_]:=-(S6Mom[k1,k2]S6Mom[k2,k3]/PTF[k1,k2,k3,k4])^2(S6Mom[k1,k2]^2+S6Mom[k3,k2]^2+S6Mom[k1,k3]^2)/120;

MOneLoop[allPP__][vecs__]/;Length@{allPP}===Length@{vecs} && AllTrue[{allPP},SameQ[#,"+"]&] && 4 <= Length@{vecs} <= 5 :=
Module[{n = Length@{vecs}},
	(-1)/2(-1)^n/960*Sum[
		halfSoftFunction[ab[[1]],n,ab[[2]]]*halfSoftFunction[ab[[2]],Complement[{vecs},Join[ab,n]],ab[[1]]]*
			(Chain[$angle,ab[[1]],{PlusM@@n,ab[[2]],PlusM@@Complement[{vecs},Join[ab,n]]},ab[[1]],$square]+
			Chain[$square,ab[[1]],{PlusM@@n,ab[[2]],PlusM@@Complement[{vecs},Join[ab,n]]},ab[[1]],$angle])^3,
	{ab,Subsets[{vecs},{2}]},{i,1,n-2-1},{n,Subsets[Complement[{vecs},ab],{i}]}]
]


halfSoftFunction[a_,{c_},b_]:=1/(SpinorAngleBracket[a,c]^2*SpinorAngleBracket[c,b]^2);
halfSoftFunction[a_,{c_,d_},b_]:=1/(SpinorAngleBracket[a,c]*SpinorAngleBracket[a,d])SpinorSquareBracket[c,d]/SpinorAngleBracket[c,d]1/(SpinorAngleBracket[c,b]*SpinorAngleBracket[d,b]);


halfSoftFunction[a_,N_List,b_]/;False :=Null;


(* ::Subsubsection::Closed:: *)
(*Backwards compatibility*)


ATree["S","++","S"][p1_,p2_,p3_]:=MSppS[p1,p2,p3];
ATree["S","--","S"][p1_,p2_,p3_]:=MSmmS[p1,p2,p3];
ATree[hel1_,hel2_,hel3_][p1_,p2_,p3_] /; AllTrue[{hel1,hel2,hel3},MemberQ[{"++","--"},#]&]:=
Module[{helYM1,helYM2,helYM3},
	{helYM1,helYM2,helYM3} = {hel1,hel2,hel3}/.{"++"->"+","--"->"-"};
	-ATree[helYM1,helYM2,helYM3][p1,p2,p3]^2
];

ATree["S","++","++","S"][p1_,p2_,p3_,p4_]:=MSppPPS[p1,p2,p3,p4];
ATree["--","++","++","--"][p1_,p2_,p3_,p4_]:=MmmPPppMM[p1,p2,p3,p4];


(* ::Subsection::Closed:: *)
(*Numerical optimization for selected amplitudes*)


ATreeN[hels__][moms__]:=ATree[hels][moms];


ATreeN[hels__][moms__] /; MemberQ[$amplitudesNumericalOptimize,{hels}] ||  Count[{hels},"S"|"s"] >2 := 
	ATree[hels][moms]/.amplitudePartsReplacementRules[ATree[hels][moms]]


(* ::Subsection::Closed:: *)
(*Default Case*)

ATree[argsA___][argsB___] /; $forceBerendsGiele === "Adaptive" :=
    ATreeBerendsGiele[argsA][argsB];


ATree[hels___][moms___] /; Message[ATree::noamp, {hels}] := Null;
ATreeGluon[hels___][moms___] /; Message[ATreeGluon::noamp, {hels}] := Null;
ATreeContact[hels___][moms___] /; Message[ATreeContact::noamp, {hels}] := Null;

MTree[hels___][moms___] /; Message[MTree::noamp, {hels}] := Null;

AOneLoop[hels___][moms___] /; Message[AOneLoop::noamp, {hels}] := Null;

MOneLoop[hels___][moms___] /; Message[MOneLoop::noamp, {hels}] := Null;


(* ::Section::Closed:: *)
(*Load Amplitudes File*)


LoadAmplitudesFile[filenameArg_ : ""]:=Module[{filename},
	If[filenameArg =!= "",
		WriteString["stdout","Looking for amplitudes file at: ", filenameArg," ..."];
		Which[
			FileExistsQ[filenameArg],
				WriteString["stdout"," Found!\n"];
				filename = filenameArg,
			FileExistsQ[FileNameJoin[{filenameArg,"amplitudeDefinitions.mx"}]],
				WriteString["stdout"," Found amplitudeDefinitions.mx in folder!\n"];
				filename = FileNameJoin[{filenameArg,"amplitudeDefinitions.mx"}],
			_,
				WriteString["stdout"," Not found!\n"];
				Abort[];
		];
		,
		WriteString["stdout","Looking for amplitudes file in default locations..."];
		filename = Select[$amplitudesFilePath,FileExistsQ[FileNameJoin[{#,"amplitudeDefinitions.mx"}]]&];
		If[filename === {},
			WriteString["stdout"," No file named amplitudeDefinitions.mx found!\n"];
			Abort[];
		,
			filename = FileNameJoin[{First@filename,"amplitudeDefinitions.mx"}];
			WriteString["stdout"," File ",filename , " found!\n"];
		];
	];
	WriteString["stdout","Loading file..."];
	Get[filename];
	WriteString["stdout"," Done!\n"];
];


(* ::Section::Closed:: *)
(*Save Amplitudes To File*)


Options[SaveAmplitudesToFile] = {Overwrite -> False};
SaveAmplitudesToFile[filenameArg_ : "", opts : OptionsPattern[]] :=
Module[{filename},
	If[filenameArg =!= "",
		filename = filenameArg;
		WriteString["stdout","Saving amplitude definitions to file: ", filename," ..."];
	,
		filename = FileNameJoin[{$amplitudesFilePath[[1]],"amplitudeDefinitions.mx"}];
		WriteString["stdout","Saving amplitude definitions to default file: ", filename," ..."];
	];
	If[FileExistsQ[filename] && ! OptionValue[Overwrite],
		WriteString["stdout"," File already exists. Rerun the command with the option Overwrite->True to overwrite it.\n"];
		Abort[];
	,
		DumpSave[filename,{A3GGG,A3SGS,A4GGGG,A4SGGS,A5SGGGSBCFWGG(*,A5SSPss*),A6SSPPss,
			A6SSPPssGluon,A6SSPPssContact,A6SSPssP,A7SSPPPss,A7SSPssPP,A7SPSsPPs,
			A7SPSPsPs,A7SPSPssP,A7SPSPPss,A7SPPSPss,A7SPPPSss,A7SPPPPPS,A7SPPPSPP,
			A7SPsPSPs,A7SPsSsPP,A7SPsSPPs,A7SSPPPssGluon,A7SSPPPssContact,A7SPPPsSs,
			A8SSPPPPssContact
		}];
		WriteString["stdout"," Done!\n"];
	];
];


(* ::Section::Closed:: *)
(*Reevaluation of BCFW Amplitudes  *)




RecomputeAmplitudes[]:=Module[{},
	WriteString["stdout","Recomputing amplitudes. This will take a second..."];
	
	ClearAll[A3GGG,A3SGS,A4GGGG,A4SGGS,A5SGGGSBCFWGG(*,A5SSPss*),A6SSPPss,
		A6SSPPssGluon,A6SSPPssContact,A6SSPssP,A7SSPPPss,A7SSPssPP,A7SPSsPPs,
		A7SPSPsPs,A7SPSPssP,A7SPSPPss,A7SPPSPss,A7SPPPSss,A7SPPPPPS,A7SPPPSPP,
		A7SPsPSPs,A7SPsSsPP,A7SPsSPPs,A7SSPPPssGluon,A7SSPPPssContact,A7SPPPsSs,
		A8SSPPPPssContact
	];
	
	Activate@Map[Function[hels,Inactive[SetDelayed][A3GGG[Sequence@@hels][k1_,k2_,k3_],A3GGGEval[hels][k1,k2,k3]]],(Tuples[Tuples[{1,2},{2}],{3}])/.TranslHel];
	
	Activate@Map[Function[hels,Inactive[SetDelayed][A3SGS[Sequence@@hels][k1_,k2_,k3_],A3SGSEval[hels][k1,k2,k3]]],(Tuples[Tuples[{1,2},{2}],{1}])/.TranslHel];
	
	Activate@Map[Function[hels,Inactive[SetDelayed][A4GGGG[Sequence@@hels][k1_,k2_,k3_,k4_],A4GGGG6D[hels][k1,k2,k3,k4]]],(Tuples[Tuples[{1,2},{2}],{4}])/.TranslHel];
	
	Activate@Map[Function[hels,Inactive[SetDelayed][A4SGGS[Sequence@@hels][k1_,k2_,k3_,k4_],A4SGGSEval[hels][k1,k2,k3,k4]]],(Tuples[Tuples[{1,2},{2}],{2}])/.TranslHel];
	
	Quiet[
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","+","+"][p1_,p2_,p3_,p4_,p5_],{1,0,0,0,1} . A5SGGGSBCFWGGAnalytics["+","+","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","+","+"][p1_,p2_,p3_,p4_,p5_],{0,1,0,0,1} . A5SGGGSBCFWGGAnalytics["-","+","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["<","+","+"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,1} . A5SGGGSBCFWGGAnalytics["<","+","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG[">","+","+"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,1} . A5SGGGSBCFWGGAnalytics[">","+","+"][p1,p2,p3,p4,p5]]];
	
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","+","-"][p1_,p2_,p3_,p4_,p5_],{1,0,0,0,1} . A5SGGGSBCFWGGAnalytics["+","+","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","+","-"][p1_,p2_,p3_,p4_,p5_],{0,1,0,0,1} . A5SGGGSBCFWGGAnalytics["-","+","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["<","+","-"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,1} . A5SGGGSBCFWGGAnalytics["<","+","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG[">","+","-"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,1} . A5SGGGSBCFWGGAnalytics[">","+","-"][p1,p2,p3,p4,p5]]];
	
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","-","+"][p1_,p2_,p3_,p4_,p5_],{0,1,0,0,1} . A5SGGGSBCFWGGAnalyticsReversed["-","-","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["<","-","+"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,1} . A5SGGGSBCFWGGAnalyticsReversed["<","-","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG[">","-","+"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,1} . A5SGGGSBCFWGGAnalyticsReversed[">","-","+"][p1,p2,p3,p4,p5]]];

	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","-","-"][p1_,p2_,p3_,p4_,p5_],{1,0,0,0,1} . A5SGGGSBCFWGGAnalyticsReversed["+","-","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","-","-"][p1_,p2_,p3_,p4_,p5_],{0,1,0,0,1} . A5SGGGSBCFWGGAnalyticsReversed["-","-","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["<","-","-"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,1} . A5SGGGSBCFWGGAnalyticsReversed["<","-","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG[">","-","-"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,1} . A5SGGGSBCFWGGAnalyticsReversed[">","-","-"][p1,p2,p3,p4,p5]]];

	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","-","+"][p1_,p2_,p3_,p4_,p5_],{0,1,0,0,1,0,0,0,1,1} . A5SGGGSBCFWGGNAdjAnalytics["+","-","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","<","+"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,0,0,0,0,1,1} . A5SGGGSBCFWGGNAdjAnalytics["+","<","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+",">","+"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,0,0,0,0,1,1} . A5SGGGSBCFWGGNAdjAnalytics["+",">","+"][p1,p2,p3,p4,p5]]];

	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+","<","-"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,0,0,1,0,1,1} . A5SGGGSBCFWGGNAdjAnalytics["+","<","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["+",">","-"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,0,0,0,1,1,1} . A5SGGGSBCFWGGNAdjAnalytics["+",">","-"][p1,p2,p3,p4,p5]]];

	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","<","+"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,0,0,1,0,1,1} . A5SGGGSBCFWGGNAdjReversedAnalytics["-","<","+"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-",">","+"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,0,0,0,1,1,1} . A5SGGGSBCFWGGNAdjReversedAnalytics["-",">","+"][p1,p2,p3,p4,p5]]];
	
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-","<","-"][p1_,p2_,p3_,p4_,p5_],{0,0,0,1,0,0,0,0,1,1} . A5SGGGSBCFWGGNAdjReversedAnalytics["-","<","-"][p1,p2,p3,p4,p5]]];
	Activate[Inactive[SetDelayed][A5SGGGSBCFWGG["-",">","-"][p1_,p2_,p3_,p4_,p5_],{0,0,1,0,0,0,0,0,1,1} . A5SGGGSBCFWGGNAdjReversedAnalytics["-",">","-"][p1,p2,p3,p4,p5]]];
	
	
(*	Activate@Inactive[SetDelayed][A5SSPss[p1_,p2_,p3_,p4_,p5_],A5SSPssBCFWSG[p1,p2,p3,p4,p5]];*)
(*	Activate@Inactive[SetDelayed][A5SSPss[p1_,p2_,p3_,p4_,p5_],A5SSPssFeynman[p1,p2,p3,p4,p5]];*)
	
	Activate@Inactive[SetDelayed][A6SSPPss[p1_,p2_,p3_,p4_,p5_,p6_],A6SSPPssBCFWSG[p1,p2,p3,p4,p5,p6]];
	Activate@Inactive[SetDelayed][A6SSPPssGluon[p1_,p2_,p3_,p4_,p5_,p6_],A6SSPPssGluonBCFWSG[p1,p2,p3,p4,p5,p6]];
	Activate@Inactive[SetDelayed][A6SSPPssContact[p1_,p2_,p3_,p4_,p5_,p6_],A6SSPPssContactBCFWSG[p1,p2,p3,p4,p5,p6]];
	Activate@Inactive[SetDelayed][A6SSPssP[p1_,p2_,p3_,p4_,p5_,p6_],A6SSPSSPBCFWSG[p1,p2,p3,p4,p5,p6]];
	
	Activate@Inactive[SetDelayed][A7SSPPPss[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SSPPPSSBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SSPPPssGluon[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SSPPPssGluonBCFWSG[p1, p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SSPPPssContact[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SSPPPssContactBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SSPssPP[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SSPSSPPBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	
	Activate@Inactive[SetDelayed][A7SPSsPPs[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPSsPPsBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPSPsPs[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPSPsPsBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPSPssP[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPSPssPBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPSPPss[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPSPPssBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPPSPss[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPPSPssBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPPPSss[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPPPPSssBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	
	Activate@Inactive[SetDelayed][A7SPPPPPS[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPPPPPSBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPPPSPP[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPPPSPPBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	
	Activate@Inactive[SetDelayed][A7SPsPSPs[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPsPSPsBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPsSsPP[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPsSsPPBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPsSPPs[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPsSPPsBCFWSG[p1,p2,p3,p4,p5,p6,p7]];
	Activate@Inactive[SetDelayed][A7SPPPsSs[p1_,p2_,p3_,p4_,p5_,p6_,p7_],A7SPPPsSsBCFWSG[p1,p2,p3,p4,p5,p6,p7]];

	Activate@Inactive[SetDelayed][A8SSPPPPssContact[p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_],A8SSPPPPssContactBCFWSG[p1,p2,p3,p4,p5,p6,p7,p8]];
	
	
	];
	WriteString["stdout"," Done!\n"];
];


(* ::Section::Closed:: *)
(*Tree Amplitude Definitions Yang-Mills*)


(* ::Subsection::Closed:: *)
(*Amplitude 3: GGG*)


A3GGGRaw[n1_,n2_,n3_,k_][a_,adot_,b_,bdot_,c_,cdot_]:=1/(S[n1,k]S[n2,k]S[n3,k])*
(AngAngInvChain[n1,k,n2][a,b]AngAngInvChain[n2,k,n1][bdot,adot]AngSquInvChain[n3,n1,k,n3][c,cdot]+
AngAngInvChain[n2,k,n3][b,c]AngAngInvChain[n3,k,n2][cdot,bdot]AngSquInvChain[n1,n2,k,n1][a,adot]+
AngAngInvChain[n3,k,n1][c,a]AngAngInvChain[n1,k,n3][adot,cdot]AngSquInvChain[n2,n3,k,n2][b,bdot]);


ClearAll[A3GGGEval];
A3GGGEval[configSymbolic_][k1_,k2_,k3_] /; Length@configSymbolic ==3 :=Module[{res,config=configSymbolic/.(Reverse/@TranslHel),k,p1,p2,p3},
KillMasses[{k}];
res=A3GGGRaw[p1,p2,p3,k][a,adot,b,bdot,c,cdot]//To4D;
res=res/.{a->config[[1,1]],adot->config[[1,2]],b->config[[2,1]],bdot->config[[2,2]],c->config[[3,1]],cdot->config[[3,2]]}//To4D;
res=SpinorReplace[res,{SpinorUndotPure[p1][$mu]->SpinorUndotPure[k][$lam],SpinorUndotPure[p2][$mu]->SpinorUndotPure[k][$lam],SpinorUndotPure[p3][$mu]->SpinorUndotPure[k][$lam],
SpinorDotPure[p1][$mu]->SpinorDotPure[k][$lam],SpinorDotPure[p2][$mu]->SpinorDotPure[k][$lam],SpinorDotPure[p3][$mu]->SpinorDotPure[k][$lam]}];
res=res/.{s[p1,k]->-Spinorsquarebracket[p1,k]Spinoranglebracket[p1,k],s[p2,k]->-Spinorsquarebracket[p2,k]Spinoranglebracket[p2,k],s[p3,k]->-Spinorsquarebracket[p3,k]Spinoranglebracket[p3,k]}//Simplify;
res=res/.{Spinoranglebracket[p2,p3] Spinorsquarebracket[p2,k]->-Spinoranglebracket[p1,p3]Spinorsquarebracket[p1,k]};
res=res/.{Spinoranglebracket[p1,p3] Spinoranglebracket[p2,k] Spinorsquarebracket[p1,k]->-Spinoranglebracket[p2,p3] Spinoranglebracket[p2,k] Spinorsquarebracket[p2,k]}//Simplify;
res=res//.{(Spinoranglebracket[p1,k] Spinorsquarebracket[p1,k]+Spinoranglebracket[p3,k] Spinorsquarebracket[p3,k])->-Spinoranglebracket[p2,k] Spinorsquarebracket[p2,k]};
res=res/.{Spinorsquarebracket[k,p2]->-Spinoranglebracket[p3,p1]/Spinoranglebracket[p3,p2]Spinorsquarebracket[k,p1],Spinorsquarebracket[k,p3]->-Spinoranglebracket[p2,p1]/Spinoranglebracket[p2,p3]Spinorsquarebracket[k,p1],Spinoranglebracket[k,p2]->-Spinorsquarebracket[p3,p1]/Spinorsquarebracket[p3,p2]Spinoranglebracket[k,p1],Spinoranglebracket[k,p3]->-Spinorsquarebracket[p2,p1]/Spinorsquarebracket[p2,p3]Spinoranglebracket[k,p1]}//Simplify;
(*res=res/.{extramass[p2]->-extramass[p1]-extramass[p3],extramasstilde[p2]->-extramasstilde[p1]-extramasstilde[p3]}//Simplify;*)
res=res/.{Spinoranglebracket[p1,p3] Spinorsquarebracket[p1,p3]+Spinoranglebracket[p2,p3] Spinorsquarebracket[p2,p3]->-Spinoranglebracket[p1,p2] Spinorsquarebracket[p1,p2],Spinoranglebracket[p1,p2] Spinorsquarebracket[p1,p2]+Spinoranglebracket[p2,p3] Spinorsquarebracket[p2,p3]->-Spinoranglebracket[p1,p3] Spinorsquarebracket[p1,p3],Spinoranglebracket[p1,p2] Spinorsquarebracket[p1,p2]+Spinoranglebracket[p1,p3] Spinorsquarebracket[p1,p3]->-Spinoranglebracket[p2,p3] Spinorsquarebracket[p2,p3]}//Expand;
res = res//.{HoldPattern[SpinorAngleBracket[i_,j_]]/;Or[i===k,j===k]:>S[i,j]/SpinorSquareBracket[j,i]};
res = res /.{S[k,p3]->S[k,p1]*S[p1,p3]/S[p2,p3],S[k,p2]->S[k,p1]S[p1,p2]/S[p2,p3]}/.{S[i_,j_]:>SpinorSquareBracket[i,j]SpinorAngleBracket[j,i]};
res /. MapThread[Rule,{{p1,p2,p3,k},{k1,k2,k3,$AmplitudesReference}}]
];



(* ::Subsection::Closed:: *)
(*Amplitude 3: SGS*)


(* ::Subsubsection::Closed:: *)
(*All Helicities*)


ClearAll[A3SGS,A3SGSEval,A3SGSRaw];
A3SGSRaw[p1_,p2_,p3_,q_][b_,bdot_]:=1/S[p2,q]*AngSquInvChain[p2,p1,q,p2][b,bdot];


A3SGSEval[configSymbolic_][k1_,k2_,k3_] /; Length@configSymbolic == 1 :=Module[{res,config=configSymbolic/.(Reverse/@TranslHel),q,p1,p2,p3},
KillMasses[{q,p2}];
res = A3SGSRaw[p1,p2,p3,q][Sequence@@Flatten[config]]//To4D;
res=SpinorReplace[res,{SpinorUndotPure[p1][$mu]->SpinorUndotPure[q][$lam],SpinorUndotPure[p2][$mu]->SpinorUndotPure[q][$lam],SpinorUndotPure[p3][$mu]->SpinorUndotPure[q][$lam],
		SpinorDotPure[p1][$mu]->SpinorDotPure[q][$lam],SpinorDotPure[p2][$mu]->SpinorDotPure[q][$lam],SpinorDotPure[p3][$mu]->SpinorDotPure[q][$lam]}];
res = res/.{S[i_,j_]:>SpinorAngleBracket[i,j]SpinorSquareBracket[j,i]};
res = Divide@(Sequence@@ToChain/@NumeratorDenominator[res]);
res = res /.{Chain[a_,b_,{Except[_UnderBar,e_]},c_,d_]:>Chain[a,b,{UnderBar[e]},c,d]};
res /. MapThread[Rule,{{p1,p2,p3,q},{k1,k2,k3,$AmplitudesReference}}]
];


(* ::Subsubsection:: *)
(*Vertices for Scalar Gluon BCFW Computations*)


ATreeSPSScalarGluonShift[extScalar_,extGluon_,sign_,shiftSMom_,shiftGMom_]:=
	-sign*Chain[$square,extGluon,{UnderBar[extScalar],UnderBar@shiftSMom},shiftGMom,
		$square]/Chain[$angle,extGluon,{UnderBar@shiftSMom},shiftGMom,$square];

ATreeSMSScalarGluonShift[extScalar_,extGluon_,sign_,shiftSMom_,shiftGMom_]:=
	sign*Chain[$angle,extGluon,{UnderBar[extScalar]},shiftGMom,$square]/SpinorSquareBracket[extGluon,shiftGMom];


(* ::Subsection::Closed:: *)
(*Amplitude 4: GGGG*)


(* ::Subsubsection:: *)
(*All Helicities*)


A4GGGGRaw[k_,l_,n_,m_][a_,adot_,b_,bdot_,c_,cdot_,d_,ddot_]:=-1/(S6Mom[k,l]S6Mom[l,n])AngAngInvariant[k,l,n,m][a,b,c,d]SquSquInvariant[k,l,n,m][adot,bdot,cdot,ddot];


NewProcess;
A4GGGG6D[configSymbolic_][k1_,k2_,k3_,k4_] /; Length@configSymbolic == 4 :=Module[{res,config=configSymbolic/.(Reverse/@TranslHel),k,p1,p2,p3,p4},
	KillMasses[{k}];
		res=A4GGGGRaw[p1,p2,p3,p4][Sequence@@Flatten@config]//To4D;
		res=res/.{a->config[[1,1]],adot->config[[1,2]],b->config[[2,1]],bdot->config[[2,2]],c->config[[3,1]],cdot->config[[3,2]]}//To4D;
		res=SpinorReplace[res,{SpinorUndotPure[p1][$mu]->SpinorUndotPure[k][$lam],SpinorUndotPure[p2][$mu]->SpinorUndotPure[k][$lam],SpinorUndotPure[p3][$mu]->SpinorUndotPure[k][$lam],SpinorUndotPure[p4][$mu]->SpinorUndotPure[k][$lam],
		SpinorDotPure[p1][$mu]->SpinorDotPure[k][$lam],SpinorDotPure[p2][$mu]->SpinorDotPure[k][$lam],SpinorDotPure[p3][$mu]->SpinorDotPure[k][$lam],SpinorDotPure[p4][$mu]->SpinorDotPure[k][$lam]}];
		res = Factor[Simplify@res]//.{Spinoranglebracket[k,i_]Spinorsquarebracket[k,i_]:>-S[k,i],Power[SpinorAngleBracket[k,i_],-1]*Power[SpinorSquareBracket[k,i_],-1]:>-Power[S[k,i],-1],
			Spinoranglebracket[k,i_] Spinoranglebracket[k,j_] Spinorsquarebracket[i_,j_]:>-Chain[$angle,k,{UnderBar[i],UnderBar[j]},k,$angle],Spinoranglebracket[i_,j_] Spinorsquarebracket[k,i_] Spinorsquarebracket[k,j_]:>-Chain[$square,k,{UnderBar[i],UnderBar[j]},k,$square]};
		res = Simplify[res];
		res /. {S->S6Mom};
		res /. MapThread[Rule,{{p1,p2,p3,p4,k},{k1,k2,k3,k4,$AmplitudesReference}}]
	];


(* ::Subsection::Closed:: *)
(*Amplitude 4: SS'S'S*)


A4SSss[p1S_,p2S_,p3s_,p4s_]:=S6Mom[p1S,p3s]/S6Mom[p1S,p2S];


(* ::Subsection::Closed:: *)
(*Amplitude 4: SGGS*)


(* ::Subsubsection::Closed:: *)
(*All Helicities*)


A4SGGSRaw[p1_,p2_,p3_,p4_][a2_,adot2_,a3_,adot3_]:=Module[{cdot,ddot,c,d},-1/(4*S[p3,p4]S[p2,p3])SumContracted[c,d][AngAngInvariant[p2,p3,p4,p4][a2,a3,c,d]*levicivita2Up[c,d]]*SumContracted[cdot,ddot][levicivita2Up[cdot,ddot]SquSquInvariant[p2,p3,p1,p1][adot2,adot3,cdot,ddot]]];


A4SGGSEval[configSymbolic_][k1_,k2_,k3_,k4_] /; Length@configSymbolic == 2 :=Module[{res,config=configSymbolic/.(Reverse/@TranslHel),q,p1,p2,p3,p4},
	KillMasses[{q,p2,p3}];
	res = A4SGGSRaw[p1,p2,p3,p4][Sequence@@Flatten[config]]//To4D;
	res=SpinorReplace[res,{SpinorUndotPure[p1][$mu]->SpinorUndotPure[q][$lam],SpinorUndotPure[p2][$mu]->SpinorUndotPure[q][$lam],SpinorUndotPure[p3][$mu]->SpinorUndotPure[q][$lam],SpinorUndotPure[p4][$mu]->SpinorUndotPure[q][$lam],
			SpinorDotPure[p1][$mu]->SpinorDotPure[q][$lam],SpinorDotPure[p2][$mu]->SpinorDotPure[q][$lam],SpinorDotPure[p3][$mu]->SpinorDotPure[q][$lam],SpinorDotPure[p4][$mu]->SpinorDotPure[q][$lam]}];
	res = res/.{HoldPattern[SpinorAngleBracket[p1,i_]SpinorSquareBracket[p1,j_]]:>-Chain[$angle,i,{UnderBar[p1]},j,$square]+Extramass[p1]Extramasstilde[p1]SpinorAngleBracket[i,q]SpinorSquareBracket[j,q]/SpinorAngleBracket[p1,q]*1/SpinorSquareBracket[p1,q],
				HoldPattern[SpinorAngleBracket[i_,p4]SpinorSquareBracket[j_,p4]]:>-Chain[$angle,i,{UnderBar[p4]},j,$square]+Extramass[p4]Extramasstilde[p4]SpinorAngleBracket[i,q]SpinorSquareBracket[j,q]/SpinorAngleBracket[p4,q]*1/SpinorSquareBracket[p4,q]};
	res = res/.{S->S6Mom};
	res /. MapThread[Rule,{{p1,p2,p3,p4,q},{k1,k2,k3,k4,$AmplitudesReference}}] //Simplify
];


(* ::Subsection::Closed:: *)
(*Amplitude 4: SGGS*)

(* ::Subsubsection::Closed:: *)
(*SPPS*)

A4SPPS[p1_,p2_,p3_,p4_]:=
		-(Extramass[p1]*Extramasstilde[p1]*SpinorSquareBracket[p2, p3]^2)/(S6Mom[p2, p3]*S6Mom[p3, p4]);

(* ::Subsubsection::Closed:: *)
(*SPSP*)


A4SPSP[p1_,p2_,p3_,p4_]:=-Extramass[p1]Extramasstilde[p1]SpinorSquareBracket[p2,p4]^2/(S6Mom[p1,p2]S6Mom[p2,p3]);


(* ::Subsection::Closed:: *)
(*Amplitude 5: SGGGS*)


(* ::Subsubsection::Closed:: *)
(*A5SPPPS Literature*)


(* ::Text:: *)
(*Results from hep-th/0504159*)


A5SPPPSBadger[p1s_,p2g_,p3g_,p4g_,p5s_]:=-(Extramass[p5s]*Extramasstilde[p5s]Chain[$square,p4g,
	{PlusM[p2g,p3g],UnderBar[p1s]},p2g,$square]*
										  1/(S6Mom[p1s,p2g]*S6Mom[p4g,p5s]*SpinorAngleBracket[p2g,
												p3g]*SpinorAngleBracket[p3g,p4g]));


(* ::Subsubsection::Closed:: *)
(*A5SPPMS Literature*)


A5SPPMSBadger[p1_,p2_,p3_,p4_,p5_]:=
1/(Chain[$square,p4,{PlusM[p2,p3],UnderBar[p1]},p2,$square])*(
	(Chain[$angle,p4,{UnderBar[p5],PlusM[p2,p3],UnderBar[p1]},p2,$square])^2/(SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]S6Mom[p1,p2]S6Mom[p4,p5])-
	Extramass[p5]Extramasstilde[p5]SpinorSquareBracket[p2,p3]^3/(SpinorSquareBracket[p3,p4]S6Mom[p1,p5])
)


(* ::Subsubsection:: *)
(*A5SPMPS Literature*)


A5SPMPSBadger[p1_,p2_,p3_,p4_,p5_]:=-1/Chain[$square,p4,{PlusM[p2,p3],UnderBar[p1]},p2,$square]*(
	-Chain[$angle,p3,{UnderBar[p1]},p2,$square]^2Chain[$angle,p3,{UnderBar[p5]},p4,$square]^2/
		(SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]S6Mom[p1,p2]S6Mom[p4,p5])+
	Extramass[p1]Extramasstilde[p1]SpinorSquareBracket[p2,p4]^4/(S6Mom[p1,p5]SpinorSquareBracket[p2,p3]SpinorSquareBracket[p3,p4]));


(* ::Subsubsection::Closed:: *)
(*BCFW SGGGS Gluon-Gluon*)


NewProcess;
A5SGGGSBCFWGGAnalytics[hel2_,"+",hel4_][k1_,k2_,k3_,k4_,k5_] := Module[{p1,p2,p3,p4,p5},
		({
			BCFWRulesMasslessChannel[p3,p4,{p2,OverHat[p3]}][ATree[hel2,"+","-"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","+",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p3,p4,{p2,OverHat[p3]}][ATree[hel2,"+","+"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","-",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p3,p4,{p2,OverHat[p3]}][-ATree[hel2,"+","<"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S",">",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p3,p4,{p2,OverHat[p3]}][-ATree[hel2,"+",">"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","<",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			 
			BCFWRulesMassiveChannel[p3,p4,{p5,OverHat[p4]}][ATree["S",hel2,"+","S"][p1,p2,OverHat[p3],OverHat[K]]*(-1)/S6Mom[p4,p5]*ATree["S",hel4,"S"][pm[OverHat[K]],OverHat[p4],p5]]
		}(*/.{UnderBar[x_]/;MemberQ[{p2,p3,p4},x]\[RuleDelayed]x}*)/.{Chain[_,i_,{j_},i_,_]:>(S6Mom[i,j]/.UnderBar->Identity),Extramass[p5]->-Extramass[p1],
			Extramasstilde[p5]->-Extramasstilde[p1],S->S6Mom})/.{Indeterminate->0,ComplexInfinity->0}/.MapThread[Rule,{{p1,p2,p3,p4,p5},{k1,k2,k3,k4,k5}}]
];

A5SGGGSBCFWGGAnalyticsReversed[hel2_,"-",hel4_][k1_,k2_,k3_,k4_,k5_] := Module[{p1,p2,p3,p4,p5},
		({
			BCFWRulesMasslessChannel[p4,p3,{p2,OverHat[p3]}][ATree[hel2,"-","-"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","+",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p3,{p2,OverHat[p3]}][ATree[hel2,"-","+"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","-",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p3,{p2,OverHat[p3]}][-ATree[hel2,"-","<"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S",">",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p3,{p2,OverHat[p3]}][-ATree[hel2,"-",">"][p2,OverHat[p3],pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","<",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			 
			BCFWRulesMassiveChannel[p4,p3,{p5,OverHat[p4]}][ATree["S",hel2,"-","S"][p1,p2,OverHat[p3],OverHat[K]]*(-1)/S6Mom[p4,p5]*ATree["S",hel4,"S"][pm[OverHat[K]],OverHat[p4],p5]]
		}(*/.{UnderBar[x_]/;MemberQ[{p2,p3,p4},x]\[RuleDelayed]x}*)/.{Chain[_,i_,{j_},i_,_]:>(S6Mom[i,j]/.UnderBar->Identity),Extramass[p5]->-Extramass[p1],
			Extramasstilde[p5]->-Extramasstilde[p1],S->S6Mom})/.{Indeterminate->0,ComplexInfinity->0}/.MapThread[Rule,{{p1,p2,p3,p4,p5},{k1,k2,k3,k4,k5}}]
];

A5SGGGSBCFWGGNAdjAnalytics["+",hel3_,hel4_][k1_,k2_,k3_,k4_,k5_] := Module[{p1,p2,p3,p4,p5},
		({
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p2]}][ATree["+",hel3,"-"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","+",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p2]}][ATree["+",hel3,"+"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","-",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p2]}][-ATree["+",hel3,"<"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S",">",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p2]}][-A3GGG["+",hel3,">"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","<",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p4]}][ATree["S","+","+","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["-",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p4]}][ATree["S","+","-","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["+",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p4]}][-ATree["S","+","<","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree[">",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p2,p4,{p3,OverHat[p4]}][-ATree["S","+",">","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["<",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			
			BCFWRulesMassiveChannel[p2,p4,{p5,OverHat[p4]}][A4SGGS["+",hel3][p1,OverHat[p2],p3,OverHat[K]](-1/S6Mom[p4,p5])ATree["S",hel4,"S"][pm[OverHat[K]],OverHat[p4],p5]],
			
			BCFWRulesMassiveChannel[p2,p4,{p1,OverHat[p2]}][ATree["S","+","S"][p1,OverHat[p2],pm[OverHat[K]]](-1/S6Mom[p2,p1])ATree["S",hel3,hel4,"S"][OverHat[K],p3,OverHat[p4],p5]]
		}(*/.{UnderBar[x_]/;MemberQ[{p2,p3,p4},x]\[RuleDelayed]x}*)/.{Chain[_,i_,{j_},i_,_]:>(S6Mom[i,j]/.UnderBar->Identity),Extramass[p5]->-Extramass[p1],
			Extramasstilde[p5]->-Extramasstilde[p1],S->S6Mom})/.{Indeterminate->0,ComplexInfinity->0}/.MapThread[Rule,{{p1,p2,p3,p4,p5},{k1,k2,k3,k4,k5}}]
];

A5SGGGSBCFWGGNAdjReversedAnalytics["-",hel3_,hel4_][k1_,k2_,k3_,k4_,k5_] := Module[{p1,p2,p3,p4,p5},
		({
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p2]}][ATree["-",hel3,"-"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","+",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p2]}][ATree["-",hel3,"+"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","-",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p2]}][-ATree["-",hel3,"<"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S",">",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p2]}][-A3GGG["-",hel3,">"][OverHat[p2],p3,pm[OverHat[K]]](-1/S6Mom[p2,p3])ATree["S","<",hel4,"S"][p1,OverHat[K],OverHat[p4],p5]],
			
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p4]}][ATree["S","-","+","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["-",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p4]}][ATree["S","-","-","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["+",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p4]}][-ATree["S","-","<","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree[">",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			BCFWRulesMasslessChannel[p4,p2,{p3,OverHat[p4]}][-ATree["S","-",">","S"][p1,OverHat[p2],OverHat[K],p5](-1/S6Mom[p4,p3])ATree["<",hel3,hel4][pm[OverHat[K]],p3,OverHat[p4]]],
			
			BCFWRulesMassiveChannel[p4,p2,{p5,OverHat[p4]}][A4SGGS["-",hel3][p1,OverHat[p2],p3,OverHat[K]](-1/S6Mom[p4,p5])ATree["S",hel4,"S"][pm[OverHat[K]],OverHat[p4],p5]],
			
			BCFWRulesMassiveChannel[p4,p2,{p1,OverHat[p2]}][ATree["S","-","S"][p1,OverHat[p2],pm[OverHat[K]]](-1/S6Mom[p2,p1])ATree["S",hel3,hel4,"S"][OverHat[K],p3,OverHat[p4],p5]]
		}(*/.{UnderBar[x_]/;MemberQ[{p2,p3,p4},x]\[RuleDelayed]x}*)/.{Chain[_,i_,{j_},i_,_]:>(S6Mom[i,j]/.UnderBar->Identity),Extramass[p5]->-Extramass[p1],
			Extramasstilde[p5]->-Extramasstilde[p1],S->S6Mom})/.{Indeterminate->0,ComplexInfinity->0}/.MapThread[Rule,{{p1,p2,p3,p4,p5},{k1,k2,k3,k4,k5}}]
];




(* ::Subsection::Closed:: *)
(*Amplitude 5: SGGSG*)


(* ::Subsubsection:: *)
(*A5SPPSP*)


A5SPPSP[p1_,p2_,p3_,p4_,p5_]:=-Extramass[p1]Extramasstilde[p1]/(SpinorAngleBracket[p2,p3])*
	1/(S6Mom[p1,p5]Chain[$square,p5,{UnderBar[p3],UnderBar[p4]},p5,$square]+S6Mom[p3,p4]Chain[$square,p5,{UnderBar[p1],UnderBar[p4]},p5,$square])*
	(
		Chain[$square,p5,{UnderBar[p2],UnderBar[p3]},p5,$square]SpinorSquareBracket[p2,p5]SpinorSquareBracket[p3,p5]Extramass[p1]Extramasstilde[p1]/
			(S6Mom[p1,p2]S6Mom[p3,p4])+
		Chain[$square,p5,{UnderBar[p1],UnderBar[p4]},p5,$square]^2SpinorSquareBracket[p2,p3]/
			(S6Mom[p1,p5]S6Mom[p4,p5])
	);


(* ::Subsection:: *)
(*Amplitude 5: SSGS'S'*)


(* ::Subsubsection::Closed:: *)
(*SS+SS Feynman Diagrams*)



A5SSPssFeynmanRulesNumerator[k1_,k2_,k3_,k4_,k5_] := {
	<|
		"Rules"->{
			{
				Chain[$square,k3,{UnderBar[k5]},$AmplitudesReference,$angle]:>
						-Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle]-
								Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle]-
								Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle]
			},
			{
				S6Mom[k1,k3]->S6Mom[k4,k5]-S6Mom[k2,k3]-S6Mom[k1,k2],
				S6Mom[k1,k4]->S6Mom[k2,k3]-S6Mom[k4,k5]-S6Mom[k1,k5],
				S6Mom[k2,k4]->S6Mom[k1,k5]-S6Mom[k2,k3]-S6Mom[k3,k4],
				S6Mom[k2,k5]->S6Mom[k3,k4]-S6Mom[k1,k5]-S6Mom[k1,k2],
				S6Mom[k3,k5]->S6Mom[k1,k2]-S6Mom[k3,k4]-S6Mom[k4,k5]
			}
		},
		"Post"->Numerator@*Together
	|>,
	<|
		"Pre"->(MultiCollect[#,{
			Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle],
			Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle],
			Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle]},
			Simplify]&),
		"Rules" -> {
			Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle] S6Mom[k2,k3]->
					-SpinorAngleBracket[k3,$AmplitudesReference]*Chain[$square,k3,{UnderBar[k2],UnderBar[k1]},k3,
						$square]+Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle]S6Mom[k3,k1],
			Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle] S6Mom[k2,k3]->
					- SpinorAngleBracket[k3,$AmplitudesReference]*Chain[$square,k3,{UnderBar[k2],UnderBar[k4]},k3,
						$square]+Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle]S6Mom[k3,k4]
		}
	|>,
	<|
		"Pre"->(MultiCollect[#,{
			Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle],
			Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle],
			Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle]},
			Simplify]&),
		"Rules" -> {
			{
				(-S6Mom[k1,k5]+S6Mom[k2,k3]+S6Mom[k3,k4])->-S6Mom[k2,k4],
				(S6Mom[k1,k2]+S6Mom[k1,k3]+S6Mom[k2,k3]-S6Mom[k4,k5])->0
			},
			{
				(S6Mom[k2,k3] S6Mom[k3,k4]+(S6Mom[k1,k5]-S6Mom[k3,k4]) S6Mom[k4,k5]+S6Mom[k1,k2] (S6Mom[k3,k4]+S6Mom[k4,k5]))->
						-S6Mom[k2,k5] S6Mom[k4,k5]+S6Mom[k3,k4] (-S6Mom[k1,k3]+S6Mom[k4,k5])
			}
		},
		"Post"->((1/2(#-ReplaceKinematics[#,{k1,k2,k3,k4,k5},{k5,k4,k3,k2,k1}])//SortChains//Simplify)&)
	|>,
	<|
		"Rules"->{
			S6Mom[k1,k2] (S6Mom[k1,k4]-S6Mom[k2,k3])+S6Mom[k1,k3] S6Mom[k3,k4]+S6Mom[k2,k3] S6Mom[k3,k5]+
					S6Mom[k2,k5] S6Mom[k4,k5]-S6Mom[k3,k4] S6Mom[k4,k5]->
					-(-S6Mom[k1,k2] S6Mom[k1,k4]+S6Mom[k1,k2] S6Mom[k3,k4]+
							2 S6Mom[k2,k3] S6Mom[k3,k4]+S6Mom[k2,k3] S6Mom[k4,k5]-
							S6Mom[k2,k5] S6Mom[k4,k5])
		},
		"Post"->Identity
	|>
};


A5SSPssFeynmanRulesDenominator[k1_,k2_,k3_,k4_,k5_] := {
	<|
		"Rules"->{
			{
				Chain[$square,k3,{UnderBar[k5]},$AmplitudesReference,$angle]:>
						-Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle]-
								Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle]-
								Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle]
			},
			{
				S6Mom[k1,k3]->S6Mom[k4,k5]-S6Mom[k2,k3]-S6Mom[k1,k2],
				S6Mom[k1,k4]->S6Mom[k2,k3]-S6Mom[k4,k5]-S6Mom[k1,k5],
				S6Mom[k2,k4]->S6Mom[k1,k5]-S6Mom[k2,k3]-S6Mom[k3,k4],
				S6Mom[k2,k5]->S6Mom[k3,k4]-S6Mom[k1,k5]-S6Mom[k1,k2],
				S6Mom[k3,k5]->S6Mom[k1,k2]-S6Mom[k3,k4]-S6Mom[k4,k5]
			}
		},
		"Post"->Denominator@*Together
	|>
};


A5SSPssFeynman[p1_,p2_,p3_,p4_,p5_] := Module[
	{res = 0,s12,s34,s23,s45,chn34,chn35,chn31,chn32,spaa3k,k1,k2,k3,k4,k5,resultUnprocessed,resultNum,resultDenom},
	s45 = S6Mom[k4,k5];
	s12 = S6Mom[k1,k2];
	s34 = S6Mom[k4,k3];
	s23 = S6Mom[k2,k3];
	chn34 = Chain[$square,k3,{UnderBar[k4]},$AmplitudesReference,$angle];
	chn35 = Chain[$square,k3,{UnderBar[k5]},$AmplitudesReference,$angle];
	chn31 = Chain[$square,k3,{UnderBar[k1]},$AmplitudesReference,$angle];
	chn32 = Chain[$square,k3,{UnderBar[k2]},$AmplitudesReference,$angle];
	spaa3k = SpinorAngleBracket[k3,$AmplitudesReference];
	res += 1/2*1/(s45)*(chn35-chn34);
	res += 1/2*1/(s12)*(chn32-chn31);
	res += 1/(s45*s23)chn32*(S6Mom[k1,k5]-S6Mom[k1,k4]);
	res += 1/(s12*s34)chn34*(S6Mom[k5,k2]-S6Mom[k5,k1]);
	res += 1/(s45*s12)(chn34+chn35)/2*(S6Mom[k5,k2]-S6Mom[k5,k1]-S6Mom[k4,k2]+S6Mom[k4,k1]);
	res += -1/(s45*s12)(chn32-chn31)/2*(S6Mom[k5,k3]-S6Mom[k4,k3]);
	res += 1/(s45*s12)(chn35-chn34)/2*(S6Mom[k3,k2]-S6Mom[k3,k1]);
	res += -1/(s34)*(chn34);
	res += 1/(s23)*(chn32);

	resultUnprocessed = 1/2 * 1/spaa3k * res;

	resultNum=ApplyAnalyticTransformationRules[resultUnprocessed,A5SSPssFeynmanRulesNumerator[k1,k2,k3,k4,k5]];
	resultDenom=ApplyAnalyticTransformationRules[resultUnprocessed,A5SSPssFeynmanRulesDenominator[k1,k2,k3,k4,k5]];

	-ReplaceKinematics[
		Divide@@{resultNum,resultDenom},
		{k1,k2,k3,k4,k5},{p1,p2,p3,p4,p5}
	]
];

A5SSPss[p1_,p2_,p3_,p4_,p5_] := -1/2*(
	S6Mom[p2, p4]*(
		Chain[$square, p3, {UnderBar[p4], UnderBar[p5]}, p3, $square]*S6Mom[p2, p3] +
		Chain[$square, p3, {UnderBar[p1], UnderBar[p2]}, p3, $square]*S6Mom[p3, p4]
	) - Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]*(
		S6Mom[p1, p2]*S6Mom[p1, p4] - S6Mom[p1, p2]*S6Mom[p3, p4] -
		2*S6Mom[p2, p3]*S6Mom[p3, p4] - S6Mom[p2, p3]*S6Mom[p4, p5] +
		S6Mom[p2, p5]*S6Mom[p4, p5])
	)/(S6Mom[p1, p2]*S6Mom[p2, p3]*S6Mom[p3, p4]*S6Mom[p4, p5]);

A5SSPssGluon[p1_,p2_,p3_,p4_,p5_] := -1/2*(
	S6Mom[p2, p4]*(
		Chain[$square, p3, {UnderBar[p4], UnderBar[p5]}, p3, $square]*S6Mom[p2, p3] +
				Chain[$square, p3, {UnderBar[p1], UnderBar[p2]}, p3, $square]*S6Mom[p3, p4]
	) - Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]*(
		S6Mom[p1, p2]*S6Mom[p1, p4] - S6Mom[p1, p2]*S6Mom[p3, p4] -
				2*S6Mom[p2, p3]*S6Mom[p3, p4] - S6Mom[p2, p3]*S6Mom[p4, p5] +
				S6Mom[p2, p5]*S6Mom[p4, p5]+S6Mom[p1, p2]*S6Mom[p4, p5])
)/(S6Mom[p1, p2]*S6Mom[p2, p3]*S6Mom[p3, p4]*S6Mom[p4, p5]);

A5SSPssContact[p1_,p2_,p3_,p4_,p5_] :=
		 -Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]/(2*S6Mom[p2, p3]*S6Mom[p3, p4]);

(* ::Subsubsection:: *)
(*BCFW Scalar-Gluon*)


A5SSPssBCFWSG[p1_,p2_,p3_,p4_,p5_]:=Module[{k1,k2,k3,k4,k5,K},
-Plus@@{ApplyBcfwScalarGluonShiftRules[
	ATree["S","S","s","s"][k5,OverHat[k1],OverHat[K],k4](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
{k2,k3},k1,k2,{k2},K],
ApplyBcfwScalarGluonShiftRules[
	ATreeSPSScalarGluonShift[k5,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4])ATree["S","-","+","S"][k4,pm@OverHat[K],OverHat[k2],k3],
{k2,k3,k4},k1,k2,{k2},K],
ApplyBcfwScalarGluonShiftRules[
	ATreeSMSScalarGluonShift[k5,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][k4,pm@OverHat[K],OverHat[k2],k3],
{k2,k3,k4},k1,k2,{k2},K]}/.
	{PlusM[x_]:>x}/.AlignExtramassRules[{k5,k1},{k3,k4}]/.
	MapThread[Rule,{{k1,k2,k3,k4,k5},{p2,p3,p4,p5,p1}}]//ReduceMomenta[{p1,p2,p3,p4,p5},{p3}]//SortChains//simplifyBCFW
];


(* ::Subsection::Closed:: *)
(*Amplitude 5: SGSS'S'*)


A5SPSss[p4_,p5_,p1_,p2_,p3_]:=Chain[$square,p5,{UnderBar[p1],UnderBar[p4]},p5,$square]/(S6Mom[p4,p5]S6Mom[p5,p1])*
	(1+1/S6Mom[p2,p3]*(2Mom4DVector[p3][Null][$up] . Mom4DVector[p4][Null][$down]+S6Mom[p5,p1]Chain[$square,p5,{UnderBar[p3],UnderBar[p4]},p5,$square]/Chain[$square,p5,{UnderBar[p1],UnderBar[p4]},p5,$square]))

A5SPSssGluon[p1_,p2_,p3_,p4_,p5_] := -1/2*((Chain[$square, p2, {UnderBar[p4], UnderBar[p3]}, p2,
	$square]	- Chain[$square, p2, {UnderBar[p5], UnderBar[p3]}, p2, $square])*S6Mom[p1, p2] +
		Chain[$square, p2, {UnderBar[p1], UnderBar[p3]}, p2, $square]*(S6Mom[p3, p4] - S6Mom[p3, p5]))/
      (S6Mom[p1, p2]*S6Mom[p2, p3]*S6Mom[p4, p5])

A5SPSssContact[p1_,p2_,p3_,p4_,p5_] := A5SSPssContact[p5,p1,p2,p3,p4];

(* ::Subsection::Closed:: *)
(*Amplitude 5: SGS'SS'*)


A5SPsSs[p1_,p2_,p3_,p4_,p5_]:=-Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]/
    (S6Mom[p1,p2]S6Mom[p2,p3]);



(* ::Subsection::Closed:: *)
(*Amplitude 6: SGGGGS*)


(* ::Subsubsection::Closed:: *)
(*A6S++++S Literature*)


A6SPPPPS[p1s_,p2g_,p3g_,p4g_,p5g_,p6s_]:=-Extramass[p6s]*Extramasstilde[p6s]/(S6Mom[p1s,p2g]*S6Mom[p1s,p2g,p3g]*S6Mom[p1s,p2g,p3g,p4g]*SpinorAngleBracket[p2g,p3g]*SpinorAngleBracket[p3g,p4g]*SpinorAngleBracket[p4g,p5g])*
	(Chain[$square,p5g,{UnderBar[p6s],PlusM[p4g,p5g],PlusM[p2g,p3g],UnderBar[p1s]},p2g,$square]);


(* ::Subsubsection::Closed:: *)
(*A6S+++-S Literature*)


(* IMPORTANT: the sign of the second term is flipped compared to the result in hep-th/0504159,
	but should be correct *)


A6SPPPMS[l1_,p1_,p2_,p3_,p4_,l2_]:=
	Plus@@{(S6Mom[l1,p1,p2]Chain[$angle,p4,{UnderBar@l2,PlusM[p1,p2,p3],UnderBar@l1},p1,$square]-Extramass[l1]Extramasstilde[l1]Chain[$angle,p4,{UnderBar@l2,p3,p2},p1,$square])^2/
		(S6Mom[l1,p1]S6Mom[l1,p1,p2]S6Mom[l1,p1,p2,p3]SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]*
			Chain[$square,p4,{UnderBar@l2,PlusM[p3,p4],PlusM[p1,p2],UnderBar@l1},p1,$square])+
	-Extramass[l1]Extramasstilde[l1]Chain[$square,p3,{PlusM[p1,p2],UnderBar@l1},p1,$square]^3/
		(S6Mom[l1,p1]SpinorAngleBracket[p1,p2]SpinorSquareBracket[p3,p4]Chain[$angle,p2,{PlusM[p3,p4],PlusM[UnderBar[l1],UnderBar[l2]],UnderBar@l1},p1,$square]*
			Chain[$square,p4,{UnderBar@l2,PlusM[p3,p4],PlusM[p1,p2],UnderBar@l1},p1,$square]),
	-Extramass[l1]Extramasstilde[l1]Chain[$angle,p4,{PlusM[p3,p2]},p1,$square]^3/
		(S6Mom[p1,p2,p3,p4]S6Mom[p2,p3,p4]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]Chain[$angle,p2,{PlusM[p3,p4],PlusM[UnderBar[l1],UnderBar[l2]],UnderBar@l1},p1,$square])};


(* ::Subsubsection::Closed:: *)
(*A6S++-+S Literature*)


(* The first two terms correspond to the {p1,p2} channel in a l1,p1 shift. Simons expressions for this term appear to be be wrong, so I added my own BCFW result *)
A6SPPMPS[l1_,p1_,p2_,p3_,p4_,l2_]:=
	Plus@@{
(*		(S6Mom[l1,p1,p2]Chain[$angle,p3,{UnderBar@l1},p1,$square]-Extramass[l1]Extramasstilde[l1]SpinorSquareBracket[p3,p2]SpinorSquareBracket[p2,p1])^2*
			Chain[$angle,p3,{UnderBar@l2},p4,$square]^2/
				(S6Mom[l1,p1]S6Mom[l1,p1,p2]S6Mom[l1,p1,p2,p3]SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]*
					Chain[$square,p4,{UnderBar@l2,PlusM[p3,p4],PlusM[p1,p2],UnderBar@l1},p1,$square]),
		Extramass[l1]Extramasstilde[l1]Chain[$square,p4,{PlusM[p1,p2],UnderBar@l1},p1,$square]^4/
			(S6Mom[l1,p1]SpinorAngleBracket[p1,p2]SpinorSquareBracket[p3,p4]Chain[$square,p3,{PlusM[p1,p2],UnderBar@l1},p1,$square]*
				Chain[$angle,p2,{PlusM[p3,p4],PlusM[UnderBar[l1],UnderBar[l2]],UnderBar[l1]},p1,$square]Chain[$square,p4,{UnderBar@l2,PlusM[p3,p4],PlusM[p1,p2],UnderBar@l1},p1,$square]),
		*)
		-(SpinorSquareBracket[p1, p2]^3*
  (-((Chain[$angle, p3, {UnderBar[l2]}, p4, $square]^2*Chain[$angle, p3, {PlusM[p1, p2, UnderBar[l1]], PlusM[p1, p2], UnderBar[l1]}, p1, $square]^2*
      Chain[$square, p1, {p2, UnderBar[l1]}, p1, $square])/(Chain[$angle, p3, {p2}, p1, $square]*S6Mom[l2, p4]*S6Mom[l1, p1, p2]*
      SpinorAngleBracket[p3, p4])) - (Chain[$square, p4, {PlusM[p1, p2], UnderBar[l1]}, p1, $square]^4*Extramass[l1]*Extramasstilde[l1])/
    (Chain[$square, p3, {PlusM[p1, p2], UnderBar[l1]}, p1, $square]*(S6Mom[l1, l2] + (Chain[$square, p1, {UnderBar[l2], UnderBar[l1]}, p1, $square]*
        S6Mom[p1, p2])/Chain[$square, p1, {p2, UnderBar[l1]}, p1, $square])*SpinorSquareBracket[p3, p4])))/
 (Chain[$square, p1, {p2, UnderBar[l1]}, p1, $square]*Chain[$square, p2, {p1, UnderBar[l1]}, p1, $square]*S6Mom[p1, p2]*
  (Chain[$square, p4, {p3, PlusM[p1, p2, UnderBar[l1]], PlusM[p1, p2], UnderBar[l1]}, p1, $square] + 
   Chain[$square, p4, {PlusM[p1, p2], PlusM[p1, p2, UnderBar[l1]], PlusM[p1, p2], UnderBar[l1]}, p1, $square] + 
   (Chain[$square, p1, {UnderBar[l1], PlusM[p1, p2, UnderBar[l1]], PlusM[p1, p2], UnderBar[l1]}, p1, $square]*S6Mom[p1, p2]*SpinorSquareBracket[p1, p4])/
    Chain[$square, p1, {UnderBar[l1], p2}, p1, $square])),
		
		-Extramass[l1]Extramasstilde[l1]Chain[$angle,p3,{PlusM[p2,p4]},p1,$square]^4/
			(S6Mom[p1,p2,p3,p4]S6Mom[p2,p3,p4]Chain[$angle,p4,{PlusM[p2,p3]},p1,$square]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]*
				Chain[$angle,p2,{PlusM[p3,p4],PlusM[UnderBar[l1],UnderBar[l2]],UnderBar[l1]},p1,$square]),
		+Extramass[l1]Extramasstilde[l1]Chain[$square,p4,{PlusM[p1,p2,p3],UnderBar@l1},p1,$square]SpinorSquareBracket[p1,p2]^3/
			(S6Mom[p1,p2,p3]Chain[$angle,p4,{PlusM[p2,p3]},p1,$square]S6Mom[l1,p1,p2,p3]SpinorSquareBracket[p2,p3]*
				Chain[$square,p3,{PlusM[p1,p2],UnderBar@l1},p1,$square])
	}


(* ::Subsection:: *)
(*Amplitude 6: SSGGS'S'*)


(* ::Subsubsection:: *)
(* ATreeContact from Scalar Currents*)

A6SSPPssContactScalarCurrents[p1_,p2_,p3_,p4_,p5_,p6_]:=
		-1/2*(Chain[$square, p3, {UnderBar[p2], UnderBar[p5]}, p4, $square]/(S6Mom[p2, p3]*S6Mom[p4, p5]) -
				((Extramass[p2]*Extramasstilde[p2])/(S6Mom[p2, p3]*S6Mom[p2, p3, p4]) +
        (Extramass[p5]*Extramasstilde[p5])/(S6Mom[p4, p5]*S6Mom[p3, p4, p5]))*
            SpinorSquareBracket[p3, p4])/SpinorAngleBracket[p3, p4];

A6SSPssPContactScalarCurrents[p1_,p2_,p3_,p4_,p5_,p6_]:=
		-1/2*Chain[$square, p3, {UnderBar[p2], UnderBar[p4]}, p3, $square]*
      Chain[$square, p6, {UnderBar[p5], UnderBar[p1]}, p6, $square]/
					(S6Mom[p2,p3]S6Mom[p3,p4]S6Mom[p5,p6]S6Mom[p6,p1]);

(* ::Subsubsection:: *)
(*BCFW Scalar Gluon*)


A6SSPPssBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6,K},
	-Plus@@simplifyBCFW[{
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","S","+","s","s"][k6,OverHat[k1],OverHat[K],k4,k5](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
		{k2,k3},k1,k2,{k2,k3},K],
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","S","s","s"][k6,OverHat[k1],OverHat[K],k5](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4],
		{k2,k3,k4},k1,k2,{k2,k3},K]/.{Extramasstilde[K]->Extramasstilde[k4], Extramass[K]->Extramass[k4]},
		ApplyBcfwScalarGluonShiftRules[
			ATreeSMSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][k5,pm@OverHat[K],OverHat[k2],k3,k4],
		{k2,k3,k4,k5},k1,k2,{k2,k3},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeSPSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","-","+","+","S"][k5,pm@OverHat[K],OverHat[k2],k3,k4],
		{k2,k3,k4,k5},k1,k2,{k2,k3},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k6,k1},{k4,k5}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6},{k2,k3}]//SortChains]/.
		MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p2,p3,p4,p5,p6,p1}}]
];

A6SSPPssContactBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6,K},
	-Plus@@({
		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","+","s","s"][k6,OverHat[k1],OverHat[K],k4,k5](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
			{k2,k3},k1,k2,{k2,k3},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","s","s"][k6,OverHat[k1],OverHat[K],k5](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4],
			{k2,k3,k4},k1,k2,{k2,k3},K]/.{Extramasstilde[K]->Extramasstilde[k4], Extramass[K]->Extramass[k4]}
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k6,k1},{k4,k5}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6},{k2,k3}]//SortChains)/.
			MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p2,p3,p4,p5,p6,p1}}]
];

A6SSPPssGluonBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6,K},
	-Plus@@({
		ApplyBcfwScalarGluonShiftRules[
			ATreeGluon["S","S","+","s","s"][k6,OverHat[k1],OverHat[K],k4,k5](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
			{k2,k3},k1,k2,{k2,k3},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeGluon["S","S","s","s"][k6,OverHat[k1],OverHat[K],k5](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4],
			{k2,k3,k4},k1,k2,{k2,k3},K]/.{Extramasstilde[K]->Extramasstilde[k4], Extramass[K]->Extramass[k4]},
		ApplyBcfwScalarGluonShiftRules[
			ATreeSMSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][k5,pm@OverHat[K],OverHat[k2],k3,k4],
			{k2,k3,k4,k5},k1,k2,{k2,k3},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeSPSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","-","+","+","S"][k5,pm@OverHat[K],OverHat[k2],k3,k4],
			{k2,k3,k4,k5},k1,k2,{k2,k3},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k6,k1},{k4,k5}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6},{k2,k3}]//SortChains)/.
			MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p2,p3,p4,p5,p6,p1}}]
];



(* ::Subsection:: *)
(*Amplitude 6: SS+S'S'+*)


(* ::Subsubsection:: *)
(*BCFW Scalar Gluon*)


A6SSPSSPBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6,K},
	-Plus@@simplifyBCFW[{
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","S","+","s","s"][OverHat[K],k4,k5,k6,OverHat@k1](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
		{k2,k3},k1,k2,{k2,k5},K],
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","+","+","S"][OverHat[k1],OverHat[K],k5,k6](1/S6Mom[k2,k3,k4])ATree["S","-","+","S"][k4,pm@OverHat@K,OverHat@k2,k3],
		{k2,k3,k4},k1,k2,{k2,k5},K],
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","-","+","S"][OverHat[k1],OverHat[K],k5,k6](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][k4,pm@OverHat@K,OverHat@k2,k3],
			{k2,k3,k4},k1,k2,{k2,k5},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeSMSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][k4,k5,pm@OverHat[K],OverHat[k2],k3],
		{k2,k3,k4,k5},k1,k2,{k2,k5},K],
		ApplyBcfwScalarGluonShiftRules[
			ATreeSPSScalarGluonShift[k6,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5])ATree["S","+","-","+","S"][k4,k5,pm@OverHat[K],OverHat[k2],k3],
		{k2,k3,k4,k5},k1,k2,{k2,k5},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k6,k1},{k3,k4}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6},{k2,k5}]//SortChains]/.
		MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p2,p3,p4,p5,p6,p1}}]
];


(* ::Subsection::Closed:: *)
(*Amplitude 6: S+SS'S'+*)


ASPSssPChannel61[p1_,p2_,p3_,p4_,p5_,p6_]:=-1/S6Mom[p1,p6]Extramass[p1]Extramasstilde[p1]SpinorSquareBracket[p6,p2]/Chain[$angle,p6,{UnderBar[p1]},p2,$square]*
	Chain[$square,p2,{UnderBar[p3],PlusM[p6,UnderBar[p1]]},p2,$square]/(S6Mom[p6,p1,p2]S6Mom[OverHat[p2],p3])*
	(1+
		(S6Mom[p5,p6,OverHat[p1]]+
			S6Mom[OverHat[p2],p3]Chain[$square,p2,{UnderBar[p5],PlusM[p6,UnderBar[p1]]},p2,$square]/Chain[$square,p2,{UnderBar[p3],PlusM[p6,UnderBar[p1]]},p2,$square])/
		S6Mom[p4,p5]
	)/.{
		S6Mom[p5,p6,OverHat[p1]]->
			S6Mom[p5,p6]+2*Mom4DVector[p5][Null][$up] . Mom4DVector[p1][Null][$down]-S6Mom[p1,p6]Chain[$square,p2,{UnderBar[p5],UnderBar[p1]},p2,$square]/Chain[$square,p2,{p6,UnderBar[p1]},p2,$square],
		S6Mom[p3,OverHat[p2]]->
			S6Mom[p2,p3]+S6Mom[p1,p6]Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]/Chain[$square,p2,{p6,UnderBar[p1]},p2,$square]
		}
		
ASPSssPChannel23Delta[p1_,p2_,p3_,p4_,p5_,p6_]:=
	S6Mom[p6,OverHat[p1]]Chain[$square,p6,{UnderBar[p4],UnderBar[p5]},p6,$square]+S6Mom[p4,p5]Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]/.
	{
		S6Mom[p6,OverHat[p1]]->
			S6Mom[p6,p1]+S6Mom[p2,p3]Chain[$square,p2,{p6,UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square],
		Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]->
			Chain[$square,p6,{UnderBar[p1],UnderBar[p5]},p6,$square]+S6Mom[p2,p3]SpinorSquareBracket[p6,p2]*
				Chain[$square,p6,{UnderBar[p5],UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]
	};

ASPSssPChannel23[p1_,p2_,p3_,p4_,p5_,p6_]:=
Module[{K},
	1/(S6Mom[p2,p3]S6Mom[p1,p2])Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]*
	(
		1/(ASPSssPChannel23Delta[p1,p2,p3,p4,p5,p6]S6Mom[p4,p5]S6Mom[p1,p2,p3])*
			(
				Extramass[p4]Extramasstilde[p4]Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[OverHat[K]]},p6,$square]^2+
				Extramass[p1]Extramasstilde[p1]Chain[$square,p6,{UnderBar[p5],UnderBar[p4]},p6,$square]^2
			)-
		Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]/(S6Mom[p5,p6]S6Mom[p6,OverHat[p1]])*
			(1+S6Mom[p5,p6,OverHat[p1]]*Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]/ASPSssPChannel23Delta[p1,p2,p3,p4,p5,p6])
	)/.
	{
		S6Mom[p6,OverHat[p1]]->
			S6Mom[p6,p1]+S6Mom[p2,p3]Chain[$square,p2,{p6,UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square],
		Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]->
			Chain[$square,p6,{UnderBar[p1],UnderBar[p5]},p6,$square]+S6Mom[p2,p3]SpinorSquareBracket[p6,p2]*
				Chain[$square,p6,{UnderBar[p5],UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square],
		S6Mom[p5,p6,OverHat[p1]]->
			2*(Mom4DVector[p5][Null][$up] . Mom4DVector[p6][Null][$down]+Mom4DVector[p5][Null][$up] . Mom4DVector[p1][Null][$down]+Mom4DVector[p1][Null][$up] . Mom4DVector[p6][Null][$down])+
				S6Mom[p2,p3]Chain[$square,p2,{PlusM[UnderBar[p5],p6],UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square],
		Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[OverHat[K]]},p6,$square]->
			Chain[$square,p6,{UnderBar[p1],PlusM[p2,UnderBar[p3]]},p6,$square]-
				SpinorSquareBracket[p2,p6]S6Mom[p2,p3]Chain[$square,p6,{PlusM[UnderBar[p1],p2,UnderBar[p3]],UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]
	}
]


A6SPSssP[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6},
	ASPSssPChannel61[k1,k2,k3,k4,k5,k6]+ASPSssPChannel23[k1,k2,k3,k4,k5,k6]/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p1,p2,p3,p4,p5,p6}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 6: SSS'++S'*)


A6sSSsPPChannel45[p1_,p2_,p3_,p4_,p5_,p6_]:=
Module[{K},
	1/S6Mom[p4,p5]Chain[$square,p5,{UnderBar[p4]},p6,$angle]/(SpinorAngleBracket[p5,p6]S6Mom[p4,p5,p6]S6Mom[p1,OverHat[p6]])*
		Chain[$square,OverHat[p6],{UnderBar[p1],UnderBar[OverHat[K]]},OverHat[p6],$square]*
		(1+(S6Mom[p3,OverHat[K]]+
			S6Mom[p1,OverHat[p6]]*Chain[$square,OverHat[p6],{UnderBar[p3],UnderBar[OverHat[K]]},OverHat[p6],$square]/
			Chain[$square,OverHat[p6],{UnderBar[p1],UnderBar[OverHat[K]]},OverHat[p6],$square])/S6Mom[p2,p3]
		)/.
	{
		S6Mom[p3,OverHat[K]]->
			2Mom4DVector[p3][Null][$up] . Mom4DVector[p4][Null][$down]+S6Mom[p3,p5]-
				S6Mom[p4,p5]Chain[$square,p5,{UnderBar[p3]},p6,$angle]/Chain[$square,p5,{UnderBar[p4]},p6,$angle],
		S6Mom[p1,OverHat[p6]]->
			S6Mom[p1,p6]+S6Mom[p4,p5]Chain[$square,p5,{UnderBar[p1]},p6,$angle]/Chain[$square,p5,{UnderBar[p4]},p6,$angle],
		Chain[$square,OverHat[p6],{UnderBar[p1],UnderBar[OverHat[K]]},OverHat[p6],$square]->
			-(Chain[$square,p6,{UnderBar[p1],PlusM[UnderBar[p2],UnderBar[p3]]},p6,$square]+
				S6Mom[p4,p5]/Chain[$square,p5,{UnderBar[p4]},p6,$angle]*
				(
					Chain[$square,p5,{UnderBar[p1],PlusM[UnderBar[p2],UnderBar[p3]]},p6,$square]+
					Chain[$square,p6,{UnderBar[p1],PlusM[UnderBar[p2],UnderBar[p3]]},p5,$square]
				)+
				S6Mom[p4,p5]^2/Chain[$square,p5,{UnderBar[p4]},p6,$angle]^2*
					Chain[$square,p5,{UnderBar[p1],PlusM[UnderBar[p2],UnderBar[p3]]},p5,$square]
			),
		Chain[$square,OverHat[p6],{UnderBar[p3],UnderBar[OverHat[K]]},OverHat[p6],$square]->
			-(Chain[$square,p6,{UnderBar[p3],PlusM[UnderBar[p2],UnderBar[p1]]},p6,$square]+
				S6Mom[p4,p5]/Chain[$square,p5,{UnderBar[p4]},p6,$angle]*
				(
					Chain[$square,p5,{UnderBar[p3],PlusM[UnderBar[p2],UnderBar[p1]]},p6,$square]+
					Chain[$square,p6,{UnderBar[p3],PlusM[UnderBar[p2],UnderBar[p1]]},p5,$square]
				)+
				S6Mom[p4,p5]^2/Chain[$square,p5,{UnderBar[p4]},p6,$angle]^2*
					Chain[$square,p5,{UnderBar[p3],PlusM[UnderBar[p2],UnderBar[p1]]},p5,$square]
			)
	}
]


A6sSSsPPChannel61[p1_,p2_,p3_,p4_,p5_,p6_]:=-1/S6Mom[p1,p6]Chain[$square,OverHat[p6],{UnderBar[p1]},p5,$angle]/SpinorAngleBracket[p6,p5]*
	Chain[$square,p5,{PlusM[UnderBar[p6],UnderBar[p1]],UnderBar[p4]},p5,$square]/(S6Mom[p4,OverHat[p5]]*S6Mom[p5,p6,p1])*
	(1+(2Mom4DVector[p3][Null][$up] . Mom4DVector[p4][Null][$down]+
		S6Mom[p5,p6,p1]*Chain[$square,p5,{UnderBar[p3],UnderBar[p4]},p5,$square]/
			Chain[$square,p5,{PlusM[UnderBar[p6],UnderBar[p1]],UnderBar[p4]},p5,$square]
	)/S6Mom[p2,p3])/.
	{
		Chain[$square,OverHat[p6],{UnderBar[p1]},p5,$angle]->
			Chain[$square,p6,{UnderBar[p1]},p5,$angle]-S6Mom[p1,p6]S6Mom[p1,p5]/Chain[$square,p5,{UnderBar[p1]},p6,$angle],
		S6Mom[p4,OverHat[p5]]->
			S6Mom[p4,p5]+S6Mom[p1,p6]Chain[$square,p5,{UnderBar[p4]},p6,$angle]/Chain[$square,p5,{UnderBar[p1]},p6,$angle]
	}


A6sSSsPP[p1_,p2_,p3_,p4_,p5_,p6_] := Module[{k1,k2,k3,k4,k5,k6},
	A6sSSsPPChannel61[k1,k2,k3,k4,k5,k6]+A6sSSsPPChannel45[k1,k2,k3,k4,k5,k6]/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p1,p2,p3,p4,p5,p6}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 6: S+SS'+S'*)


A6sSPSsPChannel56[p1_,p2_,p3_,p4_,p5_,p6_]:=Chain[$square,p6,{UnderBar[p5],UnderBar[p1]},p6,$square]/(S6Mom[p1,p6]S6Mom[p5,p6])*
	Chain[$square,p3,{UnderBar[p4],UnderBar[p2]},p3,$square]/(S6Mom[p2,p3]S6Mom[p3,p4])*
	(1+(S6Mom[p2,OverHat[p1]]+
		S6Mom[p3,p4]*Chain[$square,p3,{UnderBar[OverHat[p1]],UnderBar[p2]},p3,$square]/Chain[$square,p3,{UnderBar[p4],UnderBar[p2]},p3,$square])/S6Mom[p5,p6,p1])/.
	{
		S6Mom[p2,OverHat[p1]]->
			2Mom4DVector[p1][Null][$up] . Mom4DVector[p2][Null][$down]+
				S6Mom[p5,p6]Chain[$square,p6,{UnderBar[p2],UnderBar[p1]},p6,$square]/Chain[$square,p6,{UnderBar[p5],UnderBar[p1]},p6,$square],
		Chain[$square,p3,{UnderBar[OverHat[p1]],UnderBar[p2]},p3,$square]->
			Chain[$square,p3,{UnderBar[p1],UnderBar[p2]},p3,$square]+S6Mom[p5,p6]SpinorSquareBracket[p3,p6]*
				Chain[$square,p3,{UnderBar[p2],UnderBar[p1]},p6,$square]/Chain[$square,p6,{UnderBar[p5],UnderBar[p1]},p6,$square]
	}


A6sSPSsP[p1_,p2_,p3_,p4_,p5_,p6_] := Module[{k1,k2,k3,k4,k5,k6},
	A6sSPSsPChannel56[k1,k2,k3,k4,k5,k6]/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p1,p2,p3,p4,p5,p6}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 6: S++S'SS'*)


A6SPPsSsChannel234[p1_,p2_,p3_,p4_,p5_,p6_]:=
	Extramass[p4]Extramasstilde[p4]/(S6Mom[p2,p3,p4]S6Mom[p3,p4])*
	SpinorSquareBracket[p2,p3]^2/S6Mom[OverHat[p2],p3]/.
	{
		S6Mom[OverHat[p2],p3]->
			S6Mom[p2,p3]-S6Mom[p2,p3,p4]*
				Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]/Chain[$square,p2,{PlusM[UnderBar[p3],UnderBar[p4]],UnderBar[p1]},p2,$square]
	};


A6SPPsSsChannel23[p1_,p2_,p3_,p4_,p5_,p6_]:=
	-Chain[$square,p2,{UnderBar[p1],PlusM[p2,p3],UnderBar[p4],PlusM[UnderBar[p5],UnderBar[p6]],PlusM[p2,p3],UnderBar[p1]},p2,$square]*
	1/(S6Mom[p1,p2,p3]S6Mom[p1,p2]Chain[$angle,p3,{p2,UnderBar[p1]},p3,$angle]S6Mom[OverHat[p2],p3,p4])/.
	{
		S6Mom[OverHat[p2],p3,p4]->S6Mom[p2,p4]+S6Mom[p3,p4]-S6Mom[p2,p3]*
			Chain[$square,p2,{UnderBar[p4],UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]
	};


A6SPPsSs[p1_,p2_,p3_,p4_,p5_,p6_]:=Module[{k1,k2,k3,k4,k5,k6},
	A6SPPsSsChannel234[k1,k2,k3,k4,k5,k6]+A6SPPsSsChannel23[k1,k2,k3,k4,k5,k6]/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p1,p2,p3,p4,p5,p6}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 6: S+S'SS'+*)


A6SPsSsPChannel16[p1_,p2_,p3_,p4_,p5_,p6_]:=
	Extramass[p1]Extramasstilde[p1]/(S6Mom[OverHat[p2],p3]S6Mom[p6,p1,p2]S6Mom[p1,p6]Chain[$angle,p6,{UnderBar[p1]},p2,$square])*
	SpinorSquareBracket[p6,p2]Chain[$square,p2,{UnderBar[p3],PlusM[UnderBar[p1],p6]},p2,$square]/.
	{
		S6Mom[OverHat[p2],p3]->
			S6Mom[p2,p3]+S6Mom[p1,p6]*
				Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]/Chain[$square,p2,{p6,UnderBar[p1]},p2,$square]
	}


A6SPsSsPChannel23[p1_,p2_,p3_,p4_,p5_,p6_]:=
	Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]*Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square]*
	1/(S6Mom[p2,p3]S6Mom[p1,p2]S6Mom[p5,p6]S6Mom[OverHat[p1],p6])/.
	{
		S6Mom[OverHat[p1],p6]->S6Mom[p6,p1]+S6Mom[p2,p3]*
			Chain[$square,p2,{p6,UnderBar[p1]},p2,$square]/Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square],
		Chain[$square,p6,{UnderBar[OverHat[p1]],UnderBar[p5]},p6,$square] ->
			Chain[$square,p6,{UnderBar[p1],UnderBar[p5]},p6,$square]+S6Mom[p2,p3]*SpinorSquareBracket[p6,p2]*
			Chain[$square,p2,{UnderBar[p1],UnderBar[p5]},p6,$square]/Chain[$square,p2,{UnderBar[p1],UnderBar[p3]},p2,$square]
	}


A6SPsSsP[p1_,p2_,p3_,p4_,p5_,p6_]:=
Module[{k1,k2,k3,k4,k5,k6},
	A6SPsSsPChannel16[k1,k2,k3,k4,k5,k6]+A6SPsSsPChannel23[k1,k2,k3,k4,k5,k6]/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6},{p1,p2,p3,p4,p5,p6}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 6: S+S'S+S'*)


A6SPsSPs[p1_,p2_,p3_,p4_,p5_,p6_]:=
	Chain[$square,p2,{UnderBar[p3],UnderBar[p1]},p2,$square]Chain[$square,p5,{UnderBar[p6],UnderBar[p4]},p5,$square]/
	(S6Mom[p2,p3]S6Mom[p1,p2]S6Mom[p5,p6]S6Mom[p4,p5]);


(* ::Subsection:: *)
(*Amplitude 7: SGGGGGS*)


A7SPPPPPSBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","+","+","S"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATree["-","+","+"][pm@OverHat[K],OverHat[k2],k3],
	{k2,k3},k1,k2,{k2,k3,k4,k5,k6},K]
	
}/.{PlusM[x_]:>x}//.{Extramass[k6]->-Extramass[k1],Extramasstilde[k6]->-Extramasstilde[k1]}//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4,k5,k6}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}//.
{
	Chain[a__,{b___,PlusM[c__],PlusM[c__,d__],e___},f__]:>S6Mom[c]Chain[a,{b,e},f]+Chain[a,{b,PlusM[c],d,e},f],
	Chain[a__,{b___,c_,c_,e___},f__]:>(Extramass[c]Extramasstilde[c]/.UnderBar->Identity)Chain[a,{b,e},f],
	Chain[a__,{b___,PlusM[c_,d__]},c_,f_]:>Chain[a,{b,PlusM[d]},c,f]
}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]//SortChains//FullSimplify
]


(* ::Subsection:: *)
(*Amplitude 7: SGGGSGG*)


A7SPPPSPPBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","S","+","+"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATree["-","+","+"][pm@OverHat[K],OverHat[k2],k3],
	{k2,k3},k1,k2,{k2,k3,k4,k6,k7},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","S"][OverHat[K],k6,k7,OverHat[k1]](1/S6Mom[k6,k7,k1])ATree["S","+","+","+","S"][pm@OverHat[K],OverHat[k2],k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k3,k4,k6,k7},K]/.{Extramass[K]->Extramass[k5],Extramasstilde[K]->Extramasstilde[k5]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATreeSPSScalarGluonShift[OverHat[k1],k7,-1,k1,k2](1/S6Mom[k7,k1])ATree["S","+","+","+","S","+"][pm@OverHat[K],OverHat[k2],k3,k4,k5,k6],
	{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4,k6,k7},K]
	
}/.{PlusM[x_]:>x}//.{Extramass[k5]->-Extramass[k1],Extramasstilde[k5]->-Extramasstilde[k1]}//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4,k6,k7}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}//.
{
	Chain[a__,{b___,PlusM[c__],PlusM[c__,d__],e___},f__]:>S6Mom[c]Chain[a,{b,e},f]+Chain[a,{b,PlusM[c],d,e},f],
	Chain[a__,{b___,c_,c_,e___},f__]:>(Extramass[c]Extramasstilde[c]/.UnderBar->Identity)Chain[a,{b,e},f],
	Chain[a__,{b___,PlusM[c_,d__]},c_,f_]:>Chain[a,{b,PlusM[d]},c,f],
	Chain[a_,b_,{PlusM[b_,d__],c___},e_,f_]:>Chain[a,b,{PlusM[d],c},e,f]
}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]//SortChains
]


(* ::Subsection:: *)
(*Amplitude 7: SSGGS'S'G*)


(* ::Subsubsection:: *)
(*BCFW Scalar Gluon*)


A7SSPSSPPBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{
	ApplyBcfwScalarGluonShiftRules[
		ATreeSMSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","+","+","+","+","S"][k4,k5,k6,pm@OverHat[K],OverHat[k2],k3],
	{k2,k3,k4,k5,k6},k1,k2,{k2,k5,k6},K],

	ApplyBcfwScalarGluonShiftRules[
		ATreeSPSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","+","+","-","+","S"][k4,k5,k6,pm@OverHat[K],OverHat[k2],k3],
	{k2,k3,k4,k5,k6},k1,k2,{k2,k5,k6},K],		
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","S"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k2,k3,k4,k5])ATree["S","+","-","+","S"][k4,k5,pm@OverHat@K,OverHat@k2,k3],
	{k2,k3,k4,k5},k1,k2,{k2,k5,k6},K],	
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","-","+","S"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][k4,k5,pm@OverHat@K,OverHat@k2,k3],
	{k2,k3,k4,k5},k1,k2,{k2,k5,k6},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","+","S"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","-","+","S"][k4,pm@OverHat@K,OverHat@k2,k3],
	{k2,k3,k4},k1,k2,{k2,k5,k6},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","-","+","+","S"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][k4,pm@OverHat@K,OverHat@k2,k3],
	{k2,k3,k4},k1,k2,{k2,k5,k6},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","+","s","s"][OverHat[K],k4,k5,k6,k7,OverHat[k1]](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k5,k6},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]}
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k7,k1},{k3,k4}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k5,k6}]//SortChains//MapAt[simplifyBCFW,{{1},{2},{3},{4},{5},{6}}];
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p2,p3,p4,p5,p6,p7,p1}}]
]


(* ::Subsection:: *)
(*Amplitude 7: SSGGGS'S'*)


(* ::Subsubsection:: *)
(*BCFW Scalar Gluon*)


A7SSPPPSSBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{
	ApplyBcfwScalarGluonShiftRules[
			ATree["S","S","+","+","s","s"][k7,OverHat[k1],OverHat[K],k4,k5,k6](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
		{k2,k3},k1,k2,{k2,k3,k4},K],
		
		ApplyBcfwScalarGluonShiftRules[
			ATree["S","S","s","s"][k7,OverHat[k1],OverHat[K],k6](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
		{k2,k3,k4,k5},k1,k2,{k2,k3,k4},K],
		
		ApplyBcfwScalarGluonShiftRules[
			ATreeSMSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","+","+","+","+","S"][k6,pm@OverHat[K],OverHat[k2],k3,k4,k5],
		{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4},K],
		
		ApplyBcfwScalarGluonShiftRules[
			ATreeSPSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","-","+","+","+","S"][k6,pm@OverHat[K],OverHat[k2],k3,k4,k5],
		{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4},K]
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k7,k1},{k5,k6}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4}]//SortChains;
res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.
    MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p2,p3,p4,p5,p6,p7,p1}}]//
		MapAt[simplifyBCFW,{{2},{3},{4}}]//Apply[Plus]
];

A7SSPPPssContactBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
	res=-{
		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","+","+","s","s"][k7,OverHat[k1],OverHat[K],k4,k5,k6](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
			{k2,k3},k1,k2,{k2,k3,k4},K],

		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","s","s"][k7,OverHat[k1],OverHat[K],k6](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
			{k2,k3,k4,k5},k1,k2,{k2,k3,k4},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k7,k1},{k5,k6}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4}]//SortChains;
	res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.
			MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p2,p3,p4,p5,p6,p7,p1}}]//Apply[Plus]
];

A7SSPPPssGluonBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
	res=-{
		ApplyBcfwScalarGluonShiftRules[
			ATreeGluon["S","S","+","+","s","s"][k7,OverHat[k1],OverHat[K],k4,k5,k6](1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
			{k2,k3},k1,k2,{k2,k3,k4},K],

		ApplyBcfwScalarGluonShiftRules[
			ATreeGluon["S","S","s","s"][k7,OverHat[k1],OverHat[K],k6](1/S6Mom[k2,k3,k4,k5])ATree["S","+","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
			{k2,k3,k4,k5},k1,k2,{k2,k3,k4},K],

		ApplyBcfwScalarGluonShiftRules[
			ATreeSMSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","+","+","+","+","S"][k6,pm@OverHat[K],OverHat[k2],k3,k4,k5],
			{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4},K],

		ApplyBcfwScalarGluonShiftRules[
			ATreeSPSScalarGluonShift[k7,OverHat[K],-1,k1,k2](1/S6Mom[k2,k3,k4,k5,k6])ATree["S","-","+","+","+","S"][k6,pm@OverHat[K],OverHat[k2],k3,k4,k5],
			{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k7,k1},{k5,k6}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4}]//SortChains;
	res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.
			MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p2,p3,p4,p5,p6,p7,p1}}]//Apply[Plus]
];


(* ::Subsection:: *)
(*Amplitude 7:  SPSsPPs*)


A7SPSsPPsBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","s","+","+","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k5,k6},K]

}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k3},{k4,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k5,k6}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
]


(* ::Subsection:: *)
(*Amplitude 7: SPSPsPs*)


A7SPSPsPsBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","s","+","s"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","+","S","+"][pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4},k1,k2,{k2,k4,k6},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","s","+","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k4,k6},K]
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k3},{k5,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k4,k6}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]//simplifyBCFW
];


(* ::Subsection:: *)
(*Amplitude 7: SPSPssP*)


A7SPSPssPBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","s","s","+"][k5,k6,k7,OverHat[k1],OverHat[K],k4](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k4,k7},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","s","s"][k5,k6,k7,OverHat[k1],OverHat[K]](1/S6Mom[k2,k3,k4])ATree["S","+","S","+"][pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4},k1,k2,{k2,k4,k7},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATreeSPSScalarGluonShift[OverHat[k1],k7,-1,k1,k2](1/S6Mom[k1,k7])ATree["S","S","s","+","s","+"][k5,k6,pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4,k5,k6},k1,k2,{k2,k4,k7},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k3},{k5,k6}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k4,k7}]//SortChains//MapAt[simplifyBCFW,{2}];
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
]


(* ::Subsection:: *)
(*Amplitude 7: SPSPPss*)


A7SPSPPssBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","+","s","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k4,k5},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","s","s"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","+","S","+"][pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4},k1,k2,{k2,k4,k5},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","s","s"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k1,k6,k7])ATree["S","+","S","+","+"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k4,k5},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k3},{k6,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k4,k5}]//SortChains;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0,Extramass[k4]->0,
	Extramass[k5]->0,Extramasstilde[k4]->0, Extramasstilde[k5]->0}/.MapThread[Rule,
	{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
];


(* ::Subsection::Closed:: *)
(*Amplitude 7: SPPSPss*)


A7SPPSPssBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","S","+","s","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATree["-","+","+"][pm@OverHat[K],OverHat[k2],k3],
	{k2,k3},k1,k2,{k2,k3,k5},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","+","s","s"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4},k1,k2,{k2,k3,k5},K]//.{Extramass[K]->Extramass[k4],Extramasstilde[K]->Extramasstilde[k4]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","s","s"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k1,k6,k7])ATree["S","+","+","S","+"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k3,k5},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k4},{k6,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k5}]//SortChains;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
]


(* ::Subsection:: *)
(*Amplitude 7: SPPPSss*)


A7SPPPPSssBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","S","s","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATree["-","+","+"][pm@OverHat[K],OverHat[k2],k3],
	{k2,k3},k1,k2,{k2,k3,k4},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","S","s","s"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k1,k6,k7])ATree["S","+","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k3,k4},K]/.{Extramass[K]->Extramass[k5],Extramasstilde[K]->Extramasstilde[k5]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k5},{k6,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]//simplifyBCFW
]


(* ::Subsection:: *)
(*Amplitude 7: SPsPSPs*)


A7SPsPSPsBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","s","+","S","+","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k4,k6},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","s","S","+","s"][OverHat[k1],OverHat[K],k5,k6,k7](1/S6Mom[k2,k3,k4])ATree["S","+","S","+"][pm@OverHat@K,OverHat@k2,k3,k4],
	{k2,k3,k4},k1,k2,{k2,k4,k6},K]/.{Extramass[K]->Extramass[k3],Extramasstilde[K]->Extramasstilde[k3]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k5},{k3,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k4,k6}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]//simplifyBCFW
]


(* ::Subsection::Closed:: *)
(*Amplitude 7: SPsSsPP*)


A7SPsSsPPBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","s","S","s","+","+"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k6,k7},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","S"][OverHat[K],k6,k7,OverHat[k1]](1/S6Mom[k6,k7,k1])ATree["S","+","s","S","s"][pm@OverHat@K,OverHat@k2,k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k6,k7},K]/.{Extramass[K]->Extramass[k4],Extramasstilde[K]->Extramasstilde[k4]},
	
	ApplyBcfwScalarGluonShiftRules[
		ATreeSPSScalarGluonShift[OverHat[k1],k7,-1,k1,k2](1/S6Mom[k7,k1])ATree["S","+","s","S","s","+"][pm@OverHat[K],OverHat[k2],k3,k4,k5,k6],
	{k2,k3,k4,k5,k6},k1,k2,{k2,k6,k7},K]/.{Extramass[K]->Extramass[k4],Extramasstilde[K]->Extramasstilde[k4]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k4},{k3,k5}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k6,k7}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 7: SPsSPPs*)


A7SPsSPPsBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","s","S","+","+","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATreeSPSScalarGluonShift[k3,OverHat[k2],-1,k1,k2],
	{k2,k3},k1,k2,{k2,k5,k6},K]
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k4},{k3,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k5,k6}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
]


(* ::Subsection::Closed:: *)
(*Amplitude 7: SPPPsSs*)


A7SPPPsSsBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=Module[{k1,k2,k3,k4,k5,k6,k7,K,res},
res=-{

	ApplyBcfwScalarGluonShiftRules[
		ATree["S","+","+","s","S","s"][OverHat[k1],OverHat[K],k4,k5,k6,k7](1/S6Mom[k2,k3])ATree["-","+","+"][pm@OverHat@K,OverHat[k2],k3],
	{k2,k3},k1,k2,{k2,k3,k4},K],
	
	ApplyBcfwScalarGluonShiftRules[
		ATree["S","s","S","s"][OverHat[k1],OverHat[K],k6,k7](1/S6Mom[k6,k7,k1])ATree["S","+","+","+","S"][pm@OverHat[K],OverHat[k2],k3,k4,k5],
	{k2,k3,k4,k5},k1,k2,{k2,k3,k4},K]/.{Extramass[K]->Extramass[k5],Extramasstilde[K]->Extramasstilde[k5]}
	
}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k1,k6},{k5,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7},{k2,k3,k4}]//SortChains//simplifyBCFW;
Plus@@res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7},{p1,p2,p3,p4,p5,p6,p7}}]
];

(* ::Subsection::Closed:: *)
(*Amplitude 8: SSPPPPss*)

A8SSPPPPssContactBCFWSG[p1_,p2_,p3_,p4_,p5_,p6_,p7_,p8_]:=Module[{k1,k2,k3,k4,k5,k6,k7,k8,K,res},
	res=-{
		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","+","+","+","s","s"][k8,OverHat[k1],OverHat[K],k4,k5,k6,k7]
					(1/S6Mom[k2,k3])ATree["+","+","-"][OverHat@k2,k3,pm@OverHat@K],
			{k2,k3},k1,k2,{k2,k3,k4,k5},K],

		ApplyBcfwScalarGluonShiftRules[
			ATreeContact["S","S","s","s"][k8,OverHat[k1],OverHat[K],k7]
					(1/S6Mom[k2,k3,k4,k5,k6])
					ATree["S","+","+","+","+","S"][pm@OverHat@K,OverHat@k2,k3,k4,k5,k6],
			{k2,k3,k4,k5,k6},k1,k2,{k2,k3,k4,k5},K]
	}/.{PlusM[x_]:>x}/.AlignExtramassRules[{k8,k1},{k6,k7}]//ReduceMomenta[{k1,k2,k3,k4,k5,k6,k7,k8}, {k2,k3,k4,k5}]//SortChains;
	res//.{SpinorSquareBracket[k2,OverHat[K]]->1, S6Mom[_]:>0}/.
			MapThread[Rule,{{k1,k2,k3,k4,k5,k6,k7,k8},{p2,p3,p4,p5,p6,p7,p8,p1}}]//Apply[Plus]
];


(* ::Section::Closed:: *)
(*Tree Amplitude Definitions Gravity*)


(* ::Subsection::Closed:: *)
(*Amplitude 3: SgS*)


MSppS[k1_,k2_,k3_]:=-Chain[$square,k2,{UnderBar[k1]},$AmplitudesReference,$angle]^2/SpinorAngleBracket[k2,$AmplitudesReference]^2;


MSmmS[k1_,k2_,k3_]:=-Chain[$angle,k2,{UnderBar[k1]},$AmplitudesReference,$square]^2/SpinorSquareBracket[k2,$AmplitudesReference]^2;


(* ::Subsection:: *)
(*Amplitude 4: SggS*)


MSppPPS[k1_,k2_,k3_,k4_]:=-(Extramass[k1]Extramasstilde[k1]SpinorSquareBracket[k2,k3]/SpinorAngleBracket[k2,k3])^2(1/S6Mom[k1,k3]+1/S6Mom[k1,k2]);


(* ::Subsection:: *)
(*Amplitude 4: gggg*)


MmmPPppMM[k1_,k2_,k3_,k4_]:=SpinorSquareBracket[k2,k1]/SpinorAngleBracket[k1,k2]*SpinorAngleBracket[k4,k1]^7/
	(SpinorAngleBracket[k2,k3]SpinorAngleBracket[k3,k4]^2*SpinorAngleBracket[k2,k4]SpinorAngleBracket[k1,k3])


(* ::Subsection:: *)
(*Amplitude 4: SSss*)


MSSss[k1_,k2_,k3_,k4_]:=
	AmplitudeCoefficient["MSSssS12"]S6Mom[k1,k2]+
	AmplitudeCoefficient["MSSssS23"]S6Mom[k2,k3]+
	AmplitudeCoefficient["MSSssMu1"]Extramass[k2]Extramasstilde[k2]+
	AmplitudeCoefficient["MSSssMu2"]Extramass[k3]Extramasstilde[k3]+
	AmplitudeCoefficient["MSSssdilaton"]*2Extramass[k2]Extramasstilde[k2]Extramass[k3]Extramasstilde[k3]/S6Mom[k1,k2]+
	(2*Mom4DVector[k1][Null][$up] . Mom4DVector[k3][Null][$down])*(2*Mom4DVector[k2][Null][$up] . Mom4DVector[k3][Null][$down])/S6Mom[k1,k2];


(* ::Section::Closed:: *)
(*Amplitude Infinity: SGG...GGS*)


AnSnGS[plus__][p1_,pGl__,pn_] /; {plus} === ConstantArray["+",Length@{pGl}] := Module[{res,y,nGl = Length@{pGl},piProd},

piProd = piProduct[p1,pGl,pn];
y = Product[ S6Mom[p1,Sequence@@{pGl}[[1;;j]]],{j,1,Length@{pGl}-1}];

(* For some reason we need a factor of -1 for every gluon. Can't be bothered to search for it*)
(* I am guessing convention differences for polarization vectors... *)
(-1)^(nGl)*Extramass[p1]*Extramasstilde[p1]/(PTFOpen[pGl])*Plus@@(piProd/y)
];


piProduct[mom__]:=Module[{dummies, dummyRules,prods},

dummies = Table[Unique["p"],{i,1,Length@{mom}}];
dummyRules = MapThread[#1->#2&,{dummies,{mom}}];

prods = Reverse@Table[{S6Mom[Sequence@@dummies[[1;;k-1]]],{dummies[[k]],plus@@dummies[[1;;k-1]]}},{k,3,Length@{mom}-2}];
prods = Dot@@(Replace[prods,{a_,{b_,c_}}:>a*IdentityMatrix+b . c,{1}]);

prods = Distribute[prods]//. {Dot[a___, d_S6Mom b_, c___]:> d Dot[a, b, c],Dot[a___,IdentityMatrix,b___]:>Dot[a,b]}/.IdentityMatrix->1;

prods = First@Replace[Verbatim[List@@prods],{HoldPattern[a_*Dot[b___]] :>{a,{b}},b_/;FreeQ[{b},Dot]:>{b,{}}, HoldPattern[Dot[b__]]:>{1,{b}}},{2}];

(*Print[#1*Chain[$square,dummies[[-2]],UnderBar/@#2,dummies[[2]],$square]&@@@prods /. dummyRules];*)
prods = prods //.{{a___,{pref_,{b1___,plus[b2__],b3___}},c___}:>{a,Sequence@@Map[{pref,{b1,#,b3}}&,{b2}],c}};
prods = #1*Chain[$square,dummies[[-2]],UnderBar/@#2,dummies[[2]],$square]&@@@prods;
prods /. dummyRules
];


(* ::Section:: *)
(*One-loop Amplitude Definitions Yang-Mills*)


(* ::Subsection::Closed:: *)
(*One-Loop All-Plus O(\[Epsilon]^0) (n-Point)*)


A1LAllPlusnPt[n_] := -1/3*1/PTF[n]Sum[TrM[Sequence@@config],{config,Subsets[Range@n,{4}]}];
A1LAllPlusnPt[moms__] /; Length@{moms}>3 := -1/3*1/PTF[moms]Sum[TrM[Sequence@@config],{config,Subsets[{moms},{4}]}]


(* ::Subsection::Closed:: *)
(*One-Loop One-Minus O(\[Epsilon]^0) *)


(* ::Text:: *)
(*Negative helicity is always on the first argument.*)
(*Taken from hep-ph/0505055*)
(*For general form look at hep-ph/9312276 (though the form is pretty complicated.)*)


(* ::Subsubsection:: *)
(*4-Point*)


A1LOneMinus4Pt[p1_,p2_,p3_,p4_] := 1/3*SpinorAngleBracket[p2,p4]SpinorSquareBracket[p2,p4]^3/(SpinorSquareBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]SpinorSquareBracket[p4,p1]);
A1LOneMinus4Pt[] := A1LOneMinus4Pt[Sequence@@Range@4];


(* ::Subsubsection::Closed:: *)
(*5-Point*)


A1LOneMinus5Pt[p1_,p2_,p3_,p4_,p5_] := 1/3*1/SpinorAngleBracket[p3,p4]^2*
										(
											-SpinorSquareBracket[p2,p5]^3/(SpinorSquareBracket[p1,p2]SpinorSquareBracket[p5,p1])
											+SpinorAngleBracket[p1,p4]^3*SpinorSquareBracket[p4,p5]SpinorAngleBracket[p3,p5]/(SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p4,p5]^2)
											-SpinorAngleBracket[p1,p3]^3*SpinorSquareBracket[p3,p2]SpinorAngleBracket[p4,p2]/(SpinorAngleBracket[p1,p5]SpinorAngleBracket[p5,p4]SpinorAngleBracket[p3,p2]^2)
										);
A1LOneMinus5Pt[] := A1LOneMinus5Pt@@(Range@5);


(* ::Subsubsection::Closed:: *)
(*6-Point*)


A1LOneMinus6Pt[p1_,p2_,p3_,p4_,p5_,p6_] := 1/3*Plus@@{
	(Chain[$angle,p1,{p2},p6,$square]+Chain[$angle,p1,{p3},p6,$square])^3/
		(S6Mom[p1,p2,p3]*SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p4,p5]^2(Chain[$angle,p3,{p1},p6,$square]+Chain[$angle,p3,{p2},p6,$square])),+
	(Chain[$angle,p1,{p3},p2,$square]+Chain[$angle,p1,{p4},p2,$square])^3/
		(S6Mom[p2,p3,p4]*SpinorAngleBracket[p5,p6]SpinorAngleBracket[p6,p1]SpinorAngleBracket[p3,p4]^2(Chain[$angle,p5,{p3},p2,$square]+Chain[$angle,p5,{p4},p2,$square])),+
	SpinorSquareBracket[p2,p6]^3/(SpinorSquareBracket[p1,p2]SpinorSquareBracket[p6,p1]S6Mom[p3,p4,p5])*(
		SpinorSquareBracket[p2,p3]SpinorSquareBracket[p3,p4]/(SpinorAngleBracket[p4,p5](Chain[$angle,p5,{p3},p2,$square]+Chain[$angle,p5,{p4},p2,$square]))-
		SpinorSquareBracket[p4,p5]SpinorSquareBracket[p5,p6]/(SpinorAngleBracket[p3,p4](Chain[$angle,p3,{p1},p6,$square]+Chain[$angle,p3,{p2},p6,$square]))+
		SpinorSquareBracket[p3,p5]/(SpinorAngleBracket[p3,p4]SpinorAngleBracket[p4,p5])),-
	SpinorAngleBracket[p1,p3]^3SpinorSquareBracket[p2,p3]SpinorAngleBracket[p2,p4]/
		(SpinorAngleBracket[p2,p3]^2SpinorAngleBracket[p3,p4]^2SpinorAngleBracket[p4,p5]SpinorAngleBracket[p5,p6]SpinorAngleBracket[p6,p1]),+
	SpinorAngleBracket[p1,p5]^3SpinorSquareBracket[p5,p6]SpinorAngleBracket[p4,p6]/
		(SpinorAngleBracket[p4,p5]^2SpinorAngleBracket[p5,p6]^2SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]),-
	SpinorAngleBracket[p1,p4]^3SpinorAngleBracket[p3,p5](Chain[$angle,p1,{p2},p4,$square]+Chain[$angle,p1,{p3},p4,$square])/
		(SpinorAngleBracket[p1,p2]SpinorAngleBracket[p2,p3]SpinorAngleBracket[p3,p4]^2SpinorAngleBracket[p4,p5]^2SpinorAngleBracket[p5,p6]SpinorAngleBracket[p6,p1])
	};
	
A1LOneMinus6Pt[] := A1LOneMinus6Pt@@(Range@6);


(* ::Subsection:: *)
(*One-Loop Scalar Pair Amplitudes (Rational only)*)


(* ::Subsubsection:: *)
(*4-pt: SSPP*)


A1LSSPP[p1_,p2_,p3_,p4_]:=1/3*1/SpinorAngleBracket[p3,p4]^2(2*S6Mom[p1,p2]+S6Mom[p2,p3])


(* ::Subsubsection:: *)
(*5-pt: SSPPP*)


A1LSSPPP[k1_,k2_,k3_,k4_,k5_]:=-((Chain[$square, k3, {UnderBar[k2], PlusM[k3, k4]}, k5,
	$square]*Chain[$square, k3, {UnderBar[k2], UnderBar[k1]}, k3, $square]*SpinorSquareBracket[k3, k4])/
   (Chain[$angle, k5, {k4}, k3, $square]*S6Mom[k1, k2]*(Chain[$angle, k4, {UnderBar[k2]}, k3, $square]*S6Mom[k1, k2] - Chain[$square, k3, {UnderBar[k1], UnderBar[k2]}, k3, $square]*
      SpinorAngleBracket[k3, k4])) + (SpinorSquareBracket[k3, k4]*(-(Chain[$square, k3, {k4, UnderBar[k2]}, k3, $square]*(2*S6Mom[k1, k2] + S6Mom[k1, k5])) + 
     2*Chain[$square, k3, {UnderBar[k1], UnderBar[k2]}, k3, $square]*SpinorAngleBracket[k3, k4]*SpinorSquareBracket[k3, k4]))/
   (Chain[$angle, k5, {k4}, k3, $square]^2*S6Mom[k2, k3]*SpinorAngleBracket[k3, k4]) - 
  (Extramass[k2]*Extramasstilde[k2]*(Chain[$angle, k5, {UnderBar[k2]}, k3, $square]*S6Mom[k1, k2] - Chain[$square, k3, {UnderBar[k1], UnderBar[k2]}, k3, $square]*SpinorAngleBracket[k3, k5])*
    SpinorSquareBracket[k3, k5]^3)/(Chain[$square, k3, {UnderBar[k2], PlusM[k3, k4]}, k5, $square]*S6Mom[k1, k2]*(Chain[$angle, k4, {UnderBar[k2]}, k3, $square]*S6Mom[k1, k2] - 
     Chain[$square, k3, {UnderBar[k1], UnderBar[k2]}, k3, $square]*SpinorAngleBracket[k3, k4])*SpinorAngleBracket[k4, k5]) + 
  (Chain[$angle, k4, {UnderBar[k1]}, k5, $square]*Chain[$angle, k5, {k4}, k3, $square]*Chain[$square, k3, {UnderBar[k2], PlusM[k3, k4]}, k5, $square]^2*
     (Extramass[k2]*Extramasstilde[k2] - S6Mom[k1, k5])*SpinorAngleBracket[k3, k4]*SpinorSquareBracket[k3, k4]^2 + Chain[$angle, k4, {UnderBar[k1]}, k5, $square]^2*
     Chain[$square, k3, {UnderBar[k2], PlusM[k3, k4]}, k5, $square]*SpinorSquareBracket[k3, k4]^2*(Chain[$angle, k5, {k4}, k3, $square]*S6Mom[k1, k5]*S6Mom[k2, k3] - 
      Chain[$angle, k5, {PlusM[k5, UnderBar[k1]], PlusM[k3, k4], UnderBar[k2]}, k3, $square]*SpinorAngleBracket[k3, k4]*SpinorSquareBracket[k3, k4]) + 
    Chain[$angle, k5, {k4}, k3, $square]*Chain[$square, k3, {k4, UnderBar[k2]}, k3, $square]*Extramass[k2]*Extramasstilde[k2]*S6Mom[k1, k5]^2*SpinorAngleBracket[k3, k4]*
     SpinorSquareBracket[k3, k5]^2*SpinorSquareBracket[k4, k5])/(Chain[$angle, k5, {k4}, k3, $square]*Chain[$square, k3, {PlusM[UnderBar[k1], UnderBar[k2]], UnderBar[k1]}, k5, $square]*
    Chain[$square, k3, {UnderBar[k2], PlusM[k3, k4]}, k5, $square]*S6Mom[k1, k5]^2*SpinorAngleBracket[k3, k4]^2*SpinorAngleBracket[k4, k5]*SpinorSquareBracket[k3, k4]))/3


(* ::Section::Closed:: *)
(*Two-loop All-Plus Definitions Yang-Mills*)


(* ::Subsection::Closed:: *)
(*Two-Loop Divergent Part (n-Point)*)


TriInt[i_,j_]:=(-\[Mu]/S[i,j])^\[Epsilon];


A2LnPtDiv[n_] := -A1LAllPlusnPt[n]*1/\[Epsilon]^2*Sum[TriInt[Sequence@@ij],{ij,Partition[Range@n,2,1,{1}]}];
A2LnPtDiv[moms__] /; Length@{moms}>3  :=-A1LAllPlusnPt[moms]*1/\[Epsilon]^2*Sum[TriInt[Sequence@@ij],{ij,Partition[{moms},2,1,{1}]}];


A2LnPtDivPartial[n_Integer,i_Integer] /; 1<=i<=n := -A1LAllPlusnPt[n]*1/\[Epsilon]^2*TriInt[Sequence@@(Partition[Range@n,2,1,{1}][[i]])];
A2LnPtDivPartial[moms__,i_Integer] /; (1<=i<=Length@moms)&&Length@{moms}>3  :=-A1LAllPlusnPt[moms]*1/\[Epsilon]^2*TriInt[Sequence@@(Partition[{moms},2,1,{1}][[i]])];


(* ::Subsection:: *)
(*Two-Loop Finite PolyLog Part*)


(* ::Subsubsection::Closed:: *)
(*Two Mass Box Integral Truncated*)


F2ME[s_,t_,m2_,m4_]:=PolyLog[2,1-m2/s]+PolyLog[2,1-m2/t]+PolyLog[2,1-m4/s]+PolyLog[2,1-m4/t]-PolyLog[2,1-m2*m4/(s*t)]+1/2*Log[s/t]^2


(* ::Subsubsection:: *)
(*Finite PolyLog Part*)


K2K4ST[n_Integer,i_Integer,r_Integer,moms_List]:=Module[{ni=moms[[CyclicPermutations[n][[i]]]],K2,K4,S,T},
	K2=ni[[1;;r]];
	T=ni[[1;;r+1]];
	K4=ni[[r+2;;n-1]];
	S=Prepend[K2,ni[[n]]];
	{K2,K4,S,T}
];

K2K4ST[n_Integer,i_Integer,r_Integer]:=K2K4ST[n,i,r,Range@n];


(* i runs from 1 to n, while r runs from 1 to (n-4) *)
cri[n_Integer,i_Integer,r_Integer,moms_List] /; Length@moms == n := Module[{p4,k41,k42,ab,abc,abcd,k2,k4,s,t,res=0,ni=moms[[CyclicPermutations[n][[i]]]]},
	{k2,k4,s,t}=K2K4ST[n,i,r,moms];
	res += Sum[TrM[Sequence@@abcd],{abcd,Subsets[k4,{4}]}];
	res -= Sum[Sum[TrM[Sequence@@abc,p4],{abc,Subsets[k4,{3}]}],{p4,k4}];
	res += Sum[Sum[Chain[$angle,Last@ni,{k41,Sequence@@ab,k42},ni[[1+r]],$angle],{k41,k4},{k42,k4}],{ab,Subsets[k4,{2}]}]/SpinorAngleBracket[Last[ni],ni[[1+r]]];

	res
];

cri[n_Integer,i_Integer,r_Integer]:=cri[n,i,r,Range@n];



A2LnPtPoly[moms__] /; Length@{moms} >3 :=Module[{n=Length@{moms},i,r,K2,K4,S,T},
	-1/3*1/PTF[moms]*Sum[
		{K2,K4,S,T}=K2K4ST[n,i,r,{moms}];
		cri[n,i,r,{moms}]*F2ME[S4Mom[Sequence@@S],S4Mom[Sequence@@T],S4Mom[Sequence@@K2],S4Mom[Sequence@@K4]]/.{S4Mom[_]:>0}
	,{i,1,n},{r,1,n-4}]
];

A2LnPtPolyPartial[moms__,i_Integer] /; Length@{moms} > 3 := Module[{n=Length@{moms},r,K2,K4,S,T},
	-1/3*1/PTF[moms]*Sum[
		{K2,K4,S,T}=K2K4ST[n,i,r,{moms}];
		cri[n,i,r,{moms}]*F2ME[S4Mom[Sequence@@S],S4Mom[Sequence@@T],S4Mom[Sequence@@K2],S4Mom[Sequence@@K4]]/.{S4Mom[_]:>0}
	,{r,1,n-4}]
];


(* ::Subsection::Closed:: *)
(*Rational Part*)


(* ::Subsubsection::Closed:: *)
(*4-Point*)



(*
	R243 and R241B Extracted from hep-ph/0201161 (adding an additional factor of 4, which is
	also required to match R241)
*)

R241[a_,b_,c_,d_]:=1/9*1/PTF[a,b,c,d]*(S[a,c]^2+8S[a,b]S[b,c]);
R241[]:=R241[Sequence@@Range@4];

R243[a_,b_,c_,d_]:=
    -2/9 S[a,b]/(S[b,c]S[a,c])*1/(SpinorAngleBracket[a,b]^2SpinorAngleBracket[c,d]^2)*
        (S[a,b]^3+S[b,c]^3+S[a,c]^3+24*S[a,b]S[b,c]S[a,c]);

R241B[a_,b_,c_,d_]:=0;

(* ::Subsubsection::Closed:: *)
(*5-Point Planar*)


(* ::Text:: *)
(*Badger*)


A52LRatLitBadgerPartial[p1_,p2_,p3_,p4_,p5_]:=Module[{intBoxF3,intBoxF3l1l2,intButF3,intButF3tr,intButF3l1l2},
intBoxF3=S[p1,p3]/4;
intBoxF3l1l2=(TrM[p1,p3,p4,p5]-S[p1,p3](6*S[p4,p5]-2*S[p1,p3]-S[p3,p4]-S[p1,p5]))/36;
intButF3=1/4;
intButF3tr=(2*Chain[$square,p1,{p2,p3,p1,p4,p3,p4,p5},p1,$angle]+Chain[$square,p1,{p2,p3,p2,p4,p3,p4,p5},p1,$angle]+
			4*Chain[$square,p1,{p2,p3,p1,p5,p3,p4,p5},p1,$angle]+2*Chain[$square,p1,{p2,p3,p2,p5,p3,p4,p5},p1,$angle])/36;
intButF3l1l2=1/36*(2*S[p1,p4]+S[p2,p4]+4*S[p1,p5]+2*S[p2,p5]);
TrP[p1,p3,p4,p5]/S[p1,p3]*intBoxF3+TrP[p1,p3,p4,p5]/(S[p1,p3]S[p4,p5])*intBoxF3l1l2+
TrP[p1,p2,p4,p5]*intButF3+1/(S[p1,p2]S[p4,p5])intButF3tr+(S[p1,p5]+TrP[p1,p2,p4,p5](S[p1,p2]+S[p4,p5])/(S[p1,p2]S[p4,p5]))intButF3l1l2
]


A52LRatLitBadger[p1_,p2_,p3_,p4_,p5_]:=-(4-2)^2/PTF[{p1,p2,p3,p4,p5}]*Sum[A52LRatLitBadgerPartial[Sequence@@perm],{perm,CyclicPermutations[{p1,p2,p3,p4,p5}]}]


(* ::Text:: *)
(*Dunbar*)

(* Taken from 1603.07514*)
R5aPartial[p1_,p2_,p3_,p4_,p5_]:=2/3*TrP[p4,p5,p1,p2]^2/(S[p1,p2]S[p4,p5]);
R5a[p1_,p2_,p3_,p4_,p5_]:=Sum[R5aPartial[Sequence@@i],{i,CyclicPermutations[{p1,p2,p3,p4,p5}]}];
R5bPartial[p1_,p2_,p3_,p4_,p5_]:=(10/3*S[p1,p2]S[p2,p3]+2/3*S[p1,p2]S[p3,p4]);
R5b[p1_,p2_,p3_,p4_,p5_]:=Sum[R5bPartial[Sequence@@i],{i,CyclicPermutations[{p1,p2,p3,p4,p5}]}];

A52LRatLitDunbar[p1_,p2_,p3_,p4_,p5_]:=1/6*1/PTF[p1,p2,p3,p4,p5](R5a[p1,p2,p3,p4,p5]+R5b[p1,p2,p3,p4,p5]);
A52LRatLitDunbarPartial[p1_,p2_,p3_,p4_,p5_]:=1/6*1/PTF[p1,p2,p3,p4,p5](R5aPartial[p1,p2,p3,p4,p5]+R5bPartial[p1,p2,p3,p4,p5]);

(* Taken from 1911.06547*)
R251Partial[a_,b_,c_,d_,e_]:=
		TrP[d,e,a,b]^2/(S[d,e] S[a,b])+5 S[a,b] S[b,c]+S[a,b] S[c,d];
R251[a_,b_,c_,d_,e_]:=
    1/9 * 1/PTF[a,b,c,d,e] Sum[R251Partial[Sequence@@i],{i,CyclicPermutations[{a,b,c,d, e}]}];



(* ::Subsubsection::Closed:: *)
(*5-Point Non-Planar*)

(*Taken from 1911.06547*)
(* The sign is different, but this is just a typo in the paper. *)
(* I checked the color relations from Edison & Naculich, and they work only *)
(* when including the sign below *)
(* The color relation given in eq.(3.6) in 1911.06547 is also off *)
(* by an overall sign. *)

R253Partial[a_,b_,c_,d_,e_]:=TrM[a,c,d,e]TrM[e,c,b,a]/(S[a,e]S[c,d])+3/2 S[a,b]^2;
R253[a_,b_,c_,d_,e_]:=-(*!!!*)2/3*1/(PTF[a,b]PTF[c,d,e])*
		SumCycliclyPermutedKinematics[R253Partial@@{a,b,c,d,e},{a,b},{c,d,e}];

(* The reason for the opposite sign compared to 1911.06547 is the usage
 of \[Epsilon][p1,p2,p3,p4] there, which is defined as TrP - TrM, so -Tr5 *)
R251B[p1_,p2_,p3_,p4_,p5_]:=-2*Tr5[p1,p2,p3,p4](-1/PTF[p1,p2,p3,p4,p5]+2/PTF[p1,p3,p4,p2,p5]+2/PTF[p1,p4,p3,p2,p5]+2/PTF[p1,p4,p2,p3,p5]);


(* ::Subsubsection::Closed:: *)
(*6-Point Planar*)


(* ::Text:: *)
(*Dunbar*)


G61[p1_,p2_,p3_,p4_,p5_,p6_] := S4Mom[p3,p4]S4Mom[p4,p6]/(S4Mom[p5,p6]S4Mom[p1,p2,p3])*
								(TrM[p6,p1,p2,p5]+TrM[p6,p1,p3,p5])+S4Mom[p1,p3]S4Mom[p3,p4]/(S4Mom[p1,p2]S4Mom[p4,p5,p6])(TrM[p1,p6,p5,p2]+TrM[p1,p6,p4,p2]);
G62[p1_,p2_,p3_,p4_,p5_,p6_] := TrM[p1,p5,p6,p2]^2/(S4Mom[p1,p2]S4Mom[p6,p5])+1/2*TrM[p1,p3,p4,p6]^2/(S4Mom[p1,p6]S4Mom[p3,p4]);
G63[p1_,p2_,p3_,p4_,p5_,p6_] := S4Mom[p4,p6]/S4Mom[p1,p2,p3]TrM[p1,p6,p4,p3];
G64[p1_,p2_,p3_,p4_,p5_,p6_] := S4Mom[p4,p5,p6]/S4Mom[p1,p6]TrM[p1,p2,p5,p6];
G65[p1_,p2_,p3_,p4_,p5_,p6_] := 2*S4Mom[p1,p3]^2+S4Mom[p5,p2]^2+S4Mom[p1,p2](-3*S4Mom[p1,p3]-2*S4Mom[p1,p4]+6*S4Mom[p1,p5]+4*S4Mom[p2,p3]+S4Mom[p2,p4]+2*S4Mom[p2,p5]+
								4*S4Mom[p2,p6]+7*S4Mom[p3,p4]-S4Mom[p3,p5]-S4Mom[p4,p5]+3*S4Mom[p4,p6])+
								S4Mom[p1,p3](2*S4Mom[p1,p4]+3*S4Mom[p1,p5]-2*S4Mom[p2,p4]-S4Mom[p2,p5]+S4Mom[p3,p6]-5/2*S4Mom[p4,p6])+
								3/2*S4Mom[p1,p4]S4Mom[p2,p5]-8*TrM[p2,p3,p4,p5]+5*TrM[p6,p1,p3,p4];

A62LRatLitDunbarPartial[p1_,p2_,p3_,p4_,p5_,p6_]:=G61@@{p1,p2,p3,p4,p5,p6}+G62@@{p1,p2,p3,p4,p5,p6}+G63@@{p1,p2,p3,p4,p5,p6}+
	G64@@{p1,p2,p3,p4,p5,p6}+G65@@{p1,p2,p3,p4,p5,p6};

A62LRatLitDunbar[p1_,p2_,p3_,p4_,p5_,p6_]:=1/9*1/PTF[p1,p2,p3,p4,p5,p6]*Sum[G61@@perm+G62@@perm+G63@@perm+G64@@perm+G65@@perm,{perm,CyclicPermutations[{p1,p2,p3,p4,p5,p6}]}];

R261[p1_,p2_,p3_,p4_,p5_,p6_] := A62LRatLitDunbar[p1,p2,p3,p4,p5,p6];

(* ::Text:: *)
(*Badger*)


(* ::Text:: *)
(*Result from Simons local integrands paper*)


fR[p1_,p2_,p3_,p4_,p5_,p6_] := 2S6Mom[p2,p3]S6Mom[p3,p4]S6Mom[p4,p5]TrP[p1,p2,p5,p6]/(S6Mom[p1,p2]S6Mom[p5,p6]S6Mom[p1,p2,p3])+
								TrP[p1,p2,p3,p6]/(S6Mom[p1,p2]S6Mom[p1,p2,p3])(-4S6Mom[p2,p3]S6Mom[p3,p4]+2TrP[p1,p3,p4,p5]+6TrP[p2,p3,p4,p5])+
								1/S6Mom[p1,p2,p3](-12S6Mom[p3,p4]TrP[p1,p2,p3,p6]-S6Mom[p3,p4]TrP[p1,p2,p5,p6]+
								3(S6Mom[p1,p2]+S6Mom[p3,p4]+S6Mom[p5,p6])TrP[p1,p3,p4,p6]-
								16S6Mom[p1,p2]S6Mom[p1,p6]S6Mom[p3,p4])-
								TrP[p1,p2,p4,p5]^2/(S6Mom[p1,p2]S6Mom[p4,p5])+
								TrP[p1,p2,p5,p6]/(S6Mom[p1,p2]S6Mom[p5,p6])(-2TrP[p1,p3,p4,p5]-2TrP[p4,p3,p4,p5]-2TrP[p1,p2,p5,p6])+
								1/S6Mom[p1,p2](2(S6Mom[p1,p6]-S6Mom[p3,p4]+S6Mom[p4,p5])TrP[p1,p2,p3,p4]+2(S6Mom[p2,p3]-S6Mom[p3,p4]+S6Mom[p4,p5]-S6Mom[p1,p2,p3])TrP[p1,p2,p3,p5]+
									2(S6Mom[p2,p3]+S6Mom[p4,p5]-3S6Mom[p1,p2,p3])TrP[p1,p2,p4,p5])-
								S6Mom[p1,p2]/4(59*S6Mom[p2,p3]-8S6Mom[p3,p4]-56S6Mom[p4,p5])+
								S6Mom[p1,p2,p3]/4(-4S6Mom[p1,p2]-4S6Mom[p2,p3]+39S6Mom[p3,p4]-40S6Mom[p2,p3,p4])+
								9/4*TrP[p1,p2,p3,p4]+35/4*TrP[p1,p2,p3,p5]+15/4*TrP[p1,p2,p4,p5];

A62LRatLitBadger[p1_,p2_,p3_,p4_,p5_,p6_] := -1*4/144*1/PTF[p1,p2,p3,p4,p5,p6]*Sum[fR@@perm+fR@@(Reverse@perm),{perm,CyclicPermutations[{p1,p2,p3,p4,p5,p6}]}];


A62LRatLitBadgerPartial[p1_,p2_,p3_,p4_,p5_,p6_] := -1*4/144*1/PTF[p1,p2,p3,p4,p5,p6]*Function[perm,fR@@perm+fR@@(Reverse@perm)][{p1,p2,p3,p4,p5,p6}];


(* ::Subsubsection::Closed:: *)
(*6-Point Non-Planar*)


(* ::Text:: *)
(*Tr[xx] Tr[xxxx]*)


H163[a_,b_,c_,d_,e_,f_]:=G163[a,b,c,d,e,f]/(PTF[a,b,c,d,e,f]SpinorAngleBracket[c,d])+
	SpinorSquareBracket[c,d]/SpinorAngleBracket[c,d]^2*(SpinorAngleBracket[c,f]SpinorAngleBracket[d,b]Chain[$square,b,{f},d,$angle])/
		(SpinorAngleBracket[a,b]SpinorAngleBracket[a,f]SpinorAngleBracket[b,f]SpinorAngleBracket[d,e]SpinorAngleBracket[e,f]);

G163[a_,b_,c_,d_,e_,f_]:=S[c,e]Chain[$angle,c,{b,f},d,$angle]-S[c,f]Chain[$angle,c,{b,e},d,$angle];
G263[a_,b_,c_,d_,e_,f_]:=1/S6Mom[d,e,f](Chain[$square,d,{PlusM[d,e,f],b},a,$square]Chain[$angle,d,{f,PlusM[d,e,f]},a,$angle]+
	S6Mom[d,e]Chain[$square,f,{c,b,d},f,$angle]+Chain[$square,b,{d,f},e,$square]Chain[$angle,b,{c,PlusM[a,b,c]},e,$angle]);
G363[a_,b_,c_,d_,e_,f_]:=-S6Mom[d,f]Chain[$angle,d,{f,b},c,$angle]Chain[$square,c,{PlusM[a,b,c]},e,$angle]/(SpinorAngleBracket[d,e]S6Mom[d,e,f])-
	S6Mom[d,e]Chain[$angle,f,{d,b},c,$angle]Chain[$square,c,{d},e,$angle]/(SpinorAngleBracket[e,f]S6Mom[d,e,f]);
G463[a_,b_,c_,d_,e_,f_]:=-S6Mom[b,d]S6Mom[d,e]-Chain[$square,a,{b,d,e},a,$angle]+Chain[$square,b,{c,d,e},b,$angle]-Chain[$square,a,{b,d,f},a,$angle]+
	Chain[$square,b,{c,d,f},b,$angle]+Chain[$square,b,{c,e,f},b,$angle]-Chain[$square,b,{d,e,f},b,$angle];
G563[a_,b_,c_,d_,e_,f_]:=-4S[a,c]^2+2S[a,b]S[a,d]-2S[a,c]S[a,d]+2S[a,b]S[a,e]-2S[a,c]S[a,e]+2S[b,d]^2-2S[b,e]^2+2S[b,f]^2+
	-8S[a,c]S[c,d]+4S[b,c]S[c,d]+12S[b,d]S[c,d]+6S[c,d]^2-8S[a,c]S[c,e]+12 S[b,c]S[c,e]+16S[b,d]S[c,e]+
	+4S[b,e]S[c,e]+8S[c,d]S[c,e]+2S[c,e]^2+2S[c,f]^2-8S[a,c]S[d,e]-4S[a,d]S[d,e]-4S[b,c]S[d,e]+4S[c,d]S[d,e]+
	+4S[c,e]S[d,e]-8Chain[$square,a,{b,c,e},a,$angle]-39Chain[$square,a,{b,c,f},a,$angle]-18Chain[$square,a,{b,d,f},a,$angle]+2Chain[$square,a,{b,e,f},a,$angle]+
	-10Chain[$square,a,{c,d,f},a,$angle]-2Chain[$square,a,{c,e,f},a,$angle]-4Chain[$square,a,{d,e,f},a,$angle]+8Chain[$square,b,{c,d,e},b,$angle]-4Chain[$square,b,{c,d,f},b,$angle]+
	-4Chain[$square,b,{c,e,f},b,$angle]-4Chain[$square,b,{d,e,f},b,$angle]-4Chain[$square,c,{d,e,f},c,$angle];

R263[a_,b_,c_,d_,e_,f_]:=1/3*SumCycliclyPermutedKinematics[
	H163[a,b,c,d,e,f]-H163[a,b,c,d,f,e]+
	(G263[a,b,c,d,e,f]+G363[a,b,c,d,e,f]+G463[a,b,c,d,e,f])/(PTF[a,b,c]PTF[d,e,f])+
	1/4*G563[a,b,c,d,e,f]/PTF[a,b,c,d,e,f],{a,b},{c,d,e,f}];


(* ::Text:: *)
(*Tr[xxx]Tr[xxx]*)


R264[a_,b_,c_,d_,e_,f_]:=1/36*SumCycliclyPermutedKinematics[
	1/PTF[a,b,c]1/PTF[d,e,f]*(G164[a,b,c,d,e,f]+G264[a,b,c,d,e,f])+
	12/(PTF[c,d,e,f]SpinorAngleBracket[a,b])*(G364[a,b,c,d,e,f]+G464[a,b,c,d,e,f])
	+
	1/PTF[d,e,f]1/PTF[a,b,c]*(G164[d,e,f,a,b,c]+G264[d,e,f,a,b,c])+
	12/(PTF[f,a,b,c]SpinorAngleBracket[d,e])*(G364[d,e,f,a,b,c]+G464[d,e,f,a,b,c])
,{a,b,c},{d,e,f}];


G164[a_,b_,c_,d_,e_,f_]:=4Chain[$angle,e,{PlusM[a,b,c],a},b,$angle]*Chain[$square,e,{d,PlusM[a,b,c]},b,$square]/S6Mom[a,b,c];
G264[a_,b_,c_,d_,e_,f_]:=S[a,d]^2+106*S[a,b]S[a,d]+102 TrP[a,b,c,d]-4TrP[a,b,d,e]-4TrP[a,d,b,e];
G364[a_,b_,c_,d_,e_,f_]:=1/SpinorAngleBracket[a,b](TrP[b,a,c,d]+TrP[b,a,e,f]);
G464[a_,b_,c_,d_,e_,f_]:=Chain[$square,a,{c,d},b,$square]+Chain[$square,a,{e,f},b,$square];


(* ::Text:: *)
(*Tr[xx] Tr[xx]Tr[xx]*)


(* Color structures: *)
R2622[a_,b_,c_,d_,e_,f_]:=SumCycliclyPermutedKinematics[
	1/(PTF[a,b,c]PTF[d,e,f])(G1622[a,b,c,d,e,f]+G2622[a,b,c,d,e,f])+
	1/(PTF[a,b,e]PTF[f,c,d])(G1622[a,b,e,f,c,d]+G2622[a,b,e,f,c,d])+
	1/(PTF[c,d,a]PTF[b,e,f])(G1622[c,d,a,b,e,f]+G2622[c,d,a,b,e,f])+
	1/(PTF[e,f,a]PTF[b,c,d])(G1622[e,f,a,b,c,d]+G2622[e,f,a,b,c,d])+
	1/(PTF[c,d,e]PTF[f,a,b])(G1622[c,d,e,f,a,b]+G2622[c,d,e,f,a,b])+
	1/(PTF[e,f,c]PTF[d,a,b])(G1622[e,f,c,d,a,b]+G2622[e,f,c,d,a,b])
,{a,b},{c,d},{e,f}];


G1622[a_,b_,c_,d_,e_,f_]:=Chain[$angle,b,{PlusM[a,b,c],f},d,$angle]*Chain[$square,b,{c,PlusM[a,b,c]},d,$square]/S6Mom[a,b,c];
G2622[a_,b_,c_,d_,e_,f_]:=S[a,d]Chain[$square,e,{PlusM[b,c]},e,$angle]-S[a,c]Chain[$square,e,{PlusM[f,a]},e,$angle]-S[a,f]S[a,e]-S[a,e]S[c,d]


(* ::Text:: *)
(*N_c^0 Tr[xxxxxx]*)


R261B[p1_,p2_,p3_,p4_,p5_,p6_]:=
	-2*1/PTF[p1,p2,p3,p4,p5,p6]Sum[-Tr5@@moms,{moms,Subsets[{p1,p2,p3,p4,p5,p6},{4}]}]+
	-4*(
		Tr5[p3,p4,p5,p6]/PTF[p1,p2,p4,p5,p3,p6]+Tr5[p3,p4,p5,p6]/PTF[p1,p2,p5,p3,p4,p6]+Tr5[p3,p4,p5,p6]/PTF[p1,p2,p5,p4,p3,p6]+
		Tr5[p1,p2,p3,p4]/PTF[p1,p3,p4,p2,p5,p6]-Tr5[p1,p2,p3,p6]/PTF[p1,p3,p4,p5,p2,p6]+Tr5[p1,p2,p3,p4]/PTF[p1,p4,p2,p3,p5,p6]-Tr5[p1,p3,p4,p6]/PTF[p1,p4,p2,p5,p3,p6]+
		Tr5[p1,p2,p3,p4]/PTF[p1,p4,p3,p2,p5,p6]+Tr5[p1,p2,p4,p6]/PTF[p1,p4,p3,p5,p2,p6]-Tr5[p1,p3,p4,p6]/PTF[p1,p4,p5,p2,p3,p6]+Tr5[p1,p2,p4,p6]/PTF[p1,p4,p5,p3,p2,p6]+
		-Tr5[p1,p4,p5,p6]/PTF[p1,p5,p2,p3,p4,p6]+Tr5[p1,p3,p5,p6]/PTF[p1,p5,p2,p4,p3,p6]+Tr5[p1,p3,p5,p6]/PTF[p1,p5,p4,p2,p3,p6]-Tr5[p1,p2,p5,p6]/PTF[p1,p5,p4,p3,p2,p6]
	)


(* ::Subsubsection:: *)
(*7-Point*)


(* ::Text:: *)
(*Dunbar*)


G17[a_,b_,c_,d_,e_,f_,g_] := ((SpinorAngleBracket[g,a])/( S6Mom[a,b,c]S6Mom[e,f,g])) (((SpinorAngleBracket[c,d]SpinorSquareBracket[e,g]Chain[$square,d,{K[a,b,c]},e,$angle]Chain[$square,a,{K[a,b,c]},e,$angle] Chain[$square,c,{K[a,b,c]},f,$angle])/( SpinorAngleBracket[e,f]))-((SpinorAngleBracket[d,e]SpinorSquareBracket[c,a] Chain[$square,d,{K[e,f,g]},c,$angle]  Chain[$square,g,{K[e,f,g]},c,$angle] Chain[$square,e,{K[e,f,g]},b,$angle]  )/( SpinorAngleBracket[b,c]))+ (( SpinorAngleBracket[e,f] SpinorAngleBracket[c,d]  SpinorSquareBracket[c,a]SpinorSquareBracket[f,g] Chain[$square,e,{K[e,f,g]},a,$angle] Chain[$square,d,{K[e,f,g]},b,$angle] )/( SpinorAngleBracket[a,b])) - ((SpinorAngleBracket[b,c] SpinorAngleBracket[d,e] SpinorSquareBracket[e,g] SpinorSquareBracket[a,b] Chain[$square,c,{K[a,b,c]},g,$angle] Chain[$square,d,{K[a,b,c]},f,$angle])/( SpinorAngleBracket[f,g])));
G27[a_,b_,c_,d_,e_,f_,g_] := ((1)/( S6Mom[a,b,c]S6Mom[e,f,g])) S6Mom[c,d]S6Mom[d,e]SpinorAngleBracket[g,a]Chain[$square,g,{K[e,f,g],K[a,b,c]},a,$square];
G37[a_,b_,c_,d_,e_,f_,g_] := ((1)/( S6Mom[c,d,e])) ( S6Mom[c,e](((S6Mom[e,f]Chain[$angle,c,{K[a,b],K[f,g,a]},d,$angle])/( SpinorAngleBracket[c,d]))-((S6Mom[b,c]Chain[$angle,e,{K[f,g],K[g,a,b]},d,$angle] )/( SpinorAngleBracket[d,e])))+((SpinorAngleBracket[e,f] SpinorAngleBracket[b,c] SpinorSquareBracket[f,b]Chain[$square,c,{K[c,d,e]},g,$angle]Chain[$square,e,{K[c,d,e]},a,$angle])/( SpinorAngleBracket[g,a]))+((SpinorAngleBracket[b,c] Chain[$square,c,{K[c,d,e]},b,$angle]Chain[$square,e,{K[c,d,e]},a,$angle]Chain[$square,b,{K[f,g]},e,$angle])/( SpinorAngleBracket[a,b]))+(( SpinorAngleBracket[e,f] Chain[$square,e,{K[c,d,e]},f,$angle]Chain[$square,c,{K[c,d,e]},g,$angle]Chain[$square,f,{K[a,b]},c,$angle])/( SpinorAngleBracket[f,g])));
G47[a_,b_,c_,d_,e_,f_,g_] := ((SpinorSquareBracket[g,a])/( SpinorAngleBracket[g,a])) SpinorAngleBracket[g,e] SpinorAngleBracket[a,e]( ((SpinorSquareBracket[d,e])/( SpinorAngleBracket[d,e]))SpinorAngleBracket[d,g] SpinorAngleBracket[d,a]+((SpinorSquareBracket[e,f])/( SpinorAngleBracket[e,f]))SpinorAngleBracket[f,g]SpinorAngleBracket[f,a] );
G57[a_,b_,c_,d_,e_,f_,g_] := ((1)/( S6Mom[c,d,e]))( SpinorSquareBracket[c,e] (SpinorAngleBracket[e,f] SpinorSquareBracket[d,f] Chain[$angle,c,{K[a,b],K[f,g,a]},d,$angle]+SpinorAngleBracket[b,c] SpinorSquareBracket[d,b] Chain[$angle,e,{K[f,g],K[g,a,b]},d,$angle]  )+SpinorAngleBracket[b,c]SpinorAngleBracket[e,f](2 SpinorAngleBracket[g,a] SpinorSquareBracket[c,e]SpinorSquareBracket[f,g] SpinorSquareBracket[a,b]+  SpinorSquareBracket[b,f] Chain[$square,e,{K[a,b],K[f,g]},c,$square]));
G67[a_,b_,c_,d_,e_,f_,g_] := ((1)/(SpinorAngleBracket[g,a])) (Chain[$angle,g,{f,K[b,c]},a,$angle] S6Mom[e,f,g] -Chain[$angle,a,{b,K[e,f]},g,$angle] S6Mom[a,b,c]);
G77[a_,b_,c_,d_,e_,f_,g_] := S6Mom[b,f]^2-2S6Mom[g,a]^2-3S6Mom[d,b]S6Mom[d,f]+4S6Mom[d,a]S6Mom[d,g]-6S6Mom[a,c]S6Mom[e,g]+7(S6Mom[e,b]S6Mom[f,c]+S6Mom[e,a]S6Mom[g,c])+S6Mom[a,b]S6Mom[f,g]+3S6Mom[f,a]S6Mom[g,b]+S6Mom[c,e](S6Mom[c,f]+S6Mom[e,b]-4(S6Mom[a,b]+S6Mom[f,g]+S6Mom[g,a])+5Chain[$square,d,{K[g,a]},d,$angle])+4Chain[$square,e,{b,c,f},e,$angle]-2Chain[$square,f,{g,a,b},f,$angle]+3Chain[$square,g,{b,a,f},g,$angle]+2Chain[$square,g,{c,e,a},g,$angle];

A72LRatLitDunbarPartial[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=(G17@@{p1,p2,p3,p4,p5,p6,p7}+G27@@{p1,p2,p3,p4,p5,p6,p7}+G37@@{p1,p2,p3,p4,p5,p6,p7}+
	G47@@{p1,p2,p3,p4,p5,p6,p7}+G57@@{p1,p2,p3,p4,p5,p6,p7}+G67@@{p1,p2,p3,p4,p5,p6,p7}+G77@@{p1,p2,p3,p4,p5,p6,p7})//.{Chain[a1__,{b1___,K[c1___],d1___},f1__]:>Sum[Chain[a1,{b1,i,d1},f1],{i,{c1}}]};


A72LRatLitDunbar[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=1/9*1/PTF[p1,p2,p3,p4,p5,p6,p7]*Sum[G17@@perm+G27@@perm+G37@@perm+G47@@perm+G57@@perm+G67@@perm+G77@@perm,{perm,CyclicPermutations[{p1,p2,p3,p4,p5,p6,p7}]}]//.{Chain[a1__,{b1___,K[c1___],d1___},f1__]:>Sum[Chain[a1,{b1,i,d1},f1],{i,{c1}}]};

R271[p1_,p2_,p3_,p4_,p5_,p6_,p7_]:=A72LRatLitDunbar[p1,p2,p3,p4,p5,p6,p7];


(* ::Subsubsection:: *)
(*R2n1B*)


R2n1BS1[list_,r_Integer,i_Integer]/; r<=Length@list-4 && i>r :=list[[r+1;;i-1]];
R2n1BS2[list_,i_Integer,j_Integer]/; j>i :=list[[i+1;;j-1]];
R2n1BS3[list_,j_Integer,s_Integer]/; s>j :=list[[j+1;;s-1]];


Srsij[list_,r_Integer,i_Integer,j_Integer,s_Integer]/; r<i<j<s := 
	ShuffleProduct[R2n1BS1[list,r,i],Reverse@R2n1BS2[list,i,j],R2n1BS3[list,j,s]];


alphaSrsij[list_,r_Integer,i_Integer,j_Integer,s_Integer]/; r<i<j<s :=
	Join[list[[;;r]],list[[{j}]],#,list[[{i}]],list[[s;;]]]&/@Srsij[list,r,i,j,s];


R2n1B[momenta_List]:=R2n1B1[momenta]+R2n1B2[momenta];
R2n1B1[momenta_List]:=-2/PTF@@momenta*Sum[-Tr5@@ijkl,{ijkl,Subsets[momenta,{4}]}];
R2n1B2[momenta_List]:=Module[{n=Length@momenta},
	4*Sum[
		epsilon[momenta[[;;r]],momenta[[j]],momenta[[i]],momenta[[s;;]]](-1)^(i-j+1)
				Sum[1/PTF@@alpha,{alpha,alphaSrsij[momenta,r,i,j,s]}],
	{r,1,n-4},{s,r+4,n},{i,r+1,s-2},{j,i+1,s-1}]
];

epsilon[l1_List,p1_,p2_,l2_List]:=epsilon[#1,p1,p2,#2]&@@@Tuples[{l1,l2}]//Apply[Plus];
epsilon[p1_,p2_,p3_,p4_]:=-Tr5[p1,p2,p3,p4];


(* ::Section::Closed:: *)
(*Two-loop All-Plus Definitions Gravity*)


(* ::Subsection:: *)
(*Two-Loop Divergent Part (n-Point)*)


(* ::Subsection:: *)
(*Two-Loop Finite PolyLog Part*)


(* ::Subsection:: *)
(*Rational Part*)


(* ::Subsubsection::Closed:: *)
(*4-Point*)


MTwoLoopRemainderFreiburg[k1_,k2_,k3_,k4_]:=(SpinorSquareBracket[k1,k2]SpinorSquareBracket[k3,k4])^2/(SpinorAngleBracket[k1,k2]SpinorAngleBracket[k3,k4])^2*
	(
		(117617*S*T*U)/21600  + 
		(S*(14*S^2 + 9*S*T + 9*T^2)*Log[-S])/60 + 
		(T*(9*S^2 + 9*S*T + 14*T^2)*Log[-T])/60 + 
		((-14*S^3 - 33*S^2*T - 33*S*T^2 - 14*T^3)*Log[-U])/60 + 30*S*T*U*(cGB[mu]-2*cR3[mu])/.{S->S[k1,k2],T->S[k2,k3],U->S[k1,k3]}
	);


(*MTwoLoopRemainderFreiburg[k1_,k2_,k3_,k4_]:=(SpinorSquareBracket[k1,k2]SpinorSquareBracket[k3,k4])^2/(SpinorAngleBracket[k1,k2]SpinorAngleBracket[k3,k4])^2*
	(
		(-117617*S*T*(S + T))/21600 - (I/60)*Pi*S*(14*S^2 + 9*S*T + 9*T^2) + 
		(S*(14*S^2 + 9*S*T + 9*T^2)*Log[S])/60 + 
		(T*(9*S^2 + 9*S*T + 14*T^2)*Log[-T])/60 + 
		((-14*S^3 - 33*S^2*T - 33*S*T^2 - 14*T^3)*Log[S + T])/60 + 30*S*T*U*(cGB[mu]-2*cR3[mu])/.{S->S[k1,k2],T->S[k2,k3],U->S[k1,k3]}
	);*)


MTwoLoopRationalFreiburg[k1_,k2_,k3_,k4_]:=MTwoLoopRemainderFreiburg[k1,k2,k3,k4]/.{cGB[_]:>0,cR3[_]:>0,Log[_]:>0,Pi->0};


(* ::Chapter:: *)
(*Verifications*)


(* ::Section::Closed:: *)
(*KLT Checks*)

Options[VerifyKK4ScalarAmplitudes] = {TreeHead ->ATree};

VerifyKK4ScalarAmplitudes[n_,opts:OptionsPattern[]]:=(
	ClearKinematics[Range@n];
	GenFourScalarMomenta[{1,2},{3,4},Range[5,#],Seed->223234982,ParameterRange->100]&[n];
	Table[
		KKRelation[#,m,PassRules[VerifyKK4ScalarAmplitudes,KKRelation,opts]]&/@
		InitialKKHelicityAssignments4ScalarAmplitudes[n],
		{m,3, n-1}
	]//Activate//ToNum
);


VerifyKK2ScalarAmplitudes[n_]:=(
	ClearKinematics[Range@n];
	GenSpinors[Range@#,FourD->Range[3,#],Seed->223234982,ParameterRange->100]&[n];
	Table[KKRelation[#,m]&/@InitialKKHelicityAssignments2ScalarAmplitudes[n],{m,3,n-1}]//Activate//ToNum
);


InitialKKHelicityAssignments4ScalarAmplitudes[n_]:=
	DeleteDuplicatesBy[First]@Map[{Prepend[#/.Join[{2->"S",3->"s",4->"s"},Rule[#,"+"]&/@Range[5,n]],"S"],Prepend[#,1]}&]@
		Permutations[Range[2,n]];


InitialKKHelicityAssignments2ScalarAmplitudes[n_]:=
	DeleteDuplicatesBy[First]@Map[{Prepend[#/.Join[{2->"S"},Rule[#,"+"]&/@Range[3,n]],"S"],Prepend[#,1]}&]@
		Permutations[Range[2,n]];

Options[KKRelation] = {TreeHead ->ATree};

KKRelation[{hels_List,momenta_List},m_,opts:OptionsPattern[]]:=Module[{orderings,n=Length@hels},
	orderings=Join[{1},#,{m}]&/@ShuffleProduct[Range[n][[2;;m-1]],Reverse[Range[n][[m+1;;n]]]];
	Function[treeHead,Inactive[treeHead][Sequence@@hels][Sequence@@momenta]-
		(-1)^(n-m)Plus@@(Inactive[treeHead][Sequence@@hels[[#]]][Sequence@@momenta[[#]]]&/@orderings)][OptionValue[TreeHead]]
];

(* ::Chapter:: *)
(*Postamble*)


(* definitions for system functions *)


End[];


Protect[ "SpinorHelicity6D`Amplitudes`*" ];    (* protect exported symbols *)
Unprotect[$AmplitudesReference, $centerNGluonMax, $simplifyBCFWAmplitudes, $forceBerendsGiele];
(*LoadAmplitudesFile[];*)

EndPackage[ ];  (* end the package context *)
