(* ::Package:: *)

(* ::Title:: *)
(*Numerics*)


(* ::Chapter:: *)
(*Preamble*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`Numerics`", "SpinorHelicity6D`","SpinorHelicity6D`Unitarity`",
	"SpinorHelicity6D`Unitarity`LoopMomentum`","SpinorHelicity6D`Amplitudes`","Utils`","NumSymb`"];

Unprotect["SpinorHelicity6D`Numerics`*"];
ClearAll["SpinorHelicity6D`Numerics`*"];
ClearAll["SpinorHelicity6D`Numerics`Private`*"];


(* ::Section:: *)
(*Exports*)


Numerics::usage = "Numerics.wl is a package for generating numeric momentum configurations.";

$epsKinematic::usage = "";
$NumericsReference::usage = "";


(* ::Subsection::Closed:: *)
(*Momentum Utilities*)


CopyMomentum::usage = "";


NullMomentumQ::usage = "";


FlattenMomentum::usage = "";
FlattenMomentum::ref = "Reference momentum may not be null.";


SumMomenta::usage = "";
SumMomenta::ref = "Reference momentum may not be null.";


(* ::Subsection::Closed:: *)
(*Scalar Momenta*)


GenTwoScalarMomenta::usage = "";
GenFourScalarMomenta::usage = "";


(* ::Subsection::Closed:: *)
(*BCFW*)


GenScalarGluonShiftMomenta::usage = "";
GenGluonGluonShiftMomenta::usage = "";
zBCFW::usage = "Variable in BCFW shift numerics.";

GenAllPlusRisagerShiftMomenta::usage = "Generates shifted momenta according
	to eq. (32) of [hep-th/0501240]. Used for one-loop all-plus recursion relations.";

(* ::Subsection:: *)
(*Momenta in Special Limits*)


(* ::Subsubsection::Closed:: *)
(*Planar Kinematics*)


GenPlanarLimitKinematicsMassive::usage = ""


(* ::Subsubsection::Closed:: *)
(*GenCollinear*)

CollinearAngleToZ::usage = "";
CollinearZToAngle::usage = "";

Collinear::usage = "Collinear.wl is a package for generating collinear momentum configurations.";

GenCollinear::usage = "Generates collinear momentum configurations. Takes as argument Rule[momOrgi,{momCol1,momCol2}] and Rule[{pkIn,piIn},{pkOut,piOut}],
		       where the second rule specifies the momenta used for the recoil.";
		


(* ::Subsubsection::Closed:: *)
(*GenSoftMasslessMassive*)


GenSoftMasslessMassive::usage = "";


(* ::Subsubsection:: *)
(*Splitting Functions*)

$NumericsEpsilonDimReg::usage = "";
$NumericsEpsilonLimit::usage = "";
$RenormalizationScaleMu::usage = "";
MakeBoxes[$NumericsEpsilonDimReg,StandardForm|TraditionalForm] := InterpretationBox["\[Epsilon]",$NumericsEpsilonDimReg];
MakeBoxes[$NumericsEpsilonLimit,StandardForm|TraditionalForm] := InterpretationBox[SubscriptBox["\[Epsilon]","L"],$NumericsEpsilonLimit];
MakeBoxes[$RenormalizationScaleMu,StandardForm|TraditionalForm] := InterpretationBox["\[Mu]",$RenormalizationScaleMu];

SplittingFunctionTree::usage = "Generates tree-level splitting functions";
SplittingFunctionOneLoop::usage = "Generates one-loop splitting functions";
SoftTreeGGG::usage = "Generates tree-level soft functions";
SoftTreeSGG::usage = "Generates tree-level soft functions";
SoftTreeGGS::usage = "Generates tree-level soft functions";
SoftTreeSGS::usage = "Generates tree-level soft functions";

(* ::Chapter:: *)
(*Main*)


Begin["`Private`"];    (* begin the private context (implementation part) *)
protected = {};

raiseLowerMom4DNSpinorIndices[mat_/;Dimensions[mat]==={2,2}] := -Transpose[LeviCivitaTensor[2] . mat . LeviCivitaTensor[2]];


(* ::Section:: *)
(*Momentum Utilities*)


(* ::Subsection::Closed:: *)
(*Copy Momentum*)


(* Definition of the exported functions *)
CopyMomentum[init_,final_]:=
Module[{},
	ClearKinematics[final];
	Mom4DN[final][Null][$up] = Mom4DN[init][Null][$up];
	Mom4DN[final][$flat][$up] = Mom4DN[init][$flat][$up];
	Mom4DN[final][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[init][Null][$up]];
	Mom4DN[final][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[init][$flat][$up]];
	Mom4DVectorNFromMatrix[final][Null];
	Mom4DVectorNFromMatrix[final][$flat];
	
	SpinorDotN[final][$lam][$up] = SpinorDotN[init][$lam][$up];
	SpinorUndotN[final][$lam][$up] = SpinorUndotN[init][$lam][$up];
	SpinorDotN[final][$lam][$down] = SpinorDotN[init][$lam][$down];
	SpinorUndotN[final][$lam][$down] = SpinorUndotN[init][$lam][$down];
	SpinorDotN[final][$mu][$up] = SpinorDotN[init][$mu][$up];
	SpinorUndotN[final][$mu][$up] = SpinorUndotN[init][$mu][$up];
	SpinorDotN[final][$mu][$down] = SpinorDotN[init][$mu][$down];
	SpinorUndotN[final][$mu][$down] = SpinorUndotN[init][$mu][$down];
	
	SpinorDotN[final][$lam][i:1|2][Null]:=Indexed[SpinorDotN[final][$lam][$up],i];
	SpinorDotN[final][$lam][Null][i:1|2]:=Indexed[SpinorDotN[final][$lam][$down],i];
	SpinorUndotN[final][$lam][i:1|2][Null]:=Indexed[SpinorUndotN[final][$lam][$up],i];
	SpinorUndotN[final][$lam][Null][i:1|2]:=Indexed[SpinorUndotN[final][$lam][$down],i];
	
	SpinorDotN[final][$mu][i:1|2][Null]:=Indexed[SpinorDotN[final][$mu][$up],i];
	SpinorDotN[final][$mu][Null][i:1|2]:=Indexed[SpinorDotN[final][$mu][$down],i];
	SpinorUndotN[final][$mu][i:1|2][Null]:=Indexed[SpinorUndotN[final][$mu][$up],i];
	SpinorUndotN[final][$mu][Null][i:1|2]:=Indexed[SpinorUndotN[final][$mu][$down],i];
	
	
	Unprotect[Extramass,Extramasstilde];
	Extramass[final] = Extramass[init];
	Extramasstilde[final] = Extramasstilde[init];
	Protect[Extramass,Extramasstilde];
	ExtramassN[final] = ExtramassN[init];
	ExtramasstildeN[final] = ExtramasstildeN[init];
];


(* ::Subsection::Closed:: *)
(*NullMomentumQ*)


NullMomentumQ[mom_]:=And[ExtramassN[mom]===0,ExtramasstildeN[mom]===0,
	Mom4DVectorN[mom][Null][$up]===Mom4DVectorN[mom][$flat][$up],Mom4DVectorN[mom][Null][$down]===Mom4DVectorN[mom][$flat][$down]]


(* ::Subsection::Closed:: *)
(*Flattening Momenta*)


Options[FlattenMomentum] = {Normalization->Null,NumericalEvaluationFunction->Identity};

FlattenMomentum[momentum_,reference_,opts:OptionsPattern[]] /; NullMomentumQ[reference] :=Module[{},
	Mom4DN[momentum][$flat][$up] = Mom4DN[momentum][Null][$up]-
		Mom4DVectorN[momentum][Null][$down] . Mom4DVectorN[momentum][Null][$up]Mom4DN[reference][Null][$up]/
			(2Mom4DVectorN[momentum][Null][$up] . Mom4DVectorN[reference][Null][$down]);
	Mom4DN[momentum][$flat][$down] = raiseLowerMom4DNSpinorIndices@Mom4DN[momentum][$flat][$up];
	Mom4DVectorNFromMatrix[momentum][$flat];
	
	SpinorsFromMom4DN[momentum,Reference->reference,PassRules[FlattenMomentum,SpinorsFromMom4DN,opts]];
];

FlattenMomentum[momentum_,reference_] /; Message[FlattenMomentum::ref] := Null;



(* ::Subsection::Closed:: *)
(*New Momentum from Sum of Momenta*)


Options[SumMomenta]={Normalization->Null,NumericalEvaluationFunction->Identity};

SumMomenta[label_,momenta_List,reference_,opts:OptionsPattern[]] /; NullMomentumQ[reference]:=Module[{},
	Mom4DN[label][Null][$down] = Sum[Mom4DN[i][Null][$down],{i,momenta}];
	Mom4DN[label][Null][$up] = RaiseLowerMom4DNSpinorIndices@Mom4DN[label][Null][$down];
	Mom4DVectorNFromMatrix[label][Null];
	
	Mom4DN[label][$flat][$down] = Mom4DN[label][Null][$down] -
		(Mom4DVectorN[label][Null][$down] . Mom4DVectorN[label][Null][$up])/(2*Mom4DVectorN[label][Null][$down] . Mom4DVectorN[reference][Null][$up]) Mom4DN[reference][$flat][$down];
	Mom4DN[label][$flat][$up] = RaiseLowerMom4DNSpinorIndices@Mom4DN[label][$flat][$down];
	Mom4DVectorNFromMatrix[label][$flat];
	
	ExtramassN[label] = Sum[ExtramassN[i],{i,momenta}];
	ExtramasstildeN[label] = Sum[ExtramasstildeN[i],{i,momenta}];
	
	SpinorsFromMom4DN[label,Reference->reference, PassRules[SumMomenta,SpinorsFromMom4DN,opts]];
];
SumMomenta[momentum_,reference_,opts:OptionsPattern[]] /; Message[SumMomenta::ref] := Null;


(* ::Section::Closed:: *)
(*Scalar Momentum Generation*)


(* ::Subsection::Closed:: *)
(*One Scalar Pair Momentum Generation (with fixed mass)*)


Options[GenTwoScalarMomenta]=Join[{Mass2 -> 100}, Options[GenSpinors]];

GenTwoScalarMomenta[scalarPair1_List,masslessMomenta_List,opts:OptionsPattern[]]:=
Module[
	{l11=Unique["l"],l12=Unique["l"],l13=Unique["l"],tempMomenta},
	
	
	tempMomenta=Unique["p"]&/@Range@2;
	KillMasses[Join[tempMomenta,masslessMomenta]];
	GenSpinors[#,FourD->#,PassRules[GenTwoScalarMomenta,GenSpinors,opts]]&[Join[tempMomenta,masslessMomenta]];
	
	FixTriangleLoopMomentum6D[{{tempMomenta[[1]]},{tempMomenta[[2]]},masslessMomenta},
		LoopMomentumLabels->{l11,l12,l13},tVar->RandomInteger[{1,OptionValue[ParameterRange]}],
		MuVar->#,MuTildeVar->#,PostProcessLoopMomentum->False]&@Sqrt[OptionValue[Mass2]];
		
	MapThread[CopyMomentum,{{l11[$star],pm@l13[$star]},scalarPair1}];
	
	Unprotect[Extramass,Extramasstilde];
	Map[Unset@*Extramass,Hold@Evaluate@scalarPair1,{2}]//ReleaseHold;
	Map[Unset@*Extramasstilde,Hold@Evaluate@scalarPair1,{2}]//ReleaseHold;
	Protect[Extramass,Extramasstilde];

	ClearKinematics[
		l11[$star],l12[$star],l13[$star],
		l11[$unstar],l12[$unstar],l13[$unstar],
		Sequence@@tempMomenta
	];
];


(* ::Subsection::Closed:: *)
(*Two Scalar Pairs Momentum Generation (with fixed masses)*)


Options[GenFourScalarMomenta]=Options[GenSpinors];
GenFourScalarMomenta[scalarPair1_List,scalarPair2_List,masslessMomenta_List,opts:OptionsPattern[]]:=
Module[
	{l11=Unique["l"],l12=Unique["l"],l13=Unique["l"],
	 l21=Unique["l"],l22=Unique["l"],l23=Unique["l"],
	 tempMomenta},
	
	
	tempMomenta=Unique["p"]&/@Range@4;
	KillMasses[Join[tempMomenta,masslessMomenta]];
	GenSpinors[#,FourD->#,PassRules[GenFourScalarMomenta,GenSpinors,opts]]&[Join[tempMomenta,masslessMomenta]];
	
	FixTriangleLoopMomentum6D[{{tempMomenta[[1]]},{tempMomenta[[2]]},Join[tempMomenta[[{3,4}]],masslessMomenta]},
		LoopMomentumLabels->{l11,l12,l13},tVar->RandomInteger[{1,OptionValue[ParameterRange]}],
		MuVar->#,MuTildeVar->#,PostProcessLoopMomentum->False]&@RandomInteger[{1,OptionValue[ParameterRange]}];
		
	FixTriangleLoopMomentum6D[{{tempMomenta[[3]]},{tempMomenta[[4]]},Join[tempMomenta[[{1,2}]],masslessMomenta]},
		LoopMomentumLabels->{l21,l22,l23},tVar->RandomInteger[{1,OptionValue[ParameterRange]}],
		MuVar->#,MuTildeVar->-#,PostProcessLoopMomentum->False]&@RandomInteger[{1,OptionValue[ParameterRange]}];
	
	MapThread[CopyMomentum,{{l11[$star],pm@l13[$star],pm@l23[$star],l21[$star]},Join[scalarPair1,scalarPair2]}];

	Unprotect[Extramass,Extramasstilde];
	Map[Unset@*Extramass,Hold@Evaluate@Join[scalarPair1,scalarPair2],{2}]//ReleaseHold;
	Map[Unset@*Extramasstilde,Hold@Evaluate@Join[scalarPair1,scalarPair2],{2}]//ReleaseHold;
	Protect[Extramass,Extramasstilde];

	ClearKinematics[
		l11[$star],l12[$star],l13[$star],
		l11[$unstar],l12[$unstar],l13[$unstar],
		l21[$star],l22[$star],l23[$star],
		l21[$unstar],l22[$unstar],l23[$unstar],
		Sequence@@tempMomenta
	];
];


(* ::Section::Closed:: *)
(*BCFW*)


(* ::Subsection::Closed:: *)
(*Scalar Gluon BCFW Shift Numerics*)


(* Here we perform a scalar-gluon shift, where the shift momentum is defined by flattening the scalar momentum against the gluon one. As a convention we assume
	the momentum K to flow from the left to the right amplitude, where the left one carries the scalar, and the right one the gluon *)

Options[GenScalarGluonShiftMomenta] = {FlattenReference -> $NumericsReference};
GenScalarGluonShiftMomenta[momenta_List,shiftSMom_,shiftGMom_,KLabel_,channel_,opts:OptionsPattern[]]:=
    Module[{zVar,shiftMom,shiftMomVector},
	
	ClearKinematics[{KLabel,OverHat[KLabel],OverHat[shiftGMom],OverHat[shiftSMom],shiftMom}];

	(* Projecting with shiftGMom is essential here*)
	FlattenMomentum[shiftSMom,shiftGMom];
	
	zVar = S6MomN@@channel/ChainN[$angle,shiftSMom,{UnderBar/@PlusM@@channel},shiftGMom,$square];
	
	shiftMomVector = ChainN[$angle,shiftSMom,{PauliMatrixSpinor[#, $up]},shiftGMom,$square]&/@Range[0,3];
	Mom4DNFromVector[shiftMom][$up][shiftMomVector];
	Mom4DVectorNFromMatrix[shiftMom][Null];
	
	Mom4DN[OverHat[shiftSMom]][Null][$up] = Mom4DN[shiftSMom][Null][$up] + zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftSMom]][$flat][$up] = Mom4DN[shiftSMom][$flat][$up] + zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftSMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftSMom]][Null][$up]];
	Mom4DN[OverHat[shiftSMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftSMom]][$flat][$up]];
	Mom4DVectorNFromMatrix[OverHat[shiftSMom]][Null];
	Mom4DVectorNFromMatrix[OverHat[shiftSMom]][$flat];
	
	SpinorDotN[OverHat[shiftSMom]][$lam][$up] = SpinorDotN[shiftSMom][$lam][$up]+ zBCFW SpinorDotN[shiftGMom][$lam][$up];
	SpinorUndotN[OverHat[shiftSMom]][$lam][$up] = SpinorUndotN[shiftSMom][$lam][$up];
	SpinorDotN[OverHat[shiftSMom]][$lam][$down] = SpinorDotN[shiftSMom][$lam][$down]+ zBCFW SpinorDotN[shiftGMom][$lam][$down];
	SpinorUndotN[OverHat[shiftSMom]][$lam][$down] = SpinorUndotN[shiftSMom][$lam][$down];
	SpinorDotN[OverHat[shiftSMom]][$mu][$up] = SpinorDotN[shiftSMom][$mu][$up];
	SpinorUndotN[OverHat[shiftSMom]][$mu][$up] = SpinorUndotN[shiftSMom][$mu][$up];
	SpinorDotN[OverHat[shiftSMom]][$mu][$down] = SpinorDotN[shiftSMom][$mu][$down];
	SpinorUndotN[OverHat[shiftSMom]][$mu][$down] = SpinorUndotN[shiftSMom][$mu][$down];
	
	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[shiftSMom]] = Extramass[shiftSMom];
	Extramasstilde[OverHat[shiftSMom]] = Extramasstilde[shiftSMom];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[shiftSMom]] = ExtramassN[shiftSMom];
	ExtramasstildeN[OverHat[shiftSMom]] = ExtramasstildeN[shiftSMom];
	
	
	
	Mom4DN[OverHat[shiftGMom]][Null][$up] = Mom4DN[shiftGMom][Null][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGMom]][$flat][$up] = Mom4DN[shiftGMom][$flat][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGMom]][Null][$up]];
	Mom4DN[OverHat[shiftGMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGMom]][$flat][$up]];
	Mom4DVectorNFromMatrix[OverHat[shiftGMom]][Null];
	Mom4DVectorNFromMatrix[OverHat[shiftGMom]][$flat];
	
	SpinorDotN[OverHat[shiftGMom]][$lam][$up] = SpinorDotN[shiftGMom][$lam][$up];
	SpinorUndotN[OverHat[shiftGMom]][$lam][$up] = SpinorUndotN[shiftGMom][$lam][$up] - zBCFW SpinorUndotN[shiftSMom][$lam][$up];
	SpinorDotN[OverHat[shiftGMom]][$lam][$down] = SpinorDotN[shiftGMom][$lam][$down];
	SpinorUndotN[OverHat[shiftGMom]][$lam][$down] = SpinorUndotN[shiftGMom][$lam][$down] - zBCFW SpinorUndotN[shiftSMom][$lam][$down];
	SpinorDotN[OverHat[shiftGMom]][$mu][$up] = SpinorDotN[shiftGMom][$mu][$up];
	SpinorUndotN[OverHat[shiftGMom]][$mu][$up] = SpinorUndotN[shiftGMom][$mu][$up];
	SpinorDotN[OverHat[shiftGMom]][$mu][$down] = SpinorDotN[shiftGMom][$mu][$down];
	SpinorUndotN[OverHat[shiftGMom]][$mu][$down] = SpinorUndotN[shiftGMom][$mu][$down];
	
	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[shiftGMom]] = Extramass[shiftGMom];
	Extramasstilde[OverHat[shiftGMom]] = Extramasstilde[shiftGMom];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[shiftGMom]] = ExtramassN[shiftGMom];
	ExtramasstildeN[OverHat[shiftGMom]] = ExtramasstildeN[shiftGMom];
	
	
	
	SumMomenta[KLabel,channel,shiftGMom];
	
	Mom4DN[OverHat[KLabel]][Null][$up] = Mom4DN[KLabel][Null][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[KLabel]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[KLabel]][Null][$up]];
	Mom4DVectorNFromMatrix[OverHat[KLabel]][Null];
	
	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[KLabel]] = Extramass[KLabel];
	Extramasstilde[OverHat[KLabel]] = Extramasstilde[KLabel];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[KLabel]] = ExtramassN[KLabel];
	ExtramasstildeN[OverHat[KLabel]] = ExtramasstildeN[KLabel];
	
	(* Normalization should be irrelevant, as the phase weight of OverHat[K] is zero, and an OverHat[K] spinor cannot be
	   changed to a K one*)
	(* Okay, it turns out that the choice of reference momentum can make a difference in
	 very special cases. The one I stumbled over is a case where I want to cancel a 0/0
	 divergence first. Take the channel S_{34}, where 3 is shiftGMom. In that case
	  SpinorAngleBracket[OverHat[K],OverHat[3]]SpinorAngleBracket[OverHat[K],
	  4]/SpinorAngleBracket[4,OverHat[3]]^2\
	  Should be finite. However, when setting the reference momentum for OverHat[K] to
	   be shiftGMom, SpinorAngleBracket[4,OverHat[K]] always evaluates to 0, and
	   we cannot take the nice zBCFW->z34 limit. *)
	FlattenMomentum[OverHat[KLabel], OptionValue[FlattenReference]];
	
	ClearKinematics[shiftMom];
	
	{zBCFW->zVar}
];


(* ::Subsection::Closed:: *)
(*Gluon Gluon BCFW Shift Numerics*)


(* Here we perform a gluon-gluon shift, where we shift \tilde{\lambda} of *)
(* shiftGLatMom and \lambda of shiftGLaMom. *)
(* As a convention we assume the momentum K to flow from the amplitude of *)
(* shiftGLatMom to that of shiftGLaMom *)


GenGluonGluonShiftMomenta[momenta_List,shiftGLatMom_,shiftGLaMom_,KLabel_,channel_]:=
    Module[{zVar,shiftMom,shiftMomVector},

	ClearKinematics[{KLabel,OverHat[KLabel],OverHat[shiftGLatMom],OverHat[shiftGLaMom],shiftMom}];

	zVar = S6MomN@@channel/ChainN[$angle,shiftGLatMom,{UnderBar/@PlusM@@channel},shiftGLaMom,$square];

	shiftMomVector = ChainN[$angle,shiftGLatMom,{PauliMatrixSpinor[#, $up]},shiftGLaMom,$square]&/@Range[0,3];
	Mom4DNFromVector[shiftMom][$up][shiftMomVector];
	Mom4DVectorNFromMatrix[shiftMom][Null];

	Mom4DN[OverHat[shiftGLatMom]][Null][$up] = Mom4DN[shiftGLatMom][Null][$up] + zBCFW/2Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGLatMom]][$flat][$up] = Mom4DN[shiftGLatMom][$flat][$up] + zBCFW/2Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGLatMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGLatMom]][Null][$up]];
	Mom4DN[OverHat[shiftGLatMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGLatMom]][$flat][$up]];
	Mom4DVectorNFromMatrix[OverHat[shiftGLatMom]][Null];
	Mom4DVectorNFromMatrix[OverHat[shiftGLatMom]][$flat];

	SpinorUndotN[OverHat[shiftGLatMom]][$lam][$up] = SpinorUndotN[shiftGLatMom][$lam][$up];
	SpinorDotN[OverHat[shiftGLatMom]][$lam][$up] = SpinorDotN[shiftGLatMom][$lam][$up] + zBCFW SpinorDotN[shiftGLaMom][$lam][$up];
	SpinorUndotN[OverHat[shiftGLatMom]][$lam][$down] = SpinorUndotN[shiftGLatMom][$lam][$down];
	SpinorDotN[OverHat[shiftGLatMom]][$lam][$down] = SpinorDotN[shiftGLatMom][$lam][$down] + zBCFW SpinorDotN[shiftGLaMom][$lam][$down];
	SpinorUndotN[OverHat[shiftGLatMom]][$mu][$up] = SpinorUndotN[shiftGLatMom][$mu][$up];
	SpinorDotN[OverHat[shiftGLatMom]][$mu][$up] = SpinorDotN[shiftGLatMom][$mu][$up];
	SpinorUndotN[OverHat[shiftGLatMom]][$mu][$down] = SpinorUndotN[shiftGLatMom][$mu][$down];
	SpinorDotN[OverHat[shiftGLatMom]][$mu][$down] = SpinorDotN[shiftGLatMom][$mu][$down];

	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[shiftGLatMom]] = Extramass[shiftGLatMom];
	Extramasstilde[OverHat[shiftGLatMom]] = Extramasstilde[shiftGLatMom];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[shiftGLatMom]] = ExtramassN[shiftGLatMom];
	ExtramasstildeN[OverHat[shiftGLatMom]] = ExtramasstildeN[shiftGLatMom];



	Mom4DN[OverHat[shiftGLaMom]][Null][$up] = Mom4DN[shiftGLaMom][Null][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGLaMom]][$flat][$up] = Mom4DN[shiftGLaMom][$flat][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[shiftGLaMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGLaMom]][Null][$up]];
	Mom4DN[OverHat[shiftGLaMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[shiftGLaMom]][$flat][$up]];
	Mom4DVectorNFromMatrix[OverHat[shiftGLaMom]][Null];
	Mom4DVectorNFromMatrix[OverHat[shiftGLaMom]][$flat];

	SpinorDotN[OverHat[shiftGLaMom]][$lam][$up] = SpinorDotN[shiftGLaMom][$lam][$up];
	SpinorUndotN[OverHat[shiftGLaMom]][$lam][$up] = SpinorUndotN[shiftGLaMom][$lam][$up] - zBCFW SpinorUndotN[shiftGLatMom][$lam][$up];
	SpinorDotN[OverHat[shiftGLaMom]][$lam][$down] = SpinorDotN[shiftGLaMom][$lam][$down];
	SpinorUndotN[OverHat[shiftGLaMom]][$lam][$down] = SpinorUndotN[shiftGLaMom][$lam][$down] - zBCFW SpinorUndotN[shiftGLatMom][$lam][$down];
	SpinorDotN[OverHat[shiftGLaMom]][$mu][$up] = SpinorDotN[shiftGLaMom][$mu][$up];
	SpinorUndotN[OverHat[shiftGLaMom]][$mu][$up] = SpinorUndotN[shiftGLaMom][$mu][$up];
	SpinorDotN[OverHat[shiftGLaMom]][$mu][$down] = SpinorDotN[shiftGLaMom][$mu][$down];
	SpinorUndotN[OverHat[shiftGLaMom]][$mu][$down] = SpinorUndotN[shiftGLaMom][$mu][$down];

	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[shiftGLaMom]] = Extramass[shiftGLaMom];
	Extramasstilde[OverHat[shiftGLaMom]] = Extramasstilde[shiftGLaMom];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[shiftGLaMom]] = ExtramassN[shiftGLaMom];
	ExtramasstildeN[OverHat[shiftGLaMom]] = ExtramasstildeN[shiftGLaMom];



	SumMomenta[KLabel,channel,shiftGLaMom];

	Mom4DN[OverHat[KLabel]][Null][$up] = Mom4DN[KLabel][Null][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];
	Mom4DN[OverHat[KLabel]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[KLabel]][Null][$up]];
	Mom4DVectorNFromMatrix[OverHat[KLabel]][Null];

	Unprotect[Extramass,Extramasstilde];
	Extramass[OverHat[KLabel]] = Extramass[KLabel];
	Extramasstilde[OverHat[KLabel]] = Extramasstilde[KLabel];
	Protect[Extramass,Extramasstilde];
	ExtramassN[OverHat[KLabel]] = ExtramassN[KLabel];
	ExtramasstildeN[OverHat[KLabel]] = ExtramasstildeN[KLabel];

	(* Normalization should be irrelevant, as the phase weight of OverHat[K] is zero, and an OverHat[K] spinor cannot be
	   changed to a K one*)
	FlattenMomentum[OverHat[KLabel],shiftGLaMom];

	ClearKinematics[shiftMom];

	{zBCFW->zVar}
];

(* ::Subsection::Closed:: *)
(* Risager-like shift, as given in eq.(32) of [hep-th/0501240] *)



GenAllPlusRisagerShiftMomenta[momenta_List, {jMom_,lMom_,nMom_},KLabel_,channel_]:=
		Module[{zVar,shiftjMom,shiftlMom,shiftnMom,shiftjMomVector,shiftlMomVector,shiftnMomVector},

			ClearKinematics[{KLabel,OverHat[KLabel],OverHat[jMom],OverHat[lMom],OverHat[nMom],
				shiftjMom,shiftlMom,shiftnMom}];

(*			zVar = S6MomN@@channel/ChainN[$angle,shiftGLatMom,{UnderBar/@PlusM@@channel},shiftGLaMom,$square];*)


			(* OverHat[j] *)
			shiftjMomVector =
					(
						ChainN[$angle,jMom,{PauliMatrixSpinor[#, $up]},lMom, $square]+
							SpinorAngleBracketN[nMom,jMom]/SpinorAngleBracketN[lMom,jMom]*
           		ChainN[$angle,jMom,{PauliMatrixSpinor[#, $up]},nMom,$square]
					)&/@Range[0,3];
			Mom4DNFromVector[shiftjMom][$up][shiftjMomVector];
			Mom4DVectorNFromMatrix[shiftjMom][Null];

			Mom4DN[OverHat[jMom]][Null][$up] = Mom4DN[jMom][Null][$up] - zBCFW/2 Mom4DN[shiftjMom][Null][$up];
			Mom4DN[OverHat[jMom]][$flat][$up] = Mom4DN[jMom][$flat][$up] - zBCFW/2 Mom4DN[shiftjMom][Null][$up];
			Mom4DN[OverHat[jMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[jMom]][Null][$up]];
			Mom4DN[OverHat[jMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[jMom]][$flat][$up]];
			Mom4DVectorNFromMatrix[OverHat[jMom]][Null];
			Mom4DVectorNFromMatrix[OverHat[jMom]][$flat];

			SpinorUndotN[OverHat[jMom]][$lam][$up] = SpinorUndotN[jMom][$lam][$up];
			SpinorDotN[OverHat[jMom]][$lam][$up] = SpinorDotN[jMom][$lam][$up] - zBCFW SpinorDotN[lMom][$lam][$up] -
				zBCFW SpinorAngleBracketN[nMom, jMom]/SpinorAngleBracketN[lMom,jMom] SpinorDotN[nMom][$lam][$up];
			SpinorUndotN[OverHat[jMom]][$lam][$down] = SpinorUndotN[jMom][$lam][$down];
			SpinorDotN[OverHat[jMom]][$lam][$down] = SpinorDotN[jMom][$lam][$down] - zBCFW SpinorDotN[lMom][$lam][$down] -
					zBCFW SpinorAngleBracketN[nMom, jMom]/SpinorAngleBracketN[lMom,jMom] SpinorDotN[nMom][$lam][$down];
			SpinorUndotN[OverHat[jMom]][$mu][$up] = SpinorUndotN[jMom][$mu][$up];
			SpinorDotN[OverHat[jMom]][$mu][$up] = SpinorDotN[jMom][$mu][$up];
			SpinorUndotN[OverHat[jMom]][$mu][$down] = SpinorUndotN[jMom][$mu][$down];
			SpinorDotN[OverHat[jMom]][$mu][$down] = SpinorDotN[jMom][$mu][$down];

			Unprotect[Extramass,Extramasstilde];
			Extramass[OverHat[jMom]] = Extramass[jMom];
			Extramasstilde[OverHat[jMom]] = Extramasstilde[jMom];
			Protect[Extramass,Extramasstilde];
			ExtramassN[OverHat[jMom]] = ExtramassN[jMom];
			ExtramasstildeN[OverHat[jMom]] = ExtramasstildeN[jMom];

			(* OverHat[l] *)
			shiftlMomVector = ChainN[$angle,jMom,{PauliMatrixSpinor[#, $up]},lMom, $square]&/@Range[0,3];
			Mom4DNFromVector[shiftlMom][$up][shiftlMomVector];
			Mom4DVectorNFromMatrix[shiftlMom][Null];

			Mom4DN[OverHat[lMom]][Null][$up] = Mom4DN[lMom][Null][$up] + zBCFW/2 Mom4DN[shiftlMom][Null][$up];
			Mom4DN[OverHat[lMom]][$flat][$up] = Mom4DN[lMom][$flat][$up] + zBCFW/2 Mom4DN[shiftlMom][Null][$up];
			Mom4DN[OverHat[lMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[lMom]][Null][$up]];
			Mom4DN[OverHat[lMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[lMom]][$flat][$up]];
			Mom4DVectorNFromMatrix[OverHat[lMom]][Null];
			Mom4DVectorNFromMatrix[OverHat[lMom]][$flat];

			SpinorDotN[OverHat[lMom]][$lam][$up] = SpinorDotN[lMom][$lam][$up];
			SpinorUndotN[OverHat[lMom]][$lam][$up] = SpinorUndotN[lMom][$lam][$up] + zBCFW SpinorUndotN[jMom][$lam][$up];
			SpinorDotN[OverHat[lMom]][$lam][$down] = SpinorDotN[lMom][$lam][$down];
			SpinorUndotN[OverHat[lMom]][$lam][$down] = SpinorUndotN[lMom][$lam][$down] + zBCFW SpinorUndotN[jMom][$lam][$down];
			SpinorDotN[OverHat[lMom]][$mu][$up] = SpinorDotN[lMom][$mu][$up];
			SpinorUndotN[OverHat[lMom]][$mu][$up] = SpinorUndotN[lMom][$mu][$up];
			SpinorDotN[OverHat[lMom]][$mu][$down] = SpinorDotN[lMom][$mu][$down];
			SpinorUndotN[OverHat[lMom]][$mu][$down] = SpinorUndotN[lMom][$mu][$down];

			Unprotect[Extramass,Extramasstilde];
			Extramass[OverHat[lMom]] = Extramass[lMom];
			Extramasstilde[OverHat[lMom]] = Extramasstilde[lMom];
			Protect[Extramass,Extramasstilde];
			ExtramassN[OverHat[lMom]] = ExtramassN[lMom];
			ExtramasstildeN[OverHat[lMom]] = ExtramasstildeN[lMom];


			(* OverHat[n] *)
			shiftnMomVector = SpinorAngleBracketN[nMom,jMom]/SpinorAngleBracketN[lMom,jMom]*
				ChainN[$angle,jMom,{PauliMatrixSpinor[#, $up]},nMom, $square]&/@Range[0,3];
			Mom4DNFromVector[shiftnMom][$up][shiftnMomVector];
			Mom4DVectorNFromMatrix[shiftnMom][Null];

			Mom4DN[OverHat[nMom]][Null][$up] = Mom4DN[nMom][Null][$up] + zBCFW/2 Mom4DN[shiftnMom][Null][$up];
			Mom4DN[OverHat[nMom]][$flat][$up] = Mom4DN[nMom][$flat][$up] + zBCFW/2 Mom4DN[shiftnMom][Null][$up];
			Mom4DN[OverHat[nMom]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[nMom]][Null][$up]];
			Mom4DN[OverHat[nMom]][$flat][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[nMom]][$flat][$up]];
			Mom4DVectorNFromMatrix[OverHat[nMom]][Null];
			Mom4DVectorNFromMatrix[OverHat[nMom]][$flat];

			SpinorDotN[OverHat[nMom]][$lam][$up] = SpinorDotN[nMom][$lam][$up];
			SpinorUndotN[OverHat[nMom]][$lam][$up] = SpinorUndotN[nMom][$lam][$up] +
       zBCFW SpinorAngleBracketN[nMom,jMom]/SpinorAngleBracketN[lMom,jMom] SpinorUndotN[jMom][$lam][$up];
			SpinorDotN[OverHat[nMom]][$lam][$down] = SpinorDotN[nMom][$lam][$down];
			SpinorUndotN[OverHat[nMom]][$lam][$down] = SpinorUndotN[nMom][$lam][$down] +
        zBCFW SpinorAngleBracketN[nMom,jMom]/SpinorAngleBracketN[lMom,jMom] SpinorUndotN[jMom][$lam][$down];
			SpinorDotN[OverHat[nMom]][$mu][$up] = SpinorDotN[nMom][$mu][$up];
			SpinorUndotN[OverHat[nMom]][$mu][$up] = SpinorUndotN[nMom][$mu][$up];
			SpinorDotN[OverHat[nMom]][$mu][$down] = SpinorDotN[nMom][$mu][$down];
			SpinorUndotN[OverHat[nMom]][$mu][$down] = SpinorUndotN[nMom][$mu][$down];

			Unprotect[Extramass,Extramasstilde];
			Extramass[OverHat[nMom]] = Extramass[nMom];
			Extramasstilde[OverHat[nMom]] = Extramasstilde[nMom];
			Protect[Extramass,Extramasstilde];
			ExtramassN[OverHat[nMom]] = ExtramassN[nMom];
			ExtramasstildeN[OverHat[nMom]] = ExtramasstildeN[nMom];


(*			SumMomenta[KLabel,channel,shiftGLaMom];*)

(*			Mom4DN[OverHat[KLabel]][Null][$up] = Mom4DN[KLabel][Null][$up] - zBCFW/2 Mom4DN[shiftMom][Null][$up];*)
(*			Mom4DN[OverHat[KLabel]][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[OverHat[KLabel]][Null][$up]];*)
(*			Mom4DVectorNFromMatrix[OverHat[KLabel]][Null];*)

(*			Unprotect[Extramass,Extramasstilde];*)
(*			Extramass[OverHat[KLabel]] = Extramass[KLabel];*)
(*			Extramasstilde[OverHat[KLabel]] = Extramasstilde[KLabel];*)
(*			Protect[Extramass,Extramasstilde];*)
(*			ExtramassN[OverHat[KLabel]] = ExtramassN[KLabel];*)
(*			ExtramasstildeN[OverHat[KLabel]] = ExtramasstildeN[KLabel];*)

(*			*)(* Normalization should be irrelevant, as the phase weight of OverHat[K] is zero, and an OverHat[K] spinor cannot be*)
(*         changed to a K one*)
(*			FlattenMomentum[OverHat[KLabel],shiftGLaMom];*)

			ClearKinematics[shiftjMom,shiftlMom,shiftnMom];

			{zBCFW->zVar}
		];



(* ::Section:: *)
(*Momenta in Special Limits*)


(* ::Subsection::Closed:: *)
(*Planar Limit*)


(* Parameterized kinematics for two massive momenta with the same mass (i.e. a scalar pair),
	where the limit of vanishing $epsKinematic parametrizes a planar limit. For example,
	GenPlanarLimitKinematicsMassive[1\[Rule]Tilde[1],{3,4},2\[Rule]Tilde[2]] generates momenta Tilde[1], Tilde[2],
	such that Tilde[1]+Tilde[2] = p1+p2, Extramass[Tilde[1]]\[Equal]-Extramass[Tilde[2]],
	Extramasstilde[Tilde[1]]\[Equal]-Extramasstilde[Tilde[2]] and for example <3|Tilde[1]|4] vanishes  for $epsKinematic\[Rule]0 *)


Options[GenPlanarLimitKinematicsMassive]={LimitVariable->$epsKinematic,ParameterRange->100,Reference->Global`q,Seed->Null,FreeParameter->Null}


GenPlanarLimitKinematicsMassive[momLimit:HoldPattern[pL_->hatpL_],plane:{pp1_,pp2_},massiveAbsorber:HoldPattern[pAbs_->hatpAbs_],opts:OptionsPattern[]]:=
Module[{a,b,c,muSqHatpL,k,killEpsilon,gamma},

	(*NumericalEvaluationFunction passed to SpinorsFromMom4DN in FlattenMomentum*)
	killEpsilon = ReplaceAll[OptionValue[LimitVariable]->0];
	
	If[IntegerQ@OptionValue[Seed], 
		SeedRandom[OptionValue[Seed]]; 
	];
	
	k = OptionValue[Reference];
	
	(* Choosing dependence on plane momenta randomly *)
	If[OptionValue[FreeParameter] =!= Null,
		b = OptionValue[FreeParameter];,
		b = RandomInteger[OptionValue[ParameterRange]];
	];
	
	a = -(S6Mom[pL,pAbs]+b(S6Mom[pL,pp2]+S6Mom[pAbs,pp2])+OptionValue[LimitVariable](S6Mom[pL,k]+S6Mom[pAbs,k]))/(S6Mom[pL,pp1]+S6Mom[pAbs,pp1])//ToNum;
	
	Mom4DN[hatpL][Null][$up] = a Mom4DN[pp1][Null][$up] + b Mom4DN[pp2][Null][$up] + OptionValue[LimitVariable] Mom4DN[k][Null][$up];
	Mom4DN[hatpL][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[hatpL][Null][$up]];
	Mom4DVectorNFromMatrix[hatpL][Null];
	FlattenMomentum[hatpL,k,NumericalEvaluationFunction->killEpsilon];
	
	muSqHatpL = -(a*b*S6MomN[pp1,pp2]+OptionValue[LimitVariable](a S6MomN[pp1,k]+b S6MomN[pp2,k]));	
	ExtramassN[hatpL] = muSqHatpL;
	ExtramasstildeN[hatpL] = -1;
	
	(* Massive Absorber *)
	(*c =( 2Mom4DVectorN[hatpL][Null][$up] . (Mom4DVectorN[pL][Null][$down]+Mom4DVectorN[pAbs][Null][$down])-S6MomN[pL,pAbs])/
		(2Mom4DVectorN[k][Null][$up] . (Mom4DVectorN[pL][Null][$down]+Mom4DVectorN[pAbs][Null][$down]-Mom4DVectorN[hatpL][Null][$down]));*)

	Mom4DN[hatpAbs][Null][$up] = Mom4DN[pL][Null][$up]+Mom4DN[pAbs][Null][$up] - Mom4DN[hatpL][Null][$up](* + c Mom4DN[k][Null][$up]*);
	Mom4DN[hatpAbs][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[hatpAbs][Null][$up]];
	Mom4DVectorNFromMatrix[hatpAbs][Null];
	FlattenMomentum[hatpAbs,k,NumericalEvaluationFunction->killEpsilon];
		
	ExtramassN[hatpAbs] = -ExtramassN[hatpL];
	ExtramasstildeN[hatpAbs] = -ExtramasstildeN[hatpL];

	(*
	(* Absorbing recoil momentum, such that we again have momentum conservation *)
	
	gamma = 1 + S6MomN[pi,pL,pm@hatpL]/(S6MomN[pi,pj]+S6MomN[pj,pL,pm@hatpL]);
	
	Mom4DN[hatpj][Null][$up] = gamma Mom4DN[pj][Null][$up] ;
	Mom4DN[hatpj][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[hatpj][Null][$up]];
	Mom4DN[hatpj][$flat][$up] = Mom4DN[hatpj][Null][$up];
	Mom4DN[hatpj][$flat][$down] = Mom4DN[hatpj][Null][$down];
	Mom4DVectorNFromMatrix[hatpj][Null];
	Mom4DVectorNFromMatrix[hatpj][$flat];
	SpinorsFromMom4DN[hatpj,Normalization->pj,NumericalEvaluationFunction->killEpsilon];
	
	Mom4DN[hatpi][Null][$up] = Mom4DN[pi][Null][$up] + Mom4DN[pL][Null][$up] - Mom4DN[hatpL][Null][$up] + (1-gamma) Mom4DN[pj][Null][$up];
	Mom4DN[hatpi][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[hatpi][Null][$up]];
	Mom4DN[hatpi][$flat][$up] = Mom4DN[hatpi][Null][$up];
	Mom4DN[hatpi][$flat][$down] = Mom4DN[hatpi][Null][$down];
	Mom4DVectorNFromMatrix[hatpi][Null];
	Mom4DVectorNFromMatrix[hatpi][$flat];
	SpinorsFromMom4DN[hatpi,Normalization\[Rule]pi,NumericalEvaluationFunction->killEpsilon];
	
	KillMasses[{hatpi,hatpj}];*)
		
]


(* ::Subsection::Closed:: *)
(*GenCollinear*)

CollinearAngleToZ[theta_]:=Cos[theta]^2;
CollinearZToAngle[z_]:=ArcCos[Sqrt[z]];

Options[GenCollinear] = {Reference -> q, AngleVar -> theta, LimitVar-> $NumericsEpsilonLimit, MomentumConserving->True};

GenCollinear[Rule[momOrig_,{momCol1_,momCol2_}], opts : OptionsPattern[]] := GenCollinear[Rule[momOrig,{momCol1,momCol2}],Rule[{Null,Null},{Null,Null}],MomentumConserving->False,opts];

GenCollinear[Rule[momOrig_,{momCol1_,momCol2_}],Rule[{pk_,pi_},{pkOut_,piOut_}], opts : OptionsPattern[]]:=Module[{protected,gamma},

	Mom4DN[momCol1][Null][$up] = Mom4DN[momOrig][Null][$up]*Cos[OptionValue[AngleVar]]^2+OptionValue[LimitVar]^2*Sin[OptionValue[AngleVar]]^2*Mom4DN[OptionValue[Reference]][Null][$up]-
					OptionValue[LimitVar]Sin[OptionValue[AngleVar]]*Cos[OptionValue[AngleVar]]*
					(KroneckerProduct[SpinorUndotN[momOrig][$lam][$up],SpinorDotN[OptionValue[Reference]][$lam][$up]]+KroneckerProduct[SpinorUndotN[OptionValue[Reference]][$lam][$up],SpinorDotN[momOrig][$lam][$up]]);
	SpinorHelicity6D`Mom4DN[momCol1][$flat][$up] = Mom4DN[momCol1][Null][$up];
	SpinorHelicity6D`Mom4DN[momCol1][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momCol1][Null][$up]];
	SpinorHelicity6D`Mom4DN[momCol1][$flat][$down] = Mom4DN[momCol1][Null][$down];
	Mom4DVectorNFromMatrix[momCol1][Null];
	Mom4DVectorNFromMatrix[momCol1][$flat];
	
	SpinorHelicity6D`SpinorDotN[momCol1][$lam][$up] = SpinorDotN[momOrig][$lam][$up]*Cos[OptionValue[AngleVar]]-SpinorDotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol1][$lam][$up] = SpinorUndotN[momOrig][$lam][$up]*Cos[OptionValue[AngleVar]]-SpinorUndotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol1][$lam][$down] = SpinorDotN[momOrig][$lam][$down]*Cos[OptionValue[AngleVar]]-SpinorDotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol1][$lam][$down] = SpinorUndotN[momOrig][$lam][$down]*Cos[OptionValue[AngleVar]]-SpinorUndotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol1][$mu][$up] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol1][$mu][$up] = 0;
	SpinorHelicity6D`SpinorDotN[momCol1][$mu][$down] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol1][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	SpinorHelicity6D`Extramass[momCol1] = 0;
	SpinorHelicity6D`Extramasstilde[momCol1] = 0;
	Protect[Extramass,Extramasstilde];
	SpinorHelicity6D`ExtramassN[momCol1] = 0;
	SpinorHelicity6D`ExtramasstildeN[momCol1] = 0;



	Mom4DN[momCol2][Null][$up] = Mom4DN[momOrig][Null][$up]*Sin[OptionValue[AngleVar]]^2+OptionValue[LimitVar]^2*Cos[OptionValue[AngleVar]]^2*Mom4DN[OptionValue[Reference]][Null][$up]+
					OptionValue[LimitVar]Sin[OptionValue[AngleVar]]*Cos[OptionValue[AngleVar]]*
					(KroneckerProduct[SpinorUndotN[momOrig][$lam][$up],SpinorDotN[OptionValue[Reference]][$lam][$up]]+KroneckerProduct[SpinorUndotN[OptionValue[Reference]][$lam][$up],SpinorDotN[momOrig][$lam][$up]]);
	SpinorHelicity6D`Mom4DN[momCol2][$flat][$up] = Mom4DN[momCol2][Null][$up];
	SpinorHelicity6D`Mom4DN[momCol2][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momCol2][Null][$up]];
	SpinorHelicity6D`Mom4DN[momCol2][$flat][$down] = Mom4DN[momCol2][Null][$down];
	Mom4DVectorNFromMatrix[momCol2][Null];
	Mom4DVectorNFromMatrix[momCol2][$flat];
	
	SpinorHelicity6D`SpinorDotN[momCol2][$lam][$up] = SpinorDotN[momOrig][$lam][$up]*Sin[OptionValue[AngleVar]]+SpinorDotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol2][$lam][$up] = SpinorUndotN[momOrig][$lam][$up]*Sin[OptionValue[AngleVar]]+SpinorUndotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol2][$lam][$down] = SpinorDotN[momOrig][$lam][$down]*Sin[OptionValue[AngleVar]]+SpinorDotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol2][$lam][$down] = SpinorUndotN[momOrig][$lam][$down]*Sin[OptionValue[AngleVar]]+SpinorUndotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol2][$mu][$up] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol2][$mu][$up] = 0;
	SpinorHelicity6D`SpinorDotN[momCol2][$mu][$down] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol2][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	SpinorHelicity6D`Extramass[momCol2] = 0;
	SpinorHelicity6D`Extramasstilde[momCol2] = 0;
	Protect[Extramass,Extramasstilde];
	SpinorHelicity6D`ExtramassN[momCol2] = 0;
	SpinorHelicity6D`ExtramasstildeN[momCol2] = 0;
	

	If[OptionValue[MomentumConserving],

	   gamma = 1 + (-OptionValue[LimitVar]^2*ChainN[$angle,pi,{OptionValue[Reference]},pi,$square])/(ChainN[$angle,pi,{pk},pi,$square]-OptionValue[LimitVar]^2ChainN[$angle,pk,{OptionValue[Reference]},pk,$square]);
	
		SpinorHelicity6D`Mom4DN[pkOut][Null][$up] = gamma * Mom4DN[pk][Null][$up];
		SpinorHelicity6D`Mom4DN[pkOut][$flat][$up] = Mom4DN[pkOut][Null][$up];
		SpinorHelicity6D`Mom4DN[pkOut][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[pkOut][Null][$up]];
		SpinorHelicity6D`Mom4DN[pkOut][$flat][$down] = Mom4DN[pkOut][Null][$down];
		Mom4DVectorNFromMatrix[pkOut][Null];
		Mom4DVectorNFromMatrix[pkOut][$flat];
		SpinorsFromMom4DN[pkOut,Normalization->pk];
		Unprotect[Extramass,Extramasstilde];
		SpinorHelicity6D`Extramass[pkOut] = 0;
		SpinorHelicity6D`Extramasstilde[pkOut] = 0;
		Protect[Extramass,Extramasstilde];
		SpinorHelicity6D`ExtramassN[pkOut] = 0;
		SpinorHelicity6D`ExtramasstildeN[pkOut] = 0;
	   
		SpinorHelicity6D`Mom4DN[piOut][Null][$up] = Mom4DN[pi][Null][$up] - OptionValue[LimitVar]^2 Mom4DN[OptionValue[Reference]][Null][$up] + (1-gamma) * Mom4DN[pk][Null][$up];
		SpinorHelicity6D`Mom4DN[piOut][$flat][$up] = Mom4DN[piOut][Null][$up];
		SpinorHelicity6D`Mom4DN[piOut][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[piOut][Null][$up]];
		SpinorHelicity6D`Mom4DN[piOut][$flat][$down] = Mom4DN[piOut][Null][$down];
		Mom4DVectorNFromMatrix[piOut][Null];
		Mom4DVectorNFromMatrix[piOut][$flat];
		SpinorsFromMom4DN[piOut,Normalization->pi];
		Unprotect[Extramass,Extramasstilde];
		SpinorHelicity6D`Extramass[piOut] = 0;
		SpinorHelicity6D`Extramasstilde[piOut] = 0;
		Protect[Extramass,Extramasstilde];
		SpinorHelicity6D`ExtramassN[piOut] = 0;
		SpinorHelicity6D`ExtramasstildeN[piOut] = 0;
	];
];




(* ::Subsection::Closed:: *)
(*GenSoftMasslessMassive *)

Options[GenSoftMasslessMassive] = {LimitVariable-> \[Epsilon]Soft,Reference->$NumericsReference};

GenSoftMasslessMassive[Rule[mom1_,momSoft1_],{Rule[mom2_,momRecoil2_],Rule[mom3_,momRecoil3_]}, opts : OptionsPattern[]]:=
Module[{a, q},
	q = OptionValue[Reference];

	Mom4DN[momSoft1][Null][$up] = OptionValue[LimitVariable]^2*Mom4DN[mom1][Null][$up];
	Mom4DN[momSoft1][$flat][$up] = Mom4DN[momSoft1][Null][$up];
	Mom4DN[momSoft1][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momSoft1][Null][$up]];
	Mom4DN[momSoft1][$flat][$down] = Mom4DN[momSoft1][Null][$down];
	Mom4DVectorNFromMatrix[momSoft1][Null];
	Mom4DVectorNFromMatrix[momSoft1][$flat];
	
	SpinorDotN[momSoft1][$lam][$up] = OptionValue[LimitVariable]*SpinorDotN[mom1][$lam][$up];
	SpinorUndotN[momSoft1][$lam][$up] = OptionValue[LimitVariable]*SpinorUndotN[mom1][$lam][$up];
	SpinorDotN[momSoft1][$lam][$down] = OptionValue[LimitVariable]*SpinorDotN[mom1][$lam][$down];
	SpinorUndotN[momSoft1][$lam][$down] = OptionValue[LimitVariable]*SpinorUndotN[mom1][$lam][$down];
	{SpinorDotN[momSoft1][$lam][1][Null],SpinorDotN[momSoft1][$lam][2][Null]} = SpinorDotN[momSoft1][$lam][$up];
	{SpinorUndotN[momSoft1][$lam][1][Null],SpinorUndotN[momSoft1][$lam][2][Null]} = SpinorUndotN[momSoft1][$lam][$up];
	{SpinorDotN[momSoft1][$lam][Null][1],SpinorDotN[momSoft1][$lam][Null][2]} = SpinorDotN[momSoft1][$lam][$down];
	{SpinorUndotN[momSoft1][$lam][Null][1],SpinorUndotN[momSoft1][$lam][Null][2]} = SpinorUndotN[momSoft1][$lam][$down];
	
	SpinorDotN[momSoft1][$mu][$up] = 0;
	SpinorUndotN[momSoft1][$mu][$up] = 0;
	SpinorDotN[momSoft1][$mu][$down] = 0;
	SpinorUndotN[momSoft1][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	Extramass[momSoft1] = 0;
	Extramasstilde[momSoft1] = 0;
	Protect[Extramass,Extramasstilde];
	ExtramassN[momSoft1] = 0;
	ExtramasstildeN[momSoft1] = 0;

	(* parameter required to match the masses of the scalar pair *)
	a = -(1-OptionValue[LimitVariable]^2)S6MomN[mom2,mom1]/
     	(S6MomN[mom2,q]+S6MomN[mom3,q]+(1-OptionValue[LimitVariable]^2)S6MomN[mom1,q]);

	Mom4DN[momRecoil2][Null][$up] =
     	Mom4DN[mom2][Null][$up]+(1-OptionValue[LimitVariable]^2)*Mom4DN[mom1][Null][$up]+a*Mom4DN[q][Null][$up];
	Mom4DN[momRecoil2][$flat][$up] = Mom4DN[momRecoil2][Null][$up];
	Mom4DN[momRecoil2][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momRecoil2][Null][$up]];
	Mom4DN[momRecoil2][$flat][$down] = Mom4DN[momRecoil2][Null][$down];
	Mom4DVectorNFromMatrix[momRecoil2][Null];

	ExtramassN[momRecoil2] = ExtramassN[mom2]ExtramasstildeN[mom2]-a S6MomN[mom3,q];
	ExtramasstildeN[momRecoil2] = 1;


	Mom4DN[momRecoil3][Null][$up] = Mom4DN[mom3][Null][$up] - a Mom4DN[q][Null][$up];
	Mom4DN[momRecoil3][$flat][$up] = Mom4DN[momRecoil3][Null][$up];
	Mom4DN[momRecoil3][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momRecoil3][Null][$up]];
	Mom4DN[momRecoil3][$flat][$down] = Mom4DN[momRecoil3][Null][$down];
	Mom4DVectorNFromMatrix[momRecoil3][Null];

	ExtramassN[momRecoil3] = -ExtramassN[momRecoil2];
	ExtramasstildeN[momRecoil3] = -ExtramasstildeN[momRecoil2];
	
];




(* ::Subsection::Closed:: *)
(*Collinear Splitting Functions*)


(* ::Subsubsection::Closed:: *)
(*Tree*)


SplittingFunctionTree[z_,a_,b_]:=
	<|"-"->
	      Association[{"-","-"}->0
			 ,{"+","+"}->1/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b],
			  {"+","-"}->-z^2/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b],
			  {"-","+"}->-(1-z)^2/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b]],
	  "+"->
	      Association[{"+","+"}->0,
			  {"-","-"}->-1/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b],
			  {"-","+"}->+z^2/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b],
			  {"+","-"}->(1-z)^2/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b]
	      ]
|>;


(* ::Subsubsection::Closed:: *)
(*One-Loop*)


cGamma[eps_]:=Exp[eps*EulerGamma]/2*Gamma[1+eps]Gamma[1-eps]^2/Gamma[1-2*eps];

rSSYM[z_,s_,mu_,eps_, cutoff_ : 1]:=cGamma[eps](-mu^2/s)^eps*1/eps^2(-((1-z)/z)^eps*Pi*eps/(Sin[Pi*eps])+Sum[2*eps^(2*m-1)PolyLog[2m-1,z/(z-1)],{m,1,cutoff}])

Options[rSQCDPP] = {DeltaHV->0,NFoNC->0};
rSQCDPP[z_,s_,mu_,eps_,opts:OptionsPattern[]]:=rSSYM[z,s,mu,eps]+cGamma[eps](-mu^2/s)^eps*2*z*(1-z)/((1-2*eps)(2-2*eps)(3-2*eps))(1-eps*OptionValue[DeltaHV]-OptionValue[NFoNC]);


Options[SplittingFunctionOneLoop] = {DeltaHV->0,NFoNC->0,MuVar->$RenormalizationScaleMu,EpsVar->$NumericsEpsilonDimReg};
				   
SplittingFunctionOneLoop[z_,a_,b_,opts :OptionsPattern[]]:=
	2*<|"-"->
		Association[
			{"-","-"}->cGamma[OptionValue[EpsVar]]Sqrt[z(1-z)]Spinoranglebracket[a,b]/SpinorSquareBracket[a,b]^2*
					2/((1-2*OptionValue[EpsVar])(2-2*OptionValue[EpsVar])(3-2*OptionValue[EpsVar]))*(-OptionValue[MuVar]^2/S[a,b])^OptionValue[EpsVar]*(1-OptionValue[EpsVar]*OptionValue[DeltaHV]-OptionValue[NFoNC]),
			{"+","+"}->SplittingFunctionTree[z,a,b]["-",{"+","+"}]*rSQCDPP[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"-","+"}->SplittingFunctionTree[z,a,b]["-",{"-","+"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"+","-"}->SplittingFunctionTree[z,a,b]["-",{"+","-"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]]
		],
	    "+"->
		Association[
			{"+","+"}->-cGamma[OptionValue[EpsVar]]Sqrt[z(1-z)]SpinorSquareBracket[a,b]/SpinorAngleBracket[a,b]^2*
					2/((1-2*OptionValue[EpsVar])(2-2*OptionValue[EpsVar])(3-2*OptionValue[EpsVar]))*(-OptionValue[MuVar]^2/S[a,b])^OptionValue[EpsVar]*(1-OptionValue[EpsVar]*OptionValue[DeltaHV]-OptionValue[NFoNC]),
			{"-","-"}->SplittingFunctionTree[z,a,b]["+",{"-","-"}]*rSQCDPP[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"+","-"}->SplittingFunctionTree[z,a,b]["+",{"+","-"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"-","+"}->SplittingFunctionTree[z,a,b]["+",{"-","+"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]]
		]
|>;


(* ::Subsection:: *)
(*Soft Functions*)


SoftTreeGGG["+"][a_,s_,b_]:=SpinorAngleBracket[a,b]/(SpinorAngleBracket[a,s]SpinorAngleBracket[s,b]);
SoftTreeGGG["-"][a_,s_,b_]:=-SpinorSquareBracket[a,b]/(SpinorSquareBracket[a,s]SpinorSquareBracket[s,b]);
SoftTreeSGG["+"][a_,s_,b_]:=-ATree["S","+","S"][a,s,pm@a]1/S6Mom[a,s]/.{Global`q->b};
SoftTreeSGG["-"][a_,s_,b_]:=-ATree["S","-","S"][a,s,pm@a]1/S6Mom[a,s]/.{Global`q->b};
SoftTreeGGS["+"][a_,s_,b_]:=-ATree["S","+","S"][pm@b,s,b]1/S6Mom[b,s]/.{Global`q->a};
SoftTreeGGS["-"][a_,s_,b_]:=-ATree["S","-","S"][pm@b,s,b]1/S6Mom[b,s]/.{Global`q->a};
SoftTreeSGS["+"][a_,s_,b_]:=-ATree["S","+","S"][a,s,pm@a]1/S6Mom[a,s]-ATree["S","+","S"][pm@b,s,b]1/S6Mom[b,s];
SoftTreeSGS["-"][a_,s_,b_]:=-ATree["S","-","S"][a,s,pm@a]1/S6Mom[a,s]-ATree["S","-","S"][pm@b,s,b]1/S6Mom[b,s];


(* ::Chapter:: *)
(*Postamble*)


Protect[Evaluate[protected]];
End[];
Protect[ "SpinorHelicity6D`Numeric`*" ];    (* protect exported symbols *)

EndPackage[ ];  (* end the package context *)
