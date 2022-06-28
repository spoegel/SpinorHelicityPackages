(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: LoopMomentum *)
(* :Context: LoopMomentum` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`Unitarity`LoopMomentum`",
	     {
		     "SpinorHelicity6D`",
		     "SpinorHelicity6D`Unitarity`",
		     "Utils`",
		     "NumSymb`",
		     "RationalSeries`",
		     "SpinorHelicity6D`SpinorFunctions`"
	     }
];
(* Exported symbols added here with SymbolName::usage *)


InitLoopMomentum::usage = "";
InitLoopMomentumVariant::usage = "";
FixBoxLoopMomentum6D::usage = "";
FixTriangleLoopMomentum6D::usage = "";
FixBubbleLoopMomentum6D::usage = "";
FixBubbleTriSubLoopMomentum6D::usage = "";

LoopMomentumLabels::usage = "";

LoopMomentumPostProcessingFirstPass::usage = "Option for first pass in LoopMomentumPostprocessing.";
LoopMomentumPostProcessingSecondPass::usage  = "Option for second pass in LoopMomentumPostprocessing.";

ReduceSqrts::usage = "Option for all loop-momentum definition functions. Reduces squareroots when their argument is a perfect square."

ComputeLoopMomentumSpinors::usage = "";
PostProcessLoopMomentum::usage = "";

IntermediateSimplificationFunction::usage = "";
IntermediateSimplificationFunction = Identity;

$BubbleFudgeFactor = 1;
$LoopMomentumSqrtsDebug::usage = "Contains the unique square roots defined in the last Box, Triangle or BubbleTri loop-momentum definition.";
$LoopMomentumSqrtsDebug = {};

$LoopMomentumReference::usage = "Default massless momentum used as reference for loop-momentum definitions.";
KillMasses[{$LoopMomentumReference}];
GenSpinors[
	{$LoopMomentumReference},
	FourD->{$LoopMomentumReference},
	RandomSpinors->True,
	Seed->235,
	ParameterRange->50,
	DisplaySpinors->False
];

Begin["`Private`"];


(* ::Section:: *)
(*Helper Functions*)

Options[LoopMomentumPostProcessing] =
    {LoopMomentumPostProcessingFirstPass -> Identity,LoopMomentumPostProcessingSecondPass -> Identity(*((SimplifyNumSymbs[#];#)&)*)};


LoopMomentumPostProcessing[expr_,variables_,opts : OptionsPattern[]] :=
    OptionValue[LoopMomentumPostProcessingSecondPass][
      CoalesceNumSymb@NumSymbReplace[OptionValue[LoopMomentumPostProcessingFirstPass][expr]/.sqrt->Sqrt,variables]
    ];


NumericRootReduce[expr_?NumericQ] := RootReduce@expr;
NumericRootReduce[expr_] := expr;


(* ::Section:: *)
(*Square Root Unique*)


SqrtUnique /: Power[SqrtUnique[args___],pow_] /; EvenQ[pow] := Power[SqrtUniqueArgument[args],pow/2];


CreateSqrtUnique[argExpr_]:=Module[{label=Unique[sqrt]},
  SqrtUniqueArgument[label] = argExpr;
  SqrtUnique[label]
];


(* ::Section:: *)
(* Loop Momentum *)


(* ::Subsection:: *)
(* Initialization of Loop Momentum Variables *)


Options[InitLoopMomentum]={LabelVariants->{Rule[$minus,lm],Rule[$plus,lp]},Reference->$LoopMomentumReference,
  MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,Assumptions->{},Rules->{},ComputeLoopMomentumSpinors->True,LargeMuLimit->False};

InitLoopMomentum[label_, subtractMom_List, subtractMass_List, opts : OptionsPattern[]]:=Module[{},
  If[OptionValue[LargeMuLimit],
    Unprotect[Extramass,Extramasstilde];
    Extramass[label]=OptionValue[MuVar]^2-Sum[Extramass[i],{i,subtractMass}];
    Extramasstilde[label]=OptionValue[MuTildeVar]^2-Sum[Extramasstilde[i],{i,subtractMass}];
    Protect[Extramass,Extramasstilde];

    ExtramassN[label]=OptionValue[MuVar]^2-Sum[ExtramassN[i],{i,subtractMass}];
    ExtramasstildeN[label]=OptionValue[MuTildeVar]^2-Sum[ExtramasstildeN[i],{i,subtractMass}];
    ,
    Unprotect[Extramass,Extramasstilde];
    Extramass[label]=OptionValue[MuVar]-Sum[Extramass[i],{i,subtractMass}];
    Extramasstilde[label]=OptionValue[MuTildeVar]-Sum[Extramasstilde[i],{i,subtractMass}];
    Protect[Extramass,Extramasstilde];

    ExtramassN[label]=OptionValue[MuVar]-Sum[ExtramassN[i],{i,subtractMass}];
    ExtramasstildeN[label]=OptionValue[MuTildeVar]-Sum[ExtramasstildeN[i],{i,subtractMass}];
  ];



  Do[
    InitLoopMomentumVariant[label,variant,subtractMom,subtractMass,FilterRules[{opts, Options[InitLoopMomentum]},Options[InitLoopMomentumVariant]]];
    ,{variant,OptionValue[LabelVariants]}];
];


Options[InitLoopMomentumVariant]={MuVar->$LoopMomentumMu,MuTildeVar->$LoopMomentumMuTilde,
  Reference->$LoopMomentumReference,Assumptions->{},Rules->{},ComputeLoopMomentumSpinors->True};

InitLoopMomentumVariant[label_,variantRule_Rule,subtractMom_List,subtractMass_List,opts:OptionsPattern[]]:=Module[{},


  Mom4DN[label[variantRule[[1]]]][Null][$down]=Mom4DN[variantRule[[2]]][Null][$down]-Sum[Mom4DN[i][Null][$down],{i,subtractMom}];
  Mom4DN[label[variantRule[[1]]]][Null][$up]=RaiseLowerMom4DNSpinorIndices[Mom4DN[label[variantRule[[1]]]][Null][$down]];


  Unprotect[Extramass,Extramasstilde];
  Extramass[label[variantRule[[1]]]]=Extramass[label]-Sum[Extramass[i],{i,subtractMass}];
  Extramasstilde[label[variantRule[[1]]]]=Extramasstilde[label]-Sum[Extramasstilde[i],{i,subtractMass}];
  Protect[Extramass,Extramasstilde];

  ExtramassN[label[variantRule[[1]]]]=ExtramassN[label]-Sum[ExtramassN[i],{i,subtractMass}];
  ExtramasstildeN[label[variantRule[[1]]]]=ExtramasstildeN[label]-Sum[ExtramasstildeN[i],{i,subtractMass}];

  Mom4DVectorNFromMatrix[label[variantRule[[1]]]][Null];


  Mom4DN[label[variantRule[[1]]]][$flat][$down]=Mom4DN[label[variantRule[[1]]]][Null][$down]-
      Mom4DN[OptionValue[Reference]][$flat][$down]*ExtramassN[label[variantRule[[1]]]]*ExtramasstildeN[label[variantRule[[1]]]]*
          1/ChainN[$angle,OptionValue[Reference],{UnderBar[label[variantRule[[1]]]]},OptionValue[Reference],$square];

  Mom4DN[label[variantRule[[1]]]][$flat][$up]=RaiseLowerMom4DNSpinorIndices[Mom4DN[label[variantRule[[1]]]][$flat][$down]];


  Mom4DVectorNFromMatrix[label[variantRule[[1]]]][$flat];


  If[OptionValue[ComputeLoopMomentumSpinors],
    SpinorsFromMom4DN[label[variantRule[[1]]],FilterRules[{opts,Options[InitLoopMomentumVariant]},Options[SpinorsFromMom4DN]]];
  ];
];


(* ::Subsection:: *)
(*Compute Gamma*)


Options[BoxGamma] =
{
	ReduceSqrts -> True
};

BoxGamma[P1P4_,s1_,s4_,opts:OptionsPattern[]]:=
	Module[{gamSqrt,res,eval},
	       
	       Which[
		       OptionValue[ReduceSqrts] && PossibleZeroQ[SampleEval[s1*s4],Method->"ExactAlgebraics"],
		       gamSqrt = P1P4;
		     ,
		       True
		     ,
		       gamSqrt = Sqrt[P1P4^2 - s1*s4];
	       ];
	       
	       AppendTo[$LoopMomentumSqrtsDebug,gamSqrt];

	       res = P1P4 + gamSqrt;

	       res
	];

Options[TriangleGamma] =
{
	ReduceSqrts -> True
};

TriangleGamma[P1P3_,s1_,s3_,opts:OptionsPattern[]]:=
	Module[{gamSqrt,res,eval,P1P2,s2},

	       s2 = s1 + s3 + 2 P1P3;
	       P1P2 = -s1 - P1P3;
	       
	       Which[
		       OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[s1*s3],Method->"ExactAlgebraics"],
		       gamSqrt = P1P3;
		       res = 2*P1P3;
		     ,
		       OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[s2],Method->"ExactAlgebraics"] && SampleEval[P1P2] > 0,
		       gamSqrt = P1P2;
		       res = -s1;
		     ,
		       OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[s2],Method->"ExactAlgebraics"] && SampleEval[P1P2] < 0,
		       gamSqrt = -P1P2;
		       res = -s1-2 P1P2;
		     ,
		       True
		     ,
		       gamSqrt = Sqrt[P1P3^2 - s1*s3];
		       res = P1P3 + gamSqrt;
	       ];
	       
	       AppendTo[$LoopMomentumSqrtsDebug,gamSqrt];
	       
	       res
	];


(* ::Subsection:: *)
(*Boxes (6D) *)


Options[FixBoxLoopMomentum6D]=
{
	LoopMomentumLabels->{l1,l2,l3,l4},
	Reference->$LoopMomentumReference,
	PrecisionGoal->MachinePrecision,
	MuVar->$LoopMomentumMu,
	MuTildeVar->$LoopMomentumMuTilde,
	LargeMuLimit->False,
	Rules->{},
	Assumptions->{},
	ComputeLoopMomentumSpinors->True,
	LoopMomentumPostProcessingFirstPass -> Identity,
	Gravity->False,
	ReduceSqrts->True
};

FixBoxLoopMomentum6D[{K1Ind_List, K2Ind_List, K3Ind_List, K4Ind_List}, opts : OptionsPattern[]] :=
	Module[{gam, c, d, Delta, DeltaPrime, ooapSeries, ooamSeries, am, ap,
		a0 = 0, a1 = 0, a2 = 0, a3 = 0, Lm, Lp, P1, P2 , P3, P4, P1b, P4b, P1P4,
		s1, s2, s4, mu2P2, aSqrt, TwoMulK1, TwoMulK4, MuK1, MuK4, MK1, MK3, MK4,
		MTildeK1, MTildeK3, MTildeK4, gamSqrt, aInnerSqrt, lm, lp, simpFunc
	       },
	       
	       Set4DMassive[{P1,P2,P3,P4,Lm,Lp}];

	       $LoopMomentumSqrtsDebug = {};

	       simpFunc = IntermediateSimplificationFunction;

	       Mom4DN[P1][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K1Ind}]//simpFunc;
	       Mom4DN[P2][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K2Ind}]//simpFunc;
	       Mom4DN[P3][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K3Ind}]//simpFunc;
	       Mom4DN[P4][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K4Ind}]//simpFunc;
	       Mom4DN[P1][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K1Ind}]//simpFunc;
	       Mom4DN[P2][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K2Ind}]//simpFunc;
	       Mom4DN[P3][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K3Ind}]//simpFunc;
	       Mom4DN[P4][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K4Ind}]//simpFunc;

	       Mom4DVectorNFromMatrix[P1][Null];
	       Mom4DVectorNFromMatrix[P2][Null];
	       Mom4DVectorNFromMatrix[P3][Null];
	       Mom4DVectorNFromMatrix[P4][Null];

	       s1 = (Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P1][Null][$down])//simpFunc;
	       s2 = (Mom4DVectorN[P2][Null][$up]) . (Mom4DVectorN[P2][Null][$down])//simpFunc;
	       s4 = (Mom4DVectorN[P4][Null][$up]) . (Mom4DVectorN[P4][Null][$down])//simpFunc;
	       mu2P2 = Expand[Sum[Extramass[i],{i,K2Ind}]*Sum[Extramasstilde[i],{i,K2Ind}]]//simpFunc;

	       MK1 = Sum[Extramass[i],{i,K1Ind}]//simpFunc;
	       MTildeK1 = Sum[Extramasstilde[i],{i,K1Ind}]//simpFunc;
	       MK3 = Sum[Extramass[i],{i,K3Ind}]//simpFunc;
	       MTildeK3 = Sum[Extramasstilde[i],{i,K3Ind}]//simpFunc;
	       MK4 = Sum[Extramass[i],{i,K4Ind}]//simpFunc;
	       MTildeK4 = Sum[Extramasstilde[i],{i,K4Ind}]//simpFunc;

	       TwoMulK1 =(OptionValue[MuVar]*MTildeK1 + OptionValue[MuTildeVar]*MK1)//simpFunc;
	       TwoMulK4 =(OptionValue[MuVar]*MTildeK4 + OptionValue[MuTildeVar]*MK4)//simpFunc;
	       MuK1 = Sum[Extramasstilde[i],{i,K1Ind}]*Sum[Extramass[i],{i,K1Ind}]//simpFunc;
	       MuK4 = Sum[Extramasstilde[i],{i,K4Ind}]*Sum[Extramass[i],{i,K4Ind}]//simpFunc;

	       P1P4=(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P4][Null][$down])//simpFunc;

	       gam = BoxGamma[P1P4,s1,s4,PassRules[FixBoxLoopMomentum6D,BoxGamma,opts]]//simpFunc;

	       Mom4DN[P1b][Null][$down] = gam/(gam^2 - s1*s4) (gam*Mom4DN[P1][Null][$down] - s1*Mom4DN[P4][Null][$down])//simpFunc;
	       Mom4DN[P4b][Null][$down] = gam/(gam^2 - s1*s4) (gam*Mom4DN[P4][Null][$down] - s4*Mom4DN[P1][Null][$down])//simpFunc;
	       Mom4DN[P1b][$flat][$down] = Mom4DN[P1b][Null][$down];
	       Mom4DN[P4b][$flat][$down] = Mom4DN[P4b][Null][$down];

	       Mom4DVectorNFromMatrix[P1b][$flat];
	       Mom4DVectorNFromMatrix[P4b][$flat];

	       SpinorsFromMom4DN[P1b,FilterRules[{opts,NumericalEvaluationFunction->SampleEval,Options[FixBoxLoopMomentum6D]},Options[SpinorsFromMom4DN]]];
	       SpinorsFromMom4DN[P4b,FilterRules[{opts,NumericalEvaluationFunction->SampleEval,Options[FixBoxLoopMomentum6D]},Options[SpinorsFromMom4DN]]];

	       

	       If[OptionValue[LargeMuLimit],
		  (* LargeMuLimit === True *)
		  If[And@@(PossibleZeroQ/@ToNum[{MK1,MTildeK1,MK3,MTildeK3,MK4,MTildeK4}]) == False,
		     Print["ERROR: The option LargeMuLimit->True in FixBoxLoopMomentum6D is only supported for purely 4D external kinematics. Aborting."];
		     Abort;
		  ];

		  c = -(gam + s1)*s4/(gam^2 - s1*s4);
		  d = (gam + s4)*s1/(gam^2 - s1*s4);

		  Delta = 2 c*(Mom4DVectorN[P1b][$flat][$up]) . (Mom4DVectorN[P2][Null][$down]) + 2*d*(Mom4DVectorN[P4b][$flat][$up]) . (Mom4DVectorN[P2][Null][$down]);
		  DeltaPrime = -Delta + 2*(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P2][Null][$down]) + s2;

		  aInnerSqrt = NumericRootReduce@Sqrt[TrMN[P4b, P2, P1b, P2]]/NumericRootReduce@Sqrt[gam];
		  AppendTo[$LoopMomentumSqrtsDebug,aInnerSqrt];

		  a0 = 2*aInnerSqrt;
		  a1 = (DeltaPrime^2 - 4*c*d*TrMN[P4b, P2, P1b, P2])/(4*aInnerSqrt);
		  If[OptionValue[Gravity],
		     a2 =  -(DeltaPrime^2 - 4*c*d*TrMN[P4b, P2, P1b, P2])^2/(64*(aInnerSqrt)^3);
		     a3 = (gam^3*aInnerSqrt*(DeltaPrime^2 - 4*c*d*TrMN[P4b, P2, P1b, P2])^3)/(512*TrMN[P4b, P2, P1b, P2]^3);
		  ];


		  (* If LargeMuLimit is true, we make the implicit redefinition that MuVar is Sqrt[MuVar], etc. to
           not have to deal with square roots *)
		  aSqrt = (OptionValue[MuVar]*OptionValue[MuTildeVar]) a0 + (OptionValue[MuVar]OptionValue[MuTildeVar])^(-1) a1 +
			  (OptionValue[MuVar]*OptionValue[MuTildeVar])^(-3) a2 + (OptionValue[MuVar]OptionValue[MuTildeVar])^(-5) a3;

		  {am,ap} = 1/(2*ChainN[$angle,P1b,{P2},P4b,$square])*(DeltaPrime + {-1,+1}* aSqrt);

		  If[OptionValue[Gravity]===True,
		     {ooamSeries,ooapSeries} = (2*ChainN[$angle,P1b,{P2},P4b,$square])*Plus@@(
			     {(DeltaPrime^6 - 5*a0*a1*DeltaPrime^4*pm^2 - 3*a0^2*(-2*a1^2 + a0*a2)*DeltaPrime^2*pm^4 - a0^3*(a1^3 - 2*a0*a1*a2 + a0^2*a3)*pm^6)/(a0^7*pm^7*OptionValue[MuTildeVar]^7*OptionValue[MuVar]^7),
			      -((DeltaPrime*(DeltaPrime^4 - 4*a0*a1*DeltaPrime^2*pm^2 + a0^2*(3*a1^2 - 2*a0*a2)*pm^4))/(a0^6*pm^6*OptionValue[MuTildeVar]^6*OptionValue[MuVar]^6)),
			      (DeltaPrime^4 - 3*a0*a1*DeltaPrime^2*pm^2 + a0^2*(a1^2 - a0*a2)*pm^4)/(a0^5*pm^5*OptionValue[MuTildeVar]^5*OptionValue[MuVar]^5),
			      -((DeltaPrime^3 - 2*a0*a1*DeltaPrime*pm^2)/(a0^4*pm^4*OptionValue[MuTildeVar]^4*OptionValue[MuVar]^4)),
			      (DeltaPrime^2 - a0*a1*pm^2)/(a0^3*pm^3*OptionValue[MuTildeVar]^3*OptionValue[MuVar]^3),
			      -(DeltaPrime/(a0^2*pm^2*OptionValue[MuTildeVar]^2*OptionValue[MuVar]^2)),
			      1/(a0*pm*OptionValue[MuTildeVar]*OptionValue[MuVar])}/.{pm->{-1,+1}});
		   ,
		     {ooamSeries,ooapSeries} = (2*ChainN[$angle,P1b,{P2},P4b,$square])*Plus@@(
			     {(DeltaPrime^2 - a0*a1*pm^2)/(a0^3*pm^3*OptionValue[MuTildeVar]^3*OptionValue[MuVar]^3),
			      -(DeltaPrime/(a0^2*pm^2*OptionValue[MuTildeVar]^2*OptionValue[MuVar]^2)),
			      1/(a0*pm*OptionValue[MuTildeVar]*OptionValue[MuVar])}/.{pm->{-1,+1}});
		  ];

		  Lm = 1/2*((am*ChainN[$angle,P1b,{#},P4b,$square] + (c*d*gam - OptionValue[MuVar]^2*OptionValue[MuTildeVar]^2)/gam * ooamSeries * ChainN[$angle,P4b, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
		  Lp = 1/2*((ap*ChainN[$angle,P1b,{#},P4b,$square] + (c*d*gam - OptionValue[MuVar]^2*OptionValue[MuTildeVar]^2)/gam * ooapSeries * ChainN[$angle,P4b, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
		  Lm += c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P4b][$flat][$up];
		  Lp += c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P4b][$flat][$up];
		,
		  (* LargeMuLimit === False*)
		  {c,d}= {(s4*(MuK1 - s1 - TwoMulK1) + gam*(MuK4 - s4 + TwoMulK4))/(gam^2 - s1*s4), (gam*(-MuK1 + s1 + TwoMulK1) - s1*(MuK4 - s4 + TwoMulK4))/(gam^2 - s1*s4)};

		  Delta = 2 c*(Mom4DVectorN[P1b][$flat][$up]) . (Mom4DVectorN[P2][Null][$down]) + 2*d*(Mom4DVectorN[P4b][$flat][$up]) . (Mom4DVectorN[P2][Null][$down]);
		  DeltaPrime = -Delta + 2*(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P2][Null][$down]) - Sum[Extramass[i]Extramasstilde[j]+Extramass[j]Extramasstilde[i],{i,K1Ind},{j,K2Ind}] +
			       s2 - mu2P2 + Sum[OptionValue[MuVar]Extramasstilde[i]+Extramass[i]OptionValue[MuTildeVar],{i,K2Ind}];

		  AppendTo[$LoopMomentumSqrtsDebug,Sqrt[DeltaPrime^2 - 4 (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/gam*TrMN[P4b,P2,P1b,P2]]];
		  
		  am = 1/(2*ChainN[$angle,P1b,{P2},P4b,$square])*(DeltaPrime - Sqrt[DeltaPrime^2 - 4 (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/gam*TrMN[P4b,P2,P1b,P2]]);
		  ap = 1/(2*ChainN[$angle,P1b,{P2},P4b,$square])*(DeltaPrime + Sqrt[DeltaPrime^2 - 4 (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/gam*TrMN[P4b,P2,P1b,P2]]);

		  Lm = 1/2*((am*ChainN[$angle,P1b,{#},P4b,$square] + (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(am*gam) ChainN[$angle,P4b, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
		  Lm += + c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P4b][$flat][$up];
		  Lp = 1/2*((ap*ChainN[$angle,P1b,{#},P4b,$square] + (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(ap*gam) ChainN[$angle,P4b, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
		  Lp += + c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P4b][$flat][$up];
	       ];

	       {Lm, Lp} = LoopMomentumPostProcessing[ToNum@{Lm,Lp},{OptionValue[MuVar],OptionValue[MuTildeVar]},
						     PassRules[FixBoxLoopMomentum6D,LoopMomentumPostProcessing,opts]];


	       Mom4DNFromVector[lm][$up][Lm];
	       Mom4DNFromVector[lp][$up][Lp];
	       Mom4DVectorNFromMatrix[lm][Null];
	       Mom4DVectorNFromMatrix[lp][Null];



	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[1]],{},{},
				LabelVariants->{$plus->lp,$minus->lm},PassRules[FixBoxLoopMomentum6D,InitLoopMomentum,opts]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[2]],K1Ind,K1Ind,
				LabelVariants->{$plus->lp,$minus->lm},PassRules[FixBoxLoopMomentum6D,InitLoopMomentum,opts]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[3]],Join[K1Ind,K2Ind],Join[K1Ind,K2Ind],
				LabelVariants->{$plus->lp,$minus->lm},PassRules[FixBoxLoopMomentum6D,InitLoopMomentum,opts]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[4]],Join[K1Ind,K2Ind,K3Ind],Join[K1Ind,K2Ind,K3Ind],
				LabelVariants->{$plus->lp,$minus->lm},PassRules[FixBoxLoopMomentum6D,InitLoopMomentum,opts]];

	       ClearKinematics[{P1,P2,P3,P4,P1b,P4b,lm,lp}];
	];


(* ::Subsection:: *)
(*Triangles (6D) *)


Options[FixTriangleLoopMomentum6D]=
{
	LoopMomentumLabels->{l1,l2,l3},
	Reference->$LoopMomentumReference,
	tVar->$LoopMomentumT,
	PrecisionGoal->MachinePrecision,
	MuVar->$LoopMomentumMu,
	MuTildeVar->$LoopMomentumMuTilde,
	InitializeReversed->False,
	Rules->{},
	Assumptions->{},
	ComputeLoopMomentumSpinors->True,
	LoopMomentumPostProcessingFirstPass -> Identity,
	PostProcessLoopMomentum ->True,
	ReduceSqrts->True
};

FixTriangleLoopMomentum6D[{K1Ind_List, K2Ind_List, K3Ind_List}, opts : OptionsPattern[]] :=
	Module[{gam, c, d, LUnstar, LStar, P1, P2, P3, P1b, P3b, P1P3, s1, s2, s3,
		TwoMulK1, TwoMulK3, MuK1,MuK3, P1bNorm, P3bNorm },
	       
	       Set4DMassive[{P1,P2,P3,LUnstar,LStar}];

	       $LoopMomentumSqrtsDebug = {};
	       
	       Mom4DN[P1][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K1Ind}];
	       Mom4DN[P2][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K2Ind}];
	       Mom4DN[P3][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K3Ind}];
	       Mom4DN[P1][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K1Ind}];
	       Mom4DN[P2][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K2Ind}];
	       Mom4DN[P3][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K3Ind}];



	       Mom4DVectorNFromMatrix[P1][Null];
	       Mom4DVectorNFromMatrix[P2][Null];
	       Mom4DVectorNFromMatrix[P3][Null];

	       s1 = (Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P1][Null][$down]);
	       s2 = (Mom4DVectorN[P2][Null][$up]) . (Mom4DVectorN[P2][Null][$down]);
	       s3 = (Mom4DVectorN[P3][Null][$up]) . (Mom4DVectorN[P3][Null][$down]);

	       TwoMulK1 =(OptionValue[MuVar]*Sum[Extramasstilde[i],{i,K1Ind}] + OptionValue[MuTildeVar]*Sum[Extramass[i],{i,K1Ind}]);
	       TwoMulK3 =(OptionValue[MuVar]*Sum[Extramasstilde[i],{i,K3Ind}] + OptionValue[MuTildeVar]*Sum[Extramass[i],{i,K3Ind}]);
	       MuK1 = Sum[Extramasstilde[i],{i,K1Ind}]*Sum[Extramass[i],{i,K1Ind}];
	       MuK3 = Sum[Extramasstilde[i],{i,K3Ind}]*Sum[Extramass[i],{i,K3Ind}];

	       P1P3=(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P3][Null][$down]);

	       gam = TriangleGamma[P1P3,s1,s3,PassRules[FixTriangleLoopMomentum6D,TriangleGamma,opts]];

	       Mom4DN[P1b][Null][$down] = gam/(gam^2 - s1*s3) (gam*Mom4DN[P1][Null][$down] - s1*Mom4DN[P3][Null][$down]);
	       Mom4DN[P3b][Null][$down] = gam/(gam^2 - s1*s3) (gam*Mom4DN[P3][Null][$down] - s3*Mom4DN[P1][Null][$down]);
	       Mom4DN[P1b][$flat][$down] = Mom4DN[P1b][Null][$down];
	       Mom4DN[P3b][$flat][$down] = Mom4DN[P3b][Null][$down];

	       Mom4DVectorNFromMatrix[P1b][$flat];
	       Mom4DVectorNFromMatrix[P3b][$flat];

	       If[Length@K1Ind === 1, P1bNorm = K1Ind[[1]], P1bNorm = OptionValue[Reference]];
	       If[Length@K3Ind === 1, P3bNorm = K3Ind[[1]], P3bNorm = OptionValue[Reference]];

	       SpinorsFromMom4DN[P1b,NumericalEvaluationFunction->SampleEval,Normalization->P1bNorm,PassRules[FixTriangleLoopMomentum6D,SpinorsFromMom4DN,opts]];
	       SpinorsFromMom4DN[P3b,NumericalEvaluationFunction->SampleEval,Normalization->P3bNorm,PassRules[FixTriangleLoopMomentum6D,SpinorsFromMom4DN,opts]];

	       {c,d} = {(s3*(MuK1 - s1 - TwoMulK1) + gam*(MuK3 - s3 + TwoMulK3))/(gam^2 - s1*s3),(gam*(-MuK1 + s1 + TwoMulK1) - s1*(MuK3 - s3 + TwoMulK3))/(gam^2 - s1*s3)};

	       LUnstar = 1/2*((OptionValue[tVar]*ChainN[$angle,P1b,{#},P3b,$square] + (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(OptionValue[tVar]*gam) ChainN[$angle,P3b, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
	       LUnstar += + c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P3b][$flat][$up];
	       LStar = 1/2*((OptionValue[tVar]*ChainN[$angle,P3b,{#},P1b,$square] + (c*d*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(OptionValue[tVar]*gam) ChainN[$angle,P1b, {#}, P3b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
	       LStar += + c*Mom4DVectorN[P1b][$flat][$up] + d*Mom4DVectorN[P3b][$flat][$up];

	       {LUnstar,LStar} = Switch[OptionValue[PostProcessLoopMomentum],
					True,
					LoopMomentumPostProcessing[ToNum@{LUnstar,LStar},{OptionValue[MuVar],OptionValue[MuTildeVar],OptionValue[tVar]},
								   PassRules[FixBoxLoopMomentum6D,LoopMomentumPostProcessing,opts]],
					False,
					ToNum@{LUnstar,LStar},
					_,
					OptionValue[PostProcessLoopMomentum]@ToNum@{LUnstar,LStar}
				 ];


	       Mom4DNFromVector[lus][$up][LUnstar];
	       Mom4DNFromVector[ls][$up][LStar];
	       Mom4DVectorNFromMatrix[lus][Null];
	       Mom4DVectorNFromMatrix[ls][Null];

	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[1]],{},{},LabelVariants->{$unstar->lus,$star->ls},FilterRules[{opts,Options[FixTriangleLoopMomentum6D]},Options[InitLoopMomentum]]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[2]],K1Ind,K1Ind,LabelVariants->{$unstar->lus,$star->ls},FilterRules[{opts,Options[FixTriangleLoopMomentum6D]},Options[InitLoopMomentum]]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[3]],Join[K1Ind,K2Ind],Join[K1Ind,K2Ind],LabelVariants->{$unstar->lus,$star->ls},FilterRules[{opts,Options[FixTriangleLoopMomentum6D]},Options[InitLoopMomentum]]];

	       ClearKinematics[{P1,P2,P3,P1b,P3b,lus,ls}];
	];



(* ::Subsection:: *)
(* Bubbles (6D) *)


Options[FixBubbleLoopMomentum6D]=
{
	LoopMomentumLabels->{l1,l2},
	Reference->$LoopMomentumReference,
	tVar->$LoopMomentumT,
	yVar->$LoopMomentumY,
	PrecisionGoal->MachinePrecision,
	MuVar->$LoopMomentumMu,
	MuTildeVar->$LoopMomentumMuTilde,
	Rules->{},
	Assumptions->{},
	ComputeLoopMomentumSpinors->True,
	LoopMomentumPostProcessingFirstPass -> Identity,
	Gravity->False,
	ReduceSqrts->True
};

FixBubbleLoopMomentum6D[{K1Ind_List, K2Ind_List}, Xi_, opts : OptionsPattern[]] :=
	Module[{gam, c, L, P1, P1b, s1, TwoMulK1, MuK1, l, P1bNorm},
	       Set4DMassive[{P1,L}];

	       $LoopMomentumSqrtsDebug = {};

	       Mom4DN[P1][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K1Ind}];
	       Mom4DN[P1][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K1Ind}];

	       Mom4DVectorNFromMatrix[P1][Null];

	       s1 = (Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P1][Null][$down]);
	       TwoMulK1 =(OptionValue[MuVar]*Sum[Extramasstilde[i],{i,K1Ind}] + OptionValue[MuTildeVar]*Sum[Extramass[i],{i,K1Ind}]);
	       MuK1 = Sum[Extramasstilde[i],{i,K1Ind}]*Sum[Extramass[i],{i,K1Ind}];
	       gam = 2*(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[Xi][$flat][$down]);


	       Mom4DN[P1b][Null][$down] = Mom4DN[P1][Null][$down] - s1/gam*Mom4DN[Xi][$flat][$down];
	       Mom4DN[P1b][$flat][$down] = Mom4DN[P1b][Null][$down];

	       P1bNorm = K1Ind[[-1]];

	       Mom4DVectorNFromMatrix[P1b][$flat];
	       SpinorsFromMom4DN[P1b,NumericalEvaluationFunction->SampleEval,Normalization->P1bNorm,PassRules[FixBubbleLoopMomentum6D,SpinorsFromMom4DN,opts]];

	       c = (s1(1-OptionValue[yVar])-MuK1+TwoMulK1)/gam;

	       L = 1/2*((OptionValue[tVar]*ChainN[$angle,P1b,{#},Xi,$square] + (OptionValue[yVar]*c*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(OptionValue[tVar]*gam) ChainN[$angle,Xi, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
	       L += OptionValue[yVar]*Mom4DVectorN[P1b][$flat][$up] + c*Mom4DVectorN[Xi][$flat][$up];

	       L= LoopMomentumPostProcessing[ToNum@L,{OptionValue[MuVar],OptionValue[MuTildeVar],OptionValue[tVar],OptionValue[yVar]},
					     PassRules[FixBoxLoopMomentum6D,LoopMomentumPostProcessing,opts]];

	       Mom4DNFromVector[l][$up][L];
	       Mom4DVectorNFromMatrix[l][Null];

	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[1]],{},{},LabelVariants->{Null->l},FilterRules[{opts,Options[FixBubbleLoopMomentum6D]},Options[InitLoopMomentum]]];
	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[2]],K1Ind,K1Ind,LabelVariants->{Null->l},FilterRules[{opts,Options[FixBubbleLoopMomentum6D]},Options[InitLoopMomentum]]];

	       ClearKinematics[{P1,P1b,l}];

	       If[OptionValue[Gravity],
		  Table[BubbleBubIntegral[i,s1,OptionValue[MuVar]*OptionValue[MuTildeVar]],{i,0,4}]
		,
		  Table[BubbleBubIntegral[i,s1,OptionValue[MuVar]*OptionValue[MuTildeVar]],{i,0,2}]
	       ]

(*{1,-1/2*(s1-MuK1+TwoMulK1)/s1,1/3*(((s1-MuK1+TwoMulK1)/s1)^2-OptionValue[MuVar]*OptionValue[MuTildeVar]/s1),s1}*)
	];


BubbleBubIntegral[n_,s1_,mu2_]:=Module[{dummyMu2},
  (*Global`yCoeff[n]**)1/(n+1)Sum[Binomial[n-i,i](-dummyMu2/s1)^i,{i,0,Floor[n/2]}]/.{dummyMu2->mu2}
];


(* ::Subsection:: *)
(* Bubbles Triangle Contribution (6D) *)


Options[FixBubbleTriSubLoopMomentum6D]=
{
	LoopMomentumLabels->{l1,l2,l3},
	Reference->$LoopMomentumReference,
	tVar->$LoopMomentumT,
	PrecisionGoal->100,
	LargeTLimit->True,
	MuVar->$LoopMomentumMu,
	MuTildeVar->$LoopMomentumMuTilde,
	ComputeLoopMomentumSpinors->True,
	ReplaceSqrts->False,
	Rules->{},
	Assumptions->{},
	LoopMomentumPostProcessingFirstPass -> Identity,
	Gravity->False,
	ReduceSqrts->True
};

FixBubbleTriSubLoopMomentum6D[{K1Ind_List, K2Ind_List}, KExtraIndCompl_, Xi_, opts : OptionsPattern[]] :=
	Module[{gam, c, Lp, Lm, P1, P1b, s1, TwoMulK1, MuK1, KExtraInd,
		PExtra, PExtrab, sPExtra, TwoMulPExtra, MuPExtra, P1PExtra,
		\[Alpha], y, c0, c1, c2, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0,
		sPExtraTilde, P1bPP1b,P1bPXi,XiPP1b,XiPXi, Delta, Ti, \[Alpha]Sqrt,
		simplifyAssume, sK1MKp},

	       $LoopMomentumSqrtsDebug = {};

	       KExtraInd = Complement[Join[K1Ind,K2Ind],KExtraIndCompl];
	       Set4DMassive[{P1,Lp,Lm,PExtra}];

	       (*Mom4DN[P1][Null][$down]=N[Sum[Mom4DN[i][Null][$down],{i,K1Ind}],OptionValue[PrecisionGoal]];
      Mom4DN[P1][Null][$up]=N[Sum[Mom4DN[i][Null][$up],{i,K1Ind}],OptionValue[PrecisionGoal]];
      Mom4DN[PExtra][Null][$down]=N[Sum[Mom4DN[i][Null][$down],{i,KExtraInd}],OptionValue[PrecisionGoal]];
      Mom4DN[PExtra][Null][$up]=N[Sum[Mom4DN[i][Null][$up],{i,KExtraInd}],OptionValue[PrecisionGoal]];*)

	       Mom4DN[P1][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,K1Ind}];
	       Mom4DN[P1][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,K1Ind}];
	       Mom4DN[PExtra][Null][$down]=Sum[Mom4DN[i][Null][$down],{i,KExtraInd}];
	       Mom4DN[PExtra][Null][$up]=Sum[Mom4DN[i][Null][$up],{i,KExtraInd}];

	       Mom4DVectorNFromMatrix[P1][Null];
	       Mom4DVectorNFromMatrix[PExtra][Null];

	       s1 = (Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[P1][Null][$down]);
	       sPExtra = (Mom4DVectorN[PExtra][Null][$up]) . (Mom4DVectorN[PExtra][Null][$down]);

	       TwoMulK1 =(OptionValue[MuVar]*Sum[Extramasstilde[i],{i,K1Ind}] + OptionValue[MuTildeVar]*Sum[Extramass[i],{i,K1Ind}]);
	       MuK1 = Sum[Extramasstilde[i],{i,K1Ind}]*Sum[Extramass[i],{i,K1Ind}];

	       TwoMulPExtra = OptionValue[MuVar]*Sum[Extramasstilde[i],{i,KExtraInd}] + OptionValue[MuTildeVar]*Sum[Extramass[i],{i,KExtraInd}];
	       MuPExtra = Sum[Extramasstilde[i],{i,KExtraInd}]*Sum[Extramass[i],{i,KExtraInd}];
	       P1PExtra = (Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[PExtra][Null][$down]);

	       gam = 2*(Mom4DVectorN[P1][Null][$up]) . (Mom4DVectorN[Xi][$flat][$down]);

	       Mom4DN[P1b][Null][$down] = (Mom4DN[P1][Null][$down] - s1/gam*Mom4DN[Xi][$flat][$down]);
	       Mom4DN[P1b][$flat][$down] = Mom4DN[P1b][Null][$down];

	       Mom4DVectorNFromMatrix[P1b][$flat];
	       SpinorsFromMom4DN[P1b,FilterRules[{opts,NumericalEvaluationFunction->SampleEval,Options[FixBubbleTriSubLoopMomentum6D]},Options[SpinorsFromMom4DN]]];

	       XiPXi=ChainN[$angle,Xi,{PExtra},Xi,$square];
	       P1bPXi=ChainN[$angle,P1b,{PExtra},Xi,$square];
	       XiPP1b=ChainN[$angle,Xi,{PExtra},P1b,$square];
	       P1bPP1b=ChainN[$angle,P1b,{PExtra},P1b,$square];

	       Delta=P1PExtra^2-s1*sPExtra;

	       If[OptionValue[Gravity],
		  Ti = Table[BubbleTriangleIntegral[i,s1,sPExtra,P1PExtra,XiPP1b,gam,OptionValue[MuVar]*OptionValue[MuTildeVar]],{i,1,6}];
		,
		  Ti = Table[BubbleTriangleIntegral[i,s1,sPExtra,P1PExtra,XiPP1b,gam,OptionValue[MuVar]*OptionValue[MuTildeVar]],{i,1,3}];
	       ];
	       (*
		 Ti[[1]]=-s1*XiPP1b/(2*gam*Delta);
		 Ti[[2]]=-3*s1*XiPP1b^2/(8*gam^2*Delta^2)(-s1*sPExtra+P1PExtra*s1);
		 Ti[[3]]=-XiPP1b^3/(48*gam^3*Delta^3)(15*s1^3*sPExtra^2+30*P1PExtra*s1^3*sPExtra+11*P1PExtra^2*s1^3+4*s1^4*sPExtra+16*OptionValue[MuVar]*OptionValue[MuTildeVar]*s1^2*Delta);
		*)

	       c2 = s1*XiPP1b;

	       If[OptionValue[LargeTLimit],
		  (* LargeTLimit === True *)

			  
		  sK1MKp = s1+sPExtra+P1bPP1b+s1/gam*XiPXi;

		  
		  Which[
			  (* Check for (K')^2=0 *)
			  OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[gam^2*s1*sPExtra],Method->"ExactAlgebraics"] && SampleEval[(gam*P1bPP1b+s1*XiPXi)] < 0,
			  \[Alpha] = (gam*P1bPP1b+s1*XiPXi)^2;
			  \[Alpha]Sqrt = -(gam*P1bPP1b+s1*XiPXi);
			,
			  (* Check for (K')^2=0 *)
			  OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[gam^2*s1*sPExtra],Method->"ExactAlgebraics"] && SampleEval[(gam*P1bPP1b+s1*XiPXi)] > 0,
			  \[Alpha] = (gam*P1bPP1b+s1*XiPXi)^2;
			  \[Alpha]Sqrt = (gam*P1bPP1b+s1*XiPXi);
			,
			  (* Check for (K-K')^2=0 *)
			  OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[sK1MKp],Method->"ExactAlgebraics"] && SampleEval[gam*P1bPP1b+s1*(2*gam+XiPXi)] < 0,
			  \[Alpha] = (gam*P1bPP1b+s1*(2*gam+XiPXi))^2;
			  \[Alpha]Sqrt = -(gam*P1bPP1b+s1*(2*gam+XiPXi));
			,
			  (* Check for (K-K')^2=0 *)
			  OptionValue[ReduceSqrts] === True && PossibleZeroQ[SampleEval[sK1MKp],Method->"ExactAlgebraics"] && SampleEval[gam*P1bPP1b+s1*(-2*gam+XiPXi)] > 0,
			  \[Alpha] = (gam*P1bPP1b+s1*(2*gam+XiPXi))^2;
			  \[Alpha]Sqrt = (gam*P1bPP1b+s1*(2*gam+XiPXi));
			,
			  True
			,
			  \[Alpha] = (gam*P1bPP1b+s1*XiPXi)^2-4*gam^2*s1*sPExtra;
			  \[Alpha]Sqrt = Sqrt[\[Alpha]];
		  ];
		  
		  sPExtraTilde=+s1+TwoMulK1-MuK1;

		  AppendTo[$LoopMomentumSqrtsDebug,\[Alpha]Sqrt];

		  a1 = (gam*P1bPP1b - s1*XiPXi + {+1,-1}*\[Alpha]Sqrt)/(2*c2);
		  a2 = (-2*c2*{+1,-1}*(gam*(MuPExtra - sPExtra + TwoMulPExtra) + (MuK1 - s1 - TwoMulK1)*XiPXi) -
			(MuK1 - s1 - TwoMulK1)*XiPP1b*(gam*P1bPP1b*{+1,-1} - {+1,-1}*s1*XiPXi + \[Alpha]Sqrt))/(2*c2*\[Alpha]Sqrt);
		  a3 = ({+1,-1}*(-((MuK1 - s1 - TwoMulK1)*XiPP1b*(gam*P1bPP1b - s1*XiPXi) + 2*c2*(gam*(MuPExtra - sPExtra + TwoMulPExtra) + (MuK1 - s1 - TwoMulK1)*XiPXi))^2 +
				 XiPP1b*(-4*c2*OptionValue[MuVar]*OptionValue[MuTildeVar] + (-MuK1 + s1 + TwoMulK1)^2*XiPP1b)*\[Alpha]))/(4*c2*\[Alpha]Sqrt^3);
		  a4 = -({+1,-1}*XiPP1b^2*(gam*(-2*MuPExtra*s1 + 2*s1*sPExtra + P1bPP1b*sPExtraTilde - 2*s1*TwoMulPExtra) + s1*sPExtraTilde*XiPXi)*
			 (-(gam*(2*MuPExtra*s1 - 2*s1*sPExtra - P1bPP1b*sPExtraTilde + 2*s1*TwoMulPExtra) - s1*sPExtraTilde*XiPXi)^2 + (-4*OptionValue[MuVar]*OptionValue[MuTildeVar]*s1 + sPExtraTilde^2)*\[Alpha]))/(4*s1*\[Alpha]Sqrt^5);
		  a5 = ({+1,-1}*XiPP1b^3*(-5*(gam*(-2*MuPExtra*s1 + 2*s1*sPExtra + P1bPP1b*sPExtraTilde - 2*s1*TwoMulPExtra) + s1*sPExtraTilde*XiPXi)^4 -
					  6*(4*OptionValue[MuVar]*OptionValue[MuTildeVar]*s1 - sPExtraTilde^2)*(gam*(-2*MuPExtra*s1 + 2*s1*sPExtra + P1bPP1b*sPExtraTilde - 2*s1*TwoMulPExtra) + s1*sPExtraTilde*XiPXi)^2*\[Alpha] -
					  (-4*OptionValue[MuVar]*OptionValue[MuTildeVar]*s1 + sPExtraTilde^2)^2*\[Alpha]^2))/(16*s1*\[Alpha]Sqrt^7);

		  y = OptionValue[tVar] * a1 + a2 + 1/OptionValue[tVar] * a3 + 1/OptionValue[tVar]^2 * a4 + 1/OptionValue[tVar]^3 * a5;
		,
		  (* LargeTLimit === False *)

		  c1 = gam*P1bPP1b*OptionValue[tVar]-s1*XiPXi*OptionValue[tVar]+(s1-MuK1+TwoMulK1)*XiPP1b;
		  c0 = gam*OptionValue[tVar](sPExtra-MuPExtra-TwoMulPExtra)+OptionValue[tVar](s1-MuK1+TwoMulK1)*XiPXi+OptionValue[tVar]^2*gam*P1bPXi-OptionValue[MuVar]*OptionValue[MuTildeVar]*XiPP1b;

		  y = 1/(2*c2)*(c1+{1,-1}*Sqrt[c1^2+4*c0*c2]);
		  AppendTo[$LoopMomentumSqrtsDebug,Sqrt[c1^2+4*c0*c2]];
	       ];
	       
	       c = (s1(1-y)-MuK1+TwoMulK1)/gam;

	       Lm = 1/2*((OptionValue[tVar]*ChainN[$angle,P1b,{#},Xi,$square] + (y[[1]]*c[[1]]*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(OptionValue[tVar]*gam) ChainN[$angle,Xi, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
	       Lm += y[[1]]*Mom4DVectorN[P1b][$flat][$up] + c[[1]]*Mom4DVectorN[Xi][$flat][$up];

	       Lp = 1/2*((OptionValue[tVar]*ChainN[$angle,P1b,{#},Xi,$square] + (y[[2]]*c[[2]]*gam - OptionValue[MuVar]*OptionValue[MuTildeVar])/(OptionValue[tVar]*gam) ChainN[$angle,Xi, {#}, P1b,$square]) & /@ {PauliMatrixSpinor[0, $up], PauliMatrixSpinor[1, $up], PauliMatrixSpinor[2, $up], PauliMatrixSpinor[3, $up]});
	       Lp += y[[2]]*Mom4DVectorN[P1b][$flat][$up] + c[[2]]*Mom4DVectorN[Xi][$flat][$up];


	       {Lm,Lp} = LoopMomentumPostProcessing[ToNum@{Lm,Lp},{OptionValue[MuVar],OptionValue[MuTildeVar],OptionValue[tVar]},
						    PassRules[FixBoxLoopMomentum6D,LoopMomentumPostProcessing,opts]];


	       Mom4DNFromVector[lm][$up][Lm];
	       Mom4DNFromVector[lp][$up][Lp];
	       Mom4DVectorNFromMatrix[lm][Null];
	       Mom4DVectorNFromMatrix[lp][Null];

	       InitLoopMomentum[OptionValue[LoopMomentumLabels][[1]],{},{},LabelVariants->{$minus->lm,$plus->lp},FilterRules[{opts,Options[FixBubbleTriSubLoopMomentum6D]},Options[InitLoopMomentum]]];
	       Switch[Position[{K1Ind,K2Ind},Last@KExtraIndCompl][[1,1]],
		      1,
		      InitLoopMomentum[OptionValue[LoopMomentumLabels][[2]],KExtraIndCompl,KExtraIndCompl,LabelVariants->{$minus->lm,$plus->lp},FilterRules[{opts,Options[FixBubbleTriSubLoopMomentum6D]},Options[InitLoopMomentum]]];
		      InitLoopMomentum[OptionValue[LoopMomentumLabels][[3]],K1Ind,K1Ind,LabelVariants->{$minus->lm,$plus->lp},FilterRules[{opts,Options[FixBubbleTriSubLoopMomentum6D]},Options[InitLoopMomentum]]];,
				      2,
		      InitLoopMomentum[OptionValue[LoopMomentumLabels][[2]],K1Ind,K1Ind,LabelVariants->{$minus->lm,$plus->lp},FilterRules[{opts,Options[FixBubbleTriSubLoopMomentum6D]},Options[InitLoopMomentum]]];
		      InitLoopMomentum[OptionValue[LoopMomentumLabels][[3]],KExtraIndCompl,KExtraIndCompl,LabelVariants->{$minus->lm,$plus->lp},FilterRules[{opts,Options[FixBubbleTriSubLoopMomentum6D]},Options[InitLoopMomentum]]];,
				      _,
				      Print["Error: Last momentum of KExtraInd cannot be found in K1Ind or K2Ind."];
		      Abort[];
	       ];

	       ClearKinematics[{P1,P1b,PExtra,lm,lp}];

	       ToNum@Ti
	];



(* adapting my result to Badgers, i.e. no (-1)^n here. Not sure whether that's correct *)
BubbleTriangleIntegral[n_,s1_,s2_,p1p2_,chaink1k2chi_,gamma_,mu2_]:=(*Global`tCoeff[n]**)(-1)^(n+1)*
    chaink1k2chi^n/(gamma)^n*BubbleTrianglePVCoeff[n,s1,s2,p1p2,mu2];


BubbleTrianglePVCoeff[1,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  -s1/(2*delta)/.{delta->p1p2^2-s1*s2}
];


BubbleTrianglePVCoeff[2,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  -(3*s1^2*(p1p2 + s2))/(8*delta^2)/.{delta->p1p2^2-s1*s2}
];


BubbleTrianglePVCoeff[3,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  -(s1^2*(p1p2^2*(16*mu2 + 11*s1) + 30*p1p2*s1*s2 + s1*s2*(-16*mu2 + 4*s1 + 15*s2)))/(48*delta^3)/.{delta->p1p2^2-s1*s2}
];


BubbleTrianglePVCoeff[4,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  (-5*s1^3*(p1p2 + s2)*(2*p1p2^2*(22*mu2 + 5*s1) + 42*p1p2*s1*s2 + s1*s2*(-44*mu2 + 11*s1 + 21*s2)))/(384*delta^4)/.{delta->p1p2^2-s1*s2}
];


BubbleTrianglePVCoeff[5,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  (s1^3*(-1024*delta^2*mu2^2 + s1*(-4*mu2*(p1p2^2 - s1*s2)*(607*p1p2^2 + 1470*p1p2*s2 + 128*s1*s2 + 735*s2^2) -
      s1*(274*p1p2^4 + p1p2^2*(2310*p1p2 + 607*s1)*s2 + (4935*p1p2^2 + 1470*p1p2*s1 + 64*s1^2)*s2^2 + 105*(36*p1p2 + 7*s1)*s2^3 + 945*s2^4))))/(3840*delta^5)/.{delta->p1p2^2-s1*s2}
];


BubbleTrianglePVCoeff[6,s1_,s2_,p1p2_,mu2_]:=Module[{delta},
  (-7*s1^4*(p1p2 + s2)*(528*delta^2*mu2^2 + s1^2*(28*p1p2^4 + 8*p1p2^2*(40*p1p2 + 13*s1)*s2 + (820*p1p2^2 + 340*p1p2*s1 + 33*s1^2)*s2^2 + 10*(66*p1p2 + 17*s1)*s2^3 +
      165*s2^4) + 8*delta*mu2*s1*(52*p1p2^2 + 170*p1p2*s2 + s2*(33*s1 + 85*s2))))/(5120*delta^6)/.{delta->p1p2^2-s1*s2}
];

End[]; (* `Private` *)
Unprotect[IntermediateSimplificationFunction,$BubbleFudgeFactor];
EndPackage[]
