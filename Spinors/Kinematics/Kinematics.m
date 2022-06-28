(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: Kinematics *)
(* :Context: Kinematics` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-10-01 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Spinors`Kinematics`",
  {
    "Spinors`",
    "MomentumTwistorDefinitions`",
    "Spinors`NRational`"
  }
];
(* Exported symbols added here with SymbolName::usage *)

GenRationalMomenta::usage = "Defines rational valued momenta. Use NRational for evaluation.";

GenParametrizedMomenta::usage = "Defines momenta in terms of momentum twistor variables.
  Use NRational for evaluation.";
MomentumTwistorParameterSolution::usage = "Returns replacement rules to translate twistor
  parameters into spinorial expressions";

GenMomentaFromTwistors::usage = "Defines spinors from a provided twistor matrix.
  GenMomentaFromTwistors[labels, matrix]";
GenMomentaFromTwistors::missmatch = "Number of labels does not match number of provided twistors.";

RemoveNumSpinors::usage = "";
RemoveNumVectors::usage = "";


Begin["`Private`"];

AppendTo[$ContextPath,"Spinors`Private`"];
Unprotect/@$SpinorsFunctions;


(* :Section: *)
(* Setup scalar momenta *)
RemoveNumSpinors[a_]:=(NumSpinorList:=Evaluate[Alternatives@@DeleteCases[List@@NumSpinorList,
  Alternatives@@Flatten@{a}]]);
RemoveNumVectors[a_]:=(NumVectorList:=Evaluate[Alternatives@@DeleteCases[List@@NumVectorList,
  Alternatives@@Flatten@{a}]]);


(* :Section: *)
(* Twistor kinematics *)

Options[GenMomentaFromTwistors] = {Quiet->False};

GenMomentaFromTwistors[labelsBare_List, ZMatrix_List,opts:OptionsPattern[]] /;
    Length@labelsBare === Dimensions[ZMatrix][[2]] :=
    Module[{$blocked,n=Length@labelsBare,lambda,lambdaT,labels},

      labels = labelsBare /. {i_Integer?Positive :> Sp[i]};

      $blocked = If[OptionValue[Quiet]===True,Print,$blocked];

      {lambda,lambdaT} = LambdaLambdaTildeFromTwistor[ZMatrix];

      Block[{Print},
        RemoveNumSpinors[labels];
        RemoveNumVectors[labels];
        UndeclareLVector /@ labels;
        UndeclareSpinor /@ labels;
      ];

      Block[{#},
        MapIndexed[
          DeclareSpinorMomentum[#1,
            List /@ lambda[[First@#2]],
            {lambdaT[[First@#2]]}
          ]&,
          labels
        ];
      ]&@$blocked;
    ];

GenMomentaFromTwistors[___] /; Message[GenMomentaFromTwistors::missmatch] := Null;



(* :Section: *)
(* Rational kinematics *)

Options[GenRationalMomenta] = Join[Options[GenMomentaFromTwistors],{ParameterRange->1000}];

GenRationalMomenta[labelsBare_List,opts:OptionsPattern[]] :=
    GenMomentaFromTwistors[
      labelsBare,
      RandomInteger[{1,OptionValue[ParameterRange]},{4,Length@labelsBare}],
      FilterRules[{opts},Options[GenMomentaFromTwistors]]
    ];



(* :Section: *)
(* Parametrized kinematics *)

Options[GenParametrizedMomenta] = Options[GenMomentaFromTwistors];

GenParametrizedMomenta[labelsBare_List,opts:OptionsPattern[]] :=
    GenMomentaFromTwistors[
      labelsBare,
      MomentumTwistorMatrixParametrized[Length@labelsBare],
      FilterRules[{opts},Options[GenMomentaFromTwistors]]
    ];


MomentumTwistorParameterSolution[labelsBare_] := Module[
  {labels,n=Length@labelsBare},

  labels = labelsBare /. {i_Integer?Positive :> Sp[i]};

  MapAt[(#/.MapThread[Rule,{Sp/@(Range@n),labels}])&,
    MomentumTwistorParameterSolutionTemplate[n]/.{
      $momTwistorInv -> s,
      $momTwistorSpaaProd -> Spaa,
      $momTwistorSpbbProd -> Spbb,
      $momTwistorSpaaChain -> Spaa,
      $momTwistorSpabChain -> Spab,
      $momTwistorSpbaChain -> Spba,
      $momTwistorSpbbChain -> Spbb
    },
    {All,2}
  ]
];


  Protect/@$SpinorsFunctions;

  End[]; (* `Private` *)

  EndPackage[];