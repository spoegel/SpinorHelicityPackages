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

(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: ScalarKinematics *)
(* :Context: ScalarKinematics` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-06-10 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Spinors`ScalarKinematics`",{"Spinors`","Spinors`NRational`","Spinors`Kinematics`"}];
(* Exported symbols added here with SymbolName::usage *)

lp::usage = "";
SumM::usage = "";

ATreeSAM::usage = "";
ATreeSAMContact::usage = "";
ATreeSAMGluon::usage = "";

SetupSAMFourScalarMomentumKinematics::usage = "";
SetupSAMTwoScalarMomentumKinematics::usage = "";

Begin["`Private`"];

AppendTo[$ContextPath,"Spinors`Private`"];

Unprotect/@$SpinorsFunctions;

(* :Section: *)
(* lp *)

NRationalEval[lp[x_,y_]]:=2 NRationalEval[MP[x,y]];
NRationalEval[lp[x:Repeated[_,{3,Infinity}]]]:= 2*NRationalEval/@MP@@@Subsets[{x},{2}]//Apply[Plus];

AppendTo[$numericHeads,lp];


(* :Section: *)
(* Setup scalar momenta *)
(*RemoveNumSpinors[a_]:=(NumSpinorList:=Evaluate[Alternatives@@DeleteCases[List@@NumSpinorList,*)
(*  Alternatives@@Flatten@{a}]]);*)
(*RemoveNumVectors[a_]:=(NumVectorList:=Evaluate[Alternatives@@DeleteCases[List@@NumVectorList,*)
(*  Alternatives@@Flatten@{a}]]);*)


Options[SetupSAMFourScalarMomentumKinematics] = {Quiet->False};

SetupSAMFourScalarMomentumKinematics[scalar1Labels : {_, _}, scalar2Labels : {_, _},
  masslessLabels_List, kinematicPoint_Association, opts:OptionsPattern[]] /;
    Length@masslessLabels ===
        Length@kinematicPoint["Massless"] :=
    Module[{$blocked},

      $blocked = If[OptionValue[Quiet]===True,Print,$blocked];

      Block[{Print},
        RemoveNumSpinors[Join[scalar1Labels,scalar2Labels,masslessLabels]];
        RemoveNumVectors[Join[scalar1Labels,scalar2Labels,masslessLabels]];
        Map[UndeclareLVector[#]&, Join[scalar1Labels,scalar2Labels,masslessLabels]];
        Map[UndeclareSpinor[#]&, Join[scalar1Labels,scalar2Labels,masslessLabels]];
      ];

      Block[{#},

        MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars1", First@#2]]]&,
          scalar1Labels];
        MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars2", First@#2]]]&,
          scalar2Labels];

        MapIndexed[DeclareSpinorMomentum[#1, List /@ kinematicPoint[["MasslessSpinors", First@#2, 1]],
          {kinematicPoint[["MasslessSpinors", First@#2, 2]]}]&, masslessLabels];

      ]&@$blocked;
    ];

Options[SetupSAMTwoScalarMomentumKinematics] = {Quiet->False};

SetupSAMTwoScalarMomentumKinematics[scalarLabels : {_, _}, masslessLabels_List,
  kinematicPoint_Association, opts:OptionsPattern[]] /; Length@masslessLabels === Length@kinematicPoint["Massless"]:=
    Module[{$blocked},

      $blocked = If[OptionValue[Quiet]===True,Print,$blocked];

      Block[{Print},
        RemoveNumSpinors[Join[scalarLabels,masslessLabels]];
        RemoveNumVectors[Join[scalarLabels,masslessLabels]];
        Map[UndeclareLVector[#]&, Join[scalarLabels,masslessLabels]];
        Map[UndeclareSpinor[#]&, Join[scalarLabels,masslessLabels]];
      ];

      Block[{#},

        MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars", First@#2]]]&, scalarLabels];

        MapIndexed[DeclareSpinorMomentum[#1, List /@ kinematicPoint[["MasslessSpinors", First@#2, 1]],
          {kinematicPoint[["MasslessSpinors", First@#2, 2]]}]&, masslessLabels];

      ]&@$blocked;
    ];


(* :Section: *)
(* SumM *)

SetAttributes[SumM,Orderless];

SumM /: MakeBoxes[SumM[y__],StandardForm|TraditionalForm] := TemplateBox[ToBoxes/@{y},"SumM",
  DisplayFunction->(RowBox[{"(",TemplateSlotSequence[1,"+"],")"}]&),
  InterpretationFunction->(RowBox[{"SumM","[",TemplateSlotSequence[1,","],"]"}]&)
];

SumM /: NumSMatrixQ[SumM[x__]] := AllTrue[{x},NumSMatrixQ];
SumM /: Explicit[SumM[x__]] := Plus@@Explicit/@{x};

N[Spaa[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ],p_:$MachinePrecision]:=
    (N[UbarSpa[a],p].(Dot@@((Explicit/@{dm})/.SumM->Plus)).N[USpa[b],p])[[1,1]];
N[Spab[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ],p_:$MachinePrecision]:=
    (N[UbarSpa[a].(Dot@@((Explicit/@{dm})/.SumM->Plus)).USpb[b],p])[[1,1]];
N[Spba[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ],p_:$MachinePrecision]:=
    (N[UbarSpb[a].(Dot@@((Explicit/@{dm})/.SumM->Plus)).USpa[b],p])[[1,1]];
N[Spbb[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ],p_:$MachinePrecision]:=
    (N[UbarSpb[a].(Dot@@((Explicit/@{dm})/.SumM->Plus)).USpb[b],p])[[1,1]];


NRationalEval[Spaa[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumCLa[a].(Dot@@(NRationalEval/@(alternate[Sm2,CSm2][{dm}]/.SumM->Plus))).NumLa[b])//First//First;

NRationalEval[Spab[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumCLa[a].(Dot@@(NRationalEval/@(alternate[Sm2,CSm2][{dm}]/.SumM->Plus))).NumCLat[b])//First//First;

NRationalEval[Spba[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumLat[a].(Dot@@(NRationalEval/@(alternate[CSm2,Sm2][{dm}]/.SumM->Plus))).NumLa[b]) //First//First;

NRationalEval[Spbb[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumLat[a].(Dot@@(NRationalEval/@(alternate[CSm2,Sm2][{dm}]/.SumM->Plus))).NumCLat[b])//First//First;

$ContextPath = DeleteCases[$ContextPath,"Spinors`Private`"];

Protect/@$SpinorsFunctions;

End[]; (* `Private` *)

EndPackage[]