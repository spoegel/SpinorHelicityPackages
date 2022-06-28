(* Mathematica Source File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-05-24 *)




Begin["Spinors`"];

lp::usage = "";
SumM::usage = "";

ATreeSAM::usage = "";
ATreeSAMContact::usage = "";

NRational::usage = "";

SetupSAMFourScalarMomentumKinematics::usage = "";
SetupSAMTwoScalarMomentumKinematics::usage = "";

Begin["`Private`"];

Unprotect/@$SpinorsFunctions;

RemoveNumSpinors[a_]:=(NumSpinorList:=Evaluate[Alternatives@@DeleteCases[List@@NumSpinorList,
  Alternatives@@Flatten@{a}]]);
RemoveNumVectors[a_]:=(NumVectorList:=Evaluate[Alternatives@@DeleteCases[List@@NumVectorList,
  Alternatives@@Flatten@{a}]]);


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
(* lp *)
SetAttributes[lp,Orderless];

(* :Section: *)
(* NRationalEval *)

NRationalEval[x_?NumericQ] := NRationalEval[x] =  x;

NRationalEval[Dot[Lat[b_?NumSpinorQ],CLat[a_?NumSpinorQ]]]:=NRationalEval[(NumLat[b].NumCLat[a])[[1,1]]];

NRationalEval[d:Dot[Lat[b_?NumSpinorQ],sms__?NumSm2Q,La[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;RightSM2Order[{sms}]&&OddQ[Length[{sms}]];

NRationalEval[d:Dot[Lat[b_?NumSpinorQ],sms__?NumSm2Q,CLat[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;RightSM2Order[{sms}]&&EvenQ[Length[{sms}]];

NRationalEval[Dot[CLa[b_?NumSpinorQ],la:La[a_?NumSpinorQ]]]:=(NRationalEval[NumCLa[b].NumLa[a]])[[1,1]];

NRationalEval[d:Dot[CLa[b_?NumSpinorQ],sms__?NumSm2Q,CLat[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;RightCSM2Order[{sms}]&&OddQ[Length[{sms}]];

NRationalEval[d:Dot[CLa[b_?NumSpinorQ],sms__?NumSm2Q,La[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;RightCSM2Order[{sms}]&&EvenQ[Length[{sms}]];

NRationalEval[Dot[UbarSpa[b_?NumSpinorQ],USpa[a_?NumSpinorQ]]]:=NRationalEval[(NumCLa[b].NumLa[a])[[1,1]]];

NRationalEval[Dot[UbarSpb[b_?NumSpinorQ],USpb[a_?NumSpinorQ]]]:=NRationalEval[(NumLat[b].NumCLat[a])[[1,1]]];

NRationalEval[d:Dot[UbarSpa[b_?NumSpinorQ],sms__?NumSm4Q,USpa[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;EvenQ[Length[{sms}]];

NRationalEval[d:Dot[UbarSpa[b_?NumSpinorQ],sms__?NumSm4Q,USpb[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;OddQ[Length[{sms}]];

NRationalEval[d:Dot[UbarSpb[b_?NumSpinorQ],sms__?NumSm4Q,USpb[a_?NumSpinorQ]]]:=(NRationalEval[#]&/@d)[[1,1]]/;EvenQ[Length[{sms}]];


NRationalEval[SmBA2[b_?NumSpinorQ,a_?NumSpinorQ]]:=NRationalEval[NumCLat[b].NumCLa[a]];
NRationalEval[CSmBA2[b_?NumSpinorQ,a_?NumSpinorQ]]:=NRationalEval[NumLa[a].NumLat[b]];
NRationalEval[SmBA4[b_?NumSpinorQ,a_?NumSpinorQ]]:=NRationalEval[USpb[b].UbarSpa[a]+USpa[a].UbarSpb[b]];


Sm2/:NRationalEval[Sm2[a:(_?NumSpinorQ|_?NumVectorQ|_?NumSMatrixQ)]]:=NRationalEval[SlashM2D[a]];


CSm2/:NRationalEval[CSm2[a:(_?NumSpinorQ|_?NumVectorQ|_?NumSMatrixQ)]]:=NRationalEval[SlashM2D2[a]];


Sm4/:NRationalEval[Sm4[a:(_?NumSpinorQ|_?NumVectorQ|_?NumSMatrixQ)]]:=NRationalEval[SlashM[a]];


NRationalEval[ProjPlus]:=((DiagonalMatrix[{1,1,1,1}]+Ga5)/2);
NRationalEval[ProjMinus]:=((DiagonalMatrix[{1,1,1,1}]-Ga5)/2);


NRationalEval[Sm4[Gamma0]]:=Ga0;
NRationalEval[Sm4[Gamma1]]:=Ga1;
NRationalEval[Sm4[Gamma2]]:=Ga2;
NRationalEval[Sm4[Gamma3]]:=Ga3;
NRationalEval[Sm4[Gamma5]]:=Ga5;
NRationalEval[Sm4[Gamma0]]:=Ga0;
NRationalEval[Sm4[Gamma1]]:=Ga1;
NRationalEval[Sm4[Gamma2]]:=Ga2;
NRationalEval[Sm4[Gamma3]]:=Ga3;
NRationalEval[Sm4[Gamma5]]:=Ga5;

NRationalEval[Gamma0]:=Ga0;
NRationalEval[Gamma1]:=Ga1;
NRationalEval[Gamma2]:=Ga2;
NRationalEval[Gamma3]:=Ga3;
NRationalEval[Gamma5]:=Ga5;

NRationalEval[Sm2[Gamma0]]:=sigma0;
NRationalEval[Sm2[Gamma0]]:=sigma0;
NRationalEval[Sm2[Gamma1]]:=sigma1;
NRationalEval[Sm2[Gamma1]]:=sigma1;
NRationalEval[Sm2[Gamma2]]:=sigma2;
NRationalEval[Sm2[Gamma2]]:=sigma2;
NRationalEval[Sm2[Gamma3]]:=sigma3;
NRationalEval[Sm2[Gamma3]]:=sigma3;



NRationalEval[USpa[a_?NumSpinorQ]]:=1/Sqrt[2] List/@Flatten[{NumLa[a],NumLa[a]}];


NRationalEval[USpb[a_?NumSpinorQ]]:=1/Sqrt[2] List/@Flatten[{NumCLat[a],-NumCLat[a]}];


NRationalEval[UbarSpb[a_?NumSpinorQ]]:=1/Sqrt[2] {Flatten[{NumLat[a],-NumLat[a]}]};


NRationalEval[UbarSpa[a_?NumSpinorQ]]:=1/Sqrt[2] {Flatten[{NumCLa[a],NumCLa[a]}]};


NRationalEval[Spaa[a_?NumSpinorQ,b_?NumSpinorQ]]:=(NumCLa[a].NumLa[b])[[1,1]];
NRationalEval[Spbb[a_?NumSpinorQ,b_?NumSpinorQ]]:=(NumLat[a].NumCLat[b])[[1,1]];


(*NRationalEval[Spaa[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=*)
(*    (NRationalEval[UbarSpa[a]].(Dot@@(NRationalEval@*Explicit/@{dm})).NRationalEval[USpa[b]])[[1,1]];*)

(*NRationalEval[Spab[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=*)
(*    (NRationalEval[UbarSpa[a]].(Dot@@(NRationalEval@*Explicit/@{dm})).NRationalEval[USpb[b]])[[1,1]];*)

(*NRationalEval[Spba[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=*)
(*    (NRationalEval[UbarSpb[b]].(Dot@@ (NRationalEval@*Explicit/@{dm})).NRationalEval[USpa[b]])[[1,1]];*)

(*NRationalEval[Spbb[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=*)
(*    (NRationalEval[UbarSpb[a]].(Dot@@(NRationalEval@*Explicit/@{dm})).NRationalEval[USpb[b]])[[1,1]];*)


NRationalEval[Spaa[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumCLa[a].(Dot@@(NRationalEval/@(alternate[Sm2,CSm2][{dm}]/.SumM->Plus))).NumLa[b])//First//First;

NRationalEval[Spab[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumCLa[a].(Dot@@(NRationalEval/@(alternate[Sm2,CSm2][{dm}]/.SumM->Plus))).NumCLat[b])//First//First;

NRationalEval[Spba[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumLat[a].(Dot@@(NRationalEval/@(alternate[CSm2,Sm2][{dm}]/.SumM->Plus))).NumLa[b])//First//First;

NRationalEval[Spbb[a_?NumSpinorQ,dm__?NumSMatrixQ,b_?NumSpinorQ]]:=
    (NumLat[a].(Dot@@(NRationalEval/@(alternate[CSm2,Sm2][{dm}]/.SumM->Plus))).NumCLat[b])//First//First;


NRationalEval[s[a_?NumSpinorQ,b_?NumSpinorQ]]:=NRationalEval[2 Num4V[a].gg.Num4V[b]];


NRationalEval[s[a:(_?NumSpinorQ|_?NumVectorQ),b:(_?NumSpinorQ|_?NumVectorQ)]]:=NRationalEval[(Num4V[a]+Num4V[b]).gg.(Num4V[a]+Num4V[b])];


NRationalEval[s[a:((_?NumSpinorQ|_?NumVectorQ)..)]]:=NRationalEval[Plus@@((Num4V[#]&/@{a})).gg.Plus@@((Num4V[#]&/@{a}))];


NRationalEval[MP[a_?(NumVectorQ[#]||NumSpinorQ[#]&),b_?(NumVectorQ[#]||NumSpinorQ[#]&)]]:=NRationalEval[Num4V[a].gg.Num4V[b]];


NRationalEval[MP[a_?(NumVectorQ[#]||NumSpinorQ[#]&),b_?(NumVectorQ[#]||NumSpinorQ[#]&)]]:=NRationalEval[Num4V[a].gg.Num4V[b]];


NRationalEval[MP2[a_?(NumVectorQ[#]||NumSpinorQ[#]&)]]:=NRationalEval[Num4V[a].gg.Num4V[a]];


La/:NRationalEval[La[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=NRationalEval[SlashM2D2[m].NumCLat[a]];


Lat/:NRationalEval[Lat[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=(-1) NRationalEval[NumCLa[a].SlashM2D2[m]];


CLa/:NRationalEval[CLa[u:((m_?NumSMatrixQ)[a_?NumSpinorQ])]]:=Transpose[\[Epsilon]\[Epsilon].NRationalEval[La[u]]];


CLat/:NRationalEval[CLat[u:((m_?NumSMatrixQ)[a_?NumSpinorQ])]]:=\[Epsilon]\[Epsilon].Transpose[NRationalEval[Lat[u]]];


La/:NRationalEval[La[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=(NRationalEval[Dot@@alternate[SlashM2D,SlashM2D2][d].If[OddQ[Length[d]],NumCLat[a],NumLa[a]]]);


Lat/:NRationalEval[Lat[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=
    ((-1)^Length[d])(NRationalEval[If[OddQ[Length[d]],NumCLa[a],NumLat[a]].(Dot@@Reverse[alternate[SlashM2D,SlashM2D2][d]])]);



CLa/:NRationalEval[CLa[u:((d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ])]]:=Transpose[\[Epsilon]\[Epsilon].NRationalEval[La[u]]];


CLat/:NRationalEval[CLat[u:((d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ])]]:=\[Epsilon]\[Epsilon].Transpose[NRationalEval[Lat[u]]];


USpa/:NRationalEval[USpa[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=NRationalEval[Sm4[m]].USpb[a];


USpb/:NRationalEval[USpb[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=NRationalEval[Sm4[m]].USpa[a];


UbarSpa/:NRationalEval[UbarSpa[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=(-1)UbarSpb[a].NRationalEval[Sm4[m]];


UbarSpb/:NRationalEval[UbarSpb[(m_?NumSMatrixQ)[a_?NumSpinorQ]]]:=(-1)UbarSpa[a].NRationalEval[Sm4[m]];


USpa/:NRationalEval[USpa[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=(NRationalEval[(Sm4/@d).If[OddQ[Length[d]],USpb[a],USpa[a]]]);


USpb/:NRationalEval[USpb[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=(NRationalEval[(Sm4/@d).If[OddQ[Length[d]],USpa[a],USpb[a]]]);


UbarSpa/:NRationalEval[UbarSpa[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=(-1)^Length[d](NRationalEval[If[OddQ[Length[d]],UbarSpb[a],UbarSpa[a]].(Sm4/@Reverse[d])]);


UbarSpb/:NRationalEval[UbarSpb[(d:Dot[(_?NumSMatrixQ),(_?NumSMatrixQ)..])[a_?NumSpinorQ]]]:=(-1)^Length[d](NRationalEval[If[OddQ[Length[d]],UbarSpa[a],UbarSpb[a]].(Sm4/@Reverse[d])]);


NRationalEval[s:(Spaa[_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ],_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ]])]:=
    NRationalEval[UnCompact[s]];


NRationalEval[s:(Spab[_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ],_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ]])]:=
    NRationalEval[UnCompact[s]];


NRationalEval[s:(Spba[_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ],_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ]])]:=
    NRationalEval[UnCompact[s]];


NRationalEval[s:(Spbb[_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ],_?NumSpinorQ|(_?NumSMatrixQ)[_?NumSpinorQ]|(Dot[_?NumSMatrixQ,(_?NumSMatrixQ)..])[_?NumSpinorQ]])]:=
    NRationalEval[UnCompact[s]];

NRationalEval[La[a_?NumSpinorQ]]:=NRationalEval[NumLa[a]];
NRationalEval[Lat[a_?NumSpinorQ]]:=NRationalEval[NumLat[a]];
NRationalEval[CLa[a_?NumSpinorQ]]:=NRationalEval[NumCLa[a]];
NRationalEval[CLat[a_?NumSpinorQ]]:=NRationalEval[NumCLat[a]];

NRationalEval[lp[x_,y_]]:=2 NRational[MP[x,y]];
NRationalEval[lp[x:Repeated[_,{3,Infinity}]]]:= 2*NRationalEval/@MP@@@Subsets[{x},{2}]//Apply[Plus];

NRationalEval[Verbatim[Plus][expr___]] := Plus@@NRationalEval/@{expr};
NRationalEval[Verbatim[Times][expr___]] := Times@@NRationalEval/@{expr};
NRationalEval[Verbatim[List][expr___]] := List@@NRationalEval/@{expr};
NRationalEval[Verbatim[Power][expr___]] := Power@@NRationalEval/@{expr};

$numericHeads = {Spbb,Spaa,Spab,Spba,USpa,USpb,UbarSpa,UbarSpb,s,NumCLa,NumCLat,
  NumLa,NumLat,Num4V,CSm2,Sm2,Sm4,SmBA2,CSmBA2,SmBA4,La,Lat,CLa,CLat,MP2,MP,
  SlashM2D,lp};

$constantHeads= {Gamma0,Gamma1,Gamma2,Gamma3,Gamma5,ProjMinus,ProjPlus};

(*NRational[expr:Except[Alternatives@@Blank/@$numericHeads]] :=*)
(*Module[{uniqueObjects,uniqueObjectsEvaluated,evaluationRules},*)

(*  uniqueObjects = DeleteDuplicates@Cases[*)
(*      {expr},*)
(*      Alternatives@@Blank/@$numericHeads,*)
(*      Infinity*)
(*  ];*)
(*  uniqueObjectsEvaluated = NRationalEval/@uniqueObjects;*)
(*  evaluationRules = Dispatch@MapThread[Rule,{uniqueObjects,uniqueObjectsEvaluated}];*)

(*  expr /. evaluationRules*)
(*];*)

setupNRationalCache[]:=
Module[{},
  ClearAll[nRationalCache];

  nRationalCache[Verbatim[Plus][expr___]] :=
    nRationalCache[Plus[expr]] = Plus@@nRationalCache/@{expr};
  nRationalCache[Verbatim[Times][expr___]] :=
    nRationalCache[Times[expr]] = Times@@nRationalCache/@{expr};
  nRationalCache[Verbatim[List][expr___]] :=
    nRationalCache[List[expr]] = List@@nRationalCache/@{expr};
  nRationalCache[Verbatim[Power][expr___]] :=
    nRationalCache[Power[expr]] = Power@@nRationalCache/@{expr};

  nRationalCache[expr:Alternatives@@Join[(Verbatim[#][BlankSequence[]]&/@$numericHeads),Verbatim/@$constantHeads]]:=
    nRationalCache[expr] = NRationalEval[expr];
  nRationalCache[x_] := x;
];

NRational[expr_]:=Module[{},
  setupNRationalCache[];
  nRationalCache@expr
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


Protect/@$SpinorsFunctions;

End[];

End[];