(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: SAM *)
(* :Context: SAM` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-05-21 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`SAM`", {"SpinorHelicity6D`",
    "SpinorHelicity6D`Numerics`", "Utils`", "SpinorHelicity6D`Amplitudes`",
    "Spinors`","Spinors`NRational`","Spinors`ScalarKinematics`"
  }
];
(* Exported symbols added here with SymbolName::usage *)

SpinorHelicityToSAM::usage = "";
SAMToSpinorHelicity6D::usage = "";

GenerateSAMFourScalarMomentumKinematics::usage = "";
GenerateSAMTwoScalarMomentumKinematics::usage = "";

DefineATreeSAM::usage = "";

Begin["`Private`"];

(* ::Section:: *)
(* SpinorHelicity6D to SAM *)


(* Express S6Mom in terms of s and MP *)
(*S6MomToSAM[s6mom:S6Mom[momenta__],massless_List,scalars___List]:=*)
(*    Module[{masses,massReplacements,scalarPairs},*)

(*  massReplacements = Rule[Spinors`MP2[#2],Spinors`MP2[#1]]&@@@{scalars};*)
(*  scalarPairs = Select[{scalars},ContainsAll[{momenta},#]&];*)
(*  masses = -Sum[Spinors`MP2[momentum],{momentum,{momenta}}]*)
(*      + 2 Sum[Spinors`MP2[scalarPair[[1]]],{scalarPair,scalarPairs}];*)

(*  Spinors`s[momenta]+masses /. {*)
(*    Spinors`MP[i_, i_] /; MemberQ[massless, i] :> 0*)
(*  } /.{*)
(*    Spinors`s[i_,j_]-Spinors`MP[i_,i_]-Spinors`MP[j_,j_]:>2 Spinors`MP[i,j],*)
(*    Spinors`s[i_,j_]-Spinors`MP[i_,i_]/;MemberQ[massless,j] :>2 Spinors`MP[i,j]*)
(*  } /. massReplacements*)
(*];*)

(* Express S6Mom in terms of MP only *)
S6MomToSAM[s6mom:S6Mom[momenta__],massless_List,scalars___List]:=
    Module[{masses,massReplacements,scalarPairs,mps},

      massReplacements = Rule[Spinors`MP2[#2],Spinors`MP2[#1]]&@@@{scalars};
      scalarPairs = Select[{scalars},ContainsAll[{momenta},#]&];
      masses = + 2 Sum[Spinors`MP2[scalarPair[[1]]],{scalarPair,scalarPairs}];
      mps = Sum[2MP@@momentumPair, {momentumPair,Subsets[{momenta},{2}]}];

      mps+masses /. {
        Spinors`MP[i_, i_] /; MemberQ[massless, i] :> 0
      } /. massReplacements
    ];


spinorHelicity6DToSAMRules[massless_List,scalars___List] :=
    {
      UnderBar -> Identity,
      S6Mom[momenta__] :> S6MomToSAM[S6Mom[momenta],massless,scalars],
      S[momenta__] :> S6MomToSAM[S6Mom[momenta],massless,scalars],
      S4Mom[momenta__] :> S6MomToSAM[S6Mom[momenta],massless,scalars],
      PlusM[a___] :> SumM[a],
      Chain[$angle, a_, {b___}, c_, $angle] :> Spinors`Spaa[a, b, c],
      Chain[$angle, a_, {b___}, c_, $square] :> Spinors`Spab[a, b, c],
      Chain[$square, a_, {b___}, c_, $angle] :> Spinors`Spba[a, b, c],
      Chain[$square, a_, {b___}, c_, $square] :> Spinors`Spbb[a, b, c],
      SpinorAngleBracket -> Spinors`Spaa,
      SpinorSquareBracket -> Spinors`Spbb,
      Extramass[i_]^(p_.) * Extramasstilde[i_]^(p_.) :> Spinors`MP2[i]^p,
      Mom4DVector[i_][Null][$up] . Mom4DVector[j_][Null][$down] :> Spinors`MP[i,j],
      Mom4DVector[i_][Null][$down] . Mom4DVector[j_][Null][$up] :> Spinors`MP[i,j]
    };


SpinorHelicityToSAM[expr_,massless_List,scalars___List] := expr //.
    spinorHelicity6DToSAMRules[massless,scalars];

(* ::Section:: *)
(* SAM to SpinorHelicity6D *)


ClearAll[SAMToSpinorHelicity6D];
Options[SAMToSpinorHelicity6D] = {ScalarPairs->{},Massless->{}};

SAMToSpinorHelicity6D[opts:OptionsPattern[]] := Function[x,SAMToSpinorHelicity6D[x,opts]];

SAMToSpinorHelicity6D[expr_,opts:OptionsPattern[]]:=expr/.{
  lp[x___]->2*LorentzProduct4D[x],Spbb[x_,y_]:>SpinorSquareBracket[x,y],Spaa[x_,
    y_]:>SpinorAngleBracket[x,y],
  Spbb[a_,b__,c_]:>Chain[$square,a,UnderBar/@{b},c,$square],Spab[a_,b__,c_]:>Chain[$angle,a,UnderBar/@{b},c,$square],
  Spba[a_,b__,c_]:>Chain[$square,a,UnderBar/@{b},c,$angle],Spaa[a_,b__,c_]:>Chain[$angle,a,UnderBar/@{b},c,$angle],
  s[x___]:>SAMInvariantToSpinorHelicity6D[{x},PassRules[SAMToSpinorHelicity6D,SAMInvariantToSpinorHelicity6D,opts]]
}//.SAMToSpinorHelicity6DMomentumRules[PassRules[SAMToSpinorHelicity6D,SAMToSpinorHelicity6DMomentumRules,opts]];



Options[SAMInvariantToSpinorHelicity6D] = Options[SAMToSpinorHelicity6D];
SAMInvariantToSpinorHelicity6D[{x___},opts:OptionsPattern[]]:=
Module[{},
  Sum[Extramass[i]Extramasstilde[i],{i,{x}}]+2LorentzProduct4D[x]
];

Options[SAMToSpinorHelicity6DMomentumRules] = Options[SAMToSpinorHelicity6D];
SAMToSpinorHelicity6DMomentumRules[opts:OptionsPattern[]]:=
Module[{scalarPairs,scalars,massless},
  scalarPairs = OptionValue[ScalarPairs];
  scalars = scalarPairs//Flatten;
  massless = OptionValue[Massless];
  {
    LorentzProduct4D[x_,x_]/; MemberQ[#2,x] :> Extramass[x]Extramasstilde[x],
    S6Mom[x_] /; MemberQ[Join[#2,#3],x] :> 0,
    Splice@Table[Extramass[pair[[2]]]->-Extramass[pair[[1]]],{pair,#1}],
    Splice@Table[Extramasstilde[pair[[2]]]->-Extramasstilde[pair[[1]]],{pair,#1}],
    UnderBar[x_] /; MemberQ[#3,x] :> x,
    Extramass[x_]/; MemberQ[#3,x] :> 0,
    Extramasstilde[x_]/; MemberQ[#3,x] :> 0
  }&[scalarPairs,scalars,massless]
];


(* ::Section:: *)
(* Generate SAM Scalar Kinematics *)


Options[GenerateSAMFourScalarMomentumKinematics] =
    Join[Options[GenFourScalarMomenta], {Quiet->False, Parallel->False}];
GenerateSAMFourScalarMomentumKinematics[n_ /; n >= 4, nPoints_ /; nPoints >= 1, opts : OptionsPattern[]] :=
    Module[{labels, points, save$PMPMMetric = $PMPMMetric, k,table},
      labels = Table[Unique["p"], {i, n}];
      $PMPMMetric = False;

      table = If[
        OptionValue[Parallel]===True,
        ParallelTable,
        Table
      ];

      Block[{$Notebooks = False},
        k = 0;

        points = Function[tableHead,
          tableHead[
            GenFourScalarMomenta[labels[[{1, 2}]], labels[[{3, 4}]], labels[[5 ;;]],
              PassRules[GenerateSAMFourScalarMomentumKinematics, GenFourScalarMomenta, opts]];
            If[OptionValue[Quiet]===True, UpdateProgressBar[k,nPoints];];
            k++;
            (ClearKinematics[labels];#)&@<|
              "Scalars1" -> {Mom4DVectorN[labels[[1]]][Null][$up], Mom4DVectorN[labels[[2]]][Null][$up]},
              "Scalars2" -> {Mom4DVectorN[labels[[3]]][Null][$up], Mom4DVectorN[labels[[4]]][Null][$up]},
              "Massless" -> (Mom4DVectorN[#][Null][$up]& /@ labels[[5;;]]),
              "MasslessSpinors" -> ({SpinorUndotN[#][$lam][$up], SpinorDotN[#][$lam][$up]}& /@ labels[[5 ;;]])
            |>,
            {i, nPoints}
          ]
        ][table];
        If[OptionValue[Quiet]===True, StopProgressBar[nPoints];];

      ];

      $PMPMMetric = save$PMPMMetric;

      points
    ];

Options[GenerateSAMTwoScalarMomentumKinematics] = Options[GenTwoScalarMomenta];
GenerateSAMTwoScalarMomentumKinematics[n_ /; n >= 4, nPoints_ /; nPoints >= 1, opts : OptionsPattern[]] :=
    Module[{labels, points, save$PMPMMetric = $PMPMMetric},
      labels = Table[Unique["p"], {i, n}];
      $PMPMMetric = False;
      points = Table[
        GenTwoScalarMomenta[labels[[{1, 2}]], labels[[3 ;;]],
          PassRules[GenerateSAMTwoScalarMomentumKinematics, GenTwoScalarMomenta, opts]];
        (ClearKinematics[labels];#)&@<|
          "Scalars" -> {Mom4DVectorN[labels[[1]]][Null][$up], Mom4DVectorN[labels[[2]]][Null][$up]},
          "Massless" -> (Mom4DVectorN[#][Null][$up]& /@ labels[[3;;]]),
          "MasslessSpinors" -> ({SpinorUndotN[#][$lam][$up], SpinorDotN[#][$lam][$up]}& /@ labels[[3 ;;]])
        |>,
        {i, nPoints}
      ];

      $PMPMMetric = save$PMPMMetric;

      points
    ];

(*SetupSAMFourScalarMomentumKinematics[scalar1Labels : {_, _}, scalar2Labels : {_, _},*)
(*  masslessLabels_List, kinematicPoint_Association] /; Length@masslessLabels ===*)
(*  Length@kinematicPoint["Massless"] :=*)
(*    Module[{},*)

(*      MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars1", First@#2]]]&,*)
(*        scalar1Labels];*)
(*      MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars2", First@#2]]]&,*)
(*        scalar2Labels];*)

(*      MapIndexed[DeclareSpinorMomentum[#1, List /@ kinematicPoint[["MasslessSpinors", First@#2, 1]],*)
(*        {kinematicPoint[["MasslessSpinors", First@#2, 2]]}]&, masslessLabels];*)
(*    ];*)

(*SetupSAMTwoScalarMomentumKinematics[scalarLabels : {_, _}, masslessLabels_List,*)
(*  kinematicPoint_Association] /; Length@masslessLabels === Length@kinematicPoint["Massless"]:=*)
(*    Module[{},*)

(*      MapIndexed[DeclareLVectorMomentum[#1, kinematicPoint[["Scalars", First@#2]]]&, scalarLabels];*)

(*      MapIndexed[DeclareSpinorMomentum[#1, List /@ kinematicPoint[["MasslessSpinors", First@#2, 1]],*)
(*        {kinematicPoint[["MasslessSpinors", First@#2, 2]]}]&, masslessLabels];*)
(*    ];*)


(* ::Section:: *)
(* ATreeSAM ATreeSAMContact *)

DefineATreeSAM[] := (
  (* 4-pt *)
  ATreeSAM["S","+","+","S"][p1_,p2_,p3_,p4_]=SpinorHelicityToSAM[ATree["S","+","+","S"][p1,p2,p3,p4],{p2,p3},{p1,p4}];
  ATreeSAM["S","+","-","S"][p1_,p2_,p3_,p4_]=SpinorHelicityToSAM[ATree["S","+","-","S"][p1,p2,p3,p4],{p2,p3},{p1,p4}];
  ATreeSAM["S","+","S","+"][p1_,p2_,p3_,p4_]=SpinorHelicityToSAM[ATree["S","+","S","+"][p1,p2,p3,p4],{p2,p4},{p1,p3}];

  (* 5-pt *)
  ATreeSAM["S","+","+","+","S"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","+","+","S"][p1,p2,p3,p4,p5],{p2,p3,p4},{p1,p5}];
  ATreeSAM["S","+","+","-","S"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","+","-","S"][p1,p2,p3,p4,p5],{p2,p3,p4},{p1,p5}];
  ATreeSAM["S","+","-","+","S"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","-","+","S"][p1,p2,p3,p4,p5],{p2,p3,p4},{p1,p5}];
  ATreeSAM["S","+","+","S","+"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","+","S","+"][p1,p2,p3,p4,p5],{p2,p3,p5},{p1,p4}];

  (* 6-pt *)
  ATreeSAM["S","+","+","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","+","+","S"][p1,p2,p3,p4,p5,p6],{p2,p3,p4,p5},{p1,p6}];
  ATreeSAM["S","+","+","+","S","+"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","+","S","+"][p1,p2,p3,p4,p5,p6],{p2,p3,p4,p6},{p1,p5}];
  ATreeSAM["S","+","+","S","+","+"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","S","+","+"][p1,p2,p3,p4,p5,p6],{p2,p3,p5,p6},{p1,p4}];
  ATreeSAM["S","+","+","+","-","S"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","+","-","S"][p1,p2,p3,p4,p5,p6],{p2,p3,p4,p5},{p1,p6}];
  ATreeSAM["S","+","+","-","+","S"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","-","+","S"][p1,p2,p3,p4,p5,p6],{p2,p3,p4,p5},{p1,p6}];

  (* 7-pt *)
  ATreeSAM["S","+","+","+","+","+","S"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","+","+","+","S"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p4,p5,p6},{p1,p7}];
  ATreeSAM["S","+","+","+","+","S","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","+","+","S","+"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p4,p5,p7},{p1,p6}];
  ATreeSAM["S","+","+","+","S","+","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","+","S","+","+"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p4,p6,p7},{p1,p5}];


  (* 5-pt *)
  ATreeSAM["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","S","+","s","s"][p1,p2,p3,p4,p5],{p3},{p1,p2},{p4,p5}];
  ATreeSAM["S","+","S","s","s"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","S","s","s"][p1,p2,p3,p4,p5],{p2},{p1,p3},{p4,p5}];
  ATreeSAM["S","+","s","S","s"][p1_,p2_,p3_,p4_,p5_]=SpinorHelicityToSAM[ATree["S","+","s","S","s"][p1,p2,p3,p4,p5],{p2},{p1,p4},{p3,p5}];

  (* 6-pt *)
  ATreeSAM["S","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","S","+","+","s","s"][p1,p2,p3,p4,p5,p6],{p3,p4},{p1,p2},{p5,p6}];
  ATreeSAM["S","S","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","S","+","s","s","+"][p1,p2,p3,p4,p5,p6],{p3,p6},{p1,p2},{p4,p5}];

  ATreeSAM["S","+","+","S","s","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","S","s","s"][p1,p2,p3,p4,p5,p6],{p2,p3},{p1,p4},{p5,p6}];
  ATreeSAM["S","+","S","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","S","+","s","s"][p1,p2,p3,p4,p5,p6],{p2,p4},{p1,p3},{p5,p6}];
  ATreeSAM["S","+","S","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","S","s","+","s"][p1,p2,p3,p4,p5,p6],{p2,p5},{p1,p3},{p4,p6}];

  ATreeSAM["S","+","+","s","S","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","+","s","S","s"][p1,p2,p3,p4,p5,p6],{p2,p3},{p1,p5},{p4,p6}];
  ATreeSAM["S","+","s","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","s","+","S","s"][p1,p2,p3,p4,p5,p6],{p2,p4},{p1,p5},{p3,p6}];
  ATreeSAM["S","+","s","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_]=SpinorHelicityToSAM[ATree["S","+","s","S","+","s"][p1,p2,p3,p4,p5,p6],{p2,p5},{p1,p4},{p3,p6}];

  (* 7-pt *)
  ATreeSAM["S","S","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","S","+","+","+","s","s"][p1,p2,p3,p4,p5,p6,p7],{p3,p4,p5},{p1,p2},{p6,p7}];
  ATreeSAM["S","S","+","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","S","+","+","s","s","+"][p1,p2,p3,p4,p5,p6,p7],{p3,p4,p7},{p1,p2},{p5,p6}];

  ATreeSAM["S","+","+","+","S","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","+","S","s","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p4},{p1,p5},{p6,p7}];
  ATreeSAM["S","+","+","S","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","S","+","s","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p5},{p1,p4},{p6,p7}];
  ATreeSAM["S","+","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","S","+","+","s","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p4,p5},{p1,p3},{p6,p7}];
  ATreeSAM["S","+","+","S","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","S","s","+","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p6},{p1,p4},{p5,p7}];
  ATreeSAM["S","+","S","+","s","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","S","+","s","+","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p4,p6},{p1,p3},{p5,p7}];
  ATreeSAM["S","+","S","+","s","s","+"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","S","+","s","s","+"][p1,p2,p3,p4,p5,p6,p7],{p2,p4,p7},{p1,p3},{p5,p6}];

  ATreeSAM["S","+","+","+","s","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","+","s","S","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p4},{p1,p6},{p5,p7}];
  ATreeSAM["S","+","+","s","+","S","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","s","+","S","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p5},{p1,p6},{p4,p7}];
  ATreeSAM["S","+","+","s","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","+","s","S","+","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p3,p6},{p1,p5},{p4,p7}];
  ATreeSAM["S","+","s","+","S","+","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_]=SpinorHelicityToSAM[ATree["S","+","s","+","S","+","s"][p1,p2,p3,p4,p5,p6,p7],{p2,p4,p6},{p1,p5},{p3,p7}];

  (* Contact-term contributions only *)
  ATreeSAMContact["S","S","s","s"][p1_,p2_,p3_,p4_] = SpinorHelicityToSAM[ATreeContact["S","S","s","s"][p1,p2,p3,p4],{},{p1,p2},{p3,p4}];
  ATreeSAMContact["S","S","+","s","s"][p1_,p2_,p3_,p4_,p5_] = SpinorHelicityToSAM[ATreeContact["S","S","+","s","s"][p1,p2,p3,p4,p5],{p3},{p1,p2},{p4,p5}];
  ATreeSAMContact["S","S","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_] = SpinorHelicityToSAM[ATreeContact["S","S","+","+","s","s"][p1,p2,p3,p4,p5,p6],{p3,p4},{p1,p2},{p5,p6}];
  ATreeSAMContact["S","S","+","+","+","s","s"][p1_,p2_,p3_,p4_,p5_,p6_,p7_] = SpinorHelicityToSAM[ATreeContact["S","S","+","+","+","s","s"][p1,p2,p3,p4,p5,p6,p7],{p3,p4,p5},{p1,p2},{p6,p7}];
);

End[]; (* `Private` *)

EndPackage[]