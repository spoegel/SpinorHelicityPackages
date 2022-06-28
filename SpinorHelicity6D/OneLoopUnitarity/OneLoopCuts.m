(* Wolfram Language Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: OneLoopCuts *)
(* :Context: OneLoopCuts` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-17 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

(* For new style packages see: https://mathematica.stackexchange.com/a/176489) *)
(* Declare package context *)
Package["SpinorHelicity6D`OneLoopUnitarity`"]

(* Import other packages *)
PackageImport["SpinorHelicity6D`"]
PackageImport["SpinorHelicity6D`Unitarity`"]
PackageImport["Utils`"]

(* Exporting *)
PackageExport["GenOneLoopConfigs"]


PackageScope["oneLoopBoxConfigQ"]
oneLoopBoxPattern = {_, _, _, _};
oneLoopBoxConfigQ[x_] := MatchQ[x, oneLoopBoxPattern];
PackageScope["oneLoopTriangleConfigQ"]
oneLoopTrianglePattern = {a_, b_, c_}  /; FreeQ[{a, b, c}, $cutBubble];
oneLoopTriangleConfigQ[x_] := MatchQ[x, oneLoopTrianglePattern];
PackageScope["oneLoopBubbleConfigQ"]
oneLoopBubblePattern = {_, _};
oneLoopBubbleConfigQ[x_] := MatchQ[x, oneLoopBubblePattern];
PackageScope["oneLoopBubbleTriangleConfigQ"]
oneLoopBubbleTrianglePattern = {a_, b_, c_} /; !FreeQ[{a, b, c}, $cutBubble];
oneLoopBubbleTriangleConfigQ[x_] := MatchQ[x, oneLoopBubbleTrianglePattern];

PackageScope["oneLoopConfigQ"]
oneLoopConfigPattern = oneLoopBoxPattern | oneLoopTrianglePattern | oneLoopBubblePattern | oneLoopBubbleTrianglePattern;
oneLoopConfigQ[x_] := MatchQ[x, oneLoopConfigPattern];

PackageScope["OneLoopConfigToLoopMomentumConfig"]
OneLoopConfigToLoopMomentumConfig[config_] := Join[#["Outer"][[All, 2]], #["Inner"][[All, 2]]]& /@ config;
OneLoopConfigToLoopMomentumConfig[config : oneLoopBubbleTrianglePattern] /; Position[config, $cutBubble][[All, 1]] === {1, 3} := OneLoopConfigToLoopMomentumConfig@RotateLeft@config;

PackageScope["FuseVertices"]
FuseVertices[v1_Association, v2_Association] := <|"Outer" -> Join[v1["Outer"], v2["Outer"]], "Inner" -> Join[v2["Inner"], v1["Inner"]], "Properties" -> Union[v1["Properties"], v2["Properties"]]|>;
PackageScope["RemoveCutBubble"]
RemoveCutBubble[vertex_Association] := ReplacePart[vertex, {"Properties" -> vertex["Properties"] /. $cutBubble -> Nothing}];

PackageScope["FuseOneLoopCutBubble"]
FuseOneLoopCutBubble[config : oneLoopBubbleTrianglePattern] /; Position[config, $cutBubble][[All, 1]] === {1, 3} :=
    ReplacePart[config, {3 -> RemoveCutBubble@FuseVertices[config[[3]], config[[1]]], 1 -> Nothing}];
FuseOneLoopCutBubble[config : oneLoopBubbleTrianglePattern] /; Position[config, $cutBubble][[All, 1]] === {1, 2} :=
    ReplacePart[config, {1 -> RemoveCutBubble@FuseVertices[config[[1]], config[[2]]], 2 -> Nothing}];
FuseOneLoopCutBubble[config : oneLoopBubbleTrianglePattern] /; Position[config, $cutBubble][[All, 1]] === {2, 3} :=
    ReplacePart[config, {2 -> RemoveCutBubble@FuseVertices[config[[2]], config[[3]]], 3 -> Nothing}];

FuseOneLoopCutBubble[config_] := config




Options[GetOneLoopCutPositions] = {DoNotCut -> {}};
GetOneLoopCutPositions[momenta_List, nCuts_Integer, opts : OptionsPattern[]] :=
    DeleteCases[
      Subsets[Range@Length@momenta, {nCuts}],
      {___, x : Alternatives @@ OptionValue[DoNotCut], ___} | {a_, b_} /; (b - a ===1 ||
          b -a === (Length@momenta-1))
    ];


Options[CutPositionsToCut] = {DoNotCut -> {}};
CutPositionsToCut[momenta_List, cutPositions_List] :=
    SplitAt[RotateRight[momenta, Length[momenta] - cutPositions[[-1]]], cutPositions + Length[momenta] - cutPositions[[-1]]]

DecorateOneLoopConfigs[configBare_, doNotCut_] := <|"Outer" -> #, "Inner" -> {}, "Properties" -> {"DoNotCut" -> doNotCut}|>& /@ configBare;
RemoveDoNotCut[configs_] := MapAt[# /. HoldPattern["DoNotCut" -> _] -> Nothing &, configs, {All, All, "Properties"}];



Options[GenOneLoopConfigs] = {DoNotCut -> {}};
GenOneLoopConfigs[momenta_List, opts : OptionsPattern[]] :=
    CutOneLoopBubbles[DecorateOneLoopConfigs[#, momenta[[OptionValue[DoNotCut]]]]& /@
        Join[
          GenOneLoopConfigs[momenta, 4, PassRules[GenOneLoopConfigs, GenOneLoopConfigs, opts]],
          GenOneLoopConfigs[momenta, 3, PassRules[GenOneLoopConfigs, GenOneLoopConfigs, opts]],
          GenOneLoopConfigs[momenta, 2, PassRules[GenOneLoopConfigs, GenOneLoopConfigs, opts]]
        ]
    ];

GenOneLoopConfigs[momenta_List, nCuts_Integer, opts : OptionsPattern[]] :=
    CutPositionsToCut[momenta, #]& /@ GetOneLoopCutPositions[momenta, nCuts, PassRules[GenOneLoopConfigs, GetOneLoopCutPositions, opts]];


CutOneLoopBubbles[configs_] := CutOneLoopBubble /@ configs;

CutOneLoopBubble[config_List] /; Length@config != 2 := config;
CutOneLoopBubble[config_List] := Splice@Replace[Join[
  {config},
  Thread[{CutOneLoopBubbleVertex@config[[1]], config[[2]]}],
  Thread[{config[[1]], CutOneLoopBubbleVertex@config[[2]]}]
], {{{a__}, b_} :> {a, b}, {a_, {b__}} :> {a, b}}, 1];

(* Can/should be extended for inner cuts as well *)
CutOneLoopBubbleVertex[vertex : <|"Outer" -> pOut_, "Inner" -> pIn_, "Properties" -> props_|>] := Module[
  {doNotCut},
  If[MemberQ[props, HoldPattern["DoNotCut" -> _]],
    doNotCut = Join @@ Cases[props, HoldPattern["DoNotCut" -> x_List]][[All, 2]];,
    doNotCut = {};
  ];
  RemoveDoNotCuts[doNotCut] /@ CutOuterVertex[vertex]
];


CutOuterVertex[vertex_Association] := Module[{cuts},
  cuts = Replace[
    Transpose /@ Tuples[{Table[{#1[[i + 1 ;;]], #1[[;; i]]}, {i, 0, Length@#1}], Table[{#2[[i + 1 ;;]], #2[[;; i]]}, {i, 0, Length@#2}]}],
    {{___, {{}, _}, {_, {}}, ___} -> Nothing, {___, {_, {}}, {{}, _}, ___} -> Nothing}, 2
  ]&[vertex["Outer"], vertex["Inner"]];

  {<|"Outer" -> #[[2, 1]], "Inner" -> #[[1, 2]], "Properties" -> Append[vertex["Properties"], $cutBubble]|>,
    <|"Outer" -> #[[1, 1]], "Inner" -> #[[2, 2]], "Properties" -> Append[vertex["Properties"], $cutBubble]|>}& /@ cuts
];


RemoveDoNotCuts[doNotCut_] := Function[cut, RemoveDoNotCuts[cut, doNotCut]];
RemoveDoNotCuts[cut_, doNotCut_] /; MemberQ[doNotCut, cut[[1, "Outer", -1]]] := Nothing;
RemoveDoNotCuts[cut_, doNotCut_] := cut;
