(* Wolfram Language Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: OneLoopUnitarity *)
(* :Context: OneLoopUnitarity` *)
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
PackageImport["SpinorHelicity6D`Unitarity`LoopMomentum`"]
PackageImport["SpinorHelicity6D`Unitarity`RationalTerms`"]
PackageImport["SpinorHelicity6D`Unitarity`OneLoopIntegrals`"]
PackageImport["SpinorHelicity6D`Unitarity`ConfigurationTools`"]
PackageImport["Utils`"]
PackageImport["RationalSeries`"]
PackageImport["NumSymb`"]

(* Export functions *)
PackageExport["ComputeOneLoopConfig"]


(* ::Section::Closed:: *)
(*Loop-Momenta*)


(* ::Subsection::Closed:: *)
(*Prototypes*)


Options[CreateLoopMomentumAssoc] = {LoopMomentumPostProcessingFirstPass -> Identity, Gravity -> False, LargeTLimit -> True, LargeMuLimit -> True, ComputeLoopMomentumSpinors -> False, BubbleReference->$bubbleReference};


CreateLoopMomentumAssoc[config_?oneLoopBoxConfigQ,opts:OptionsPattern[]]:=
    CreateBoxLoopMomentumAssoc[
      OneLoopConfigToLoopMomentumConfig[config],
      PassRules[CreateLoopMomentumAssoc,CreateBoxLoopMomentumAssoc,opts]
    ];
CreateLoopMomentumAssoc[config_?oneLoopTriangleConfigQ,opts:OptionsPattern[]]:=
    CreateTriangleLoopMomentumAssoc[
      OneLoopConfigToLoopMomentumConfig[config],
      PassRules[CreateLoopMomentumAssoc,CreateTriangleLoopMomentumAssoc,opts]
    ];
CreateLoopMomentumAssoc[config_?oneLoopBubbleConfigQ,opts:OptionsPattern[]]:=
    CreateBubbleLoopMomentumAssoc[
      OneLoopConfigToLoopMomentumConfig[config],
      PassRules[CreateLoopMomentumAssoc,CreateBubbleLoopMomentumAssoc,opts]
    ];
CreateLoopMomentumAssoc[config_?oneLoopBubbleTriangleConfigQ,opts:OptionsPattern[]]:=
    CreateBubbleTriangleLoopMomentumAssoc[
      OneLoopConfigToLoopMomentumConfig@FuseOneLoopCutBubble@config,
      OneLoopConfigToLoopMomentumConfig@config,
      PassRules[CreateLoopMomentumAssoc,CreateBubbleTriangleLoopMomentumAssoc,opts]
    ];


ClearLoopMomentumDefinitions[loopMomAssoc_Association]:=Join@@(Through[loopMomAssoc["Momenta"]@#]&/@loopMomAssoc["Variations"])//ClearKinematics


(* ::Subsection::Closed:: *)
(*Implementation*)


(* ::Subsubsection::Closed:: *)
(*Box*)


Options[CreateBoxLoopMomentumAssoc] = Join[Options[CreateLoopMomentumAssoc], {LoopName -> ""}];

CreateBoxLoopMomentumAssoc[config_?oneLoopBoxConfigQ,opts:OptionsPattern[]] := Module[
  {assocOut=<||>,l1=Unique["l" <> OptionValue[LoopName] <> "1$"],l2=Unique["l" <> OptionValue[LoopName] <> "2$"],l3=Unique["l3" <> OptionValue[LoopName] <> "$"],l4=Unique["l" <> OptionValue[LoopName] <> "4$"],
    Mu=Unique["Mu" <> OptionValue[LoopName] <> "$"],MuTilde=Unique["MuTilde" <> OptionValue[LoopName] <> "$"]},

  FixBoxLoopMomentum6D[config,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l1,l2,l3,l4}, LargeMuLimit->True,
    PassRules[CreateBoxLoopMomentumAssoc,FixBoxLoopMomentum6D,opts]];

  assocOut["Type"] = "Box";
  assocOut["Momenta"] = {l1,l2,l3,l4};
  assocOut["Variations"] = {$plus,$minus};
  assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde|>;

  assocOut
];


(* ::Subsubsection::Closed:: *)
(*Triangle*)


Options[CreateTriangleLoopMomentumAssoc] = Join[Options[CreateLoopMomentumAssoc], {LoopName -> ""}];

CreateTriangleLoopMomentumAssoc[config_,opts:OptionsPattern[]] /; Length@config == 3 := Module[
  {assocOut=<||>,l1=Unique["l" <> OptionValue[LoopName] <> "1$"],l2=Unique["l" <> OptionValue[LoopName] <> "2$"],l3=Unique["l3" <> OptionValue[LoopName] <> "$"],
    Mu=Unique["Mu" <> OptionValue[LoopName] <> "$"],MuTilde=Unique["MuTilde" <> OptionValue[LoopName] <> "$"],t=Unique["t" <> OptionValue[LoopName] <> "$"]},

  FixTriangleLoopMomentum6D[config,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l1,l2,l3},
    PassRules[CreateTriangleLoopMomentumAssoc,FixTriangleLoopMomentum6D,opts]];

  assocOut["Type"] = "Triangle";
  assocOut["Momenta"] = {l1,l2,l3};
  assocOut["Variations"] = {$star,$unstar};
  assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;

  assocOut
];


(* ::Subsubsection::Closed:: *)
(*Bubble*)


Options[CreateBubbleLoopMomentumAssoc] = Join[Options[CreateLoopMomentumAssoc], {LoopName -> ""}];

CreateBubbleLoopMomentumAssoc[config_,opts:OptionsPattern[]] /; Length@config == 2 := Module[
  {assocOut=<||>,yIntegrals,l1=Unique["l" <> OptionValue[LoopName] <> "1$"],l2=Unique["l" <> OptionValue[LoopName] <> "2$"],
    Mu=Unique["Mu" <> OptionValue[LoopName] <> "$"],MuTilde=Unique["MuTilde" <> OptionValue[LoopName] <> "$"],t=Unique["t" <> OptionValue[LoopName] <> "$"],y=Unique["y" <> OptionValue[LoopName] <> "$"]},

  yIntegrals = FixBubbleLoopMomentum6D[config,OptionValue[BubbleReference],yVar->y,tVar->t,MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l1,l2},
    PassRules[CreateBubbleLoopMomentumAssoc,FixBubbleLoopMomentum6D,opts]];

  assocOut["Type"] = "Bubble";
  assocOut["Momenta"] = {l1,l2};
  assocOut["Variations"] = {Null};
  assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t,"y"->y|>;
  assocOut["yIntegrals"] = ToNum@yIntegrals;

  assocOut
];


(* ::Subsubsection::Closed:: *)
(*BubbleTriangle*)

Options[CreateBubbleTriangleLoopMomentumAssoc] = Join[Options[CreateLoopMomentumAssoc], {LoopName -> ""}];

CreateBubbleTriangleLoopMomentumAssoc[baseConfig_,subConfig_,opts:OptionsPattern[]] /; Length@subConfig == 3 := Module[
  {assocOut=<||>,tIntegrals,splittingMomenta,l1=Unique["l" <> OptionValue[LoopName] <> "1$"],l2=Unique["l" <> OptionValue[LoopName] <> "2$"],l3=Unique["l" <> OptionValue[LoopName] <> "3$"],
    Mu=Unique["Mu" <> OptionValue[LoopName] <> "$"],MuTilde=Unique["MuTilde" <> OptionValue[LoopName] <> "$"],t=Unique["t" <> OptionValue[LoopName] <> "$"],y=Unique["y" <> OptionValue[LoopName] <> "$"]},

  splittingMomenta = GetSplittingMomenta[baseConfig,subConfig];

  tIntegrals = FixBubbleTriSubLoopMomentum6D[baseConfig,splittingMomenta,OptionValue[BubbleReference],tVar->t,
    MuVar->Mu,MuTildeVar->MuTilde,LoopMomentumLabels->{l1,l2,l3},
    PassRules[CreateBubbleTriangleLoopMomentumAssoc,FixBubbleTriSubLoopMomentum6D,opts]];

  assocOut["Type"] = "BubbleTriangle";
  assocOut["Momenta"] = {l1,l2,l3};
  assocOut["Variations"] = {$plus,$minus};
  assocOut["Parameters"] = <|"Mu"->Mu,"MuTilde"->MuTilde,"t"->t|>;
  assocOut["tIntegrals"] = ToNum@tIntegrals;

  assocOut
];


(* ::Section::Closed:: *)
(*Rational Term*)


ComputeOneLoopRationalTerm[config_?oneLoopBoxConfigQ,trees_,loopMom_,opts:OptionsPattern[]] :=
    GetBoxRationalTerm[ToNum @ Activate @ trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
      LargeMuLimit->True, PassRules[ComputeOneLoopRationalTerm,GetBoxRationalTerm,opts]];


ComputeOneLoopRationalTerm[config_?oneLoopTriangleConfigQ,trees_,loopMom_,opts:OptionsPattern[]] :=
    GetTriangleRationalTerm[ToNum @ Activate @ trees,MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
      tVar->loopMom["Parameters"]["t"], PassRules[ComputeOneLoopRationalTerm,GetTriangleRationalTerm,opts]];


ComputeOneLoopRationalTerm[config_?oneLoopBubbleConfigQ,trees_,loopMom_,opts:OptionsPattern[]] :=
    GetBubbleBubRationalTerm[ToNum @ Activate @ trees,loopMom["yIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
      tVar->loopMom["Parameters"]["t"],yVar->loopMom["Parameters"]["y"],PassRules[ComputeOneLoopRationalTerm,GetBubbleBubRationalTerm,opts]];


ComputeOneLoopRationalTerm[config_?oneLoopBubbleTriangleConfigQ,trees_,loopMom_,opts:OptionsPattern[]] :=
    GetBubbleTriRationalTerm[ToNum @ Activate @ trees,loopMom["tIntegrals"],MuVar->loopMom["Parameters"]["Mu"],MuTildeVar->loopMom["Parameters"]["MuTilde"],
      tVar->loopMom["Parameters"]["t"],PassRules[ComputeOneLoopRationalTerm,GetBubbleTriRationalTerm,opts]];


(* ::Section::Closed:: *)
(*Integrals*)


(* ::Subsection::Closed:: *)
(*SecondLoopConfigIntegral*)


Options[OneLoopConfigIntegrals]={Gravity->False};


OneLoopConfigIntegrals[config_?oneLoopBoxConfigQ,opts:OptionsPattern[]]:=
    If[OptionValue[Gravity]===True,
      Through[{Apply[IntegralBoxMu4],Apply[IntegralBoxMu6],Apply[IntegralBoxMu8]}[OneLoopConfigToLoopMomentumConfig@config]],
      Through[{Apply[IntegralBoxMu4]}[OneLoopConfigToLoopMomentumConfig@config]]
    ];


OneLoopConfigIntegrals[config_?oneLoopTriangleConfigQ,opts:OptionsPattern[]] :=
    If[OptionValue[Gravity]===True,
      Through[{Apply[IntegralTriangleMu2],Apply[IntegralTriangleMu4],Apply[IntegralTriangleMu6]}[OneLoopConfigToLoopMomentumConfig@config]],
      Through[{Apply[IntegralTriangleMu2]}[OneLoopConfigToLoopMomentumConfig@config]]
    ];


OneLoopConfigIntegrals[config_?oneLoopBubbleConfigQ,opts:OptionsPattern[]] :=
    If[OptionValue[Gravity]===True,
      Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[OneLoopConfigToLoopMomentumConfig@config]],
      Through[{Apply[IntegralBubbleMu2]}[OneLoopConfigToLoopMomentumConfig@config]]
    ];


OneLoopConfigIntegrals[config_?oneLoopBubbleTriangleConfigQ,opts:OptionsPattern[]] :=
    If[OptionValue[Gravity]===True,
      Through[{Apply[IntegralBubbleMu2],Apply[IntegralBubbleMu4]}[OneLoopConfigToLoopMomentumConfig@FuseOneLoopCutBubble@config]],
      Through[{Apply[IntegralBubbleMu2]}[OneLoopConfigToLoopMomentumConfig@FuseOneLoopCutBubble@config]]
    ];


(* ::Section::Closed:: *)
(*Putting Everything Together: ComputeOneLoopConfig*)

Options[ComputeOneLoopConfig] = {};
ComputeOneLoopConfig[config_?oneLoopConfigQ,opts:OptionsPattern[]]:= Module[{loopMomentum,trees,rational},
  loopMomentum = CreateLoopMomentumAssoc[config,PassRules[ComputeOneLoopConfig,CreateLoopMomentumAssoc,opts]];

  trees = ConstructLoopAmplitudes[config,#,ClosedLoop->True,PassRules[ComputeOneLoopConfig,ConstructLoopAmplitudes,opts]]&/@
      (Through[loopMomentum["Momenta"][#]]&/@loopMomentum["Variations"]);

  rational = Sum[ComputeOneLoopRationalTerm[config,treeVariation,loopMomentum,PassRules[ComputeOneLoopConfig,ComputeOneLoopRationalTerm,opts]],{treeVariation,trees}];

  ClearLoopMomentumDefinitions[loopMomentum];

  2*(rational . OneLoopConfigIntegrals[config,PassRules[ComputeOneLoopConfig,OneLoopConfigIntegrals,opts]])//ToNum
]
