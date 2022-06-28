(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: Unitarity *)
(* :Context: Unitarity` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`Unitarity`"];

Unitarity::usage = "Unitarity is a package collecting 4-dimensional unitarity technology based on the SpinorHelicity6D package.";

$star::usage = "";
$unstar::usage = "";
$plus::usage = "";
$minus::usage = "";

LargeMuLimit::usage = "Option for FixBoxLoopMomentum6D and GetBoxRationalTerm. If True, computation will use expansion
	of the sqrt in the box loop momentum. In this case, the parameters MuVar and MuTildeVar will be interpreted as the
	square roots of the respective parameters, as to avoid unnecessary complexity. This is also taken into account by
	GetBoxRationalTerm.";

LargeTLimit::usage = "";


MuVar::usage = "";
MuTildeVar::usage = "";
tVar::usage = "";
yVar::usage = "";

$LoopMomentumMu::usage = "";
$LoopMomentumMuTilde::usage = "";
$LoopMomentumT::usage = "";
$LoopMomentumY::usage = "";

MakeBoxes[$LoopMomentumMu,StandardForm|TraditionalForm] := InterpretationBox["\[Mu]",$LoopMomentumMu];
MakeBoxes[$LoopMomentumMuTilde,StandardForm|TraditionalForm] := InterpretationBox[OverscriptBox["\[Mu]","~"],$LoopMomentumMuTilde];
MakeBoxes[$LoopMomentumT,StandardForm|TraditionalForm] := InterpretationBox["t",$LoopMomentumT];
MakeBoxes[$LoopMomentumY,StandardForm|TraditionalForm] := InterpretationBox["y",$LoopMomentumY];

BubbleReference::usage = "";
$bubbleReference::usage = "";

Gravity::usage = "";

$replacedSquareRoots = "";

SqrtUnique::usage = "Unique square roots that were replaced in box and bubble-triangle loop momenta.";
SqrtUniqueArgument::usage = "Returns argument of a replaced square root.";

$cutBubble::usage = "";

Begin["`Private`"];
End[]; (* `Private` *)

EndPackage[];
