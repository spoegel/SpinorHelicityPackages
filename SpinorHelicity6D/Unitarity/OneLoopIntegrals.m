(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: OneLoopIntegrals *)
(* :Context: OneLoopIntegrals` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`Unitarity`OneLoopIntegrals`",
  "SpinorHelicity6D`","SpinorHelicity6D`Unitarity`","Utils`"];
(* Exported symbols added here with SymbolName::usage *)


IntegralBoxMu4::usage = "";
IntegralBoxMu6::usage = "";
IntegralBoxMu8::usage = "";


IntegralTriangleMu2::usage = "";
IntegralTriangleMu4::usage = "";
IntegralTriangleMu6::usage = "";


IntegralBubbleMu2::usage = "";
IntegralBubbleMu4::usage = "";
IntegralBubbleMu6::usage = "";

Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Mu2 Integrals*)


(* ::Subsubsection::Closed:: *)
(*Boxes*)


IntegralBoxMu8[k1_List,k2_List,k3_List,k4_List]:=-(*Global`BoxIntegralCoeff[8]**)1/840*(2*Sum[(S6Mom@@corner)^2,{corner,{k1,k2,k3,k4}}]+
    2*Sum[S6Mom@@cornerPair[[1]]*S6Mom@@cornerPair[[2]],{cornerPair,Partition[{k1,k2,k3,k4},2,1,1]}]+
    2 (S6Mom@@Join[k1,k2]^2+S6Mom@@Join[k2,k3]^2)+
    2 (S6Mom@@Join[k1,k2]+S6Mom@@Join[k2,k3])*Sum[(S6Mom@@corner),{corner,{k1,k2,k3,k4}}]+
    S6Mom@@k1*S6Mom@@k3+S6Mom@@k2*S6Mom@@k4+
    S6Mom@@Join[k1,k2]*S6Mom@@Join[k2,k3]);


IntegralBoxMu6[k1_List,k2_List,k3_List,k4_List]:=-(*Global`BoxIntegralCoeff[6]*)1/60*(Sum[(S6Mom@@corner),{corner,{k1,k2,k3,k4}}]+S6Mom@@Join[k1,k2]+S6Mom@@Join[k2,k3]);


IntegralBoxMu4[k1_List,k2_List,k3_List,k4_List]:=-(*Global`BoxIntegralCoeff[4]*)1/6;


(* ::Subsubsection::Closed:: *)
(*Triangles*)


IntegralTriangleMu6[k1_List,k2_List,k3_List]:=-(*Global`TriangleIntegralCoeff[6]**)1/180*(Sum[(S6Mom@@corner)^2,{corner,{k1,k2}}]+
    S6Mom@@Join[k1,k2]*Sum[(S6Mom@@corner),{corner,{k1,k2}}]+
    S6Mom@@k1*S6Mom@@k2+
    S6Mom@@Join[k1,k2]^2);


IntegralTriangleMu4[k1_List,k2_List,k3_List]:=-(*Global`TriangleIntegralCoeff[4]*)1/24*(Sum[(S6Mom@@corner),{corner,{k1,k2}}]+S6Mom@@Join[k1,k2]);


IntegralTriangleMu2[k1_List,k2_List,k3_List]:=-(*Global`TriangleIntegralCoeff[2]*)1/2;


(* ::Subsubsection::Closed:: *)
(*Bubbles*)


IntegralBubbleMu6[k1_List,k2_List]:=-(*Global`BubbleIntegralCoeff[6]**)1/420*S6Mom@@k1^3;


IntegralBubbleMu4[k1_List,k2_List]:=-(*Global`BubbleIntegralCoeff[4]*)1/60*S6Mom@@k1^2;


IntegralBubbleMu2[k1_List,k2_List]:=-(*Global`BubbleIntegralCoeff[2]**)1/6*S6Mom@@k1;



End[]; (* `Private` *)

EndPackage[]