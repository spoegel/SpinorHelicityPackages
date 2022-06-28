(* ::Package:: *)

(* :Title: SpinorFunctions.m -- a package collecting 6D *)

(* :Context: SpinorHelicity6D` *)

(* :Summary:

 *)

(* :Copyright: \[Copyright] 2020 by Sebastian Poegel *)

(* :Package Version: 0.1 *)

(* :Mathematica Version: 12.0 *)

(* :History:
   0.1 Copying functions from notebook
*)

(* :Keywords: *)

(* :Sources:
*)

(* :Warnings:

*)

(* :Limitations:

*)

(* :Discussion:
   <description of algorithm, information for experts>
*)

(* :Requirements:
   SpinorHelicity6D.wl handling 6D Spinors, written by M. Huber and S. de Angelis
*)

(* :Todo:
   -Rename AngAngInvChain into Chain6DN, and define Chain6D, which looks nice in a notebook and stays unevaluated. Then add compat layer.
*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`SpinorFunctions`", "SpinorHelicity6D`","Utils`"];

Unprotect["SpinorHelicity6D`SpinorFunctions`*"];
ClearAll["SpinorHelicity6D`SpinorFunctions`*"];
ClearAll["SpinorHelicity6D`SpinorFunctions`Private`*"];

(* usage messages for the exported functions and the context itself *)

SpinorFunctions::usage = "SpinorFunctions.wl is a package collecting functions of 4D and 6D spinors.";

PTF::usage = "Generates Parke-Taylor Factor. Variants: PTF[_List], or PTF[_Integer].";
PTFn::usage = "(DEPRECATED, use PTF) Generates Parke-Taylor Factor. Variants: PTFn[_List], or PTFn[_Integer].";
PTFOpen::usage = "";
AntiPTF::usage = "Generates anti Parke-Taylor Factor. Variants: AntiPTF[_List], or AntiPTF[_Integer].";

AngAngInvChain::usage = "Evaluates an Angle-Angle chain of 6D spinors in terms of 6D spinor products.";
AngSquInvChain::usage = "Evaluates an Angle-Square chain of 6D spinors in terms of 6D spinor products.";
SquAngInvChain::usage = "Evaluates an Square-Angle chain of 6D spinors in terms of 6D spinor products.";


TrM::usage = "Unevaluated trace minus of four 4D momenta.";
TrP::usage = "Unevaluated trace plus of four 4D momenta.";
Tr5::usage = "Unevaluated trace 5 of four 4D momenta.";
TrMN::usage = "Evaluated trace minus of four 4D momenta. Evaluation via Mom4DN[Null], i.e. momenta do not need to be massless.";
TrPN::usage = "Evaluated trace plus of four 4D momenta.  Evaluation via Mom4DN[Null], i.e. momenta do not need to be massless.";
Tr5N::usage = "Evaluated trace 5 of four 4D momenta.  Evaluation via Mom4DN[Null], i.e. momenta do not need to be massless.";

TracesToSpinors::usage = "Expands traces in terms of spinor products. Must be used with care, as there are no checks for masslessness.";
ChainToTrace::usage = "Changes spinor chains to trace objects where possible.";

SpinorHelicity6DtoSAM::usage = "Change notation to that of S@M by D. Maitre and P. Mastrolia. Use only if all momenta are four-dimensional.";

ReplaceKinematics::usage = "";
SumCycliclyPermutedKinematics::usage = "";
SumPermutedKinematics::usage = "Sums kinematic expression over all permutations of elements of
list. When passed multiple lists of momentum labels, the summation is carried out over
permutations of each list separately.";

basisS::usage = "";
basisS::gramDegenerate = "The basis of products of `2` Mandelstam for `1` momenta contains linear
dependencies due to gram relations. Implement an exception for this case. Sorry...";

allS::usage = "";

InversePhase::usage = "";
InversePhaseGravity::usage = "";


(* error messages for the exported objects *)

AngAngInvChain::argxMomenta = "No or an odd number of slashed momenta provided.";
AngSquInvChain::argxMomenta = "Odd number of slashed momenta provided.";
SquAngInvChain::argxMomenta = "Odd number of slashed momenta provided.";


Begin["`Private`"];    (* begin the private context (implementation part) *)

protected = Unprotect[SpinorHelicity6D`Chain,SpinorHelicity6D`Mom4DVector];

(* definition of auxiliary functions and local (static) variables *)


(* ::Chapter:: *)
(* Definition of the exported functions *)


(* ::Section:: *)
(* Parke Taylor Factor *)


PTF[n_] := Times @@ (SpinorAngleBracket[#[[1]], #[[2]]] & /@ Partition[Range[n], 2, 1, {1}]);
PTF[n_List] := Times @@ (SpinorAngleBracket[#[[1]], #[[2]]] & /@ Partition[n, 2, 1, {1}]);
PTF[moms__] /; Length@{moms} > 1 := PTF[{moms}];

PTFOpen[n_] := Times @@ (SpinorAngleBracket[#[[1]], #[[2]]] & /@ Partition[Range[n], 2, 1]);
PTFOpen[n_List] := Times @@ (SpinorAngleBracket[#[[1]], #[[2]]] & /@ Partition[n, 2, 1]);
PTFOpen[moms__] /; Length@{moms} > 1 := PTFOpen[{moms}];

PTFn[x_] := PTF[x];


AntiPTF[x___]:=PTF[x]/.{SpinorAngleBracket->SpinorSquareBracket};



(* ::Section::Closed:: *)
(* 6D Spinor Chains *)


AngAngInvChain[n_,ks__,m_][a_,bdot_]/; And[EvenQ[Length@List[ks]],Message[AngAngInvChain::argxMomenta]] := Null;
AngAngInvChain[n_,m_][a_,bdot_] /; Message[AngAngInvChain::argxMomenta] := Null;
AngAngInvChain[n_,ks__,m_][a_,bdot_]:=Module[{dotVar1=Unique[dot],dotVar2=Unique[dot]},
SumContracted[dotVar1,dotVar2][AngSquInvChain[n,Sequence@@Most[List[ks]],Last@List[ks]][a,dotVar1]levicivita2up[dotVar1,dotVar2]SquAngInvariant[Last@List[ks],m][dotVar2,bdot]]
]
AngSquInvChain[n_,ks__,m_][a_,bdot_] /; And[OddQ[Length@List[ks]],Message[AngSquInvChain::argxMomenta]] := Null;
AngSquInvChain[n_,m_][a_,bdot_] := AngSquInvariant[n,m][a,bdot];
AngSquInvChain[n_,ks__,m_][a_,bdot_]:=Module[{dotVar1=Unique[dot],dotVar2=Unique[dot]},
SumContracted[dotVar1,dotVar2][AngAngInvChain[n,Sequence@@Most[List[ks]],Last@List[ks]][a,dotVar1]levicivita2up[dotVar1,dotVar2]AngSquInvariant[Last@List[ks],m][dotVar2,bdot]]
]

SquAngInvChain[n_,ks__,m_][a_,bdot_] /; And[OddQ[Length@List[ks]],Message[SquAngInvChain::argxMomenta]] := Null;
SquAngInvChain[n_,m_][a_,bdot_]:=AngSquInvariant[n,m][a,bdot];
SquAngInvChain[n_,ks__,m_][adot_,b_]:=Module[{dotVar1=Unique[dot],dotVar2=Unique[dot]},
SumContracted[dotVar1,dotVar2][SquAngInvariant[n,First@List[ks]][adot,dotVar1]levicivita2up[dotVar1,dotVar2]AngAngInvChain[First@List[ks],Sequence@@Rest[List[ks]],m][dotVar2,b]]
]



(* ::Section:: *)
(* 4D Traces *)


TrM[args___] /; And[Length[List[args]] != 4, Message[TrM::argrx, "TrM", Length@{args}, 4]] := Null;
TrP[args___] /; And[Length[List[args]] != 4, Message[TrP::argrx, "TrP", Length@{args}, 4]] := Null;
Tr5[args___] /; And[Length[List[args]] != 4, Message[Tr5::argrx, "Tr5", Length@{args}, 4]] := Null;

TrMN[args___] /; And[Length[List[args]] != 4, Message[TrMN::argrx, "TrMN", Length@{args}, 4]] := Null;
TrPN[args___] /; And[Length[List[args]] != 4, Message[TrPN::argrx, "TrPN", Length@{args}, 4]] := Null;
Tr5N[args___] /; And[Length[List[args]] != 4, Message[Tr5N::argrx, "Tr5N", Length@{args}, 4]] := Null;
TrMN[i1_,i2_,i3_,i4_] := Tr[Dot@@Table[Mom4DN[Replace[i,UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,i]||Head[i]===UnderBar,Null,$flat]][If[i===i1||i===i3,$down,$up]],{i,{i1,i2,i3,i4}}]];
TrPN[i1_,i2_,i3_,i4_] := Tr[Dot@@Table[Mom4DN[Replace[i,UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,i]||Head[i]===UnderBar,Null,$flat]][If[i===i2||i===i4,$down,$up]],{i,{i1,i2,i3,i4}}]];
Tr5N[i1_,i2_,i3_,i4_] := TrMN[i1,i2,i3,i4]-TrPN[i1,i2,i3,i4];


TracesToSpinors[expr_] := expr/.{Tr5[x__] :> TrM[x]-TrP[x]}/.{HoldPattern[TrM[sp__]]:>Module[{parts=Partition[List[sp],2,1,{1}]},Times@@(SpinorAngleBracket[#[[1]],#[[2]]]&/@(parts[[1;; ;;2]]))*Times@@(SpinorSquareBracket[#[[1]],#[[2]]]&/@(parts[[2;; ;;2]]))],
HoldPattern[TrP[sp__]]:>Module[{parts=Partition[List[sp],2,1,{1}]},Times@@(SpinorSquareBracket[#[[1]],#[[2]]]&/@(parts[[1;; ;;2]]))*Times@@(SpinorAngleBracket[#[[1]],#[[2]]]&/@(parts[[2;; ;;2]]))]};

ChainToTrace[expr_] := expr /.{Chain[$angle,i_,{x__},i_,$square]/;Length@{x}==3 :>TrM[i,x],Chain[$square,i_,{x__},i_,$angle]/;Length@{x}==3 :>TrP[i,x]};


SpinorHelicity6DtoSAM[expr_] := expr/.{SpinorHelicity6D`SpinorAngleBracket->Spinors`Spaa,SpinorHelicity6D`SpinorSquareBracket->Spinors`Spbb,
	SpinorHelicity6D`Chain[$angle,a_,{b__},c_,$angle]:>Spinors`Spaa[a,b,c],SpinorHelicity6D`Chain[$angle,a_,{b__},c_,$square]:>Spinors`Spab[a,b,c],
	SpinorHelicity6D`Chain[$square,a_,{b__},c_,$angle]:>Spinors`Spba[a,b,c],SpinorHelicity6D`Chain[$square,a_,{b__},c_,$square]:>Spinors`Spbb[a,b,c],S4Mom->Spinors`s,S->Spinors`s};



(* ::Section:: *)
(*Kinematic Cyclic Permutation*)

Options[ReplaceKinematics] = {ExtraHeads->{}};
ReplaceKinematics[expr_,momOrig_List,momRepl_List,opts:OptionsPattern[]]/;Length@momOrig===Length@momRepl	:=
    Module[{labelRules},
	labelRules = MapThread[#1->#2&,{momOrig,momRepl}];
	expr /. {obj: _S|_S6Mom|_SpinorAngleBracket|_SpinorSquareBracket|_Chain|_Tr5|_TrM|_TrP|
		_AngAngInvChain|_AngSquInvChain|_SquAngInvChain|_PTF|_Mom4D|_Mom4DVector|_Extramass
			|_Extramasstilde|_LorentzProduct4D :> (obj/.labelRules)}/.
			{obj:(Blank/@Alternatives@@OptionValue[ExtraHeads]) :> (obj/.labelRules)}
];

Options[SumCycliclyPermutedKinematics] = Options[ReplaceKinematics];
SumCycliclyPermutedKinematics[expr_,mom_List,moms__List,opts:OptionsPattern[]]:=
    SumCycliclyPermutedKinematics[SumCycliclyPermutedKinematics[expr,moms,opts],mom,opts];
SumCycliclyPermutedKinematics[expr_,mom_List,opts:OptionsPattern[]]:=
	Sum[ReplaceKinematics[
		expr,mom,momReplace,
		PassRules[SumCycliclyPermutedKinematics,ReplaceKinematics,opts]
	],
	{momReplace,CyclicPermutations[mom]}
];

Options[SumPermutedKinematics] = Options[ReplaceKinematics];
SumPermutedKinematics[expr_,mom_List,moms__List,opts:OptionsPattern[]]:=
		SumCycliclyPermutedKinematics[SumPermutedKinematics[expr,moms,opts],mom,opts];
SumPermutedKinematics[expr_,mom_List,opts:OptionsPattern[]]:=
	Sum[
		ReplaceKinematics[expr,mom,momReplace,
			PassRules[SumPermutedKinematics,ReplaceKinematics,opts]
		],
		{momReplace,Permutations[mom]}
	];


basisS[n_]:=
    S6Mom @@@ Join[Partition[Range @ n, 2, 1, 1],Sequence@@ Table[Partition[Range @ (n-1), i, 1],{i,3,n-3}]];
basisS[n_,0]:={1};
basisS[n_,1]:=basisS[n];

basisS[6,5]:=
    Delete[SymmetricTensorProductElements[basisS[6],5],{63}];
basisS[6,6]:=
    Delete[SymmetricTensorProductElements[basisS[6],6],{{63}, {183}, {211}, {232}, {237}, {241}, {242}, {243}, {244}}];
basisS[n_,i_] /; ((n == 6 && i <= 6) || (n >= 7 && i <= 4) || n <= 5) :=
    SymmetricTensorProductElements[basisS[n],i];
basisS[n_,i_]  /; Message[basisS::gramDegenerate,n, i] :=
    Null;

allS[n_] := S6Mom @@@ Subsets[Range @ n, {2, n-3}];


InversePhase[plus_List, minus_List] /; Length@plus > 1 && Length@minus > 1 :=
    PTF @@ plus*AntiPTF @@ minus;
InversePhase[plus_List, {minus_}] /; Length@plus > 2 :=
	PTF@@plus*Chain[$square,minus,plus[[{1,2}]],minus,$square];
InversePhase[{plus_}, minus_List] /; Length@minus > 2 :=
		AntiPTF@@minus*Chain[$angle,plus,minus[[{1,2}]],plus,$angle];

InversePhase[plus_List, {}] /; Length@plus >= 3  :=
		PTF @@ plus;
InversePhase[{}, minus_List] /; Length@minus >= 3 :=
		AntiPTF @@ minus;

InversePhaseGravity[args___]:=InversePhase[args]^2;

(* ::Chapter:: *)
(*Postamble*)


(* definitions for system functions *)

Protect[Evaluate[protected]];

End[];


Protect[ "SpinorHelicity6D`SpinorFunctions`*" ];    (* protect exported symbols *)

EndPackage[ ];  (* end the package context *)


(* ::Section:: *)
(**)
