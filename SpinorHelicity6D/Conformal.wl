(* ::Package:: *)

(* :Title: Conformal.m -- a package for conformal operations on spinor expressions *)

(* :Context: SpinorHelicity6D`Conformal` *)

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

*)

(* :Requirements:
   SpinorHelicity6D.wl handling 6D Spinors, written by M. Huber and S. de Angelis
	SpinorFunctions.wl for definition of trace objectss
*)

(* :Todo:
*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`Conformal`", "SpinorHelicity6D`", "SpinorHelicity6D`SpinorFunctions`"];

Unprotect["SpinorHelicity6D`Conformal`*"];
ClearAll["SpinorHelicity6D`Conformal`*"];
ClearAll["SpinorHelicity6D`Conformal`Private`*"];

(* usage messages for the exported functions and the context itself *)

Conformal::usage = "Conformal.wl is a package for conformal operations on spinor expressions.";


La::usage = "Argument for D to take (angle)spinor derivatives. Use as La[i,j], where the derivative for \[Lambda]_i is contracted with \[Lambda]_j.";
Lat::usage = "Argument for D to take (square)spinor derivatives. Use as Lat[i,j], where the derivative for Tilde[\[Lambda]]_i is contracted with Tilde[\[Lambda]]_j.";

KSp::usage = "Special conformal generator in terms of spinor derivatives.";
KSpCrossTerm::usage = "Computes crossterm when acting with the special conformal generator on a product of two expressions.";

(* error messages for the exported objects *)


Begin["`Private`"];    (* begin the private context (implementation part) *)

(* definition of auxiliary functions and local (static) variables *)


(* ::Chapter:: *)
(* Definition of the exported functions *)


SetOptions[D, NonConstants -> {SpinorHelicity6D`SpinorAngleBracket,SpinorHelicity6D`SpinorSquareBracket,SpinorHelicity6D`S,
										SpinorHelicity6D`S4Mom,SpinorHelicity6D`S6Mom,SpinorHelicity6D`SpinorFunctions`TrM,
										SpinorHelicity6D`SpinorFunctions`TrP,SpinorHelicity6D`SpinorFunctions`Tr5,SpinorHelicity6D`Chain}];


La /: D[HoldPattern[SpinorAngleBracket[i_,k_]],La[i_,laRep_],___] := SpinorAngleBracket[laRep,k];
La /: D[HoldPattern[SpinorAngleBracket[j_,i_]],La[i_,laRep_],___] := SpinorAngleBracket[j,laRep];
La /: D[HoldPattern[SpinorAngleBracket[j_,k_]],La[i_,laRep_],___] := 0;
La /: D[HoldPattern[SpinorSquareBracket[__]],La[i_,laRep_],___] := 0;

La /: D[S[i_,k_],La[i_,laRep_],___] := SpinorAngleBracket[laRep,k]SpinorSquareBracket[k,i];
La /: D[S[j_,i_],La[i_,laRep_],___] := SpinorAngleBracket[j,laRep]SpinorSquareBracket[i,j];
La /: D[S[j_,k_],La[i_,laRep_],___] := 0;

La /: D[HoldPattern[S4Mom[k__]],La[i_,laRep_],args___] /; Length@{k} > 2 := Sum[D[S4Mom[Sequence@@subset],La[i,laRep],args],{subset,Subsets[{k},{2}]}];
La /: D[S4Mom[i_,k_],La[i_,laRep_],___] := SpinorAngleBracket[laRep,k]SpinorSquareBracket[k,i];
La /: D[S4Mom[j_,i_],La[i_,laRep_],___] := SpinorAngleBracket[j,laRep]SpinorSquareBracket[i,j];
La /: D[S4Mom[j_,k_],La[i_,laRep_],___] := 0;

La /: D[S6Mom[i_,k_],La[i_,laRep_],___] := SpinorAngleBracket[laRep,k]SpinorSquareBracket[k,i];
La /: D[S6Mom[j_,i_],La[i_,laRep_],___] := SpinorAngleBracket[j,laRep]SpinorSquareBracket[i,j];
La /: D[S6Mom[j_,k_],La[i_,laRep_],___] := 0;

La /: D[HoldPattern[TrM[moms__]],La[i_,laRep_],args___] := D[TracesToSpinors[TrM[moms]],La[i,laRep],args];
La /: D[HoldPattern[TrP[moms__]],La[i_,laRep_],args___] := D[TracesToSpinors[TrP[moms]],La[i,laRep],args];
La /: D[HoldPattern[Tr5[moms__]],La[i_,laRep_],args___] := D[TrM[moms]-TrP[moms],La[i,laRep],args];

La /: D[HoldPattern[Chain[chainargs__]],La[i_,laRep_],args___] := D[ChainToSpinor[Chain[chainargs]],La[i,laRep],args];




Lat /: D[HoldPattern[SpinorSquareBracket[i_,k_]],Lat[i_,latRep_],___] := SpinorSquareBracket[latRep,k];
Lat /: D[HoldPattern[SpinorSquareBracket[j_,i_]],Lat[i_,latRep_],___] := SpinorSquareBracket[j,latRep];
Lat /: D[HoldPattern[SpinorSquareBracket[j_,k_]],Lat[i_,latRep_],___] := 0;
Lat /: D[HoldPattern[SpinorAngleBracket[__]],Lat[i_,latRep_],___] := 0;

Lat /: D[S[i_,k_],Lat[i_,latRep_],___] := SpinorAngleBracket[i,k]SpinorSquareBracket[k,latRep];
Lat /: D[S[j_,i_],Lat[i_,latRep_],___] := SpinorAngleBracket[j,i]SpinorSquareBracket[latRep,j];
Lat /: D[S[j_,k_],Lat[i_,latRep_],___] := 0;

Lat /: D[HoldPattern[S4Mom[k__]],Lat[i_,latRep_],args___] /; Length@{k} > 2 := Sum[D[S4Mom[Sequence@@subset],Lat[i,latRep],args], {subset,Subsets[{k},{2}]}];
Lat /: D[S4Mom[i_,k_],Lat[i_,latRep_],___] := SpinorAngleBracket[i,k]SpinorSquareBracket[k,latRep];
Lat /: D[S4Mom[j_,i_],Lat[i_,latRep_],___] := SpinorAngleBracket[j,i]SpinorSquareBracket[latRep,j];
Lat /: D[S4Mom[j_,k_],Lat[i_,latRep_],___] := 0;

Lat /: D[S6Mom[i_,k_],Lat[i_,latRep_],___] := SpinorAngleBracket[i,k]SpinorSquareBracket[k,latRep];
Lat /: D[S6Mom[j_,i_],Lat[i_,latRep_],___] := SpinorAngleBracket[j,i]SpinorSquareBracket[latRep,j];
Lat /: D[S4Mom[j_,k_],Lat[i_,latRep_],___] := 0;

Lat /: D[HoldPattern[Tr5[moms__]],Lat[i_,latRep_],args___] := D[TrM[moms]-TrP[moms],Lat[i,latRep],args];
Lat /: D[HoldPattern[TrM[moms__]],Lat[i_,latRep_],args___] := D[TracesToSpinors[TrM[moms]],Lat[i,latRep],args];
Lat /: D[HoldPattern[TrP[moms__]],Lat[i_,latRep_],args___] := D[TracesToSpinors[TrP[moms]],Lat[i,latRep],args];

Lat /: D[HoldPattern[Chain[chainargs__]],Lat[i_,latRep_],args___] := D[ChainToSpinor[Chain[chainargs]],Lat[i,latRep],args];



KSp[k_,mom_] := Function[expr,KSp[expr,k,mom]];
KSp[expr_,k_,n_Integer] := KSp[expr,k,Range@n];
KSp[expr_,k_,mom_List] := Sum[D[D[expr,La[i,k]],Lat[i,k]],{i,mom}];

KSpCrossTerm[k_,mom_] := Function[{expr1,expr2},KSpCrossTerm[expr1,expr2,k,mom]];
KSpCrossTerm[expr1_,expr2_,k_,n_Integer] := KSpCrossTerm[expr1,expr2,k,Range@n];
KSpCrossTerm[expr1_,expr2_,k_,mom_List] := Sum[D[expr1,La[i,k]]*D[expr2,Lat[i,k]]+D[expr2,La[i,k]]*D[expr1,Lat[i,k]],{i,mom}];

(* definitions for system functions *)

Protect[Evaluate[protected]];

End[];


Protect[ "SpinorHelicity6D`Conformal`*" ];    (* protect exported symbols *)

EndPackage[ ];  (* end the package context *)
