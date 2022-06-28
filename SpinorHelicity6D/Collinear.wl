(* ::Package:: *)

(* :Title: Collinear.m -- a package for generating collinear momentum configurations *)

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

*)

(* :Requirements:
   SpinorHelicity6D.wl handling 6D Spinors, written by M. Huber and S. de Angelis
	SpinorFunctions.wl for definition of trace objectss
*)

(* :Todo:
*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`Collinear`", "SpinorHelicity6D`"];

Unprotect["SpinorHelicity6D`Collinear`*"];
ClearAll["SpinorHelicity6D`Collinear`*"];
ClearAll["SpinorHelicity6D`Collinear`Private`*"];

(* usage messages for the exported functions and the context itself *)

Collinear::usage = "Collinear.wl is a package for generating collinear momentum configurations.";

GenCollinear::usage = "Generates collinear momentum configurations. Takes as argument Rule[momOrgi,{momCol1,momCol2}] and Rule[{pkIn,piIn},{pkOut,piOut}],
		       where the second rule specifies the momenta used for the recoil.";
		    
GenSoftMasslessMassive::usage = "";

SplittingFunctionTree::usage = "Generates tree-level splitting functions";
SplittingFunctionOneLoop::usage = "Generates one-loop splitting functions";
SoftFunctionTree::usage = "Generates tree-level soft functions";

(* error messages for the exported objects *)


Begin["`Private`"];    (* begin the private context (implementation part) *)

protected = {};

(* definition of auxiliary functions and local (static) variables *)
raiseLowerMom4DNSpinorIndices[mat_/;Dimensions[mat]==={2,2}] := -Transpose[LeviCivitaTensor[2] . mat . LeviCivitaTensor[2]];

cGamma[eps_]:=Exp[eps*EulerGamma]/2*Gamma[1+eps]Gamma[1-eps]^2/Gamma[1-2*eps];

ClearAll[rSSYM];
rSSYM[z_,s_,mu_,eps_, cutoff_ : 1]:=cGamma[eps](-mu^2/s)^eps*1/eps^2(-((1-z)/z)^eps*Pi*eps/(Sin[Pi*eps])+Sum[2*eps^(2*m-1)PolyLog[2m-1,z/(z-1)],{m,1,cutoff}])

Options[rSQCDPP] = {DeltaHV->0,NFoNC->0};
rSQCDPP[z_,s_,mu_,eps_,opts:OptionsPattern[]]:=rSSYM[z,s,mu,eps]+cGamma[eps](-mu^2/s)^eps*2*z*(1-z)/((1-2*eps)(2-2*eps)(3-2*eps))(1-eps*OptionValue[DeltaHV]-OptionValue[NFoNC]);

(* Definition of the exported functions *)



Options[GenCollinear] = {Reference -> q, AngleVar -> ArcCos[Sqrt[z]], LimitVar-> \[Epsilon]Coli, MomentumConserving->True};

GenCollinear[Rule[momOrig_,{momCol1_,momCol2_}], opts : OptionsPattern[]] := GenCollinear[Rule[momOrig,{momCol1,momCol2}],Rule[{Null,Null},{Null,Null}],MomentumConserving->False,opts];

GenCollinear[Rule[momOrig_,{momCol1_,momCol2_}],Rule[{pk_,pi_},{pkOut_,piOut_}], opts : OptionsPattern[]]:=Module[{protected,gamma},

	Mom4DN[momCol1][Null][$up] = Mom4DN[momOrig][Null][$up]*Cos[OptionValue[AngleVar]]^2+OptionValue[LimitVar]^2*Sin[OptionValue[AngleVar]]^2*Mom4DN[OptionValue[Reference]][Null][$up]-
					OptionValue[LimitVar]Sin[OptionValue[AngleVar]]*Cos[OptionValue[AngleVar]]*
					(KroneckerProduct[SpinorUndotN[momOrig][$lam][$up],SpinorDotN[OptionValue[Reference]][$lam][$up]]+KroneckerProduct[SpinorUndotN[OptionValue[Reference]][$lam][$up],SpinorDotN[momOrig][$lam][$up]]);
	SpinorHelicity6D`Mom4DN[momCol1][$flat][$up] = Mom4DN[momCol1][Null][$up];
	SpinorHelicity6D`Mom4DN[momCol1][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momCol1][Null][$up]];
	SpinorHelicity6D`Mom4DN[momCol1][$flat][$down] = Mom4DN[momCol1][Null][$down];
	Mom4DVectorNFromMatrix[momCol1][Null];
	Mom4DVectorNFromMatrix[momCol1][$flat];
	
	SpinorHelicity6D`SpinorDotN[momCol1][$lam][$up] = SpinorDotN[momOrig][$lam][$up]*Cos[OptionValue[AngleVar]]-SpinorDotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol1][$lam][$up] = SpinorUndotN[momOrig][$lam][$up]*Cos[OptionValue[AngleVar]]-SpinorUndotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol1][$lam][$down] = SpinorDotN[momOrig][$lam][$down]*Cos[OptionValue[AngleVar]]-SpinorDotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol1][$lam][$down] = SpinorUndotN[momOrig][$lam][$down]*Cos[OptionValue[AngleVar]]-SpinorUndotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Sin[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol1][$mu][$up] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol1][$mu][$up] = 0;
	SpinorHelicity6D`SpinorDotN[momCol1][$mu][$down] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol1][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	SpinorHelicity6D`Extramass[momCol1] = 0;
	SpinorHelicity6D`Extramasstilde[momCol1] = 0;
	Protect[Extramass,Extramasstilde];
	SpinorHelicity6D`ExtramassN[momCol1] = 0;
	SpinorHelicity6D`ExtramasstildeN[momCol1] = 0;



	Mom4DN[momCol2][Null][$up] = Mom4DN[momOrig][Null][$up]*Sin[OptionValue[AngleVar]]^2+OptionValue[LimitVar]^2*Cos[OptionValue[AngleVar]]^2*Mom4DN[OptionValue[Reference]][Null][$up]+
					OptionValue[LimitVar]Sin[OptionValue[AngleVar]]*Cos[OptionValue[AngleVar]]*
					(KroneckerProduct[SpinorUndotN[momOrig][$lam][$up],SpinorDotN[OptionValue[Reference]][$lam][$up]]+KroneckerProduct[SpinorUndotN[OptionValue[Reference]][$lam][$up],SpinorDotN[momOrig][$lam][$up]]);
	SpinorHelicity6D`Mom4DN[momCol2][$flat][$up] = Mom4DN[momCol2][Null][$up];
	SpinorHelicity6D`Mom4DN[momCol2][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momCol2][Null][$up]];
	SpinorHelicity6D`Mom4DN[momCol2][$flat][$down] = Mom4DN[momCol2][Null][$down];
	Mom4DVectorNFromMatrix[momCol2][Null];
	Mom4DVectorNFromMatrix[momCol2][$flat];
	
	SpinorHelicity6D`SpinorDotN[momCol2][$lam][$up] = SpinorDotN[momOrig][$lam][$up]*Sin[OptionValue[AngleVar]]+SpinorDotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol2][$lam][$up] = SpinorUndotN[momOrig][$lam][$up]*Sin[OptionValue[AngleVar]]+SpinorUndotN[OptionValue[Reference]][$lam][$up]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol2][$lam][$down] = SpinorDotN[momOrig][$lam][$down]*Sin[OptionValue[AngleVar]]+SpinorDotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorUndotN[momCol2][$lam][$down] = SpinorUndotN[momOrig][$lam][$down]*Sin[OptionValue[AngleVar]]+SpinorUndotN[OptionValue[Reference]][$lam][$down]*OptionValue[LimitVar]Cos[OptionValue[AngleVar]];
	SpinorHelicity6D`SpinorDotN[momCol2][$mu][$up] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol2][$mu][$up] = 0;
	SpinorHelicity6D`SpinorDotN[momCol2][$mu][$down] = 0;
	SpinorHelicity6D`SpinorUndotN[momCol2][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	SpinorHelicity6D`Extramass[momCol2] = 0;
	SpinorHelicity6D`Extramasstilde[momCol2] = 0;
	Protect[Extramass,Extramasstilde];
	SpinorHelicity6D`ExtramassN[momCol2] = 0;
	SpinorHelicity6D`ExtramasstildeN[momCol2] = 0;
	

	If[OptionValue[MomentumConserving],

	   gamma = 1 + (-OptionValue[LimitVar]^2*ChainN[$angle,pi,{OptionValue[Reference]},pi,$square])/(ChainN[$angle,pi,{pk},pi,$square]-OptionValue[LimitVar]^2ChainN[$angle,pk,{OptionValue[Reference]},pk,$square]);
	
		SpinorHelicity6D`Mom4DN[pkOut][Null][$up] = gamma * Mom4DN[pk][Null][$up];
		SpinorHelicity6D`Mom4DN[pkOut][$flat][$up] = Mom4DN[pkOut][Null][$up];
		SpinorHelicity6D`Mom4DN[pkOut][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[pkOut][Null][$up]];
		SpinorHelicity6D`Mom4DN[pkOut][$flat][$down] = Mom4DN[pkOut][Null][$down];
		Mom4DVectorNFromMatrix[pkOut][Null];
		Mom4DVectorNFromMatrix[pkOut][$flat];
		SpinorsFromMom4DN[pkOut,Normalization->pk];
		Unprotect[Extramass,Extramasstilde];
		SpinorHelicity6D`Extramass[pkOut] = 0;
		SpinorHelicity6D`Extramasstilde[pkOut] = 0;
		Protect[Extramass,Extramasstilde];
		SpinorHelicity6D`ExtramassN[pkOut] = 0;
		SpinorHelicity6D`ExtramasstildeN[pkOut] = 0;
	   
		SpinorHelicity6D`Mom4DN[piOut][Null][$up] = Mom4DN[pi][Null][$up] - OptionValue[LimitVar]^2 Mom4DN[OptionValue[Reference]][Null][$up] + (1-gamma) * Mom4DN[pk][Null][$up];
		SpinorHelicity6D`Mom4DN[piOut][$flat][$up] = Mom4DN[piOut][Null][$up];
		SpinorHelicity6D`Mom4DN[piOut][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[piOut][Null][$up]];
		SpinorHelicity6D`Mom4DN[piOut][$flat][$down] = Mom4DN[piOut][Null][$down];
		Mom4DVectorNFromMatrix[piOut][Null];
		Mom4DVectorNFromMatrix[piOut][$flat];
		SpinorsFromMom4DN[piOut,Normalization->pi];
		Unprotect[Extramass,Extramasstilde];
		SpinorHelicity6D`Extramass[piOut] = 0;
		SpinorHelicity6D`Extramasstilde[piOut] = 0;
		Protect[Extramass,Extramasstilde];
		SpinorHelicity6D`ExtramassN[piOut] = 0;
		SpinorHelicity6D`ExtramasstildeN[piOut] = 0;
	];
];




Options[GenSoftMasslessMassive] = {LimitVar-> \[Epsilon]Soft};

GenSoftMasslessMassive[Rule[mom1_,momSoft1_],Rule[mom2_,momRecoil2_], opts : OptionsPattern[]]:=
Module[{protected},

	Mom4DN[momSoft1][Null][$up] = OptionValue[LimitVar]^2*Mom4DN[mom1][Null][$up];
	Mom4DN[momSoft1][$flat][$up] = Mom4DN[momSoft1][Null][$up];
	Mom4DN[momSoft1][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momSoft1][Null][$up]];
	Mom4DN[momSoft1][$flat][$down] = Mom4DN[momSoft1][Null][$down];
	Mom4DVectorNFromMatrix[momSoft1][Null];
	Mom4DVectorNFromMatrix[momSoft1][$flat];
	
	SpinorDotN[momSoft1][$lam][$up] = OptionValue[LimitVar]*SpinorDotN[mom1][$lam][$up];
	SpinorUndotN[momSoft1][$lam][$up] = OptionValue[LimitVar]*SpinorUndotN[mom1][$lam][$up];
	SpinorDotN[momSoft1][$lam][$down] = OptionValue[LimitVar]*SpinorDotN[mom1][$lam][$down];
	SpinorUndotN[momSoft1][$lam][$down] = OptionValue[LimitVar]*SpinorUndotN[mom1][$lam][$down];
	{SpinorDotN[momSoft1][$lam][1][Null],SpinorDotN[momSoft1][$lam][2][Null]} = SpinorDotN[momSoft1][$lam][$up];
	{SpinorUndotN[momSoft1][$lam][1][Null],SpinorUndotN[momSoft1][$lam][2][Null]} = SpinorUndotN[momSoft1][$lam][$up];
	{SpinorDotN[momSoft1][$lam][Null][1],SpinorDotN[momSoft1][$lam][Null][2]} = SpinorDotN[momSoft1][$lam][$down];
	{SpinorUndotN[momSoft1][$lam][Null][1],SpinorUndotN[momSoft1][$lam][Null][2]} = SpinorUndotN[momSoft1][$lam][$down];
	
	SpinorDotN[momSoft1][$mu][$up] = 0;
	SpinorUndotN[momSoft1][$mu][$up] = 0;
	SpinorDotN[momSoft1][$mu][$down] = 0;
	SpinorUndotN[momSoft1][$mu][$down] = 0;
	Unprotect[Extramass,Extramasstilde];
	Extramass[momSoft1] = 0;
	Extramasstilde[momSoft1] = 0;
	Protect[Extramass,Extramasstilde];
	ExtramassN[momSoft1] = 0;
	ExtramasstildeN[momSoft1] = 0;

	Mom4DN[momRecoil2][Null][$up] = Mom4DN[mom2][Null][$up]+(1-OptionValue[LimitVar]^2)*Mom4DN[mom1][Null][$up];
	Mom4DN[momRecoil2][$flat][$up] = Mom4DN[momRecoil2][Null][$up];
	Mom4DN[momRecoil2][Null][$down] = raiseLowerMom4DNSpinorIndices[Mom4DN[momRecoil2][Null][$up]];
	Mom4DN[momRecoil2][$flat][$down] = Mom4DN[momRecoil2][Null][$down];
	Mom4DVectorNFromMatrix[momRecoil2][Null];
	Mom4DVectorNFromMatrix[momRecoil2][$flat];
	
	Unprotect[Extramass,Extramasstilde];
	Extramass[momRecoil2] = Extramass[mom2];
	Extramasstilde[momRecoil2] = Extramasstilde[mom2];
	Protect[Extramass,Extramasstilde];
	ExtramassN[momRecoil2] = ExtramassN[mom2];
	ExtramasstildeN[momRecoil2] = ExtramasstildeN[mom2];
	
];




SplittingFunctionTree[z_,a_,b_]:=
	<|"-"->
	      Association[{"-","-"}->0
			 ,{"+","+"}->1/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b],
			  {"+","-"}->-z^2/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b],
			  {"-","+"}->-(1-z)^2/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b]],
	  "+"->
	      Association[{"+","+"}->0,
			  {"-","-"}->-1/Sqrt[z(1-z)]1/SpinorSquareBracket[a,b],
			  {"-","+"}->+z^2/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b],
			  {"+","-"}->(1-z)^2/Sqrt[z(1-z)]1/SpinorAngleBracket[a,b]
	      ]
|>;


Options[SplittingFunctionOneLoop] = {DeltaHV->0,NFoNC->0,MuVar->Global`\[Mu],EpsVar->Global`\[Epsilon]};
				   
SplittingFunctionOneLoop[z_,a_,b_,opts :OptionsPattern[]]:=
	2*<|"-"->
		Association[
			{"-","-"}->cGamma[OptionValue[EpsVar]]Sqrt[z(1-z)]Spinoranglebracket[a,b]/SpinorSquareBracket[a,b]^2*
					2/((1-2*OptionValue[EpsVar])(2-2*OptionValue[EpsVar])(3-2*OptionValue[EpsVar]))*(-OptionValue[MuVar]^2/S[a,b])^OptionValue[EpsVar]*(1-OptionValue[EpsVar]*OptionValue[DeltaHV]-OptionValue[NFoNC]),
			{"+","+"}->SplittingFunctionTree[z,a,b]["-",{"+","+"}]*rSQCDPP[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"-","+"}->SplittingFunctionTree[z,a,b]["-",{"-","+"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"+","-"}->SplittingFunctionTree[z,a,b]["-",{"+","-"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]]
		],
	    "+"->
		Association[
			{"+","+"}->-cGamma[OptionValue[EpsVar]]Sqrt[z(1-z)]SpinorSquareBracket[a,b]/SpinorAngleBracket[a,b]^2*
					2/((1-2*OptionValue[EpsVar])(2-2*OptionValue[EpsVar])(3-2*OptionValue[EpsVar]))*(-OptionValue[MuVar]^2/S[a,b])^OptionValue[EpsVar]*(1-OptionValue[EpsVar]*OptionValue[DeltaHV]-OptionValue[NFoNC]),
			{"-","-"}->SplittingFunctionTree[z,a,b]["+",{"-","-"}]*rSQCDPP[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"+","-"}->SplittingFunctionTree[z,a,b]["+",{"+","-"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]],
			{"-","+"}->SplittingFunctionTree[z,a,b]["+",{"-","+"}]*rSSYM[z,S[a,b],OptionValue[MuVar],OptionValue[EpsVar]]
		]
|>;


(* definitions for system functions *)

Protect[Evaluate[protected]];

End[];


Protect[ "SpinorHelicity6D`Collinear`*" ];    (* protect exported symbols *)

EndPackage[ ];  (* end the package context *)
