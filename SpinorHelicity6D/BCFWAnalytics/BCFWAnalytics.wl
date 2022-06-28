(* ::Package:: *)

(* :Title: BCFW.m -- a package for performing BCFW recursion computations of amplitudes, both analytically and numerically *)

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


(* ::Chapter:: *)
(*Preamble*)


(* set up the package context, including public imports *)

BeginPackage["SpinorHelicity6D`BCFWAnalytics`", "SpinorHelicity6D`SpinorFunctions`","SpinorHelicity6D`"];

Unprotect["SpinorHelicity6D`BCFWAnalytics`*"];
ClearAll["SpinorHelicity6D`BCFWAnalytics`*"];
ClearAll["SpinorHelicity6D`BCFWAnalytics`Private`*"];


(* ::Section::Closed:: *)
(*Old Routines*)


(* usage messages for the exported functions and the context itself *)

BCFWSpinorRules::usage = "Analytic rules for spinors in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWChainRules::usage = "Analytic rules for spinor chains in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWz::usage = "Returns solution for shift parameter z in a standard BCFW, i.e. massless, two-momenta shift.";

BCFWRulesMasslessChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massless factorization channel.";
BCFWRulesMassiveChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massive factorization channel.";



BCFWGluonScalarSpinorRules::usage = "Analytic rules for spinors in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWGluonScalarChainRules::usage = "Analytic rules for spinor chains in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWGluonScalarz::usage = "Returns solution for shift parameter z in a standard BCFW, i.e. massless, two-momenta shift.";

BCFWGluonScalarRulesMasslessChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massless factorization channel.";
BCFWGluonScalarRulesMassiveChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massive factorization channel.";


BCFWScalarGluonSpinorRules::usage = "Analytic rules for spinors in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWScalarGluonChainRules::usage = "Analytic rules for spinor chains in a standard BCFW, i.e. massless, two-momenta shift.";
BCFWScalarGluonz::usage = "Returns solution for shift parameter z in a standard BCFW, i.e. massless, two-momenta shift.";

BCFWScalarGluonRulesMasslessChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massless factorization channel.";
BCFWScalarGluonRulesMassiveChannel::usage = "Applies all appropriate BCFW rules for a massless two-momenta shift with a massive factorization channel.";


MassiveKHatRules::usage = "Analytic rules to replace a massive channel momentum Overhat[K].";
MasslessKHatRules::usage = "Analytic rules to replace a massless channel momentum Overhat[K], as well as its spinors in terms of their projection on the spinors of two reference momenta.";


qBCFW = Global`q;


(* error messages for the exported objects *)


(* ::Section::Closed:: *)
(*New Routines*)


ApplyBcfwScalarGluonShiftRules::usage = "";
SortChains::usage = "";
AlignExtramassRules::usage = "";
ReduceMomenta::usage = "";


(* ::Chapter:: *)
(*Main*)


Begin["`Private`"];    (* begin the private context (implementation part) *)
protected = {};
Unprotect[protected];


(* ::Section::Closed:: *)
(*Old Routines*)


(* definition of auxiliary functions and local (static) variables *)

(* Definition of the exported functions *)

(* Generic rules for KHat *)

MassiveKHatRules[momsK_List]:={
			Chain[a__,{b___,UnderBar[OverHat[K]],c___},d__] :> Sum[Chain[a,{b,UnderBar[ki],c},d],{ki,momsK}],
			S6Mom[OverHat[K],x__]:>S6Mom[Sequence@@momsK,x],
			S6Mom[pm[OverHat[K]],x__]:>S6Mom[Sequence@@(pm /@ momsK),x],
			S4Mom[OverHat[K],x__]:>S4Mom[Sequence@@momsK,x],
			S4Mom[pm[OverHat[K]],x__]:>S4Mom[Sequence@@(pm /@ momsK),x],
			Extramass[OverHat[K]] -> Sum[Extramass[ki],{ki,momsK}],
			Extramasstilde[OverHat[K]] -> Sum[Extramasstilde[ki],{ki,momsK}]
		};

MasslessKHatRules[k1_,k2_,momsK_List]:=Module[{a1,a2,b1,b2},
			a1 = -Sum[Chain[$angle,k2,{UnderBar[ki]},k1,$square],{ki,momsK}]/S[k1,k2];
			a2 =  Sum[Chain[$angle,k1,{UnderBar[ki]},k1,$square],{ki,momsK}]/S[k1,k2];
			b1 = -Sum[Chain[$angle,k2,{UnderBar[ki]},k2,$square],{ki,momsK}]/Sum[Chain[$angle,k2,{UnderBar[ki]},k1,$square],{ki,momsK}];
			b2 = 1;
			
			{HoldPattern[Chain[a__,{b___,UnderBar[OverHat[K]]|OverHat[K],c___},d__]] :> Sum[Chain[a,{b,UnderBar[ki],c},d],{ki,momsK}],
			
			 HoldPattern[SpinorSquareBracket[OverHat[K],x_]]  :> b1*SpinorSquareBracket[k1,x] + b2*SpinorSquareBracket[k2,x],
		 	HoldPattern[SpinorSquareBracket[x_,OverHat[K]]] :> b1*SpinorSquareBracket[x,k1] + b2*SpinorSquareBracket[x,k2],
	 		HoldPattern[SpinorAngleBracket[OverHat[K],x_]]  :> a1*SpinorAngleBracket[k1,x] + a2*SpinorAngleBracket[k2,x],
 			HoldPattern[SpinorAngleBracket[x_,OverHat[K]]]  :> a1*SpinorAngleBracket[x,k1] + a2*SpinorAngleBracket[x,k2],
		 	HoldPattern[Chain[$angle,OverHat[K],{x__},c__]] :> a1*Chain[$angle,k1,{x},c] + a2*Chain[$angle,k2,{x},c],
		 	HoldPattern[Chain[$square,OverHat[K],{x__},c__]] :> b1*Chain[$square,k1,{x},c] + b2*Chain[$square,k2,{x},c],
		 	HoldPattern[Chain[a__,{x__},OverHat[K],$angle]] :> a1*Chain[a,{x},k1,$angle] + a2*Chain[a,{x},k2,$angle],
		 	HoldPattern[Chain[a__,{x__},OverHat[K],$square]] :> b1*Chain[a,{x},k1,$square] + b2*Chain[a,{x},k2,$square],
	 		S6Mom[OverHat[K],x__]:>S6Mom[Sequence@@momsK,x],
			S6Mom[pm[OverHat[K]],x__]:>S6Mom[Sequence@@(pm /@ momsK),x],
			S4Mom[OverHat[K],x__]:>S4Mom[Sequence@@momsK,x],
			S4Mom[pm[OverHat[K]],x__]:>S4Mom[Sequence@@(pm /@ momsK),x],
			Extramass[OverHat[K]] -> Sum[Extramass[ki],{ki,momsK}],
			Extramasstilde[OverHat[K]] -> Sum[Extramasstilde[ki],{ki,momsK}]
			}
		];



(****************************************************)
(* Shift on two massless legs, i.e. the standard BCFW shift *)
(****************************************************)

BCFWSpinorRules[ang_,squ_]:={
				HoldPattern[SpinorSquareBracket[OverHat[ang],x_]] :> SpinorSquareBracket[ang,x],
				HoldPattern[SpinorSquareBracket[x_,OverHat[ang]]] :> SpinorSquareBracket[x,ang],
				HoldPattern[SpinorAngleBracket[OverHat[squ],x_]] :> SpinorAngleBracket[squ,x],
				HoldPattern[SpinorAngleBracket[x_,OverHat[squ]]] :> SpinorAngleBracket[x,squ],
				
				HoldPattern[SpinorAngleBracket[OverHat[ang],x_]] :> SpinorAngleBracket[ang,x]-Global`z*SpinorAngleBracket[squ,x],
				HoldPattern[SpinorAngleBracket[x_,OverHat[ang]]] :> SpinorAngleBracket[x,ang]-Global`z*SpinorAngleBracket[x,squ],
				HoldPattern[SpinorSquareBracket[OverHat[squ],x_]] :> SpinorSquareBracket[squ,x]+Global`z*SpinorSquareBracket[ang,x],
				HoldPattern[SpinorSquareBracket[x_,OverHat[squ]]] :> SpinorSquareBracket[x,squ]+Global`z*SpinorSquareBracket[x,ang],
				
				HoldPattern[Chain[$angle,OverHat[ang],{x__},c__]] :> Chain[$angle,ang,{x},c]-Global`z*Chain[$angle,squ,{x},c],
				HoldPattern[Chain[$square,OverHat[ang],{x__},c__]] :> Chain[$square,ang,{x},c],
				HoldPattern[Chain[a__,{x__},OverHat[ang],$angle]] :> Chain[a,{x},ang,$angle]-Global`z*Chain[a,{x},squ,$angle],
				HoldPattern[Chain[a__,{x__},OverHat[ang],$square]] :> Chain[a,{x},ang,$square],
				
				HoldPattern[Chain[$square,OverHat[squ],{x__},c__]] :> Chain[$square,squ,{x},c]+Global`z*Chain[$square,ang,{x},c],
				HoldPattern[Chain[$angle,OverHat[squ],{x__},c__]] :> Chain[$angle,squ,{x},c],
				HoldPattern[Chain[a__,{x__},OverHat[squ],$square]] :> Chain[a,{x},squ,$square]+Global`z*Chain[a,{x},ang,$square],
				HoldPattern[Chain[a__,{x__},OverHat[squ],$angle]] :> Chain[a,{x},squ,$angle],
				
				Extramass[OverHat[x_]] :> Extramass[x],
				Extramasstilde[OverHat[x_]] :> Extramasstilde[x]
				
			};


BCFWChainRules[ang_,squ_]:={
				HoldPattern[Chain[$angle,a_,{b___,OverHat[ang]|UnderBar[OverHat[ang]],c___},d__]] :> Chain[$angle,a,{b,ang,c},d]-Global`z*If[EvenQ@Length@{b},
					Chain[$angle,a,{b},squ,$angle]Chain[$square,ang,{c},d]
				, 
					Chain[$angle,a,{b},ang,$square]Chain[$angle,squ,{c},d]
				],
				HoldPattern[Chain[$square,a_,{b___,OverHat[ang]|UnderBar[OverHat[ang]],c___},d__]] :> Chain[$square,a,{b,ang,c},d]-Global`z*If[OddQ@Length@{b},
					Chain[$square,a,{b},squ,$angle]Chain[$square,ang,{c},d]
					, 
					Chain[$square,a,{b},ang,$square]Chain[$angle,squ,{c},d]
				],
				HoldPattern[Chain[$angle,a_,{b___,OverHat[squ]|UnderBar[OverHat[squ]],c___},d__]] :> Chain[$angle,a,{b,squ,c},d]+Global`z*If[EvenQ@Length@{b},
					Chain[$angle,a,{b},squ,$angle]Chain[$square,ang,{c},d]
					, 
					Chain[$angle,a,{b},ang,$square]Chain[$angle,squ,{c},d]
				],
				HoldPattern[Chain[$square,a_,{b___,OverHat[squ]|UnderBar[OverHat[squ]],c___},d__]] :> Chain[$square,a,{b,squ,c},d]+Global`z*If[OddQ@Length@{b},
					Chain[$square,a,{b},squ,$angle]Chain[$square,ang,{c},d]
					, 
					Chain[$square,a,{b},ang,$square]Chain[$angle,squ,{c},d]
				],
				S6Mom[i___,OverHat[ang],j___,OverHat[squ],k___]:>S6Mom[i,ang,j,squ,k],
				S6Mom[OverHat[ang],j__] :> S6Mom[ang,j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S6Mom[OverHat[squ],j__] :> S6Mom[squ,j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S6Mom[pm[OverHat[ang]],j__] :> S6Mom[pm[ang],j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S6Mom[pm[OverHat[squ]],j__] :> S6Mom[pm[squ],j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S4Mom[i___,OverHat[ang],j___,OverHat[squ],k___]:>S4Mom[i,ang,j,squ,k],
				S4Mom[OverHat[ang],j__] :> S4Mom[ang,j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S4Mom[OverHat[squ],j__] :> S4Mom[squ,j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S4Mom[pm[OverHat[ang]],j__] :> S4Mom[pm[ang],j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				S4Mom[pm[OverHat[squ]],j__] :> S4Mom[pm[squ],j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k]},squ,$angle],{k,{j}}],
				Chain[$angle,i_,{},j_,$angle]:>SpinorAngleBracket[i,j],
				Chain[$square,i_,{},j_,$square]:>SpinorSquareBracket[i,j]
			};



BCFWz[ang_,squ_,channel_List] /; MemberQ[channel,ang] := +S6Mom@@channel/Sum[Chain[$angle,squ,{UnderBar[ki]},ang,$square],{ki,Select[channel,# =!= ang &]}];
BCFWz[ang_,squ_,channel_List] /; MemberQ[channel,squ] := -S6Mom@@channel/Sum[Chain[$angle,squ,{UnderBar[ki]},ang,$square],{ki,Select[channel,# =!= squ &]}];


(* BCFWRulesMassless and BCFWRulesMassive take hatted momenta in the list momsK defining the factorization channel *)
BCFWRulesMasslessChannel[expr_,ang_,squ_,momsK_List] := BCFWRulesMasslessChannel[ang,squ,momsK][expr];
BCFWRulesMassiveChannel[expr_,ang_,squ_,momsK_List] := BCFWRulesMassiveChannel[ang,squ,momsK][expr];
BCFWRulesMasslessChannel[ang_,squ_,momsK_List] := Function[expr,expr//.MasslessKHatRules[ang,squ,momsK]//.BCFWChainRules[ang,squ]//.BCFWSpinorRules[ang,squ]/.{Global`z->BCFWz[ang,squ,momsK/.{OverHat->Identity}]}];
BCFWRulesMassiveChannel[ang_,squ_,momsK_List] := Function[expr,expr//.MassiveKHatRules[momsK]//.BCFWChainRules[ang,squ]//.BCFWSpinorRules[ang,squ]/.{Global`z->BCFWz[ang,squ,momsK/.{OverHat->Identity}]}];



(****************************************************)
(* Shift on massless and massive leg. For the       *)
(* massless momentum the angle spinor is shifted    *)
(****************************************************)

BCFWGluonScalarSpinorRules[ang_,mass_]:={
				HoldPattern[SpinorSquareBracket[OverHat[ang],x_]] :> SpinorSquareBracket[ang,x],
				HoldPattern[SpinorSquareBracket[x_,OverHat[ang]]] :> SpinorSquareBracket[x,ang],
				
				HoldPattern[SpinorAngleBracket[OverHat[ang],x_]] :> SpinorAngleBracket[ang,x]-Global`z*Chain[$square,ang,{UnderBar[mass]},x,$angle],
				HoldPattern[SpinorAngleBracket[x_,OverHat[ang]]] :> SpinorAngleBracket[x,ang]+Global`z*Chain[$angle,x,{UnderBar[mass]},ang,$square],
								
				HoldPattern[Chain[$angle,OverHat[ang],{x__},c__]] :> Chain[$angle,ang,{x},c]-Global`z*Chain[$square,ang,{UnderBar[mass],x},c],
				HoldPattern[Chain[$square,OverHat[ang],{x__},c__]] :> Chain[$square,ang,{x},c],
				HoldPattern[Chain[a__,{x__},OverHat[ang],$angle]] :> Chain[a,{x},ang,$angle]+Global`z*Chain[a,{x,UnderBar[mass]},ang,$square],
				HoldPattern[Chain[a__,{x__},OverHat[ang],$square]] :> Chain[a,{x},ang,$square],
				
				Extramass[OverHat[x_]] :> Extramass[x],
				Extramasstilde[OverHat[x_]] :> Extramasstilde[x]
				
			};


BCFWGluonScalarChainRules[ang_,mass_]:={
	HoldPattern[Chain[$angle,a_,{b___,OverHat[ang]|UnderBar[OverHat[ang]],c___},d__]] :> Chain[$angle,a,{b,ang,c},d]+Global`z*
		   If[EvenQ@Length@{b},
		      Chain[$angle,a,{b,UnderBar[mass]},ang,$square]Chain[$square,ang,{c},d]
		    , 
		      -Chain[$angle,a,{b},ang,$square]Chain[$square,ang,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$square,a_,{b___,OverHat[ang]|UnderBar[OverHat[ang]],c___},d__]] :> Chain[$square,a,{b,ang,c},d]+Global`z*
		   If[OddQ@Length@{b},
		      Chain[$square,a,{b,UnderBar[mass]},ang,$square]Chain[$square,ang,{c},d]
		    , 
		      -Chain[$square,a,{b},ang,$square]Chain[$square,ang,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$angle,a_,{b___,OverHat[mass]|UnderBar[OverHat[mass]],c___},d__]] :> Chain[$angle,a,{b,UnderBar[mass],c},d]-Global`z*
		   If[EvenQ@Length@{b},
		      Chain[$angle,a,{b,UnderBar[mass]},ang,$square]Chain[$square,ang,{c},d]
		    , 
		      -Chain[$angle,a,{b},ang,$square]Chain[$square,ang,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$square,a_,{b___,OverHat[mass]|UnderBar[OverHat[mass]],c___},d__]] :> Chain[$square,a,{b,UnderBar[mass],c},d]-Global`z*
		   If[OddQ@Length@{b},
		      Chain[$square,a,{b,UnderBar[mass]},ang,$square]Chain[$square,ang,{c},d]
		    , 
		      -Chain[$square,a,{b},ang,$square]Chain[$square,ang,{UnderBar[mass],c},d]
		   ],
	S6Mom[i___,OverHat[ang],j___,OverHat[mass],k___]:>S6Mom[i,ang,j,mass,k],
	S6Mom[OverHat[ang],j__] :> S6Mom[ang,j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S6Mom[OverHat[mass],j__] :> S6Mom[mass,j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S6Mom[pm[OverHat[ang]],j__] :> S6Mom[pm[ang],j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S6Mom[pm[OverHat[mass]],j__] :> S6Mom[pm[mass],j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S4Mom[i___,OverHat[ang],j___,OverHat[mass],k___]:>S4Mom[i,ang,j,mass,k],
	S4Mom[OverHat[ang],j__] :> S4Mom[ang,j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S4Mom[OverHat[mass],j__] :> S4Mom[mass,j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S4Mom[pm[OverHat[ang]],j__] :> S4Mom[pm[ang],j]-Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	S4Mom[pm[OverHat[mass]],j__] :> S4Mom[pm[mass],j]+Global`z*Sum[Chain[$square,ang,{UnderBar[k],UnderBar[mass]},ang,$square],{k,{j}}],
	Chain[$angle,i_,{},j_,$angle]:>SpinorAngleBracket[i,j],
	Chain[$square,i_,{},j_,$square]:>SpinorSquareBracket[i,j]
	};


BCFWGluonScalarz[ang_,mass_,channel_List] /; MemberQ[channel,ang] := +S6Mom@@channel/Sum[Chain[$square,ang,{UnderBar[mass],UnderBar[ki]},ang,$square],{ki,Select[channel,# =!= ang &]}];
BCFWGluonScalarz[ang_,mass_,channel_List] /; MemberQ[channel,mass] := -S6Mom@@channel/Sum[Chain[$square,ang,{UnderBar[mass],UnderBar[ki]},ang,$square],{ki,Select[channel,# =!= mass &]}];


(* BCFWRulesMassless and BCFWRulesMassive take hatted momenta in the list momsK defining the factorization channel *)
BCFWGluonScalarRulesMasslessChannel[expr_,ang_,mass_,momsK_List] := BCFWGluonScalarRulesMasslessChannel[ang,mass,momsK][expr];
BCFWGluonScalarRulesMassiveChannel[expr_,ang_,mass_,momsK_List] := BCFWGluonScalarRulesMassiveChannel[ang,mass,momsK][expr];
BCFWGluonScalarRulesMasslessChannel[ang_,mass_,momsK_List] := Function[expr,expr//.MasslessKHatRules[ang,qBCFW,momsK]//.BCFWGluonScalarChainRules[ang,mass]//.BCFWGluonScalarSpinorRules[ang,mass]/.{Global`z->BCFWGluonScalarz[ang,mass,momsK/.{OverHat->Identity}]}];
BCFWGluonScalarRulesMasslessChannel[ang_,mass_,momsK_List,{pSp1_,pSp2_}] := Function[expr,expr//.MasslessKHatRules[pSp1,pSp2,momsK]//.BCFWGluonScalarChainRules[ang,mass]//.BCFWGluonScalarSpinorRules[ang,mass]/.{Global`z->BCFWGluonScalarz[ang,mass,momsK/.{OverHat->Identity}]}];
BCFWGluonScalarRulesMassiveChannel[ang_,mass_,momsK_List] := Function[expr,expr//.MassiveKHatRules[momsK]//.BCFWGluonScalarChainRules[ang,mass]//.BCFWGluonScalarSpinorRules[ang,mass]/.{Global`z->BCFWGluonScalarz[ang,mass,momsK/.{OverHat->Identity}]}];







(****************************************************)
(* Shift on massive and massless leg. For the       *)
(* massless momentum the square spinor is shifted    *)
(****************************************************)

BCFWScalarGluonSpinorRules[mass_,squ_]:={
				HoldPattern[SpinorAngleBracket[OverHat[squ],x_]] :> SpinorAngleBracket[squ,x],
				HoldPattern[SpinorAngleBracket[x_,OverHat[squ]]] :> SpinorAngleBracket[x,squ],
				
				HoldPattern[SpinorSquareBracket[OverHat[squ],x_]] :> SpinorSquareBracket[squ,x]+Global`z*Chain[$angle,squ,{UnderBar[mass]},x,$square],
				HoldPattern[SpinorSquareBracket[x_,OverHat[squ]]] :> SpinorSquareBracket[x,squ]-Global`z*Chain[$square,x,{UnderBar[mass]},squ,$angle],
								
				HoldPattern[Chain[$square,OverHat[squ],{x__},c__]] :> Chain[$square,squ,{x},c]+Global`z*Chain[$angle,squ,{UnderBar[mass],x},c],
				HoldPattern[Chain[$angle,OverHat[squ],{x__},c__]] :> Chain[$angle,squ,{x},c],
				HoldPattern[Chain[a__,{x__},OverHat[squ],$square]] :> Chain[a,{x},squ,$square]-Global`z*Chain[a,{x,UnderBar[mass]},squ,$angle],
				HoldPattern[Chain[a__,{x__},OverHat[squ],$angle]] :> Chain[a,{x},squ,$angle],
				
				Extramass[OverHat[x_]] :> Extramass[x],
				Extramasstilde[OverHat[x_]] :> Extramasstilde[x]
				
			};


BCFWScalarGluonChainRules[mass_,squ_]:={
	HoldPattern[Chain[$angle,a_,{b___,UnderBar[OverHat[mass]],c___},d__]] :> Chain[$angle,a,{b,UnderBar[mass],c},d]+Global`z*
		   If[OddQ@Length@{b},
		      Chain[$angle,a,{b,UnderBar[mass]},squ,$angle]Chain[$angle,squ,{c},d]
		    , 
		      -Chain[$angle,a,{b},squ,$angle]Chain[$angle,squ,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$square,a_,{b___,UnderBar[OverHat[mass]],c___},d__]] :> Chain[$square,a,{b,UnderBar[mass],c},d]+Global`z*
		   If[EvenQ@Length@{b},
		      Chain[$square,a,{b,UnderBar[mass]},squ,$angle]Chain[$angle,squ,{c},d]
		    , 
		      -Chain[$square,a,{b},squ,$angle]Chain[$angle,squ,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$angle,a_,{b___,OverHat[squ]|UnderBar[OverHat[squ]],c___},d__]] :> Chain[$angle,a,{b,UnderBar[squ],c},d]-Global`z*
		   If[OddQ@Length@{b},
		      Chain[$angle,a,{b,UnderBar[mass]},squ,$angle]Chain[$angle,squ,{c},d]
		    , 
		      -Chain[$angle,a,{b},squ,$angle]Chain[$angle,squ,{UnderBar[mass],c},d]
		   ],
	HoldPattern[Chain[$square,a_,{b___,OverHat[squ]|UnderBar[OverHat[squ]],c___},d__]] :> Chain[$square,a,{b,UnderBar[squ],c},d]-Global`z*
		   If[OddQ@Length@{b},
		      Chain[$square,a,{b,UnderBar[mass]},squ,$angle]Chain[$angle,squ,{c},d]
		    , 
		      -Chain[$square,a,{b},squ,$angle]Chain[$angle,squ,{UnderBar[mass],c},d]
		   ],
	S6Mom[i___,OverHat[ang],j___,OverHat[squ],k___]:>S6Mom[i,ang,j,squ,k],
	S6Mom[OverHat[mass],j__] :> S6Mom[mass,j]+Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S6Mom[OverHat[squ],j__] :> S6Mom[squ,j]-Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S6Mom[pm[OverHat[ang]],j__] :> S6Mom[pm[ang],j]-Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S6Mom[pm[OverHat[squ]],j__] :> S6Mom[pm[squ],j]+Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S4Mom[i___,OverHat[ang],j___,OverHat[squ],k___]:>S4Mom[i,ang,j,squ,k],
	S4Mom[OverHat[mass],j__] :> S4Mom[mass,j]+Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S4Mom[OverHat[squ],j__] :> S4Mom[squ,j]-Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S4Mom[pm[OverHat[ang]],j__] :> S4Mom[pm[ang],j]-Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	S4Mom[pm[OverHat[squ]],j__] :> S4Mom[pm[squ],j]+Global`z*Sum[Chain[$angle,squ,{UnderBar[k],UnderBar[mass]},squ,$angle],{k,{j}}],
	Chain[$angle,i_,{},j_,$angle]:>SpinorAngleBracket[i,j],
	Chain[$square,i_,{},j_,$square]:>SpinorSquareBracket[i,j]
	};


BCFWScalarGluonz[mass_,squ_,channel_List] /; MemberQ[channel,mass] := +S6Mom@@channel/Sum[Chain[$angle,squ,{UnderBar[mass],UnderBar[ki]},squ,$angle],{ki,Select[channel,# =!= mass &]}];
BCFWScalarGluonz[mass_,squ_,channel_List] /; MemberQ[channel,squ] := -S6Mom@@channel/Sum[Chain[$angle,squ,{UnderBar[mass],UnderBar[ki]},squ,$angle],{ki,Select[channel,# =!= squ &]}];


(* BCFWRulesMassless and BCFWRulesMassive take hatted momenta in the list momsK defining the factorization channel *)
BCFWScalarGluonRulesMasslessChannel[expr_,mass_,squ_,momsK_List] := BCFWScalarGluonRulesMasslessChannel[mass,squ,momsK][expr];
BCFWScalarGluonRulesMassiveChannel[expr_,mass_,squ_,momsK_List] := BCFWScalarGluonRulesMassiveChannel[mass,squ,momsK][expr];
BCFWScalarGluonRulesMasslessChannel[mass_,squ_,momsK_List] := Function[expr,expr//.MasslessKHatRules[Global`q,squ,momsK]//.BCFWScalarGluonChainRules[mass,squ]//.BCFWScalarGluonSpinorRules[mass,squ]/.{Global`z->BCFWScalarGluonz[mass,squ,momsK/.{OverHat->Identity}]}/.{UnderBar[squ]->squ}];
BCFWScalarGluonRulesMasslessChannel[mass_,squ_,momsK_List,{projSp1_,projSp2_}] := Function[expr,expr//.MasslessKHatRules[projSp1,projSp2,momsK]//.BCFWScalarGluonChainRules[mass,squ]//.BCFWScalarGluonSpinorRules[mass,squ]/.{Global`z->BCFWScalarGluonz[mass,squ,momsK/.{OverHat->Identity}]}/.{UnderBar[squ]->squ}];
BCFWScalarGluonRulesMassiveChannel[mass_,squ_,momsK_List] := Function[expr,expr//.MassiveKHatRules[momsK]//.BCFWScalarGluonChainRules[mass,squ]//.BCFWScalarGluonSpinorRules[mass,squ]/.{Global`z->BCFWScalarGluonz[mass,squ,momsK/.{OverHat->Identity}]}/.{UnderBar[squ]->squ}];





(* ::Section:: *)
(*New Routines*)


(* ::Subsection:: *)
(*Analytic BCFW Rules*)


(* Momentum always flows from left vertex to the rigth vertex *)
ApplyBcfwScalarGluonShiftRules[expr_,channel_List,shiftSMom_,shiftGMom_,masslessMomenta_List,K_]:=Module[{result},
	result=expr//.
	{
		Mom4DVector[k1_][Null][$up].Mom4DVector[k2_][Null][$down]:>1/2*S6Mom[k1,k2],
(*		S6Mom[a___,OverHat[K],b___,OverHat[shiftSMom],c___]:>S6Mom[a,b,c,Sequence@@channel,shiftSMom],*)
		Chain[$square,OverHat[K],{UnderBar@OverHat[shiftSMom],restmoms___},rest__]:>Chain[$square,OverHat[K],{PlusM[UnderBar[shiftSMom],Sequence@@UnderBar/@channel],restmoms},rest],
		Chain[most__,{mostmoms___,UnderBar@OverHat[shiftSMom]},OverHat[K],$square]:>Chain[most,{mostmoms,PlusM[UnderBar[shiftSMom],Sequence@@UnderBar/@channel]},OverHat[K],$square]
	};
	result = result //.
	{
		HoldPattern[SpinorSquareBracket[k_,OverHat[K]]]:>
			Chain[$square,k,{UnderBar@K,UnderBar@shiftSMom},shiftGMom,$square]/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square],
		HoldPattern[SpinorSquareBracket[OverHat[K],k_]]:>
			-Chain[$square,k,{UnderBar@K,UnderBar@shiftSMom},shiftGMom,$square]/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square],
		Chain[$square,OverHat[K],moms_List,rest___]:>1/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square]Chain[$square,shiftGMom,Join[{UnderBar[shiftSMom],UnderBar[K]},moms],rest],
		Chain[most__,moms_List,OverHat[K],$square]:>1/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square]Chain[most,Join[moms,{UnderBar[K],UnderBar[shiftSMom]}],shiftGMom,$square],
		HoldPattern[SpinorSquareBracket[k_,OverHat[K]]]:>Chain[$square,k,{UnderBar@K,UnderBar@shiftSMom},shiftGMom,$square]/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square],
		HoldPattern[SpinorSquareBracket[OverHat[K],k_]]:>-Chain[$square,k,{UnderBar@K,UnderBar@shiftSMom},shiftGMom,$square]/Chain[$angle,OverHat[K],{UnderBar[shiftSMom]},shiftGMom,$square]
	};
	result = result //.
	{
		Chain[$angle,OverHat[K],moms_List,rest___]:>1/SpinorSquareBracket[shiftGMom,OverHat[K]]Chain[$square,shiftGMom,Prepend[moms,UnderBar[K]],rest],
		Chain[most__,moms_List,OverHat[K],$angle]:>1/SpinorSquareBracket[OverHat[K],shiftGMom]Chain[most,Append[moms,UnderBar@K],shiftGMom,$square],
		
		HoldPattern[SpinorAngleBracket[k_,OverHat[K]]]:>Chain[$angle,k,{UnderBar@K},shiftGMom,$square]/SpinorSquareBracket[OverHat[K],shiftGMom],
		HoldPattern[SpinorAngleBracket[OverHat[K],k_]]:>-Chain[$angle,k,{UnderBar@K},shiftGMom,$square]/SpinorSquareBracket[OverHat[K],shiftGMom],
		
		Chain[$angle,OverHat@shiftGMom,moms_List,rest__]:>Chain[$angle,shiftGMom,moms,rest]-
			S6Mom@@channel*Chain[$square,shiftGMom,Prepend[moms,UnderBar@shiftSMom],rest]/
				Chain[$square,shiftGMom,{UnderBar@shiftSMom,PlusM@@UnderBar/@channel},shiftGMom,$square],
		Chain[most__,moms_List,OverHat@shiftGMom,$angle]:>Chain[most,moms,shiftGMom,$angle]-
			S6Mom@@channel*Chain[most,Append[moms,UnderBar@shiftSMom],shiftGMom,$square]/
				Chain[$square,shiftGMom,{PlusM@@UnderBar/@channel,UnderBar@shiftSMom},shiftGMom,$square],
				
		HoldPattern[SpinorAngleBracket[OverHat@shiftGMom,b_]]:>SpinorAngleBracket[shiftGMom,b]-S6Mom@@channel*
			Chain[$square,shiftGMom,{UnderBar@shiftSMom},b,$angle]/Chain[$square,shiftGMom,{UnderBar@shiftSMom,UnderBar/@PlusM@@channel},shiftGMom,$square],
		HoldPattern[SpinorAngleBracket[b_,OverHat@shiftGMom]]:>-(SpinorAngleBracket[shiftGMom,b]-S6Mom@@channel*
			Chain[$square,shiftGMom,{UnderBar@shiftSMom},b,$angle]/Chain[$square,shiftGMom,{UnderBar@shiftSMom,UnderBar/@PlusM@@channel},shiftGMom,$square])
	};
	result = result //.
	{
		HoldPattern[SpinorSquareBracket[a_,OverHat[shiftGMom]]]:>SpinorSquareBracket[a,shiftGMom],
		HoldPattern[SpinorSquareBracket[OverHat[shiftGMom],a_]]:>SpinorSquareBracket[shiftGMom,a],
		Chain[$square,OverHat@shiftGMom,rest__]:>Chain[$square,shiftGMom,rest],
		Chain[most__,OverHat@shiftGMom,$square]:>Chain[most,shiftGMom,$square]
	};
	result = result //.
	{
		S6Mom[rest___,OverHat[K],OverHat[shiftSMom]]:>S6Mom[rest,Sequence@@channel,shiftSMom],
		S6Mom[rest___,pm@OverHat[K],OverHat[shiftGMom]]:>S6Mom[rest,Sequence@@pm/@channel,shiftGMom]
	};
	result = result //.
	{
		S6Mom[rest__,OverHat[K]]:>S6Mom[rest,Sequence@@channel]-S6Mom@@channel*
			Chain[$square,shiftGMom,{UnderBar/@PlusM@@Join[{rest},channel],UnderBar@shiftSMom},shiftGMom,$square]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square],
		S6Mom[rest__,pm@OverHat[K]]:>S6Mom[Sequence@@pm/@{rest},OverHat[K]],
		S6Mom[rest__,OverHat[shiftGMom]]:>S6Mom[rest,shiftGMom]-S6Mom@@channel*
			Chain[$square,shiftGMom,{UnderBar/@PlusM[rest],UnderBar@shiftSMom},shiftGMom,$square]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square],
		S6Mom[rest__,OverHat[shiftSMom]]:>S6Mom[rest,shiftSMom]+S6Mom@@channel*
			Chain[$square,shiftGMom,{UnderBar/@PlusM[rest],UnderBar@shiftSMom},shiftGMom,$square]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square]
	};
	result = result //.
	{
		Chain[first__,{a___,PlusM[b__,OverHat[x_]],c___},rest__]:>Chain[first,{a,PlusM[b],c},rest]+Chain[first,{a,OverHat[x],c},rest],
		Chain[first__,{a___,PlusM[b__,UnderBar@OverHat[x_]],c___},rest__]:>Chain[first,{a,PlusM[b],c},rest]+Chain[first,{a,UnderBar@OverHat[x],c},rest],
		Chain[first__,{a___,PlusM[b__,pm@OverHat[x_]],c___},rest__]:>Chain[first,{a,PlusM[b],c},rest]-Chain[first,{a,OverHat[x],c},rest],
		Chain[first__,{a___,PlusM[b__,pm@UnderBar@OverHat[x_]],c___},rest__]:>Chain[first,{a,PlusM[b],c},rest]-Chain[first,{a,UnderBar@OverHat[x],c},rest],
(*		Chain[first__,{a___,OverHat[K],c___},rest__]:>Chain[first,{a,PlusM@@UnderBar/@channel,c},rest],*)
		UnderBar[K]:>PlusM@@UnderBar/@channel,Sequence@@(Rule[UnderBar@#,#]&/@masslessMomenta),PlusM[x_]:>x
	};
	result = result //.
	{
		Chain[type_,first_,{a___,OverHat[K]|UnderBar@OverHat[K],b___},rest__]/; (type===$square && EvenQ@Length@{a})|| (type===$angle && OddQ@Length@{a}) :>
			Chain[type,first,{a,PlusM@@UnderBar/@channel,b},rest]-S6Mom@@channel Chain[type,first,{a},shiftGMom,$square]*
				Chain[$square,shiftGMom,{UnderBar@shiftSMom,b},rest]/Chain[$square,shiftGMom,{UnderBar@shiftSMom,UnderBar/@PlusM@@channel},shiftGMom,$square],
		Chain[type_,first_,{a___,OverHat[K]|UnderBar@OverHat[K],b___},rest__]/; (type===$square && OddQ@Length@{a})|| (type===$angle && EvenQ@Length@{a}) :>
			Chain[type,first,{a,PlusM@@UnderBar/@channel,b},rest]-S6Mom@@channel Chain[type,first,{a,UnderBar@shiftSMom},shiftGMom,$square]*
				Chain[$square,shiftGMom,{b},rest]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square],
	
		Chain[type_,first_,{a___,OverHat[shiftGMom],b___},rest__]/; (type===$square && EvenQ@Length@{a})|| (type===$angle && OddQ@Length@{a}) :>
			Chain[type,first,{a,shiftGMom,b},rest]-S6Mom@@channel Chain[type,first,{a},shiftGMom,$square]*
				Chain[$square,shiftGMom,{UnderBar@shiftSMom,b},rest]/Chain[$square,shiftGMom,{UnderBar@shiftSMom,UnderBar/@PlusM@@channel},shiftGMom,$square],
		Chain[type_,first_,{a___,OverHat[shiftGMom],b___},rest__]/; (type===$square && OddQ@Length@{a})|| (type===$angle && EvenQ@Length@{a}) :>
			Chain[type,first,{a,shiftGMom,b},rest]-S6Mom@@channel Chain[type,first,{a,UnderBar@shiftSMom},shiftGMom,$square]*
				Chain[$square,shiftGMom,{b},rest]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square],
				
		Chain[type_,first_,{a___,UnderBar@OverHat[shiftSMom],b___},rest__]/; (type===$square && EvenQ@Length@{a})|| (type===$angle && OddQ@Length@{a}) :>
			Chain[type,first,{a,UnderBar@shiftSMom,b},rest]+S6Mom@@channel Chain[type,first,{a},shiftGMom,$square]*
				Chain[$square,shiftGMom,{UnderBar@shiftSMom,b},rest]/Chain[$square,shiftGMom,{UnderBar@shiftSMom,UnderBar/@PlusM@@channel},shiftGMom,$square],
		Chain[type_,first_,{a___,UnderBar@OverHat[shiftSMom],b___},rest__]/; (type===$square && OddQ@Length@{a})|| (type===$angle && EvenQ@Length@{a}) :>
			Chain[type,first,{a,UnderBar@shiftSMom,b},rest]+S6Mom@@channel Chain[type,first,{a,UnderBar@shiftSMom},shiftGMom,$square]*
				Chain[$square,shiftGMom,{b},rest]/Chain[$square,shiftGMom,{UnderBar/@PlusM@@channel,UnderBar@shiftSMom},shiftGMom,$square]
	};
	result = result //.
	{
		UnderBar[K]:>PlusM@@UnderBar/@channel,Sequence@@(Rule[UnderBar@#,#]&/@masslessMomenta),PlusM[x_]:>x,
		Chain[a_,k:Alternatives@@masslessMomenta,{PlusM[b___,k:Alternatives@@masslessMomenta,c___],d___},e___]:> Chain[a,k,{PlusM[b,c],d},e],
		Chain[a__,{b___,PlusM[c___,k:Alternatives@@masslessMomenta,d___]},k:Alternatives@@masslessMomenta,e_]:>Chain[a,{b,PlusM[c,d]},k,e],
		Chain[a_,b_,{PlusM[c___,d_,e___],d_},b_,a_]:>Chain[a,b,{PlusM[c,e],d},b,a],
		Chain[a_,b_,{d_,PlusM[c___,d_,e___]},b_,a_]:>Chain[a,b,{d,PlusM[c,e]},b,a]
	};
	result = result //.
	{
		PlusM[x_]:>x,Chain[first__,{a___,b_,b_,c___},rest__]/; Head[b]=!=PlusM :> (Extramass[b]Extramasstilde[b]//.{UnderBar->Identity,OverHat->Identity})Chain[first,{a,c},rest],
		UnderBar[K]:>PlusM@@UnderBar/@channel,Sequence@@(Rule[UnderBar@#,#]&/@masslessMomenta),PlusM[x_]:>x,
		Extramass[OverHat[x_]]:>Extramass[x],Extramasstilde[OverHat[x_]]:>Extramasstilde[x],
		SpinorSquareBracket[shiftGMom,OverHat[K]]->1
	}(*;
	result //.{
		Chain[a__,{b___,PlusM[c__],PlusM[c__,d__],e___},f__]\[RuleDelayed](S6Mom[c]/.UnderBar\[Rule]Identity)Chain[a,{b,e},f]+Chain[a,{b,PlusM[c],d,e},f],
		Chain[a__,{b___,c:Except[_PlusM],c:Except[_PlusM],e___},f__]\[RuleDelayed](Extramass[c]Extramasstilde[c]/.UnderBar\[Rule]Identity)Chain[a,{b,e},f],
		Chain[a__,{b___,PlusM[c_,d__]},c_,f_]\[RuleDelayed]Chain[a,{b,PlusM[d]},c,f]
		(* These rules break stuff, though I do not know why... *)
	}*)
]


(* ::Subsection::Closed:: *)
(*Simplification Functions*)


SortChains[expr_]:=expr//.
	{
		Chain[$square,a_,moms_List,b_,$angle]:>Chain[$angle,b,Reverse@moms,a,$square]
	}//.
	{
		Chain[a_,b_,moms_List,c_,a_] /; b=!=c && Order[b,c]===-1:>-Chain[a,c,Reverse@moms,b,a],
		Chain[a_,b_,moms_List,b_,a_] /; Order[moms,Reverse@moms]===-1 :>-Chain[a,b,Reverse@moms,b,a]
	};
	
AlignExtramassRules[{p1_,p2_},{p3_,p4_}]:={Extramass[p2]->-Extramass[p1],Extramass[p4]->-Extramass[p3],Extramasstilde[p2]->-Extramasstilde[p1],Extramasstilde[p4]->-Extramasstilde[p3]}


ReduceMomenta[allMomenta_List,masslessMomenta_List]:=Function[expr,
	Module[{n=Length@allMomenta},
		expr//.
		{
			S6Mom[moms__] /; ContainsAll[allMomenta,{moms}] && Length[{moms}]>n/2 :> S6Mom@@Complement[allMomenta,{moms}],
			Chain[first__,{a___,PlusM[moms__],b___},last__] /; ContainsAll[allMomenta,{moms}/.UnderBar->Identity] && Length[{moms}]>n/2 :> 
				-Chain[first,{a,PlusM@@UnderBar/@Complement[allMomenta,{moms}/.UnderBar->Identity],b},last],
			Sequence@@(Rule[UnderBar@#,#]&/@masslessMomenta)
		}
	]
]


(* ::Chapter:: *)
(*Postamble*)


(* definitions for system functions *)

Protect[Evaluate[protected]];

End[];


Protect[ "SpinorHelicity6D`BCFWAnalytics`*" ];    (* protect exported symbols *)
Unprotect[qBCFW];
EndPackage[ ];  (* end the package context *)
