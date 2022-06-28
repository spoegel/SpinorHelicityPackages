(* ::Package:: *)

(* ::Title:: *)
(*S@M Scalar Kinematics Example*)


(* ::Chapter:: *)
(*Setup*)


SetDirectory[NotebookDirectory[]];


<<Spinors`


<<"setupScalarKinematics.m"


<<"ATreeSAM.m";


(* ::Chapter:: *)
(*Two Scalars*)


twoScalarFivePointKinematics = <<"twoScalar5PtPoints.m";


SetupSAMTwoScalarMomentumKinematics[{p1,p5},{p2,p3,p4},twoScalarFivePointKinematics[[1]],Quiet->True]


ATreeSAM["S","+","+","-","S"][p1,p2,p3,p4,p5]
%//N
%%//NRational
%//N


(* ::Chapter:: *)
(*Four Scalars*)


fourScalarSevenPointKinematics = <<"fourScalar7PtPoints.m";


SetupSAMFourScalarMomentumKinematics[{p1,p2},{p6,p7},{p3,p4,p5},fourScalarSevenPointKinematics[[1]],Quiet->True]


ATreeSAM["S","S","+","+","+","s","s"][p1,p2,p3,p4,p5,p6,p7];
%//N
%%//NRational
%//N


(* ::Chapter:: *)
(*S6Mom to s/MP*)


SetupSAMFourScalarMomentumKinematics[{p1,p2},{p6,p7},{p3,p4,p5},fourScalarSevenPointKinematics[[1]],Quiet->True]


(* ::Text:: *)
(*These are*)


(*S6Mom[p1,p2]*) 
s[p1,p2]
2MP[p1,p1]+2MP[p1,p2]

(*S6Mom[p1,p3]*) 
2MP[p1,p3]

(*S6Mom[p1,p2,p3]*) 
s[p1,p2,p3]
2MP[p1,p1]+2MP[p1,p2]+2MP[p1,p3]+2MP[p2,p3]

(*S6Mom[p3,p4]*) 
s[p3,p4]
2MP[p3,p4]

(*S6Mom[p1,p3,p4]*) 
s[p1,p3,p4]-MP[p1,p1]
2MP[p1,p3]+2MP[p1,p4]+2MP[p3,p4]

(*S6Mom[p1,p7]: Here we are suppressing \[Mu]_{17} terms, as these should not be present in four-dimensional scalar amplitudes, and also end up being irrelevant to rational contributions *) 
2MP[p1,p7]

(*S6Mom[p1,p2,p6]: Again supressing \[Mu]_{16},\[Mu]_{26} cross-terms *) 
s[p1,p2,p6]-MP[p6,p6]
2MP[p1,p1]+2MP[p1,p2]+2MP[p1,p6]+2MP[p2,p6]


Subsets[{p1,p2,p3,p4},{2}]


(* ::Chapter:: *)
(*Test*)


SetupSAMFourScalarMomentumKinematics[{p1,p2},{p6,p7},{p3,p4,p5},fourScalarSevenPointKinematics[[1]],Quiet->True]
test1 = s[p1,p2] //NRational


SetupSAMFourScalarMomentumKinematics[{p1,p2},{p6,p7},{p3,p4,p5},fourScalarSevenPointKinematics[[2]],Quiet->True]
test2 = s[p1,p2] //NRational


tests = Array[(
	SetupSAMFourScalarMomentumKinematics[{p1,p2},{p6,p7},{p3,p4,p5},fourScalarSevenPointKinematics[[#]],Quiet->True];
	s[p1,p2]//NRational)&,
	2
]
