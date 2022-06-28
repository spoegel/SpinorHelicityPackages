(* ::Package:: *)

(* ::Title:: *)
(*Berends-Giele*)


BeginPackage["SpinorHelicity6D`BerendsGiele`",{"SpinorHelicity6D`"}];

Unprotect["SpinorHelicity6D`BerendsGiele`*"];
ClearAll["SpinorHelicity6D`BerendsGiele`*"];
ClearAll["SpinorHelicity6D`BerendsGiele`*"];

Jg::usage = "";
Js::usage = "";
JS::usage = "";

ATreeBerendsGiele::usage = "";
ATreeBerendsGieleGluon::usage = "";
ATreeBerendsGieleContact::usage = "";

$BerendsGieleDefaultRef::usage = "Default reference momentum used for polarization vectors in Berends-Giele recursion.";
KillMasses[{$BerendsGieleDefaultRef}];
GenSpinors[
	{$BerendsGieleDefaultRef},
	FourD->{$BerendsGieleDefaultRef},
	RandomSpinors->True,
	Seed->12349867,
	ParameterRange->50,
	DisplaySpinors->False
];


$BerendsGieleRefs::usage = "Association that allows to defines specific reference momenta to be used for polarization vectors. To assign reference momentum x to gluon with label y, append y->x to this association. $BerendsGieleRefs also contains a Key \"Default\" which specifies the default reference momentum to be used otherwise.";
$BerendsGieleRefs = <|"Default"->$BerendsGieleDefaultRef|>;

$BerendsGieleSignGGG::usage = "Specifies sign of cubic gluon vertex.";
$BerendsGieleSignP::usage = "Specifies sign of positive helicity gluon polarization vectors.";
$BerendsGieleSignM::usage = "Specifies sign of negative helicity gluon polarization vectors."

$BerendsGieleSignP = -1;
$BerendsGieleSignM = +1;
$BerendsGieleSignGGG = 1;

$BerendsGieleSimplificationFunction = Identity; (* Used to be Together, maybe for a good reason? *)

$BerendsGieleSSSSVertexFactor = +1;

Begin["`Private`"];

EtaBerendsGiele = DiagonalMatrix[Metric];

UpdateEtaBerendsGiele[] := EtaBerendsGiele = DiagonalMatrix[Metric];

JS[args___] /; Head[JSCache[args]] =!= JSCache := JSCache[args];

JS[{{"S",_}}] :=1;
JS[{{_,_}}] := 0;
JS[helicityMomentumArgs_List] /; Length@helicityMomentumArgs >1 := JSCache[helicityMomentumArgs] =
Module[{init = ("init" === helicityMomentumArgs[[1]]),helicitiesMomenta = helicityMomentumArgs/.{"init"->Nothing},
	propagator,helicities,momenta,n,Sg,gS,Sgg,ggS,gSg,ssS,Sss,sSs},
	n=Length@helicitiesMomenta;
	{helicities,momenta} = Transpose@helicitiesMomenta;
	propagator = Identity@If[init,1,I/S6MomN@@momenta];

	gS = $BerendsGieleSignGGG*Sum[Identity[I/Sqrt[2]*JS[helicitiesMomenta[[i+1;;]]]*
			Dot[-Sum[Mom4DVectorN[k2][Null][$down],{k2,momenta[[1;;i]]}] - 2*Sum[Mom4DVectorN[k3][Null][$down],{k3,momenta[[i+1;;]]}],Jg[helicitiesMomenta[[1;;i]]]]],
		{i,1,n-1}];

	Sg = $BerendsGieleSignGGG*Sum[Identity[I/Sqrt[2]*JS[helicitiesMomenta[[1;;i]]]*
			Dot[2*Sum[Mom4DVectorN[k2][Null][$down],{k2,momenta[[1;;i]]}] + Sum[Mom4DVectorN[k3][Null][$down],{k3,momenta[[i+1;;]]}],Jg[helicitiesMomenta[[i+1;;]]]]],
		{i,1,n-1}];

	ggS = Sum[Identity[I/2 JS[helicitiesMomenta[[i2+1;;]]]*Dot[Jg[helicitiesMomenta[[1;;i1]]],EtaBerendsGiele,
		Jg[helicitiesMomenta[[i1+1;;i2]]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	Sgg = Sum[Identity[I/2 JS[helicitiesMomenta[[1;;i1]]]*Dot[Jg[helicitiesMomenta[[i1+1;;i2]]],EtaBerendsGiele,Jg[helicitiesMomenta[[i2+1;;]]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	gSg = Sum[Identity[-I JS[helicitiesMomenta[[i1+1;;i2]]]*Dot[Jg[helicitiesMomenta[[1;;i1]]],
		EtaBerendsGiele,Jg[helicitiesMomenta[[i2+1;;]]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	ssS = Sum[Identity[-I/2 *Js[helicitiesMomenta[[1;;i1]]]*Js[helicitiesMomenta[[i1+1;;i2]]]*JS[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	Sss = Sum[Identity[-I/2 *JS[helicitiesMomenta[[1;;i1]]]*Js[helicitiesMomenta[[i1+1;;i2]]]*Js[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	sSs = Sum[Identity[I *Js[helicitiesMomenta[[1;;i1]]]*JS[helicitiesMomenta[[i1+1;;
			i2]]]*Js[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	Identity@propagator*Plus@@{gS,Sg,ggS,Sgg,gSg,$BerendsGieleSSSSVertexFactor*
     	ssS,$BerendsGieleSSSSVertexFactor*Sss,$BerendsGieleSSSSVertexFactor*sSs}
];

Js[args___] /; Head[JsCache[args]] =!= JsCache := JsCache[args];

(*Js[{{"s",_}}] := Series[x*t1*Mu1*MuTilde1*t2*Mu2*MuTilde2,{t1,Infinity,5},{Mu1,Infinity,5},{MuTilde1,Infinity,5},{t2,Infinity,5},{Mu2,Infinity,5},{MuTilde2,Infinity,5}];*)
Js[{{"s",_}}] := 1;
Js[{{_,_}}] := 0;
Js[helicityMomentumArgs_List] /; Length@helicityMomentumArgs >1 := JsCache[helicityMomentumArgs] =
Module[{init = ("init" === helicityMomentumArgs[[1]]),
	helicitiesMomenta = helicityMomentumArgs/.{"init"->Nothing},propagator,helicities,momenta,n,sg,
	gs,sgg,ggs,gsg,SSs,sSS,SsS},

	n=Length@helicitiesMomenta;
	{helicities,momenta} = Transpose@helicitiesMomenta;
	propagator = Identity@If[init,1,I/S6MomN@@momenta];

	gs = $BerendsGieleSignGGG*Sum[Identity[I/Sqrt[2]*Js[helicitiesMomenta[[i+1;;]]]*
			Dot[-Sum[Mom4DVectorN[k2][Null][$down],{k2,momenta[[1;;i]]}] - 2*Sum[Mom4DVectorN[k3][Null][$down],{k3,momenta[[i+1;;]]}],Jg[helicitiesMomenta[[1;;i]]]]],
		{i,1,n-1}];

	sg = $BerendsGieleSignGGG*Sum[Identity[I/Sqrt[2]*Js[helicitiesMomenta[[1;;i]]]*
			Dot[2*Sum[Mom4DVectorN[k2][Null][$down],{k2,momenta[[1;;i]]}] + Sum[Mom4DVectorN[k3][Null][$down],{k3,momenta[[i+1;;]]}],Jg[helicitiesMomenta[[i+1;;]]]]],
		{i,1,n-1}];

	ggs = Sum[Identity[I/2 Js[helicitiesMomenta[[i2+1;;]]]*Dot[Jg[helicitiesMomenta[[1;;i1]]],EtaBerendsGiele,
		Jg[helicitiesMomenta[[i1+1;;i2]]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	sgg = Sum[Identity[I/2 Js[helicitiesMomenta[[1;;i1]]]*Dot[Jg[helicitiesMomenta[[i1+1;;i2]]],EtaBerendsGiele,Jg[helicitiesMomenta[[i2+1;;]]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	gsg = Sum[Identity[-I Js[helicitiesMomenta[[i1+1;;i2]]]*Dot[Jg[helicitiesMomenta[[1;;i1]]],
		EtaBerendsGiele,Jg[helicitiesMomenta[[i2+1;;]]]]], {i1,1,n-2},{i2,i1+1,n-1}];

	SSs = Sum[Identity[-I/2 *JS[helicitiesMomenta[[1;;i1]]]*JS[helicitiesMomenta[[i1+1;;i2]]]*Js[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	sSS = Sum[Identity[-I/2 *Js[helicitiesMomenta[[1;;i1]]]*JS[helicitiesMomenta[[i1+1;;i2]]]*JS[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	SsS = Sum[Identity[I *JS[helicitiesMomenta[[1;;i1]]]*Js[helicitiesMomenta[[i1+1;;
			i2]]]*JS[helicitiesMomenta[[i2+1;;]]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	Identity@propagator*Plus@@{gs,sg,ggs,sgg,gsg,$BerendsGieleSSSSVertexFactor* sSS,
		$BerendsGieleSSSSVertexFactor*SSs,$BerendsGieleSSSSVertexFactor*SsS}
];


(* ::Text:: *)
(*Jg is always assumed to have index $up*)

Jg[args___] /; Head[JgCache[args]] =!= JgCache := JgCache[args];

Jg[{{"+",p_}}] /; KeyExistsQ[$BerendsGieleRefs,p] := JgCache[helicityMomentumArgs] =
	$BerendsGieleSignP /(Sqrt[2]*SpinorAngleBracketN[p,$BerendsGieleRefs[p]])*
			(ChainN[$square,p,{PauliMatrixSpinor[#, $up]},$BerendsGieleRefs[p],$angle]&/@ {0,1,2,3});

Jg[{{"+",p_}}] := JgCache[helicityMomentumArgs] =
	$BerendsGieleSignP /(Sqrt[2]*SpinorAngleBracketN[p,$BerendsGieleRefs["Default"]])*
     (ChainN[$square,p,{PauliMatrixSpinor[#, $up]},$BerendsGieleRefs["Default"],$angle]&/@ {0,1,2,3});

Jg[{{"-",p_}}] /; KeyExistsQ[$BerendsGieleRefs,p] := JgCache[helicityMomentumArgs] =
	$BerendsGieleSignM /(Sqrt[2]*SpinorSquareBracketN[p,$BerendsGieleRefs[p]])*
			(ChainN[$angle,p,{PauliMatrixSpinor[#, $up]},$BerendsGieleRefs[p],$square]&/@{0,1,2,3});

Jg[{{"-",p_}}] := JgCache[helicityMomentumArgs] =
	$BerendsGieleSignM /(Sqrt[2]*SpinorSquareBracketN[p,$BerendsGieleRefs["Default"]])*
     (ChainN[$angle,p,{PauliMatrixSpinor[#, $up]},$BerendsGieleRefs["Default"],$square]&/@{0,1,2,3});

Jg[{{"0",p__}}] := JgCache[helicityMomentumArgs] =
    Sqrt[S6MomN[p]]/(Sum[Mom4DVectorN[i][Null][$down],{i,{p}}] .
		Mom4DVectorN[$BerendsGieleRefs["Default"]][Null][$up])Mom4DVectorN[$BerendsGieleRefs["Default"]][Null][$up];

Jg[{{_,_}}] := {0,0,0,0};

Jg[helicityMomentumArgs_List] /; Length@helicityMomentumArgs >1 := JgCache[helicityMomentumArgs] =
Module[{init = ("init" === helicityMomentumArgs[[1]]),helicitiesMomenta = helicityMomentumArgs/.{"init"->Nothing},
	propagator,helicities,momenta,n,gg,ggg,ssg,gss,ss,SSg,gSS,SS,SgS,sgs},
	n=Length@helicitiesMomenta;
	{helicities,momenta} = Transpose@helicitiesMomenta;
	propagator = $BerendsGieleSimplificationFunction@If[init,1,-I/S6MomN@@momenta];

	ss = $BerendsGieleSignGGG*Sum[I/Sqrt[2]*Js[helicitiesMomenta[[1;;i]]]*Js[helicitiesMomenta[[i+1;;]]]*
			$BerendsGieleSimplificationFunction[(Sum[Mom4DVectorN[k3][Null][$up],{k3,momenta[[i+1;;]]}] -
       Sum[Mom4DVectorN[k2][Null][$up],{k2,momenta[[1;;i]]}])],
		{i,1,n-1}];

	SS = $BerendsGieleSignGGG*Sum[I/Sqrt[2]*JS[helicitiesMomenta[[1;;i]]]*JS[helicitiesMomenta[[i+1;;]]]*
			$BerendsGieleSimplificationFunction[(Sum[Mom4DVectorN[k3][Null][$up],{k3,momenta[[i+1;;]]}] - Sum[Mom4DVectorN[k2][Null][$up],{k2,momenta[[1;;i]]}])],
		{i,1,n-1}];

	gg = $BerendsGieleSignGGG*Sum[VJggg[
		$BerendsGieleSimplificationFunction[Sum[Mom4DVectorN[k2][Null][$up],{k2,momenta[[1;;i]]}]],Jg[helicitiesMomenta[[1;;i]]],
		$BerendsGieleSimplificationFunction[Sum[Mom4DVectorN[k3][Null][$up],{k3,momenta[[i+1;;]]}]],Jg[helicitiesMomenta[[i+1;;]]]
		],
		{i,1,n-1}];
	
	gss = Sum[I/2 Jg[helicitiesMomenta[[1;;i1]]]*Js[helicitiesMomenta[[i1+1;;i2]]]*Js[helicitiesMomenta[[i2+1;;]]],
		{i1,1,n-2},{i2,i1+1,n-1}];
		
	ssg = Sum[I/2 Jg[helicitiesMomenta[[i2+1;;]]]*Js[helicitiesMomenta[[1;;i1]]]*Js[helicitiesMomenta[[i1+1;;i2]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	sgs = Sum[-I Js[helicitiesMomenta[[i2+1;;]]]*Js[helicitiesMomenta[[1;;i1]]]*Jg[helicitiesMomenta[[i1+1;;i2]]],
		{i1,1,n-2},{i2,i1+1,n-1}];
		
	gSS = Sum[I/2 Jg[helicitiesMomenta[[1;;i1]]]*JS[helicitiesMomenta[[i1+1;;i2]]]*JS[helicitiesMomenta[[i2+1;;]]],
		{i1,1,n-2},{i2,i1+1,n-1}];
		
	SSg = Sum[I/2 Jg[helicitiesMomenta[[i2+1;;]]]*JS[helicitiesMomenta[[1;;i1]]]*JS[helicitiesMomenta[[i1+1;;i2]]],
		{i1,1,n-2},{i2,i1+1,n-1}];

	SgS = Sum[-I JS[helicitiesMomenta[[i2+1;;]]]*JS[helicitiesMomenta[[1;;i1]]]*Jg[helicitiesMomenta[[i1+1;;i2]]],
		{i1,1,n-2},{i2,i1+1,n-1}];
		
	ggg = Sum[VJgggg[
			Jg[helicitiesMomenta[[1;;i1]]],
			Jg[helicitiesMomenta[[i1+1;;i2]]],
			Jg[helicitiesMomenta[[i2+1;;]]]
		],
		{i1,1,n-2},{i2,i1+1,n-1}];


	Identity[propagator]*Plus@@{gg,ggg,ssg,gss,ss,SSg,gSS,SS,SgS,sgs}
];


(* ::Text:: *)
(*VJggg has Lorentz index $up*)


VJggg[p2_,Jg2_,p3_,Jg3_]:= I/Sqrt[2](
	Dot[Jg2,EtaBerendsGiele,Jg3]*(p2-p3) +
	Jg2*Dot[Jg3,EtaBerendsGiele,-p3-2*p2]+
	Jg3*Dot[Jg2,EtaBerendsGiele,2*p3+p2]);


VJgggg[Jg2_,Jg3_,Jg4_]:= 
	I*Jg3*Dot[Jg2,EtaBerendsGiele,Jg4]-
	I/2(
		Jg2*Dot[Jg3,EtaBerendsGiele,Jg4]+
		Jg4*Dot[Jg2,EtaBerendsGiele,Jg3]
	);



ATreeBerendsGiele["S",particles__][_,momenta__]:=
		(
			UpdateEtaBerendsGiele[];
			ClearBerendsGieleCurrentCache[];
			(-I) * JS[{"init", Sequence @@ Transpose@{{particles}, {momenta}}}]
		);
ATreeBerendsGiele["s",particles__][_,momenta__]:=
		(
			UpdateEtaBerendsGiele[];
			ClearBerendsGieleCurrentCache[];
			(-I)*Js[{"init",Sequence@@Transpose@{{particles},{momenta}}}]
		);
ATreeBerendsGiele["-",particles__][momentum1_,momenta__]:=
		(
			UpdateEtaBerendsGiele[];
			ClearBerendsGieleCurrentCache[];
			(-I)*Dot[
				Jg[{{"-",momentum1}}],
				DiagonalMatrix[Metric],
				Jg[{"init",Sequence@@Transpose@{{particles},{momenta}}}]
			]
		);
ATreeBerendsGiele["+",particles__][momentum1_,momenta__]:=
		(
			UpdateEtaBerendsGiele[];
			ClearBerendsGieleCurrentCache[];
			(-I)*Dot[
				Jg[{{"+",momentum1}}],
				DiagonalMatrix[Metric],
				Jg[{"init",Sequence@@Transpose@{{particles}, {momenta}}}]
			]
		);

ATreeBerendsGieleGluon[args1___][args2___]:=
Module[{saveContact = $BerendsGieleSSSSVertexFactor, res},
	$BerendsGieleSSSSVertexFactor = 0;
	res = ATreeBerendsGiele[args1][args2];
	$BerendsGieleSSSSVertexFactor = saveContact;
	res
];

ATreeBerendsGieleContact[args1___][args2___]:=
Module[{saveContact = $BerendsGieleSSSSVertexFactor, resFull, resGluon},
	resFull = ATreeBerendsGiele[args1][args2];
	$BerendsGieleSSSSVertexFactor = 0;
	resGluon = ATreeBerendsGiele[args1][args2];
	$BerendsGieleSSSSVertexFactor = saveContact;
	resFull - resGluon
];

ClearBerendsGieleCurrentCache[] := ClearAll[JgCache,JsCache,JSCache];

(* ::Chapter:: *)
(*Postamble*)


End[];
Protect["SpinorHelicity6D`BerendsGiele`*"];
Unprotect[
	$BerendsGieleSignGGG,
	$BerendsGieleSignM,
	$BerendsGieleSignP,
	$BerendsGieleSSSSVertexFactor,
	$BerendsGieleRefs,
	$BerendsGieleSimplificationFunction
];
EndPackage[];
