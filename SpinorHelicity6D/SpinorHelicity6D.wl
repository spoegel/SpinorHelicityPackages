(*

    This file is part of SpinorHelicity.

    SpinorHelicity is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SpinorHelicityAddons is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SpinorHelicityAddons.  If not, see <http://www.gnu.org/licenses/>.

*)

(* ::Package:: *)

(* ::Title:: *)
(*SpinorHelicity6D*)


BeginPackage["SpinorHelicity6D`","Utils`"];


Unprotect["SpinorHelicity6D`*"];
Unprotect["SpinorHelicity6D`Private`*"];
ClearAll["SpinorHelicity6D`*"];
ClearAll["SpinorHelicity6D`Private`*"];


(* ::Section::Closed:: *)
(*Messages*)


$Shortcuts::usage="This is a list comprising all the available shortcuts and their output"


MinusSignQ::usage="Tells you if its argument is negative"


overbar::usage="Auxiliary function for OverBar."


SpinorDot::usage="This is the generic tilde 4d spinor. It serves both for the \[Lambda] and the \[Mu], it takes four arguments SpinorDot[momlabel][type][upper][lower].
The first labels the momentum, the second is either $lam or $mu and distinguishes between \[Lambda] and \[Mu], the third is for the upper index and the fourth is for the lower index. If no index
of a certain kind is present the argument to be given is Null. Do not worry, you are not supposed to use this functions. The shortcuts directly give you what you are seeking."


SpinorUndot::usage="This is the generic 4d spinor. It serves both for the \[Lambda] and the \[Mu], it takes four arguents SpinorUndot[momlabel][type][upper][lower].
The first labels the momentum, the second is either $lam or $mu and distinguishes between \[Lambda] and \[Mu], the firsd is for the upper index and the fourth is for the lower index. If no index
of a certain kind is present the argument to be given is Null. Do not worry, you are not supposed to use this functions. The shortcuts directly give you what you are seeking."


Spinorundot::usage="Auxiliary function for SpinorUndot"


$lam::usage="Protected label for the \[Lambda] spinors"


$mu::usage="Protected label for the \[Mu] spinors"


Spinordot::usage="Auxiliary function for SpinorDot"


SpinorDot6D::usage="SpinorDot6D[momlabel][Loretzindex][Null][lgindex] is the 6D spinor with momentum momlabel, Lorentz index Lorentzindex and the downstairs little group index lgindex"


SpinorUndot6D::usage="SpinorUndot6D[momlabel][Loretzindex][Null][lgindex] is the 6D tilde spinor with momentum momlabel, Lorentz index Lorentzindex and the downstairs little group index lgindex"


SpinorUndotPure::usage="The equivaqlent of SpinorUndot but without index structure"
SpinorDotPure::usage="The equivalent of SpinorDot but without index structure"


extramass::usage="Mass given by the extra dimensional part of the momenta."
extramasstilde::usage="auxiliary function encorporating all the properties of the mass \[CapitalMu] tilde"
Extramass::usage="auxiliary function for extramass"
Extramasstilde::usage="auxiliary function for extramasstilde"


KillMasses::usage="KillMasses[{a,b,c,...}] sets \[CapitalMu][a],\[CapitalMu][b],\[CapitalMu][c] and so on to zero (as well as \[CapitalMu] tilde)"
Momenta4D::usage="This a variable which tells you which momenta have been specified to be 4 dimensional. It takes no argument of course."
Set4DMassive::usage="Set4DMassive[{a,b,c,...}] sets \[CapitalMu][a],\[CapitalMu][b],\[CapitalMu][c] and so on to zero (as well as \[CapitalMu] tilde). In addition, some functions will treat these momenta differently
in simplifications and numerical evaluations, as they do not allow decomposition into 4D spinors."
Momenta4DMassive::usage="This a variable which tells you which momenta have been specified to be 4-dimensional massive."
NewProcess::usage="NewProcess deletes all the downvalues of any function fixed during the computation of a given process. For example masses of particles that were previously set to zero (being 4D) are restore using this function. Also conservation relations stored in AllIdentities6D are cleared."
ClearDownValues::usage="ClearDownValues[f] clears the downvalues of the function f given as an argument"


SquAngInvariant::usage="SquAngInvariant[n,m][adot,b] is the invariant [n_adot m_b > with adot and b little group indices and n,m labels of the particles. The associated shortcut is esc + sai + esc"
AngSquInvariant::usage="AngSquInvariant[n,m][a,bdot] is the invariant < n_a m_bdot ] with a and bdot littlegroup indices and n,m particle labels. The associated shortcut is esc + asi + esc"
AngAngInvariant::usage="AngAngInvariant[m,n,l,k][a,b,c,d] is the invariant < m_a n_b l_c k_d > where m,n,l,k are particle labels and a,b,c,d are little group indices. The shotcut is esc + aa + esc"
SquSquInvariant::usage="BraBrainvariant[m,n,l,k][a,b,c,d] is the invariant [ m_a n_b l_c k_d ] where m,n,l,k are particle labels and a,b,c,d are the little group indices"


To4D::usage="When applied to a function of the six-dimensional invariants, To4D breaks them down into 4 dimensional invariants, provided the little group indices assume values 1 or 2."
To4DAlways::usage="To4DAlways[True] will change global settings so that six-dimensional invariants are always automatically broken down into four-dimensional invariants. To4DAlways[False] undoes this global setting."


SpinorAngleBracket::usage = "SpinorAngleBracket[\*StyleBox[\"x\", \"TI\"], \*StyleBox[\"y\", \"TI\"]] represents the contraction of undotted spinors \!\(\*SubscriptBox[\(\[Lambda]\), \(x\)]\) and \!\(\*SubscriptBox[\(\[Lambda]\), \(y\)]\)."
SpinorSquareBracket::usage = "SpinorAngleBracket[\*StyleBox[\"x\", \"TI\"], \*StyleBox[\"y\", \"TI\"]] represents the contraction of undotted spinors" (*SubscriptBox[OverscriptBox[\"\[Lambda]\", \"~\"], #]*)
Spinoranglebracket::usage="Auxiliary function for SpinorAngleBracket"
Spinorsquarebracket::usage="Auxiliary function for SpinorSquareBracket"


levicivita2up::usage="levicivita2[a,b] is the element (a,b) (with a,b integers) of the two dimensional Levi-Civita tensor with upper indices and following the conventions in 1010.0494"
levicivita2down::usage="levicivita2[a,b] is the element (a,b) (with a,b integers) of the two dimensional Levi-Civita tensor with lower indices and following the conventions in 1010.0494"
levicivita2Up::usage="This is an auxiliary function needed in order to define the good contractions properties for the Levi Civita tensor with the spinors. The function levicivita2up[a,b] assumes this value when a and b are not both integers."
levicivita2Down::usage="This is an auxiliary function needed in order to define the good contractions properties for the Levi Civita tensor with the spinors. The function levicivita2down[a,b] assumes this value when a and b are not both integers."


$contractvariable::usage="Internal variable used in some contractions."
$contractvariable1::usage="Internal variable used in some contractions."
$contractvariable2::usage="Internal variable used in some contractions."
$contractvariable3::usage="Internal variable used in some contractions."
$contractvariable4::usage="Internal variable used in some contractions."


Momenta::usage = "Momenta[\*StyleBox[\"expr\", \"TI\"]] extracts the name of the momenta appearing in brackes in \*StyleBox[\"expr\", \"TI\"]."
AllMomenta::usage="Same as Momenta but distinguishes between the \[Lambda] from the \[Mu]"


SpinorReplace::usage="SpinorReplace[exp,subs] applies the replacements subs in exp. Notice that subs must be given in terms of pure spinors without index structure"
SpinorReplaceSequential::usage="SpinorReplaceSequential[exp,subs] works like SpinorReplace but does the replacements sequentially, one after another in the given order."
SubCounter::usage="Auxiliary function used inside SpinorReplace"
$crep::usage="Auxiliary global variable needed in SubCounter"


ConvenientMu::usage="ConvenientMu[exp] return a convenient substitution of the first encountered \[Mu] spinor in terms of other spinors, in order to kill as much terms in exp as possible."


MuReplace::usage="MuReplace[exp] replaces the \[Mu] and \[Mu] tilde in exp with the most convenient ones in order to kill as much terms in exp as possible. It admits the option DisplayReplacements->True to show the replacements."
DisplayReplacements::usage="DisplayReplacements is an option for MuReplace. If set to True the replacements chosen by MuReplace are shown along with the result. Default is False."
AllowMutoMu::usage="AllowMutoMu is an option for MuReplace. If set to True replecements of \[Mu] with another \[Mu] become allowed. Default is False."
GlobalReplacements::usage="GlobalReplacements is an option for MuReplace. If set to True it makes the replecemnts globally defined, and stores the in FixedSpinors. Default is False."


FixedSpinors::usage="FixedSpinors is a global variable in which all the arbitrary constraints imposed on the spinors are stored."
ClearSpinors::usage="ClearSpinors[] clears all the arbitrary constraints enforced on the spinors."
FixSpinors::usage="FixSpinors[{\[Mu][a]->\[Lambda][b]}] fixes \[Mu][a] to be equal to \[Lambda][b]. FixSpinors needs a list of pure spinor replacements as input. In order to clear the given definitions use ClearSpinors."


SchoutenSimplify::usage = "SchoutenSimplify[\*StyleBox[\"expr\", \"TI\"]] simplifies the expression \*StyleBox[\"expr\", \"TI\"] by using Schouten identities."


Mom4D::usage="Mom4D[label][type][undotted][dotted] produces expression of the 4D massive momentum label with indices undotted and dotted, which are upper or lower indices based on type, which should be either $up or $down"
Mom4DN::usage="Mom4DN[label][type] produces numerical expression of the 4D massive momentum label in matric form."
Mom4DVector::usage="Mom4DVectorN[label][type][lindex] produces numerical expression of the 4D massive momentum label with Lorentz index lindex. The argument type should be either $up or $down, depending on position of Lorentz index."
Mom4DVectorN::usage = "";

LorentzProduct4D::usage = "";

Mom4DVectorNFromMatrix::usage="Creates a Mom4DVector from previously generated spinors."
Mom4DNFromVector::usage="";
Mom4DNFromVectorInPlace::usage = "";
SpinorsFromMom4DN::usage = "Generate spinors from the Mom4DN[][$flat] matrix of the provided momentum name.";
$flat::usage=""
$up::usage="Variable needed to specify upper or lower indices";
$down::usage="Variable needed to specify upper or lower indices";
pm::usage = "Momentum label marker that signals a sign flip on a momentum";

$PMPMMetric = False;
Metric::usage = "";
PauliMatrixSpinor::usage = "Paulimatrices defined to match the $up/$down convention in this package.";

Infycheck::usage="Infycheck[exp] returns True if the expression is divergent and False otherwise.";

SafeDivide::usage="SafeDivide[x,y] returns zero if PossibleZeroQ[x]===True, otherwise x/y.";


S6::usage="S6[p1,p2] is the six dimensional mandelstam invariant s=(p1+p2)^2 expressed in terms of the four dimensional pieces";
S6Mom::usage = "Six dimensional Mandelstam invariant using Mom4DVector";
S6MomN::usage = "Six dimensional Mandelstam invariant using Mom4DVector, numerically evaluated.";
S6N::usage="S6N[p1,p2] is the six dimensional mandelstam invariant s=(p1+p2)^2, numerically evaluated";
S6many::usage="S6many[p1,p2,...,pn] is the six dimensional mandelstam invariant s=(p1+p2+...+pn)^2 in terms of the 4d objects";
S6manyN::usage="S6many[p1,p2,...,pn] is the six dimensional mandelstam invariant s=(p1+p2+...+pn)^2, numerically evaluated";

S4Mom::usage = "Four dimensional Mandelstam invariant using Mom4DVector";
S4MomN::usage = "Four dimensional Mandelstam invariant using Mom4DVector, numerically evaluated";

CompleteDenominators::usage="CompleteDenominators[exp] returns exp over a single denominator where spinor products have been completed to Mandelstam invariants.";
CompleteMandelstam::usage="CompleteMandelstam[exp] returns exp where as many spinor products as possible have been replaced with the corresponding four-dimensional Mandelstam invariants. CompleteMandelstam does only gather existing spinor products without adding new terms."
ToChain::usage="ToChain[exp] returns exp where angle and square brackets are gathered together in chains of invariants.";
Chain::usage="Chain[type1,first,momList,last,type2] is the invariant obtained chaining together angle and square brackest. Type1 and type2 assume values $angle or $bracket and represent the type of bracket with which the invariants start or end. First is the first momentum in the invariant, last is the last and momList are all the intermediate ones."
ChainN::usage="ChainN[type1,first,momList,last,type2] is the invariant obtained chaining together angle and square brackest, numerically evaluated.";
chain::usage="chain is an auxiliary function for Chain. It has all the contraction properties one would expect from Chain.";
ReverseChain::usage = "";

PlusM::usage = "Construct for summing momenta in Chains.";
ExpandPlusM::usage = ""

S::usage="Global variable used to Label six-dimensional Mandelstam invariants."
S4::usage="Global variable used to Label four-dimensional Mandelstam invariants."
ChainToSpinor::usage="ChainToSpinor[exp] transforms chain objects into contracted spinors in exp."


$angle::usage="Global variable of SpinorHelicity6D, labels angle brackets."
$square::usage="Global variable of SpinorHelicity6D, labels square brackets."


mp::usage="mp[p_a,p_b] is the scalar product fo momenta p_a and p_b. mp[] has the attribute Orderless."
eps::usage="eps[p_a,p_b,p_c,p_d] is the Levi-Civita tensor contracted with the four momenta pa_,p_b,p_c and p_d."
TrG::usage="TrG[mom_List] is the trace over the list of slashed momenta mom."
TrG5::usage="TrG5[mom_List] is the trace over the slashed momenta of the list mom but has also a Gamma 5 matrix inserted in the first position of the trace."
ToTrace::usage="ToTrace[exp] converts all closed chains in exp into traces and evaluates them. It admits the Option KillEpsilon whose default is False. If KillEpsilon is set to True the terms in the traces containing an epsilon tensor are set to zero."
ScalProdToS::usage="ScalProdToS[exp] convert all the four-dimensional scalar products in exp to six dimensional mandelstam invariants S taking into account possible masses."


Relabel::usage="Relabel[exp_,{LabelRep->{lab1->lab2,...},ScalarProduct->funtion}] allows to perform relabelings of the momenta inside the functions S[], S4[] and mp[]. These relabelings are specified by the optional argument LabelRep. Furthermore it allows to replace the definition of the scalar product given in SpinorHelicity6D with any given function (also optional feature). This is intended to be used when interfacing the results of SpinorHelicity6D with other Mathematica package. For example by setting ScalarProduct->sp will allow direct input of the results into LiteRed for integral reductions."
LabelRep::usage="Option for the function Relabel. LabelRep specifies the label and functional replacements to be performed in a given expression. For example LabelRep->{1->p1,2->p2,S[1,2]->s12,...}."
ScalarProduct::usage="Option for the function Relabel which specifies if the scalar product has to be redefined to a new function. For example ScalarProduct->sp sets the scalar product to be the function sp (used in LiteRed for example)."


SumContracted::usage="SumContracted[index__][exp_] performs the sum over the contracted little group indices index in the expression exp."
StateSumParallel::usage="StateSumParallel is an option for SumContracted. SumContracted[a1,...,an,{StateSumParallel->True},exp] evaluates the state sum in parallel on all the avaialble parallel kernels. If no parallel kernel was generated so far, it launches the maximum number of parallel kernels and proceeds with the computation."
DisplayKernel::usage="DisplayKernel is an option for SumContracted, if set to True the kernel ID of the parallel kernels is displayed in the final result, allowing to backtrack which kernel evaulated what."
CompleteToMassive::usage="CompleteToMassive[exp_,reps_] restores the \[Mu] dependence in exp after this has been removed by choosing specific values for them. It uses the replacements reps to reinsert the \[Mu] backwards. The replacemnt argument is optional, if no replacement is given the lements of FixedSpinors are used by default."


SpinorUndotN::usage="SpinorUndotN[label][type][indexpos] is the numeric equivalent of SpinorUndot. The variable type is either $lam or $mu and the index position is either $up or $down. This function is used to store the numeric values of the spinors generated by GenSpinors"


SpinorDotN::usage="SpinorUndotN[label][type][indexpos] is the numeric equivalent of SpinorDot. The variable type is either $lam or $mu and the index position is either $up or $down. This function is used to store the numeric values of the spinors generated by GenSpinors"


SpinorAngleBracketN::usage="SpinorAngleBracketN[x,y] returns the numeric value of the angle bracket of the two spinors x and y."


SpinorSquareBracketN::usage="SpinorSQuareBracketN[x,y] returns the numeric value of the square bracket of the two spinors x and y."


RaiseLowerMom4DNSpinorIndices::usage = "";


GenSpinors::usage="GenSpinors[{x1,x2,...,xn},{Options}] generates numerical values for the spinors corresponding to on-shell, conserved, complex kinematics. The spinors are labelled by x1,x2,...,xn and the function allows for the following options: {Dimension,DisplaySpinors,Parametric,ParameterName,ParameterRange,SpecialKinematics,Seed}. For further details type ?OptionName."
GenSpinors::NotPrime="GenSpinors was called with the FiniteField option, but the value of ParameterRange is not a prime number."
GenSpinors::UnsupportedDimension="GenSpinors currently only supports four and six dimensional spinors."
GenSpinors::FiniteFieldParametric="GenSpinors currently does not support parametrized kinematics in finite fields."
GenSpinors::FixedReferenceParametric="GenSpinors does not support currently support parametrized kinematics with fixed refernce spinors for the 6D momenta."
GenSpinors::FixedReferenceBadSyntax="Reference spinors are specified in the wrong format. Required: Rule[label,{\[Xi],\[Eta],\[Xi]t,\[Eta]t}]."
Dimension::usage="Dimension is an option for GenSpinors. It specifies the dimension of the generated kinematics. Default is 6, allowed are 6 and 4."
DisplaySpinors::usage="DisplaySpinors is an option for GenSpinors. If set to True the generated kinematics is displayed. Default is False."
Parametric::usage="Parametric is an option for GenSpinors. If set to True the kinematics is generated in terms of a minimal set of independent variables (3n-10 in 4 dimensions and 5n-15 is six-dimensions) instead of numbers. Default is False."
ParameterName::usage="ParameterName is an option for GenSpinors. It allows to choose a custom label for the independent variables in terms of which the kinematics is defined if Parametric->True. Default is $par."
ParameterRange::usage="ParameterRange is an option for GenSpinors. It allows to set the range in which the independent variables are chosen when numerical values of the kinmeatics are required. Default is 10000."
SpecialKinematics::usage="SpecialKinematics is an option for GenSpinors. It allows to to chose between rational (Rational), finite field (FiniteField) and real (False) kinematics. Default is Rational. "
Seed::usage="Seed is an option for GenSpinors. It allows to set a seed for SeedRandom so that the randomly generated kinematics can be rederived if necessary. Default is False, but any integer number is admitted."
RandomSpinors::usage="RandomSpinors is an option for GenSpinors. It allows to generate a random set of spinors, which lead to on-shell momenta but these do not satisfy momentum conservation"
$par::usage="Protected symbol. It is the default name of the variables in GenSpinors."
FiniteField::usage="Protected symbol. It is the option name for the SpecialKinematics option of GenSpinors to generate finite field valued spinors."
FourD::usage = "Option for GenSpinors to designate momenta that are four-dimensional"

ClearKinematics::usage="ClearKinematics clears all the so far generated and stored numerical values for the kinematics."


ExtramassN::usage="ExtramassN[label] is the numerical equivalent of Extramass[label]."
ExtramasstildeN::usage="ExtramasstildeN[label] is the numerical equivalent of Extramasstilde[label]."


ToNum::usage="TuNum[exp] return numeric value of exp. It requires some kind of numerical kinematics to be generated first using GenSpinors."


Cyclic::usage="Cyclic[exp,cycle_List,OptionsPattern[]] allows to generate a cyclic permutation of the labels in a given expression. The input cycle must be of the form {x1,x2,...,xn} which means x1 is replaced with x2, x2 with x3 and so on, or {{x1,x2,...},...,{...xn}} where the replacements are done in the subcycles. The option SumAll is allowed in OptionPattern."


SumAll::usage="SumAll is an option for Cyclic, which if set to True allows to not only generate the cyclic permutations but also to sum over them. Default is false"


MDelta::usage="MDelta[dim][up1,up2][down1,down2] is the Kronecker delta in dimension dim with upper indices up1 and up2 and lower indices down1 down2. Not present indices are replaced by Null"


MDeltaBox::usage="Just the boxing function for MDelta"
$MDimension::usage="Global variable setting the dimension of MDelta. Default is 4."


Antisymmetrize::usage="Antisymmetrize[exp,{A1,...,An},{B1,...,Bm},...] returns exp antisymmetrized on the indices in the index lists. Be careful, the indices in the list need to be in the same order as in exp for all the signs to be correct in the output."


Fstrength::usage="Fstrength[momlabel][A,B][C,D][lg,lgdot] is the field strength of mmentum momlabel, with A,B Lorentz indices transforming in the fundamental representation and C,D in th eantifundamental, and lg,lgdot little group indices."
FstrengthBox::usage="FstrengthBox is just the boxing function for Fstrength."


ContractReplace::usage="ContractReplace[exp] performs simplifications in exp by recognizing terms which are contracted with an \[Epsilon] tensor and differ by a simple relabelling of indices. Pay attention, it works only on <1_a,2_b] invariants and could be rather inefficient in terms of time...Will hopefully be inproved in later versions of the package."


Fst::usage="Fst[mom][mu,nu][Null,Null] is the field strength tensor with standard Lorentz Indices mu nu and momentum label mom"
Eta::usage="Eta[mu,nu][$up/$down] is the flat metric tensor with upper/lower indices mu nu"
Polar::usage="Polar[mom,refmom][mu][Null] is the polarization vector for the particle with momentum mom and reference momentum refmom, with upper index mu"
PolarPure::usage="PolarPure[mom,refmom] is the polarization vector for the particle with momentum mom and reference momentum refmom with Lorentz indices stripped off"
Mom::usage="Mom[mom][mu][Null] represents the momentum for the particle with momentum mom and upper Lorentz index mu"
MomPure::usage="MomPure[mom] represents the momentum of the particle with momentum mom with Lorentz index stripped off"


FeynCToSpinor::usage="FeynCToSpinor[exp,momenta_List] converst exp from FeynClac notation to SPinorHelicity6D notation. It requires a list of the momenta in exp to work properly, also it is just a first version covering parts of the FeynCalc functions. Any suggestions for improvements are welcome."


MomReplace::usage="MomReplace[exp,replacements_List] is specifically designed to apply replacements of momenta which are contracted into chains or into Mandelstam invariants. Itadmits the option NumeratorOnly"


NumeratorOnly::usage="NumeratorOnly is an option for MomReplace. It allows to perform replacements only in the numerator if set to True (default), or on the whole expression if set to False."


ToFile::usage="ToFile[exp,filename] writes to the file filename.txt located in the documents directory all the frontend oprinted output appearing in exp, just returning the output of exp."


SpinorPalette::usage="Opens the palette associated to the package SpinorHelicity6D."


ClearSubValues::usage="ClearSubValues[f] clears all the SubValues of f leaving downvalues and upvalues intact."


(* ::Section:: *)
(*Private: 6D spinor helicity definitions*)


Begin["`Private`"]


(*We define a private variable needed for the package to decide whether to run shortcuts and the palette or not. This is related to the availability of a frontend.*)


frontend=If[TrueQ[$FrontEnd==Null],0,1];


(* ::Subsection::Closed:: *)
(*Infycheck*)


Attributes[Infycheck]={HoldAll};
Infycheck[x_]:=TrueQ[Quiet[Check[x,$Failed,{PowerMod::ninv,Power::infy,Infinity::indet,Power::indet}],{PowerMod::ninv,Power::infy,Infinity::indet,Power::indet}]==$Failed];


(* ::Subsection::Closed:: *)
(*SafeDivide*)


SafeDivide[x_,y_]:=If[PossibleZeroQ[x],0,x/y];


(* ::Subsection::Closed:: *)
(*overbar*)


overbar[Times[-1,x_]]:=-1*OverBar[x];
overbar[x_Integer]:=If[x>=0,OverBar[x],-OverBar[-x]];
overbar[n_]:=OverBar[n];


(* ::Subsection::Closed:: *)
(*MinuSignQ*)


(*MinusSignQ[x_]:=If[StringContainsQ[ToString[x],"-"],True,False];*)
MinusSignQ[-x_]:=True;
MinusSignQ[x_?Negative]:=True;
MinusSignQ[x_]:=False


(* ::Subsection::Closed:: *)
(*Extra mass*)


ExtramassBox[x_]:=TemplateBox[{x},"extramass",
DisplayFunction->(RowBox[{"\[CapitalMu]","[",#,"]"}]&),
InterpretationFunction->(RowBox[{"Extramass","[",#,"]"}]&)
];

extramass[x_]:=Module[{input},
input=x/.OverBar->overbar;
If[MinusSignQ[input],Return[-Extramass[-input]],Return[Extramass[input]]];
];

Extramass /: MakeBoxes[Extramass[x_],StandardForm|TraditionalForm]:=ExtramassBox[ToBoxes[x]];

Extramass[UnderBar[x_]] := Extramass[x];
ExtramassN[UnderBar[x_]] := ExtramassN[x];

(* ::Subsection::Closed:: *)
(*Extra Mass tilde*)


ExtramasstildeBox[x_]:=TemplateBox[{x},"extramasstilde",
DisplayFunction->(RowBox[{OverscriptBox["\[CapitalMu]","~"],"[",#,"]"}]&),
InterpretationFunction->(RowBox[{"Extramasstilde","[",#,"]"}]&)
];

extramasstilde[x_]:=Module[{input},
input=x/.OverBar->overbar;
If[MinusSignQ[input],Return[-Extramasstilde[-input]],Return[Extramasstilde[input]]];
];

Extramasstilde /: MakeBoxes[Extramasstilde[x_],StandardForm|TraditionalForm]:=ExtramasstildeBox[ToBoxes[x]];

Extramasstilde[UnderBar[x_]] := Extramasstilde[x];
ExtramasstildeN[UnderBar[x_]] := ExtramasstildeN[x];


(* ::Subsection::Closed:: *)
(*KillMasses*)


KillMasses[x:Except[_List,_]]:=KillMasses[{x}];

KillMasses[x_List]:=(
Unprotect[Extramass,Extramasstilde,Momenta4D];
Do[
	Extramass[i]=0;
	Extramasstilde[i]=0;
	ExtramassN[i]=0;
	ExtramasstildeN[i]=0;,
{i,x}];
Momenta4D=x;
Protect[Extramass,Extramasstilde,Momenta4D];
);


(* ::Subsection::Closed:: *)
(*Momenta4D*)


Momenta4D={};


(* ::Subsection::Closed:: *)
(*Set4DMassive*)


Set4DMassive[x_List]:=(
KillMasses[x];
Momenta4DMassive=Join[x,Momenta4DMassive];
);

Set4DMassive[x_]:=Set4DMassive[{x}];



(* ::Subsection::Closed:: *)
(*Momenta4D*)


Momenta4DMassive={};


(* ::Subsection::Closed:: *)
(*ClearDownValues*)


ClearDownValues[f_]:=DownValues[f]=DeleteCases[DownValues[f],_?(FreeQ[First[#],Pattern]&)];


(* ::Subsection::Closed:: *)
(*NewProcess*)


NewProcess:=(Unprotect[Extramass,Extramasstilde,Momenta4D];ClearDownValues[Extramass];ClearDownValues[Extramasstilde];Momenta4D={};Momenta4DMassive={};Protect[Extramass,Extramasstilde,Momenta4D];);



(* ::Subsection::Closed:: *)
(*Pauli Matrices*)


(*Extending Pauli Matrices *)
(* The first index specifies the position of the Lorentz index, whereas the second one is used for the position of the SU(2) indices*)

	PauliMatrixSpinor[2, $up][$down] /; $PMPMMetric := I*PauliMatrix[2];
	PauliMatrixSpinor[2, $down][$down] /; $PMPMMetric := I*PauliMatrix[2];
	PauliMatrixSpinor[2, $up][$up] /; $PMPMMetric := -I*PauliMatrix[2];
	PauliMatrixSpinor[2, $down][$up] /; $PMPMMetric := -I*PauliMatrix[2];
	

	PauliMatrixSpinor[k_, $up][$down] /; MemberQ[Range[0, 3], k] := PauliMatrix[k];
	PauliMatrixSpinor[k_, $down][$up] /; MemberQ[Range[0, 3], k] := PauliMatrix[k];
	PauliMatrixSpinor[k_, $down][$down] /; MemberQ[Range[0, 3], k] := PauliMatrix[k] * {1, -1, -1, -1}[[k + 1]];
	PauliMatrixSpinor[k_, $up][$up] /; MemberQ[Range[0, 3], k] := PauliMatrix[k] * {1, -1, -1, -1}[[k + 1]];



(* ::Subsection::Closed:: *)
(* Negative Momenta *)


pm /: UnderBar[pm[x_]] := pm[UnderBar[x]];

pm /: OverHat[pm[x_]] := pm[OverHat[x]];
pm /: pm[pm[x_]] := x;

Chain[a_,b_,{c___,pm[d_],e___},f_,g_] := -Chain[a,b,{c,d,e},f,g];
ChainN[a_,b_,{c___,pm[d_],e___},f_,g_] := -ChainN[a,b,{c,d,e},f,g];

pm /: Chain[a_,pm[b_],c_List,d__] := I*Chain[a,b,c,d];
pm /: Chain[a__,b_List,pm[c_],d_] := I*Chain[a,b,c,d];

pm /: ChainN[a_,pm[b_],c_List,d__] := I*ChainN[a,b,c,d];
pm /: ChainN[a__,b_List,pm[c_],d_] := I*ChainN[a,b,c,d];

pm /: SpinorHelicity6D`SpinorAngleBracket[pm[x_],y_]:=I*SpinorHelicity6D`SpinorAngleBracket[x,y];
pm /: SpinorHelicity6D`SpinorAngleBracket[x_,pm[y_]]:=I*SpinorHelicity6D`SpinorAngleBracket[x,y];
pm /: SpinorHelicity6D`SpinorSquareBracket[pm[x_],y_]:=I*SpinorHelicity6D`SpinorSquareBracket[x,y];
pm /: SpinorHelicity6D`SpinorSquareBracket[x_,pm[y_]]:=I*SpinorHelicity6D`SpinorSquareBracket[x,y];

SpinorHelicity6D`SpinorDotN[pm[x_]][y_][z_][w_] := I*(SpinorHelicity6D`SpinorDotN[x][y][z][w]);
SpinorHelicity6D`SpinorUndotN[pm[x_]][y_][z_][w_] := I*(SpinorHelicity6D`SpinorUndotN[x][y][z][w]);

SpinorHelicity6D`SpinorDotN[pm[x_]][y_][z_] := I*(SpinorHelicity6D`SpinorDotN[x][y][z]);
SpinorHelicity6D`SpinorUndotN[pm[x_]][y_][z_] := I*(SpinorHelicity6D`SpinorUndotN[x][y][z]);


pm /: SpinorHelicity6D`Extramass[pm[x_]]:=-SpinorHelicity6D`Extramass[x];
pm /: SpinorHelicity6D`Extramasstilde[pm[x_]]:=-SpinorHelicity6D`Extramasstilde[x];
pm /: SpinorHelicity6D`ExtramassN[pm[x_]]:=-SpinorHelicity6D`ExtramassN[x];
pm /: SpinorHelicity6D`ExtramasstildeN[pm[x_]]:=-SpinorHelicity6D`ExtramasstildeN[x];

Mom4DN[pm[x_]][y_][z_] := -Mom4DN[x][y][z];
Mom4D[pm[x_]][y_][z_] := -Mom4D[x][y][z];

pm /: S6Mom[a___,x_,b___,pm[x_],d___] := S6Mom[a,b,d];

pm /: S6Mom[moms:pm[_]..] := S6Mom@@Replace[{moms},pm->Identity,{2},Heads->True]


(* ::Subsection::Closed:: *)
(*SpinorUndot*)


(*Define all the boxing functions for the different cases, \[Lambda], \[Mu] with both up and down indices*)
SpinorLaUpBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorLaUp",
DisplayFunction->(RowBox[{SuperscriptBox["\[Lambda]",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinorundot","[",#1,"]","[","$lam","]","[",#2,"]","[","Null","]"}]&)
];
SpinorLaDownBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorLaDown",
DisplayFunction->(RowBox[{SubscriptBox["\[Lambda]",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinorundot","[",#1,"]","[","$lam","]","[","Null","]","[",#2,"]"}]&)
];
SpinorMuUpBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorMuUp",
DisplayFunction->(RowBox[{SuperscriptBox["\[Mu]",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinorundot","[",#1,"]","[","$mu","]","[",#2,"]","[","Null","]"}]&)
];
SpinorMuDownBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorMuDown",
DisplayFunction->(RowBox[{SubscriptBox["\[Mu]",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinorundot","[",#1,"]","[","$mu","]","[","Null","]","[",#2,"]"}]&)
];


(*Define the property of Spinorundot to get an i in front for negative arguments of the momentum label. Pay attention to the difference between SpinorUndot and Spinorundot!*)
Spinorundot[momlabel_][type_][upper_][lower_]:=If[MinusSignQ[momlabel],I*SpinorUndot[-momlabel][type][upper][lower],SpinorUndot[momlabel][type][upper][lower]];


(*Define the action of Makeboxes on our functions*)
SpinorUndot /: MakeBoxes[SpinorUndot[mom_][$lam][upper_][Null],StandardForm|TraditionalForm]:=SpinorLaUpBox[ToBoxes[mom],ToBoxes[upper]];
SpinorUndot /: MakeBoxes[SpinorUndot[mom_][$mu][upper_][Null],StandardForm|TraditionalForm]:=SpinorMuUpBox[ToBoxes[mom],ToBoxes[upper]];
SpinorUndot /: MakeBoxes[SpinorUndot[mom_][$lam][Null][lower_],StandardForm|TraditionalForm]:=SpinorLaDownBox[ToBoxes[mom],ToBoxes[lower]];
SpinorUndot /: MakeBoxes[SpinorUndot[mom_][$mu][Null][lower_],StandardForm|TraditionalForm]:=SpinorMuDownBox[ToBoxes[mom],ToBoxes[lower]];


(*Define all the shortcuts*)
If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lu" -> SpinorLaUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "ld" -> SpinorLaDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "muu" -> SpinorMuUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "mud" -> SpinorMuDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
];


(*Finally define the contraction properties of the Undotted spinors among themselves when repeated indices are encoutered:*)
SpinorUndot /: Times[SpinorUndot[a_][$mu][d_][Null],SpinorUndot[b_][$mu][Null][d_]]:=SpinorAngleBracket[OverBar[a],OverBar[b]];
SpinorUndot /: Times[SpinorUndot[a_][$lam][d_][Null],SpinorUndot[b_][$mu][Null][d_]]:=SpinorAngleBracket[a,OverBar[b]];
SpinorUndot /: Times[SpinorUndot[a_][$mu][d_][Null],SpinorUndot[b_][$lam][Null][d_]]:=SpinorAngleBracket[OverBar[a],b];
SpinorUndot /: Times[SpinorUndot[a_][$lam][d_][Null],SpinorUndot[b_][$lam][Null][d_]]:=SpinorAngleBracket[a,b];
SpinorUndot /: Times[levicivita2Down[a_,b_],SpinorUndot[mom_][type_][b_][Null]]:=SpinorUndot[mom][type][Null][a];
SpinorUndot /: Times[levicivita2Down[a_,b_],SpinorUndot[mom_][type_][a_][Null]]:=-SpinorUndot[mom][type][Null][b];
SpinorUndot /: Times[levicivita2Up[a_,b_],SpinorUndot[mom_][type_][Null][b_]]:=SpinorUndot[mom][type][a][Null];
SpinorUndot /: Times[levicivita2Up[a_,b_],SpinorUndot[mom_][type_][Null][a_]]:=-SpinorUndot[mom][type][b][Null];


(* ::Subsection::Closed:: *)
(*SpinorUndotPure*)


SpinorUndotPureLBox[label_]:=TemplateBox[{label},"SpinorUndotPureL",
DisplayFunction->(RowBox[{"\[Lambda]","[",#,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorUndotPure","[",#,"]","[","$lam","]"}]&)
];
SpinorUndotPureMBox[label_]:=TemplateBox[{label},"SpinorUndotPureM",
DisplayFunction->(RowBox[{"\[Mu]","[",#,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorUndotPure","[",#,"]","[","$mu","]"}]&)
];

SpinorUndotPure /: MakeBoxes[SpinorUndotPure[label_][$lam],TraditionalForm|StandardForm]:=SpinorUndotPureLBox[ToBoxes[label]];
SpinorUndotPure /: MakeBoxes[SpinorUndotPure[label_][$mu],TraditionalForm|StandardForm]:=SpinorUndotPureMBox[ToBoxes[label]];

If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lp" -> SpinorUndotPureLBox["\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "mp" -> SpinorUndotPureMBox["\[SelectionPlaceholder]"]]];
];


(* ::Subsection::Closed:: *)
(*SpinorUndot6D*)


SpinorUndot6D[momlabel_][index_][Null][1]:={-extramass[momlabel]/Spinoranglebracket[momlabel,overbar[momlabel]]*Spinorundot[momlabel][$mu][Null][index],Spinordot[momlabel][$lam][index][Null]};
SpinorUndot6D[momlabel_][index_][Null][2]:={Spinorundot[momlabel][$lam][Null][index],-extramasstilde[momlabel]/Spinorsquarebracket[momlabel,overbar[momlabel]]*Spinordot[momlabel][$mu][index][Null]};

(*Define boxing functions for the SpinorUndot6D*)

SpinorUndot6DBox[momlabel_,index_,lgindex_]:=TemplateBox[{momlabel,index,lgindex},"SpinorUndot6D",
DisplayFunction->(RowBox[{SubscriptBox[SuperscriptBox["\[CapitalLambda]",#2],#3],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorUndot6D","[",#1,"]","[",#2,"]","[","Null","]","[",#3,"]"}]&)
];

SpinorUndot6DBox2[mom_,A_,bdot_]:=TemplateBox[{mom,A,bdot},"SpinorUndot6D2",
DisplayFunction->(RowBox[{SuperscriptBox["\[CapitalLambda]",RowBox[{#2,#3}]],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorUndot6D","[",#1,"]","[",#2,"]","[",#3,"]","[","Null","]"}]&)
];

(*Define the action of MakeBoxes on SpinorUndot6D*)

SpinorUndot6D /: MakeBoxes[SpinorUndot6D[mom_][index_][Null][lgindex_],StandardForm|TraditionalForm]:=SpinorUndot6DBox[ToBoxes[mom],ToBoxes[index],ToBoxes[lgindex]];
SpinorUndot6D /: MakeBoxes[SpinorUndot6D[mom_][index_][lgindex_][Null],StandardForm|TraditionalForm]:=SpinorUndot6DBox2[ToBoxes[mom],ToBoxes[index],ToBoxes[lgindex]];

(*Define shortcut*)

If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "Ld" -> SpinorUndot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
];
(*Define contraction properties with the epsilon tensor*)

(*SpinorUndot6D /: Times[levicivita2Up[a_,b_],SpinorUndot6D[mom_][A_][Null][b_]]:=SpinorUndot6D[mom][A][a][Null];
SpinorUndot6D /: Times[levicivita2Up[b_,a_],SpinorUndot6D[mom_][A_][Null][b_]]:=-SpinorUndot6D[mom][A][a][Null];*)


(* ::Subsection::Closed:: *)
(*SpinorDot*)


(*For an explaination of the following functions see the section SpinorUndot, these work just in the same way*)
SpinorLatUpBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorLatUp",
DisplayFunction->(RowBox[{SuperscriptBox[OverscriptBox["\[Lambda]","~"],#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinordot","[",#1,"]","[","$lam","]","[",#2,"]","[","Null","]"}]&)
];
SpinorLatDownBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorLatDown",
DisplayFunction->(RowBox[{SubscriptBox[OverscriptBox["\[Lambda]","~"],#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinordot","[",#1,"]","[","$lam","]","[","Null","]","[",#2,"]"}]&)
];
SpinorMutUpBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorMutUp",
DisplayFunction->(RowBox[{SuperscriptBox[OverscriptBox["\[Mu]","~"],#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinordot","[",#1,"]","[","$mu","]","[",#2,"]","[","Null","]"}]&)
];
SpinorMutDownBox[mom_,index_]:=TemplateBox[{mom,index},"SpinorMutDown",
DisplayFunction->(RowBox[{SubscriptBox[OverscriptBox["\[Mu]","~"],#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Spinordot","[",#1,"]","[","$mu","]","[","Null","]","[",#2,"]"}]&)
];


Spinordot[momlabel_][type_][upper_][lower_]:=If[MinusSignQ[momlabel],I*SpinorDot[-momlabel][type][upper][lower],SpinorDot[momlabel][type][upper][lower]];


SpinorDot /: MakeBoxes[SpinorDot[mom_][$lam][upper_][Null],StandardForm|TraditionalForm]:=SpinorLatUpBox[ToBoxes[mom],ToBoxes[upper]];
SpinorDot /: MakeBoxes[SpinorDot[mom_][$mu][upper_][Null],StandardForm|TraditionalForm]:=SpinorMutUpBox[ToBoxes[mom],ToBoxes[upper]];
SpinorDot /: MakeBoxes[SpinorDot[mom_][$lam][Null][lower_],StandardForm|TraditionalForm]:=SpinorLatDownBox[ToBoxes[mom],ToBoxes[lower]];
SpinorDot /: MakeBoxes[SpinorDot[mom_][$mu][Null][lower_],StandardForm|TraditionalForm]:=SpinorMutDownBox[ToBoxes[mom],ToBoxes[lower]];


If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "ltu" -> SpinorLatUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "ltd" -> SpinorLatDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "mtu" -> SpinorMutUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "mtd" -> SpinorMutDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
];


(*Finally define the contraction properties of the Undotted spinors among themselves when repeated indices are encoutered:*)
SpinorDot /: Times[SpinorDot[a_][$mu][d_][Null],SpinorDot[b_][$mu][Null][d_]]:=SpinorSquareBracket[OverBar[b],OverBar[a]];
SpinorDot /: Times[SpinorDot[a_][$lam][d_][Null],SpinorDot[b_][$mu][Null][d_]]:=SpinorSquareBracket[OverBar[b],a];
SpinorDot /: Times[SpinorDot[a_][$mu][d_][Null],SpinorDot[b_][$lam][Null][d_]]:=SpinorSquareBracket[b,OverBar[a]];
SpinorDot /: Times[SpinorDot[a_][$lam][d_][Null],SpinorDot[b_][$lam][Null][d_]]:=SpinorSquareBracket[b,a];
SpinorDot /: Times[levicivita2Down[a_,b_],SpinorDot[mom_][type_][b_][Null]]:=SpinorDot[mom][type][Null][a];
SpinorDot /: Times[levicivita2Down[a_,b_],SpinorDot[mom_][type_][a_][Null]]:=-SpinorDot[mom][type][Null][b];
SpinorDot /: Times[levicivita2Up[a_,b_],SpinorDot[mom_][type_][Null][b_]]:=SpinorDot[mom][type][a][Null];
SpinorDot /: Times[levicivita2Up[a_,b_],SpinorDot[mom_][type_][Null][a_]]:=-SpinorDot[mom][type][b][Null];


(* ::Subsection::Closed:: *)
(*SpinorDotPure*)


SpinorDotPureLBox[label_]:=TemplateBox[{label},"SpinorDotPureL",
DisplayFunction->(RowBox[{OverscriptBox["\[Lambda]","~"],"[",#,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorDotPure","[",#,"]","[","$lam","]"}]&)
];
SpinorDotPureMBox[label_]:=TemplateBox[{label},"SpinorDotPureM",
DisplayFunction->(RowBox[{OverscriptBox["\[Mu]","~"],"[",#,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorDotPure","[",#,"]","[","$mu","]"}]&)
];

SpinorDotPure /: MakeBoxes[SpinorDotPure[label_][$lam],TraditionalForm|StandardForm]:=SpinorDotPureLBox[ToBoxes[label]];
SpinorDotPure /: MakeBoxes[SpinorDotPure[label_][$mu],TraditionalForm|StandardForm]:=SpinorDotPureMBox[ToBoxes[label]];

If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "ltp" -> SpinorDotPureLBox["\[SelectionPlaceholder]"]]];
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "mtp" -> SpinorDotPureMBox["\[SelectionPlaceholder]"]]];
];


(* ::Subsection::Closed:: *)
(*SpinorDot6D*)


SpinorDot6D[momlabel_][index_][Null][1]:={extramasstilde[momlabel]/Spinoranglebracket[momlabel,overbar[momlabel]]*Spinorundot[momlabel][$mu][index][Null],-Spinordot[momlabel][$lam][Null][index]};
SpinorDot6D[momlabel_][index_][Null][2]:={Spinorundot[momlabel][$lam][index][Null],-extramass[momlabel]/Spinorsquarebracket[momlabel,overbar[momlabel]]*Spinordot[momlabel][$mu][Null][index]};

(*Define boxing functions*)

SpinorDot6DBox[momlabel_,index_,lgindex_]:=TemplateBox[{momlabel,index,lgindex},"SpinorDot6D",
DisplayFunction->(RowBox[{SubscriptBox[SubscriptBox[OverscriptBox["\[CapitalLambda]","~"],#2],#3],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorDot6D","[",#1,"]","[",#2,"]","[","Null","]","[",#3,"]"}]&)
];

SpinorDot6DBox2[momlabel_,index_,lgindex_]:=TemplateBox[{momlabel,index,lgindex},"SpinorDot6D2",
DisplayFunction->(RowBox[{SuperscriptBox[SubscriptBox[OverscriptBox["\[CapitalLambda]","~"],#2],#3],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"SpinorDot6D","[",#1,"]","[",#2,"]","[",#3,"]","[","Null","]"}]&)
];

(*Define action of MakeBoxes*)

SpinorDot6D /: MakeBoxes[SpinorDot6D[mom_][index_][Null][lgindex_],StandardForm|TraditionalForm]:=SpinorDot6DBox[ToBoxes[mom],ToBoxes[index],ToBoxes[lgindex]];
SpinorDot6D /: MakeBoxes[SpinorDot6D[mom_][index_][lgindex_][Null],StandardForm|TraditionalForm]:=SpinorDot6DBox2[ToBoxes[mom],ToBoxes[index],ToBoxes[lgindex]];

(*Define Shortcut*)

If[frontend==1,
SetOptions[EvaluationNotebook[],InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "Ltd" -> SpinorDot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
];

(*Define contraction with SPinorUndot6D and epsilon tensor*)

SpinorDot6D /: Times[SpinorDot6D[label1_][L_][Null][lg1_],SpinorUndot6D[label2_][L_][Null][lg2_]]:=AngSquInvariant[label2,label1][lg2,lg1];

(*SpinorDot6D /: Times[levicivita2Up[adot_,bdot_],SpinorDot6D[mom_][A_][Null][adot_]]:=-SpinorDot6D[mom][A][bdot][Null];
SpinorDot6D /: Times[levicivita2Up[bdot_,adot_],SpinorDot6D[mom_][A_][Null][adot_]]:=SpinorDot6D[mom][A][bdot][Null];*)


(* ::Subsection::Closed:: *)
(*Levi-Civita Tensor up*)


(*levicivita2upBox[a_,b_]:=TemplateBox[{a,b},"levicivita2up",
DisplayFunction->(SuperscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2up","[",#1,",",#2,"]"}]&)
];

levicivita2up /: MakeBoxes[levicivita2up[a_,b_],TraditionalForm|StandardForm]:=levicivita2upBox[ToBoxes[a],ToBoxes[b]];
(*Boxing of levicivita2Up, actually I just realized that no auxiliary function is needed...*)
levicivita2Up /: MakeBoxes[levicivita2Up[a_,b_],TraditionalForm|StandardForm]:=TemplateBox[{ToBoxes[a],ToBoxes[b]},"levicivita2Up",
DisplayFunction->(SuperscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2Up","[",#1,",",#2,"]"}]&)];*)


(*levicivita2up[a_Integer,b_Integer]:=LeviCivitaTensor[2][[a,b]];
levicivita2up[a_,b_]:=levicivita2Up[a,b];
levicivita2Up[a_, b_] /; (a == b) := 0
levicivita2Up[a_, b_] /; \[Not]OrderedQ[{a,b}] := -levicivita2Up[b, a];

SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lcup" -> levicivita2upBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]]*)


levicivita2upBox[a_,b_]:=TemplateBox[{a,b},"levicivita2up",
DisplayFunction->(SuperscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2up","[",#1,",",#2,"]"}]&)
];

levicivita2up /: MakeBoxes[levicivita2up[a_,b_],TraditionalForm|StandardForm]:=levicivita2upBox[ToBoxes[a],ToBoxes[b]];
(*Boxing of levicivita2Up, actually I just realized that no auxiliary function is needed...*)
levicivita2Up /: MakeBoxes[levicivita2Up[a_,b_],TraditionalForm|StandardForm]:=TemplateBox[{ToBoxes[a],ToBoxes[b]},"levicivita2Up",
DisplayFunction->(SuperscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2Up","[",#1,",",#2,"]"}]&)];


levicivita2up[a_,b_]:=levicivita2Up[a,b];
levicivita2Up[a_, b_] /; (a == b) := 0;
levicivita2Up[a_, b_] /; !OrderedQ[{a,b}] := -levicivita2Up[b, a];
levicivita2Up[a_Integer,b_Integer]:=LeviCivitaTensor[2][[a,b]];

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lcup" -> levicivita2upBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Levi-Civita Tensor down*)


(*levicivita2downBox[a_,b_]:=TemplateBox[{a,b},"levicivita2down",
DisplayFunction->(SubscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2down","[",#1,",",#2,"]"}]&)
];

levicivita2down /: MakeBoxes[levicivita2down[a_,b_],TraditionalForm|StandardForm]:=levicivita2downBox[ToBoxes[a],ToBoxes[b]];
(*Boxing of levicivita2Down, actually I just realized that no auxiliary function is needed...*)
levicivita2Down /: MakeBoxes[levicivita2Down[a_,b_],TraditionalForm|StandardForm]:=TemplateBox[{ToBoxes[a],ToBoxes[b]},"levicivita2Down",
DisplayFunction->(SubscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2Down","[",#1,",",#2,"]"}]&)];*)


(*levicivita2down[a_Integer,b_Integer]:=-LeviCivitaTensor[2][[a,b]];
levicivita2down[a_,b_]:=levicivita2Down[a,b];
levicivita2Down[a_, b_] /; (a == b) := 0
levicivita2Down[a_, b_] /; \[Not]OrderedQ[{a,b}] := -levicivita2Down[b, a];

SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lcd" -> levicivita2downBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]]*)


levicivita2downBox[a_,b_]:=TemplateBox[{a,b},"levicivita2down",
DisplayFunction->(SubscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2down","[",#1,",",#2,"]"}]&)
];

levicivita2down /: MakeBoxes[levicivita2down[a_,b_],TraditionalForm|StandardForm]:=levicivita2downBox[ToBoxes[a],ToBoxes[b]];
(*Boxing of levicivita2Down, actually I just realized that no auxiliary function is needed...*)
levicivita2Down /: MakeBoxes[levicivita2Down[a_,b_],TraditionalForm|StandardForm]:=TemplateBox[{ToBoxes[a],ToBoxes[b]},"levicivita2Down",
DisplayFunction->(SubscriptBox["\[Epsilon]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"levicivita2Down","[",#1,",",#2,"]"}]&)];


levicivita2down[a_,b_]:=levicivita2Down[a,b];
levicivita2Down[a_, b_] /; (a == b) := 0
levicivita2Down[a_, b_] /; !OrderedQ[{a,b}] := -levicivita2Down[b, a];
levicivita2Down[a_Integer,b_Integer]:=-LeviCivitaTensor[2][[a,b]];

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "lcd" -> levicivita2downBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Spinor Angle Bracket*)


SpinorAngleBracket[a_, b_] /; (a == b) := 0
SpinorAngleBracket[a_, b_] /; !OrderedQ[{a,b}] := -SpinorAngleBracket[b, a];

SpinorAngleBracketBox[a_, b_] :=
    TemplateBox[{a, b}, "SpinorAngleBracket",
        DisplayFunction -> (RowBox[{"\[LeftAngleBracket]",RowBox[{#1,"\[MediumSpace]",#2}],"\[RightAngleBracket]"}]&),
        InterpretationFunction -> (RowBox[{"Spinoranglebracket","[",RowBox[{#1,",",#2}],"]"}]&)]

Spinoranglebracket[x_,y_]:=Module[{input1,input2,bool1,bool2},
input1=x/.OverBar->overbar;
bool1=MinusSignQ[input1];
input2=y/.OverBar->overbar;
bool2=MinusSignQ[input2];
Which[bool1&&bool2,
Return[-SpinorAngleBracket[-input1,-input2]],
bool1,
Return[I*SpinorAngleBracket[-input1,input2]],
bool2,
Return[I*SpinorAngleBracket[input1,-input2]],
True,
Return[SpinorAngleBracket[input1,input2]]
];
];

SpinorAngleBracket /: MakeBoxes[SpinorAngleBracket[a_, b_], StandardForm | TraditionalForm] := SpinorAngleBracketBox[ToBoxes[a], ToBoxes[b]]
If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates @ Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "ab" -> SpinorAngleBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Spinor Square Bracket*)


SpinorSquareBracket[a_, b_] /; (a == b) := 0
SpinorSquareBracket[a_, b_] /; \[Not]OrderedQ[{a, b}] := -SpinorSquareBracket[b, a]

SpinorSquareBracketBox[a_, b_] :=
    TemplateBox[{a, b}, "SpinorSquareBracket",
        DisplayFunction -> (RowBox[{"[",RowBox[{#1,"\[MediumSpace]",#2}],"]"}]&),
        InterpretationFunction -> (RowBox[{"Spinorsquarebracket","[",RowBox[{#1,",",#2}],"]"}]&)]

Spinorsquarebracket[x_,y_]:=Module[{input1,input2,bool1,bool2},
input1=x/.OverBar->overbar;
bool1=MinusSignQ[input1];
input2=y/.OverBar->overbar;
bool2=MinusSignQ[input2];
Which[bool1&&bool2,
Return[-SpinorSquareBracket[-input1,-input2]],
bool1,
Return[I*SpinorSquareBracket[-input1,input2]],
bool2,
Return[I*SpinorSquareBracket[input1,-input2]],
True,
Return[SpinorSquareBracket[input1,input2]]
];
];

SpinorSquareBracket /: MakeBoxes[SpinorSquareBracket[a_, b_], StandardForm | TraditionalForm] := SpinorSquareBracketBox[ToBoxes[a], ToBoxes[b]]

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates @ Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "sb" -> SpinorSquareBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Square-angle invariant*)


SquAngInvariant[a_, b_][x_,y_] /; (a == b) := 0;

(*SquAngInvariant[n_,m_][a_Integer,b_Integer]:=SpinorDot6D[n][$contractvariable][Null][a].SpinorUndot6D[m][$contractvariable][Null][b];
*)
SquAngInvariantBox[n_,m_,adot_,b_]:=TemplateBox[{n,m,adot,b},"SquAngInvariant",
DisplayFunction->(RowBox[{"[",SubscriptBox[#1,#3], "," ,SubscriptBox[#2,#4],"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"SquAngInvariant","[",#1,",",#2,"]","[",#3,",",#4,"]"}]&)
];

SquAngInvariant /: MakeBoxes[SquAngInvariant[n_,m_][adot_,b_],StandardForm|TraditionalForm]:=SquAngInvariantBox[ToBoxes[n],ToBoxes[m],ToBoxes[adot],ToBoxes[b]];

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "sai" -> SquAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Angle-square invariant*)


AngSquInvariant[a_, b_][x_,y_] /; (a == b) := 0;

(*AngSquInvariant[n_,m_][a_Integer,b_Integer]:=SpinorUndot6D[n][$contractvariable][Null][a].SpinorDot6D[m][$contractvariable][Null][b];
*)
AngSquInvariantBox[n_,m_,adot_,b_]:=TemplateBox[{n,m,adot,b},"AngSquInvariant",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",SubscriptBox[#1,#3], "," ,SubscriptBox[#2,#4],"]"}]&),
InterpretationFunction->(RowBox[{"AngSquInvariant","[",#1,",",#2,"]","[",#3,",",#4,"]"}]&)
];

AngSquInvariant /: MakeBoxes[AngSquInvariant[n_,m_][adot_,b_],StandardForm|TraditionalForm]:=AngSquInvariantBox[ToBoxes[n],ToBoxes[m],ToBoxes[adot],ToBoxes[b]];

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "asi" -> AngSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Angle-angle invariant*)


AngAngInvariantBox[m_,n_,l_,k_,a_,b_,c_,d_]:=TemplateBox[{m,n,l,k,a,b,c,d},"AngAngInvariant",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",SubscriptBox[#1,#5],",",SubscriptBox[#2,#6],",",SubscriptBox[#3,#7],",",SubscriptBox[#4,#8],"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"AngAngInvariant","[",#1,",",#2,",",#3,",",#4,"]","[",#5,",",#6,",",#7,",",#8,"]"}]&)
];

AngAngInvariant /: MakeBoxes[AngAngInvariant[m_,n_,l_,k_][a_,b_,c_,d_],TraditionalForm|StandardForm]:=AngAngInvariantBox[ToBoxes[m],ToBoxes[n],ToBoxes[l],ToBoxes[k],ToBoxes[a],ToBoxes[b],ToBoxes[c],ToBoxes[d]];

(*AngAngInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[k][$contractvariable3][Null][c]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[l][$contractvariable3][Null][d]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[j][$contractvariable3][Null][b]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[k][$contractvariable1][Null][c]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[l][$contractvariable1][Null][d]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[j][$contractvariable1][Null][b]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
*)
If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "aa" -> AngAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Square square invariant*)


SquSquInvariantBox[m_,n_,l_,k_,a_,b_,c_,d_]:=TemplateBox[{m,n,l,k,a,b,c,d},"SquSquInvariant",
DisplayFunction->(RowBox[{"[",SubscriptBox[#1,#5],",",SubscriptBox[#2,#6],",",SubscriptBox[#3,#7],",",SubscriptBox[#4,#8],"]"}]&),
InterpretationFunction->(RowBox[{"SquSquInvariant","[",#1,",",#2,",",#3,",",#4,"]","[",#5,",",#6,",",#7,",",#8,"]"}]&)
];

SquSquInvariant /: MakeBoxes[SquSquInvariant[m_,n_,l_,k_][a_,b_,c_,d_],StandardForm|TraditionalForm]:=SquSquInvariantBox[ToBoxes[m],ToBoxes[n],ToBoxes[l],ToBoxes[k],ToBoxes[a],ToBoxes[b],ToBoxes[c],ToBoxes[d]];

(*SquSquInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[k][$contractvariable3][Null][c]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[l][$contractvariable3][Null][d]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[j][$contractvariable3][Null][b]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[k][$contractvariable1][Null][c]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[l][$contractvariable1][Null][d]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[j][$contractvariable1][Null][b]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
*)
If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "bb" -> SquSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*To4D*)


To4D[exp_]:=Block[{out,SquAngInvariant,AngSquInvariant,AngAngInvariant,SquSquInvariant},
SquAngInvariant[n_,m_][a_Integer,b_Integer]:=SpinorDot6D[n][$contractvariable][Null][a] . SpinorUndot6D[m][$contractvariable][Null][b];
AngSquInvariant[n_,m_][a_Integer,b_Integer]:=SpinorUndot6D[n][$contractvariable][Null][a] . SpinorDot6D[m][$contractvariable][Null][b];
AngAngInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[k][$contractvariable3][Null][c]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[l][$contractvariable3][Null][d]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[j][$contractvariable3][Null][b]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[k][$contractvariable1][Null][c]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[l][$contractvariable1][Null][d]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[j][$contractvariable1][Null][b]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
SquSquInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[k][$contractvariable3][Null][c]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[l][$contractvariable3][Null][d]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[j][$contractvariable3][Null][b]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[k][$contractvariable1][Null][c]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[l][$contractvariable1][Null][d]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[j][$contractvariable1][Null][b]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
out=exp;
Return[out];
];


(* ::Subsection::Closed:: *)
(*To4DAlways*)


To4DAlways[bool_?BooleanQ]:=(
Unprotect[AngSquInvariant,SquAngInvariant,AngAngInvariant,SquSquInvariant];
If[bool,
(*Define the breackdown properties of the six-dimensional spinors*)
	SquAngInvariant[n_,m_][a_Integer,b_Integer]:=SpinorDot6D[n][$contractvariable][Null][a] . SpinorUndot6D[m][$contractvariable][Null][b];
	AngSquInvariant[n_,m_][a_Integer,b_Integer]:=SpinorUndot6D[n][$contractvariable][Null][a] . SpinorDot6D[m][$contractvariable][Null][b];
	AngAngInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[k][$contractvariable3][Null][c]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[l][$contractvariable3][Null][d]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[i][$contractvariable1][Null][a]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[j][$contractvariable3][Null][b]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[k][$contractvariable1][Null][c]]First[SpinorUndot6D[l][$contractvariable2][Null][d]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[l][$contractvariable1][Null][d]]First[SpinorUndot6D[j][$contractvariable2][Null][b]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+First[SpinorUndot6D[j][$contractvariable1][Null][b]]First[SpinorUndot6D[k][$contractvariable2][Null][c]]Last[SpinorUndot6D[i][$contractvariable3][Null][a]]Last[SpinorUndot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
	SquSquInvariant[i_,j_,k_,l_][a_Integer,b_Integer,c_Integer,d_Integer]:=
Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[k][$contractvariable3][Null][c]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[l][$contractvariable3][Null][d]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[i][$contractvariable1][Null][a]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[j][$contractvariable3][Null][b]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[k][$contractvariable1][Null][c]]Last[SpinorDot6D[l][$contractvariable2][Null][d]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[j][$contractvariable4][Null][b]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[l][$contractvariable1][Null][d]]Last[SpinorDot6D[j][$contractvariable2][Null][b]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[k][$contractvariable4][Null][c]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4]+Last[SpinorDot6D[j][$contractvariable1][Null][b]]Last[SpinorDot6D[k][$contractvariable2][Null][c]]First[SpinorDot6D[i][$contractvariable3][Null][a]]First[SpinorDot6D[l][$contractvariable4][Null][d]]*levicivita2up[$contractvariable1,$contractvariable2]levicivita2down[$contractvariable3,$contractvariable4];
,
ClearSubValues[AngSquInvariant];
ClearSubValues[SquAngInvariant];
ClearSubValues[AngAngInvariant];
ClearSubValues[SquSquInvariant];
];
Protect[AngSquInvariant,SquAngInvariant,AngAngInvariant,SquSquInvariant];);


(* ::Subsection::Closed:: *)
(*SpinorReplace*)


$crep=1;
SubCounter:=($crep=$crep+1;\[Alpha]*IntegerPart[$crep/2]);
SpinorReplace[exp_,subs_List]:=Block[{locexp,locsubs,$a,$b},
locexp=exp/.{
SpinorAngleBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]-> HoldForm[SpinorUndot[x][$mu][SubCounter][Null]*SpinorUndot[y][$mu][Null][SubCounter]],SpinorSquareBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]->HoldForm[SpinorDot[x][$mu][Null][SubCounter]*SpinorDot[y][$mu][SubCounter][Null]],
SpinorAngleBracket[x_,y_]/;TrueQ[Head[y]==OverBar]-> HoldForm[SpinorUndot[x][$lam][SubCounter][Null]*SpinorUndot[y][$mu][Null][SubCounter]],SpinorSquareBracket[x_,y_]/;TrueQ[Head[y]==OverBar]->HoldForm[SpinorDot[x][$lam][Null][SubCounter]*SpinorDot[y][$mu][SubCounter][Null]],SpinorAngleBracket[x_,y_]-> HoldForm[SpinorUndot[x][$lam][SubCounter][Null]*SpinorUndot[y][$lam][Null][SubCounter]],SpinorSquareBracket[x_,y_]->HoldForm[SpinorDot[x][$lam][Null][SubCounter]*SpinorDot[y][$lam][SubCounter][Null]]};
locexp=locexp/.{SpinorUndot[OverBar[x_]][$mu][y_][z_]->SpinorUndot[x][$mu][y][z],SpinorDot[OverBar[x_]][$mu][y_][z_]->SpinorDot[x][$mu][y][z]};
locsubs=subs/.{SpinorUndotPure->SpinorUndot,SpinorDotPure->SpinorDot};
(*Do[
locsubs[[i,1]]=locsubs[[i,1]][$a_][$b_];
locsubs[[i,2]]=locsubs[[i,2]][$a][$b];
,{i,Length[subs]}];*)
locexp=locexp/.locsubs;
$crep=1;
Return[locexp//ReleaseHold];
];


(* ::Subsection::Closed:: *)
(*Momenta*)


Momenta[expr_] := Join[Cases[{expr}, HoldPattern[SpinorAngleBracket[a_, b_] | SpinorSquareBracket[a_, b_]] :> Sequence[a, b], \[Infinity]],
Cases[{expr}, HoldPattern[extramass[a_] | extramasstilde[a_]|\[CapitalMu][a_]|OverTilde[\[CapitalMu]][a_]|SpinorUndot[a_][_][_][_]|SpinorDot[a_][_][_][_]] :> a, \[Infinity]]]/.OverBar[x_]->x // DeleteDuplicates


(* ::Subsection::Closed:: *)
(*AllMomenta*)


AllMomenta[expr_] := Join[Cases[{expr}, HoldPattern[SpinorAngleBracket[a_, b_] | SpinorSquareBracket[a_, b_]] :> Sequence[a, b], \[Infinity]],
Cases[{expr}, HoldPattern[extramass[a_] | extramasstilde[a_]|\[CapitalMu][a_]|OverTilde[\[CapitalMu]][a_]|SpinorUndot[a_][_][_][_]|SpinorDot[a_][_][_][_]] :> a, \[Infinity]]] // DeleteDuplicates


(* ::Subsection::Closed:: *)
(*ConvenientMu*)


(*ConvenientMu[exp_,mu_]:=Block[{momlist,mu2,count,tabu,localexp,locmomlist,$SquarePowerCount,$AnglePowerCount,anglecounts,squarecounts,angletarget,squaretarget,anglereplace,squarereplace,locexp},
(*Extract all the momenta appearing in expression*)
momlist=AllMomenta[exp];
(*remove mu from the list*)
momlist=DeleteCases[momlist,x_/;x==mu];

(*tabu is the list of values it cannot assume since it would set the denominator to zero.*)
(*tabu["angle"]={};
tabu["square"]={};*)
tabu["angle"]={mu}/.OverBar[x_]->x;
tabu["square"]={mu}/.OverBar[x_]->x;
localexp=Denominator[Together[exp]];
localexp=localexp//.{SpinorAngleBracket[x_,y_]/;TrueQ[x==mu]==False&&TrueQ[y==mu]==False->1,SpinorSquareBracket[x_,y_]->1};
debugPrint["localexp= ",localexp];
AppendTo[tabu["angle"],DeleteCases[AllMomenta[localexp],x_/;x==mu]];
tabu["angle"]=tabu["angle"]//Flatten;
debugPrint["tabu[angle]= ",tabu["angle"]];

(*Now determine the tabu list in the square bracket*)
localexp=Denominator[Together[exp]];
localexp=localexp//.{SpinorSquareBracket[x_,y_]/;TrueQ[x==mu]==False&&TrueQ[y==mu]==False->1,SpinorAngleBracket[x_,y_]->1};
debugPrint["localexp= ",localexp];
AppendTo[tabu["square"],DeleteCases[AllMomenta[localexp],x_/;x==mu]];
tabu["square"]=tabu["square"]//Flatten;
debugPrint["tabu[square]= ",tabu["square"]];

(*Now pick local expression and remove all the tabu contractions in order not to count them*)
locexp=exp//.{SpinorSquareBracket[x_,y_]/;MemberQ[tabu["square"],x]||MemberQ[tabu["square"],y]->1,SpinorAngleBracket[x_,y_]/;MemberQ[tabu["angle"],x]||MemberQ[tabu["angle"],y]->1};
debugPrint["localexp= ",localexp];

(*count how many times the chosen barred momentum appears in expression contracted with each single one of the other momenta*)
(*counting for the angle bracket*)
Do[
(*now count*)
localexp=locexp//.{SpinorAngleBracket[mu,i]->$AnglePowerCount,SpinorAngleBracket[i,mu]->$AnglePowerCount};
debugPrint["localexp= ",localexp];
count[i,"angle"]=Expand[D[localexp,$AnglePowerCount]];
debugPrint["count[",i,",angle]= ",count[i,"angle"]];
(*If the resulting expression is not a sum then the immediate counting could go wrong, so we have to test the Head*)
Which[TrueQ[Head[count[i,"angle"]]==Plus],
debugPrint["sum in ",i];
count[i,"angle"]=Length[count[i,"angle"]],
TrueQ[count[i,"angle"]==0],
debugPrint["zero in ",i],
True,
debugPrint["not a sum in ",i];
count[i,"angle"]=1;
];
debugPrint["count[",i,",angle]= ",count[i,"angle"]];
,{i,momlist}];

(*Now we just need to identify the momentum with the highest counter*)
debugPrint["momlist= ",momlist];
anglecounts=Table[count[i,"angle"],{i,momlist}];
debugPrint["anglecounts= ",anglecounts];
angletarget=Position[anglecounts,Max[anglecounts]]//Flatten//First;
angletarget=momlist[[angletarget]];
debugPrint["Selected substitution angle= ",angletarget];

(*Now recast the replacement into a form amiable for later use*)
mu2=mu/.OverBar[x_]->x;
If[TrueQ[Head[angletarget]==OverBar],
angletarget=angletarget/.OverBar[x_]->x;
anglereplace={SpinorUndotPure[mu2][$mu]->SpinorUndotPure[angletarget][$mu]},
anglereplace={SpinorUndotPure[mu2][$mu]->SpinorUndotPure[angletarget][$lam]};
];

(*Repeat the procedure for the squares, upon performing the angle substitutions*)
(*Pick local expression and remove all the tabu contractions in order not to count them*)
locexp=SpinorReplace[exp,anglereplace];
locexp=locexp//.{SpinorSquareBracket[x_,y_]/;MemberQ[tabu["square"],x]||MemberQ[tabu["square"],y]->1,SpinorAngleBracket[x_,y_]/;MemberQ[tabu["angle"],x]||MemberQ[tabu["angle"],y]->1};
debugPrint["locexp= ",locexp];

(*Counting the square brackets*)
Do[
(*now count*)
localexp=locexp//.{SpinorSquareBracket[mu,i]->$SquarePowerCount,SpinorSquareBracket[i,mu]->$SquarePowerCount};
(*count[i,"square"]=Exponent[localexp,$SquarePowerCount];*)
debugPrint["localexp= ",localexp];
count[i,"square"]=Expand[D[localexp,$SquarePowerCount]];
debugPrint["count[",i,",square]= ",count[i,"square"]];
(*If the resulting expression is not a sum then the immediate counting could go wrong, so we have to test the Head*)
Which[TrueQ[Head[count[i,"square"]]==Plus],
debugPrint["sum in ",i];
count[i,"square"]=Length[count[i,"square"]],
TrueQ[count[i,"square"]==0],
debugPrint["zero in ",i],
True,
debugPrint["not a sum in ",i];
count[i,"square"]=1;
];
debugPrint["count[",i,",square]= ",count[i,"square"]];
,{i,momlist}];

(*Now we just need to identify the momentum with the highest counter*)
debugPrint["momlist= ",momlist];
squarecounts=Table[count[i,"square"],{i,momlist}];
debugPrint["squarecounts= ",squarecounts];
squaretarget=Position[squarecounts,Max[squarecounts]]//Flatten//First;
squaretarget=momlist[[squaretarget]];
debugPrint["Selected subsitution square= ",squaretarget];
If[TrueQ[Head[squaretarget]==OverBar],
squaretarget=squaretarget/.OverBar[x_]->x;
squarereplace={SpinorDotPure[mu2][$mu]->SpinorDotPure[squaretarget][$mu]},
squarereplace={SpinorDotPure[mu2][$mu]->SpinorDotPure[squaretarget][$lam]};
];
Return[Join[anglereplace,squarereplace]];
];*)


ConvenientMu[exp_,mu_,tabulist_,bool_]:=Module[{momlist,mu2,count,tabu,localexp,locmomlist,$SquarePowerCount,$AnglePowerCount,anglecounts,squarecounts,angletarget,squaretarget,anglereplace,squarereplace,locexp,localtabulist,angletargetlist,squaretargetlist},
(*Extract all the momenta appearing in expression*)
momlist=AllMomenta[exp];
(*remove mu from the list*)
momlist=DeleteCases[momlist,x_/;x==mu];

(*tabu is the list of values it cannot assume since it would set the denominator to zero.*)
(*tabu["angle"]={};
tabu["square"]={};*)
tabu["angle"]={mu}/.OverBar[x_]->x;
tabu["square"]={mu}/.OverBar[x_]->x;
localexp=Denominator[Together[exp]];
localexp=localexp//.{SpinorAngleBracket[x_,y_]/;TrueQ[x==mu]==False&&TrueQ[y==mu]==False->1,SpinorSquareBracket[x_,y_]->1};
debugPrint["localexp= ",localexp];
AppendTo[tabu["angle"],DeleteCases[AllMomenta[localexp],x_/;x==mu]];
tabu["angle"]=tabu["angle"]//Flatten;
debugPrint["tabu[angle]= ",tabu["angle"]];

(*Now determine the tabu list in the square bracket*)
localexp=Denominator[Together[exp]];
localexp=localexp//.{SpinorSquareBracket[x_,y_]/;TrueQ[x==mu]==False&&TrueQ[y==mu]==False->1,SpinorAngleBracket[x_,y_]->1};
debugPrint["localexp= ",localexp];
AppendTo[tabu["square"],DeleteCases[AllMomenta[localexp],x_/;x==mu]];
tabu["square"]=tabu["square"]//Flatten;
debugPrint["tabu[square]= ",tabu["square"]];

(*Now pick local expression and remove all the tabu contractions in order not to count them*)
locexp=exp//.{SpinorSquareBracket[x_,y_]/;MemberQ[tabu["square"],x]||MemberQ[tabu["square"],y]->1,SpinorAngleBracket[x_,y_]/;MemberQ[tabu["angle"],x]||MemberQ[tabu["angle"],y]->1};
If[!bool,
(*We want the \[Mu] to be replaced by a \[Lambda] not another \[Mu], so we set to one all the contractions of two \[Mu]*)
locexp=locexp/.{SpinorAngleBracket[x_,y_]/;TrueQ[Head[x]==OverBar]&&TrueQ[Head[y]==OverBar]->1,SpinorSquareBracket[x_,y_]/;TrueQ[Head[x]==OverBar]&&TrueQ[Head[y]==OverBar]->1};
];
debugPrint["locexp= ",locexp];

(*count how many times the chosen barred momentum appears in expression contracted with each single one of the other momenta*)
(*counting for the angle bracket*)
Do[
(*now count*)
localexp=locexp//.{SpinorAngleBracket[mu,i]->$AnglePowerCount,SpinorAngleBracket[i,mu]->$AnglePowerCount};
debugPrint["localexp= ",localexp];
count[i,"angle"]=Expand[D[localexp,$AnglePowerCount]];
debugPrint["count[",i,",angle]= ",count[i,"angle"]];
(*If the resulting expression is not a sum then the immediate counting could go wrong, so we have to test the Head*)
Which[TrueQ[Head[count[i,"angle"]]==Plus],
debugPrint["sum in ",i];
count[i,"angle"]=Length[count[i,"angle"]],
TrueQ[count[i,"angle"]==0],
debugPrint["zero in ",i],
True,
debugPrint["not a sum in ",i];
count[i,"angle"]=1;
];
debugPrint["count[",i,",angle]= ",count[i,"angle"]];
,{i,momlist}];

(*Now we just need to identify the momentum with the highest counter*)
debugPrint["momlist= ",momlist];
anglecounts=Table[count[i,"angle"],{i,momlist}];
debugPrint["anglecounts= ",anglecounts];
(*Take into account the fact that there might be either no allowed substitution or no mu to be replaced, which would result in anglecounts being a list of zeros*)
(*Also account for the fact that if we already replaced \[Mu][i]\[Rule]\[Mu][j] we cannot replace \[Mu][j]\[Rule]\[Lambda][i], this is what tabulist is needed for.*)
If[Max[anglecounts]==0,
angletarget=mu,
angletargetlist=Position[anglecounts,Max[anglecounts]]//Flatten;
angletargetlist=Table[momlist[[i]],{i,angletargetlist}];
debugPrint["angletargetlist= ",angletargetlist];
If[bool,
localtabulist=Table[i[[1]],{i,tabulist}]/.{SpinorDotPure[x_][$mu]:>0,SpinorUndotPure[x_][$mu]:>x}//DeleteDuplicates;
localtabulist=Join[DeleteCases[localtabulist,0],OverBar/@DeleteCases[localtabulist,0]];
debugPrint["localtabulist for angles= ",localtabulist];
angletargetlist=DeleteCases[angletargetlist,x_/;MemberQ[localtabulist,x]];
];
If[Length[angletargetlist]>0,
angletarget=First[angletargetlist],
angletarget=mu;
];
(*angletarget=momlist[[angletarget]];*)
];
debugPrint["Selected substitution angle= ",angletarget];

(*Now recast the replacement into a form amiable for later use*)
mu2=mu/.OverBar[x_]->x;
If[TrueQ[Head[angletarget]==OverBar],
angletarget=angletarget/.OverBar[x_]->x;
anglereplace={SpinorUndotPure[mu2][$mu]->SpinorUndotPure[angletarget][$mu]},
anglereplace={SpinorUndotPure[mu2][$mu]->SpinorUndotPure[angletarget][$lam]};
];

(*Repeat the procedure for the squares, upon performing the angle substitutions*)
(*Pick local expression and remove all the tabu contractions in order not to count them*)
locexp=SpinorReplace[exp,anglereplace];
locexp=locexp//.{SpinorSquareBracket[x_,y_]/;MemberQ[tabu["square"],x]||MemberQ[tabu["square"],y]->1,SpinorAngleBracket[x_,y_]/;MemberQ[tabu["angle"],x]||MemberQ[tabu["angle"],y]->1};
debugPrint["locexp= ",locexp];

(*Counting the square brackets*)
Do[
(*now count*)
localexp=locexp//.{SpinorSquareBracket[mu,i]->$SquarePowerCount,SpinorSquareBracket[i,mu]->$SquarePowerCount};
(*count[i,"square"]=Exponent[localexp,$SquarePowerCount];*)
debugPrint["localexp= ",localexp];
count[i,"square"]=Expand[D[localexp,$SquarePowerCount]];
debugPrint["count[",i,",square]= ",count[i,"square"]];
(*If the resulting expression is not a sum then the immediate counting could go wrong, so we have to test the Head*)
Which[TrueQ[Head[count[i,"square"]]==Plus],
debugPrint["sum in ",i];
count[i,"square"]=Length[count[i,"square"]],
TrueQ[count[i,"square"]==0],
debugPrint["zero in ",i],
True,
debugPrint["not a sum in ",i];
count[i,"square"]=1;
];
debugPrint["count[",i,",square]= ",count[i,"square"]];
,{i,momlist}];

(*Now we just need to identify the momentum with the highest counter*)
debugPrint["momlist= ",momlist];
squarecounts=Table[count[i,"square"],{i,momlist}];
debugPrint["squarecounts= ",squarecounts];
(*Take into account the fact that there might be either no allowed substitution or no mu to be replaced, which would result in anglecounts being a list of zeros*)
If[Max[squarecounts]==0,
squaretarget=mu,
squaretargetlist=Position[squarecounts,Max[squarecounts]]//Flatten;
squaretargetlist=Table[momlist[[i]],{i,squaretargetlist}];
debugPrint["squaretargetlist= ",squaretargetlist];
If[bool,
localtabulist=Table[i[[1]],{i,tabulist}]/.{SpinorDotPure[x_][$mu]:>x,SpinorUndotPure[x_][$mu]:>0}//DeleteDuplicates;
localtabulist=Join[DeleteCases[localtabulist,0],OverBar/@DeleteCases[localtabulist,0]];
debugPrint["localtabulist for squares= ",localtabulist];
squaretargetlist=DeleteCases[squaretargetlist,x_/;MemberQ[localtabulist,x]];
];
If[Length[squaretargetlist]>0,
squaretarget=First[squaretargetlist],
squaretarget=mu;
];

];
debugPrint["Selected subsitution square= ",squaretarget];
If[TrueQ[Head[squaretarget]==OverBar],
squaretarget=squaretarget/.OverBar[x_]->x;
squarereplace={SpinorDotPure[mu2][$mu]->SpinorDotPure[squaretarget][$mu]},
squarereplace={SpinorDotPure[mu2][$mu]->SpinorDotPure[squaretarget][$lam]};
];
Return[Join[anglereplace,squarereplace]];
];


(* ::Subsection::Closed:: *)
(*SpinorReplaceSequential*)


SpinorReplaceSequential[exp_,subs_]:=Catch[Block[{locexp},
locexp=exp;
Do[If[Infycheck[locexp=SpinorReplace[locexp,{subs[[i]]}]],Throw["The substitution "<>ToString[subs[[i]]]<>" is not allowed"]],{i,Length[subs]}];
Return[locexp];
]];


(* ::Subsection::Closed:: *)
(*MuReplace*)


(*Options[MuReplace]={DisplayReplacements->False};
MuReplace[exp_,OptionsPattern[]]:=Block[{localexp,momlist,replacements,totalreplacements,mu},
localexp=exp;
momlist:=AllMomenta[localexp];
totalreplacements={};
While[MemberQ[momlist,x_/;TrueQ[Head[x]==OverBar]],
mu=First[Select[momlist,TrueQ[Head[#]==OverBar]&,1]];
replacements=ConvenientMu[localexp,mu];
totalreplacements=Join[totalreplacements,replacements];
localexp=SpinorReplace[localexp,replacements];
debugPrint["localexp= ",localexp];
];
Which[OptionValue[DisplayReplacements],
Return[{localexp,totalreplacements}],
OptionValue[DisplayReplacements]==False,
Return[localexp],
True,
Return["Unknown option assigned for DisplayReplacements, please choose among True or False"];
];
];*)


Options[MuReplace]={DisplayReplacements->False,AllowMutoMu->False,GlobalReplacements->False};
MuReplace[exp_,OptionsPattern[]]:=Module[{localexp,momlist,replacements,totalreplacements,mu,mu2,exceptions,candidates,candidatetest,scores,$count,tabulist},
localexp=exp;
momlist:=AllMomenta[localexp];
totalreplacements={};
exceptions={};
While[Length[candidates=DeleteCases[Select[momlist,TrueQ[Head[#]==OverBar]&],x_/;MemberQ[exceptions,x]]]>0,
(*Now among the candidates we want to choose the one which appears more often*)
scores={};
Do[
candidatetest=localexp/.{SpinorAngleBracket[x_,y_]/;TrueQ[x==i]||TrueQ[y==i]->$count,SpinorSquareBracket[x_,y_]/;TrueQ[x==i]||TrueQ[y==i]->$count};
candidatetest=D[candidatetest,$count]//Expand;
AppendTo[scores,Length[candidatetest]];
,{i,candidates}];
debugPrint["The candidates ",candidates," scored ",scores];
debugPrint["Chosen ",candidates[[Position[scores,Max[scores],1]//Flatten//First]]];
mu=candidates[[Position[scores,Max[scores],1]//Flatten//First]];
tabulist={};
If[OptionValue[AllowMutoMu],
mu2=mu/.OverBar[x_]:>x;
tabulist=DeleteCases[totalreplacements,x_/;FreeQ[x[[2]],SpinorDotPure[mu2][$mu]]&&FreeQ[x[[2]],SpinorUndotPure[mu2][$mu]]];
debugPrint["candidate= ",mu];
debugPrint["tabulist= ",tabulist];
];
replacements=ConvenientMu[localexp,mu,tabulist,OptionValue[AllowMutoMu]];
(*If ConvenientMu did not find an appropriate replacement for mu we leave it as it is, so we need to add it to the list of exceptions*)
If[!FreeQ[Table[{i[[1]],i[[2]]},{i,replacements}],{x_,x_}],
AppendTo[exceptions,mu];
];
totalreplacements=Join[totalreplacements,replacements];
localexp=SpinorReplace[localexp,replacements];
debugPrint["localexp= ",localexp];
];
(*Store in global variable fixed spinors if the given option is set to True. We only want to store the real replacements not \[Mu] going to itrself, so we have to remove such elements from our list before saving it.*)
If[OptionValue[GlobalReplacements],
replacements=Table[{i[[1]],i[[2]]},{i,totalreplacements}];
replacements=DeleteCases[Table[{i[[1]],i[[2]]},{i,replacements}],{x_,x_}];
replacements=Table[{i[[1]]->i[[2]]},{i,replacements}]//Flatten;
FixSpinors[replacements];
];
Which[OptionValue[DisplayReplacements],
Return[{localexp,totalreplacements}],
OptionValue[DisplayReplacements]==False,
Return[localexp],
True,
Return["Unknown option assigned for DisplayReplacements, please choose among True or False"];
];
];


(* ::Subsection::Closed:: *)
(*ClearSpinors*)


FixedSpinors={};
ClearSpinors[]:=(Unprotect[SpinorUndot,SpinorDot,FixedSpinors];SubValues[SpinorUndot]=DeleteCases[SubValues[SpinorUndot],_?(FreeQ[First[#],HoldPattern@SpinorUndot[Pattern]]&)];
SubValues[SpinorDot]=DeleteCases[SubValues[SpinorDot],_?(FreeQ[First[#],HoldPattern@SpinorDot[Pattern]]&)];
FixedSpinors={};
Protect[SpinorUndot,SpinorDot,FixedSpinors];)


(* ::Subsection::Closed:: *)
(*ClearSpinorsPrivate*)


FixedSpinors={};
ClearSpinorsPrivate[]:=(Unprotect[SpinorUndot,SpinorDot];SubValues[SpinorUndot]=DeleteCases[SubValues[SpinorUndot],_?(FreeQ[First[#],HoldPattern@SpinorUndot[Pattern]]&)];
SubValues[SpinorDot]=DeleteCases[SubValues[SpinorDot],_?(FreeQ[First[#],HoldPattern@SpinorDot[Pattern]]&)];
Protect[SpinorUndot,SpinorDot];)


(* ::Subsection::Closed:: *)
(*FixSpinors*)


FixSpinors[x__List]:=Module[{len,reps,duplicates,clean,auxlist,auxUndot,auxDot},
len=Length[x];
reps=x;
duplicates=Table[reps[[j,1]],{j,len}];

(*Before assigning the spinor equalities we add the new definitions to the list of fixed spinors, but to do so we need first to remove possible old overlapping definitions*)

Unprotect[FixedSpinors];
Do[FixedSpinors=DeleteCases[FixedSpinors,y_/;MemberQ[duplicates,y[[1]]]],{j,len}];
FixedSpinors=Join[FixedSpinors,x]//Sort;
Protect[FixedSpinors];

(*Now we clear all the definitions of the spinors*)
ClearSpinorsPrivate[];

(*Now fix the spinors to the desired values, also the ones which were already fixed*)
Unprotect[SpinorUndot,SpinorDot];
reps=FixedSpinors/.{SpinorUndotPure->SpinorUndot,SpinorDotPure->SpinorDot};
reps=reps/.Rule->Set;
Protect[SpinorUndot,SpinorDot];
];


(* ::Subsection::Closed:: *)
(*SchoutenSimplify*)


(* Rules for Schouten simplification *)
SchoutenRules =
    {k_. SpinorAngleBracket[a_,c_] SpinorAngleBracket[b_,d_] - k_. SpinorAngleBracket[a_,d_] SpinorAngleBracket[b_,c_] :> k SpinorAngleBracket[a,b] SpinorAngleBracket[c,d],
     k_. SpinorAngleBracket[a_,b_] SpinorAngleBracket[c_,d_] + k_. SpinorAngleBracket[a_,d_] SpinorAngleBracket[b_,c_] :> k SpinorAngleBracket[a,c] SpinorAngleBracket[b,d],
     k_. SpinorAngleBracket[a_,c_] SpinorAngleBracket[b_,d_] - k_. SpinorAngleBracket[a_,b_] SpinorAngleBracket[c_,d_] :> k SpinorAngleBracket[a,d] SpinorAngleBracket[b,c],
     k_. SpinorSquareBracket[a_,c_] SpinorSquareBracket[b_,d_] - k_. SpinorSquareBracket[a_,d_] SpinorSquareBracket[b_,c_] :> k SpinorSquareBracket[a,b] SpinorSquareBracket[c,d],
     k_. SpinorSquareBracket[a_,b_] SpinorSquareBracket[c_,d_] + k_. SpinorSquareBracket[a_,d_] SpinorSquareBracket[b_,c_] :> k SpinorSquareBracket[a,c] SpinorSquareBracket[b,d],
     k_. SpinorSquareBracket[a_,c_] SpinorSquareBracket[b_,d_] - k_. SpinorSquareBracket[a_,b_] SpinorSquareBracket[c_,d_] :> k SpinorSquareBracket[a,d] SpinorSquareBracket[b,c]}

(* Generates all possible Schouten identities associated with a given set of momenta. *)
SchoutenIdentities[momenta_List] :=
    ({a,b,c,d} \[Function] Sequence @@ {SpinorAngleBracket[a,b] SpinorAngleBracket[c,d] +
        SpinorAngleBracket[b,c] SpinorAngleBracket[a,d] + SpinorAngleBracket[c,a] SpinorAngleBracket[b,d] == 0,
      SpinorSquareBracket[a,b] SpinorSquareBracket[c,d] + SpinorSquareBracket[b,c] SpinorSquareBracket[a,d] + SpinorSquareBracket[c,a] SpinorSquareBracket[b,d] == 0}) @@@ Subsets[momenta,{4}]

(* SchoutenSimplify*)
SchoutenSimplify[expr_] :=
    Simplify[expr,
        Assumptions -> SchoutenIdentities[Momenta[expr]],
        TransformationFunctions -> {Automatic, (e \[Function] e /. SchoutenRules)}]


(* ::Subsection::Closed:: *)
(*Mom4D*)


Protect[$flat];

Mom4D[l1_][Null][$up][a_][b_]/;!MemberQ[Momenta4DMassive,l1]:=Spinorundot[l1][$lam][a][Null]Spinordot[l1][$lam][b][Null]+SafeDivide[ExtramassN[l1]*ExtramasstildeN[l1],Spinoranglebracket[l1,OverBar[l1]]Spinorsquarebracket[OverBar[l1],l1]]Spinorundot[l1][$mu][a][Null]Spinordot[l1][$mu][b][Null];
Mom4D[l1_][Null][$down][a_][b_]/;!MemberQ[Momenta4DMassive,l1]:=Spinorundot[l1][$lam][Null][a]Spinordot[l1][$lam][Null][b]+SafeDivide[ExtramassN[l1]*ExtramasstildeN[l1],Spinoranglebracket[l1,OverBar[l1]]Spinorsquarebracket[OverBar[l1],l1]]Spinorundot[l1][$mu][Null][a]Spinordot[l1][$mu][Null][b];
Mom4D[l1_][$flat][$up][a_][b_]/;!MemberQ[Momenta4DMassive,l1]:=Spinorundot[l1][$lam][a][Null]Spinordot[l1][$lam][b][Null];
Mom4D[l1_][$flat][$down][a_][b_]/;!MemberQ[Momenta4DMassive,l1]:=Spinorundot[l1][$lam][Null][a]Spinordot[l1][$lam][Null][b];

(*The switch in indices a and b in the $down case is due to the fact that the first index of Mom4D is always the undotted, while the second is the dotted one.
  Since this does not work for multiplying these matrices, the indices are swapped.*)
Mom4DN[l1_][type_][$up]/;!MemberQ[Momenta4DMassive,l1]:=Table[ToNumEval[Mom4D[l1][type][$up][a][b]],
  {a,1,2},{b,1,2}];
Mom4DN[l1_][type_][$down]/;!MemberQ[Momenta4DMassive,l1]:=Table[ToNumEval[Mom4D[l1][type][$down][b][a]],{a,1,2},{b,1,2}];

Mom4DNFromVector[label_][lIdxPos:$up|$down][mom_List]:=
Module[{idxPauli = First@Complement[{$up,$down},{lIdxPos}]},
      Mom4DN[label][Null][$down]=Sum[mom[[i]]*PauliMatrixSpinor[i-1,idxPauli][$down],{i,1,4}];
      Mom4DN[label][Null][$up]=Sum[mom[[i]]*PauliMatrixSpinor[i-1,idxPauli][$up],{i,1,4}];
];


Mom4DNFromVectorInPlace[label_]:=
Module[{mom = Mom4DVectorN[label][Null][$up]},
      Mom4DN[label][Null][$down]=Simplify@Sum[mom[[i]]*PauliMatrixSpinor[i-1,$down][$up],{i,1,4}];
      Mom4DN[label][Null][$up]=Simplify@Sum[mom[[i]]*PauliMatrixSpinor[i-1,$down][$up],{i,1,4}];
];


Mom4DVector[pm[x_]][y_][z_]:=(-1)*SpinorHelicity6D`Mom4DVector[x][y][z];
Mom4DVectorN[pm[x_]][y_][z_]:=(-1)*Mom4DVectorN[x][y][z];


LorentzProduct4D /: MakeBoxes[LorentzProduct4D[x__],StandardForm|TraditionalForm] :=
    TemplateBox[ToBoxes/@{x},"LorentzProduct4D",
  DisplayFunction->(RowBox[{"(",TemplateSlotSequence[1,"\[CenterDot]"],")"}]&),
  InterpretationFunction->(RowBox[{"LorentzProduct4D","[",TemplateSlotSequence[1,","],"]"}]&)
];

LorentzProduct4DN[x_,y_] := Mom4DVectorN[x][Null][$up].Mom4DVectorN[y][Null][$down];
LorentzProduct4DN[x:Repeated[_,{3,Infinity}]] := Plus@@LorentzProduct4DN@@@Subsets[{x},{2}];




Mom4DVectorNFromMatrix[label_][type_] /; $PMPMMetric :=
Module[{comps=Evaluate[Mom4DN[label][type][$down]]},
   Mom4DVectorN[label][type][$up] = 1/2*{comps[[1,1]]+comps[[2,2]],-comps[[1,2]]-comps[[2,1]],-comps[[2,1]]+comps[[1,2]],comps[[2,2]]-comps[[1,1]]};
   Mom4DVectorN[label][type][$down] = {1,-1,1,-1}*Mom4DVectorN[label][type][$up];
];

Mom4DVectorNFromMatrix[label_][type_]:=
Module[{comps=Evaluate[Mom4DN[label][type][$down]]},
   Mom4DVectorN[label][type][$up] = 1/2*{comps[[1,1]]+comps[[2,2]],-comps[[1,2]]-comps[[2,1]],I*comps[[2,1]]-I*comps[[1,2]],comps[[2,2]]-comps[[1,1]]};
   Mom4DVectorN[label][type][$down] = {1,-1,-1,-1}*Mom4DVectorN[label][type][$up];
];

Options[SpinorsFromMom4DN] = {Normalization->Null,Reference->{0,0,0,0},Rules->{},Assumptions->{},NumericalEvaluationFunction->Identity};

SpinorsFromMom4DN[label_, opts: OptionsPattern[]]:=
Module[{a,b,c,d,mom4d},

	mom4d = Mom4DN[label][$flat][$down];
	
	If[And@@(PossibleZeroQ[OptionValue[NumericalEvaluationFunction][#],Assumptions->OptionValue[Assumptions],Method->"ExactAlgebraics"]&/@(mom4d[[All,1]]/.OptionValue[Rules])),
		b = If[OptionValue[Normalization] === Null, 1, (SpinorUndotN[OptionValue[Normalization]][$lam][$down])[[2]] ];
		
		
		c = mom4d[[1,2]]/b;
		d = mom4d[[2,2]]/b;
		
		a = Which[!PossibleZeroQ[OptionValue[NumericalEvaluationFunction][d/.OptionValue[Rules]],Assumptions->OptionValue[Assumptions]],
					mom4d[[2,1]]/d, 
				  !PossibleZeroQ[OptionValue[NumericalEvaluationFunction][c/.OptionValue[Rules]],Assumptions->OptionValue[Assumptions]],
					mom4d[[1,1]]/c
			];
	,
		a = If[OptionValue[Normalization] === Null, 1, (SpinorUndotN[OptionValue[Normalization]][$lam][$down])[[1]] ];
	
		c = mom4d[[1,1]]/a;
		d = mom4d[[2,1]]/a;
		
		b = Which[!PossibleZeroQ[OptionValue[NumericalEvaluationFunction][d]/.OptionValue[Rules],Assumptions->OptionValue[Assumptions],Method->"ExactAlgebraics"], 
					mom4d[[2,2]]/d, 
				  !PossibleZeroQ[OptionValue[NumericalEvaluationFunction][c]/.OptionValue[Rules],Assumptions->OptionValue[Assumptions],Method->"ExactAlgebraics"], 
					mom4d[[1,2]]/c
			];
	];
	
	SpinorDotN[label][$lam][$down]={c,d};
	SpinorDotN[label][$lam][$up]={d,-c};
	SpinorUndotN[label][$lam][$down]={a,b};
	SpinorUndotN[label][$lam][$up]={b,-a};
	
	SpinorDotN[label][$lam][i:1|2][Null]:=Indexed[SpinorDotN[label][$lam][$up],i];
	SpinorDotN[label][$lam][Null][i:1|2]:=Indexed[SpinorDotN[label][$lam][$down],i];
	SpinorUndotN[label][$lam][i:1|2][Null]:=Indexed[SpinorUndotN[label][$lam][$up],i];
	SpinorUndotN[label][$lam][Null][i:1|2]:=Indexed[SpinorUndotN[label][$lam][$down],i];
	
	If[Head[OptionValue[Reference]]===List,
		SpinorDotN[label][$mu][$down]={-1,1}*OptionValue[Reference][[{4,3}]];
		SpinorDotN[label][$mu][$up]=OptionValue[Reference][[{3,4}]];
		SpinorUndotN[label][$mu][$down]=OptionValue[Reference][[{1,2}]];
		SpinorUndotN[label][$mu][$up]={1,-1}*OptionValue[Reference][[{2,1}]];
	,
		SpinorDotN[label][$mu][$down]=SpinorDotN[OptionValue[Reference]][$lam][$down];
		SpinorDotN[label][$mu][$up]=SpinorDotN[OptionValue[Reference]][$lam][$up];
		SpinorUndotN[label][$mu][$down]=SpinorUndotN[OptionValue[Reference]][$lam][$down];
		SpinorUndotN[label][$mu][$up]=SpinorUndotN[OptionValue[Reference]][$lam][$up];
	];
	
	SpinorDotN[label][$mu][i:1|2][Null]:=Indexed[SpinorDotN[label][$mu][$up],i];
	SpinorDotN[label][$mu][Null][i:1|2]:=Indexed[SpinorDotN[label][$mu][$down],i];
	SpinorUndotN[label][$mu][i:1|2][Null]:=Indexed[SpinorUndotN[label][$mu][$up],i];
	SpinorUndotN[label][$mu][Null][i:1|2]:=Indexed[SpinorUndotN[label][$mu][$down],i];
];


(* ::Subsection::Closed:: *)
(*S6*)


S6[x_,y_]:=-extramass[x]extramasstilde[y]-extramass[y]extramasstilde[x]+Mom4D[x][Null][$up][$contractvariable1][$contractvariable2]Mom4D[y][Null][$down][$contractvariable1][$contractvariable2]//Expand;

S6N[x_,y_]:=ToNumEval[S6[x,y]];

S6N[x__]:=ToNumEval[S6Mom[x]];


(* ::Subsection::Closed:: *)
(*RaiseLowerMom4DNSpinorIndices*)


(* The minus is for one of the LevitCivitas, and the result is transposed, as the dotted and undotted indices switch places in Mom4DN for easy matrix multiplication.
   As a reference, Mom4DN[$down] has indices {dotted, undotted}, and Mom4DN[$up] has {undotted, dotted}. *)
RaiseLowerMom4DNSpinorIndices[mat_/;Dimensions[mat]==={2,2}] := -Transpose[LeviCivitaTensor[2] . mat . LeviCivitaTensor[2]];



(* ::Subsection::Closed:: *)
(*Metric*)


Metric := If[$PMPMMetric,{1,-1,1,-1},{1,-1,-1,-1}];


(* ::Subsection::Closed:: *)
(*S6Mom, S6MomN*)


SetAttributes[S6Mom,Orderless];

S6Mom[a___,UnderBar[c_],b___] := S6Mom[a,c,b];


S6MomN[moms__] := Sum[Mom4DVectorN[i][Null][$up],{i,{moms}}] . Sum[Mom4DVectorN[i][Null][$down],{i,{moms}}] -
                  Sum[ExtramassN[i],{i,{moms}}]*Sum[ExtramasstildeN[i],{i,{moms}}];


(* ::Subsection::Closed:: *)
(*S*)


SetAttributes[S,Orderless];


(* ::Subsection:: *)
(*S4*)


SetAttributes[S4,Orderless];

S4[x_,y_]:=Mom4D[x][Null][$up][$contractvariable1][$contractvariable2]Mom4D[y][Null][$down][$contractvariable1][$contractvariable2]//Expand;

S4N[x_,y_]:=ToNumEval[S4[x,y]];


(* ::Subsection:: *)
(*S4Mom*)


SetAttributes[S4Mom,Orderless];

S4MomN[moms__] := Sum[ExtramassN[i]*ExtramasstildeN[i],{i,{moms}}]+ Sum[2*Mom4DVectorN[i[[1]]][Null][$up] . Mom4DVectorN[i[[2]]][Null][$down],{i,Subsets[{moms},{2}]}];



(* ::Subsection:: *)
(*S6many*)


S6many[labels__]:=Block[{locvar,locexp},
locvar={labels};
locvar=DeleteDuplicates[Sort/@Tuples[locvar,2]];
locvar=DeleteCases[locvar,x_/;Length[DeleteDuplicates[x]]==1];
debugPrint["locvar= ",locvar];
locexp=Plus@@Apply[S6[#1,#2]&,locvar,{1}];
Return[locexp];
];

S6manyN[x_,y_]:=ToNumEval[S6many[x,y]];


(* ::Subsection:: *)
(*S4many*)


S4many[labels__]:=Block[{locvar,locexp},
locvar={labels};
locvar=DeleteDuplicates[Sort/@Tuples[locvar,2]];
locvar=DeleteCases[locvar,x_/;Length[DeleteDuplicates[x]]==1];
debugPrint["locvar= ",locvar];
locexp=Plus@@Apply[S4[#1,#2]&,locvar,{1}];
Return[locexp];
];

S4manyN[x_,y_]:=ToNumEval[S4many[x,y]];


(* ::Subsection::Closed:: *)
(*CompleteDenominators*)


(*CompleteDenominators[exp_]:=Module[{localexp,numtot,dentot,numden},
localexp=Together[exp];
numtot=Numerator[localexp];
dentot=Denominator[localexp];
dentot=dentot/.{SpinorAngleBracket[x_,y_]:>S[x,y]/SpinorSquareBracket[y,x]};
numden=Numerator[dentot]/.{SpinorSquareBracket[x_,y_]:>S[x,y]/SpinorAngleBracket[y,x]};
dentot=numden/Denominator[dentot];
Return[numtot/dentot];
];*)


(* ::Subsection::Closed:: *)
(*CompleteDenominators*)


CompleteDenominators[exp_]:=Module[{localexp,numtot,dentot,numden},
localexp=Together[exp];
numtot=Numerator[localexp];
dentot=Denominator[localexp];
dentot=dentot/.{SpinorAngleBracket[x_,y_]:>S4[x,y]/SpinorSquareBracket[y,x]};
numden=Numerator[dentot]/.{SpinorSquareBracket[x_,y_]:>S4[x,y]/SpinorAngleBracket[y,x]};
dentot=numden/Denominator[dentot];
Return[numtot/dentot];
];


(* ::Subsection::Closed:: *)
(*CompleteMandelstam*)


(*CompleteMandelstam[exp_]:=Module[{localexp,numtot,dentot,numnum,dennum,denden,numden,numi},
localexp=Together[exp];
numtot={Numerator[localexp]}/.{Plus->List}//Flatten;
dentot=Denominator[localexp];

Do[
numi=numtot[[i]];
numi=numi/.{SpinorAngleBracket[x_,y_]:>S[x,y]/SpinorSquareBracket[y,x]};
numnum=Numerator[numi];
dennum=Denominator[numi]/.{SpinorSquareBracket[x_,y_]:>S[x,y]/SpinorAngleBracket[y,x]};
numi=numnum/dennum;
numtot[[i]]=numi;
,{i,Length[numtot]}];
numtot=Plus@@numtot;

dentot=dentot/.{SpinorAngleBracket[x_,y_]:>S[x,y]/SpinorSquareBracket[y,x]};
numden=Numerator[dentot];
denden=Denominator[dentot]/.{SpinorSquareBracket[x_,y_]:>S[x,y]/SpinorAngleBracket[y,x]};
dentot=numden/denden;

Return[numtot/dentot];
];*)


(* ::Subsection::Closed:: *)
(*CompleteMandelstam*)


CompleteMandelstam[exp_]:=Module[{localexp,numtot,dentot,numnum,dennum,denden,numden,numi},
localexp=Together[exp];
numtot={Numerator[localexp]}/.{Plus->List}//Flatten;
dentot=Denominator[localexp];

Do[
   numi=numtot[[i]];
   numi=numi/.{SpinorAngleBracket[x_,y_]:>S4[x,y]/SpinorSquareBracket[y,x]};
   numnum=Numerator[numi];
   dennum=Denominator[numi]/.{SpinorSquareBracket[x_,y_]:>S4[x,y]/SpinorAngleBracket[y,x]};
   numi=numnum/dennum;
   numtot[[i]]=numi;
,{i,Length[numtot]}];
numtot=Plus@@numtot;

dentot=dentot/.{SpinorAngleBracket[x_,y_]:>S4[x,y]/SpinorSquareBracket[y,x]};
numden=Numerator[dentot];
denden=Denominator[dentot]/.{SpinorSquareBracket[x_,y_]:>S4[x,y]/SpinorAngleBracket[y,x]};
dentot=numden/denden;

Return[numtot/dentot];
];


(* ::Subsection::Closed:: *)
(*ScalProdToS*)


ScalProdToS[exp_]:=exp/.{mp[i_,j_]:>S[i,j]/2+(extramass[i]extramasstilde[j]+extramass[j]extramasstilde[i])/2};


(* ::Subsection::Closed:: *)
(*chain*)


(*Contraction properties of the chains*)
(*SquareAngle to AngleSquare*)
chain[$square,x_,k_List,y_,$angle]:=(-1)^(Length[k]+1)chain[$angle,y,Reverse[k],x,$square];
(*AngleAngle with SquareSquare*)
chain /: Times[chain[$angle,x_,k_List,y_,$angle],chain[$square,y_,q_List,z_,$square]]:=chain[$angle,x,Join[k,{y},q],z,$square];
chain /: Times[chain[$angle,x_,k_List,y_,$angle],chain[$square,z_,q_List,y_,$square]]:=(-1)^(Length[q]+1)chain[$angle,x,Join[k,{y},Reverse[q]],z,$square];
chain /: Times[chain[$angle,y_,k_List,x_,$angle],chain[$square,y_,q_List,z_,$square]]:=(-1)^(Length[k]+1)chain[$angle,x,Join[Reverse[k],{y},q],z,$square];
chain /: Times[chain[$angle,x_,k_List,y_,$angle],chain[$square,z_,q_List,x_,$square]]:=chain[$square,z,Join[q,{x},k],y,$angle];
(*AngleSquare with SquareSquare*)
(*chain /: Times[chain[$square,x_,k_List,y_,$angle],chain[$square,y_,q_List,z_,$square]]:=chain[$square,x,Join[k,{y},q],z,$square];
chain /: Times[chain[$square,x_,k_List,y_,$angle],chain[$square,z_,q_List,y_,$square]]:=(-1)^(Length[q]+1)chain[$square,x,Join[k,{y},q],z,$square];*)
chain /: Times[chain[$angle,y_,k_List,x_,$square],chain[$square,y_,q_List,z_,$square]]:=(-1)^(Length[k]+1)chain[$square,x,Join[Reverse[k],{y},q],z,$square];
chain /: Times[chain[$angle,y_,k_List,x_,$square],chain[$square,z_,q_List,y_,$square]]:=chain[$square,z,Join[q,{y},k],x,$square];
(*AngleSquare with AngleAngle*)
chain /: Times[chain[$angle,x_,k_List,y_,$square],chain[$angle,y_,q_List,z_,$angle]]:=chain[$angle,x,Join[k,{y},q],z,$angle];
chain /: Times[chain[$angle,x_,k_List,y_,$square],chain[$angle,z_,q_List,y_,$angle]]:=(-1)^(Length[q]+1)chain[$angle,x,Join[k,{y},Reverse[q]],z,$angle];
(*chain /: Times[chain[$square,y_,k_List,x_,$angle],chain[$angle,y_,q_List,z_,$angle]]:=(-1)^(Length[k]+1)chain[$angle,x,Join[k,{y},q],z,$angle];
chain /: Times[chain[$square,y_,k_List,x_,$angle],chain[$angle,z_,q_List,y_,$angle]]:=chain[$angle,z,Join[q,{y},k],x,$angle];*)
(*AngleSquare with AngleSquare*)
chain /: Times[chain[$angle,x_,k_List,y_,$square],chain[$angle,y_,q_List,z_,$square]]:=chain[$angle,x,Join[k,{y},q],z,$square];
(*Vanishing chains*)
chain[type_,x_,{x_,y___},z_,type2_]:=0;
chain[type_,x_,{y___,z_},z_,type2_]:=0;
chain[type_,x_,{y___,z_,z_,k___},l_,type2_] /; Head[z]=!=UnderBar :=0;

(*Display of the chains*)
AngleSquarechainBox[beginning_,args_,end_]:=TemplateBox[{beginning,args,end},"AngleSquarechain",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",#1,#2,#3,"]"}]&),
InterpretationFunction->(RowBox[{"chain","[","$angle",",",#1,",",#2,",",#3,",","$square","]"}]&)
];
SquareAnglechainBox[beginning_,args_,end_]:=TemplateBox[{beginning,args,end},"SquareAnglechain",
DisplayFunction->(RowBox[{"[",#1,#2,#3,"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"chain","[","$square",",",#1,",",#2,",",#3,",","$angle","]"}]&)
];
AngleAnglechainBox[beginning_,args_,end_]:=TemplateBox[{beginning,args,end},"AngleAnglechain",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",#1,#2,#3,"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"chain","[","$angle",",",#1,",",#2,",",#3,",","$angle","]"}]&)
];
SquareSquarechainBox[beginning_,args_,end_]:=TemplateBox[{beginning,args,end},"SquareSquarechain",
DisplayFunction->(RowBox[{"[",#1,#2,#3,"]"}]&),
InterpretationFunction->(RowBox[{"chain","[","$square",",",#1,",",#2,",",#3,",","$square","]"}]&)
];
chain /: MakeBoxes[chain[$angle,x_,y_List,z_,$square],StandardForm|TraditionalForm] /;OddQ[Length[y]+2]:=AngleSquarechainBox[ToBoxes[x],ToBoxes[y],ToBoxes[z]];
chain /: MakeBoxes[chain[$angle,x_,y_List,z_,$angle],StandardForm|TraditionalForm] /;EvenQ[Length[y]+2]:=AngleAnglechainBox[ToBoxes[x],ToBoxes[y],ToBoxes[z]];
chain /: MakeBoxes[chain[$square,x_,y_List,z_,$angle],StandardForm|TraditionalForm] /;OddQ[Length[y]+2]:=SquareAnglechainBox[ToBoxes[x],ToBoxes[y],ToBoxes[z]];
chain /: MakeBoxes[chain[$square,x_,y_List,z_,$square],StandardForm|TraditionalForm] /;EvenQ[Length[y]+2]:=SquareSquarechainBox[ToBoxes[x],ToBoxes[y],ToBoxes[z]];


(* ::Subsection::Closed:: *)
(*Chain*)


(*Chain with negative momentum*)
Chain[a_,b_,{c___,pm[d_],e___},f_,g_] := -Chain[a,b,{c,d,e},f,g];

Chain[$angle,i_,{},j_,$angle] := SpinorAngleBracket[i,j];
Chain[$square,i_,{},j_,$square] := SpinorSquareBracket[i,j];

(*Display of the chains*)
AngleSquareChainBox[beginning_,args_,end_]:=TemplateBox[{beginning,end,Sequence@@args},"AngleSquareChain",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",#1,"|",TemplateSlotSequence[3,"\[ThinSpace]"],"|",#2,"]"}]&),
InterpretationFunction->(RowBox[{"Chain","[","$angle",",",#1,",","{",TemplateSlotSequence[3,","],"}",",",#2,",","$square","]"}]&)
];
SquareAngleChainBox[beginning_,args_,end_]:=TemplateBox[{beginning,end,Sequence@@args},"SquareAngleChain",
DisplayFunction->(RowBox[{"[",#1,"|",TemplateSlotSequence[3,"\[ThinSpace]"],"|",#2,"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"Chain","[","$square",",",#1,",","{",TemplateSlotSequence[3,","],"}",",",#2,",","$angle","]"}]&)
];
AngleAngleChainBox[beginning_,args_,end_]:=TemplateBox[{beginning,end,Sequence@@args},"AngleAngleChain",
DisplayFunction->(RowBox[{"\[LeftAngleBracket]",#1,"|",TemplateSlotSequence[3,"\[ThinSpace]"],"|",#2,"\[RightAngleBracket]"}]&),
InterpretationFunction->(RowBox[{"Chain","[","$angle",",",#1,",","{",TemplateSlotSequence[3,","],"}",",",#2,",","$angle","]"}]&)
];
SquareSquareChainBox[beginning_,args_,end_]:=TemplateBox[{beginning,end,Sequence@@args},"SquareSquareChain",
DisplayFunction->(RowBox[{"[",#1,"|",TemplateSlotSequence[3,"\[ThinSpace]"],"|",#2,"]"}]&),
  InterpretationFunction->(RowBox[{"Chain","[","$square",",",#1,",","{",TemplateSlotSequence[3,","],"}",",",#2,",","$square","]"}]&)
];
Chain /: MakeBoxes[Chain[$angle,x_,y_List,z_,$square],StandardForm|TraditionalForm] /;OddQ[Length[y]+2]:=AngleSquareChainBox[ToBoxes[x],ToBoxes/@y,ToBoxes[z]];
Chain /: MakeBoxes[Chain[$angle,x_,y_List,z_,$angle],StandardForm|TraditionalForm] /;EvenQ[Length[y]+2]:=AngleAngleChainBox[ToBoxes[x],ToBoxes/@y,ToBoxes[z]];
Chain /: MakeBoxes[Chain[$square,x_,y_List,z_,$angle],StandardForm|TraditionalForm] /;OddQ[Length[y]+2]:=SquareAngleChainBox[ToBoxes[x],ToBoxes/@y,ToBoxes[z]];
Chain /: MakeBoxes[Chain[$square,x_,y_List,z_,$square],StandardForm|TraditionalForm] /;EvenQ[Length[y]+2]:=SquareSquareChainBox[ToBoxes[x],ToBoxes/@y,ToBoxes[z]];
(*Vanishing chains*)
Chain[type_,x_,{x_,y___},z_,type2_]:=0;
Chain[type_,x_,{y___,z_},z_,type2_]:=0;
Chain[type_,x_,{y___,z_,z_,k___},l_,type2_]/; !MatchQ[z,_UnderBar|_PlusM] :=0;

Chain[a_,b_,{PlusM[b_,c___],e___},f_,g_]:=Chain[a,b,{PlusM[c],e},f,g];
Chain[a_,b_,{c___,PlusM[],e___},f_,g_]:=Chain[a,b,{c,e},f,g];
Chain[a_,b_,{c___,PlusM[h_],e___},f_,g_]:=Chain[a,b,{c,h,e},f,g];
Chain[a_,b_,{c_,c_},b_,a_]:=0;
Chain[$square,b_,{c_,c_},d_,$square]:= (S6Mom[c]+Extramass[c]Extramasstilde[c]) SpinorSquareBracket[b,d];
Chain[$angle,b_,{c_,c_},d_,$angle]:=(S6Mom[c]+Extramass[c]Extramasstilde[c]) SpinorAngleBracket[b,d];
Chain[b__,{a___,c_,c_,e___},d__] := (S6Mom[c]+Extramass[c]Extramasstilde[c])Chain[b,{a,e},d];

ReverseChain[ch:Chain[type1_,p1_,moms_List,p2_,type2_]]:=(-1)^(1+Length@moms)Chain[type2,p2,Reverse@moms,p1,type1];


SetAttributes[PlusM,Orderless];
PlusM[rest___,pm[x_],x_]:=PlusM[rest];

PlusM /: MakeBoxes[PlusM[y__],StandardForm|TraditionalForm] := TemplateBox[ToBoxes/@{y},"PlusM",
DisplayFunction->(RowBox[{"(",TemplateSlotSequence[1,"+"],")"}]&),
InterpretationFunction->(RowBox[{"PlusM","[",TemplateSlotSequence[1,","],"]"}]&)
];

PlusM /: S6Mom[PlusM[x___]] := S6Mom[x];
PlusM /: Extramass[PlusM[x___]] := Sum[Extramass[i],{i,{x}}];
PlusM /: Extramasstilde[PlusM[x___]] := Sum[Extramasstilde[i],{i,{x}}];
PlusM /: UnderBar[PlusM[x___]] := PlusM@@UnderBar/@{x};


(* ::Subsection::Closed:: *)
(*ChainN*)


ExpandPlusM[expr_] := expr //. Dispatch[Chain[a___,{b___,PlusM[c___],d___},e___]:>Sum[Chain[a,{b,cEl,d},e],{cEl,{c}}]];

(*Evaluation of ChainN*)
ChainN[$angle,a_,b_,c_,$angle]:=-Dot@@{SpinorUndotN[a][$lam][$down],
Sequence@@Table[Switch[b[[i]],
		PauliMatrixSpinor[_,_],
			b[[i]][If[OddQ[i],$up,$down]],
		_PlusM,
			Sum[
				Mom4DN[mom/.UnderBar->Identity][If[MemberQ[Momenta4DMassive,mom]||Head[mom/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$up,$down]]
			,{mom,List@@(b[[i]])}],
		_,
			Mom4DN[Replace[b[[i]],UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,b[[i]]]||Head[b[[i]]/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$up,$down]]
    ]
,{i,1,Length@b}],
SpinorUndotN[c][$lam][$up]};

ChainN[$square,a_,b_,c_,$square]:=-Dot@@{SpinorDotN[a][$lam][$up],
Sequence@@Table[Switch[b[[i]],
		PauliMatrixSpinor[_,_],
			b[[i]][If[OddQ[i],$down,$up]],
		_PlusM,
			Sum[
				Mom4DN[mom/.UnderBar->Identity][If[MemberQ[Momenta4DMassive,mom]||Head[mom/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$down,$up]]
			,{mom,List@@(b[[i]])}],
		_,
        Mom4DN[Replace[b[[i]],UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,b[[i]]]||Head[b[[i]]/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$down,$up]]
     ]
,{i,1,Length@b}],
SpinorDotN[c][$lam][$down]};

ChainN[$angle,a_,b_,c_,$square]:=Dot@@{SpinorUndotN[a][$lam][$down],
Sequence@@Table[Switch[b[[i]],
		PauliMatrixSpinor[_,_],
			b[[i]][If[OddQ[i],$up,$down]],
		_PlusM,
			Sum[
				Mom4DN[mom/.UnderBar->Identity][If[MemberQ[Momenta4DMassive,mom]||Head[mom/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$up,$down]]
			,{mom,List@@(b[[i]])}],
		_,
			Mom4DN[Replace[b[[i]],UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,b[[i]]]||Head[b[[i]]/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$up,$down]]
    ]
,{i,1,Length@b}],
SpinorDotN[c][$lam][$down]};

ChainN[$square,a_,b_,c_,$angle]:=Dot@@{SpinorDotN[a][$lam][$up],
Sequence@@Table[Switch[b[[i]],
		PauliMatrixSpinor[_,_],
			b[[i]][If[OddQ[i],$down,$up]],
		_PlusM,
			Sum[
				Mom4DN[mom/.UnderBar->Identity][If[MemberQ[Momenta4DMassive,mom]||Head[mom/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$down,$up]]
			,{mom,List@@(b[[i]])}],
		_,
			Mom4DN[Replace[b[[i]],UnderBar->Identity,1,Heads->True]][If[MemberQ[Momenta4DMassive,b[[i]]]||Head[b[[i]]/.pm->Identity]===UnderBar,Null,$flat]][If[OddQ[i],$down,$up]]
    ]
,{i,1,Length@b}],
SpinorUndotN[c][$lam][$up]};

(*Vanishing chains*)
ChainN[type_,x_,{x_,y___},z_,type2_]:=0;
ChainN[type_,x_,{y___,z_},z_,type2_]:=0;
ChainN[type_,x_,{y___,z_,z_,k___},l_,type2_]/; !MatchQ[z,_Underbar|_PlusM] :=0;



(* ::Subsection::Closed:: *)
(*ToChain*)


(*ToChain[exp_]:=Block[{localexp,MyPower,counter},
localexp=exp/.Power->MyPower;
(*First we need to group things properly in order to apply the contraction properties, since Mathematica is not able to recognize things inside powers as contractible with things outside that power...*)
If[FreeQ[localexp,MyPower]==False,
localexp=localexp/.{A_*MyPower[B_,x_]/;FreeQ[A,MyPower]:>MyPower[A,1]*MyPower[B,x]};
localexp=localexp//.MyPower[A_,x_]*MyPower[B_,x_]:>MyPower[A*B,x];
localexp=localexp/.MyPower[A_,x_?Negative]:>1/MyPower[A,-x];
localexp=localexp//.{MyPower[x_*A_,n_]*MyPower[B_,m_]/;n>m&&FreeQ[B,x]:>MyPower[x*A,n-m]MyPower[x*A*B,m]};
localexp=localexp//.{MyPower[x_,n_]*MyPower[B_,m_]/;n>m&&FreeQ[B,x]:>MyPower[x,n-m]MyPower[x*B,m]};
];

(*Now introduce the Chains and let them contract*)
localexp=localexp/.{SpinorAngleBracket[x_,y_]:>chain[$angle,x,{},y,$angle],SpinorSquareBracket[x_,y_]:>chain[$square,x,{},y,$square]};

(*Now replace the chain with Chain, where the first has all the contraction properties whereas the second does not and thus expressions stay as they are*)

localexp=localexp/.{chain->Chain};

(*Back to mathematica's Power*)
localexp=localexp/.{MyPower->Power};

(*Replace the single chains back with angle and square brackets*)
localexp=localexp/.{Chain[$angle,x_,{},y_,$angle]:>SpinorAngleBracket[x,y],Chain[$square,x_,{},y_,$square]:>SpinorSquareBracket[x,y]};

Return[localexp];
];*)


(* ::Subsection::Closed:: *)
(*ToChain*)


ToChain[exp_]:=Block[{localexp,MyPower,counter},
localexp=exp/.Power->MyPower;
(*First we need to group things properly in order to apply the contraction properties, since Mathematica is not able to recognize things inside powers as contractible with things outside that power...*)
If[FreeQ[localexp,MyPower]==False,
localexp=localexp/.{A_*MyPower[B_,x_]/;FreeQ[A,MyPower]:>MyPower[A,1]*MyPower[B,x]};
localexp=localexp/.MyPower[A_,x_?Negative]:>MyPower[MyPower[A,-x],-1];
localexp=localexp//.MyPower[A_,x_]*MyPower[B_,x_]:>MyPower[A*B,x];
localexp=localexp//.{MyPower[x_*A_,n_]*MyPower[B_,m_]/;n>m&&FreeQ[B,x]:>MyPower[x*A,n-m]MyPower[x*A*B,m]};
localexp=localexp//.{MyPower[x_,n_]*MyPower[B_,m_]/;n>m&&FreeQ[B,x]:>MyPower[x,n-m]MyPower[x*B,m]};
];

(*Now introduce the Chains and let them contract*)
localexp=localexp/.{SpinorAngleBracket[x_,y_]:>chain[$angle,x,{},y,$angle],SpinorSquareBracket[x_,y_]:>chain[$square,x,{},y,$square],Chain->chain};

(*Back to mathematica's Power*)
localexp=localexp/.{MyPower->Power};

(*Now replace the chain with Chain, where the first has all the contraction properties whereas the second does not and thus expressions stay as they are*)
localexp=localexp/.{chain->Chain};

(*Replace the single chains back with angle and square brackets*)
localexp=localexp/.{Chain[$angle,x_,{},y_,$angle]:>SpinorAngleBracket[x,y],Chain[$square,x_,{},y_,$square]:>SpinorSquareBracket[x,y]};

Return[localexp];
];


(* ::Subsection::Closed:: *)
(*ChainToSpinor*)


allFreeQ[expr_,form_List]:=And@@(FreeQ[expr,#]&/@form);

ChainToSpinor[exp_]:=Module[{localexp},
localexp=exp;
localexp=localexp/.{Chain[$angle,a_,b_,c_,$square]/;allFreeQ[b,Momenta4DMassive]:>Spinoranglebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,2,Length[b],2}]Spinorsquarebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,1,Length[b]-1,2}],
Chain[$square,a_,b_,c_,$angle]/;allFreeQ[b,Momenta4DMassive]:>Spinorsquarebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,1,Length[b]-1,2}]Spinoranglebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,2,Length[b],2}],
Chain[$angle,a_,b_,c_,$angle]/;allFreeQ[b,Momenta4DMassive]:>Spinoranglebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,2,Length[b]-1,2}]Spinoranglebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,1,Length[b],2}],
Chain[$square,a_,b_,c_,$square]/;allFreeQ[b,Momenta4DMassive]:>Spinorsquarebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,1,Length[b],2}]Spinorsquarebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,2,Length[b]-1,2}],
chain[$angle,a_,b_,c_,$square]/;allFreeQ[b,Momenta4DMassive]:>Spinoranglebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,2,Length[b],2}]Spinorsquarebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,1,Length[b]-1,2}],
chain[$angle,a_,b_,c_,$angle]/;allFreeQ[b,Momenta4DMassive]:>Spinoranglebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,2,Length[b]-1,2}]Spinoranglebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,1,Length[b],2}],
chain[$square,a_,b_,c_,$square]/;allFreeQ[b,Momenta4DMassive]:>Spinorsquarebracket[a,b[[1]]]Product[Spinoranglebracket[b[[i]],b[[i+1]]],{i,1,Length[b],2}]Spinorsquarebracket[Last[b],c]Product[Spinorsquarebracket[b[[i]],b[[i+1]]],{i,2,Length[b]-1,2}]};
Return[localexp];
];


(* ::Subsection::Closed:: *)
(*mp*)


mpBox[x_,y_]:=TemplateBox[{x,y},"ScalarProduct",
DisplayFunction->(RowBox[{"(",#1,"\[CenterDot]",#2,")"}]&),
InterpretationFunction->(RowBox[{"mp","[",#1,",",#2,"]"}]&)];
mp /: MakeBoxes[mp[x_,y_],StandardForm|TraditionalForm]:=mpBox[ToBoxes[x],ToBoxes[y]];
SetAttributes[mp,{Orderless,Protected}];


(* ::Subsection::Closed:: *)
(*eps*)


epsBox[a_,b_,c_,d_]:=TemplateBox[{a,b,c,d},"eps",
DisplayFunction->(RowBox[{"\[Epsilon]","[",#1,",",#2,",",#3,",",#4,"]"}]&),
InterpretationFunction->(RowBox[{"eps","[",#1,",",#2,",",#3,",",#4,"]"}]&)
];
eps /: MakeBoxes[eps[a_,b_,c_,d_],StandardForm|TraditionalForm]:=epsBox[ToBoxes[a],ToBoxes[b],ToBoxes[c],ToBoxes[d]];

(*Contraction with twice the same vector vanishes*)
eps[x___,y_,z___,y_,k___]:=0;

SetAttributes[eps,Protected];


(* ::Subsection::Closed:: *)
(*TrG*)


TrG[x_List]/;OddQ[Length[x]]:=0;
TrG[{}]:=4;
TrG[x_List]:=Sum[(-1)^i*mp[x[[1]],x[[i]]]TrG[Delete[x,{{1},{i}}]],{i,2,Length[x]}];


(* ::Subsection::Closed:: *)
(*TrG5*)


TrG5[x_List]/;Length[x]<4:=0;
TrG5[x_List]/;OddQ[Length[x]]:=0;
TrG5[x_List]:=mp[x[[-3]],x[[-2]]]*TrG5[Delete[x,{{-3},{-2}}]]+mp[x[[-2]],x[[-1]]]*TrG5[x[[;;-3]]]-mp[x[[-3]],x[[-1]]]*TrG5[Delete[x,{{-3},{-1}}]]-I*Sum[(-1)^i*eps[x[[-i]],x[[-3]],x[[-2]],x[[-1]]]*TrG[Delete[x[[;;-4]],{-(i-3)}]],{i,4,Length[x]}];


(* ::Subsection::Closed:: *)
(*ToTrace*)


Options[ToTrace]={KillEpsilon->False}
{KillEpsilon->False}
ToTrace[exp_,OptionsPattern[]]:=Block[{eps,localexp},
Which[OptionValue[KillEpsilon],
eps[x__]:=0,
OptionValue[KillEpsilon]==False,
Null,
True,
Print["Undefined value for the option KillEpsilon. Only True or False are allowed. Proceed assuming default value False."];
];
(*Convert the chains to traces*)
localexp=exp/.Chain[$angle,a_,b_List,a_,$square]:>(TrG[Join[{a},b]]-TrG5[Join[{a},b]])/2;

(*Return output*)
Return[localexp];
];

SetAttributes[ToTrace,Protected];


(* ::Subsection::Closed:: *)
(*SumContracted*)


Options[SumContracted]={StateSumParallel->False,DisplayKernel->False};

SumContracted[x__,OptionsPattern[]][exp_]:=Block[{out,vars,local},
vars=Sequence@@Table[{i,2},{i,{x}}];
If[TrueQ[OptionValue[StateSumParallel]],
(*Parallel sum over states*)
If[$KernelCount==0,
LaunchKernels[];
];
(*Prepare objects to be evaluated in parallel*)
If[TrueQ[OptionValue[DisplayKernel]],
(*Display which kernel computed what in the final result*)
local=Table[ParallelSubmit[{x},Labeled[Framed[exp],$KernelID]],Evaluate[vars]]//Flatten,
(*Without displaying the Kernel*)
local=Table[ParallelSubmit[{x},exp],Evaluate[vars]]//Flatten;
];

(*Evaluate*)
out=WaitAll[local];
out=Plus@@out,
out=Sum[exp,Evaluate[vars]];
];
Return[out];
];


(* ::Subsection::Closed:: *)
(*CompleteToMassive*)


(*CompleteToMassive[exp_,reps_:FixedSpinors]:=Block[{$crep2,locvar,SubCounter2,localexp,locreps,moms,SpinorUndot,SpinorDot,out,MyPower,numer,firstrep,MomLeft,MomRight,MomMassive},
(*Auxiliary function needed for the replacements later on*)
$crep2=2;
SubCounter2:=locvar*IntegerPart[$crep2++/2];

(*Pick only the relevant replacements, i.e. those concerning the \[Mu]s*)
locreps=Select[reps,(!FreeQ[#,SpinorUndotPure[_][$mu]->SpinorUndotPure[_][$lam]]||!FreeQ[#,SpinorDotPure[_][$mu]->SpinorDotPure[_][$lam]])&];

(*Extract the momenta whose \[Mu]s need to be restored*)
moms=Table[i[[1]],{i,locreps}];
moms=DeleteDuplicates[moms//.{SpinorUndotPure[x_][_]:>x,SpinorDotPure[x_][_]:>x}];

(*Convert the whole expression into single spinors, making also use of the auxiliary function SubCounter defined along with SpinorReplace. Before doing this we split powers into single objects to get a fully correct replacement.*)
localexp=exp//.Power[A_,n_?Positive]:>Product[MyPower[A,i],{i,n}];

localexp=localexp/.{
SpinorAngleBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]:> SpinorUndot[x][$mu][SubCounter2][Null]*SpinorUndot[y][$mu][Null][SubCounter2],SpinorSquareBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]:>SpinorDot[x][$mu][Null][SubCounter2]*SpinorDot[y][$mu][SubCounter2][Null],
SpinorAngleBracket[x_,y_]/;TrueQ[Head[y]==OverBar]:> SpinorUndot[x][$lam][SubCounter2][Null]*SpinorUndot[y][$mu][Null][SubCounter2],SpinorSquareBracket[x_,y_]/;TrueQ[Head[y]==OverBar]:>SpinorDot[x][$lam][Null][SubCounter2]*SpinorDot[y][$mu][SubCounter2][Null],SpinorAngleBracket[x_,y_]:> SpinorUndot[x][$lam][SubCounter2][Null]*SpinorUndot[y][$lam][Null][SubCounter2],SpinorSquareBracket[x_,y_]:>SpinorDot[x][$lam][Null][SubCounter2]*SpinorDot[y][$lam][SubCounter2][Null]};
localexp=localexp/.{SpinorUndot[OverBar[x_]][$mu][y_][z_]:>SpinorUndot[x][$mu][y][z],SpinorDot[OverBar[x_]][$mu][y_][z_]:>SpinorDot[x][$mu][y][z]};
(*Now we can safely remove MyPower*)
localexp=localexp//.MyPower[A_,x_]:>A;

(*Now replace the massless momenta |p>[p| with a massive momentum plus the \[Mu] piece. The massive momentum cannot be written any more in terms of spinors so will appear only inside a chain object.*)

(*localexp=localexp//Expand;*)
(*Now we define a bunch of objects. These will allow us to do the replecements in a more efficient way*)
MomLeft /: Times[MomLeft[mom_][x_][y_],MomRight[mom_][z_][k_]]:=MomMassive[{mom},$angle,$square][x,z][y,k]+extramass[mom]*extramasstilde[mom]*SpinorUndot[mom][$mu][x][y]*SpinorDot[mom][$mu][z][k]/(SpinorAngleBracket[mom,overbar[mom]]SpinorSquareBracket[mom,overbar[mom]]);
(*Now we define the contraction properties of MomMassive with itself*)
MomMassive /: Times[MomMassive[x_List,$angle,$square][a_,adot_][Null,bdot_],MomMassive[y_List,$angle,$square][Null,cdot_][a_,ddot_]]:=MomMassive[Join[x,y],$square,$square][adot,cdot][bdot,ddot];
MomMassive /: Times[MomMassive[x_List,$angle,$square][a_,Null][b_,bdot_],MomMassive[y_List,$angle,$square][c_,bdot_][d_,Null]]:=MomMassive[Join[x,y],$angle,$angle][a,c][b,d];
MomMassive /: Times[MomMassive[x_List,$angle,$angle][a1_,a2_][b1_,Null],MomMassive[y_List,$angle,type_][Null,cdot_][a2_,ddot_]]:=MomMassive[Join[x,y],$angle,type][a1,cdot][b1,ddot];
MomMassive /: Times[MomMassive[x_List,$angle,$angle][Null,a2_][c_,b2_],MomMassive[y_List,type_,$angle][cdot_,c_][ddot_,Null]]:=MomMassive[Join[y,x],type,$angle][cdot,a2][ddot,b2];
MomMassive /: Times[MomMassive[x_List,$square,$square][a1_,a2_][b1_,Null],MomMassive[y_List,$square,type_][Null,cdot_][a2_,ddot_]]:=-MomMassive[Join[x,y],$square,type][a1,cdot][b1,ddot];
MomMassive /: Times[MomMassive[x_List,$square,$square][Null,a2_][c_,b2_],MomMassive[y_List,type_,$square][cdot_,c_][ddot_,Null]]:=-MomMassive[Join[y,x],type,$square][cdot,a2][ddot,b2];

Do[
localexp=localexp/.{SpinorUndot[i][$lam][x_][y_]:>MomLeft[i][x][y],SpinorDot[i][$lam][x_][y_]:>MomRight[i][x][y]};
localexp=localexp/.{MomLeft[i][x_][y_]:>SpinorUndot[i][$lam][x][y],MomRight[i][x_][y_]:>SpinorDot[i][$lam][x][y]};
,{i,moms}];


(*Now we replace the resulting \[Mu]s with the chosen \[Lambda]s again*)
localexp=SpinorReplaceSequential[localexp,locreps];

(*Finally we define the contraction properties of MomMassive with angle and square spinors*)

MomMassive /: Times[MomMassive[x_List,$angle,$angle][a1_,a2_][Null,Null],SpinorUndot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][Null][a2_]]:=-chain[$angle,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$angle,$angle][Null,Null][a1_,a2_],SpinorUndot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][a2_][Null]]:=-chain[$angle,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$angle,$angle][a1_,Null][Null,a2_],SpinorUndot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][a2_][Null]]:=chain[$angle,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$angle,$angle][Null,a2_][a1_,Null],SpinorUndot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][Null][a2_]]:=chain[$angle,lab,x,lab2,$angle];

MomMassive /: Times[MomMassive[x_List,$square,$square][a1_,a2_][Null,Null],SpinorDot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][Null][a2_]]:=-chain[$square,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$square,$square][Null,Null][a1_,a2_],SpinorDot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][a2_][Null]]:=-chain[$square,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$square,$square][a1_,Null][Null,a2_],SpinorDot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][a2_][Null]]:=chain[$square,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$square,$square][Null,a2_][a1_,Null],SpinorDot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][Null][a2_]]:=chain[$square,lab,x,lab2,$square];

MomMassive /: Times[MomMassive[x_List,$square,$angle][a1_,a2_][Null,Null],SpinorDot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][Null][a2_]]:=+chain[$square,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$square,$angle][Null,Null][a1_,a2_],SpinorDot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][a2_][Null]]:=+chain[$square,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$square,$angle][a1_,Null][Null,a2_],SpinorDot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][a2_][Null]]:=-chain[$square,lab,x,lab2,$angle];
MomMassive /: Times[MomMassive[x_List,$square,$angle][Null,a2_][a1_,Null],SpinorDot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][Null][a2_]]:=-chain[$square,lab,x,lab2,$angle];

MomMassive /: Times[MomMassive[x_List,$angle,$square][a1_,a2_][Null,Null],SpinorUndot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][Null][a2_]]:=+chain[$angle,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$angle,$square][Null,Null][a1_,a2_],SpinorUndot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][a2_][Null]]:=+chain[$angle,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$angle,$square][a1_,Null][Null,a2_],SpinorUndot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][a2_][Null]]:=-chain[$angle,lab,x,lab2,$square];
MomMassive /: Times[MomMassive[x_List,$angle,$square][Null,a2_][a1_,Null],SpinorUndot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][Null][a2_]]:=-chain[$angle,lab,x,lab2,$square];

(*Just to make sure all the properties have been applied*)
out=localexp//Expand;

(*Return the output*)
Return[out];
];*)


(* ::Subsection::Closed:: *)
(*CompleteToMassive*)


(*MomMassive is an auxiliary function for CompleteToMassive*)

(*Now we define the contraction properties of MomMassive with itself*)
(*MomMassive /: Times[MomMassive[l1_List,type1_,$angle][a1_,A_][a2_,Null],MomMassive[l2_List,type2_,$angle][b1_,Null][b2_,A_]]:=(-1)^(Length[l2]+1)MomMassive[Join[l1,Reverse[l2]],type1,type2][a1,b1][a2,b2];*)
(*MomMassive /: Times[MomMassive[l1_List,$angle,type1_][A_,a1_][Null,a2_],MomMassive[l2_List,$angle,type2_][Null,b1_][A_,b2_]]:=(-1)^(Length[l1]+1)MomMassive[Join[Reverse[l1],l2],type1,type2][a1,b1][a2,b2];*)
(*MomMassive /:  Times[MomMassive[l1_List,type1_,$angle][a1_,A_][a2_,Null],MomMassive[l2_List,$angle,type2_][Null,b1_][A_,b2_]]:=MomMassive[Join[l1,l2],type1,type2][a1,b1][a2,b2];*)
(*MomMassive /: Times[MomMassive[l1_List,$square,type1_][Null,a1_][A_,a2_],MomMassive[l2_List,$square,type2_][A_,b1_][Null,b2_]]:=(-1)^(Length[l1]+1)MomMassive[Join[Reverse[l1],l2],type1,type2][a1,b1][a2,b2];*)
(*MomMassive /: Times[MomMassive[l1_List,type1_,$square][a1_,Null][a2_,A_],MomMassive[l2_List,type2_,$square][b1_,A_][b2_,Null]]:=(-1)^(Length[l2]+1)MomMassive[Join[l1,Reverse[l2]],type1,type2][a1,b1][a2,b2];*)
(*MomMassive /: Times[MomMassive[l1_List,type1_,$square][a1_,Null][a2_,A_],MomMassive[l2_List,$square,type2_][A_,b1_][Null,b2_]]:=MomMassive[Join[l1,l2],type1,type2][a1,b1][a2,b2];*)
(*Finally we define the contraction properties of MomMassive with angle and square spinors*)

(*MomMassive /: Times[MomMassive[x_List,$angle,$angle][a1_,a2_][Null,Null],SpinorUndot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][Null][a2_]]:=-chain[$angle,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$angle][Null,Null][a1_,a2_],SpinorUndot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][a2_][Null]]:=-chain[$angle,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$angle][a1_,Null][Null,a2_],SpinorUndot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][a2_][Null]]:=chain[$angle,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$angle][Null,a2_][a1_,Null],SpinorUndot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][Null][a2_]]:=chain[$angle,lab,x,lab2,$angle];*)

(*MomMassive /: Times[MomMassive[x_List,$square,$square][a1_,a2_][Null,Null],SpinorDot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][Null][a2_]]:=-chain[$square,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$square][Null,Null][a1_,a2_],SpinorDot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][a2_][Null]]:=-chain[$square,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$square][a1_,Null][Null,a2_],SpinorDot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][a2_][Null]]:=chain[$square,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$square][Null,a2_][a1_,Null],SpinorDot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][Null][a2_]]:=chain[$square,lab,x,lab2,$square];*)

(*MomMassive /: Times[MomMassive[x_List,$square,$angle][a1_,a2_][Null,Null],SpinorDot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][Null][a2_]]:=+chain[$square,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$angle][Null,Null][a1_,a2_],SpinorDot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][a2_][Null]]:=+chain[$square,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$angle][a1_,Null][Null,a2_],SpinorDot[lab_][$lam][Null][a1_],SpinorUndot[lab2_][$lam][a2_][Null]]:=-chain[$square,lab,x,lab2,$angle];*)
(*MomMassive /: Times[MomMassive[x_List,$square,$angle][Null,a2_][a1_,Null],SpinorDot[lab_][$lam][a1_][Null],SpinorUndot[lab2_][$lam][Null][a2_]]:=-chain[$square,lab,x,lab2,$angle];*)

(*MomMassive /: Times[MomMassive[x_List,$angle,$square][a1_,a2_][Null,Null],SpinorUndot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][Null][a2_]]:=+chain[$angle,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$square][Null,Null][a1_,a2_],SpinorUndot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][a2_][Null]]:=+chain[$angle,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$square][a1_,Null][Null,a2_],SpinorUndot[lab_][$lam][Null][a1_],SpinorDot[lab2_][$lam][a2_][Null]]:=-chain[$angle,lab,x,lab2,$square];*)
(*MomMassive /: Times[MomMassive[x_List,$angle,$square][Null,a2_][a1_,Null],SpinorUndot[lab_][$lam][a1_][Null],SpinorDot[lab2_][$lam][Null][a2_]]:=-chain[$angle,lab,x,lab2,$square];*)

(*CompleteToMassive[exp_,reps_:FixedSpinors]:=Block[{$crep2,locvar,SubCounter2,localexp,locreps,moms,numtot,dentot,out,MyPower,SpinorUndot2,SpinorDot2,numer,firstrep,MomLeft,MomRight},*)

(*Auxiliary function needed for the replacements later on*)
(*$crep2=2;*)
(*SubCounter2:=locvar*IntegerPart[$crep2++/2];*)

(*Pick only the relevant replacements, i.e. those concerning the \[Mu]s*)
(*locreps=Select[reps,(!FreeQ[#,SpinorUndotPure[_][$mu]->SpinorUndotPure[_][$lam]]||!FreeQ[#,SpinorDotPure[_][$mu]->SpinorDotPure[_][$lam]])&];*)

(*Extract the momenta whose \[Mu]s need to be restored*)
(*moms=Table[i[[1]],{i,locreps}];*)
(*moms=DeleteDuplicates[moms//.{SpinorUndotPure[x_][_]:>x,SpinorDotPure[x_][_]:>x}];*)

(*Now we are going to separate numerator and denominator*)

(*localexp=Together[exp];*)
(*numtot={Numerator[localexp]//Expand//CompleteMandelstam}/.{Plus->List}//Flatten;*)
(*dentot=Denominator[localexp];*)

(*Convert the whole expression into single spinors, making also use of the auxiliary function SubCounter defined along with SpinorReplace. Before doing this we split powers into single objects to get a fully correct replacement.*)
(*numtot=numtot//.Power[A_,n_?Positive]:>Product[MyPower[A,i],{i,n}];*)

(*numtot=numtot/.{*)
(*SpinorAngleBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]:> SpinorUndot2[x][$mu][SubCounter2][Null]*SpinorUndot2[y][$mu][Null][SubCounter2],SpinorSquareBracket[x_,y_]/;TrueQ[Head[x]==OverBar&&Head[y]==OverBar]:>SpinorDot2[x][$mu][Null][SubCounter2]*SpinorDot2[y][$mu][SubCounter2][Null],*)
(*SpinorAngleBracket[x_,y_]/;TrueQ[Head[y]==OverBar]:> SpinorUndot2[x][$lam][SubCounter2][Null]*SpinorUndot2[y][$mu][Null][SubCounter2],SpinorSquareBracket[x_,y_]/;TrueQ[Head[y]==OverBar]:>SpinorDot2[x][$lam][Null][SubCounter2]*SpinorDot2[y][$mu][SubCounter2][Null],SpinorAngleBracket[x_,y_]:> SpinorUndot2[x][$lam][SubCounter2][Null]*SpinorUndot2[y][$lam][Null][SubCounter2],SpinorSquareBracket[x_,y_]:>SpinorDot2[x][$lam][Null][SubCounter2]*SpinorDot2[y][$lam][SubCounter2][Null]};*)
(*numtot=numtot/.{SpinorUndot2[OverBar[x_]][$mu][y_][z_]:>SpinorUndot2[x][$mu][y][z],SpinorDot2[OverBar[x_]][$mu][y_][z_]:>SpinorDot2[x][$mu][y][z]};*)
(*Now we can safely remove MyPower*)
(*numtot=numtot//.MyPower[A_,x_]:>A;*)

(*Now replace the massless momenta |p>[p| with a massive momentum plus the \[Mu] piece. The massive momentum cannot be written any more in terms of spinors so will appear only inside a chain object.*)

(*Now we define a bunch of objects. These will allow us to do the replecements in a more efficient way*)
(*MomLeft /: Times[MomLeft[mom_][x_][y_],MomRight[mom_][z_][k_]]:=MomMassive[{mom},$angle,$square][x,z][y,k]+extramass[mom]*extramasstilde[mom]*SpinorUndot[mom][$mu][x][y]*SpinorDot[mom][$mu][z][k]/(SpinorAngleBracket[mom,overbar[mom]]SpinorSquareBracket[mom,overbar[mom]]);*)

(*Do[*)
(*numtot=numtot/.{SpinorUndot2[i][$lam][x_][y_]:>MomLeft[i][x][y],SpinorDot2[i][$lam][x_][y_]:>MomRight[i][x][y]};*)
(*numtot=numtot/.{MomLeft[i][x_][y_]:>SpinorUndot2[i][$lam][x][y],MomRight[i][x_][y_]:>SpinorDot2[i][$lam][x][y]};*)
(*,{i,moms}];*)


(*Now we replace the four-dimensional Mandelstam invariants with six-dimensional ones*)

(*numtot=numtot/.{S4[x1_,x2_]:>S[x1,x2]+extramass[x2]*extramasstilde[x1]+extramass[x1]*extramasstilde[x2]-(extramass[x1]*extramasstilde[x1] Spinoranglebracket[x2,\!\(\*OverscriptBox[\(x1\), \(_\)]\)] Spinorsquarebracket[x2,\!\(\*OverscriptBox[\(x1\), \(_\)]\)])/(Spinoranglebracket[x1,\!\(\*OverscriptBox[\(x1\), \(_\)]\)] Spinorsquarebracket[x1,\!\(\*OverscriptBox[\(x1\), \(_\)]\)])-(extramass[x2]*extramasstilde[x2] Spinoranglebracket[x1,\!\(\*OverscriptBox[\(x2\), \(_\)]\)] Spinorsquarebracket[x1,\!\(\*OverscriptBox[\(x2\), \(_\)]\)])/(Spinoranglebracket[x2,\!\(\*OverscriptBox[\(x2\), \(_\)]\)] Spinorsquarebracket[x2,\!\(\*OverscriptBox[\(x2\), \(_\)]\)])+(extramass[x1]*extramasstilde[x1]extramass[x2]*extramasstilde[x2] Spinoranglebracket[\!\(\*OverscriptBox[\(x1\), \(_\)]\),\!\(\*OverscriptBox[\(x2\), \(_\)]\)] Spinorsquarebracket[\!\(\*OverscriptBox[\(x1\), \(_\)]\),\!\(\*OverscriptBox[\(x2\), \(_\)]\)])/(Spinoranglebracket[x1,\!\(\*OverscriptBox[\(x1\), \(_\)]\)] Spinoranglebracket[x2,\!\(\*OverscriptBox[\(x2\), \(_\)]\)] Spinorsquarebracket[x1,\!\(\*OverscriptBox[\(x1\), \(_\)]\)] Spinorsquarebracket[x2,\!\(\*OverscriptBox[\(x2\), \(_\)]\)])};*)

(*Now we replace the resulting \[Mu]s with the chosen \[Lambda]s again*)
(*numtot=SpinorReplaceSequential[numtot,reps];*)

(*Just to make sure all the properties have been applied*)
(*numtot=numtot/.{SpinorUndot2->SpinorUndot,SpinorDot2->SpinorDot}//Expand;*)

(*Make numtot a sum again*)
(*If[Length[numtot]==1,*)
(*numtot=First[numtot],*)
(*numtot=Plus@@numtot;*)
(*];*)

(*Return the output*)
(*Return[numtot/dentot];*)
(*];*)


(* ::Subsection::Closed:: *)
(*Relabel*)


Options[Relabel]={LabelRep->{},ScalarProduct->True};

Relabel[exp_,OptionsPattern[]]:=Module[{localexp,funrep,labrep,labrepneg,sp},
(*Divide the replacements in two cathegories: the ines acting on labesl and the ones acting on the entire functions S[__] or mp[_,_]*)

funrep=Select[OptionValue[LabelRep],(!FreeQ[#,S[__]]||!FreeQ[#,mp[_,_]]||!FreeQ[#,S4[__]])&];

(*Label replacements are harder to deal with since they may be simple numbers,so they need to be separated and treated more carefully. In order to do these replacemenst properly we transfor the labels of all the functions into strings. This allows us to distinguish for example the label 2 from the number 2 in 2*l1 where l1 is the label. This has the drawback that we also need to duplicate the replacements by introducing a set of negative labels as well. This seems to be the safest way of doing things.*)

labrep=Select[OptionValue[LabelRep],(FreeQ[#,S[__]]&&FreeQ[#,mp[_,_]]&&FreeQ[#,S4[__]])&];
labrepneg=Table[-i[[1]]->-i[[2]],{i,labrep}];
labrep=Table[ToString[i[[1]]]->ToString[i[[2]]],{i,labrep}];
labrepneg=Table[ToString[i[[1]]]->ToString[i[[2]]],{i,labrepneg}];

(*Convert all the labels to strings*)
localexp=exp/.{S[x__]:>S[Sequence@@ToString/@{x}],S4[x__]:>S4[Sequence@@ToString/@{x}],mp[x__]:>mp[Sequence@@ToString/@{x}]};
(*Apply the label replacements*)
(*localexp=localexp/.f_[x__]\[RuleDelayed]f[StringReplace[#,Join[labrep,labrepneg]]&/@{x}];*)
localexp=localexp/.Join[labrep,labrepneg];
(*Convert labels back to symbols*)
localexp=localexp/.{S[x__]:>S[Sequence@@ToExpression/@{x}],S4[x__]:>S4[Sequence@@ToExpression/@{x}],mp[x__]:>mp[Sequence@@ToExpression/@{x}]};
(*Apply the functional replacements*)
localexp=localexp/.funrep;

(*Finally convert the scalar product if a new one has been defined*)

If[!TrueQ[OptionValue[ScalarProduct]],
sp=OptionValue[ScalarProduct];
localexp=localexp/.{mp[x_,y_]:>sp[x,y],S4[x_,y_]:>sp[x+y,x+y],S[x_,y_]:>sp[x+y,x+y]};
];

Return[localexp];
];


(* ::Subsection::Closed:: *)
(*SpinorAngleBracketN*)


SpinorAngleBracketN[x_,y_]:=SpinorUndotN[x][$lam][$up] . SpinorUndotN[y][$lam][$down];
SpinorAngleBracketN[x_,OverBar[y_]]:=SpinorUndotN[x][$lam][$up] . SpinorUndotN[y][$mu][$down];
SpinorAngleBracketN[OverBar[x_],y_]:=SpinorUndotN[x][$mu][$up] . SpinorUndotN[y][$lam][$down];
SpinorAngleBracketN[OverBar[x_],OverBar[y_]]:=SpinorUndotN[x][$mu][$up] . SpinorUndotN[y][$mu][$down];


(* ::Subsection::Closed:: *)
(*SpinorSquareBracketN*)


SpinorSquareBracketN[x_,y_]:=SpinorDotN[x][$lam][$down] . SpinorDotN[y][$lam][$up];
SpinorSquareBracketN[x_,OverBar[y_]]:=SpinorDotN[x][$lam][$down] . SpinorDotN[y][$mu][$up];
SpinorSquareBracketN[OverBar[x_],y_]:=SpinorDotN[x][$mu][$down] . SpinorDotN[y][$lam][$up];
SpinorSquareBracketN[OverBar[x_],OverBar[y_]]:=SpinorDotN[x][$mu][$down] . SpinorDotN[y][$mu][$up];


(* ::Subsection::Closed:: *)
(*GenSpinors*)


Options[GenSpinors]={Dimension->6,DisplaySpinors->False,Parametric->False,ParameterName->$par,ParameterRange->200,
	Modulus->NextPrime[10000,-1],SpecialKinematics->Rational,Seed->False,RandomSpinors->False,FourD->{},FixedReference->False};
SetAttributes[$par,Protected];
SetAttributes[FiniteField,Protected];

GenSpinors[labels_List,opts:OptionsPattern[]]:=Module[
{lambda,lambdaout,num,ra,rs,la,ls,\[Xi],\[Eta],\[Xi]t,\[Eta]t,mass,masst,mur,mul,mutr,mutl,parn,usedpar,system,sol,dim,prime,FreeIndcs,FixIndcs,solvVars},
	
	
	
	(*First we check whether the give combination of options work together.*)
	If[OptionValue[FixedReference]=!=False,
	   If[TrueQ[OptionValue[Parametric]],
	      Message[GenSpinors::FixedReferenceParametric];
	      Abort[];
	   ];
	   If[!MatchQ[OptionValue[FixedReference],Rule[_,List[_?NumericQ,_?NumericQ,_?NumericQ,_?NumericQ]]],
	      Message[GenSpinors::FixedReferenceBadSyntax];
	      Abort[];
	   ];
	];
	
	
	
	(*Next we have to clear the definitions of SpinorUndotN and SpinorDotN*)
	ClearKinematics[labels];
	
	If[OptionValue[FixedReference]=!=False,
	   SubValues[SpinorUndotN]=DeleteCases[SubValues[SpinorUndotN],_?(!FreeQ[First[#],With[{j=OptionValue[FixedReference][[1]]},HoldPattern[SpinorUndotN[j]]]]&)];
	   SubValues[SpinorDotN]=DeleteCases[SubValues[SpinorDotN],_?(!FreeQ[First[#],With[{j=OptionValue[FixedReference][[1]]},HoldPattern[SpinorDotN[j]]]]&)];
	   SubValues[ExtramassN]=DeleteCases[SubValues[ExtramassN],_?(!FreeQ[First[#],With[{j=OptionValue[FixedReference][[1]]},HoldPattern[ExtramassN[j]]]]&)];
	   SubValues[ExtramasstildeN]=DeleteCases[SubValues[ExtramasstildeN],_?(!FreeQ[First[#],With[{j=OptionValue[FixedReference][[1]]},HoldPattern[ExtramasstildeN[j]]]]&)];
	];
	
		If[OptionValue[Dimension]==6 && (Length[labels] - Length@Intersection[labels,OptionValue[FourD]]) < 2,
			dim = 4;
		,
			dim = OptionValue[Dimension];
		];
		
	(*Fix the seed of the random generator if required using the Seed option*)
	If[TrueQ[Head[OptionValue[Seed]]==Integer],
		SeedRandom[OptionValue[Seed]];
	];
	
	(*Four dimensions*)
	Switch[dim,
		4,
			genSpinors4D[labels,FilterRules[{opts,Options[GenSpinors]},Options[genSpinors4D]]];,
		6,
			genSpinors6D[labels,FilterRules[{opts,Options[GenSpinors]},Options[genSpinors6D]]];,
		_,
		Message[GenSpinors::UnsupportedDimension];
	];
];


(* ::Subsubsection:: *)
(*GenSpinors 4D*)


Options[genSpinors4D]={DisplaySpinors->False,Parametric->False,ParameterName->$par,ParameterRange->200,Modulus->NextPrime[10000,-1],
	SpecialKinematics->Rational,Seed->False,RandomSpinors->False,FixedReference->False};

genSpinors4D[labels_List,opts:OptionsPattern[]]:=Module[
{lambda,lambdaout,num,ra,rs,la,ls,\[Xi],\[Eta],\[Xi]t,\[Eta]t,mass,masst,mur,mul,mutr,mutl,parn,usedpar,system,sol,dim,prime,FreeIndcs,FixIndcs,solvVars},

	If[OptionValue[SpecialKinematics]===FiniteField,
	   If[PrimeQ[OptionValue[Modulus]],
	      prime = Abs[OptionValue[Modulus]];
		,
	      Message[GenSpinors::NotPrime];
	      Abort[];
	   ];
	   If[TrueQ[OptionValue[Parametric]],
	      Message[GenSpinors::FiniteFieldParametric];
	      Abort[];
	   ];
	];
	
	(*Definition of the spinors in terms of the components. In the six-dimensionla case there will be 2n of tese spinors, where the first n refer to \[Lambda] and the secodn n to \[Lambda]' which is redefinition of the \[Mu] encoding also the masses, see paper to be written.*)
	ra[i_]:={\[Xi][i],\[Eta][i]}(*|i>*);
	rs[i_]:={\[Xi]t[i],\[Eta]t[i]}(*|i]*);
	la[i_]:={\[Eta][i],-\[Xi][i]}(*<i|*);
	ls[i_]:={-\[Eta]t[i],\[Xi]t[i]}(*[i|*);

	num=Length[labels];
	
	(*Fix n+6 components to arbitrary numbers, distinguishing between rational, finite field, and real kinematics*)
         Switch[OptionValue[SpecialKinematics],
            Rational,
               Do[\[Xi][i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,num}];
               Do[\[Xi]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}];
               	\[Eta]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}];
               	\[Eta][i+2]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,2}],
            FiniteField,
               Do[\[Xi][i]=RandomInteger[{1,prime}],{i,num}];
               Do[\[Xi]t[i]=RandomInteger[{1,prime}];
               	\[Eta]t[i]=RandomInteger[{1,prime}];
               	\[Eta][i+2]=RandomInteger[{1,prime}],{i,2}],
            _,
               Do[\[Xi][i]=RandomReal[OptionValue[ParameterRange]],{i,num}];
               Do[\[Xi]t[i]=RandomReal[OptionValue[ParameterRange]];
               	\[Eta]t[i]=RandomReal[OptionValue[ParameterRange]];
               	\[Eta][i+2]=RandomReal[OptionValue[ParameterRange]],{i,2}];
         ];

         (*We check if we want numeric values or parametric values*)
         If[!TrueQ[OptionValue[Parametric]],
            (*Generate numerical values for all the remaining components except 4, which will solve momentum conservation*)
            Switch[OptionValue[SpecialKinematics],
               Rational,
                  Do[\[Eta][i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,5,num}];
                  Do[\[Xi]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}];
                  	\[Eta]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,4,num}],
               FiniteField,
                  Do[\[Eta][i]=RandomInteger[{1,prime}],{i,5,num}];
                  Do[\[Xi]t[i]=RandomInteger[{1,prime}];
                  	\[Eta]t[i]=RandomInteger[{1,prime}],{i,4,num}],
               _,
                  Do[\[Eta][i]=RandomReal[OptionValue[ParameterRange]],{i,5,num}];
                  Do[\[Xi]t[i]=RandomReal[OptionValue[ParameterRange]];
                  	\[Eta]t[i]=RandomReal[OptionValue[ParameterRange]],{i,4,num}];
            ],
            (*Assign a parameter name to the remaining components*)
            parn=OptionValue[ParameterName];
            Do[\[Xi]t[i]=parn[i-3],{i,4,num}];
            usedpar=num-3;
            Do[\[Eta]t[i]=parn[i+usedpar-3],{i,4,num}];
            usedpar=usedpar+(num-3);
            Do[\[Eta][i]=parn[i+usedpar-4],{i,5,num}];
         ];

         If[OptionValue[FixedReference] =!= False,
            (*The +1 definitions will be used later to set the global variables (as they are guaranteed not to be overwritten by having 4D momenta)*)
            \[Xi][num+1]=OptionValue[FixedReference][[2,1]];
            \[Eta][num+1]=OptionValue[FixedReference][[2,2]];
            \[Xi]t[num+1]=OptionValue[FixedReference][[2,3]];
            \[Eta]t[num+1]=OptionValue[FixedReference][[2,4]];
         ];
         
         (*Solve momentum conservation in other 4 components, or set to arbitrary numbers if RandomSpinors is set to true*)
         If[TrueQ[OptionValue[RandomSpinors]],
            {\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}=Table[RandomInteger[{1,OptionValue[ParameterRange]}],{i,4}],
            If[OptionValue[SpecialKinematics]===FiniteField,
               sol=Solve[Sum[\[Xi][i]\[Xi]t[i],{i,num}]==0&&Sum[\[Xi][i]\[Eta]t[i],{i,num}]==0&&Sum[\[Eta][i]\[Xi]t[i],{i,num}]==0&&Sum[\[Eta][i]\[Eta]t[i],{i,num}]==0,
                   {\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]},Modulus->prime];
               {\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}=Mod[{\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}/.First[sol],prime];
           ,
               sol=Solve[Sum[\[Xi][i]\[Xi]t[i],{i,num}]==0&&Sum[\[Xi][i]\[Eta]t[i],{i,num}]==0&&Sum[\[Eta][i]\[Xi]t[i],{i,num}]==0&&Sum[\[Eta][i]\[Eta]t[i],{i,num}]==0,
                   {\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}];
               {\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}={\[Eta][1],\[Eta][2],\[Xi]t[3],\[Eta]t[3]}/.First[sol];
            ];
         ];

         (*Relate generated kinematics to the spinors:*)
         If[OptionValue[SpecialKinematics]===FiniteField,
            Do[
               (*\[Lambda] spinors, (mod prime)*)
               SpinorUndotN[labels[[i]]][$lam][$up]=Mod[la[i],prime];
               SpinorUndotN[labels[[i]]][$lam][$down]=Mod[ra[i],prime];
               SpinorDotN[labels[[i]]][$lam][$down]=Mod[ls[i],prime];
               SpinorDotN[labels[[i]]][$lam][$up]=Mod[rs[i],prime];
            ,{i,num}];
            If[OptionValue[FixedReference]=!=False,
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up]=Mod[la[num+1],prime];
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down]=Mod[ra[num+1],prime];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down]=Mod[ls[num+1],prime];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up]=Mod[rs[num+1],prime];
            ];
         ,
            Do[
               (*\[Lambda] spinors*)
               SpinorUndotN[labels[[i]]][$lam][$up]=la[i];
               SpinorUndotN[labels[[i]]][$lam][$down]=ra[i];
               SpinorDotN[labels[[i]]][$lam][$down]=ls[i];
               SpinorDotN[labels[[i]]][$lam][$up]=rs[i];
            ,{i,num}];
            If[OptionValue[FixedReference]=!=False,
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up]=la[num+1];
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down]=ra[num+1];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down]=ls[num+1];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up]=rs[num+1];

               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;
               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;

               ExtramassN[OptionValue[FixedReference][[1]]]=0;
               ExtramasstildeN[OptionValue[FixedReference][[1]]]=0;
            ];
         ];

         Activate @ Table[{
            (*\[Lambda] Spinor components*)
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$lam][j:1|2][Null],Indexed[SpinorUndotN[labels[[i]]][$lam][$up],{j}]],
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$lam][Null][j:1|2],Indexed[SpinorUndotN[labels[[i]]][$lam][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$lam][Null][j:1|2],Indexed[SpinorDotN[labels[[i]]][$lam][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$lam][j:1|2][Null],Indexed[SpinorDotN[labels[[i]]][$lam][$up],{j}]]}
         ,{i,num}];
         If[OptionValue[FixedReference]=!=False,
            Activate @ {
               Inactive[SetDelayed][SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][j:1|2][Null],Indexed[SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up],{j}]],
               Inactive[SetDelayed][SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][Null][j:1|2],Indexed[SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down],{j}]],
               Inactive[SetDelayed][SpinorDotN[OptionValue[FixedReference][[1]]][$lam][Null][j:1|2],Indexed[SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down],{j}]],
               Inactive[SetDelayed][SpinorDotN[OptionValue[FixedReference][[1]]][$lam][j:1|2][Null],Indexed[SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up],{j}]]
            };
            SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][_]:=0;
            SpinorDotN[OptionValue[FixedReference][[1]]][$mu][_]:=0;
         ];

	(*In case dim=4 because of n-1 or n 4D momenta being specified in FourD:*)
	Do[
		(*Setting \[Mu] spinors to zero*)
		SpinorUndotN[labels[[i]]][$mu][$up]={0,0};
		SpinorUndotN[labels[[i]]][$mu][$down]={0,0};
		SpinorDotN[labels[[i]]][$mu][$down]={0,0};
		SpinorDotN[labels[[i]]][$mu][$up]={0,0};
		(*Setting masses to zero*)
		ExtramassN[labels[[i]]]=0;
		ExtramasstildeN[labels[[i]]]=0;
	,{i,num}];
	Do[
		(*Setting \[Mu] Spinor components set to zero*)
		SpinorUndotN[labels[[i]]][$mu][j][Null]=0;
		SpinorUndotN[labels[[i]]][$mu][Null][j]=0;
		SpinorDotN[labels[[i]]][$mu][j][Null]=0;
		SpinorDotN[labels[[i]]][$mu][Null][j]=0;
	,{i,num},{j,{1,2}}];

		Do[Mom4DVectorNFromMatrix[lbl][Null],{lbl,labels}];
		Do[Mom4DVectorNFromMatrix[lbl][$flat],{lbl,labels}];

         (*Print the computed spinors if required*)

	If[TrueQ[OptionValue[DisplaySpinors]],
		Print["Output reads: {|\[Lambda]>,|\[Lambda]],|\[Mu]>,|\[Mu]]}"];
		Print[Table[{ra[i],rs[i],{0,0},{0,0}},{i,num}]];
	,
		Null
	]
];


(* ::Subsubsection:: *)
(*GenSpinors 6D*)


Options[genSpinors6D]={DisplaySpinors->False,Parametric->False,ParameterName->$par,ParameterRange->200,Modulus->NextPrime[10000,-1],
	SpecialKinematics->Rational,Seed->False,RandomSpinors->False,FixedReference->False,FourD->{}};

genSpinors6D[labels_List,opts:OptionsPattern[]]:=Module[
{lambda,lambdaout,num,ra,rs,la,ls,\[Xi],\[Eta],\[Xi]t,\[Eta]t,mass,masst,mur,mul,mutr,mutl,parn,usedpar,system,sol,dim,prime,
	freeIndices,fixedIndices,fourDIndices,solvVars},

	If[OptionValue[SpecialKinematics]===FiniteField,
	   If[PrimeQ[OptionValue[Modulus]],
	      prime = Abs[OptionValue[Modulus]];
		,
	      Message[GenSpinors::NotPrime];
	      Abort[];
	   ];
	   If[TrueQ[OptionValue[Parametric]],
	      Message[GenSpinors::FiniteFieldParametric];
	      Abort[];
	   ];
	];

	(*Definition of the spinors in terms of the components. In the six-dimensionla case there will be 2n of tese spinors, where the first n refer to \[Lambda] and the secodn n to \[Lambda]' which is redefinition of the \[Mu] encoding also the masses, see paper to be written.*)
	ra[i_]:={\[Xi][i],\[Eta][i]}(*|i>*);
	rs[i_]:={\[Xi]t[i],\[Eta]t[i]}(*|i]*);
	la[i_]:={\[Eta][i],-\[Xi][i]}(*<i|*);
	ls[i_]:={-\[Eta]t[i],\[Xi]t[i]}(*[i|*);

	num=Length[labels];
	
	fourDIndices = DeleteDuplicates[Position[labels,#][[1,1]]&/@OptionValue[FourD]];
	
         (*If we want six-dimensional objects we have to generate also the \[Mu] and the masses. We will do this starting from scratch, which should be much easier than recycle momentum twistor variables*)

         (*Determining the indices of the momenta for which components have to be randomly chosen, and for which the system below can be solved.
           This depends on the 4D momenta: the vanishing \[Xi](t),\[Eta](t)[num+i] can multiply the variables for which we want to solve. Therefore
           we choose our variables such that at least two of the \[Xi]t we solve for belong to 6D momenta (of which there are always at least two).*)
         (*freeIndices=Join[OptionValue[FourD],Complement[Range[1,num],OptionValue[FourD]][[3;;]]];*)
         freeIndices=Join[fourDIndices,Complement[Range[1,num],fourDIndices][[3;;]]];
         fixedIndices=Complement[Range[1,num],freeIndices];

         (*Fix 3n+9 arbitrary components, differentiating between rational, finite field, and floating point kinematics*)
         Switch[OptionValue[SpecialKinematics],
            FiniteField,
               (*3n*)
               Do[\[Xi][i]=RandomInteger[{1,prime}];
                  \[Eta][i+num]=RandomInteger[{1,prime}];
                  \[Xi]t[i+num]=RandomInteger[{1,prime}],
               {i,num}];
               (*The other 9 arbitrary ones*)
               \[Eta]t[1+num]=RandomInteger[{1,prime}];
               \[Eta]t[2+num]=RandomInteger[{1,prime}];
               \[Eta]t[3+num]=RandomInteger[{1,prime}];
               \[Eta]t[4+num]=RandomInteger[{1,prime}];
               \[Xi][1+num]=RandomInteger[{1,prime}];
               \[Xi][2+num]=RandomInteger[{1,prime}];
               \[Xi][3+num]=RandomInteger[{1,prime}];
               \[Xi][4+num]=RandomInteger[{1,prime}];
               \[Xi]t[freeIndices[[2]]]=RandomInteger[{1,prime}];,
            Rational | _,
               (*3n*)
               Do[\[Xi][i]=RandomInteger[{1,OptionValue[ParameterRange]}];
                  \[Eta][i+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
                  \[Xi]t[i+num]=RandomInteger[{1,OptionValue[ParameterRange]}],
               {i,num}];
               (*The other 9 arbitrary ones*)
               \[Eta]t[1+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Eta]t[2+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Eta]t[3+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Eta]t[4+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Xi][1+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Xi][2+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Xi][3+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Xi][4+num]=RandomInteger[{1,OptionValue[ParameterRange]}];
               \[Xi]t[freeIndices[[2]]]=RandomInteger[{1,OptionValue[ParameterRange]}];
         ];

         (*Now we have to set the free-parameters to either be equal to the chosen variable name or generate numeric values for them, depending on the OptionValue of Parametric*)
         If[TrueQ[OptionValue[Parametric]],
            (*Parametric expression*)
            (*In order to write less lets introduce this local variable*)
            parn=OptionValue[ParameterName];
            usedpar=1;
            Do[\[Xi]t[i]=parn[usedpar];
               usedpar++;,{i,freeIndices[[3;;]]}];
            (* debugPrint[Table[\[Xi]t[i]=parn[i-Length@freeIndices+2],{i,freeIndices[[3;;]]}; *)
            Do[\[Eta]t[i]=parn[usedpar];
               usedpar++;,{i,Complement[Range[1,num],{freeIndices[[1]]}]}];
            (* debugPrint[Table[\[Eta]t[i]=parn[i-1+usedpar],{i,2,num}]]; *)
            Do[\[Xi][i+num]=parn[usedpar];
               usedpar++;,{i,5,num}];
            (* debugPrint[Table[\[Xi][i+num]=parn[i-4+usedpar],{i,5,num}]]; *)
            Do[\[Eta][i]=parn[usedpar];
               usedpar++;,{i,freeIndices}];
            (* debugPrint[Table[\[Eta][i+num]=parn[i-2+usedpar],{i,3,num}]]; *)
            Do[\[Eta]t[i+num]=parn[usedpar];
               usedpar++,{i,5,num}];
            (* debugPrint[Table[\[Eta]t[i+num]=parn[i-4+usedpar],{i,5,num}]]; *)
            (* debugPrint[usedpar==(5*num-15)] *)
            ,
            (*Numerical values*)
            Switch[OptionValue[SpecialKinematics],
               FiniteField,
                  (*In the case finite field kinematics is required*)
                  Do[\[Xi]t[i]=RandomInteger[{1,prime}],{i,freeIndices[[3;;]]}];
                  Do[\[Eta]t[i]=RandomInteger[{1,prime}],{i,Complement[Range[1,num],{freeIndices[[1]]}]}];
                  Do[\[Xi][i+num]=RandomInteger[{1,prime}],{i,5,num}];
                  Do[\[Eta][i]=RandomInteger[{1,prime}],{i,freeIndices}];
                  Do[\[Eta]t[i+num]=RandomInteger[{1,prime}],{i,5,num}],
               Rational | _,
                  (*In the case rational kinematics is required*)
                  Do[\[Xi]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,freeIndices[[3;;]]}];
                  Do[\[Eta]t[i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,Complement[Range[1,num],{freeIndices[[1]]}]}];
                  Do[\[Xi][i+num]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,5,num}];
                  Do[\[Eta][i]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,freeIndices}];
                  Do[\[Eta]t[i+num]=RandomInteger[{1,OptionValue[ParameterRange]}],{i,5,num}]
            ];
         ];

      If[OptionValue[FixedReference] =!= False,
         (*The +1 definitions will be used later to set the global variables (as they are guaranteed not to be overwritten by having 4D momenta)*)
         Do[\[Xi][num+i]=OptionValue[FixedReference][[2,1]],{i,1,num+1}];
         Do[\[Eta][num+i]=OptionValue[FixedReference][[2,2]],{i,1,num+1}];
         Do[\[Xi]t[num+i]=OptionValue[FixedReference][[2,3]],{i,1,num+1}];
         Do[\[Eta]t[num+i]=OptionValue[FixedReference][[2,4]],{i,1,num+1}];
      ];

	Do[
		Module[{labelIdx=Position[labels,label][[1,1]]},
			\[Xi][num+labelIdx]=0;
			\[Xi]t[num+labelIdx]=0;
			\[Eta][num+labelIdx]=0;
			\[Eta]t[num+labelIdx]=0;
		];
	,{label,OptionValue[FourD]}];

         (*Now the last 6 variables fixed by momentum conservation. If RandomSpinors is set to True momentum conservation is not imposed and the last variables are fixed randomly.*)

         solvVars={\[Xi]t[fixedIndices[[1]]],\[Xi]t[fixedIndices[[2]]],\[Xi]t[freeIndices[[1]]],\[Eta][fixedIndices[[1]]],\[Eta][fixedIndices[[2]]],\[Eta]t[freeIndices[[1]]]};

         If[TrueQ[OptionValue[RandomSpinors]],
            Evaluate[solvVars]=Table[RandomInteger[{1,OptionValue[ParameterRange]}],{i,6}],
            system={Sum[\[Xi][i]\[Xi]t[i]+\[Xi][i+num]\[Xi]t[i+num],{i,num}]==0,
                    Sum[\[Xi][i]\[Eta]t[i]+\[Xi][i+num]\[Eta]t[i+num],{i,num}]==0,
                    Sum[\[Eta][i]\[Xi]t[i]+\[Eta][i+num]\[Xi]t[i+num],{i,num}]==0,
                    Sum[\[Eta][i]\[Eta]t[i]+\[Eta][i+num]\[Eta]t[i+num],{i,num}]==0,
                    Sum[\[Xi][i]\[Eta][i+num]-\[Xi][i+num]\[Eta][i],{i,num}]==0,
                    Sum[\[Xi]t[i]\[Eta]t[i+num]-\[Xi]t[i+num]\[Eta]t[i],{i,num}]==0};
            Switch[OptionValue[SpecialKinematics],
               FiniteField,
                  sol=Solve[system,solvVars,Modulus->prime]//Flatten;
                  Evaluate[solvVars]=Mod[solvVars/.sol,prime];,
               Rational,
                  sol=Solve[system,solvVars]//Flatten;
                  Evaluate[solvVars]=solvVars/.sol;,
               _,
                  sol=Solve[system,solvVars]//Flatten;
                  Evaluate[solvVars]=N[solvVars/.sol];
            ];
         ];


         (*Now we have to write the masses and the \[Mu] in terms of these variables. For the relations among these see the paper to come...*)


         (*We choose the \[Mu] to be equal to \[Lambda]' which is an allowed choice, so nothing has to be done there.*)

         (*Now we just need to set the global variables equal to te locally generated ones:*)

         If[OptionValue[SpecialKinematics]===FiniteField,
            Do[
               (*\[Lambda] spinors*)
               SpinorUndotN[labels[[i]]][$lam][$up]=Mod[la[i],prime];
               SpinorUndotN[labels[[i]]][$lam][$down]=Mod[ra[i],prime];
               SpinorDotN[labels[[i]]][$lam][$down]=Mod[ls[i],prime];
               SpinorDotN[labels[[i]]][$lam][$up]=Mod[rs[i],prime];
               (*\[Mu] spinors*)
               SpinorUndotN[labels[[i]]][$mu][$up]=Mod[la[i+num],prime];
               SpinorUndotN[labels[[i]]][$mu][$down]=Mod[ra[i+num],prime];
               SpinorDotN[labels[[i]]][$mu][$down]=Mod[ls[i+num],prime];
               SpinorDotN[labels[[i]]][$mu][$up]=Mod[rs[i+num],prime];
               (*Masses:*)
               ExtramassN[labels[[i]]]=Mod[la[i] . ra[i+num],prime];
               ExtramasstildeN[labels[[i]]]=Mod[ls[i+num] . rs[i],prime];
            ,{i,num}];
            (*Fixed reference spinors*)
            If[OptionValue[FixedReference]=!=False,
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up]=Mod[la[2*num+1],prime];
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down]=Mod[ra[2*num+1],prime];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down]=Mod[ls[2*num+1],prime];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up]=Mod[rs[2*num+1],prime];

               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;
               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;

               ExtramassN[OptionValue[FixedReference][[1]]]=0;
               ExtramasstildeN[OptionValue[FixedReference][[1]]]=0;
            ];,
            Do[
               (*\[Lambda] spinors*)
               SpinorUndotN[labels[[i]]][$lam][$up]=la[i];
               SpinorUndotN[labels[[i]]][$lam][$down]=ra[i];
               SpinorDotN[labels[[i]]][$lam][$down]=ls[i];
               SpinorDotN[labels[[i]]][$lam][$up]=rs[i];
               (*\[Mu] spinors*)
               SpinorUndotN[labels[[i]]][$mu][$up]=la[i+num];
               SpinorUndotN[labels[[i]]][$mu][$down]=ra[i+num];
               SpinorDotN[labels[[i]]][$mu][$down]=ls[i+num];
               SpinorDotN[labels[[i]]][$mu][$up]=rs[i+num];
               (*Masses:*)
               ExtramassN[labels[[i]]]=la[i] . ra[i+num];
               ExtramasstildeN[labels[[i]]]=ls[i+num] . rs[i];
            ,{i,num}];
            (*Fixed reference spinors*)
            If[OptionValue[FixedReference]=!=False,
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up]=la[2*num+1];
               SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down]=ra[2*num+1];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down]=ls[2*num+1];
               SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up]=rs[2*num+1];

               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;
               SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$down]=0;
               SpinorDotN[OptionValue[FixedReference][[1]]][$mu][$up]=0;

               ExtramassN[OptionValue[FixedReference][[1]]]=0;
               ExtramasstildeN[OptionValue[FixedReference][[1]]]=0;
            ];
         ];

         Activate @ Table[{
            (*\[Lambda] Spinor components*)
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$lam][j:1|2][Null],Indexed[SpinorUndotN[labels[[i]]][$lam][$up],{j}]],
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$lam][Null][j:1|2],Indexed[SpinorUndotN[labels[[i]]][$lam][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$lam][Null][j:1|2],Indexed[SpinorDotN[labels[[i]]][$lam][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$lam][j:1|2][Null],Indexed[SpinorDotN[labels[[i]]][$lam][$up],{j}]],
            (*\[Mu] Spinor components*)
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$mu][j:1|2][Null],Indexed[SpinorUndotN[labels[[i]]][$mu][$up],{j}]],
            Inactive[SetDelayed][SpinorUndotN[labels[[i]]][$mu][Null][j:1|2],Indexed[SpinorUndotN[labels[[i]]][$mu][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$mu][Null][j:1|2],Indexed[SpinorDotN[labels[[i]]][$mu][$down],{j}]],
            Inactive[SetDelayed][SpinorDotN[labels[[i]]][$mu][j:1|2][Null],Indexed[SpinorDotN[labels[[i]]][$mu][$up],{j}]]}
         ,{i,num}];
         (*Fixed reference spinor components*)
         If[OptionValue[FixedReference]=!=False,
            Activate @ {
               Inactive[SetDelayed][SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][j:1|2][Null],Indexed[SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$up],{j}]],
               Inactive[SetDelayed][SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][Null][j:1|2],Indexed[SpinorUndotN[OptionValue[FixedReference][[1]]][$lam][$down],{j}]],
               Inactive[SetDelayed][SpinorDotN[OptionValue[FixedReference][[1]]][$lam][Null][j:1|2],Indexed[SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$down],{j}]],
               Inactive[SetDelayed][SpinorDotN[OptionValue[FixedReference][[1]]][$lam][j:1|2][Null],Indexed[SpinorDotN[OptionValue[FixedReference][[1]]][$lam][$up],{j}]]
            };
            SpinorUndotN[OptionValue[FixedReference][[1]]][$mu][_]:=0;
            SpinorDotN[OptionValue[FixedReference][[1]]][$mu][_]:=0;
         ];

         Do[Mom4DVectorNFromMatrix[lbl][Null],{lbl,labels}];
         Do[Mom4DVectorNFromMatrix[lbl][$flat],{lbl,labels}];


         (*Print the computed spinors if required*)

         If[TrueQ[OptionValue[DisplaySpinors]],
            Print["Output reads: {|\[Lambda]>,|\[Lambda]],|\[Mu]>,|\[Mu]]}"];
            Return[Table[{ra[i],rs[i],ra[i+num],rs[i+num]},{i,num}]],
            Return[Null];
         ];
];


(* ::Subsection::Closed:: *)
(*ClearKinematics*)


ClearSubValues[label_,function_] := SubValues[function] = DeleteCases[SubValues[function],subV_?(!FreeQ[First[#],HoldPattern[function[label]]]&)];
ClearSubValues[function_] := SubValues[function] = DeleteCases[SubValues[function],subV_?(!FreeQ[First[#],HoldPattern[function[_]]]&)];


ClearKinematics[labels_List] := ClearKinematics /@ labels;
ClearKinematics[labels__] := ClearKinematics /@ {labels};

ClearKinematics[label_]:=(
	ClearSubValues[label,SpinorUndotN];
	ClearSubValues[label,SpinorDotN];
	ClearSubValues[label,Mom4DN];
	ClearSubValues[label,Mom4DVectorN];
	ClearSubValues[label,ExtramassN];
	ClearSubValues[label,ExtramasstildeN];
	
	Unprotect[Extramass,Extramasstilde];
	ClearSubValues[label,Extramass];
	ClearSubValues[label,Extramasstilde];
	Protect[Extramass,Extramasstilde];
);


ClearKinematics[]:=(
	ClearSubValues[SpinorUndotN];
	ClearSubValues[SpinorDotN];
  ClearSubValues[label,Mom4DN];
  ClearSubValues[label,Mom4DVectorN];
	ClearSubValues[ExtramassN];
	ClearSubValues[ExtramasstildeN];
	
	Unprotect[Extramass,Extramasstilde];
	ClearSubValues[Extramass];
	ClearSubValues[Extramasstilde];
	Protect[Extramass,Extramasstilde];
);


(* ::Subsection:: *)
(*ToNum*)

$ToNumRules = {S -> S6N, S6Mom -> S6MomN, S4Mom -> S4MomN, LorentzProduct4D -> LorentzProduct4DN,
  SpinorAngleBracket->SpinorAngleBracketN,SpinorSquareBracket->SpinorSquareBracketN,SpinorDot->SpinorDotN,SpinorUndot->SpinorUndotN,
  Chain->ChainN,SpinorHelicity6D`SpinorFunctions`TrM->SpinorHelicity6D`SpinorFunctions`TrMN,SpinorHelicity6D`SpinorFunctions`TrP->SpinorHelicity6D`SpinorFunctions`TrPN,
  SpinorHelicity6D`SpinorFunctions`Tr5->SpinorHelicity6D`SpinorFunctions`Tr5N,
  Mom4DVector->Mom4DVectorN, Extramass->ExtramassN,Extramasstilde->ExtramasstildeN,Mom4D -> Mom4DN
};

$ToNumHeads = $ToNumRules[[All,1]]

Table[
  (ToNumEval[x:headPattern@#[[1]]]:= x/.#) &@functionRule,
  {functionRule,$ToNumRules}
];

ToNumEval[Verbatim[Plus][expr___]] := Plus@@ToNumEval/@{expr};
ToNumEval[Verbatim[Times][expr___]] := Times@@ToNumEval/@{expr};
ToNumEval[Verbatim[List][expr___]] := List@@ToNumEval/@{expr};
ToNumEval[Verbatim[Power][expr___]] := Power@@ToNumEval/@{expr};
ToNumEval[Verbatim[Dot][expr___]] := Dot@@ToNumEval/@{expr};

ToNumEval[x_]:=x;

(*applyToNumRules[expr_] :=expr//.{S -> S6, S6Mom -> S6MomN, S4Mom -> S4MomN, LorentzProduct4D -> LorentzProduct4DN}//.{*)
(*      SpinorAngleBracket->SpinorAngleBracketN,SpinorSquareBracket->SpinorSquareBracketN,SpinorDot->SpinorDotN,SpinorUndot->SpinorUndotN,*)
(*      Chain->ChainN,SpinorHelicity6D`SpinorFunctions`TrM->SpinorHelicity6D`SpinorFunctions`TrMN,SpinorHelicity6D`SpinorFunctions`TrP->SpinorHelicity6D`SpinorFunctions`TrPN,*)
(*      SpinorHelicity6D`SpinorFunctions`Tr5->SpinorHelicity6D`SpinorFunctions`Tr5N, Mom4DVector->Mom4DVectorN, Extramass->ExtramassN,Extramasstilde->ExtramasstildeN*)
(*};*)
                          
(*ToNum[expr_] := Module[{uniqueObjects,uniqueObjectsEvaluated,evaluationRules},*)

(*	uniqueObjects = DeleteDuplicates@Cases[{expr},*)
(*		headPattern[S]|headPattern[S6Mom]|headPattern[S4Mom]|headPattern[SpinorAngleBracket]|headPattern[SpinorSquareBracket]|headPattern[SpinorDot]|headPattern[SpinorUndot]|*)
(*		headPattern[Chain]|headPattern[SpinorHelicity6D`SpinorFunctions`TrM]|headPattern[SpinorHelicity6D`SpinorFunctions`TrP]|headPattern[SpinorHelicity6D`SpinorFunctions`Tr5]|*)
(*		headPattern[Mom4DVector]|headPattern[Extramass]|headPattern[Extramasstilde], *)
(*	Infinity];*)
(*	uniqueObjectsEvaluated = applyToNumRules@uniqueObjects;*)
(*	evaluationRules = Dispatch@MapThread[Rule,{uniqueObjects,uniqueObjectsEvaluated}];*)
(*	*)
(*	expr /. evaluationRules*)
(*]*)

setupToNumCache[]:=
    Module[{},
      ClearAll[ToNumCache];

      ToNumCache[Verbatim[Plus][expr___]] :=
          ToNumCache[Plus[expr]] = Plus@@ToNumCache/@{expr};
      ToNumCache[Verbatim[Times][expr___]] :=
          ToNumCache[Times[expr]] = Times@@ToNumCache/@{expr};
      ToNumCache[Verbatim[List][expr___]] :=
          ToNumCache[List[expr]] = List@@ToNumCache/@{expr};
      ToNumCache[Verbatim[Power][expr___]] :=
          ToNumCache[Power[expr]] = Power@@ToNumCache/@{expr};
      ToNumCache[Verbatim[Dot][expr___]] :=
          ToNumCache[Dot[expr]] = Dot@@ToNumCache/@{expr};
      ToNumCache[Verbatim[Log][args___,expr_]] :=
          ToNumCache[Log[args,expr]] = Log[args,ToNumCache@expr];
      ToNumCache[Verbatim[PolyLog][args__,expr_]] :=
          ToNumCache[PolyLog[args,expr]] = PolyLog[args,ToNumCache@expr];


      ToNumCache[expr:Alternatives@@headPattern/@$ToNumHeads]:=
          ToNumCache[expr] = ToNumEval[expr];
      ToNumCache[x_] := x;
    ];

ToNum[expr_]:=Module[{},
  setupToNumCache[];
  ToNumCache@expr
];


(* ::Subsection::Closed:: *)
(*Cyclic*)


(*Auxiliary function cyclic for Cyclic*)

cyclic[exp_,cycle_List]:=Module[{localexp,reps,localcycle},
(*Consider the possibility of many cycles*)
If[TrueQ[Head[cycle[[1]]]==List],
localcycle=cycle,
localcycle={cycle};
];
(*Transform cycle_List into a list of replacements*)
reps={};
Do[
AppendTo[reps,Join[Table[{subcycle[[i]]->subcycle[[i+1]],ToExpression["r"<>ToString[subcycle[[i]]]]->ToExpression["r"<>ToString[subcycle[[i+1]]]]},{i,Length[subcycle]-1}],{subcycle[[-1]]->subcycle[[1]],ToExpression["r"<>ToString[subcycle[[-1]]]]->ToExpression["r"<>ToString[subcycle[[1]]]]}]//Flatten];
,{subcycle,localcycle}];
reps=reps//Flatten;
(*Apply replacements to the functions in the package.*)
localexp=exp/.{MomPure[x_]:>(MomPure[x]/.reps),PolarPure[x__]:>(PolarPure[x]/.reps),S[x__]:>(S[x]/.reps),S4[x__]:>(S4[x]/.reps),Mom[x_]:>(Mom[x]/.reps),Polar[x__]:>(Polar[x]/.reps),SpinorAngleBracket[x__]:>(SpinorAngleBracket[x]/.reps),SpinorSquareBracket[x__]:>(SpinorSquareBracket[x]/.reps),AngSquInvariant[x_]:>(AngSquInvariant[x]/.reps),SquAngInvariant[x__]:>(SquAngInvariant[x]/.reps),AngAngInvariant[x__]:>(AngAngInvariant[x]/.reps),SquSquInvariant[x__]:>(SquSquInvariant[x]/.reps),SpinorUndot[x__]:>(SpinorUndot[x]/.reps),SpinorDot[x__]:>(SpinorDot[x]/.reps),SpinorUndot6D[x__]:>(SpinorUndot6D[x]/.reps),SpinorDot6D[x__]:>(SpinorDot6D[x]/.reps),SpinorDotPure[x__]:>(SpinorDotPure[x]/.reps),SpinorUndotPure[x__]:>(SpinorUndotPure[x]/.reps)};
Return[localexp];
];

(*Options*)
Options[Cyclic]={SumAll->False};

(*The actual function*)
Cyclic[exp_,cycle_List,OptionsPattern[]]:=Module[{len,addenda,out},
If[TrueQ[OptionValue[SumAll]]&&!TrueQ[Head[cycle[[1]]]==List],
len=Length[cycle];
addenda[1]=exp;
Do[addenda[i]=cyclic[addenda[i-1],cycle],{i,2,len}];
out=Sum[addenda[i],{i,len}];
,
out=cyclic[exp,cycle];
];
Return[out];
];


(* ::Subsection::Closed:: *)
(*FeynCToSpinor*)


FeynCToSpinor[exp_,momenta_List]:=Block[{Eta,locexp,Sreps,momrep,proprep1,proprep2,Pair},
Sreps=Table[ToExpression["s"<>ToString[i]<>ToString[j]]->S[i,j],{i,momenta},{j,momenta}]//Flatten;
locexp=exp/.Sreps;
momrep={Pair[LorentzIndex[m_,___],Momentum[p1_,___]]:>Mom[p1][m][Null],Pair[LorentzIndex[\[Mu]_,___],LorentzIndex[\[Nu]_,___]]:>Eta[\[Mu],\[Nu]][$up]};
locexp=locexp/.momrep;
proprep1={FeynAmpDenominator[x__]:>1/(Times[x])};
proprep2={PropagatorDenominator[x_,m_]:>mp[x,x]-mp[m,m]};
locexp=locexp/.proprep1;
locexp=locexp/.proprep2;
proprep1={mp[0,0]->0,Momentum[x_,_]:>MomPure[x]};
locexp=locexp/.proprep1;
Return[locexp];
];


(* ::Subsection::Closed:: *)
(*MomReplace*)


MomReplace[exp_,reps_List,OptionsPattern[]]:=Block[{locexp,localden,count1,count2,\[Mu],\[Nu],locreps,subcounter,Eta2,loc,loopmom,len,backrep},

(*Set up: start by creating a list of th emomenta to be replaced*)
loopmom=Table[i[[1]],{i,reps}];
loopmom=loopmom/.{Mom[x_][_][_]:>x};
(*Definition of counters needed later on in the replacemnets*)
count1=1;
count2=2;
subcounter:=IntegerPart[count2++/2];
(*Definition of the replecements*)
locreps={S[i_,j_]:>2*Mom[i][\[Mu][subcounter]][Null]Eta2[\[Mu][subcounter],\[Mu][subcounter]][$down]Mom[j][\[Mu][subcounter]][Null],chain[type1_,a_,{b___,x_,c___},d_,type2_]/;MemberQ[loopmom,x]:>chain[type1,a,{b,\[Nu][loc=count1++],c},d,type2]Mom[x][\[Nu][loc]][Null],Chain[type1_,a_,{b___,x_,c___},d_,type2_]/;MemberQ[loopmom,x]:>Chain[type1,a,locexp=Plus@@Table[locexp[[i]]/localden[[i]],{i,len}];{b,\[Nu][loc=count1++],c},d,type2]Mom[x][\[Nu][loc]][Null]};
backrep={chain[type1_,a_,{b___,\[Mu]_,c___},d_,type2_]Mom[x_][\[Mu]_][Null]:>chain[type1,a,{b,x,c},d,type2],Chain[type1_,a_,{b___,\[Mu]_,c___},d_,type2_]Mom[x_][\[Mu]_][Null]:>Chain[type1,a,{b,x,c},d,type2]};

(*Create a list of numerators and denominators where the replacemnt has to be performed. this separation is needed in order to handle the replacements in the denominator safely*)
locexp=exp//Expand;
If[TrueQ[Head[locexp]==Plus],
locexp=List@@locexp,
locexp={locexp}//Flatten[#,1]&;
];
debugPrint[locexp];
len=Length[locexp];
localden=Table[Denominator[i],{i,locexp}];
debugPrint[localden];
locexp=Table[Numerator[i],{i,locexp}];
debugPrint[locexp];

(*Apply replecemnts in the numerator*)
locexp=locexp/.locreps;
locexp=locexp/.reps;
locexp=locexp/.{Eta2->Eta}//Expand;
debugPrint[locexp];
locexp=locexp//.backrep;
debugPrint[locexp];
(*Differentiate among the possibilities of performing the replacement in the numerator only or in the denominator as well*)
If[TrueQ[OptionValue[NumeratorOnly]],
locexp=Plus@@Table[locexp[[i]]/localden[[i]],{i,len}],
localden=localden/.locreps;
localden=localden/.reps;
localden=localden/.{Eta2->Eta}//Expand;
debugPrint[localden];
localden=localden//.backrep;
debugPrint[localden];
locexp=Plus@@Table[locexp[[i]]/localden[[i]],{i,len}];
];
Return[locexp];
];


(* ::Section::Closed:: *)
(*Symbolic calculus*)


(* ::Subsection::Closed:: *)
(*Fstrength*)


Fstrength[label_][A_,B_][C_,D_][lg1_Symbol,lg2_Symbol]:=2*Antisymmetrize[MDelta[4][A ][C]SpinorUndot6D[label ][B ][Null][lg1]SpinorDot6D[label ][D ][Null][lg2],{A,B},{C,D}];

(*Boxing function*)
FstrengthBox[mom_,AA_,BB_,CC_,DD_,lg1_,lg2_]:=TemplateBox[{mom,AA,BB,CC,DD,lg1,lg2},"Fstrength",
DisplayFunction->(RowBox[{SubscriptBox[SuperscriptBox["F",RowBox[{#2,#3}]],RowBox[{#4,#5}]],"[",SubscriptBox[#1,RowBox[{#6,#7}]],"]"}]&),
InterpretationFunction->(RowBox[{"Fstrength[",#1,"][",#2,",",#3,"][",#4,",",#5,"][",#6,",",#7,"]"}]&)
];

Fstrength /: MakeBoxes[Fstrength[mom_][AA_,BB_][CC_,DD_][lg1_Symbol,lg2_Symbol],StandardForm|TraditionalForm]:=FstrengthBox[ToBoxes[mom],ToBoxes[AA],ToBoxes[BB],ToBoxes[CC],ToBoxes[DD],ToBoxes[lg1],ToBoxes[lg2]];

If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "fst" -> FstrengthBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*MDelta*)


(*Properties of MDelta*)
(*Trace*)
MDelta[dim_Integer][up1_][down1_] /;TrueQ[up1==down1]:=dim;

(*Kronecker delta for integers*)
MDelta[dim_Integer][up1_Integer][down1_Integer]:=If[TrueQ[up1==down1],1,0];

(*Contraction of two deltas*)
MDelta /: Times[MDelta[dim_][up1_][down1_],MDelta[dim_][down1_][down2_]]:=MDelta[dim][up1][down2];

(*Define contractions with SpinorDot*)

MDelta /: Times[MDelta[4][x_][y_],SpinorDot6D[label_][y_][lg1_][lg2_]]:=SpinorDot6D[label][x][lg1][lg2];
MDelta /: Times[MDelta[4][x_][y_],SpinorDot6D[label_][x_][lg1_][lg2_]]:=SpinorDot6D[label][y][lg1][lg2];
MDelta /: Times[MDelta[4][x_][y_],SpinorUndot6D[label_][y_][lg1_][lg2_]]:=SpinorUndot6D[label][x][lg1][lg2];
MDelta /: Times[MDelta[4][x_][y_],SpinorUndot6D[label_][x_][lg1_][lg2_]]:=SpinorUndot6D[label][y][lg1][lg2];

(*Proper visualization through custom box function*)
MDeltaBox[dim_,up_,down_]:=TemplateBox[{up,down,dim},"MDelta",
DisplayFunction->(RowBox[{SubscriptBox[SuperscriptBox["\[Delta]",#1],#2]}]&),
InterpretationFunction->(RowBox[{"MDelta","[",#3,"][",#1,"][",#2,"]"}]&)
];

(*Define the action of MakeBoxes on MDelta*)
MDelta /: MakeBoxes[MDelta[dim_][up_][down_],StandardForm|TraditionalForm]:=MDeltaBox[ToBoxes[dim],ToBoxes[up],ToBoxes[down]];

(*Definion of the default dimension of the delta*)
$MDimension=4;

(*Shortcuts*)
If[frontend==1,
SetOptions[EvaluationNotebook[],
    InputAliases -> DeleteDuplicates@Append[InputAliases /. Options[EvaluationNotebook[], InputAliases], "md" -> MDeltaBox[$MDimension,"\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]];
    ];


(* ::Subsection::Closed:: *)
(*Antisymmetrize*)


SetAttributes[Antisymmetrize,HoldFirst];

Antisymmetrize[exp_,indices__List]:=Block[{canonicalsign,len,combos,localexp,output,localindices,nsubs},
(*Start by counting the number of different lists given as input*)
localindices={indices};
nsubs=Length[localindices];

output=HoldForm[exp];

(*Now we loop over every single list of antisymmetrized indices*)
Do[
(*Since mathematica has its own canonical ordering we need to take into account that what we started from may already have a minus sign. So we simply use the starting lists sign as normalization*)
canonicalsign=Signature[localindices[[k]]];

(*Now generate all the possible reshuffelings of the indices in list*)
len=Length[localindices[[k]]];
combos=Select[DeleteDuplicates[DeleteDuplicates/@Tuples[Table[localindices[[k]],{j,len}]]],(Length[#]==len)&];
debugPrint[combos];

(*Now do the replacements in the expression*)
localexp=output;
output={};
Do[AppendTo[output,(Signature[combos[[i]]]/canonicalsign)*localexp/.Table[localindices[[k,j]]->combos[[i,j]],{j,len}]],{i,Length[combos]}];
output=(1/(len!))*Plus@@output//Expand;
,{k,nsubs}];

(*Before returning the result release the hold on the expression*)
output=ReleaseHold[output];
Return[output];
];


(* ::Subsection::Closed:: *)
(*ContractReplace*)


ContractReplace[exp_]:=Block[{locexp,angsqurep,squangrep,angsqurepdot,squangrepdot,backangsqu,backsquang},
(*Just an enourmous and probably unefficient list of replacements coming with appropriate minus signs...*)
angsqurep={Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;OrderedQ[{x,z}]&&x!=z:>AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot],Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;OrderedQ[{x,z}]==False:>-AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot],Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&OrderedQ[{y,k}]&&y!=k:>AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&OrderedQ[{y,k}]==False:>-AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&TrueQ[y==k]&&OrderedQ[{adot,bdot}]:>AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&TrueQ[y==k]&&OrderedQ[{adot,bdot}]==False:>-AngSquInvariant[x,y][{a,b},adot]*AngSquInvariant[z,k][{a,b},bdot]};
angsqurepdot={Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;OrderedQ[{y,k}]&&y!=k:>AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;OrderedQ[{y,k}]==False:>-AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&OrderedQ[{x,z}]&&x!=z:>AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&OrderedQ[{x,z}]==False:>-AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&TrueQ[x==z]&&OrderedQ[{a,b}]:>AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}],
Times[AngSquInvariant[x_,y_][a_,adot_],AngSquInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&TrueQ[x==z]&&OrderedQ[{a,b}]==False:>-AngSquInvariant[x,y][a,{adot,bdot}]*AngSquInvariant[z,k][b,{adot,bdot}]};
squangrep={Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;OrderedQ[{x,z}]&&x!=z:>SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot],Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;OrderedQ[{x,z}]==False:>-SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot],Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&OrderedQ[{y,k}]&&y!=k:>SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&OrderedQ[{y,k}]==False:>-SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&TrueQ[y==k]&&OrderedQ[{adot,bdot}]:>SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[a_,b_]]/;TrueQ[x==z]&&TrueQ[y==k]&&OrderedQ[{adot,bdot}]==False:>-SquAngInvariant[x,y][{a,b},adot]*SquAngInvariant[z,k][{a,b},bdot]};
squangrepdot={Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;OrderedQ[{y,k}]&&y!=k:>SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;OrderedQ[{y,k}]==False:>-SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&OrderedQ[{x,z}]&&x!=z:>SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&OrderedQ[{x,z}]==False:>-SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&TrueQ[x==z]&&OrderedQ[{a,b}]:>SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}],
Times[SquAngInvariant[x_,y_][a_,adot_],SquAngInvariant[z_,k_][b_,bdot_],levicivita2Up[adot_,bdot_]]/;TrueQ[y==k]&&TrueQ[x==z]&&OrderedQ[{a,b}]==False:>-SquAngInvariant[x,y][a,{adot,bdot}]*SquAngInvariant[z,k][b,{adot,bdot}]};
locexp=exp//.Join[angsqurep,squangrep,angsqurepdot,squangrepdot];
(*Now the expression will be simplified and then we backsubstitute*)
backangsqu={Times[AngSquInvariant[x_,y_][{a_,b_},adot_],AngSquInvariant[z_,k_][{a_,b_},bdot_]]:>AngSquInvariant[x,y][a,adot]*AngSquInvariant[z,k][b,bdot]*levicivita2Up[a,b],Times[AngSquInvariant[x_,y_][a_,{adot_,bdot_}],AngSquInvariant[z_,k_][b_,{adot_,bdot_}]]:>AngSquInvariant[x,y][a,adot]*AngSquInvariant[z,k][b,bdot]*levicivita2Up[adot,bdot]};
backsquang={Times[SquAngInvariant[x_,y_][{a_,b_},adot_],SquAngInvariant[z_,k_][{a_,b_},bdot_]]:>SquAngInvariant[x,y][a,adot]*SquAngInvariant[z,k][b,bdot]*levicivita2Up[a,b],Times[SquAngInvariant[x_,y_][a_,{adot_,bdot_}],SquAngInvariant[z_,k_][b_,{adot_,bdot_}]]:>SquAngInvariant[x,y][a,adot]*SquAngInvariant[z,k][b,bdot]*levicivita2Up[adot,bdot]};
locexp=Expand[Simplify[locexp]];
locexp=locexp//.Join[backsquang,backangsqu];
Return[locexp];

];


(* ::Subsection::Closed:: *)
(*Computing F with Lorentz indices*)


FstBoxup[mom_,mu_,nu_]:=TemplateBox[{mom,mu,nu},"Fstup",
DisplayFunction->(RowBox[{SuperscriptBox["F",RowBox[{#2,#3}]],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Fst","[",#1,"][",#2,",",#3,"]","[","Null",",","Null","]"}]&)
];
FstBoxdown[mom_,mu_,nu_]:=TemplateBox[{mom,mu,nu},"Fstdown",
DisplayFunction->(RowBox[{SubscriptBox["F",RowBox[{#2,#3}]],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Fst","[",#1,"]","[","Null",",","Null","]","[",#2,",",#3,"]"}]&)
];
Fst /: MakeBoxes[Fst[mom_][mu_,nu_][Null,Null],StandardForm|TraditionalForm]:=FstBoxup[ToBoxes[mom],ToBoxes[mu],ToBoxes[nu]];
Fst /: MakeBoxes[Fst[mom_][Null,Null][mu_,nu_],StandardForm|TraditionalForm]:=FstBoxdown[ToBoxes[mom],ToBoxes[mu],ToBoxes[nu]];
Fst[mom_][mu_,nu_][Null,Null]:=Mom[mom][mu][Null]Polar[mom,ToExpression["ref"<>ToString[mom]]][nu][Null]-Mom[mom][nu][Null]Polar[mom,ToExpression["ref"<>ToString[mom]]][mu][Null];

EtaBoxup[mu_,nu_]:=TemplateBox[{mu,nu},"etaup",
DisplayFunction->(SuperscriptBox["\[Eta]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"Eta","[",#1,",",#2,"]","[","$up","]"}]&)
];
EtaBoxdown[mu_,nu_]:=TemplateBox[{mu,nu},"etadown",
DisplayFunction->(SubscriptBox["\[Eta]",RowBox[{#1,#2}]]&),
InterpretationFunction->(RowBox[{"Eta","[",#1,",",#2,"]","[","$down","]"}]&)
];
Eta /: MakeBoxes[Eta[mu_,nu_][$up],StandardForm|TraditionalForm]:=EtaBoxup[ToBoxes[mu],ToBoxes[nu]];
Eta /: MakeBoxes[Eta[mu_,nu_][$down],StandardForm|TraditionalForm]:=EtaBoxdown[ToBoxes[mu],ToBoxes[nu]];
SetAttributes[Eta,Orderless];

PolarBoxup[mom_,refmom_,mu_]:=TemplateBox[{mom,refmom,mu},"Polarup",
DisplayFunction->(RowBox[{SuperscriptBox["\[CurlyEpsilon]",#3],"[",#1,",",#2,"]"}]&),
InterpretationFunction->(RowBox[{"Polar","[",#1,",",#2,"][",#3,"]","[","Null","]"}]&)
];
PolarBoxdown[mom_,refmom_,mu_]:=TemplateBox[{mom,refmom,mu},"Polardown",
DisplayFunction->(RowBox[{SubscriptBox["\[CurlyEpsilon]",#3],"[",#1,",",#2,"]"}]&),
InterpretationFunction->(RowBox[{"Polar","[",#1,",",#2,"]","[","Null","]","[",#3,"]"}]&)
];
PolarPureBox[mom_,refmom_]:=TemplateBox[{mom,refmom},"PolarPure",
DisplayFunction->(RowBox[{"\[CurlyEpsilon]","[",#1,",",#2,"]"}]&),
InterpretationFunction->(RowBox[{"PolarPure","[",#1,",",#2,"]"}]&)
];
Polar /: MakeBoxes[Polar[mom_,refmom_][mu_][Null],TraditionalForm|StandardForm]:=PolarBoxup[ToBoxes[mom],ToBoxes[refmom],ToBoxes[mu]];
Polar /: MakeBoxes[Polar[mom_,refmom_][Null][mu_],TraditionalForm|StandardForm]:=PolarBoxdown[ToBoxes[mom],ToBoxes[refmom],ToBoxes[mu]];
PolarPure /: MakeBoxes[PolarPure[mom_,refmom_],StandardForm|TraditionalForm]:=PolarPureBox[ToBoxes[mom],ToBoxes[refmom]];

(*MomBoxup[mom_,mu_]:=TemplateBox[{mom,mu},"Momup",
DisplayFunction->(RowBox[{SuperscriptBox["p",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Mom","[",#1,"][",#2,"]","[","Null","]"}]&)
];
MomBoxdown[mom_,mu_]:=TemplateBox[{mom,mu},"Momdown",
DisplayFunction->(RowBox[{SubscriptBox["p",#2],"[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"Mom","[",#1,"]","[","Null","]","[",#2,"]"}]&)
];
MomPureBox[mom_]:=TemplateBox[{mom},"MomPure",
DisplayFunction->(RowBox[{"p","[",#1,"]"}]&),
InterpretationFunction->(RowBox[{"MomPure","[",#1,"]"}]&)
];*)

MomBoxup[mom_,mu_]:=TemplateBox[{mom,mu},"Momup",
DisplayFunction->(RowBox[{SuperscriptBox[#1,#2]}]&),
InterpretationFunction->(RowBox[{"Mom","[",#1,"][",#2,"]","[","Null","]"}]&)
];
MomBoxdown[mom_,mu_]:=TemplateBox[{mom,mu},"Momdown",
DisplayFunction->(RowBox[{SubscriptBox[#1,#2]}]&),
InterpretationFunction->(RowBox[{"Mom","[",#1,"]","[","Null","]","[",#2,"]"}]&)
];
MomPureBox[mom_]:=TemplateBox[{mom},"MomPure",
DisplayFunction->(RowBox[{#1}]&),
InterpretationFunction->(RowBox[{"MomPure","[",#1,"]"}]&)
];


Mom /: MakeBoxes[Mom[mom_][mu_][Null],StandardForm|TraditionalForm]:=MomBoxup[ToBoxes[mom],ToBoxes[mu]];
Mom /: MakeBoxes[Mom[mom_][Null][mu_],StandardForm|TraditionalForm]:=MomBoxdown[ToBoxes[mom],ToBoxes[mu]];
MomPure /: MakeBoxes[MomPure[mom_],StandardForm|TraditionalForm]:=MomPureBox[ToBoxes[mom]];

Eta /: Times[Eta[mu_,nu_][$down],f_[lab__][nu_][Null]]:=f[lab][Null][mu];
Eta /: Times[Eta[mu_,nu_][$up],f_[lab__][Null][nu_]]:=f[lab][mu][Null];
Mom /: Times[Mom[lab1_][mu_][Null],f_[lab2__][Null][mu_]]:=mp[MomPure[lab1],ToExpression[ToString[f]<>"Pure"][lab2]];
Polar /: Times[Polar[lab1_,lab2_][mu_][Null],f_[lab3__][Null][mu_]]:=mp[PolarPure[lab1,lab2],ToExpression[ToString[f]<>"Pure"][lab3]];
PolarPure /: mp[PolarPure[x_,ref_],PolarPure[y_,ref_]]:=0;
PolarPure /: mp[PolarPure[x_,ref_],MomPure[x_]]:=0;
PolarPure /: mp[PolarPure[x_,ref_],MomPure[ref_]]:=0;


(* ::Subsection::Closed:: *)
(*ToFile*)


SetAttributes[ToFile,HoldFirst];
ToFile[exp_,filename_]:=Block[{Print,olddir,stream,out},
olddir=Directory[];
SetDirectory[$UserDocumentsDirectory];
stream=OpenWrite[ToString[filename]<>".txt"];
Print[x_]:=Write[stream,x];
out=exp//ReleaseHold;
Close[stream];
SetDirectory[olddir];
Return[out];
];


(* ::Section::Closed:: *)
(*Shortcuts*)


(* ::Subsection:: *)
(*$Shortcuts*)


If[frontend==1,
$Shortcuts={RawBoxes[RowBox[{SuperscriptBox["\[Lambda]","\[Alpha]"],"[p]"}]] -> "esc + lu + esc",
			RawBoxes[RowBox[{SubscriptBox["\[Lambda]","\[Alpha]"],"[p]"}]] -> "esc + ld + esc",
			RawBoxes[RowBox[{SuperscriptBox[OverscriptBox["\[Lambda]","~"],"\[Alpha]"],"[p]"}]] -> "esc + ltu + esc",
			RawBoxes[RowBox[{SubscriptBox[OverscriptBox["\[Lambda]","~"],"\[Alpha]"],"[p]"}]] -> "esc + ltd + esc",
			RawBoxes[RowBox[{SuperscriptBox["\[Mu]","\[Alpha]"],"[p]"}]] -> "esc + muu + esc",
			RawBoxes[RowBox[{SubscriptBox["\[Mu]","\[Alpha]"],"[p]"}]] -> "esc + mud + esc",
			RawBoxes[RowBox[{SuperscriptBox[OverscriptBox["\[Mu]","~"],"\[Alpha]"],"[p]"}]] -> "esc + mtu + esc",
			RawBoxes[RowBox[{SubscriptBox[OverscriptBox["\[Mu]","~"],"\[Alpha]"],"[p]"}]] -> "esc + mtd + esc",
			RawBoxes[RowBox[{"\[Lambda][p]"}]]-> "esc + lp + esc",
			RawBoxes[RowBox[{"\[Mu][p]"}]]-> "esc + mp + esc",
			RawBoxes[RowBox[{OverscriptBox["\[Lambda]","~"],"[p]"}]]-> "esc + ltp + esc",
			RawBoxes[RowBox[{OverscriptBox["\[Mu]","~"],"[p]"}]]-> "esc + mtp + esc",
			RawBoxes[RowBox[{SubscriptBox[SuperscriptBox["\[CapitalLambda]","A"],"a"],"[p]"}]]-> "esc + Ld + esc",
			RawBoxes[RowBox[{SubscriptBox[SubscriptBox[OverscriptBox["\[CapitalLambda]","~"],"A"],"a"],"[p]"}]]-> "esc + Ltd + esc",
			RawBoxes[SuperscriptBox["\[Epsilon]",RowBox[{"a","b"}]]]->"esc + lcup + esc",
			RawBoxes[SubscriptBox["\[Epsilon]",RowBox[{"a","b"}]]]->"esc + lcd + esc",
			RawBoxes[SubscriptBox[SuperscriptBox["\[Delta]","A"],"B"]]-> "esc + md + esc",
			RawBoxes[RowBox[{"\[LeftAngleBracket]",RowBox[{"1 , 2"}],"\[RightAngleBracket]"}]] -> "esc + ab + esc",
			RawBoxes[RowBox[{"[",RowBox[{"1 , 2"}],"]"}]]-> "esc + sb + esc",
			RawBoxes[RowBox[{"\[LeftAngleBracket]",RowBox[{SubscriptBox["1","a"],",",SubscriptBox["2","b"]}],"]"}]] -> "esc + asi + esc",
			RawBoxes[RowBox[{"[",RowBox[{SubscriptBox["1","a"],",",SubscriptBox["2","b"]}],"\[RightAngleBracket]"}]] -> "esc + sai + esc",
			RawBoxes[RowBox[{"\[LeftAngleBracket]",SubscriptBox["1","a"],",",SubscriptBox["2","b"],",",SubscriptBox["3","c"],",",SubscriptBox["4","d"],"\[RightAngleBracket]"}]]-> "esc + aa + esc",
			RawBoxes[RowBox[{"[",SubscriptBox["1","a"],",",SubscriptBox["2","b"],",",SubscriptBox["3","c"],",",SubscriptBox["4","d"],"]"}]]-> "esc + bb + esc",
			RawBoxes[RowBox[{SubscriptBox[SuperscriptBox["F",RowBox[{"A","B"}]],RowBox[{"C","D"}]],"[",SubscriptBox["p",RowBox[{"a","b"}]],"]"}]] -> "esc + fst + esc"}//MatrixForm;
			];


(* ::Section::Closed:: *)
(*Attributes*)


SetAttributes[{overbar,MinusSignQ,SpinorDot,SpinorUndot,SpinorLaDownBox,SpinorLaUpBox,SpinorLatDownBox,
SpinorLatUpBox,SpinorMuDownBox,SpinorMuUpBox,SpinorMutDownBox,SpinorMutUpBox,Spinordot,Spinorundot,$lam,$mu,
extramass,extramasstilde,extramasstildeBox,extramassBox,ExtramassBox,ExtramasstildeBox,Extramasstilde,Extramass,
KillMasses,Momenta4D,SpinorAngleBracket,SpinorAngleBracketBox,SpinorSquareBracket,SpinorSquareBracketBox,
NewProcess,ClearDownValues,levicivita2up,levicivita2down,levicivita2Up,levicivita2Down,levicivita2up,levicivita2down,
SpinorDot6D,SpinorDot6DBox,SpinorUndot6D,SpinorUndot6DBox,AngAngInvariant,AngAngInvariantBox,SquSquInvariant,SquSquInvariantBox,
SquAngInvariant,SquAngInvariantBox,AngSquInvariant,AngSquInvariantBox,Momenta,AllMomenta,
SpinorReplace,SubCounter,ConvenientMu,SchoutenSimplify,Mom4D,Mom4DVector,S6,S6many,S6N,S6manyN,MDelta,MDeltaBox,Fstrength,FstrengthBox,\[Lambda],\[Mu],\[CapitalLambda],\[CapitalMu],\[Epsilon],\[Delta],FixedSpinors,ClearSpinors,FixSpinors,CompleteDenominators,CompleteMandelstam,ToChain,Chain,chain,$angle,$square,S},Protected]


(* ::Section::Closed:: *)
(*Create a palette*)


(*
If[frontend==1,
SpinorPalette:=CreatePalette[Grid[Join[Partition[PasteButton[RawBoxes[#]]&/@{SpinorLaUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLaDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLatUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLatDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMuUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMuDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMutUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMutDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorUndotPureLBox["\[SelectionPlaceholder]"],SpinorDotPureLBox["\[SelectionPlaceholder]"],SpinorUndotPureMBox["\[SelectionPlaceholder]"],SpinorDotPureMBox["\[SelectionPlaceholder]"],SpinorUndot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorDot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],levicivita2upBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],levicivita2downBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorAngleBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"],SpinorSquareBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"],SquAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],AngSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]},4],{{PasteButton[RawBoxes[AngAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]],PasteButton[RawBoxes[SquSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]],PasteButton[RawBoxes[RowBox[{"\[CapitalMu][","\[SelectionPlaceholder]","]"}]],Defer[extramass[\[SelectionPlaceholder]]]],PasteButton[RawBoxes[RowBox[{OverscriptBox["\[CapitalMu]","~"],"[","\[SelectionPlaceholder]","]"}]],Defer[extramasstilde[\[SelectionPlaceholder]]]]}}],Spacings->Automatic],WindowTitle->"SpinorHelicity6D"];
];


If[frontend==1,
CreatePalette[Grid[Join[Partition[PasteButton[RawBoxes[#]]&/@{SpinorLaUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLaDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLatUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorLatDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMuUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMuDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMutUpBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorMutDownBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorUndotPureLBox["\[SelectionPlaceholder]"],SpinorDotPureLBox["\[SelectionPlaceholder]"],SpinorUndotPureMBox["\[SelectionPlaceholder]"],SpinorDotPureMBox["\[SelectionPlaceholder]"],SpinorUndot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorDot6DBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],levicivita2upBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],levicivita2downBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]"],SpinorAngleBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"],SpinorSquareBracketBox["\[SelectionPlaceholder]", "\[Placeholder]"],SquAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"],AngSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]},4],{{PasteButton[RawBoxes[AngAngInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]],PasteButton[RawBoxes[SquSquInvariantBox["\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]","\[SelectionPlaceholder]"]]],PasteButton[RawBoxes[RowBox[{"\[CapitalMu][","\[SelectionPlaceholder]","]"}]],Defer[extramass[\[SelectionPlaceholder]]]],PasteButton[RawBoxes[RowBox[{OverscriptBox["\[CapitalMu]","~"],"[","\[SelectionPlaceholder]","]"}]],Defer[extramasstilde[\[SelectionPlaceholder]]]]}}],Spacings->Automatic],WindowTitle->"SpinorHelicity6D"];
];
 *)

Print["===============SpinorHelicity6D================"];
Print["Authors: Manuel Accettulli Huber (QMUL)"];
Print["         Stefano De Angelis (QMUL)"];
Print["Please report any bug to:"];
Print["m.accettullihuber@qmul.ac.uk"];
Print["Version 1.1 , last update 06/09/2019"];
Print["==============================================="];


If[TrueQ[$MachineID=="6239-87290-05914"],Speak["Ubi maior minor cessat"]];
(*If[TrueQ[$MachineID=="5113-13572-95048"],Speak["Tricche tracche bombe a mano"]];*)


End[]

EndPackage[]
