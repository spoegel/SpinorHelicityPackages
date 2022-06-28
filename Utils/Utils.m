
(* ::Package:: *)

(* ::Title:: *)
(* Package: Utils`*)


(* ::Chapter:: *)
(* Preamble*)


BeginPackage["Utils`"];

Unprotect["Utils`*"];
ClearAll["Utils`*"];
ClearAll["Utils`Private`*"];

StartProgressBar::usage = "Starts the progress indicator.

Example:

   StartProgressBar[Dynamic[k], 100];
   For[k = 1, k <= 100, k++,
       UpdateProgressBar[k, 100];
       DoSomething[];
      ];
   StopProgressBar[100];
";


(* ::Section:: *)
(*Exports*)


UpdateProgressBar::usage = "Updates the progress indicator.";

StopProgressBar::usage = "Stops the progress bar.";

PrimaryPermutations::usage = "removes elements that are linked by reversal and cyclicity";
PrimaryPermutationsNoReversal::usage = "removes elements that are linked by cyclic permutations";

IntersectionBasis::usage = "Finds a basis of the intersection of a list of generating sets of subspaces";

RationalReconstruct::usage = "";

CyclicPermutations::usage = "Generate all cyclic permutations of a list.";
CyPerms::usage = "(DEPRECATED, use CyclicPermutations) Generate all cyclic permutations of a list.";

BellList::usage = "";

AllFreeQ::usage = "AllFreeQ[expr,List[...]] returns True, if expr is FreeQ of all elements in List.";

ResidueDiscreteFourier::usage = "Computes residues numerically via discrete fourier projection.";

ExtractSymbolName::usage = "Strips all contexts and uniqueness postfixes from symbols.";

PassRules::usage = "Filters rules for passing options from one function to another.";

ShuffleProduct::usage = "Implementation of the shuffle product of two (or several) lists.";

AndPattern::usage = "";

shead::usage = "";
headPattern::usage = "";


heavySymbols::usage = "";

printTestResults::usage = "";

NumericSameQ::usage = "Behaves like SameQ, with the additional condition that all arguments are
numeric.";


SplitAt::usage = "";

SubsetSplit::usage = "";


RotateToFirst::usage = "";

GetPivots::usage = "";

GetBasis::usage = "";
GetBasisPositions::usage = "";

GetQuotientBasis::usage = "";
GetQuotientBasisPositions::usage = "";

CycleCount::usage = "";

PickCycles::usage = "";

SortedTally::usage = "";

RotateTo::usage = "";

EqualIntegerPartitions::usage = "";
EqualListPartition::usage = "";

GenerateOptimizedEvaluationFunctionFF::usage = "";
GenerateOptimizedEvaluationFunction::usage = "";

IncreasingIntegerSequence::usage = "";

IntegerNullspace::usage = "";

RandomIntegerMatrixRankDeficient::usage = "";

SymmetricTensorProductElements::usage = "SymmetricTensorProductElements[l,r]: Returns the unique elements in the symmetrized rank `r` tensor product of list `l`.";
SymmetricTensorProductSize::usage = "SymmetricTensorProductSize[l,r]: Returns the number of unique elements in the symmetrized rank `r` tensor product of list `l`.";
SymmetricTensorProductSize::usage = "SymmetricTensorProductSize[n,i]: Returns the number of unique elements in the symmetrized rank `r` tensor product of `n` elements.";


DistributeList::usage = "DistributeList[list, n]: Generates all possible ways of splitting list into n smaller lists.";

(* ::Chapter:: *)
(* Private*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(* Progress Bar *)

(*
  Source: FlexibleSUSY
  Licensed under GPL-3.0
  Link: https://github.com/FlexibleSUSY/FlexibleSUSY
*)

StartProgressBar[dyn:Dynamic[x_], total_, len_:50] :=
    If[$Notebooks,
       PrintTemporary @ Row[{ProgressIndicator[dyn, {0, total}], " Total: ", total}]
      ];

UpdateProgressBar[x_, total_, len_:50] :=
  If[!$Notebooks,
    WriteString["stdout", "[" <> StringJoin[
        Join[Table[".",{i,Round[len*x/total]}],
             Table[" ",{i,Round[len*(1-x/total)]}]]
       ] <> "] " <> ToString[x] <> "/" <> ToString[total] <> "\r"];
   ];

StopProgressBar[total_, len_:50] :=
    If[!$Notebooks,
       WriteString["stdout", "[" <> StringJoin[
           Table[".",{i,Round[len]}]
       ] <> "] " <> ToString[total] <> "/" <> ToString[total] <> "\n"];
      ];


(* ::Section::Closed:: *)
(* Intersection Basis *)

(*

  Author: Wile E.
  Modified by: Sebastian Poegel
  Licensed under CC BY-SA 4.0
  Link: https://mathematica.stackexchange.com/q/64165/73380
*)

IntersectionBasis[prime_Integer] := {};
IntersectionBasis[{},prime_Integer] := {};
IntersectionBasis[{}, __,prime_Integer] := {};
IntersectionBasis[__, {},prime_Integer] := {};
IntersectionBasis[l1_,prime_Integer] := IntersectionBasis[l1, l1,prime];
IntersectionBasis[l1_, l2_, l3__, prime_Integer] :=
  IntersectionBasis[IntersectionBasis[l1, l2,prime], l3,prime];
IntersectionBasis[l1_, l2_,prime_Integer] :=
  Catch[With[{ker = NullSpace[Transpose[Join[l1, l2]],Modulus->prime]},
    If[ker === {}, Throw[{}],
      DeleteCases[RowReduce[ker[[All, 1 ;; Length[l1]]] . l1 ,Modulus->prime], {0 ..}]]]];



(* ::Section::Closed:: *)
(* Cyclic Permutations *)
(*
   Author: ciao
   Modified by: Sebastian Poegel
   Licensed under CC BY-SA 3.0
   Link: https://mathematica.stackexchange.com/a/120593/73380

*)

CyclicPermutations[{}] := {{}};
CyclicPermutations[n_Integer] := HankelMatrix[#,RotateRight@#]& @ Range[n];
CyclicPermutations[x_List] := Map[x[[#]] &, HankelMatrix[#, RotateRight@#] &@Range[Length@x], {2}];
CyPerms[x_] := CyclicPermutations[x];



(* ::Section::Closed:: *)
(* Primary Permutations *)

(*
  Author: Mr. Wizard
  Licensed under CC BY-SA 3.0
  Link: https://mathematica.stackexchange.com/a/28875/73380
*)

(* Deletes Cyclic Duplicates. Taken from mathematica stackexchange. *)
PrimaryPermutations[a_List] := Module[{f1, f2},
   f1 = RotateLeft[#, Ordering[#, 1] - 1] &;
   f2 = (# ~Extract~ Ordering[#, 1] &[f1 /@ {#, Reverse@#}] )&;
   DeleteDuplicatesBy[a,f2]
];


(* Deletes Cyclic Duplicates. Taken from mathematica stackexchange. *)
PrimaryPermutationsNoReversal[a_List] := Module[{f1, f2},
   f1 = RotateLeft[#, Ordering[#, 1] - 1] &;
   f2 = (# ~Extract~ Ordering[#, 1] &[f1 /@ {#, #}] )&;
   DeleteDuplicatesBy[a,f2]
];


(* ::Section::Closed:: *)
(*BellList*)

(*
  Author: Robert M. Dickau
  Link: https://web.archive.org/web/20190331234613/http://mathforum.org/advanced/robertd/bell.html
*)


BellList[1] = {{{1}}}; (* the basic case *)

BellList[n_Integer?Positive] := BellList[n] = (* do some caching *)
    Flatten[
      Join[
        Map[ReplaceList[#,{S__}->{S,{n}}]&, BellList[n-1]],

        Map[ReplaceList[#,{b___,{S__},a___}->{b,{S,n},a}]&, BellList[n-1]]],
      1];


(* ::Section::Closed:: *)
(* AllFreeQ *)


AllFreeQ[expr_,form_List]:=And@@(FreeQ[expr,#]&/@form);


(* ::Section::Closed:: *)
(*Rational Reconstruct*)

(* Author: Sebastian Poegel *)

SetAttributes[RationalReconstruct,Listable];

RationalReconstruct[a_, n_Integer] := If[a === 0, 0, Divide @@ (Reverse@LatticeReduce[{{1, a}, {0, n}}][[1]])];
RationalReconstruct[n_Integer] := Function[a, RationalReconstruct[a, n]];



(* ::Section::Closed:: *)
(* Residue Discrete Fourier *)

(* Author: Sebastian Poegel *)

Options[ResidueDiscreteFourier]={PrecisionGoal->MachinePrecision};

ResidueDiscreteFourier[expr_,{var_,p_,order_,t0_},OptionsPattern[]] := Module[{},
   1/(2*p+1)*SetPrecision[Sum[Divide@@(({expr,var^order})/.{var->t0*Exp[I*2*Pi*j/(2*p+1)]}),{j,-p,p}],OptionValue[PrecisionGoal]]
];


(* ::Section::Closed:: *)
(* Extract Symbol Name *)


(* Removes contexts, and strips unique postfixes *)
(*
  Author: Mr.Wizard
  Licensed under CC BY-SA 3.0
  Link: https://mathematica.stackexchange.com/a/9099
*)

ExtractSymbolName[expr_]:=Module[{T, SR = StringReplace[#, a__ ~~ "$" ~~ DigitCharacter .. :> a] &},
  expr
    /. s_Symbol :> T @ MakeExpression @ SR @ SymbolName @ Unevaluated @ s
    /. T@_@x___ :> x
];



(* ::Section::Closed:: *)
(*PassRules*)

(* Author: Sebastian Poegel *)

PassRules[fctFrom_,fctTo_,opts__] := FilterRules[Join[Flatten@{opts},Options[fctFrom]],Options[fctTo]];
PassRules[fctFrom_,fctTo_] := FilterRules[Options[fctFrom],Options[fctTo]];
PassRules[fct_,fct_] := {};


(* ::Section::Closed:: *)
(*Shuffle Product*)

(* Author: Sebastian Poegel *)

ShuffleProduct[list1_]:={list1};
ShuffleProduct[list1_List,list2_List,listsRest__]:=ShuffleProduct[#,listsRest]&/@ShuffleProduct[list1,list2]//Apply[Join];
ShuffleProduct[list1_List,list2_List]:=Fold[ShuffleProductInsertUpToDepth,{{list1,Length@list1}},list2][[All,1]];
ShuffleProductInsertUpToDepth[listsAndDepths_,newElement_]:=ShuffleProductInsertElementToSingleList[newElement]@@@listsAndDepths//Apply[Join];
ShuffleProductInsertElementToSingleList[newElement_]:=Function[{list,maxDepth},{Insert[list,newElement,-(#+1)],#}&/@Range[0,maxDepth]];


(* ::Section::Closed:: *)
(* AndPattern *)

(* Author: Sebastian Poegel *)

(* Returns a meta pattern that requires matching all patterns given as arguments *)

AndPattern[patts___] := Except[Alternatives @@ Map[Except, {patts}]];

(* ::Section::Closed:: *)
(*SplitAt*)
(* Author: Sebastian Poegel *)


SplitAt[{},___]:={};
SplitAt[list_List,{position_Integer,positionsRest__Integer}]/;!OrderedQ[{position,positionsRest}]:=
	SplitAt[list,Sort[{position,positionsRest}]];


SplitAt[list_List,{lastPos__Integer}] /; Length@list === lastPos :=
    {list, Splice@ConstantArray[{},Length@{lastPos}]};

SplitAt[list_List,{position_Integer,positionsRest__Integer}] /; OrderedQ[{position,positionsRest}]:=
	{list[[;;position]],Splice@SplitAt[list[[position+1;;]],{positionsRest}-position]};

	
SplitAt[list_List,{position_Integer}]:=Replace[{list[[;;position]],list[[position+1;;]]},{{a__},{}}:>{{a}}];
SplitAt[list_List,position_Integer]:=SplitAt[list,{position}];


(* ::Section:: *)
(*shead & headPattern *)


(* ::Text:: *)
(*
  Author: Leonid Shifrin
  Licensed under CC BY-SA 3.0
  Link: https://stackoverflow.com/a/5961851
*)


SetAttributes[shead, HoldAllComplete];
shead[expr_] := Scan[Return, Unevaluated[expr], {-1}, Heads -> True];


headPattern[head_] := _?(Function[Null, shead[#] === head, HoldFirst]);


(* ::Section:: *)
(*Memory Use Analysis*)

(*
  Author: Leonid Shifrin
  Licensed under CC BY-SA 4.0
  Link: https://mathematica.stackexchange.com/a/6179/73380
*)


$globalProperties =
    {OwnValues, DownValues, SubValues, UpValues, NValues, 
     FormatValues, Options, DefaultValues, Attributes, Messages};


SetAttributes[getDefinitions, HoldAllComplete];
getDefinitions[s_Symbol] :=
    Flatten@Through[
        Map[Function[prop, Function[sym, prop[sym], HoldAll]],
           $globalProperties
        ][Unevaluated[s]]
      ];

symbolMemoryUsage[sname_String] :=
   ToExpression[sname, InputForm, 
      Function[s, ByteCount[getDefinitions[s]], HoldAllComplete]
   ];

heavySymbols[context_, sizeLim_: 10^6] :=
   Pick[#, UnitStep[# - sizeLim] &@Map[symbolMemoryUsage, #], 1] &@
        Names[context <> "*"];

(* ::Section:: *)
(*Print Test Results*)

(*
  Author: Pinti
  Licensed under CC BY-SA 4.0
  Link: https://mathematica.stackexchange.com/a/176573/73380
*)


getTestResults[tr_TestReportObject]:=Module[
  {fields,results,abbreviations},

  (* Add other querries to this list. *)
  fields={"TestIndex","Outcome","AbsoluteTimeUsed","MemoryUsed","TestID"};
  abbreviations={"TestIndex"->"Idx","AbsoluteTimeUsed"->"Time [s]"};

  results=ReplaceRepeated[
    Outer[#1[#2]&,Values[tr["TestResults"]],fields],
    {x_Quantity:>QuantityMagnitude[x],x_Real:>Round[x,0.001]}
  ];

  Join[{fields/.abbreviations},results]
];

printTestResults[tr_TestReportObject]:=Module[
  {list,indx,time,noTests},

  list=getTestResults[tr];
  indx=MapIndexed[
    If[
      MemberQ[#1,"Failure"|"MessagesFailure"|"Error"],
      First[#2],
      Nothing
    ]&,
    list
  ];
  time=Round[QuantityMagnitude[tr["TimeElapsed"]],0.01];
  noTests=Length[tr["TestResults"]];

  Print[noTests," tests run in ",time," seconds."];

  If[
    TrueQ@tr["AllTestsSucceeded"],
    Print["All tests succeeded!"],
    Print[tr["TestsFailedCount"]," tests failed!"];
  ];

  Print@Grid[list,
    Alignment->Left,
    Dividers->{None,{2->True}},
    Background->{None,Thread[indx->Pink]}
  ];
  tr
];

(* ::Section:: *)
(* NumericSameQ *)
(* Author: Sebastian Poegel *)

NumericSameQ[x__] := And@@NumericQ/@{x} && SameQ[x];

(* ::Section:: *)
(* SubsetSplit *)
(* Author: Sebastian Poegel *)

SubsetSplit[list_List,opts___]:={#,Complement[list,#]}&/@Subsets[list,opts];

(* ::Section:: *)
(* RotateToFirst *)
(* Author: Sebastian Poegel *)

RotateToFirst[list_List] := RotateLeft[list, Ordering[list, 1] - 1];

(* ::Section:: *)
(* Quotient Basis and Tools *)
(* Author: Sebastian Poegel *)

GetPivots[matrixRR_]:=FirstPosition[#,1,Nothing]&/@matrixRR//Flatten;

Options[GetBasis] = {Modulus->0};
GetBasis[matrix_,opts:OptionsPattern[]]:=Extract[matrix,
  GetBasisPositions[matrix,PassRules[GetBasis,GetBasisPositions,opts]]
];

Options[GetBasisPositions] = {Modulus->0};
GetBasisPositions[matrix_,opts:OptionsPattern[]]:=
  FirstPosition[#, 1, Nothing] & /@
      RowReduce[Transpose@matrix,PassRules[GetBasis,RowReduce,opts]];

Options[GetQuotientBasisPositions] = {Modulus->0};
GetQuotientBasisPositions[basisSubspace_,basisSpace_,opts:OptionsPattern[]]:=
    Module[{rr},
      rr=RowReduce[Transpose@Join[basisSubspace,basisSpace],PassRules[GetQuotientBasisPositions,RowReduce,opts]];
      (GetPivots@rr-Length@basisSubspace)//.{_?(!Positive[#]&):>Nothing}
    ];

Options[GetQuotientBasis] = Options[GetQuotientBasisPositions];
GetQuotientBasis[basisSubspace_,basisSpace_,opts:OptionsPattern[]]:=
  Extract[
    basisSpace,
    List/@GetQuotientBasisPositions[
      basisSubspace,
      basisSpace,
      PassRules[GetQuotientBasis,GetQuotientBasisPositions,opts]
    ]
  ];

(* ::Section:: *)
(* Cycle functions *)
(* Author: Sebastian Poegel *)

PickCycles[n_,c_]:=Select[Permutations[Range@n],CycleCount[#]==c &];

CycleCount[perm_]:=(Length@@PermutationCycles[perm]+Length[perm]-Length@Flatten[PermutationCycles[perm][[1]]]);

(* ::Section:: *)
(* SortedTally *)
(* Author: Sebastian Poegel *)

SortedTally[args_]:=SortBy[First]@Tally[args];

(* ::Section:: *)
(* RotateTo *)
(* Author: Sebastian Poegel *)

RotateTo[x_List, el_] := RotateLeft[x, FirstPosition[x, el] - 1];
RotateTo[el_] := RotateTo[#, el] &;

(* ::Section:: *)
(* EqualIntegerPartitions *)
(* Author: Sebastian Poegel *)
EqualIntegerPartitions[n_Integer,partitions_Integer]:=
    ConstantArray[Quotient[n,partitions],partitions]+Join[ConstantArray[1,Mod[n,partitions]],
      ConstantArray[0,partitions-Mod[n,partitions]]];

(* ::Section:: *)
(* EqualListPartition *)
(* Author: Sebastian Poegel *)
EqualListPartition[list_List,n_Integer]:=
    SplitAt[list,FoldList[Plus,EqualIntegerPartitions[Length@list,n]]];


(* ::Section:: *)
(* GenerateOptimizedEvaluationFunction *)
(* Author: Sebastian Poegel *)

GenerateOptimizedEvaluationFunctionFF[basis_,parameters_List,prime_Integer]:=
Module[{basisMod,pars,C$},
  basisMod = Quiet[
    basis //.{
      Power[a_,x_]:>PowerMod[a,x,prime],
      HoldPattern[Rational[a_,b_]]:>a*PowerMod[b,-1,prime]
    }/.MapThread[Rule,{parameters,Table[Indexed[pars,i],{i,Length@parameters}]}],
    {PowerMod::pmod}
  ];

  Experimental`OptimizeExpression[basisMod,
    "OptimizationSymbol"->C$,
    "OptimizationLevel"->2
  ]/.{
    Experimental`OptimizedExpression[x__]:>
        Function @@ Hold[{pars},PolynomialMod[x,prime]]
  }
];

GenerateOptimizedEvaluationFunction[expr_,parameters_List]:=
    Module[{pars,exprReparametrized},
      exprReparametrized =
          expr /. MapThread[Rule,{parameters,Table[Indexed[pars,i],{i, Length@parameters}]}];

      Experimental`OptimizeExpression[
        exprReparametrized,
        "OptimizationSymbol"->C$,
        "OptimizationLevel"->2
      ]/.{
        Experimental`OptimizedExpression[x__]:>
            Function @@ Hold[{pars},x]
      }
    ];

(* ::Section:: *)
(* IncreasingIntegerSequence *)
(* Author: Sebastian Poegel *)
IncreasingIntegerSequence[length_Integer,maxInt_Integer] :=
    IncreasingIntegerSequence[length,1,maxInt];
IncreasingIntegerSequence[length_Integer,minInt_Integer,maxInt_Integer]:=Module[{vars,bounds},
  vars=Table[Unique["var"],{length}];
  bounds = Prepend[{#2,#1,maxInt}&@@@Partition[vars,2,1],{First@vars,1,maxInt}];
  Flatten[Table@@{vars,Sequence@@bounds},length-1]
];


(* ::Section:: *)
(* IntegerNullspace *)
(* Author: Sebastian Poegel *)

IntegerNullspace[list_List]/;Length@list>=2:=
    BuildIntegerNullspace[{},list,Quiet[FindIntegerNullVector[list,30],FindIntegerNullVector::norel]];
IntegerNullspace[___]:={};

BuildIntegerNullspace[nvs_List,list_List,newnv_List]:=
    Join[
      nvs,
      {newnv},
      Insert[#,0,FirstPosition[newnv,Except[0,_Integer]]]&/@
          IntegerNullspace[Delete[list,FirstPosition[newnv,Except[0,_Integer]]]]
    ]
BuildIntegerNullspace[nvs_List,list_List,_]:=nvs;

(* ::Section:: *)
(* IntegerNullspace *)
(* Author: Sebastian Poegel *)
RandomIntegerMatrixRankDeficient[n_Integer,rank_Integer,p_Integer]:=
    Mod[RandomInteger[{1,p},{n, rank}] . RandomInteger[{1,p},{rank,n}],p];

(* ::Section:: *)
(* SymmetricTensorProductElements *)
(* Author: Sebastian Poegel *)
SymmetricTensorProductElements[elements_List,tensorRank_Integer]:=
    Times@@@Map[elements[[#]]&,IncreasingIntegerSequence[tensorRank,Length@elements],{-1}];

SymmetricTensorProductSize[n_Integer,rank_Integer]:=1/rank!*Pochhammer[n,rank];
SymmetricTensorProductSize[els_List,rank_Integer]:=SymmetricTensorProductSize[Length@els,rank];


(* ::Section:: *)
(* DistributeList *)
(* Author: Sebastian Poegel *)
DistributeList[{},nBoxes_Integer]:=
    {ConstantArray[{},nBoxes]};

DistributeList[list_List,nBoxes_Integer]:=
    SplitAt[list,#]&/@(IncreasingIntegerSequence[nBoxes-1,1+Length@list]-1)

(* ::Chapter:: *)
(* Closing*)


End[];


EndPackage[];
