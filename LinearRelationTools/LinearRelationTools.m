(*

    This file is part of SpinorHelicityAddons.

    SpinorHelicityAddons is free software: you can redistribute it and/or modify
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

(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: SpinorHelicity6D.LinearRelationTools *)
(* :Context: SpinorHelicity6D.LinearRelationTools` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-09-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage[
  "LinearRelationTools`",
  {
    "MomentumTwistorDefinitions`",
    "Utils`"
  }
];
(* Exported symbols added here with SymbolName::usage *)

GenerateNullSpaceFFMatrix::usage = "";
GenerateNullSpaceFFMatrixNTimes::usage = "";
GenerateTensoredNullSpaceFFMatrixNTimes::usage = "";
GenerateTensoredNullSpaceFFMatrix::usage = "";

GetOptimizedSubbase::usage = "";

FFEvaluation::usage = "";

DeTensorRelations::usage = "";

TwistorRandomEval::usage = "";

TwistorRandomEvalNTimes::usage = "";

FindIndependentRelations::usage = "";

$LinearRelationsPrime::usage = "";
$LinearRelationsPrime = NextPrime[2^24,-1];

SymmetricTensorProduct::usage = "SymmetricTensorProduct[basis,i]: Symbol used for
GenerateTensoredNullSpaceFFMatrixNTimes, which
represents the rank i symmetric tensor product of basis.";

Begin["`Private`"];

$linearRelationsCachedExpression;

(* ::Section:: *)
(* Generating nullSpace matrices *)

GenerateNullSpaceFFMatrix[basis_]:=
    Module[{tempbasis=basis},
      DistributeDefinitions[tempbasis];
      Map[Mod[#,$LinearRelationsPrime]&,
        ParallelTable[
            TwistorRandomEval[tempbasis]/.
                {Power[x_,y_]:>PowerMod[x,y,$LinearRelationsPrime],Rational[x_,y_]:>x*PowerMod[y,-1,$LinearRelationsPrime]},
          {Length@basis}, DistributedContexts->None, Method->"CoarsestGrained"],
      -1]
    ];

GenerateNullSpaceFFMatrixNTimes[basis_,n_Integer]:=
    Module[{tempbasis=basis},
      DistributeDefinitions[tempbasis];

      Map[Mod[#,$LinearRelationsPrime]&,
        ParallelTable[

          TwistorRandomEval[tempbasis]/.
              {Power[x_,y_]:>PowerMod[x,y,$LinearRelationsPrime],Rational[x_,y_]:>x*PowerMod[y,-1,$LinearRelationsPrime]},
          {n}, DistributedContexts->None, Method->"CoarsestGrained"],
        -1]
    ];

(*GenerateTensoredNullSpaceFFMatrixNTimes[{bases__},n_Integer]:=*)
(*    Module[{tempbases={bases}},*)
(*      DistributeDefinitions[tempbases];*)

(*      Map[Mod[#,$LinearRelationsPrime]&,*)
(*        ParallelTable[*)
(*          Times@@@Tuples[TwistorRandomEval[tempbases]]/.*)
(*              {Power[x_,y_]:>PowerMod[x,y,$LinearRelationsPrime],Rational[x_,y_]:>x*PowerMod[y,-1,$LinearRelationsPrime]},*)
(*          {n}, DistributedContexts->None, Method->"CoarsestGrained"],*)
(*        -1]*)
(*    ];*)

$linearRelationsCacheDirectory = FileNameJoin[
  {
    HomeDirectory[],
    ".cache",
    "Mathematica",
    "Packages",
    "LinearRelationTools"
  }
];

optimizedExpressionCachedQ[expr_] /; FileExistsQ @
    FileNameJoin[{$linearRelationsCacheDirectory, ToString@Hash@expr}] := True;

optimizedExpressionCachedQ[___] := False;



loadCachedOptimizedExpression[expr_] /; optimizedExpressionCachedQ[expr] := Module[{hash, file},
  hash = Hash@expr;
  file = FileNameJoin[{$linearRelationsCacheDirectory, ToString@hash}];

  Clear[$linearRelationsCachedExpression];
  file//Get;
  (Clear[$linearRelationsCachedExpression]; #) & [$linearRelationsCachedExpression]
];

loadCachedOptimizedExpression[___] := $Failed;

cacheOptimizedExpression[Rule[expr_,exprOpt_]] := Module[{hash, file},
  hash = Hash@expr;
  file = FileNameJoin[{$linearRelationsCacheDirectory, ToString@hash}];

  If[!DirectoryQ@DirectoryName@file,
    CreateDirectory@DirectoryName@file
  ];

  Clear[$linearRelationsCachedExpression];
  $linearRelationsCachedExpression = exprOpt;
  DumpSave[file,$linearRelationsCachedExpression];
  Clear[$linearRelationsCachedExpression];
];





tensoredBasisLength[bases_,extras_]:=
    Length@extras + Times@@Replace[bases,{list:{___}:>Length@list,
      SymmetricTensorProduct[args___]:>SymmetricTensorProductSize[args]},{1}];

Options[GenerateTensoredNullSpaceFFMatrixNTimes] =
    {Indent->0, RebuildCache -> False, ExtraElements -> {}};
GenerateTensoredNullSpaceFFMatrixNTimes[bases:{({__}|_SymmetricTensorProduct)..},n_Integer,
  nMultiplicity_Integer, opts:OptionsPattern[]]:=
    Module[{basesOptimized,extraElementsOptimized},

      basesOptimized = GetOptimizedSubbase[
        #,
        nMultiplicity,
        PassRules[
          GenerateTensoredNullSpaceFFMatrixNTimes,
          GetOptimizedSubbase,
          opts
        ]
      ]& /@ bases;

      extraElementsOptimized = GetOptimizedSubbase[
        OptionValue[ExtraElements],
        nMultiplicity,
        PassRules[
          GenerateTensoredNullSpaceFFMatrixNTimes,
          GetOptimizedSubbase,
          opts
        ]
      ];

      Print[StringRepeat["\t",OptionValue[Indent]],"Distributing optimized expressions..."];
      DistributeDefinitions[basesOptimized,extraElementsOptimized];
      DistributeDefinitions[SymmetricTensorProductElements];

      Print[StringRepeat["\t",OptionValue[Indent]],"Computing matrix..."];

      ParallelTable[
        Mod[#,$LinearRelationsPrime]&/@ (
          Join[extraElementsOptimized @ #,
            Times @@@ Tuples[Through[basesOptimized@#]/.SymmetricTensorProduct->SymmetricTensorProductElements]
          ]& @ RandomInteger[{0,$LinearRelationsPrime -1},Length @ MomentumTwistorParameters[nMultiplicity]]
        ),
        {n},
        DistributedContexts->None,
        Method->"CoarsestGrained"
      ]
    ];

Options[GenerateTensoredNullSpaceFFMatrix] =
    Options[GenerateTensoredNullSpaceFFMatrixNTimes];
GenerateTensoredNullSpaceFFMatrix[bases:{({__}|_SymmetricTensorProduct)..},
  nMultiplicity_Integer,  opts:OptionsPattern[]]:=
    GenerateTensoredNullSpaceFFMatrixNTimes[
      bases,
      tensoredBasisLength[
        bases,
        OptionValue[ExtraElements]
      ],
      nMultiplicity,
      opts
    ];





Options[GetOptimizedSubbase] = Options[GenerateTensoredNullSpaceFFMatrixNTimes];
GetOptimizedSubbase[subbase_,nMultiplicity_Integer,opts:OptionsPattern[]]:=
    If[optimizedExpressionCachedQ@subbase && !OptionValue[RebuildCache],
      Print[StringRepeat["\t",OptionValue[Indent]],"Found cached optimized expression, loading..."];
      loadCachedOptimizedExpression[subbase]
      ,
      Print[StringRepeat["\t",OptionValue[Indent]],"Generating optimized expression..."];
      (cacheOptimizedExpression[subbase->#]; #) & @
          GenerateOptimizedEvaluationFunctionFF[
            subbase,
            MomentumTwistorParameters[nMultiplicity],
            $LinearRelationsPrime
          ]
    ];


(* ::Section:: *)
(* FFEvaluation *)

FFEvaluation[expr_]:=Mod[
    expr/. {
      Power[x_,y_]:>PowerMod[x,y,$LinearRelationsPrime],
      Rational[x_,y_]:>x*PowerMod[y,-1,$LinearRelationsPrime]
    },
    $LinearRelationsPrime
];

(* ::Section:: *)
(* DeTensorRelations *)


DeTensorRelations[nullspace_List,basisAmps_,basisS_]:=
    Module[{expandedRelations,collected},
      expandedRelations = nullspace . Times@@@Tuples[{basisAmps,basisS}];
      ParallelMap[
        (Function[collected,
          Table[
            Coefficient[
              collected,
              amp
            ],
            {amp,basisAmps}
          ]
        ]@Collect[#, basisAmps])&,
        expandedRelations,
        DistributedContexts->None,
        Method->"CoarsestGrained"
      ]
    ];


(* ::Section:: *)
(* TwistorRandomEval *)


setupTwistorRandomEvalCache[]:=
    Module[{},
      ClearAll[TwistorRandomEvalCache];

      TwistorRandomEvalCache[Verbatim[Rational][a_,b_]] :=
          TwistorRandomEvalCache[Rational[a,b]] = a*PowerMod[b,-1,$LinearRelationsPrime];
      TwistorRandomEvalCache[Verbatim[Plus][expr___]] :=
          TwistorRandomEvalCache[Plus[expr]] = Plus@@TwistorRandomEvalCache/@{expr};
      TwistorRandomEvalCache[Verbatim[Times][expr___]] :=
          TwistorRandomEvalCache[Times[expr]] = Times@@TwistorRandomEvalCache/@{expr};
      TwistorRandomEvalCache[Verbatim[List][expr___]] :=
          TwistorRandomEvalCache[List[expr]] = List@@TwistorRandomEvalCache/@{expr};
      TwistorRandomEvalCache[Verbatim[Power][expr_,pow_]] :=
          TwistorRandomEvalCache[Power[expr,pow]] = PowerMod[TwistorRandomEvalCache@expr,pow, $LinearRelationsPrime];
      TwistorRandomEvalCache[Verbatim[Dot][expr___]] :=
          TwistorRandomEvalCache[Dot[expr]] = Dot@@TwistorRandomEvalCache/@{expr};

      TwistorRandomEvalCache[x_aPar] := TwistorRandomEvalCache[x] = RandomInteger[{1,$LinearRelationsPrime-1}];
      TwistorRandomEvalCache[x_bPar] := TwistorRandomEvalCache[x] = RandomInteger[{1,$LinearRelationsPrime-1}];
      TwistorRandomEvalCache[x_cPar] := TwistorRandomEvalCache[x] = RandomInteger[{1,$LinearRelationsPrime-1}];

      TwistorRandomEvalCache[x_] := x;
    ];

TwistorRandomEval[expr_]:=Module[{},
  setupTwistorRandomEvalCache[];
  Mod[TwistorRandomEvalCache@expr,$LinearRelationsPrime]
];

TwistorRandomEvalNTimes[basis_,n_Integer]:=
    Table[
      TwistorRandomEval[basis]/.
          {Power[x_,y_]:>PowerMod[x,y,$LinearRelationsPrime],Rational[x_,y_]:>x*PowerMod[y,-1,
            $LinearRelationsPrime]},
      {n}
    ];


(* ::Section:: *)
(* AmpToPrimary *)

AmpToPrimary[a_[args__]] /; Ordering[RotateToFirst /@ {{args}, Reverse@{args}}, 1]==={2}:= (-1)^(Length@{args})a@@RotateToFirst@Reverse@{args}
AmpToPrimary[a_[args__]] :=a@@RotateToFirst@{args};


(* ::Section:: *)
(* FindIndependentRelations *)

Options[FindIndependentRelations]={Modulus->$LinearRelationsPrime};
FindIndependentRelations[matrix_?NumericMatrixQ,opts:OptionsPattern[]]:=
    Fold[
      AddIndependentRelation[
        matrix,
        #1,
        #2,
        PassRules[FindIndependentRelations,AddIndependentRelation,opts]
      ]&,
      {},
      Range@Length@matrix
    ]

Options[AddIndependentRelation]=Options[FindIndependentRelations];
AddIndependentRelation[matrix_?NumericMatrixQ,indepRelations_List,newRelationIdx_Integer,opts:OptionsPattern[]]:=
    If[
      MatrixRankFFLAS[
        matrix[[Append[indepRelations,newRelationIdx]]],
        OptionValue[Modulus]
      ] > Length@indepRelations,
      Append[indepRelations,newRelationIdx],
      indepRelations
    ]


End[]; (* `Private` *)

EndPackage[]