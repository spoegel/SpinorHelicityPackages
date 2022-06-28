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

(* :Title: MomentumTwistorDefinitions *)
(* :Context: MomentumTwistorDefinitions` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-08-18 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["MomentumTwistorDefinitions`"];
(* Exported symbols added here with SymbolName::usage *)

LambdaTildeFromTwistor::usage = "Generates the LambdaTilde (bracket spinors) from twistor
variables.";

LambdaLambdaTildeFromTwistor::usage = "Generates both angle and bracket spinors
  from twistor matrix. Returns list of two lists, the first containing angle and the second
  the bracket spinors. In either case the associated SU(2) index is down.
  The spinor products are defined as <ab> = lambda[up].lambda[down], and
  [ab] = lambdatilde[down].lambdatilde[up].";


aPar::usage = "Parameters of the first row of the twistor matrix in the general parametrization.";
bPar::usage = "Parameters of the third row of the twistor matrix in the general parametrization.";
cPar::usage = "Parameters of the fourth row of the twistor matrix in the general parametrization.";

aPar /: HoldPattern[NumericQ[aPar[a_]]]:=True;
bPar /: HoldPattern[NumericQ[bPar[a_]]]:=True;
cPar /: HoldPattern[NumericQ[cPar[a_]]]:=True;


MomentumTwistorParameters::usage = "List of parameters aPar, bPar and cPar appearing
  in general parametrization of twistor matrix.";
MomentumTwistorMatrixParametrized::usage = "General parametrization of twistor matrix.";
MomentumTwistorParameterSolutionTemplate::usage = "Solution of aPar, bPar, cPar in terms of
  spinors and invariants. Solution given as templates, which have to be replace by either
  S@M or SpinorHelicity6D definitions.";


$momTwistorSpaaProd::usage = "";
$momTwistorSpbbProd::usage = "";

$momTwistorSpaaChain::usage = "";
$momTwistorSpbbChain::usage = "";
$momTwistorSpabChain::usage = "";
$momTwistorSpbaChain::usage = "";

$momTwistorInv::usage = "";


Begin["`Private`"];


(* ::Section::Closed:: *)
(*LambdaTilde from Twistor Variables*)


Options[LambdaTildeFromTwistor]={Modulus->0};

LambdaTildeFromTwistor[twistorMatrix_,opts:OptionsPattern[]] :=
    Module[{n = (Dimensions@twistorMatrix)[[2]], lbdAdjProds, lbdNonAdjProds,
      transpTwist = Transpose@twistorMatrix, padTranspTwistors, lambda, mu},

      padTranspTwistors = Join[{Last@transpTwist}, transpTwist, {First@transpTwist}];
      lambda = padTranspTwistors[[All, 1 ;; 2]];
      lbdAdjProds = Table[Det[{lambda[[j]], lambda[[j + 1]]}], {j, 1, n + 1}];
      lbdNonAdjProds = Table[Det[{lambda[[j - 1]], lambda[[j + 1]]}], {j, 2, n + 1}];
      mu = padTranspTwistors[[All, 3 ;; 4]];

      If[OptionValue[Modulus]===0,
        Table[Plus @@ {
          -lbdAdjProds[[i]] mu[[i - 1]],
          lbdNonAdjProds[[i - 1]] mu[[i]],
          -lbdAdjProds[[i - 1]] mu[[i + 1]]
        }*Power[lbdAdjProds[[i]] lbdAdjProds[[i - 1]], -1],
          {i, 2, n + 1}
        ]
        ,
        Table[Mod[
          Plus @@ Mod[
            {
              -lbdAdjProds[[i]] mu[[i - 1]],
              lbdNonAdjProds[[i - 1]] mu[[i]],
              -lbdAdjProds[[i - 1]] mu[[i + 1]]
            }, OptionValue[Modulus]
          ]*
              PowerMod[
                Mod[
                  lbdAdjProds[[i]] lbdAdjProds[[i - 1]],
                  OptionValue[Modulus]
                ], -1, OptionValue[Modulus]
              ], OptionValue[Modulus]
        ], {i, 2, n + 1}
        ]
      ]
    ];


(* ::Section:: *)
(*Spinors From Twistors*)

Options[GenerateSpinorsFromTwistorMatrix] = {Modulus->0};

LambdaLambdaTildeFromTwistor[twistorMatrix_List,opts:OptionsPattern[]]:=
Module[{lambda,lambdaTilde},
  lambda = (Transpose@twistorMatrix)[[All, 1 ;; 2]];

  lambdaTilde = LambdaTildeFromTwistor[
    twistorMatrix,
    FilterRules[
      {
        opts,
        Options[LambdaLambdaTildeFromTwistor]
      },
      Options[LambdaTildeFromTwistor]
    ]
  ];

  {lambda,lambdaTilde}
];

(* ::Section:: *)
(*Generic Parametrization*)


y[k_]:=Sum[Product[1/aPar[j],{j,1,i}],{i,1,k}]//Simplify;


bTilde[k_,n_]:=bTilde[k-1,n]+aPar[n-k](bTilde[k-1,n]-bTilde[k-2,n])+bPar[k]/.{aPar[i_]/;i>n-2 :>0};
bTilde[0,n_]:=1;
bTilde[k_,n_] /; k<0 := 0;


cTilde[k_,n_]:= cTilde[k-1,n]+aPar[n-k+1](cTilde[k-1,n]-cTilde[k-2,n])+bPar[k-1]/bPar[n-4](cPar[k]-1)/.{aPar[i_]/;i>n-2 :>0,bPar[0]->1};
cTilde[k_,n_]/; k <=0 := 0;


MomentumTwistorParameters[n_]:=Join[
  Table[aPar[i],{i,1,n-2}],
  Table[bPar[i],{i,1,n-4}],
  Table[cPar[i],{i,1,n-4}]
];


MomentumTwistorMatrixParametrized[n_]:={
  Join[{1,0},Table[y[i],{i,1,n-2}]],
  Join[{0},ConstantArray[1,n-1]],
  Join[{0,0,0,bPar[n-4]/aPar[2]},Table[bTilde[i,n],{i,n-5,0,-1}]],
  Join[{0,0,1,1},Table[cTilde[i,n],{i,n-4,1,-1}]]
}/.bPar[0]->aPar[2];

MomentumTwistorParameterSolutionTemplate[n_]:=Join[
  {aPar[1]->$momTwistorInv[1,2]},
  Table[aPar[i]->-$momTwistorSpaaProd[i,i+1]$momTwistorSpaaProd[i+2,1]/
      ($momTwistorSpaaProd[1, i]$momTwistorSpaaProd[i+1,i+2]),{i,2,n-2}],
  {bPar[n-4]->$momTwistorInv[2,3]/$momTwistorInv[1,2]},
  Table[bPar[i]->$momTwistorSpabChain[n-i,n-i+1,2]/$momTwistorSpabChain[n-i,1,2],{i,n-5,1,-1}],
  {cPar[1]->-$momTwistorInv[1,3]/$momTwistorInv[1,2]},
  Table[cPar[i]->-$momTwistorSpabChain[1,3,n-i+2]/$momTwistorSpabChain[1,2,n-i+2],{i,2,n-4}]
];



End[]; (* `Private` *)

EndPackage[]