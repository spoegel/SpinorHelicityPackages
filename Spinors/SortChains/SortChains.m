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

(* :Title: SortChain *)
(* :Context: SortChain` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-06-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Spinors`SortChains`",{"Spinors`"}];
(* Exported symbols added here with SymbolName::usage *)

$SortChains::usage = "If True, spinor chains are automatically sorted according to the usual S@M
rules. Otherwise, chains are left as they are.";
$SortChains = True;

Unprotect[Spaa,Spbb,Spba];


DownValues[Spaa] = DownValues[Spaa]/.{
  Verbatim[HoldPattern[Spaa[(Global`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Global`a_)?SpinorQ]] :> -Spaa[Global`a, Sequence @@ Reverse[{Spinors`Private`p}], Global`a] /;  !OrderedQ[{{Spinors`Private`p}, Reverse[{Spinors`Private`p}]}]]->
      (HoldPattern[Spaa[(Spinors`Private`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`a_)?SpinorQ]] :> -Spaa[Spinors`Private`a,
        Sequence @@ Reverse[{Spinors`Private`p}], Spinors`Private`a] /;  $SortChains&&!OrderedQ[{{Spinors`Private`p}, Reverse[{Spinors`Private`p}]}]),
  Verbatim[HoldPattern[Spaa[(Global`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ,
    (Spinors`Private`b_)?SpinorQ]] :> -Spaa[Spinors`Private`b, Sequence @@
        Reverse[{Spinors`Private`p}], Global`a] /;  !OrderedQ[{Global`a, Spinors`Private`b}]]->
      (HoldPattern[Spaa[(Spinors`Private`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`b_)
          ?SpinorQ]] :> -Spaa[Spinors`Private`b, Sequence @@ Reverse[{Spinors`Private`p}], Spinors`Private`a] /;  $SortChains&&!OrderedQ[{Spinors`Private`a, Spinors`Private`b}])
};


DownValues[Spba] = DownValues[Spba]/.{
  Verbatim[HoldPattern[Spba[(Global`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`b_)?SpinorQ]] :> Spab[Spinors`Private`b, Sequence @@ Reverse[{Spinors`Private`p}], Global`a]]->
      HoldPattern[Spba[(Spinors`Private`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`b_)
          ?SpinorQ]] :> Spab[Spinors`Private`b, Sequence @@ Reverse[{Spinors`Private`p}], Spinors`Private`a] /;$SortChains
};


DownValues[Spbb] = DownValues[Spbb]/.{
  Verbatim[HoldPattern[Spbb[(Global`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Global`a_)?SpinorQ]] :>
      -Spbb[Global`a, Sequence @@ Reverse[{Spinors`Private`p}], Global`a] /;
          !OrderedQ[{Reverse[{Spinors`Private`p}], {Spinors`Private`p}}]]->
      HoldPattern[Spbb[(Spinors`Private`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`a_)?SpinorQ]] :> -Spbb[Spinors`Private`a,
        Sequence @@ Reverse[{Spinors`Private`p}], Spinors`Private`a] /;  $SortChains && !OrderedQ[{Reverse[{Spinors`Private`p}], {Spinors`Private`p}}]
  ,
  Verbatim[HoldPattern[Spbb[(Global`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`b_)?SpinorQ]] :> -Spbb[Spinors`Private`b, Sequence @@ Reverse[{Spinors`Private`p}], Global`a] /; OrderedQ[{Global`a, Spinors`Private`b}] && Global`a =!= Spinors`Private`b]->
      HoldPattern[Spbb[(Spinors`Private`a_)?SpinorQ, (Spinors`Private`p___)?SMatrixQ, (Spinors`Private`b_)
          ?SpinorQ]] :> -Spbb[Spinors`Private`b, Sequence @@ Reverse[{Spinors`Private`p}], Spinors`Private`a] /; $SortChains && OrderedQ[{Spinors`Private`a, Spinors`Private`b}] && Spinors`Private`a =!= Spinors`Private`b
};


Protect[Spaa,Spbb,Spba];



EndPackage[]