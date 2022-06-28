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

(* Mathematica Source File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-06-04 *)

BeginPackage["FiniteField`"];

NFF::usage = "";

$prime = NextPrime[2^31,-1];

Begin["`Private`"];

NFF[expr_]:=Mod[expr/.{
    HoldPattern[Complex[x_,y_]] :> Complex[NFF@x,NFF@y]
  }//.{
    Power[x_,y_] :> PowerMod[x,y,$prime],
   HoldPattern[Rational[x_,y_]] :> x * PowerMod[y,-1,$prime]
  },
  $prime
];

End[];

EndPackage[];