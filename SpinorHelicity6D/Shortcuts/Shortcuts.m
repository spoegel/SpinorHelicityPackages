(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: Shortcuts *)
(* :Context: Shortcuts` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-11-10 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Shortcuts`"];
(* Exported symbols added here with SymbolName::usage *)

Spaa::usage = "Shortcut for SpinorAnglebracket as well as
  Chain[$angle,_,{__},_,$angle], imitating the notation of S@M";
Spbb::usage = "Shortcut for SpinorSquarebracket as well as
  Chain[$square,_,{__},_,$square], imitating the notation of S@M";
Spab::usage = "Shortcut for Chain[$angle,_,{__},_,$square], imitating the notation of S@M";
Spba::usage = "Shortcut for Chain[$square,_,{__},_,$angle], imitating the notation of S@M";

Begin["`Private`"];

Spaa[a_,b__,c_]:=Chain[$angle,a,{b},c,$angle];
Spab[a_,b__,c_]:=Chain[$angle,a,{b},c,$square];
Spba[a_,b__,c_]:=Chain[$square,a,{b},c,$angle];
Spbb[a_,b__,c_]:=Chain[$square,a,{b},c,$square];

Spaa[a_,b_]:=SpinorAngleBracket[a,b];
Spbb[a_,b_]:=SpinorSquareBracket[a,b];

End[]; (* `Private` *)



EndPackage[]