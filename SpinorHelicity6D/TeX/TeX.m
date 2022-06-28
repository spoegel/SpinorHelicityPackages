(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: TeX *)
(* :Context: SpinorHelicity6D`TeX` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-05-18 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`TeX`",
  {"SpinorHelicity6D`","SpinorHelicity6D`SpinorFunctions`","TeXUtilities`"}];
(* Exported symbols added here with SymbolName::usage *)

Begin["`Private`"];

removeUnderBar[x___]:={x}/.UnderBar->Identity//Apply[Sequence];

Unprotect[{SpinorAngleBracket,SpinorSquareBracket,S,Chain,Extramass,
  Extramasstilde,Tr5,TrM,TrP}];

Format[HoldPattern[SpinorAngleBracket[x__]],TeXForm] :=
    TeXDelimited["\\spaa{",x,"}",
      "DelimSeparator"->"","BodySeparator"->" "];
Format[HoldPattern[SpinorSquareBracket[x__]],TeXForm] :=
    TeXDelimited["\\spbb{",x,"}",
      "DelimSeparator"->"","BodySeparator"->" "];
Format[HoldPattern[S[x__]],TeXForm] :=
    TeXDelimited["s_{",x,"}",
      "DelimSeparator"->"","BodySeparator"->" "];
Format[HoldPattern[S6Mom[x__]],TeXForm] :=
    TeXDelimited["s_{",x,"}",
      "DelimSeparator"->"","BodySeparator"->" "];

Format[HoldPattern[Chain[$angle,a_,{b___},c_,$angle]],TeXForm] :=
    TeXDelimited["\\spaa{",a,"|",removeUnderBar@b,"|",c,"}",
      "DelimSeparator"->"", "BodySeparator"->" "];

Format[HoldPattern[Chain[$angle,a_,{b___},c_,$square]],TeXForm] :=
    TeXDelimited["\\spab{",a,"|",removeUnderBar@b,"|",c,"}",
      "DelimSeparator"->"", "BodySeparator"->" "];

Format[HoldPattern[Chain[$square,a_,{b___},c_,$angle]],TeXForm] :=
    TeXDelimited["\\spba{",a,"|",removeUnderBar@b,"|",c,"}",
      "DelimSeparator"->"", "BodySeparator"->" "];

Format[HoldPattern[Chain[$square,a_,{b___},c_,$square]],TeXForm] :=
    TeXDelimited["\\spbb{",a,"|",removeUnderBar@b,"|",c,"}",
      "DelimSeparator"->"", "BodySeparator"->" "];

Format[HoldPattern[Extramass[x___]],TeXForm] :=
    TeXDelimited["\\mu_{",x,"}",
      "DelimSeparator"->"", "BodySeparator"->" "];
Format[HoldPattern[Extramasstilde[x___]],TeXForm] :=
    TeXDelimited["\\tilde{\\mu}_{",x,"}",
        "DelimSeparator"->"", "BodySeparator"->" "];

Format[HoldPattern[TrM[x___]],TeXForm] :=
    TeXDelimited["\\trM (",x,")",
        "DelimSeparator"->"", "BodySeparator"->" "];
Format[HoldPattern[TrP[x___]],TeXForm] :=
    TeXDelimited["\\trP (",x,")",
        "DelimSeparator"->"", "BodySeparator"->" "];
Format[HoldPattern[Tr5[x___]],TeXForm] :=
    TeXDelimited["\\trF (",x,")",
        "DelimSeparator"->"", "BodySeparator"->" "];

(* Setting formatting precedence as described in *)
(* https://mathematica.stackexchange.com/a/153908/73380 *)
(FormatValues@# = Join[Lookup[#,False,{}], Lookup[#,True,{}]]&@
    GroupBy[FormatValues@#, FreeQ@TeXForm]) & /@
        {SpinorAngleBracket,SpinorSquareBracket,S,Chain,Extramass,
          Extramasstilde,Tr5,TrM,TrP};

Protect[{SpinorAngleBracket,SpinorSquareBracket,S,Chain,Extramass,
  Extramasstilde,Tr5,TrM,TrP}];

End[]; (* `Private` *)

EndPackage[]