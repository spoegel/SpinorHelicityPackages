#!/usr/bin/env wolframscript

(* Wolfram Script *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-06-10 *)

<<SpinorHelicity6D`
<<SpinorHelicity6D`Amplitudes`

$PMPMMetric = True;

KillMasses[{kBubRef,q,$bubbleReference}];
GenSpinors[{kBubRef,q,$bubbleReference},FourD->{kBubRef,q,$bubbleReference},RandomSpinors->True,Seed->235,ParameterRange->50,DisplaySpinors->False];


RecomputeAmplitudes[];
SaveAmplitudesToFile["",Overwrite->True];