#!/usr/bin/env wolframscript

(* Wolfram Script *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-05-11 *)

Block[{Print},
  << SpinorHelicity6D`
  << SpinorHelicity6D`Amplitudes`;
  << SpinorHelicity6D`Numerics`;
  << SpinorHelicity6D`Unitarity`;
  << SpinorHelicity6D`BerendsGiele`;
];

Block[{},
  LoadAmplitudesFile[];
];

Needs["Utils`"];
SetOptions[$Output, FormatType -> OutputForm];

logDir = FileNameJoin[{"test_logs",DateString[{"Year", "-", "Month", "-", "Day", "-", "Hour", "-", "Minute", "-", "Second" }] <> "_Amplitudes_tests"}];
CreateDirectory[logDir,CreateIntermediateDirectories -> True];
AppendTo[$Output, OpenWrite[FileNameJoin@{logDir,"output.log"},FormatType -> OutputForm]];

$testReports = TestReport/@{
  FileNameJoin[{DirectoryName[FindFile["SpinorHelicity6D`Amplitudes`"]],"tests","FourScalarAmplitudesKK.wlt"}],
  FileNameJoin[{DirectoryName[FindFile["SpinorHelicity6D`Amplitudes`"]],"tests","TwoScalarAmplitudesKK.wlt"}],
  FileNameJoin[{DirectoryName[FindFile["SpinorHelicity6D`Amplitudes`"]],"tests","FourScalarGluonContact.wlt"}],
  FileNameJoin[{DirectoryName[FindFile["SpinorHelicity6D`Amplitudes`"]],"tests","FourScalarGluonBerendsGiele.wlt"}],
  FileNameJoin[{DirectoryName[FindFile["SpinorHelicity6D`Amplitudes`"]],"tests","FourScalarBerendsGiele.wlt"}]};

DumpSave[FileNameJoin[{logDir,"testReports.mx"}],$testReports];

printTestResults/@$testReports
