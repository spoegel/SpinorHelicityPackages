#!/usr/bin/env wolframscript

(* Wolfram Script *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-05-12 *)

Needs["Utils`"];
SetOptions[$Output, FormatType -> OutputForm];

logDir = FileNameJoin[{"test_logs",DateString[{"Year", "-", "Month", "-", "Day", "-", "Hour", "-", "Minute", "-", "Second" }] <> "_TwoLoopAllPlusRationals_tests"}];
CreateDirectory[logDir,CreateIntermediateDirectories -> True];
AppendTo[$Output, OpenWrite[FileNameJoin@{logDir,"output.log"},FormatType -> OutputForm]];

$testReports = TestReport/@{
  "verify_4pt.mt"
  ,
  "verify_5pt.mt"
  ,
  "verify_6pt.mt"
};

DumpSave[FileNameJoin[{logDir,"testReports.mx"}],$testReports];

printTestResults/@$testReports