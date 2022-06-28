(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

BeginTestSection["verify_6pt.mt"]

Block[{Print},
  Needs["SpinorHelicity6D`TwoLoopAllPlusRationals`"];
];

$PMPMMetric = True;

Block[{WriteString},
  LoadAmplitudesFile[];
];

$multiplicity = 6;

NewProcess[];
KillMasses[Range@$multiplicity];
GenSpinors[Range@$multiplicity, FourD->Range@$multiplicity, Seed->1328135534,
  ParameterRange->50, Parametric->False, DisplaySpinors->False];


KillMasses[{q,$bubbleReference}];
GenSpinors[{q,$bubbleReference}, FourD->{q,$bubbleReference},
  RandomSpinors->True, Seed->24092444, ParameterRange->50];

If[$MachineName === "qft1",
  nKernels = 64,
  nKernels = 4
];

If[$KernelCount === 0,
  LaunchKernels[nKernels];
];

DistributeDefinitions["Utils`*"];
DistributeDefinitions["RationalSeries`*"];
DistributeDefinitions["NumSymb`*"];
DistributeDefinitions["SpinorHelicity6D`*"];

EvaluateConfigurations[configurations_List] :=
  4*ComputeConfigurations[configurations,
    BubbleReference->$bubbleReference, Recompute->True,
    ComputeLoopMomentumSpinors->False,PostProcessing->True,
    LogDirectory->FileNameJoin[{Directory[],"logs"}],DataDirectory -> FileNameJoin@{Directory[],"data"},
    LoopMomentumPostProcessingFirstPass->RootReduce,
    ShuffleConfigs->False];

Configurations[tr1_List, tr2_List, tr3_List]:=
  Join[
    Join@@(GenerateComputationTotalConfigs@@@
        Tuples[{CyclicPermutations[tr1],CyclicPermutations[tr2],CyclicPermutations[tr3]}]),
    Join@@(GenerateComputationTotalConfigs@@@
        Tuples[{CyclicPermutations[tr2],CyclicPermutations[tr3],CyclicPermutations[tr1]}]),
    Join@@(GenerateComputationTotalConfigs@@@
        Tuples[{CyclicPermutations[tr3],CyclicPermutations[tr1],CyclicPermutations[tr2]}])
  ];


VerificationTest[(* 1 *)
  configurations = Configurations[{1,2,3,4,5,6},{},{}];
  EvaluateConfigurations[configurations]
  ,
  R261[1,2,3,4,5,6]//ToNum,
  TestID->"R_6:1"
];

VerificationTest[(* 2 *)
  configurations = Configurations[{1,2},{3,4,5,6},{}];
  EvaluateConfigurations[configurations]
  ,
  R263[1,2,3,4,5,6]//ToNum,
  TestID->"R_6:3"
];

VerificationTest[(* 3 *)
  configurations = Configurations[{1,2,3},{4,5,6},{}];
  EvaluateConfigurations[configurations]
  ,
  R264[1,2,3,4,5,6]//ToNum,
  TestID->"R_6:4"
];

VerificationTest[(* 4 *)
  configurations = Configurations[{1,2},{3,4},{5,6}];
  EvaluateConfigurations[configurations]
  ,
  R2622[1,2,3,4,5,6]//ToNum,
  TestID->"R_6:22"
];

VerificationTest[(* 5 *)
  configurationsSingle =
      GenerateComputeSubLeadingSingleTraceConfigurations@Range@$multiplicity;

  configurations =
      Join@@Table[
        MapAt[
          ReplaceAll[MapThread[Rule,{Range@$multiplicity,perm}]],
          configurationsSingle,
          {All,2}
        ],
      {perm,CyclicPermutations[Range@$multiplicity]}];

  EvaluateConfigurations[configurations]
  ,
  R261B[1,2,3,4,5,6]//ToNum,
  TestID->"R_6:1B"
];

EndTestSection[]
