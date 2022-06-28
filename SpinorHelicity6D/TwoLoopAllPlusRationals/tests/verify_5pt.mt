(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

BeginTestSection["verify_5pt.mt"]

Block[{Print},
  Needs["SpinorHelicity6D`TwoLoopAllPlusRationals`"];
];

$PMPMMetric = True;

Block[{WriteString},
  LoadAmplitudesFile[];
];

$multiplicity = 5;

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
    LoopMomentumPostProcessingFirstPass->
        Simplify@*ReplaceAll[
          {
            Power[x_,exp_]/;exp>=2 :> Expand[Power[x,exp]],
            Power[x_,exp_]/;exp<=-2 :> 1/Expand[Power[x,-exp]]
          }]@*Simplify,
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
  configurations = Configurations[{1,2,3,4,5},{},{}];
  EvaluateConfigurations[configurations]
  ,
  -A52LRatLitBadger[1,2,3,4,5]//ToNum,
  TestID->"R_5:1"
];

VerificationTest[(* 2 *)
  configurations = Configurations[{1,2},{3,4,5},{}];
  EvaluateConfigurations[configurations]
  ,
  -R253[1,2,3,4,5]//ToNum,
  TestID->"R_5:3"
];

VerificationTest[(* 3 *)
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
  -R251B[1,2,3,4,5]//ToNum,
  TestID->"R_5:1B"
];

EndTestSection[]
