(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: ConfigurationGeneration *)
(* :Context: ConfigurationGeneration` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2022-05-16 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2022 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`TwoLoopAllPlusRationals`Configurations`",
  {
    "Utils`",
    "SpinorHelicity6D`TwoLoopAllPlusRationals`Helpers`"
  }
];
(* Exported symbols added here with SymbolName::usage *)

GeneratingCuts::usage = "GeneratingCuts[Tr1, Tr2, Tr3]: Generates the set AllLabels^{Tr1,Tr2,Tr3} for Tr1, Tr2, Tr3 being distinct.";

CutLabelToTex::usage = "CutLabelToTex[label]: Returns TeX form of a cut label.";

CutLabelToConfig::usage = "CutLabelToConfig[label]: Returns configuration in the old format that is associated to label.";

ConfigToCutLabel::usage = "ConfigToCutLabel[config]: Returns cut label associated to configuration config in the old format.";

Config1BTo1BCutLabel::usage = "Config1BTo1BCutLabel[config]: Returns 1B cut label of 1B config.";

CutLabel1BToConfig1B::usage = "CutLabel1BToConfig1B[label]: Returns 1B config of 1B cut label.";

Generate1BLabels::usage = "Generate1BLabels[n]: Generates the 1B generating set for momenta 1,...,n.";

Begin["`Private`"];

(*Creating labels for different cut classes*)
GeneratingCuts[tO_List,tL_List,tR_List]:=
    Join@@Table[GeneratingCutsClass[class,tO,tL,tR],{class,Tuples[Range@3,{2}]}];

GeneratingCutsClass[class_List,tO_List,tL_List,tR_List]:=
    Select[
      BuildCutLabel@@@Tuples[{
        DistributeList[tO,class[[1]]+class[[2]]+2],
        DistributeList[Reverse@tL,class[[1]]+1],
        DistributeList[Reverse@tR,class[[2]]+1]
      }],
      !EmptyAmplitudeLabelQ[#] && !MasslessBubbleLabelQ[#] &
    ];
BuildCutLabel[outerDist_List,leftDist_List,rightDist_List]:=Module[
  {class={Length@leftDist,Length@rightDist}-1},
  <|
    "Outer"->SplitAt[outerDist,{class[[1]],class[[1]]+1,class[[1]]+class[[2]]+1}],
    "Left"->SplitAt[RotateRight[leftDist],1],
    "Right"->SplitAt[RotateRight[rightDist],1]
  |>
];

(*Functions to remove massless bubbles*)
LeftBubbleLabelQ[label_Association]:=Length[label["Outer"][[1]]]===1;
RightBubbleLabelQ[label_Association]:=Length[label["Outer"][[3]]]===1;
LeftMasslessBubbleLabelQ[label_Association]:=
    LeftBubbleLabelQ[label] && Length@Flatten[{label["Outer"][[1]],label["Left"][[2]]}] < 2;
RightMasslessBubbleLabelQ[label_Association]:=
    RightBubbleLabelQ[label] && Length@Flatten[{label["Outer"][[3]],label["Right"][[2]]}] < 2;
MasslessBubbleLabelQ[label_Association]:=
    LeftMasslessBubbleLabelQ[label] || RightMasslessBubbleLabelQ[label];

(*Functions to remove labels with empty amplitudes*)
LeftEmptyAmplitudeLabelQ[label_Association]:=
    Or @@ (Length@Flatten[#] < 1&/@
        Transpose[{label["Outer"][[1]],label["Left"][[2]]}]);

RightEmptyAmplitudeLabelQ[label_Association]:=
    Or @@ (Length@Flatten[#] < 1&/@
        Transpose[{label["Outer"][[3]],label["Right"][[2]]}]);
EmptyAmplitudeLabelQ[label_Association]:=
    LeftEmptyAmplitudeLabelQ[label] || RightEmptyAmplitudeLabelQ[label];


(*Function for turning Association into Latex form*)
CutLabelToTex[label_Association]:=
    "(" <>
        CutLabelOuterTeX[label["Outer"]] <>
        CutLabelInnerTeX[label["Left"],"Left"] <>
        CutLabelInnerTeX[label["Right"],"Right"] <>
    ")";


CutLabelOuterTeX[label_]:=
    StringRiffle[
      StringTake[StringReplace[#,{"{"->"(","}"->")"}],{2,-2}]&/@ToString/@label,
      ";"
    ];

CutLabelInnerTeX[label_,side:"Left"|"Right"]:=
    StringRiffle[
      StringTake[StringReplace[#,{"{"->"(","}"->")"}],{2,-2}]&/@ToString/@label,
      {"|_"<>StringTake[side,1],";",""}
    ];

(*Function to turn (subleading triple trace cut labels) into the old cut format*)
CutLabelToConfig[label_Association]:=<|
  "FirstLoop"->LRLoopConfig[label["Outer"][[1]],label["Left"][[2]]],
  "SecondLoop"->LRLoopConfig[label["Outer"][[3]],label["Right"][[2]]],
  "Center"->CenterConfig[label["Left"][[1]],label["Outer"][[2]],label["Right"][[1]],label["Outer"][[4]]]
|>;

CutLabel1BToConfig1B[label_List]:=<|
  "FirstLoop"->LRLoopConfig[label[[1,"Loop"]],label[[3,"Loop"]]],
  "SecondLoop"->LRLoopConfig[label[[4,"Loop"]],label[[2,"Loop"]]],
  "Center"->CenterConfig[{label[[3,"Central"]]},{label[[2,"Central"]]},{label[[1,"Central"]]},{label[[4,"Central"]]},"Twisted"->"Top"]
|>;


LRLoopConfig[outer_List,inner_List]:=
    <|
      "Outer"->#[[1]],
      "Inner"->Reverse@#[[2]],
      "Properties"->{}
    |>&/@Transpose@{outer,inner};

Options[CenterConfig]={"Twisted"->Null};
CenterConfig[left_List,top_List,right_List,bottom_List,opts:OptionsPattern[]]:=<|
  "Top"->First@top,
  "Bottom"->First@bottom,
  "FirstInner"->First@left,
  "SecondInner"->First@right,
  "Twisted"->OptionValue["Twisted"],"Properties"->{}
|>;

ConfigToCutLabel[config_Association]:=<|
  "Outer"->{config[["FirstLoop",All,"Outer"]],{config[["Center","Top"]]},config[["SecondLoop",All,"Outer"]],{config[["Center","Bottom"]]}},
  "Left"->{{config[["Center","FirstInner"]]},Reverse/@config[["FirstLoop",All,"Inner"]]},
  "Right"->{{config[["Center","SecondInner"]]},Reverse/@config[["SecondLoop",All,"Inner"]]}
|>

(* ::Section::Closed:: *)
(*1B Configs to 1B Cut-Labels*)

Config1BTo1BCutLabel[config_Association] /; config["Center","Twisted"] === "Top":=Module[
  {
    E1,E2,E3,E4,C1,C2,C3,C4
  },
  E1=config[["FirstLoop",All,"Outer"]];
  E2=Reverse@config[["SecondLoop",All,"Inner"]];
  E3=Reverse@config[["FirstLoop",All,"Inner"]];
  E4=config[["SecondLoop",All,"Outer"]];

  C1=config[["Center","SecondInner"]];
  C2=config[["Center","Top"]];
  C3=config[["Center","FirstInner"]];
  C4=config[["Center","Bottom"]];

  {
    <|
      "Loop"->E1,
      "Central"->C1
    |>,
    <|
      "Loop"->E2,
      "Central"->C2
    |>,
    <|
      "Loop"->E3,
      "Central"->C3
    |>,
    <|
      "Loop"->E4,
      "Central"->C4
    |>
  }
];


Config1BTo1BCutLabel[config_Association] /; config["Center","Twisted"] === "Bottom":=Module[
  {
    E1,E2,E3,E4,C1,C2,C3,C4
  },
  E1=config[["FirstLoop",All,"Outer"]];
  E2=config[["SecondLoop",All,"Outer"]];
  E3=Reverse@config[["FirstLoop",All,"Inner"]];
  E4=Reverse@config[["SecondLoop",All,"Inner"]];

  C1=config[["Center","Top"]];
  C2=config[["Center","FirstInner"]];
  C3=config[["Center","Bottom"]];
  C4=config[["Center","SecondInner"]];

  {
    <|
      "Loop"->E1,
      "Central"->C1
    |>,
    <|
      "Loop"->E2,
      "Central"->C2
    |>,
    <|
      "Loop"->E3,
      "Central"->C3
    |>,
    <|
      "Loop"->E4,
      "Central"->C4
    |>
  }
]

(* ::Section::Closed:: *)
(* 1B Cut Label Generation *)

Options[Generate1BLabels] = {TrianglesOnly -> False};

Generate1BLabels[n_Integer,opts:OptionsPattern[]]:=
    Delete1BDuplicateReps[{},
      Config1BTo1BCutLabel/@Flatten[
        DistributeSubleadingSingleTrace[#,Range@n] & /@
            DistributeTemplateFreeToFreeInnerFreeOuter /@
                CreateTwistedCentralVertex /@
                    FilterDoubleTriangleTemplates[OptionValue[TrianglesOnly]] @
                        CreateTotalConfigTemplates[n]
      ]
    ];

Generate1BRepresentations[label1B_List]:=
    Table[
      RotateLeft[label1B,rot]/.MapThread[
        Rule,
        {
          Momenta1BLabel[label1B],
          RotateRight[Momenta1BLabel[label1B],Length@Momenta1BLabel[label1B[[;;rot]]]]
        }
      ],
      {rot,0,3}
    ];

Momenta1BLabel[label_List]:=Flatten[Values/@label];

Delete1BDuplicateReps[passed_,{}]:=Rotate1BToCanonical/@passed;
Delete1BDuplicateReps[passed_List,rest_List]:=
    Delete1BDuplicateReps[
      Append[passed,First@rest],
      DeleteCases[Rest@rest,Alternatives@@Generate1BRepresentations[First@rest]]
    ];

Rotate1BToCanonical[label_]:=RotateLeft[label,-1+First@First@Position[Values/@label,1]];

End[]; (* `Private` *)

EndPackage[]