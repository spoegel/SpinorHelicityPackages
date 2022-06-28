(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: ConfigurationTools *)
(* :Context: ConfigurationTools` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-04-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`Unitarity`ConfigurationTools`",
              "SpinorHelicity6D`","SpinorHelicity6D`Unitarity`",
              "SpinorHelicity6D`Amplitudes`","SpinorHelicity6D`BerendsGiele`","Utils`"];
(* Exported symbols added here with SymbolName::usage *)

ConstructLoopAmplitudes::usage = "";
ConstructCornerAmplitude::usage = "";

ClosedLoop::usage = "Option for ConstructLoopAmplitudes to signal that the constructed trees should close on themselves
                     in a loop. A counter example where this is not wanted are one-loop squared configurations with
                     a central vertex.";

GenConfigBare::usage = "";
GenConfigsBare::usage = "";
GenConfigPartitions::usage = "";
GenSubConfigs::usage ="";
GetSplittingMomenta::usage = "";



Begin["`Private`"];

(* ::Section:: *)
(*Tree Amplitudes*)


Options[ConstructLoopAmplitudes] = {Gravity->False, ClosedLoop -> False};

ConstructLoopAmplitudes[corners_List,loopMom_List,opts:OptionsPattern[]]:=Module[{loopMomentaPairs},
  loopMomentaPairs = If[OptionValue[ClosedLoop] === True, Partition[loopMom,{2},1,1], Partition[loopMom,{2},1]];
  Product[
    ConstructCornerAmplitude[corner[[1]],Sequence@@corner[[2]],PassRules[ConstructLoopAmplitudes,ConstructCornerAmplitude,opts]],
    {corner,Transpose[{corners,loopMomentaPairs}]}]
];


Options[ConstructCornerAmplitude] = {Gravity->False};
ConstructCornerAmplitude[momenta_Association,lMomIn_,lMomOut_,opts:OptionsPattern[]] /; !VectorQ[momenta["Inner"],ListQ] || !VectorQ[momenta["Outer"],ListQ] :=
    Module[{TreeHead},
      TreeHead = Which[
        OptionValue[Gravity]===True,
          Inactive@MTree,
        ContainsAll[momenta["Properties"],{"OneLoop"}],
          Inactive@AOneLoop,
        ContainsAll[momenta["Properties"],{"OneLoop","Gluon"}],
          Inactive@AOneLoopGluon,
        ContainsAll[momenta["Properties"],{"OneLoop","Contact"}],
          Inactive@AOneLoopGluon,
        ContainsAll[momenta["Properties"],{"Gluon"}],
          Inactive@ATreeGluon,
        ContainsAll[momenta["Properties"],{"Contact"}],
          Inactive@ATreeContact,
        ContainsAll[momenta["Properties"],{"BerendsGiele"}],
          Inactive@ATreeBerendsGiele,
        True,
          Inactive@ATree
      ];

      TreeHead[
        "S",
        Sequence@@ConstantArray["+",Length@momenta["Outer"]],
        "S",
        Sequence@@ConstantArray["+",Length@momenta["Inner"]]
      ][
        pm[lMomIn],
        Sequence@@momenta["Outer"],
        lMomOut,
        Sequence@@momenta["Inner"]
      ]
    ];

(* For the one-loop case, where external particles are specified *)
ConstructCornerAmplitude[momenta_Association,lMomIn_,lMomOut_,opts:OptionsPattern[]] /; VectorQ[momenta["Inner"],ListQ] && VectorQ[momenta["Outer"],ListQ] :=
    Module[{TreeHead,labelsOuter,typesOuter,labelsInner,typesInner},

      {typesInner,labelsInner} = Transpose[momenta["Inner"]]/.{{}->{{},{}}};
      {typesOuter,labelsOuter} = Transpose[momenta["Outer"]]/.{{}->{{},{}}};

      TreeHead = Which[
        OptionValue[Gravity]===True,
         Inactive@MTree,
        ContainsAll[momenta["Properties"],{"OneLoop"}],
          Inactive@AOneLoop,
        ContainsAll[momenta["Properties"],{"OneLoop","Gluon"}],
          Inactive@AOneLoopGluon,
        ContainsAll[momenta["Properties"],{"OneLoop","Contact"}],
         Inactive@AOneLoopGluon,
        ContainsAll[momenta["Properties"],{"Gluon"}],
          Inactive@ATreeGluon,
        ContainsAll[momenta["Properties"],{"Contact"}],
         Inactive@ATreeContact,
        True,
          Inactive@ATree
      ];

      TreeHead["S",Sequence@@typesOuter,"S",Sequence@@typesInner][
        pm[lMomIn],Sequence@@labelsOuter,lMomOut,Sequence@@labelsInner]
    ];

(* ::Section:: *)
(*Configuration Combinatorics*)


GenConfigBare[nBlocks_,momenta_List,blocks_,start_]/;DuplicateFreeQ[momenta]:=Module[{valid={{start}},newChains,n=Length@momenta},
  For[i=2,i<=nBlocks,i++,
    newChains=Function[chain,
      Append[chain,#]&/@Select[blocks,First[#]=== Part[momenta, Mod[Position[momenta,Last[Last[chain]]][[1,1]]+1,n,1]] && (Length[#]+Length[Flatten@chain]<=n)&]
    ]/@valid;
    valid=Flatten[newChains,1];
  ];
  valid
];


GenConfigsBare[nBlocks_,n_]:=Module[{configs,res},
  If[nBlocks===2,
    configs=Join@@(Partition[Range@n,#,1,1]&/@Range[2,n-nBlocks+1]);,
    configs=Join@@(Partition[Range@n,#,1,1]&/@Range[1,n-nBlocks+1]);
  ];
  res=Flatten[Map[GenConfigBare[nBlocks,Range@n,configs,#]&,configs],1];
  res=PrimaryPermutations[res];
  res=Select[res,Length[Flatten@#]===n&];
  res
];

GenConfigsBare[nBlocks_,lbls_List]:=Module[{configs,res},
  If[nBlocks===2,
    configs=Join@@(Partition[lbls,#,1,1]&/@Range[2,Length@lbls-nBlocks+1]);,
    configs=Join@@(Partition[lbls,#,1,1]&/@Range[1,Length@lbls-nBlocks+1]);
  ];
  res=Flatten[Map[GenConfigBare[nBlocks,lbls,configs,#]&,configs],1];
  res=PrimaryPermutations[res];
  res=Select[res,Length[Flatten@#]===Length@lbls &];
  res
];


(* Given a list of elements, this function returns all possible ways to cut the list into nblocks pieces *)
(* E.g.: GenConfigPartitions[3,{a,b,c,d,e}] \[Rule] {{{a},{b},{c,d,e}},{{a},{b,c},{d,e}},{{a},{b,c,d},{e}},{{a,b},{c},{d,e}},{{a,b},{c,d},{e}},{{a,b,c},{d},{e}}} *)
GenConfigPartitions[nBlocks_,lbls_List]:=
    Module[{configs=Join@@(Partition[lbls,#,1,1]&/@Range[1,Length@lbls-nBlocks+1]),res,startConfigs},
  startConfigs=Cases[configs,{First@lbls,___}];

  res=Flatten[Map[GenConfigBare[nBlocks,lbls,configs,#]&,startConfigs],1];
  res=PrimaryPermutations[res];
  res=Select[res,Length[Flatten@#]===Length@lbls &];
  res
];


(* Given a configuration with nBlocks pieces, this functions return all possible ways to introduce nSplittings cuts to obtain nBlocks+nSplittings pieces *)
(* E.g.: GenSubConfigs[{{a,b},{c,d,e},{f,g}},1] \[Rule] {{{a,b},{c,d,e},{f},{g}},{{a,b},{c},{d,e},{f,g}},{{a,b},{c,d},{e},{f,g}},{{a},{b},{c,d,e},{f,g}}} *)
(* Dear god, this thing is a nightmare.... *)
GenSubConfigs[initConfig_,nSplittings_]:=Module[{configConfigs,newSplittings,validSplittings,tempSplitting,block},
  configConfigs=Function[confVar,Join[Sequence@@GenConfigPartitions[1+#,confVar]&/@Range[0,Min[nSplittings,Length@confVar-1]]]]/@initConfig;

  validSplittings=configConfigs[[1]];
  For[block=2,block<=Length@initConfig,block++,
    newSplittings=Function[config,
      Join[config,#]&/@Select[configConfigs[[block]],Length[#]-1 <= nSplittings - (Length[config]-block+1)&]
    ]/@validSplittings;
    validSplittings=Flatten[newSplittings,1];
  ];
  Select[validSplittings,Length[#] == Length[initConfig]+nSplittings&]
];


(* Returns the momenta needed to compute the box subtraction denominators*)
(* E.g.: GetSplittingMomenta[{{1,2},{3,4}},{{1,2},{3},{4}}] \[Rule] {1,2,3} *)
GetSplittingMomenta[tri_,box_]:=Module[{boxSeq={},triSeq={},i=1},
  While[i<=Min[Length@tri,Length@box] && Length@boxSeq == Length@triSeq,
    boxSeq=Join[boxSeq,box[[i]]];
    triSeq=Join[triSeq,tri[[i]]];
    i++;
  ];
  First@MinimalBy[{boxSeq,triSeq},Length]
];

(* ::Chapter:: *)
(*Postamble*)
End[]; (* `Private` *)

EndPackage[]

