(* Mathematica Package *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

(* :Title: DrawConfigs *)
(* :Context: DrawConfigs` *)
(* :Author: Sebastian Poegel *)
(* :Date: 2021-11-09 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.1 *)
(* :Copyright: (c) 2021 Sebastian Poegel *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinorHelicity6D`TwoLoopAllPlusRationals`DrawConfigs`",
  {
    "SpinorHelicity6D`TwoLoopAllPlusRationals`",
    "SpinorHelicity6D`Unitarity`",
    "Utils`"
  }
];
(* Exported symbols added here with SymbolName::usage *)

(* ::Subsection::Closed:: *)
(*Drawing Total Configurations*)


DrawTotalConfig::usage = "";
DrawTotalConfigData::usage = "";

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Drawing Non-Planar Total Configurations*)


Options[DrawTotalConfig]={InternalLengthBox->0.5,InternalLengthTriangle->0.7,InternalLengthBubble->0.7,
  ExternalLength->0.2,VertexSize->10,BubbleLooseness->0.6,ExternalSpreadAngle->Pi/10,PlotRange->{{-1.4,1.4},{-1,1}}};

DrawTotalConfig[totalConfig_Association,opts:OptionsPattern[]]:=
    Module[{centVertexLabel, verts,vertCoords,edges,vertLabels,edgeStyles,vertexsizes},

      vertexsizes={Rule[centVertexLabel,1]};


      {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles} =
          DrawTotalConfigData[totalConfig, PassRules[DrawTotalConfig,DrawTotalConfigData,opts]];

      GraphPlot[
        Graph[verts,edges,VertexCoordinates->vertCoords,VertexLabels->vertLabels,EdgeStyle->edgeStyles,VertexSize->vertexsizes],
        VertexShapeFunction ->  Function[{pt, v, size}, Disk[pt, If[size === {0.,0.}, Offset[0],Offset[OptionValue[VertexSize]]]]],
        MultiedgeStyle->OptionValue[BubbleLooseness],PlotRange->OptionValue[PlotRange]
      ]
    ];


Options[DrawTotalConfigData]=Options[DrawTotalConfig];

DrawTotalConfigData[totalConfig_Association,opts:OptionsPattern[]]:=
    Module[{centVertexLabel=Unique["vert"], verts,vertCoords={{0,0}},edges={},vertLabels={},
      edgeStyles={},
      vertexsizes},
      verts={centVertexLabel};
      vertexsizes={Rule[centVertexLabel,1]};

      Which[
        TopTwistedConfigQ[totalConfig],
        vertLabels={centVertexLabel->Placed["T",Center]},
        BottomTwistedConfigQ[totalConfig],
        vertLabels={centVertexLabel->Placed["B",Center]},
        _,
        Null
      ];

      {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles} =
          Join[
            GenerateLoopSubGraph[
              totalConfig["FirstLoop"],
              Lookup[totalConfig["Center"],{"Top","FirstInner"}],
              centVertexLabel,
              PassRules[DrawTotalConfig,GenerateLoopSubGraph,opts]
            ],
            {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles},
            2];
      {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles} =
          Join[
            GenerateLoopSubGraphFlipped[
              totalConfig["SecondLoop"],
              Lookup[totalConfig["Center"],{"Bottom","SecondInner"}],
              centVertexLabel,
              PassRules[DrawTotalConfig,GenerateLoopSubGraph,opts]
            ],
            {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles},
            2];
      {verts,vertCoords,edges,vertLabels,vertexsizes,edgeStyles}
    ];

Options[GenerateLoopSubGraph]=Options[DrawTotalConfig];
GenerateLoopSubGraph[loop_List,{centerTop_,centerFirstInner_},centVertexLabel_,opts:OptionsPattern[]]:=Join[
  AttachExt[centerTop,centVertexLabel,{0,0},{0,OptionValue[ExternalLength]},PassRules[GenerateLoopSubGraph,AttachExt,opts]],
  AttachExt[centerFirstInner,centVertexLabel,{0,0},{-OptionValue[ExternalLength],0},PassRules[GenerateLoopSubGraph,AttachExt,opts]],
  Switch[Length@loop,
    3,GenBoxSubGraph[loop,centVertexLabel,{0,0},PassRules[GenerateLoopSubGraph,GenBoxSubGraph,opts]],
    2,GenTriangleSubGraph[loop,centVertexLabel,{0,0},PassRules[GenerateLoopSubGraph,GenTriangleSubGraph,opts]],
    1,GenBubbleSubGraph[loop,centVertexLabel,{0,0},PassRules[GenerateLoopSubGraph,GenBubbleSubGraph,opts]],
    _,Print["Error: Unsupported number of corners in loop."]; {{},{},{},{},{}}
  ],
  2];

Options[GenerateLoopSubGraphFlipped]=Options[DrawTotalConfig];
GenerateLoopSubGraphFlipped[loop_List,{centerBottom_,centerSecondInner_},centVertexLabel_,opts:OptionsPattern[]] :=
    MapAt[
      -1*#&,
      GenerateLoopSubGraph[loop,
        {centerBottom,centerSecondInner},
        centVertexLabel,
        PassRules[GenerateLoopSubGraphFlipped,GenerateLoopSubGraph,opts]
      ],
    {2}
    ]/.{Before->After,After->Before,Above->Below,Below->Above};


Options[GenBoxSubGraph]=Options[DrawTotalConfig];
GenBoxSubGraph[corners_,centVertexLabel_,centCoord_,opts:OptionsPattern[]] /; Length@corners == 3 :=
    Module[{cornersGraph,cornerVerts = Table[Unique["vert"],{i,3}],cornerCoords,cornerEdges,
      cutBubbleEdgeStyles = {},vertexsizes},

      cornerCoords = {centCoord+{-OptionValue[InternalLengthBox],-OptionValue[InternalLengthBox]},
        centCoord+{-2*OptionValue[InternalLengthBox],0},
        centCoord+{-OptionValue[InternalLengthBox],+OptionValue[InternalLengthBox]}};

      cornerEdges = UndirectedEdge@@@Partition[Append[cornerVerts,centVertexLabel],{2},1,1];

      cutBubbleEdgeStyles = CreateCutBubbleEdgeStyle[corners,cornerEdges];

      vertexsizes = Rule[#,1]&/@cornerVerts;

      cornersGraph = {
        cornerVerts,
        cornerCoords,
        cornerEdges,
        {},
        vertexsizes,
        cutBubbleEdgeStyles
      };

      Join[
        cornersGraph,
        AttachExt[corners[[1]]["Outer"],cornerVerts[[1]],cornerCoords[[1]],{0,-OptionValue[ExternalLength]},PassRules[GenBoxSubGraph,AttachExt,opts]],
        AttachExt[corners[[1]]["Inner"],cornerVerts[[1]],cornerCoords[[1]],{0,+OptionValue[ExternalLength]},PassRules[GenBoxSubGraph,AttachExt,opts]],
        AttachExt[corners[[2]]["Outer"],cornerVerts[[2]],cornerCoords[[2]],{-OptionValue[ExternalLength],0},PassRules[GenBoxSubGraph,AttachExt,opts]],
        AttachExt[corners[[2]]["Inner"],cornerVerts[[2]],cornerCoords[[2]],{+OptionValue[ExternalLength],0},PassRules[GenBoxSubGraph,AttachExt,opts]],
        AttachExt[corners[[3]]["Outer"],cornerVerts[[3]],cornerCoords[[3]],{0,+OptionValue[ExternalLength]},PassRules[GenBoxSubGraph,AttachExt,opts]],
        AttachExt[corners[[3]]["Inner"],cornerVerts[[3]],cornerCoords[[3]],{0,-OptionValue[ExternalLength]},PassRules[GenBoxSubGraph,AttachExt,opts]],
        2]
    ];


Options[GenTriangleSubGraph]=Options[DrawTotalConfig];
GenTriangleSubGraph[corners_,centVertexLabel_,centCoord_,opts:OptionsPattern[]] /; Length@corners == 2 :=
    Module[{cornersGraph,cornerVerts = Table[Unique["vert"],{i,2}],cornerCoords,cornerEdges,cutBubbleEdgeStyles={},vertexsizes},
      cornerCoords = {centCoord+{-OptionValue[InternalLengthTriangle],-OptionValue[InternalLengthTriangle]/Sqrt[3]},
        centCoord+{-OptionValue[InternalLengthTriangle],+OptionValue[InternalLengthTriangle]/Sqrt[3]}};

      cornerEdges = UndirectedEdge@@@Partition[Append[cornerVerts,centVertexLabel],{2},1,1];

      cutBubbleEdgeStyles = CreateCutBubbleEdgeStyle[corners,cornerEdges];

      vertexsizes = Rule[#,1]&/@cornerVerts;

      cornersGraph = {
        cornerVerts,
        cornerCoords,
        cornerEdges,
        {},
        vertexsizes,
        cutBubbleEdgeStyles
      };
      Join[
        cornersGraph,
        AttachExt[corners[[1]]["Outer"],cornerVerts[[1]],cornerCoords[[1]],OptionValue[ExternalLength]*Normalize[{-Tan[Pi/6],-1}],PassRules[GenTriangleSubGraph,AttachExt,opts]],
        AttachExt[corners[[1]]["Inner"],cornerVerts[[1]],cornerCoords[[1]],OptionValue[ExternalLength]*Normalize[{Tan[Pi/6],1}],PassRules[GenTriangleSubGraph,AttachExt,opts]],
        AttachExt[corners[[2]]["Outer"],cornerVerts[[2]],cornerCoords[[2]],OptionValue[ExternalLength]*Normalize[{-Tan[Pi/6],+1}],PassRules[GenTriangleSubGraph,AttachExt,opts]],
        AttachExt[corners[[2]]["Inner"],cornerVerts[[2]],cornerCoords[[2]],OptionValue[ExternalLength]*Normalize[{Tan[Pi/6],-1}],PassRules[GenTriangleSubGraph,AttachExt,opts]],
        2]
    ];


Options[GenBubbleSubGraph]=Options[DrawTotalConfig];
GenBubbleSubGraph[corners_,centVertexLabel_,centCoord_,opts:OptionsPattern[]] /; Length@corners == 1 :=
    Module[{cornersGraph,cornerVerts = Table[Unique["vert"],{i,1}],cornerCoords,cornerEdges,
      cutBubbleEdgeStyles={},vertexsizes},
      cornerCoords = {centCoord+{-OptionValue[InternalLengthBubble],0}};

      cornerEdges = UndirectedEdge@@@Partition[Append[cornerVerts,centVertexLabel],{2},1,1];

      cutBubbleEdgeStyles = CreateCutBubbleEdgeStyle[corners,cornerEdges];

      vertexsizes = Rule[#,1]&/@cornerVerts;

      cornersGraph = {
        cornerVerts,
        cornerCoords,
        cornerEdges,
        {},
        vertexsizes,
        cutBubbleEdgeStyles
      };

      Join[
        cornersGraph,
        AttachExt[corners[[1]]["Outer"],cornerVerts[[1]],cornerCoords[[1]],{-OptionValue[ExternalLength],0},PassRules[GenBubbleSubGraph,AttachExt,opts]],
        AttachExt[corners[[1]]["Inner"],cornerVerts[[1]],cornerCoords[[1]],{+OptionValue[ExternalLength],0},PassRules[GenBubbleSubGraph,AttachExt,opts]],
        2]
    ];





CreateCutBubbleEdgeStyle[corners_,cornerEdges_] := Switch[Count[corners,$cutBubble,Infinity],
  0, {},
  1, Switch[Position[corners,$cutBubble][[1,1]],
    1,{Last@cornerEdges -> Dashed},
    Length@corners,{cornerEdges[[-2]] -> Dashed},
    _, {}
  ],
  2, {First@Take[cornerEdges,{Position[corners,$cutBubble][[1,1]]}] -> Dashed}
];


Options[AttachExt]=Options[DrawTotalConfig];
AttachExt[{},___] := {{},{},{},{},{}};

(*AttachExt[extMomenta_,corner_,cornerCoord_,offset_,opts:OptionsPattern[]] /; Length@extMomenta >1 :=
Module[{maxAng = OptionValue[ExternalSpreadAngle],uniqueMomenta = Unique[Symbol["p"<>ToString[#]]]&/@extMomenta},
	{
		uniqueMomenta,
		Table[cornerCoord+offset.RotationMatrix[ang],{ang,-maxAng,maxAng,2*maxAng/(-1+Length@extMomenta)}],
		Map[UndirectedEdge[corner,#]&,uniqueMomenta],
		MapThread[Rule,{uniqueMomenta,Placed[ToString[#],LabelPlacement[offset]]&/@extMomenta}],
		Rule[#,0]&/@uniqueMomenta,
		{}
	}
];
*)

AttachExt[extMomenta_,corner_,cornerCoord_,offset_,opts:OptionsPattern[]] /; Length@extMomenta >1 :=
    Module[{angle = OptionValue[ExternalSpreadAngle],uniqueMomenta = Unique[Symbol["p"<>ToString[#]]]&/@extMomenta,n=Length@extMomenta},
      {
        uniqueMomenta,
        Table[cornerCoord+offset . RotationMatrix[ang],{ang,-(n-1)*angle/2,(n-1)*angle/2,angle}],
        Map[UndirectedEdge[corner,#]&,uniqueMomenta],
        MapThread[Rule,{uniqueMomenta,Placed[ToString[#],LabelPlacement[offset]]&/@extMomenta}],
        Rule[#,0]&/@uniqueMomenta,
        {}
      }
    ];

AttachExt[extMomenta_,corner_,cornerCoord_,offset_,opts:OptionsPattern[]] /; Length@extMomenta == 1 :=
    Module[{uniqueMomenta = Unique[Symbol["p"<>ToString[#]]]&/@extMomenta},
      {
        uniqueMomenta,
        {cornerCoord+offset},
        Map[UndirectedEdge[corner,#]&,uniqueMomenta],
        MapThread[Rule,{uniqueMomenta,Placed[ToString[#],LabelPlacement[offset]]&/@extMomenta}],
        Rule[#,0]&/@uniqueMomenta,
        {}
      }
    ];


LabelPlacement[offset_]:=
    Switch[Sign[offset/.(Rule[#,1]&/@Variables[offset])//N],
      {-1,-1},Below,
      {-1,0},Before,
      {-1,1},Above,
      {1,1},Above,
      {0,1},Above,
      {0,-1},Below,
      {1,1},Above,
      {1,0},After,
      {1,-1},Below,
      _,Center
    ];



End[]; (* `Private` *)

EndPackage[]