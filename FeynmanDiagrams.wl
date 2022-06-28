(* ::Package:: *)

(* ::Title:: *)
(*Feynman Diagram Generation*)


(* ::Section:: *)
(*Helpers*)


distributeLoops[nLoops_Integer,nVerts_Integer] := Flatten[Permutations/@(PadLeft[#,nVerts]&/@IntegerPartitions[nLoops,nVerts]),{1,2}]


mergeCurrents[currents__]:=(Join@@@Transpose[#])&/@Tuples[{currents}];


ScalesLoopFreeQ[graph_]:=Module[{adjacencyMatrix,external,externalConnected,externalConnectedBubble},
	adjacencyMatrix = Normal@AdjacencyMatrix[graph];
	external = Position[adjacencyMatrix,{0...,1,0...},1];
	externalConnected = FirstPosition[#,1]&/@Extract[adjacencyMatrix,external];
	externalConnectedBubble = Cases[Extract[adjacencyMatrix,externalConnected],{___,2,___}];
	AllTrue[externalConnectedBubble,Count[#,1]>=2&]
]


(* ::Section:: *)
(*Currents*)


(* ::Subsection:: *)
(*Gluon Current*)


CurrentJg[{},_,_] := {};

CurrentJg[{{_,_,___}},_,_] := {};
CurrentJg[{{__}},nLoops_,_]/; nLoops > 0 := {};
CurrentJg[{{"g",label_,___}},0,prevVertexIndex_] := Module[{endVertexIndex=Unique[]},
	{
		{{endVertexIndex},{UndirectedEdge[endVertexIndex,prevVertexIndex]},{UndirectedEdge[endVertexIndex,prevVertexIndex]->Thick},{endVertexIndex->ToString@label}}
	}
];

CurrentJg[particles_List,nLoops_Integer,prevVertexIndex_]/; Length@particles >1 := 
Join[
	currentJggg[particles,nLoops,prevVertexIndex],
	currentJgggg[particles,nLoops,prevVertexIndex]
	currentJggss[particles,nLoops,prevVertexIndex],
	currentJgssg[particles,nLoops,prevVertexIndex],
	currentJggSS[particles,nLoops,prevVertexIndex],
	currentJgSSg[particles,nLoops,prevVertexIndex],
	currentJgSS[particles,nLoops,prevVertexIndex],
	currentJgss[particles,nLoops,prevVertexIndex]
	(*currentJgggLoop[particles,nLoops,prevVertexIndex],
	currentJggggLoop[particles,nLoops,prevVertexIndex]*)(*,
	currentJgssLoop[particles,nLoops,prevVertexIndex],
	currentJgSSLoop[particles,nLoops,prevVertexIndex],
	currentJggssLoop[particles,nLoops,prevVertexIndex],
	currentJgssgLoop[particles,nLoops,prevVertexIndex],
	currentJggSSLoop[particles,nLoops,prevVertexIndex],
	currentJgSSgLoop[particles,nLoops,prevVertexIndex]*)
]


(* ::Subsubsection:: *)
(*currentJggg*)


currentJggg[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJg[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJg[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgggg*)


ClearAll[currentJgggg]
currentJgggg[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJgggg[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJg[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJg[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgssg*)


ClearAll[currentJgssg]
currentJgssg[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJgssg[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJs[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJs[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJg[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJggss*)


ClearAll[currentJggss]
currentJggss[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJggss[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJs[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJs[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgSSg*)


ClearAll[currentJgSSg]
currentJgSSg[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJgSSg[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJS[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJS[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJg[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJggSS*)


ClearAll[currentJggSS]
currentJggSS[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJggSS[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJS[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJS[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgSS*)


currentJgSS[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJS[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJS[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgss*)


currentJgss[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJs[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJs[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection:: *)
(*currentJgggLoop*)


currentJgggLoop[_,nLoops_Integer,_]/;nLoops <= 0 || nLoops == 2 := {};

currentJgggLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,vertexIndex=Unique[],graphExtension,protoDiagrams},

	protoDiagrams=CurrentJg[Prepend[particles,{"g",Unique[internalParticle],vertexIndex}],nLoops-1,vertexIndex];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection:: *)
(*currentJggggLoop*)


currentJggggLoop[_,nLoops_,_]/;nLoops <= 0 || nLoops == 2:= {};

currentJggggLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i]],loopAssignment[[1]],vertexIndex],
				CurrentJg[Prepend[particles[[i+1;;]],{"g",Unique[internalParticle],vertexIndex}],loopAssignment[[2]]-1,vertexIndex]
			],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgssLoop*)


currentJgssLoop[_,nLoops_Integer,_]/;nLoops <= 0 || nLoops == 2 := {};

currentJgssLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,vertexIndex=Unique[],graphExtension,protoDiagrams},

	protoDiagrams=CurrentJs[Prepend[particles,{"s",Unique[internalParticle],vertexIndex}],nLoops-1,vertexIndex];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgSSLoop*)


currentJgSSLoop[_,nLoops_Integer,_]/;nLoops <= 0  || nLoops == 1:= {};

currentJgSSLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,vertexIndex=Unique[],graphExtension,protoDiagrams},

	protoDiagrams=CurrentJS[Prepend[particles,{"S",Unique[internalParticle],vertexIndex}],nLoops-1,vertexIndex];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJggssLoop*)


currentJggssLoop[_,nLoops_,_]/;nLoops <= 0 || nLoops == 2:= {};

currentJggssLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i]],loopAssignment[[1]],vertexIndex],
				CurrentJs[Prepend[particles[[i+1;;]],{"s",Unique[internalParticle],vertexIndex}],loopAssignment[[2]]-1,vertexIndex]
			],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection:: *)
(*currentJgssgLoop*)


currentJgssgLoop[_,nLoops_,_]/;nLoops <= 0 || nLoops == 2:= {};

currentJgssgLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[i+1;;]],loopAssignment[[2]],vertexIndex],
				CurrentJs[Prepend[particles[[;;i]],{"s",Unique[internalParticle],vertexIndex}],loopAssignment[[1]]-1,vertexIndex]
			],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJggSSLoop*)


currentJggSSLoop[_,nLoops_,_]/;nLoops <= 0 || nLoops == 1:= {};

currentJggSSLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i]],loopAssignment[[1]],vertexIndex],
				CurrentJS[Prepend[particles[[i+1;;]],{"S",Unique[internalParticle],vertexIndex}],loopAssignment[[2]]-1,vertexIndex]
			],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsubsection::Closed:: *)
(*currentJgSSgLoop*)


currentJgSSgLoop[_,nLoops_,_]/; nLoops <= 0 || nLoops == 1 := {};

currentJgSSgLoop[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[i+1;;]],loopAssignment[[2]],vertexIndex],
				CurrentJS[Prepend[particles[[;;i]],{"S",Unique[internalParticle],vertexIndex}],loopAssignment[[1]]-1,vertexIndex]
			],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Thick},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsection::Closed:: *)
(*Scalar S Current*)


CurrentJS[{},_,_] := {};

CurrentJS[_,nLoops_,_]/; nLoops < 0 := {};
CurrentJS[{{_,_}},nLoops_,_]/; nLoops > 0 := {};
CurrentJS[{{_,_,___}},_,_] := {};
CurrentJS[{{"S",label_,___}},0,prevVertexIndex_] := Module[{endVertexIndex=Unique[]},
	{
		{{endVertexIndex},{UndirectedEdge[endVertexIndex,prevVertexIndex]},{UndirectedEdge[endVertexIndex,prevVertexIndex]->Dashed},{endVertexIndex->ToString@label}}
	}
];

CurrentJS[{{"S",label_,loopVertexIndex_}},0,prevVertexIndex_] := 
{
	{{},{UndirectedEdge[loopVertexIndex,prevVertexIndex]},{UndirectedEdge[loopVertexIndex,prevVertexIndex]->Dashed},{}}
};

CurrentJS[particles_List,nLoops_Integer,prevVertexIndex_]/; Length@particles >1 :=
Join[
	currentJSSg[particles,nLoops,prevVertexIndex],
	currentJSgS[particles,nLoops,prevVertexIndex],
	currentJSSgg[particles,nLoops,prevVertexIndex],
	currentJSggS[particles,nLoops,prevVertexIndex],
	currentJSSss[particles,nLoops,prevVertexIndex],
	currentJSssS[particles,nLoops,prevVertexIndex]
];


currentJSSg[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJS[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJg[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


currentJSgS[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJg[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJS[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJSSgg]
currentJSSgg[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJSSgg[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJS[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJg[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJg[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJSggS]
currentJSggS[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJSggS[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJg[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJS[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJSSss]
currentJSSss[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJSSss[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJS[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJs[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJs[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJSssS]
currentJSssS[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJSssS[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJs[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJs[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJS[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dashed},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Subsection:: *)
(*Scalar s Current*)


CurrentJs[{},_,_] := {};

CurrentJs[_,nLoops_,_]/; nLoops < 0 := {};
CurrentJs[{{_,_}},nLoops_,_]/; nLoops > 0 := {};
CurrentJs[{{_,_,___}},_,_] := {};
CurrentJs[{{"s",label_}},0,prevVertexIndex_] := Module[{endVertexIndex=Unique[]},
	{
		{{endVertexIndex},{UndirectedEdge[endVertexIndex,prevVertexIndex]},{UndirectedEdge[endVertexIndex,prevVertexIndex]->Dashed},{endVertexIndex->ToString@label}}
	}
];

CurrentJs[{{"s",label_,loopVertexIndex_}},0,prevVertexIndex_] := 
{
	{{},{UndirectedEdge[loopVertexIndex,prevVertexIndex]},{UndirectedEdge[loopVertexIndex,prevVertexIndex]->Dotted},{}}
};

CurrentJs[particles_List,nLoops_Integer,prevVertexIndex_]/; Length@particles >1 :=
Join[
	currentJssg[particles,nLoops,prevVertexIndex],
	currentJsgs[particles,nLoops,prevVertexIndex],
	currentJssgg[particles,nLoops,prevVertexIndex],
	currentJsggs[particles,nLoops,prevVertexIndex],
	currentJssSS[particles,nLoops,prevVertexIndex],
	currentJsSSs[particles,nLoops,prevVertexIndex]
];


currentJssg[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJs[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJg[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


currentJsgs[particles_,nLoops_Integer,prevVertexIndex_]:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,2];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[CurrentJg[particles[[i+1;;]],loopAssignment[[1]],vertexIndex],CurrentJs[particles[[1;;i]],loopAssignment[[2]],vertexIndex]],
		{i,1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3}];
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJssgg]
currentJssgg[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJssgg[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJs[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJg[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJg[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJsggs]
currentJsggs[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJsggs[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJg[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJg[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJs[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJssSS]
currentJssSS[particles_,nLoops_Integer,prevVertexIndex_] := {};


currentJssSS[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJs[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJS[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJS[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


ClearAll[currentJsSSs]
currentJsSSs[particles_,nLoops_Integer,prevVertexIndex_] := {};

currentJsSSs[particles_,nLoops_Integer,prevVertexIndex_] /; Length@particles >= 3:=Module[{n,loopAssignments,vertexIndex=Unique[],protoDiagrams,graphExtension},
	loopAssignments = distributeLoops[nLoops,3];
	n=Length@particles;
	protoDiagrams=Flatten[
		Table[
			mergeCurrents[
				CurrentJS[particles[[1;;i1]],loopAssignment[[1]],vertexIndex],
				CurrentJS[particles[[i1+1;;i2]],loopAssignment[[2]],vertexIndex],
				CurrentJs[particles[[i2+1;;]],loopAssignment[[3]],vertexIndex]
			],
		{i1,1,n-2},{i2,i1+1,n-1},{loopAssignment,loopAssignments}],
	{1,2,3,4}];
	
	graphExtension = {
		{vertexIndex},
		{UndirectedEdge[prevVertexIndex,vertexIndex]},
		{UndirectedEdge[prevVertexIndex,vertexIndex]->Dotted},
		{}
	};
	
	Map[Join[graphExtension,#,2]&][protoDiagrams]
];


(* ::Section:: *)
(*Main Feynman Diagrams Function*)


ClearAll[GenerateFeynmanDiagrams]


GenerateFeynmanDiagrams[particles_,nloops_Integer] /; particles[[1,1]] === "g" :=Module[{firstParticleIndex=Unique[],firstVertexIndex=Unique[]},
	Map[Join[{{firstParticleIndex},{},{},{firstParticleIndex -> ToString[particles[[1,2]]]}},#,2]&][CurrentJg[Rest@particles,nloops,firstParticleIndex]]
]


GenerateFeynmanDiagrams[{{"g",4},{"g",5},{"S",2},{"S",3}},0]
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@%


CurrentJs[{{"s",1},{"g",2}},-1,1] 


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4}},2]//Length
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4}},2],ScalesLoopFreeQ[#[[2]]]&]


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
		Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",connect}},1],ScalesLoopFreeQ[#[[2]]]&]


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
		Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5}},2],ScalesLoopFreeQ[#[[2]]]&]


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
	Select[GenerateFeynmanDiagrams[{{"g",in},{"g",2},{"g",3},{"g",4},{"g",out}},2],ScalesLoopFreeQ[#[[2]]]&]


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4}},1]


GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4}},0]


%//Length
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5},{"g",6}},0]
%//Length
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5},{"g",6},{"g",7}},0];
%//Length


$SystemWordLength


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5},{"g",6},{"g",7},{"g",8}},0];
%//Length


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
	Select[GenerateFeynmanDiagrams[{{"g","in"},{"g",2},{"g",3},{"g",4},{"g","out"}},2],ScalesLoopFreeQ[#[[2]]]&]//Length


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
	Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5}},0],ScalesLoopFreeQ[#[[2]]]&]//Length
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
	Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5}},1],ScalesLoopFreeQ[#[[2]]]&]//Length
Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@
	Select[GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5}},1],ScalesLoopFreeQ[#[[2]]]&]//Length


Graph[#[[2]],VertexLabels->#[[4]],EdgeStyle->#[[3]]]&/@GenerateFeynmanDiagrams[{{"g",1},{"g",2},{"g",3},{"g",4},{"g",5},{"g",6},{"g",7},{"g",8}},0]//Length
