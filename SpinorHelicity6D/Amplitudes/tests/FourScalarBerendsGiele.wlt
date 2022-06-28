(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)


BeginTestSection["FourScalarGluonContact"]


VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATree["S","S","s","s"][1,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","S","s","s"][1,2,3,4]//RootReduce
	}
	,
	TestID->"A4SSss"
];

VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATree["S","s","S","s"][1,3,2,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","s","S","s"][1,3,2,4]//RootReduce
	}
	,
	TestID->"A4SsSs"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","S","+","s","s"][1,2,3,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","S","+","s","s"][1,2,3,4,5]//RootReduce
	}
	,
	TestID->"A5SS+ss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","+","S","s","s"][1,3,2,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","+","S","s","s"][1,3,2,4,5]//RootReduce
	}
	,
	TestID->"A5S+Sss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","S","s","+","s"][1,2,4,3,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","S","s","+","s"][1,2,4,3,5]//RootReduce
	}
	,
	TestID->"A5SSs+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	Echo@ToNum@{
		ATree["S","s","S","+","s"][1,3,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","s","S","+","s"][1,3,2,3,4]//RootReduce
	}
	,
	TestID->"A5SsS+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{7,8},{3,4,5,6},Seed->1234];
	ToNum@{
		ATree["S","s","+","+","+","+","S","s"][1,2,3,4,5,6,7,8],
		ATree["S","+","+","+","+","s","S","s"][2,3,4,5,6,7,8,1],
		ATree["S","s","S","+","+","+","+","s"][8,1,2,3,4,5,6,7],
		ATree["S","s","S","s","+","+","+","+"][7,8,1,2,3,4,5,6]
	}
	,
	GenFourScalarMomenta[{1,2},{7,8},{3,4,5,6},Seed->1234];
	ToNum@{
		ATreeBerendsGiele["S","s","+","+","+","+","S","s"][1,2,3,4,5,6,7,8]//RootReduce,
		ATreeBerendsGiele["S","+","+","+","+","s","S","s"][2,3,4,5,6,7,8,1]//RootReduce,
		ATreeBerendsGiele["S","s","S","+","+","+","+","s"][8,1,2,3,4,5,6,7]//RootReduce,
		ATreeBerendsGiele["S","s","S","s","+","+","+","+"][7,8,1,2,3,4,5,6]//RootReduce
	}
	,
	TestID->"A8Ss++++Ss Family"
];



EndTestSection[]
