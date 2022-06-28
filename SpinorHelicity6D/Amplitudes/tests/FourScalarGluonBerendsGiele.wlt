(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

BeginTestSection["FourScalarGluonContact"]


VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","s","s"][1,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","S","s","s"][1,2,3,4]//RootReduce
	}
	,
	TestID->"A4SSss"
];

VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeGluon["S","s","S","s"][1,3,2,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","s","S","s"][1,3,2,4]//RootReduce
	}
	,
	TestID->"A4SsSs"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","+","s","s"][1,2,3,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","S","+","s","s"][1,2,3,4,5]//RootReduce
	}
	,
	TestID->"A5SS+ss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","+","S","s","s"][1,3,2,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","+","S","s","s"][1,3,2,4,5]//RootReduce
	}
	,
	TestID->"A5S+Sss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","s","+","s"][1,2,4,3,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","S","s","+","s"][1,2,4,3,5]//RootReduce
	}
	,
	TestID->"A5SSs+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","s","S","+","s"][1,3,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeBerendsGieleGluon["S","s","S","+","s"][1,3,2,3,4]//RootReduce
	}
	,
	TestID->"A5SsS+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3,6},Seed->1234];
	ToNum@{
		ATreeContact["S","S","+","s","s","+"][1,2,3,4,5,6]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3,6},Seed->1234];
	ToNum@{
		ATreeBerendsGieleContact["S","S","+","s","s","+"][1,2,3,4,5,6]//RootReduce
	}
	,
	TestID->"A6SS+ss+"
];


EndTestSection[]