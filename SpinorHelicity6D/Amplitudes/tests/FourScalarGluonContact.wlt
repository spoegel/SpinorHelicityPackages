(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)



BeginTestSection["FourScalarGluonContact"]


VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","s","s"][1,2,3,4]+ATreeContact["S","S","s","s"][1,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATree["S","S","s","s"][1,2,3,4]
	}
	,
	TestID->"A4SSss"
];

VerificationTest[(* 1 *)
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATreeGluon["S","s","S","s"][1,3,2,4]+ATreeContact["S","s","S","s"][1,3,2,4]
	}
	,
	GenFourScalarMomenta[{1,2},{3,4},{},Seed->1234];
	ToNum@{
		ATree["S","s","S","s"][1,3,2,4]
	}
	,
	TestID->"A4SsSs"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","+","s","s"][1,2,3,4,5]+ATreeContact["S","S","+","s","s"][1,2,3,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","S","+","s","s"][1,2,3,4,5]
	}
	,
	TestID->"A5SS+ss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","+","S","s","s"][1,3,2,4,5]+ATreeContact["S","+","S","s","s"][1,3,2,4,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","+","S","s","s"][1,3,2,4,5]
	}
	,
	TestID->"A5S+Sss"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","s","+","s"][1,2,4,3,5]+ATreeContact["S","S","s","+","s"][1,2,4,3,5]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","S","s","+","s"][1,2,4,3,5]
	}
	,
	TestID->"A5SSs+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATreeGluon["S","s","S","+","s"][1,3,2,3,4]+ATreeContact["S","s","S","+","s"][1,3,2,3,4]
	}
	,
	GenFourScalarMomenta[{1,2},{4,5},{3},Seed->1234];
	ToNum@{
		ATree["S","s","S","+","s"][1,3,2,3,4]
	}
	,
	TestID->"A5SsS+s"
];

VerificationTest[(* 2 *)
	GenFourScalarMomenta[{1,2},{5,6},{3,4},Seed->1234];
	ToNum@{
		ATreeGluon["S","S","+","+","s","s"][1,2,3,4,5,6]+
      ATreeContact["S","S","+","+","s","s"][1,2,3,4,5,6]
	}
	,
	GenFourScalarMomenta[{1,2},{5,6},{3,4},Seed->1234];
	ToNum@{
		ATree["S","S","+","+","s","s"][1,2,3,4,5,6]
	}
	,
	TestID->"A6SS++ss"
];

EndTestSection[]
