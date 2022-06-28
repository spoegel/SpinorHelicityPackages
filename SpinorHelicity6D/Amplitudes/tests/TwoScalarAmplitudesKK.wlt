(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)



BeginTestSection["TwoScalarAmplitudesKK"]

VerificationTest[(* 1 *)
	VerifyKK2ScalarAmplitudes[4]
	,
	List[
		List[0, 0, 0]
	]
	,
	TestID->"Four-Point Two Scalars"
];

VerificationTest[(* 2 *)
	VerifyKK2ScalarAmplitudes[5]
	,
	List[
		List[0, 0, 0, 0],
		List[0, 0, 0, 0]
	]
	,
	TestID->"Five-Point Two Scalars"
];

VerificationTest[(* 3 *)
	VerifyKK2ScalarAmplitudes[6]
	,
	List[
		List[0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0]
	]
	,
	TestID->"Six-Point Two Scalars"
];

VerificationTest[(* 4 *)
	VerifyKK2ScalarAmplitudes[7]
	,
	List[
		List[0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Seven-Point Two Scalars"
];

EndTestSection[]
