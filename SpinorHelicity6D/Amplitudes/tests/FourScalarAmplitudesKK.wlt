(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)



BeginTestSection["FourScalarAmplitudesKK"]

VerificationTest[(* 1 *)
	VerifyKK4ScalarAmplitudes[4]
	,
	List[
		List[0, 0, 0]
	]
	,
	TestID->"Four-Point Four Scalars"
];

VerificationTest[(* 2 *)
	VerifyKK4ScalarAmplitudes[5]
	,
	List[
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Five-Point Four Scalars"
];

VerificationTest[(* 3 *)
	VerifyKK4ScalarAmplitudes[6]
	,
	List[
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Six-Point Four Scalars"
];

VerificationTest[ (* 4 *)
	VerifyKK4ScalarAmplitudes[7]
	,
	List[
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Seven-Point Four Scalars"
];

VerificationTest[(* 5 *)
	VerifyKK4ScalarAmplitudes[4,TreeHead->ATreeGluon]
	,
	List[
		List[0, 0, 0]
	]
	,
	TestID->"Four-Point Four Scalars Gluon Exchange"
];

VerificationTest[(* 5 *)
	VerifyKK4ScalarAmplitudes[4,TreeHead->ATreeContact]
	,
	List[
		List[0, 0, 0]
	]
	,
	TestID->"Four-Point Four Scalars Contact Term"
];

VerificationTest[(* 6 *)
	VerifyKK4ScalarAmplitudes[5,TreeHead->ATreeGluon]
	,
	List[
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Five-Point Four Scalars Gluon Exchange"
];

VerificationTest[(* 7 *)
	VerifyKK4ScalarAmplitudes[5,TreeHead->ATreeContact]
	,
	List[
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	]
	,
	TestID->"Five-Point Four Scalars Contact Term"
];

EndTestSection[]
