(* Mathematica Test File *)
(* Created by the Wolfram Language Plugin for IntelliJ, see http://wlplugin.halirutan.de/ *)

<<SpinorHelicity6D`Amplitudes`;
<<SpinorHelicity6D`Numerics`;

Block[{WriteString},LoadAmplitudesFile[];];
KillMasses[{kBubRef,q,$bubbleReference}];
GenSpinors[{q},FourD->{q},RandomSpinors->True,Seed->235,ParameterRange->50,DisplaySpinors->False];

BeginTestSection["CollinearSplittingFunctions"]

VerificationTest[(* 1 *)
  NewProcess[];
  KillMasses[{Tilde@1,Tilde@2,3,4,56}];
  GenSpinors[{Tilde@1,Tilde@2,3,4,56},FourD->{Tilde@1,Tilde@2,3,4,56},Seed->1234,ParameterRange->100];
  GenCollinear[56->{5,6},{Tilde@1,Tilde@2}->{1,2},AngleVar->ArcSin[3/5],Reference->q];


  SeriesCoefficient[
    SeriesCoefficient[A52LRatLitBadger[1,2,3,4,56]*SplittingFunctionTree[Cos[ArcSin[3/5]]^2,5,6]["-"][{"+","+"}]+
      AOneLoop["+","+","+","+","-"][1,2,3,4,56]*SplittingFunctionOneLoop[Cos[ArcSin[3/5]]^2,5,6]["+"][{"+","+"}]+
      AOneLoop["+","+","+","+","+"][1,2,3,4,56]*SplittingFunctionOneLoop[Cos[ArcSin[3/5]]^2,5,6]["-"][{"+","+"}]/.{Log[_]:>0}//ToNum,
    {$NumericsEpsilonDimReg,0,0}]/.{Log[_]:>0,Pi->0},
  {$NumericsEpsilonLimit,0,-1}]
  ,
  SeriesCoefficient[ToNum[A62LRatLitBadger[1,2,3,4,5,6]],{$NumericsEpsilonLimit,0,-1}],
  TestID->"6-pt Two-Loop All-Plus Rationals"
]

EndTestSection[]
