ATreeSAM["S", "+", "+", "S"][p1_, p2_, p3_, p4_] = 
    -1/4*(MP[p1, p1]*Spbb[p2, p3]^2)/(MP[p2, p3]*MP[p3, p4])
 
ATreeSAM["S", "+", "-", "S"][p1_, p2_, p3_, p4_] = 
    -1/4*(Spab[p3, p1, p2]*Spab[p3, p4, p2])/(MP[p2, p3]*MP[p3, p4])
 
ATreeSAM["S", "+", "S", "+"][p1_, p2_, p3_, p4_] = 
    -1/4*(MP[p1, p1]*Spbb[p2, p4]^2)/(MP[p1, p2]*MP[p2, p3])
 
ATreeSAM["S", "+", "+", "+", "S"][p1_, p2_, p3_, p4_, p5_] = 
    (MP[p5, p5]*Spbb[p4, SumM[p2, p3], p1, p2])/(4*MP[p1, p2]*MP[p4, p5]*
      Spaa[p2, p3]*Spaa[p3, p4])
 
ATreeSAM["S", "+", "+", "-", "S"][p1_, p2_, p3_, p4_, p5_] = 
    -((Spab[p4, p5, SumM[p2, p3], p1, p2]^2/(4*MP[p1, p2]*MP[p4, p5]*
         Spaa[p2, p3]*Spaa[p3, p4]) - (MP[p5, p5]*Spbb[p2, p3]^3)/
        ((2*MP[p1, p1] + 2*MP[p1, p5])*Spbb[p3, p4]))/
      Spbb[p4, SumM[p2, p3], p1, p2])
 
ATreeSAM["S", "+", "-", "+", "S"][p1_, p2_, p3_, p4_, p5_] = 
    (-1/4*(Spab[p3, p1, p2]^2*Spab[p3, p5, p4]^2)/(MP[p1, p2]*MP[p4, p5]*
         Spaa[p2, p3]*Spaa[p3, p4]) + (MP[p1, p1]*Spbb[p2, p4]^4)/
       ((2*MP[p1, p1] + 2*MP[p1, p5])*Spbb[p2, p3]*Spbb[p3, p4]))/
     Spbb[p4, SumM[p2, p3], p1, p2]
 
ATreeSAM["S", "+", "+", "S", "+"][p1_, p2_, p3_, p4_, p5_] = 
    (MP[p4, p4]*Spbb[p3, SumM[p2, p5], p1, p2])/(4*MP[p1, p2]*MP[p3, p4]*
       Spaa[p2, p5]*Spaa[p3, p5]) + 
     (MP[p4, p4]*Spbb[p3, SumM[p2, p5], p1, p5])/(4*MP[p1, p5]*MP[p3, p4]*
       Spaa[p2, p3]*Spaa[p2, p5]) - 
     (MP[p4, p4]*Spbb[p5, SumM[p2, p3], p1, p2])/(4*MP[p1, p2]*MP[p4, p5]*
       Spaa[p2, p3]*Spaa[p3, p5])
 
ATreeSAM["S", "+", "+", "+", "+", "S"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    -1/2*(MP[p6, p6]*Spbb[p5, p6, SumM[p4, p5], SumM[p2, p3], p1, p2])/
      (MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p4] + 2*MP[p2, p3] + 
        2*MP[p2, p4] + 2*MP[p3, p4])*Spaa[p2, p3]*Spaa[p3, p4]*Spaa[p4, p5])
 
ATreeSAM["S", "+", "+", "+", "S", "+"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    -1/2*(MP[p5, p5]*Spbb[p4, p5, SumM[p3, p4], SumM[p2, p6], p1, p2])/
       (MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
        (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
         2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p6]*Spaa[p3, p4]*
        Spaa[p3, p6]) - (MP[p5, p5]*Spbb[p4, p5, SumM[p3, p4], SumM[p2, p6], 
        p1, p6])/(2*MP[p1, p6]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
        2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p3]*Spaa[p2, p6]*
       Spaa[p3, p4]) - (MP[p5, p5]*Spbb[p4, p5, SumM[p4, p6], SumM[p2, p3], 
        p1, p2])/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
        2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p3]*Spaa[p3, p6]*
       Spaa[p4, p6]) + (MP[p5, p5]*Spbb[p6, p5, SumM[p4, p6], SumM[p2, p3], 
        p1, p2])/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p4] + 2*MP[p2, p3] + 
        2*MP[p2, p4] + 2*MP[p3, p4])*Spaa[p2, p3]*Spaa[p3, p4]*Spaa[p4, p6])
 
ATreeSAM["S", "+", "+", "S", "+", "+"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    -1/2*(MP[p4, p4]*Spbb[p3, p4, SumM[p2, p3], SumM[p5, p6], p1, p6])/
       (MP[p1, p6]*(2*MP[p1, p5] + 2*MP[p1, p6] + 2*MP[p5, p6])*
        (2*MP[p1, p2] + 2*MP[p1, p5] + 2*MP[p1, p6] + 2*MP[p2, p5] + 
         2*MP[p2, p6] + 2*MP[p5, p6])*Spaa[p2, p3]*Spaa[p2, p5]*
        Spaa[p5, p6]) - (MP[p4, p4]*Spbb[p3, p4, SumM[p3, p5], SumM[p2, p6], 
        p1, p2])/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
       (2*MP[p1, p2] + 2*MP[p1, p5] + 2*MP[p1, p6] + 2*MP[p2, p5] + 
        2*MP[p2, p6] + 2*MP[p5, p6])*Spaa[p2, p6]*Spaa[p3, p5]*
       Spaa[p5, p6]) - (MP[p4, p4]*Spbb[p3, p4, SumM[p3, p5], SumM[p2, p6], 
        p1, p6])/(2*MP[p1, p6]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
       (2*MP[p1, p2] + 2*MP[p1, p5] + 2*MP[p1, p6] + 2*MP[p2, p5] + 
        2*MP[p2, p6] + 2*MP[p5, p6])*Spaa[p2, p5]*Spaa[p2, p6]*
       Spaa[p3, p5]) + (MP[p4, p4]*Spbb[p5, p4, SumM[p3, p5], SumM[p2, p6], 
        p1, p2])/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
        2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p6]*Spaa[p3, p5]*
       Spaa[p3, p6]) + (MP[p4, p4]*Spbb[p5, p4, SumM[p3, p5], SumM[p2, p6], 
        p1, p6])/(2*MP[p1, p6]*(2*MP[p1, p2] + 2*MP[p1, p6] + 2*MP[p2, p6])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
        2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p3]*Spaa[p2, p6]*
       Spaa[p3, p5]) + (MP[p4, p4]*Spbb[p5, p4, SumM[p5, p6], SumM[p2, p3], 
        p1, p2])/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p1, p6] + 2*MP[p2, p3] + 
        2*MP[p2, p6] + 2*MP[p3, p6])*Spaa[p2, p3]*Spaa[p3, p6]*Spaa[p5, p6])
 
ATreeSAM["S", "+", "+", "+", "-", "S"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    -((MP[p1, p1]*Spab[p5, SumM[p3, p4], p2]^3)/
       ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*(2*MP[p2, p3] + 
         2*MP[p2, p4] + 2*MP[p2, p5] + 2*MP[p3, p4] + 2*MP[p3, p5] + 
         2*MP[p4, p5])*Spaa[p3, p4]*Spaa[p4, p5]*Spab[p3, SumM[p4, p5], 
         SumM[p1, p6], p1, p2])) + 
     (-(MP[p1, p1]*Spab[p5, p6, p4, p3, p2]) + 
        (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p5, p6, 
          SumM[p2, p3, p4], p1, p2])^2/(2*MP[p1, p2]*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*(2*MP[p1, p2] + 
        2*MP[p1, p3] + 2*MP[p1, p4] + 2*MP[p2, p3] + 2*MP[p2, p4] + 
        2*MP[p3, p4])*Spaa[p2, p3]*Spaa[p3, p4]*Spaa[p4, p5]*
       Spbb[p5, p6, SumM[p4, p5], SumM[p2, p3], p1, p2]) - 
     (MP[p1, p1]*Spbb[p4, SumM[p2, p3], p1, p2]^3)/(2*MP[p1, p2]*Spaa[p2, p3]*
       Spab[p3, SumM[p4, p5], SumM[p1, p6], p1, p2]*Spbb[p4, p5]*
       Spbb[p5, p6, SumM[p4, p5], SumM[p2, p3], p1, p2])
 
ATreeSAM["S", "+", "+", "-", "+", "S"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    -((MP[p1, p1]*Spab[p4, SumM[p3, p5], p2]^4)/
       ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*(2*MP[p2, p3] + 
         2*MP[p2, p4] + 2*MP[p2, p5] + 2*MP[p3, p4] + 2*MP[p3, p5] + 
         2*MP[p4, p5])*Spaa[p3, p4]*Spaa[p4, p5]*Spab[p5, SumM[p3, p4], p2]*
        Spab[p3, SumM[p4, p5], SumM[p1, p6], p1, p2])) + 
     (MP[p1, p1]*Spbb[p2, p3]^3*Spbb[p5, SumM[p2, p3, p4], p1, p2])/
      ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*(2*MP[p1, p2] + 
        2*MP[p1, p3] + 2*MP[p1, p4] + 2*MP[p2, p3] + 2*MP[p2, p4] + 
        2*MP[p3, p4])*Spab[p5, SumM[p3, p4], p2]*Spbb[p3, p4]*
       Spbb[p4, SumM[p2, p3], p1, p2]) - 
     (Spbb[p2, p3]^3*(-1/2*(Spab[p4, p6, p5]^2*Spab[p4, SumM[p1, p2, p3], 
             SumM[p2, p3], p1, p2]^2*Spbb[p2, p3, p1, p2])/
          ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*MP[p5, p6]*
           Spaa[p4, p5]*Spab[p4, p3, p2]) - 
        (MP[p1, p1]*Spbb[p5, SumM[p2, p3], p1, p2]^4)/(Spbb[p4, p5]*
          (2*MP[p1, p1] + 2*MP[p1, p6] + (2*MP[p2, p3]*Spbb[p2, p6, p1, p2])/
            Spbb[p2, p3, p1, p2])*Spbb[p4, SumM[p2, p3], p1, p2])))/
      (2*MP[p2, p3]*Spbb[p2, p3, p1, p2]*Spbb[p3, p2, p1, p2]*
       ((2*MP[p2, p3]*Spbb[p2, p5]*Spbb[p2, p1, SumM[p1, p2, p3], 
           SumM[p2, p3], p1, p2])/Spbb[p2, p1, p3, p2] + 
        Spbb[p5, p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2] + 
        Spbb[p5, SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], p1, p2]))
 
ATreeSAM["S", "+", "+", "+", "+", "+", "S"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p7, p7]*Spbb[p2, p3]^3*
      (2*MP[p2, p3]*(MP[p1, p1]*Spbb[p2, SumM[p5, p6], p7, p6] + 
         Spbb[p2, p1, p4, SumM[p5, p6], p7, p6] + Spbb[p2, p1, SumM[p2, p3], 
          SumM[p5, p6], p7, p6]) + Spbb[p2, p1, SumM[p2, p3], p1, p4, 
        SumM[p5, p6], p7, p6] + Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], 
        SumM[p5, p6], p7, p6]))/(4*MP[p2, p3]*(2*MP[p1, p2] + 2*MP[p1, p3] + 
       2*MP[p2, p3])*MP[p6, p7]*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*
      Spaa[p4, p5]*Spaa[p5, p6]*Spab[p4, p3, p2]*Spbb[p2, p1, p2, p3])
 
ATreeSAM["S", "+", "+", "+", "+", "S", "+"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p6, p6]*Spbb[p2, p3]^3*
       (2*MP[p2, p3]*(MP[p1, p1]*Spbb[p2, SumM[p4, p5], p6, p5] + 
          Spbb[p2, p1, p7, SumM[p4, p5], p6, p5] + Spbb[p2, p1, SumM[p2, p3], 
           SumM[p4, p5], p6, p5]) + Spbb[p2, p1, SumM[p2, p3], p1, p7, 
         SumM[p4, p5], p6, p5] + Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], 
         SumM[p4, p5], p6, p5]))/(4*MP[p2, p3]*(2*MP[p1, p2] + 2*MP[p1, p3] + 
        2*MP[p2, p3])*MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
       Spaa[p4, p5]*Spaa[p4, p7]*Spab[p7, p3, p2]*Spbb[p2, p1, p2, p3]) + 
     (MP[p6, p6]*Spbb[p2, p3]^3*(2*MP[p2, p3]*
         (MP[p1, p1]*Spbb[p2, SumM[p5, p7], p6, p5] + Spbb[p2, p1, p4, 
           SumM[p5, p7], p6, p5] + Spbb[p2, p1, SumM[p2, p3], SumM[p5, p7], 
           p6, p5]) + Spbb[p2, p1, SumM[p2, p3], p1, p4, SumM[p5, p7], p6, 
         p5] + Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], SumM[p5, p7], p6, 
         p5]))/(4*MP[p2, p3]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       MP[p5, p6]*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*Spaa[p4, p7]*
       Spaa[p5, p7]*Spab[p4, p3, p2]*Spbb[p2, p1, p2, p3]) - 
     (MP[p6, p6]*Spbb[p2, p3]^3*(2*MP[p2, p3]*
         (MP[p1, p1]*Spbb[p2, SumM[p5, p7], p6, p7] + Spbb[p2, p1, p4, 
           SumM[p5, p7], p6, p7] + Spbb[p2, p1, SumM[p2, p3], SumM[p5, p7], 
           p6, p7]) + Spbb[p2, p1, SumM[p2, p3], p1, p4, SumM[p5, p7], p6, 
         p7] + Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], SumM[p5, p7], p6, 
         p7]))/(4*MP[p2, p3]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       MP[p6, p7]*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*Spaa[p4, p5]*
       Spaa[p5, p7]*Spab[p4, p3, p2]*Spbb[p2, p1, p2, p3]) - 
     (MP[p6, p6]*Spbb[p2, p7]^3*(2*MP[p2, p7]*
         (MP[p1, p1]*Spbb[p2, SumM[p4, p5], p6, p5] + Spbb[p2, p1, p3, 
           SumM[p4, p5], p6, p5] + Spbb[p2, p1, SumM[p2, p7], SumM[p4, p5], 
           p6, p5]) + Spbb[p2, p1, SumM[p2, p7], p1, p3, SumM[p4, p5], p6, 
         p5] + Spbb[p2, p1, SumM[p2, p7], p1, SumM[p2, p7], SumM[p4, p5], p6, 
         p5]))/(4*MP[p2, p7]*(2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
       MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p3, p4]*
       Spaa[p4, p5]*Spab[p3, p7, p2]*Spbb[p2, p1, p2, p7]) + 
     (MP[p6, p6]*Spbb[p2, p7]^3*(2*MP[p2, p7]*
         (MP[p1, p1]*Spbb[p7, SumM[p4, p5], p6, p5] + Spbb[p7, p1, p3, 
           SumM[p4, p5], p6, p5] + Spbb[p7, p1, SumM[p2, p7], SumM[p4, p5], 
           p6, p5]) + Spbb[p7, p1, SumM[p2, p7], p1, p3, SumM[p4, p5], p6, 
         p5] + Spbb[p7, p1, SumM[p2, p7], p1, SumM[p2, p7], SumM[p4, p5], p6, 
         p5]))/(4*MP[p2, p7]*(2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
       MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p3, p4]*
       Spaa[p4, p5]*Spab[p3, p2, p7]*Spbb[p7, p1, p7, p2])
 
ATreeSAM["S", "+", "+", "+", "S", "+", "+"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/4*(MP[p1, p1]^2*Spbb[p6, p7]^2*
        ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*Spbb[p2, p4] - 
         Spbb[p2, SumM[p3, p4, p5], p2, p4] - Spbb[p2, SumM[p3, p4, p5], p3, 
          p4]))/(MP[p4, p5]*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
        MP[p6, p7]*(2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*Spaa[p3, p4]*
        (Spaa[p2, p3] + ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*
           Spab[p3, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])*
        (2*MP[p1, p7] + ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*
           Spbb[p2, p7, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])) - 
     (MP[p1, p1]^2*Spbb[p2, p7]*
       (-1/2*(2*MP[p1, p7]*Spbb[p2, SumM[p4, p6], p5, p4] + 
           Spbb[p2, SumM[p1, p7], p2, SumM[p4, p6], p5, p4] + 
           Spbb[p2, SumM[p1, p7], p3, SumM[p4, p6], p5, p4])/
          ((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*MP[p4, p5]*
           (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p3, p6]*
           Spaa[p4, p6]*(-Spaa[p2, p3] + (2*MP[p1, p7]*Spab[p3, p1, p2])/
             Spbb[p2, SumM[p1, p7], p1, p2])) + 
        (2*MP[p1, p7]*Spbb[p2, SumM[p4, p6], p5, p6] + Spbb[p2, SumM[p1, p7], 
           p2, SumM[p4, p6], p5, p6] + Spbb[p2, SumM[p1, p7], p3, 
           SumM[p4, p6], p5, p6])/(2*(2*MP[p1, p2] + 2*MP[p1, p7] + 
           2*MP[p2, p7])*MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 
           2*MP[p5, p6])*Spaa[p3, p4]*Spaa[p4, p6]*(-Spaa[p2, p3] + 
           (2*MP[p1, p7]*Spab[p3, p1, p2])/Spbb[p2, SumM[p1, p7], p1, p2])) - 
        (2*MP[p1, p7]*Spbb[p2, SumM[p3, p4], p5, p4] + Spbb[p2, SumM[p1, p7], 
           p2, SumM[p3, p4], p5, p4] + Spbb[p2, SumM[p1, p7], p6, 
           SumM[p3, p4], p5, p4])/(2*(2*MP[p1, p2] + 2*MP[p1, p7] + 
           2*MP[p2, p7])*MP[p4, p5]*(2*MP[p3, p4] + 2*MP[p3, p5] + 
           2*MP[p4, p5])*Spaa[p3, p4]*Spaa[p3, p6]*(-Spaa[p2, p6] + 
           (2*MP[p1, p7]*Spab[p6, p1, p2])/Spbb[p2, SumM[p1, p7], p1, p2])) - 
        (Spbb[p2, SumM[p1, p7], p1, p2]^2*
          (2*MP[p1, p7]*(Spbb[p2, p1, SumM[p1, p7], p6]*Spbb[p2, SumM[p3, 
                p4], p5, p4] + Spbb[p2, p6]*(Spbb[p2, p1, p2, SumM[p3, p4], 
                p5, p4] + Spbb[p2, p1, p6, SumM[p3, p4], p5, p4])) + 
           Spbb[p2, SumM[p1, p7], p1, p2]*(Spbb[p4, p5, SumM[p3, p4], p2, 
              SumM[p1, p7], p6] + Spbb[p4, p5, SumM[p3, p4], p6, 
              SumM[p1, p7], p6])))/(2*MP[p4, p5]*(2*MP[p3, p4] + 
           2*MP[p3, p5] + 2*MP[p4, p5])*Spaa[p3, p4]*
          (2*MP[p1, p7]*Spab[p3, p1, p2] - Spaa[p2, p3]*
            Spbb[p2, SumM[p1, p7], p1, p2])*(2*MP[p1, p7]*Spab[p6, p1, p2] - 
           Spaa[p2, p6]*Spbb[p2, SumM[p1, p7], p1, p2])*
          ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*
            Spbb[p2, SumM[p1, p7], p1, p2] + 2*MP[p1, p7]*
            Spbb[p2, SumM[p3, p4, p5], p1, p2]))))/
      (2*MP[p1, p7]*Spab[p7, p1, p2]) + (MP[p1, p1]*Spbb[p2, p3]^3*
       Spbb[p2, p3, p1, p2]*((-4*MP[p2, p3]^2*Spbb[p2, p7]*
           Spbb[p2, p1, p5, p4]*Spbb[p2, SumM[p6, p7], p1, p2] + 
          2*MP[p2, p3]*Spbb[p2, p3, p1, p2]*
           (-(Spbb[p2, p1, p5, p4]*Spbb[p2, SumM[p6, p7], p1, p7]) + 
            Spbb[p2, p7]*(Spbb[p2, p1, SumM[p6, p7], p4, p5, p4] + 
              Spbb[p2, p1, SumM[p6, p7], SumM[p2, p3], p5, p4])) - 
          Spbb[p2, p3, p1, p2]^2*(Spbb[p4, p5, p4, SumM[p6, p7], p1, p7] + 
            Spbb[p4, p5, SumM[p2, p3], SumM[p6, p7], p1, p7]))/
         (2*MP[p4, p5]*Spaa[p6, p7]*Spab[p4, p3, p2]*Spab[p6, p3, p2]*
          (2*MP[p1, p7]*Spbb[p2, p3, p1, p2] + 2*MP[p2, p3]*
            Spbb[p2, p7, p1, p2])*((2*MP[p1, p6] + 2*MP[p1, p7] + 
             2*MP[p6, p7])*Spbb[p2, p3, p1, p2] + 2*MP[p2, p3]*
            Spbb[p2, SumM[p6, p7], p1, p2])) + 
        (2*MP[p2, p3]*Spbb[p2, p7]*(MP[p1, p1]*Spbb[p2, SumM[p4, p6], p5, 
              p4] + Spbb[p2, p1, p7, SumM[p4, p6], p5, p4] + 
            Spbb[p2, p1, SumM[p2, p3], SumM[p4, p6], p5, p4]) - 
          Spbb[p2, p3, p1, p2]*(Spbb[p4, p5, SumM[p4, p6], p7, p1, p7] + 
            Spbb[p4, p5, SumM[p4, p6], SumM[p2, p3], p1, p7]))/
         (2*MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
          Spaa[p4, p6]*Spab[p6, p3, p2]*Spab[p7, p3, p2]*
          (2*MP[p1, p7]*Spbb[p2, p3, p1, p2] + 2*MP[p2, p3]*
            Spbb[p2, p7, p1, p2])) + (-2*MP[p2, p3]*Spbb[p2, p7]*
           (MP[p1, p1]*Spbb[p2, SumM[p4, p6], p5, p6] + Spbb[p2, p1, p7, 
             SumM[p4, p6], p5, p6] + Spbb[p2, p1, SumM[p2, p3], SumM[p4, p6], 
             p5, p6]) + Spbb[p2, p3, p1, p2]*(Spbb[p6, p5, SumM[p4, p6], p7, 
             p1, p7] + Spbb[p6, p5, SumM[p4, p6], SumM[p2, p3], p1, p7]))/
         (2*MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
          Spaa[p4, p6]*Spab[p4, p3, p2]*Spab[p7, p3, p2]*
          (2*MP[p1, p7]*Spbb[p2, p3, p1, p2] + 2*MP[p2, p3]*
            Spbb[p2, p7, p1, p2])) + (2*MP[p1, p1]*MP[p2, p3]*
           Spbb[p2, p3, p1, p2]*Spbb[p2, SumM[p4, p6], p5, p4] + 
          Spbb[p2, p3, p1, p2]*(2*MP[p2, p3]*Spbb[p2, p1, p7, SumM[p4, p6], 
              p5, p4] + 2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p6], 
              p5, p4] + Spbb[p2, p1, SumM[p2, p3], p1, p7, SumM[p4, p6], p5, 
             p4] + Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], SumM[p4, p6], 
             p5, p4]))/(2*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
          Spaa[p4, p6]*Spaa[p6, p7]*Spab[p7, p3, p2]*Spbb[p2, p3, p1, p2]^
           2) + (-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p3, p1, p2]*
           Spbb[p2, SumM[p4, p6], p5, p6] - Spbb[p2, p3, p1, p2]*
           (2*MP[p2, p3]*Spbb[p2, p1, p7, SumM[p4, p6], p5, p6] + 
            2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p6], p5, p6] + 
            Spbb[p2, p1, SumM[p2, p3], p1, p7, SumM[p4, p6], p5, p6] + 
            Spbb[p2, p1, SumM[p2, p3], p1, SumM[p2, p3], SumM[p4, p6], p5, 
             p6]))/(2*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*MP[p5, p6]*
          (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p4, p6]*
          Spaa[p4, p7]*Spab[p7, p3, p2]*Spbb[p2, p3, p1, p2]^2) + 
        (-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p3, p1, p2]*Spbb[p2, SumM[p6, p7], 
            p5, p6] - Spbb[p2, p3, p1, p2]*(2*MP[p2, p3]*Spbb[p2, p1, p4, 
              SumM[p6, p7], p5, p6] + 2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], 
              SumM[p6, p7], p5, p6] + Spbb[p2, p1, SumM[p2, p3], p1, p4, 
             SumM[p6, p7], p5, p6] + Spbb[p2, p1, SumM[p2, p3], p1, 
             SumM[p2, p3], SumM[p6, p7], p5, p6]))/
         (2*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*MP[p5, p6]*
          (2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*Spaa[p4, p7]*
          Spaa[p6, p7]*Spab[p4, p3, p2]*Spbb[p2, p3, p1, p2]^2)))/
      (2*MP[p2, p3]*Spbb[p2, p1, p2, p3])
 
ATreeSAM["S", "S", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_] = 
    (-((2*(2*MP[p1, p1] + 2*MP[p1, p2])*MP[p1, p4] - 
         2*(2*MP[p1, p1] + 2*MP[p1, p2])*MP[p3, p4] - 
         8*MP[p2, p3]*MP[p3, p4] - 2*MP[p2, p3]*(2*MP[p4, p4] + 
           2*MP[p4, p5]) + 2*MP[p2, p5]*(2*MP[p4, p4] + 2*MP[p4, p5]))*
        Spbb[p3, p2, p4, p3]) + 2*MP[p2, p4]*
       (2*MP[p3, p4]*Spbb[p3, p1, p2, p3] + 2*MP[p2, p3]*
         Spbb[p3, p4, p5, p3]))/(8*(2*MP[p1, p1] + 2*MP[p1, p2])*MP[p2, p3]*
      MP[p3, p4]*(2*MP[p4, p4] + 2*MP[p4, p5]))
 
ATreeSAM["S", "+", "S", "s", "s"][p1_, p2_, p3_, p4_, p5_] = 
    -1/4*(Spbb[p2, p3, p1, p2]*
       (1 + (2*MP[p1, p5] + (2*MP[p2, p3]*Spbb[p2, p5, p1, p2])/
           Spbb[p2, p3, p1, p2])/(2*MP[p4, p4] + 2*MP[p4, p5])))/
      (MP[p1, p2]*MP[p2, p3])
 
ATreeSAM["S", "+", "s", "S", "s"][p1_, p2_, p3_, p4_, p5_] = 
    Spbb[p2, p3, p1, p2]/(4*MP[p1, p2]*MP[p2, p3])
 
ATreeSAM["S", "S", "+", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (MP[p5, p5]*Spbb[p3, p4]^2*(2*MP[p2, p6] - 
        ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
          Spbb[p3, SumM[p1, p4, p5], p2, p3])/Spbb[p3, SumM[p4, p5], p2, 
          p3]))/(2*MP[p4, p5]*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
       (2*MP[p1, p1] + 2*MP[p1, p2] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 
           2*MP[p4, p5])*Spbb[p3, p2, p1, p3])/Spbb[p3, SumM[p4, p5], p2, 
          p3])*(2*MP[p3, p4] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
          Spbb[p3, p4, p2, p3])/Spbb[p3, SumM[p4, p5], p2, p3])) - 
     (MP[p5, p5]*Spbb[p3, SumM[p4, p5, p6], p1, p3]*
       Spbb[p3, SumM[p4, p5, p6], p2, p3]*(Spbb[p3, p2, SumM[p1, p2], p6, p3, 
         p4] - Spbb[p3, p2, SumM[p1, p2], p6, SumM[p3, p5, p6], p4]))/
      (2*(2*MP[p1, p1] + 2*MP[p1, p2])*MP[p4, p5]*(-2*MP[p1, p1] - 
        2*MP[p1, p2] + Spab[p3, SumM[p4, p5, p6], p3])*
       ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, SumM[p4, p5], p2, p3] - 
        (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
         Spbb[p3, SumM[p4, p5, p6], p2, p3])*((2*MP[p1, p1] + 2*MP[p1, p2])*
         Spab[p4, p2, p3] + Spaa[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, 
          p3])) - (Spbb[p3, p2, SumM[p1, p2], p1, p2, p3]*
       ((MP[p5, p5]*Spbb[p3, p4]^3)/(2*MP[p5, p5] + 2*MP[p5, p6]) + 
        ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, p2, p5, p4]*
            Spbb[p3, SumM[p4, p5], p6, p3] + Spbb[p3, SumM[p4, p5, p6], p2, 
             p3]*(Spbb[p3, SumM[p4, p5, p6], p6, p3, p5, p4] + 
             Spbb[p3, SumM[p4, p5, p6], p6, p4, p5, p4]))^2/
         (2*MP[p4, p5]*(-2*MP[p1, p1] - 2*MP[p1, p2] + 
           Spab[p3, SumM[p4, p5, p6], p3])*((2*MP[p1, p1] + 2*MP[p1, p2])*
            Spbb[p3, SumM[p4, p5], p2, p3] - (2*MP[p3, p4] + 2*MP[p3, p5] + 
             2*MP[p4, p5])*Spbb[p3, SumM[p4, p5, p6], p2, p3])*
          ((2*MP[p1, p1] + 2*MP[p1, p2])*Spab[p4, p2, p3] + 
           Spaa[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3]))))/
      ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, SumM[p4, p5, p6], p2, p3]*
       ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, p2, p5, p4] + 
        Spbb[p3, p2, SumM[p1, p2], p3, p5, p4] + Spbb[p3, p2, SumM[p1, p2], 
         p4, p5, p4])) + (Spbb[p3, p4]^3*
       ((-2*MP[p3, p4]*Spbb[p3, p2, p5, p3] + 2*MP[p2, p5]*
           Spbb[p3, p4, p2, p3])*((2*MP[p2, p3] + 2*MP[p2, p4] + 
            2*MP[p3, p4])*Spbb[p3, p4, p2, p3]*Spbb[p3, p2, SumM[p3, p4], p5, 
            p6, SumM[p3, p4], p2, p3] + 
          (-((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*Spbb[p3, p4, p2, 
               p3]) + 2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])*
           Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, SumM[p3, p4], p2, 
            p3]) + (((2*MP[p1, p1] + 2*MP[p1, p2])*(-2*MP[p1, p5] + 
              2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5]) - 
            2*MP[p2, p6]*(2*MP[p5, p5] + 2*MP[p5, p6]) + 
            (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
             (2*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5]) + 2*MP[p5, p5] + 
              2*MP[p5, p6]))*Spbb[p3, p4, p2, p3]^2 + 4*MP[p3, p4]^2*
           Spbb[p3, p2, p1, p3]*Spbb[p3, SumM[p4, p5], p2, p3] + 
          2*MP[p3, p4]*Spbb[p3, p4, p2, p3]*((2*MP[p1, p5] - 2*MP[p3, p4] - 
              2*MP[p3, p5] - 2*MP[p4, p5])*Spbb[p3, p2, p1, p3] + 
            (2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p3, p2, p6, p3] - 
            (2*MP[p1, p1] + 2*MP[p1, p2] + 2*(2*MP[p2, p3] + 2*MP[p2, p4] + 
                2*MP[p3, p4]))*Spbb[p3, SumM[p4, p5], p2, p3]))*
         Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, SumM[p3, p4], p2, 
          p3]))/(4*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
       (2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p3, p2, p3, p4]*
       (2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
         Spbb[p3, p4, p2, p3])*Spbb[p3, p4, p2, p3]^3*
       (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] - 
        (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p4, p2, p3]))
 
ATreeSAM["S", "S", "+", "s", "s", "+"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (MP[p1, p1]*Spbb[p3, p2, SumM[p3, p4, p5], p6]^2*Spbb[p3, p4, p5, p3]^2)/
      (2*MP[p1, p6]*(2*MP[p4, p4] + 2*MP[p4, p5])*(2*MP[p3, p4] + 
        2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*
       ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*
         Spbb[p3, p2, p4, p3] + 2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])*
       (-((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, SumM[p4, p5], p2, p3]) + 
        (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*
         Spbb[p3, SumM[p4, p5, p6], p2, p3])) + 
     (MP[p4, p4]*Spbb[p3, SumM[p4, p5], p1, p6]*
       ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*
         Spbb[p3, p6] + Spbb[p3, SumM[p4, p5], p2, p6]))/
      (2*MP[p1, p6]*(2*MP[p4, p4] + 2*MP[p4, p5])*(2*MP[p3, p4] + 
        2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*(2*MP[p3, p4] + 
        ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5])*
          Spbb[p3, p2, p4, p3])/Spbb[p3, SumM[p4, p5], p2, p3])*
       (2*MP[p1, p1] + 2*MP[p1, p2] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 
           2*MP[p4, p4] + 2*MP[p4, p5])*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
         Spbb[p3, SumM[p4, p5], p2, p3])) + 
     (MP[p4, p4]*Spbb[p3, SumM[p4, p5, p6], p1, p3]*(Spbb[p3, p6, p5, p6] - 
        Spbb[p3, SumM[p4, p5, p6], p5, p6]))/(2*(2*MP[p1, p1] + 2*MP[p1, p2])*
       MP[p5, p6]*(-2*MP[p1, p1] - 2*MP[p1, p2] + Spab[p3, SumM[p4, p5, p6], 
         p3])*Spab[p6, SumM[p4, p5], p3]*(2*MP[p3, p4] + 
        ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, p2, p4, p3])/
         Spbb[p3, SumM[p4, p5, p6], p2, p3])) - 
     (Spbb[p3, p2, p4, p3]*(-((-8*MP[p1, p6]*MP[p5, p6] - 
           2*MP[p5, p6]*(2*MP[p1, p1] + 2*MP[p1, p2] + 
             (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/Spbb[p3, p2, p4, p3]) + 
           (2*MP[p1, p1] + 2*MP[p1, p2] + (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
              Spbb[p3, p2, p4, p3])*(2*MP[p2, p5] + (2*MP[p3, p4]*Spbb[p3, 
                p2, p5, p3])/Spbb[p3, p2, p4, p3]) - 2*MP[p1, p6]*
            (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5] + 
             (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p2, p4, 
               p3]) + (2*MP[p1, p3] + 2*MP[p1, p4] + 2*MP[p3, p4] + 
             (2*MP[p3, p4]*Spbb[p3, SumM[p1, p4], p2, p3])/Spbb[p3, p2, p4, 
               p3])*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 
             2*MP[p4, p5] + (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/
              Spbb[p3, p2, p4, p3]))*Spbb[p6, p5, p1, p6]) + 
        2*MP[p1, p5]*((2*MP[p3, p4]*Spbb[p3, p6]*
            (2*MP[p5, p6]*Spbb[p3, p2, p1, p6] + 2*MP[p1, p6]*
              Spbb[p3, p2, p5, p6]))/Spbb[p3, p2, p4, p3] - 
          2*MP[p5, p6]*Spbb[p6, p2, p1, p6] + 2*MP[p1, p6]*
           Spbb[p6, SumM[p3, p4], p5, p6])))/(16*MP[p1, p6]*MP[p3, p4]*
       MP[p5, p6]*Spab[p3, p2, p3]*(2*MP[p1, p1] + 2*MP[p1, p2] + 
        (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/Spbb[p3, p2, p4, p3])*
       (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p4] + 2*MP[p4, p5] + 
        (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p2, p4, 
          p3])) - 
     ((-((MP[p4, p4]*Spbb[p3, p6]^4*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
          ((2*MP[p4, p4] + 2*MP[p4, p5])*Spbb[p3, p2, SumM[p3, p4, p5], 
            p6])) - (Spbb[p3, SumM[p5, p6], p4, p3]^2*
          Spbb[p3, SumM[p4, p5, p6], p5, p6]^2)/(2*MP[p5, p6]*
          (-2*MP[p1, p1] - 2*MP[p1, p2] + Spab[p3, SumM[p4, p5, p6], p3])*
          Spab[p6, SumM[p4, p5], p3]*(2*MP[p3, p4] + 
           ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, p2, p4, p3])/
            Spbb[p3, SumM[p4, p5, p6], p2, p3])))*Spbb[p3, p2, SumM[p1, p2], 
        p1, p2, p3])/((2*MP[p1, p1] + 2*MP[p1, p2])*
       Spbb[p3, SumM[p4, p5, p6], p2, p3]^2*(Spbb[p3, p6, p5, p6] - 
        Spbb[p3, SumM[p4, p5, p6], p5, p6]))
 
ATreeSAM["S", "+", "+", "S", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (((-4*MP[p2, p4]*MP[p3, p4])/Spba[p2, p4, p3] + Spba[p3, p4, p2])*
       (1 + (2*MP[p1, p6] + ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
            Spbb[p2, p6, p1, p2])/Spbb[p2, SumM[p3, p4], p1, p2])/
         (2*MP[p5, p5] + 2*MP[p5, p6]))*Spbb[p2, SumM[p3, p4], p1, p2])/
      (2*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spaa[p2, p3]*
       (2*MP[p1, p2] + (2*MP[p3, p4]*Spba[p2, p1, p3])/Spba[p2, p4, p3])) + 
     (Spba[p2, p1, p3]*((-4*MP[p1, p2]^2*Spbb[p2, p4, SumM[p5, p6], p2])/
         Spba[p2, p1, p3]^2 - (2*MP[p1, p2]*(Spbb[p2, p4, SumM[p5, p6], p3] + 
           Spbb[p3, p4, SumM[p5, p6], p2]))/Spba[p2, p1, p3] - 
        Spbb[p3, p4, SumM[p5, p6], p3])*
       (1 + (2*MP[p1, p6] + 2*MP[p2, p6] - (2*MP[p1, p2]*Spba[p2, p6, p3])/
           Spba[p2, p1, p3] + ((2*MP[p3, p4] + (2*MP[p1, p2]*Spba[p2, p4, 
                p3])/Spba[p2, p1, p3])*((-4*MP[p1, p2]^2*Spbb[p2, p6, 
                SumM[p4, p5], p2])/Spba[p2, p1, p3]^2 - 
             (2*MP[p1, p2]*(Spbb[p2, p6, SumM[p4, p5], p3] + Spbb[p3, p6, 
                 SumM[p4, p5], p2]))/Spba[p2, p1, p3] - Spbb[p3, p6, 
              SumM[p4, p5], p3]))/((-4*MP[p1, p2]^2*Spbb[p2, p4, SumM[p5, 
                p6], p2])/Spba[p2, p1, p3]^2 - (2*MP[p1, p2]*(Spbb[p2, p4, 
                SumM[p5, p6], p3] + Spbb[p3, p4, SumM[p5, p6], p2]))/
             Spba[p2, p1, p3] - Spbb[p3, p4, SumM[p5, p6], p3]))/
         (2*MP[p5, p5] + 2*MP[p5, p6])))/(2*MP[p1, p2]*
       (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spaa[p2, p3]*
       (2*MP[p3, p4] + (2*MP[p1, p2]*Spba[p2, p4, p3])/Spba[p2, p1, p3]))
 
ATreeSAM["S", "+", "S", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (MP[p3, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p3, p4], p2]*
       (1 + (2*MP[p3, p5] + 2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p2, p5, p3, 
             p2])/Spbb[p2, p4, p3, p2] + 
          ((2*MP[p1, p2] + (2*MP[p3, p4]*Spbb[p2, p1, p3, p2])/
              Spbb[p2, p4, p3, p2])*Spbb[p2, p5, SumM[p3, p4], p2])/
           Spbb[p2, p1, SumM[p3, p4], p2])/(2*MP[p5, p5] + 2*MP[p5, p6])))/
      (2*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
       Spab[p4, p3, p2]*(2*MP[p1, p2] + (2*MP[p3, p4]*Spbb[p2, p1, p3, p2])/
         Spbb[p2, p4, p3, p2])) + (Spbb[p2, p1, p3, p2]*
       (-1/2*((Spbb[p4, p3, p5, p4] - (2*MP[p1, p2]*Spbb[p2, p4]*
              Spbb[p4, p5, p3, p2])/Spbb[p2, p1, p3, p2])*
           (1 + ((2*(MP[p3, p4] + MP[p3, p5] + MP[p4, p5]) + (2*MP[p1, p2]*
                 Spbb[p2, SumM[p4, p5], p3, p2])/Spbb[p2, p1, p3, p2])*
              (Spbb[p4, p3, p5, p4] - (2*MP[p1, p2]*Spbb[p2, p4]*Spbb[p4, p5, 
                  p3, p2])/Spbb[p2, p1, p3, p2]))/((2*MP[p5, p5] + 
                2*MP[p5, p6])*(Spbb[p4, p3, p5, p4] - (2*MP[p1, p2]*
                  Spbb[p2, p4]*Spbb[p4, p5, p3, p2])/Spbb[p2, p1, p3, p2]) + 
              (2*MP[p3, p4] + (2*MP[p1, p2]*Spbb[p2, p4, p3, p2])/
                 Spbb[p2, p1, p3, p2])*Spbb[p4, p6, p5, p4])))/
          (MP[p4, p5]*(2*MP[p3, p4] + (2*MP[p1, p2]*Spbb[p2, p4, p3, p2])/
             Spbb[p2, p1, p3, p2])) + (MP[p3, p3]*Spbb[p4, p5, p6, p4]^2 + 
          MP[p6, p6]*(Spbb[p4, p3, SumM[p1, p2], p4] - 
             (2*MP[p1, p2]*Spbb[p2, p4]*Spbb[p4, SumM[p1, p2, p3], p3, p2])/
              Spbb[p2, p1, p3, p2])^2)/((2*MP[p1, p1] + 2*MP[p1, p2] + 
           2*MP[p1, p3] + 2*MP[p2, p3])*(2*MP[p5, p5] + 2*MP[p5, p6])*
          ((2*MP[p5, p5] + 2*MP[p5, p6])*(Spbb[p4, p3, p5, p4] - 
             (2*MP[p1, p2]*Spbb[p2, p4]*Spbb[p4, p5, p3, p2])/
              Spbb[p2, p1, p3, p2]) + (2*MP[p3, p4] + 
             (2*MP[p1, p2]*Spbb[p2, p4, p3, p2])/Spbb[p2, p1, p3, p2])*
            Spbb[p4, p6, p5, p4]))))/(4*MP[p1, p2]*MP[p2, p3])
 
ATreeSAM["S", "+", "S", "s", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (Spbb[p2, p3, p1, p2]*(1 + (2*MP[p1, p6] + 
         (2*MP[p2, p3]*(Spbb[p2, p6, p1, p2] + (2*MP[p4, p5]*Spbb[p2, p5]*
              Spbb[p2, p1, p6, p5])/Spbb[p5, p4, p6, p5]))/
          Spbb[p2, p3, p1, p2] + (2*MP[p4, p5]*Spbb[p5, p1, p6, p5])/
          Spbb[p5, p4, p6, p5])/(2*MP[p4, p4] + 2*MP[p4, p5] + 2*MP[p4, p6] + 
         2*MP[p5, p6]))*Spbb[p5, p4, p6, p5])/(16*MP[p1, p2]*MP[p2, p3]*
      MP[p4, p5]*MP[p5, p6])
 
ATreeSAM["S", "+", "+", "s", "S", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (MP[p4, p4]*Spbb[p2, p3]^2)/(2*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 
        2*MP[p3, p4])*(2*MP[p2, p3] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4])*Spbb[p2, p3, p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
          p2])) - Spbb[p2, p1, SumM[p2, p3], p4, SumM[p5, p6], SumM[p2, p3], 
       p1, p2]/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
       Spaa[p3, p2, p1, p3]*(2*MP[p2, p4] + 2*MP[p3, p4] - 
        (2*MP[p2, p3]*Spbb[p2, p4, p1, p2])/Spbb[p2, p3, p1, p2]))
 
ATreeSAM["S", "+", "s", "+", "S", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    ((Spbb[p2, p3, p1, p2] + (2*MP[p4, p5]*Spbb[p2, p4]*Spbb[p4, p3, p1, p2])/
         Spbb[p4, p3, p5, p4])*Spbb[p4, p5, p3, p4])/(8*MP[p1, p2]*MP[p3, p4]*
       MP[p4, p5]*(2*MP[p2, p3] + (2*MP[p4, p5]*Spbb[p4, p2, p3, p4])/
         Spbb[p4, p5, p3, p4])) + (MP[p3, p3]*Spbb[p2, p4]*
       Spbb[p4, p5, SumM[p2, p3], p4])/(2*MP[p2, p3]*
       (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spab[p2, p3, p4]*
       (2*MP[p4, p5] + (2*MP[p2, p3]*Spbb[p4, p5, p3, p4])/
         Spbb[p4, p2, p3, p4]))
 
ATreeSAM["S", "+", "s", "S", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_] = 
    (Spbb[p2, p3, p1, p2]*Spbb[p5, p6, p4, p5])/(16*MP[p1, p2]*MP[p2, p3]*
      MP[p4, p5]*MP[p5, p6])
 
ATreeSAM["S", "S", "+", "+", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p6, p6]*(2*MP[p2, p7] + 
        ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
          Spbb[p3, SumM[p2, p7], p2, p3])/Spbb[p3, SumM[p4, p5, p6], p2, p3])*
       ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
         Spbb[p3, p5] - Spbb[p3, SumM[p4, p5, p6], p3, p5] - 
        Spbb[p3, SumM[p4, p5, p6], p4, p5]))/
      (2*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
       MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p4, p5]*
       (Spaa[p3, p4] + ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
           2*MP[p2, p7])*Spab[p4, p2, p3])/Spbb[p3, SumM[p4, p5, p6], p2, 
          p3])*(2*MP[p1, p1] + 2*MP[p1, p2] - 
        ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
          Spbb[p3, p2, p1, p3])/Spbb[p3, SumM[p4, p5, p6], p2, p3])) - 
     (MP[p6, p6]*Spbb[p3, SumM[p1, p2], p1, p3]*Spbb[p3, SumM[p1, p2], p2, 
        p3]*(Spbb[p3, p2, SumM[p1, p2], p7, p3, SumM[p4, p5], p6, p5] + 
        Spbb[p3, p2, SumM[p1, p2], p7, SumM[p1, p2], SumM[p4, p5], p6, p5]))/
      (2*(2*MP[p1, p1] + 2*MP[p1, p2])*MP[p5, p6]*(2*MP[p4, p5] + 
        2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p4, p5]*(2*MP[p1, p1] + 
        2*MP[p1, p2] + Spab[p3, SumM[p1, p2], p3])*
       ((2*MP[p1, p1] + 2*MP[p1, p2])*Spab[p4, p2, p3] - 
        Spaa[p3, p4]*Spbb[p3, SumM[p1, p2], p2, p3])*
       ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
         Spbb[p3, SumM[p1, p2], p2, p3] + (2*MP[p1, p1] + 2*MP[p1, p2])*
         Spbb[p3, SumM[p4, p5, p6], p2, p3])) - 
     (Spbb[p3, p2, SumM[p1, p2], p1, p2, p3]*
       (-((MP[p6, p6]*((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, p5] + 
             Spbb[p3, SumM[p1, p2], p3, p5] + Spbb[p3, SumM[p1, p2], p4, p5])^
            3)/((2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p5, p6] + 2*MP[p5, p7] + 
            2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p1] + 2*MP[p1, p2] + 
            Spab[p3, SumM[p1, p2], p3])*(Spab[p4, p3, SumM[p6, p7], p6, p5] + 
            Spab[p4, SumM[p1, p2], SumM[p6, p7], p6, p5])*
           ((2*MP[p1, p1] + 2*MP[p1, p2])*Spab[p4, p2, p3] - 
            Spaa[p3, p4]*Spbb[p3, SumM[p1, p2], p2, p3]))) - 
        (MP[p6, p6]*Spbb[p3, SumM[p4, p5], p6, p5]^3)/(2*MP[p5, p6]*
          Spaa[p4, p5]*(Spab[p4, p3, SumM[p6, p7], p6, p5] + 
           Spab[p4, SumM[p1, p2], SumM[p6, p7], p6, p5])*
          (Spbb[p3, p2, SumM[p1, p2], p7, p3, SumM[p4, p5], p6, p5] + 
           Spbb[p3, p2, SumM[p1, p2], p7, SumM[p1, p2], SumM[p4, p5], p6, 
            p5])) - ((2*MP[p1, p1] + 2*MP[p1, p2])*
            (-(MP[p6, p6]*Spbb[p3, p2, p4, p5]) + (2*MP[p4, p5] + 2*
                MP[p4, p6] + 2*MP[p5, p6])*Spbb[p3, p2, p6, p5])*
            Spbb[p3, SumM[p4, p5, p6], p7, p3] + Spbb[p3, SumM[p1, p2], p2, 
             p3]*(-(MP[p6, p6]*Spbb[p3, SumM[p1, p2], p7, p3, p4, p5]) + 
             (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*(Spbb[p3, 
                SumM[p1, p2], p7, p3, p6, p5] + Spbb[p3, SumM[p1, p2], p7, 
                SumM[p4, p5], p6, p5])))^2/(2*MP[p5, p6]*
          (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*Spaa[p4, p5]*
          (2*MP[p1, p1] + 2*MP[p1, p2] + Spab[p3, SumM[p1, p2], p3])*
          ((2*MP[p1, p1] + 2*MP[p1, p2])*Spab[p4, p2, p3] - 
           Spaa[p3, p4]*Spbb[p3, SumM[p1, p2], p2, p3])*
          ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
            Spbb[p3, SumM[p1, p2], p2, p3] + (2*MP[p1, p1] + 2*MP[p1, p2])*
            Spbb[p3, SumM[p4, p5, p6], p2, p3])*(Spbb[p3, p2, SumM[p1, p2], 
            p7, p3, SumM[p4, p5], p6, p5] + Spbb[p3, p2, SumM[p1, p2], p7, 
            SumM[p1, p2], SumM[p4, p5], p6, p5]))))/
      ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p3, SumM[p1, p2], p2, p3]) + 
     (Spbb[p3, p4]^3*Spbb[p3, p4, p2, p3]*
       ((MP[p6, p6]*Spbb[p3, p2, SumM[p3, p4], p5]^2*(2*MP[p2, p7] - 
           (2*MP[p3, p4]*Spbb[p3, p2, p7, p3])/Spbb[p3, p4, p2, p3] + 
           ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
              (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/Spbb[p3, p4, 
                p2, p3])*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
              SumM[p1, p5, p6], SumM[p3, p4], p2, p3])/Spbb[p3, p2, 
             SumM[p3, p4], SumM[p5, p6], SumM[p2, p3, p4], SumM[p3, p4], p2, 
             p3]))/(2*MP[p5, p6]*Spbb[p3, p4, p2, p3]^2*(2*MP[p1, p1] + 
           2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
           (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
            Spbb[p3, p4, p2, p3])*(2*MP[p3, p4] + 2*MP[p3, p5] + 
           2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/
            Spbb[p3, p4, p2, p3] - ((2*MP[p1, p1] + 2*MP[p1, p2] + 
              2*MP[p1, p7] + 2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, 
                 SumM[p4, p5, p6], p2, p3])/Spbb[p3, p4, p2, p3])*
             Spbb[p3, p2, SumM[p3, p4], p5, SumM[p2, p3, p4], SumM[p3, p4], 
              p2, p3])/Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6], 
             SumM[p2, p3, p4], SumM[p3, p4], p2, p3])*(2*MP[p1, p1] + 
           2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
            Spbb[p3, p4, p2, p3] - ((2*MP[p1, p1] + 2*MP[p1, p2] + 
              2*MP[p1, p7] + 2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, 
                 SumM[p4, p5, p6], p2, p3])/Spbb[p3, p4, p2, p3])*
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, SumM[p3, p4], 
              p2, p3])/Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6], 
             SumM[p2, p3, p4], SumM[p3, p4], p2, p3])) + 
        (MP[p6, p6]*(((2*MP[p3, p4]*Spbb[p3, p7, SumM[p3, p4], p5]*Spbb[p3, 
                p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
              Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p5]*Spbb[p3, p2, 
                SumM[p3, p4], SumM[p2, p3, p4], p1, p7, p2, p3])/
              Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p1, p7, SumM[p3, p4], p5] - 
             (2*MP[p3, p4]*Spbb[p3, p5]*((-2*MP[p3, p4]*Spbb[p3, p2, p7, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
                 Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p2, p7, p2, p3]))/Spbb[p3, p4, p2, p3] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p7, 
              SumM[p3, p4], p5])/Spbb[p3, p4, p2, p3] - 
           ((2*MP[p3, p4]*Spbb[p3, p7, SumM[p3, p4], p5]*Spbb[p3, p2, 
                SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, 
               p3] + (2*MP[p3, p4]*Spbb[p3, p7, SumM[p6, p7], p5]*Spbb[p3, 
                p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
              Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p5]*Spbb[p3, p2, 
                SumM[p3, p4], SumM[p2, p3, p4], p1, p7, p2, p3])/
              Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p1, p7, SumM[p3, p4], p5] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, p7, 
              SumM[p6, p7], p5] - (2*MP[p3, p4]*Spbb[p3, p5]*(
                (-2*MP[p3, p4]*Spbb[p3, p2, p7, p3]*Spbb[p3, p2, SumM[p3, 
                    p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p7, p2, 
                 p3]))/Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p2, p7, SumM[p3, p4], p5] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p7, 
              SumM[p6, p7], p5])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, 
           SumM[p3, p4], SumM[p2, p3, p4], SumM[p5, p6, p7], SumM[p3, p4], 
           p2, p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6, p7], p1, 
           SumM[p3, p4], p2, p3])/(2*MP[p5, p6]*(2*MP[p1, p1] + 
           2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
            Spbb[p3, p4, p2, p3])*Spbb[p3, p4, p2, p3]^4*
          (-2*MP[p1, p1] - 2*MP[p1, p2] + (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
            Spbb[p3, p4, p2, p3] + Spbb[p3, p4, SumM[p5, p6, p7], 
             SumM[p3, p4], p2, p3]/Spbb[p3, p4, p2, p3])*
          ((Spab[p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
             (2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, 
                 p3])/Spbb[p3, p4, p2, p3]))/Spbb[p3, p4, p2, p3] - 
           (Spab[p5, p4, p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
              SumM[p5, p6, p7], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3]^
             2)*(((2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, 
                 p1, p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
              SumM[p5, p6], SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/
            Spbb[p3, p4, p2, p3]^2 + ((2*MP[p1, p1] + 2*MP[p1, p2] + 
              2*MP[p1, p7] + 2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, 
                 SumM[p4, p5, p6], p2, p3])/Spbb[p3, p4, p2, p3])*
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p5, p6, p7], 
              SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3]^2)) + 
        (((MP[p6, p6]*Spbb[p3, p2, SumM[p3, p4], p5]^3)/
            ((2*MP[p6, p6] + 2*MP[p6, p7])*Spbb[p3, p4, p2, p3]^3) + 
           (((2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, 
                    p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p6, p5]*Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p5, p6], p7, SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                 p3]^3 - (Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                 SumM[p5, p6, p7], SumM[p3, p4], p2, p3]*
                (Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6, p7], p7, p5, p6, p5]/
                  Spbb[p3, p4, p2, p3] + ((2*MP[p3, p4]*Spbb[p3, p2, p6, p5]*
                     Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6, p7], p7, p3])/
                    Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                    SumM[p5, p6, p7], p7, SumM[p3, p4], p6, p5])/Spbb[p3, p4, 
                   p2, p3]))/Spbb[p3, p4, p2, p3]^2)^2/(2*MP[p5, p6]*
             (-2*MP[p1, p1] - 2*MP[p1, p2] + (2*MP[p3, p4]*Spbb[p3, p2, p1, 
                 p3])/Spbb[p3, p4, p2, p3] + Spbb[p3, p4, SumM[p5, p6, p7], 
                SumM[p3, p4], p2, p3]/Spbb[p3, p4, p2, p3])*
             ((Spab[p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
                (2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, 
                    p3])/Spbb[p3, p4, p2, p3]))/Spbb[p3, p4, p2, p3] - 
              (Spab[p5, p4, p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                 SumM[p5, p6, p7], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                 p3]^2)*(((2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*
                   Spbb[p3, p2, p1, p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p5, p6], SumM[p2, p3, p4], SumM[p3, p4], 
                 p2, p3])/Spbb[p3, p4, p2, p3]^2 + ((2*MP[p1, p1] + 
                 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
                 (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
                  Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], SumM[p5, p6, p7], SumM[p3, p4], p2, p3])/
               Spbb[p3, p4, p2, p3]^2)))*
          ((-2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, 
              p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
            Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
            SumM[p2, p3, p4], p2, p1, SumM[p2, p3, p4], SumM[p3, p4], p2, 
            p3]))/((2*MP[p1, p1] + 2*MP[p1, p2] - 
           (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/Spbb[p3, p4, p2, p3])*
          (((2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, 
                 p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p6, p5])/Spbb[p3, p4, p2, p3] + 
           ((2*MP[p3, p4]*Spbb[p3, p5, p6, p5]*Spbb[p3, p2, SumM[p3, p4], 
                SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, p5, p6, p5] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p5, p6, p5])/
            Spbb[p3, p4, p2, p3] + ((2*MP[p3, p4]*Spbb[p3, p2, p6, p5]*Spbb[
                p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p1, p3])/
              Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*Spbb[p3, p2, p6, p5]*Spbb[
                p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
              Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*Spbb[p3, p4, p6, p5]*Spbb[
                p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/
              Spbb[p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p1, SumM[p3, p4], p6, p5] + 
             Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, SumM[p3, p4], 
              p6, p5])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
           SumM[p2, p3, p4], SumM[p5, p6, p7], SumM[p3, p4], p2, p3]) + 
        (Spbb[p3, p2, SumM[p3, p4], p5]^3*Spbb[p3, p4, p2, p3]^4*
          (((((-((2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p2, p7] - (2*MP[p3, p4]*
                      Spbb[p3, p2, p7, p3])/Spbb[p3, p4, p2, p3])) + 
                 (2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, 
                      p1, p3])/Spbb[p3, p4, p2, p3])*(2*MP[p1, p1] + 
                   2*MP[p1, p2] - 2*MP[p1, p6] + 2*MP[p1, p7] + 
                   2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], 
                      p2, p3])/Spbb[p3, p4, p2, p3]) + (2*MP[p1, p6] + 
                   2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
                  (2*MP[p6, p6] + 2*MP[p6, p7] + 2*(2*MP[p1, p1] + 
                     2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
                     (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
                      Spbb[p3, p4, p2, p3])))*Spbb[p3, p2, SumM[p3, p4], p5, 
                  SumM[p2, p3, p4], SumM[p3, p4], p2, p3]^2)/Spbb[p3, p4, p2, 
                 p3]^4 + ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] - 
                  (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p4, 
                    p2, p3])^2*Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6], 
                 SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p1, SumM[p3, p4], p2, p3])/
               Spbb[p3, p4, p2, p3]^4 + ((2*MP[p3, p4] + 2*MP[p3, p5] + 
                 2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/
                  Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], p5, 
                 SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
                (-(((2*MP[p1, p1] + 2*MP[p1, p2] + 2*(2*MP[p1, p6] + 
                       2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7]) - 
                     (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/Spbb[p3, p4, p2, 
                       p3])*Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6], 
                     SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/Spbb[p3, p4, 
                     p2, p3]^2) + ((-2*MP[p1, p1] - 2*MP[p1, p2] + 
                    2*MP[p1, p6] - 2*MP[p1, p7] - 2*MP[p2, p7] + 
                    (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
                     Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
                    SumM[p2, p3, p4], p1, SumM[p3, p4], p2, p3])/
                  Spbb[p3, p4, p2, p3]^2 + ((2*MP[p6, p6] + 2*MP[p6, p7])*
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p7, 
                    SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3]^2))/
               Spbb[p3, p4, p2, p3]^2)*((-2*MP[p3, p4]*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                ((2*MP[p3, p4]*Spbb[p3, p2, p6, p3]*Spbb[p3, p5, SumM[p2, p3, 
                     p4], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                 Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p2, p6, 
                  p3]))/Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p5, p6, p2, p3])/Spbb[p3, 
                p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p4, SumM[p2, p3, p4], 
                 SumM[p3, p4], p2, p3]*Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p6, p2, p3])/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, 
                 p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, 
                 SumM[p3, p4], p6, p3])/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, 
                 p3]*((-2*MP[p3, p4]*Spbb[p3, p2, p6, p3]*Spbb[p3, p2, 
                    SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                   p3] + (2*MP[p3, p4]*Spbb[p3, p2, p6, p3]*Spbb[p3, p4, 
                    SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                   p3] + Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                  SumM[p3, p4], p2, p6, p3]))/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, 
                 p3]*((-2*MP[p3, p4]*Spbb[p3, p5, p6, p3]*Spbb[p3, p2, 
                    SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                   p3] + Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                  SumM[p3, p4], p5, p6, p3]))/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
                 p3]*((-2*MP[p3, p4]*Spbb[p3, p2, p6, p3]*Spbb[p3, p2, 
                    SumM[p3, p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, 
                   p3] - Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                  SumM[p3, p4], p6, p2, p3]))/Spbb[p3, p4, p2, p3] + 
              (2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, 
                 p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                 SumM[p3, p4], p6, p2, p3])/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p4, SumM[p2, p3, p4], SumM[p3, p4], p2, 
                 p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                 SumM[p3, p4], p6, p2, p3])/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
                 p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                 SumM[p3, p4], p6, p2, p3])/Spbb[p3, p4, p2, p3] - 
              (2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, 
                 p3]*Spbb[p3, p5, p6, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
                 p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*Spbb[p3, p5, p6, 
                 SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/Spbb[
                p3, p4, p2, p3] + Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, 
                p4], p5, p2, p6, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
               p3] + Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p2, p6, 
               SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] - 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p6, p2, SumM[
                p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] - 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p6, p5, SumM[
                p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, SumM[p3, p4], 
               p6, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, SumM[p3, p4], 
               p6, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, 
               p6, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
              Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p5, 
               p6, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3]))/
            Spbb[p3, p4, p2, p3]^2 + 
           (((2*MP[p2, p6] - (2*MP[p3, p4]*Spbb[p3, p2, p6, p3])/Spbb[p3, p4, 
                  p2, p3])*Spbb[p3, p2, SumM[p3, p4], p5, SumM[p2, p3, p4], 
                SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3]^2 - 
             ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] - 
                (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p4, 
                  p2, p3])*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p6, 
                SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3]^2)*
            (((-(((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
                    2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], 
                       p2, p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, 
                     p4], p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/
                  Spbb[p3, p4, p2, p3]^2) + ((2*MP[p3, p4] + 2*MP[p3, p5] + 
                   2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, 
                      p3])/Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
                   SumM[p5, p6], SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/
                 Spbb[p3, p4, p2, p3]^2)*((2*MP[p3, p4]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p5, p1, p2, p3])/
                 Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p4, 
                   SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p5, p1, p2, p3])/
                 Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  ((2*MP[p3, p4]*Spbb[p3, p2, p1, p3]*Spbb[p3, p5, SumM[p2, 
                       p3, p4], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p2, p1, 
                    p3]))/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p5, SumM[p3, p4], p1, p3])/
                 Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p5, 
                   SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
                  ((-2*MP[p3, p4]*Spbb[p3, p2, p1, p3]*Spbb[p3, p2, SumM[p3, 
                       p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] - 
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                    p1, p2, p3]))/Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                   p1, p2, p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*
                  Spbb[p3, p4, SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                   p1, p2, p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*
                  Spbb[p3, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                   p1, p2, p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  ((-2*MP[p3, p4]*Spbb[p3, p2, p1, p3]*Spbb[p3, p2, SumM[p3, 
                       p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                   (2*MP[p3, p4]*Spbb[p3, p2, p1, p3]*Spbb[p3, p4, SumM[p2, 
                       p3, p4], SumM[p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                    p2, p1, p3]))/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  ((-2*MP[p3, p4]*Spbb[p3, p5, p1, p3]*Spbb[p3, p2, SumM[p3, 
                       p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                    p5, p1, p3]))/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  Spbb[p3, p5, p1, p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
                   p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*Spbb[p3, p5, p1, 
                   SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/
                 Spbb[p3, p4, p2, p3] - Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p1, p2, SumM[p3, p4], SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] - Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p1, p5, SumM[p3, p4], SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p2, p1, p5, SumM[p2, p3, p4], 
                 SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p2, p1, SumM[p3, p4], SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, SumM[p3, p4], p1, p5, SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, SumM[p3, p4], p1, SumM[p3, p4], 
                 SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p1, 
                 SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
                Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                 p5, p1, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, 
                 p3]))/Spbb[p3, p4, p2, p3]^2 + ((2*MP[p1, p6] + 
                2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*Spbb[p3, p2, 
                SumM[p3, p4], p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3]*(
                (-2*MP[p3, p4]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], 
                   p2, p3]*Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, 
                   p6, p7, p3])/Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, p7, p6, 
                   p3])/Spbb[p3, p4, p2, p3] - (2*MP[p3, p4]*Spbb[p3, p2, 
                   SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  ((-2*MP[p3, p4]*Spbb[p3, p6, p7, p3]*Spbb[p3, p2, SumM[p3, 
                       p4], SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
                   Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                    p6, p7, p3]))/Spbb[p3, p4, p2, p3] + (2*MP[p3, p4]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p2, p3]*
                  Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], 
                   p7, p6, p3])/Spbb[p3, p4, p2, p3] + Spbb[p3, p2, 
                 SumM[p3, p4], SumM[p2, p3, p4], p5, p6, p7, p5, SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p6, p7, SumM[p3, p4], SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] - Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], p5, p7, p6, SumM[p3, p4], SumM[p2, p3, 
                  p4], SumM[p3, p4], p2, p3] + Spbb[p3, p2, SumM[p3, p4], 
                 SumM[p2, p3, p4], SumM[p3, p4], p6, p7, SumM[p3, p4], 
                 SumM[p2, p3, p4], SumM[p3, p4], p2, p3]))/Spbb[p3, p4, p2, 
                p3]^4)))/(2*(2*MP[p6, p6] + 2*MP[p6, p7])*
          (2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
          (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] - 
           (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p4, p2, 
             p3])*((-2*MP[p3, p4]*Spbb[p3, p5]*Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p2, p3])/Spbb[p3, p4, p2, p3] + 
           Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p5])*
          Spbb[p3, p2, SumM[p3, p4], p5, SumM[p2, p3, p4], SumM[p3, p4], p2, 
            p3]^3*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
           2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5, p6], p2, p3])/
            Spbb[p3, p4, p2, p3] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 
              2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/
               Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], SumM[p5, p6], 
              SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/Spbb[p3, p2, 
             SumM[p3, p4], p5, SumM[p2, p3, p4], SumM[p3, p4], p2, p3])*
          (2*MP[p1, p1] + 2*MP[p1, p2] - (2*MP[p3, p4]*Spbb[p3, p2, p1, p3])/
            Spbb[p3, p4, p2, p3] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 
              2*MP[p4, p5] - (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/
               Spbb[p3, p4, p2, p3])*Spbb[p3, p2, SumM[p3, p4], 
              SumM[p2, p3, p4], p1, SumM[p3, p4], p2, p3])/
            Spbb[p3, p2, SumM[p3, p4], p5, SumM[p2, p3, p4], SumM[p3, p4], 
             p2, p3]))))/(2*MP[p3, p4]*Spbb[p3, p2, p3, p4])
 
ATreeSAM["S", "S", "+", "+", "s", "s", "+"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p1, p1]*(Spbb[p7, p3, p2, p3] - Spbb[p7, SumM[p1, p2, p3], 
         p2, p3])*Spbb[p7, SumM[p1, p2, p3], p5, p4]*
       (-((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
          Spbb[p4, p7]) + Spbb[p7, SumM[p1, p2, p3], p6, p4]))/
      (4*MP[p2, p3]*MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 
        2*MP[p5, p6])*Spab[p3, SumM[p1, p2], p7]*(-2*MP[p4, p5] - 
        2*MP[p4, p6] - 2*MP[p5, p5] - 2*MP[p5, p6] + 
        Spab[p7, SumM[p1, p2, p3], p7])*(2*MP[p1, p7] + 
        ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
          Spbb[p7, p6, p1, p7])/Spbb[p7, SumM[p1, p2, p3], p6, p7])*
       (2*MP[p5, p5] + 2*MP[p5, p6] + ((2*MP[p4, p5] + 2*MP[p4, p6] + 
           2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p7, SumM[p5, p6], p6, p7])/
         Spbb[p7, SumM[p1, p2, p3], p6, p7])) + 
     (MP[p5, p5]*Spbb[p7, p6, SumM[p5, p6], p4]^2*
       (-1/2*(Spbb[p7, SumM[p2, p3], p1, p7]^2*Spbb[p7, SumM[p1, p2, p3], p2, 
             p3]^2)/(MP[p2, p3]*Spab[p3, SumM[p1, p2], p7]*
           (-2*MP[p4, p5] - 2*MP[p4, p6] - 2*MP[p5, p5] - 2*MP[p5, p6] + 
            Spab[p7, SumM[p1, p2, p3], p7])*(2*MP[p1, p7] + 
            ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
              Spbb[p7, p6, p1, p7])/Spbb[p7, SumM[p1, p2, p3], p6, p7])) - 
        (MP[p1, p1]*Spbb[p3, p7]^4*Spbb[p7, SumM[p1, p2, p3], p6, p7])/
         ((2*MP[p1, p1] + 2*MP[p1, p2])*Spbb[p7, p6, SumM[p1, p2, p7], p3])))/
      (2*MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 
        2*MP[p5, p6])*(Spbb[p7, p3, p2, p3] - Spbb[p7, SumM[p1, p2, p3], p2, 
         p3])*(2*MP[p5, p5] + 2*MP[p5, p6] + 
        ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
          Spbb[p7, SumM[p5, p6], p6, p7])/Spbb[p7, SumM[p1, p2, p3], p6, p7])*
       Spbb[p7, SumM[p1, p2, p3], p6, p7]^2) - 
     (MP[p1, p1]*Spbb[p7, SumM[p5, p6], p5, p7]*
       (Spbb[p7, p1, p7, SumM[p3, p4], p2, p3] + Spbb[p7, p1, SumM[p5, p6], 
         SumM[p3, p4], p2, p3]))/(2*MP[p2, p3]*(2*MP[p2, p3] + 2*MP[p2, p4] + 
        2*MP[p3, p4])*(2*MP[p5, p5] + 2*MP[p5, p6])*Spaa[p3, p4]*
       Spab[p4, SumM[p1, p2, p3], p7]*(2*MP[p5, p5] + 2*MP[p5, p6] + 
        Spab[p7, SumM[p5, p6], p7])*(2*MP[p1, p7] - 
        ((2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p7, p6, p1, p7])/
         Spbb[p7, SumM[p5, p6], p6, p7])) + 
     (MP[p5, p5]*Spbb[p7, p1, p2, p7]^2*
       ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
         Spbb[p4, p7]*Spbb[p7, p6, SumM[p1, p2, p7], SumM[p3, p4, p5], p6, 
          p7] + Spbb[p7, SumM[p1, p2], p6, p7]*
         (Spbb[p7, p6, SumM[p1, p2, p7], SumM[p3, p4, p5], p3, p4] + 
          Spbb[p7, p6, SumM[p1, p2, p7], SumM[p3, p4, p5], SumM[p1, p2, p7], 
           p4])))/(2*(2*MP[p1, p1] + 2*MP[p1, p2])*(2*MP[p1, p1] + 
        2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*MP[p4, p5]*
       (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*Spaa[p3, p4]*
       Spab[p3, SumM[p1, p2], p7]*Spbb[p7, SumM[p1, p2], p6, p7]*
       ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
         Spbb[p7, p6, p1, p7] + 2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, 
          p7])) - (MP[p1, p1]*Spbb[p7, SumM[p1, p2], p6, p7]*
       (-((MP[p5, p5]*Spbb[p3, p4]^3*Spbb[p7, SumM[p1, p2], p6, p7]^2)/
          (Spbb[p7, p6, SumM[p1, p2, p7], p3]*
           (-((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
              Spbb[p7, p6, p5, p7]) + (2*MP[p5, p5] + 2*MP[p5, p6])*
             Spbb[p7, SumM[p1, p2], p6, p7]))) + 
        ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
            Spbb[p7, SumM[p3, p4], p5, p4] + Spbb[p7, SumM[p1, p2], p6, 
            SumM[p3, p4], p5, p4])^2/(2*MP[p4, p5]*(2*MP[p3, p4] + 
           2*MP[p3, p5] + 2*MP[p4, p5])*Spaa[p3, p4]*Spab[p3, SumM[p1, p2], 
           p7])))/((2*MP[p1, p1] + 2*MP[p1, p2])*(2*MP[p1, p1] + 
        2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*(2*MP[p1, p7] + 
        ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
          Spbb[p7, p6, p1, p7])/Spbb[p7, SumM[p1, p2], p6, p7])*
       Spbb[p7, p6, SumM[p1, p2, p7], SumM[p3, p4], p5, p4]) - 
     (Spbb[p7, p6, SumM[p5, p6], p5, p6, p7]*
       (-((MP[p1, p1]*Spbb[p7, SumM[p5, p6], p6, p7]*
           (-((2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p3, p7]) + 
             Spbb[p7, SumM[p5, p6], p4, p3] + Spbb[p7, SumM[p5, p6], p7, p3])^
            4)/((2*MP[p1, p1] + 2*MP[p1, p2])*(2*MP[p1, p1] + 2*MP[p1, p2] + 
            2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p4, SumM[p1, p2, p3], p7]*
           (2*MP[p5, p5] + 2*MP[p5, p6] + Spab[p7, SumM[p5, p6], p7])*
           (Spab[p4, p7, SumM[p1, p2], p2, p3] + Spab[p4, SumM[p5, p6], 
             SumM[p1, p2], p2, p3])*((2*MP[p5, p5] + 2*MP[p5, p6])*
             (-(Spab[p7, p6, p7]*Spbb[p3, p7]) + Spbb[p7, p6, p4, p3] + 
              Spbb[p7, p6, SumM[p3, p5, p6], p3]) - Spab[p7, p4, p3]*
             Spbb[p7, SumM[p5, p6], p6, p7] + Spab[p7, SumM[p1, p2, p4], p3]*
             Spbb[p7, SumM[p5, p6], p6, p7]))) + 
        (MP[p1, p1]*Spbb[p3, p4]^3*(Spbb[p7, SumM[p3, p4], p2, p3] + 
           Spbb[p7, SumM[p5, p6], p2, p3])*Spbb[p7, SumM[p5, p6], p6, p7]^5)/
         (Spbb[p7, p6, SumM[p5, p6], p4]*((2*MP[p5, p5] + 2*MP[p5, p6])*
            Spbb[p7, p6, p1, p7] - 2*MP[p1, p7]*Spbb[p7, SumM[p5, p6], p6, 
             p7])*((2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p7, SumM[p1, p2], p6, 
             p7] + (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
             2*MP[p2, p7])*Spbb[p7, SumM[p5, p6], p6, p7])*
          ((2*MP[p5, p5] + 2*MP[p5, p6])*(-(Spab[p7, p6, p7]*Spbb[p3, p7]) + 
             Spbb[p7, p6, p4, p3] + Spbb[p7, p6, SumM[p3, p5, p6], p3]) - 
           Spab[p7, p4, p3]*Spbb[p7, SumM[p5, p6], p6, p7] + 
           Spab[p7, SumM[p1, p2, p4], p3]*Spbb[p7, SumM[p5, p6], p6, p7])*
          Spbb[p7, p6, SumM[p5, p6], SumM[p3, p4], p2, p3]) + 
        (Spbb[p3, p4]^3*(-((MP[p1, p1]*Spbb[p7, SumM[p3, p4], p2, p3]^4*
              Spbb[p7, SumM[p5, p6], p6, p7])/((2*MP[p1, p1] + 2*MP[p1, p2] + 
               (2*MP[p3, p4]*Spbb[p3, p1, p2, p3])/Spbb[p3, p4, p2, p3])*
              Spbb[p7, p6, SumM[p5, p6], SumM[p3, p4], p2, p3])) + 
           (Spbb[p3, p4, p2, p3]*Spbb[p7, SumM[p2, p3, p4], p1, p7]^2*
             Spbb[p7, SumM[p5, p6], SumM[p2, p3, p4], SumM[p3, p4], p2, p3]^
              2)/((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
             (2*MP[p5, p5] + 2*MP[p5, p6] + Spab[p7, SumM[p5, p6], p7])*
             Spbb[p7, SumM[p5, p6], p4, p3]*(2*MP[p1, p7] - 
              ((2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p7, p6, p1, p7])/Spbb[p7, 
                SumM[p5, p6], p6, p7]))))/(2*MP[p3, p4]*Spbb[p3, p2, p3, p4]*
          (2*MP[p3, p4]*Spbb[p3, p7]*Spbb[p3, p2, SumM[p3, p4], 
             SumM[p2, p3, p4], p2, p3] + Spbb[p3, p4, p2, p3]*
            Spbb[p7, SumM[p3, p4], SumM[p2, p3, p4], SumM[p3, p4], p2, p3] + 
           Spbb[p3, p4, p2, p3]*Spbb[p7, SumM[p5, p6], SumM[p2, p3, p4], 
             SumM[p3, p4], p2, p3]))))/((2*MP[p5, p5] + 2*MP[p5, p6])*
       Spbb[p7, SumM[p5, p6], p6, p7]^2) - 
     (Spbb[p7, p6, p1, p7]*((MP[p5, p5]*Spbb[p3, p4]^2*(2*MP[p2, p6] - 
           ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
             (Spbb[p3, SumM[p1, p7], p2, p3] + Spbb[p3, SumM[p4, p5], p2, 
               p3] - (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/Spbb[
                p7, p6, p1, p7]))/Spbb[p3, SumM[p4, p5], p2, p3] + 
           (2*MP[p1, p7]*Spbb[p7, p6, p2, p7])/Spbb[p7, p6, p1, p7]))/
         (2*MP[p4, p5]*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
          (2*MP[p3, p4] - ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
             Spbb[p3, p4, p2, p3])/Spbb[p3, SumM[p4, p5], p2, p3])*
          (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
           ((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
             (-Spbb[p3, SumM[p1, p7], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
                Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7]))/
            Spbb[p3, SumM[p4, p5], p2, p3] + 
           (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, 
             p7])) - (MP[p5, p5]*(Spbb[p3, p6, p2, p3] + 
           Spbb[p3, SumM[p4, p5], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
             Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7])*
          (-Spbb[p3, SumM[p1, p7], p6, p3] - Spbb[p3, SumM[p1, p7], 
            SumM[p4, p5], p3] - (2*MP[p1, p7]*MP[p5, p5]*Spbb[p3, p7]^2)/
            Spbb[p7, p6, p1, p7] + (2*MP[p1, p7]*Spbb[p3, p7]*
             Spbb[p7, p6, SumM[p1, p7], p3])/Spbb[p7, p6, p1, p7] + 
           (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, SumM[p4, p5], p3])/
            Spbb[p7, p6, p1, p7])*(-(MP[p5, p5]*(Spbb[p3, p2, SumM[p1, p7], 
               p4] + (2*MP[p1, p7]*Spbb[p4, p7]*Spbb[p7, p6, p2, p3])/Spbb[
                p7, p6, p1, p7])) - (2*MP[p1, p7]*Spbb[p7, p1, p2, p3]*
             Spbb[p7, p6, p3, p4])/Spbb[p7, p6, p1, p7] - 
           (2*MP[p1, p7]*Spbb[p7, p6, p2, p3]*Spbb[p7, p6, p3, p4])/
            Spbb[p7, p6, p1, p7] + MP[p1, p1]*(Spbb[p3, p6, p3, p4] + 
             (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p3, p4])/
              Spbb[p7, p6, p1, p7]) + (2*MP[p1, p7]*Spbb[p7, p1, p2, p3]*
             Spbb[p7, p6, SumM[p3, p5], p4])/Spbb[p7, p6, p1, p7] + 
           (2*MP[p1, p7]*Spbb[p7, p6, p2, p3]*Spbb[p7, p6, SumM[p3, p5], p4])/
            Spbb[p7, p6, p1, p7] - MP[p1, p1]*(MP[p5, p5]*Spbb[p3, p4] + 
             Spbb[p3, p6, SumM[p3, p5], p4] + (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[
                p7, p6, SumM[p3, p5], p4])/Spbb[p7, p6, p1, p7]) + 
           Spbb[p3, p2, SumM[p1, p7], p6, p3, p4] - Spbb[p3, p2, 
            SumM[p1, p7], p6, SumM[p3, p5], p4]))/(2*MP[p4, p5]*
          (-2*MP[p1, p1] - 2*MP[p1, p2] - 2*MP[p1, p7] - 2*MP[p2, p7] + 
           Spab[p3, p6, p3] + Spab[p3, SumM[p4, p5], p3] + 
           (2*MP[p1, p7]*Spab[p3, p6, p7]*Spbb[p3, p7])/Spbb[p7, p6, p1, 
             p7] - (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/
            Spbb[p7, p6, p1, p7])*(2*MP[p1, p1] + 2*MP[p1, p2] + 
           2*MP[p1, p7] + 2*MP[p2, p7] + (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], 
              p6, p7])/Spbb[p7, p6, p1, p7])*
          (Spaa[p3, p4]*(Spbb[p3, p6, p2, p3] + Spbb[p3, SumM[p4, p5], p2, 
              p3] + (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/
              Spbb[p7, p6, p1, p7]) + Spab[p4, p2, p3]*(2*MP[p1, p1] + 
             2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] + 
             (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, 
               p7]))*(-((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
             (Spbb[p3, p6, p2, p3] + Spbb[p3, SumM[p4, p5], p2, p3] + 
              (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/Spbb[p7, p6, 
                p1, p7])) + Spbb[p3, SumM[p4, p5], p2, p3]*
            (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] + 
             (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, 
               p7]))) - (MP[p1, p1]*(Spbb[p3, SumM[p1, p7], p2, p3] - 
           (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/
            Spbb[p7, p6, p1, p7])*((MP[p5, p5]*Spbb[p3, p4]^3)/
            (2*MP[p5, p5] + 2*MP[p5, p6] + (2*MP[p1, p7]*Spbb[p7, p6, p5, 
                p7])/Spbb[p7, p6, p1, p7]) + 
           (Spbb[p3, p2, p5, p4]*(Spbb[p3, SumM[p4, p5], p6, p3] - 
                (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, SumM[p4, p5], p3])/
                 Spbb[p7, p6, p1, p7])*(2*MP[p1, p1] + 2*MP[p1, p2] + 
                2*MP[p1, p7] + 2*MP[p2, p7] + (2*MP[p1, p7]*Spbb[p7, 
                   SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, p7]) + 
              (Spbb[p3, p6, p2, p3] + Spbb[p3, SumM[p4, p5], p2, p3] + 
                (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/Spbb[p7, p6, 
                  p1, p7])*(MP[p5, p5]*Spbb[p3, p4, p5, p4] + 
                (2*MP[p1, p7]*Spbb[p7, p3, p5, p4]*Spbb[p7, p6, SumM[p4, p5], 
                   p3])/Spbb[p7, p6, p1, p7] + (2*MP[p1, p7]*Spbb[p7, p4, p5, 
                   p4]*Spbb[p7, p6, SumM[p4, p5], p3])/Spbb[p7, p6, p1, p7] + 
                Spbb[p3, SumM[p4, p5], p6, p3, p5, p4] + Spbb[p3, 
                 SumM[p4, p5], p6, p4, p5, p4]))^2/(2*MP[p4, p5]*
             (-2*MP[p1, p1] - 2*MP[p1, p2] - 2*MP[p1, p7] - 2*MP[p2, p7] + 
              Spab[p3, p6, p3] + Spab[p3, SumM[p4, p5], p3] + 
              (2*MP[p1, p7]*Spab[p3, p6, p7]*Spbb[p3, p7])/Spbb[p7, p6, p1, 
                p7] - (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, 
                p6, p1, p7])*(Spaa[p3, p4]*(Spbb[p3, p6, p2, p3] + 
                Spbb[p3, SumM[p4, p5], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
                  Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7]) + 
              Spab[p4, p2, p3]*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
                2*MP[p2, p7] + (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/
                 Spbb[p7, p6, p1, p7]))*(-((2*MP[p3, p4] + 2*MP[p3, p5] + 
                 2*MP[p4, p5])*(Spbb[p3, p6, p2, p3] + Spbb[p3, SumM[p4, p5], 
                  p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/
                  Spbb[p7, p6, p1, p7])) + Spbb[p3, SumM[p4, p5], p2, p3]*(
                2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] + 
                (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, p6, 
                  p1, p7])))))/((Spbb[p3, p6, p2, p3] + 
           Spbb[p3, SumM[p4, p5], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
             Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7])*
          (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] + 
           (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, 
             p7])*(MP[p1, p1]*Spbb[p3, p4, p5, p4] - 
           (2*MP[p1, p7]*Spbb[p7, p3, p5, p4]*Spbb[p7, p6, p2, p3])/
            Spbb[p7, p6, p1, p7] - (2*MP[p1, p7]*Spbb[p7, p4, p5, p4]*
             Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7] + 
           Spbb[p3, p2, p5, p4]*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 
             2*MP[p2, p7] + (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, p7])/
              Spbb[p7, p6, p1, p7]) + Spbb[p3, p2, SumM[p1, p7], p3, p5, 
            p4] + Spbb[p3, p2, SumM[p1, p7], p4, p5, p4])) + 
        (Spbb[p3, p4]^3*((-2*MP[p3, p4]*Spbb[p3, p2, p5, p3] + 
             2*MP[p2, p5]*Spbb[p3, p4, p2, p3])*
            ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p3, p4, p2, 
               p3]*((2*MP[p1, p7]*Spbb[p7, SumM[p3, p4], p2, p3]*Spbb[p7, p6, 
                  p5, SumM[p3, p4], p2, p3])/Spbb[p7, p6, p1, p7] - Spbb[p3, 
                p2, SumM[p3, p4], p6, p5, SumM[p3, p4], p2, p3]) + 
             (-((2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*Spbb[p3, p4, p2, 
                  p3]) + 2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])*
              ((-2*MP[p1, p7]*Spbb[p7, SumM[p3, p4], p2, p3]*Spbb[p7, p6, 
                  SumM[p2, p3, p4], SumM[p3, p4], p2, p3])/Spbb[p7, p6, p1, 
                 p7] - Spbb[p3, p2, SumM[p3, p4], SumM[p1, p7], SumM[p2, p3, 
                 p4], SumM[p3, p4], p2, p3])) + 
           (4*MP[p3, p4]^2*Spbb[p3, SumM[p4, p5], p2, p3]*
              (-Spbb[p3, SumM[p1, p7], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
                 Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7]) + 
             Spbb[p3, p4, p2, p3]^2*(-((2*MP[p2, p6] + (2*MP[p1, p7]*
                    Spbb[p7, p6, p2, p7])/Spbb[p7, p6, p1, p7])*
                 (2*MP[p5, p5] + 2*MP[p5, p6] + (2*MP[p1, p7]*Spbb[p7, p6, 
                     p5, p7])/Spbb[p7, p6, p1, p7])) + (2*MP[p2, p3] + 
                 2*MP[p2, p4] + 2*MP[p3, p4])*(2*(2*MP[p3, p4] + 
                   2*MP[p3, p5] + 2*MP[p4, p5]) + 2*MP[p5, p5] + 
                 2*MP[p5, p6] + (2*MP[p1, p7]*Spbb[p7, p6, p5, p7])/
                  Spbb[p7, p6, p1, p7]) + (2*MP[p1, p1] + 2*MP[p1, p2] + 
                 2*MP[p1, p7] + 2*MP[p2, p7] + (2*MP[p1, p7]*Spbb[p7, 
                    SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, p7])*
                (-2*MP[p1, p5] - 2*MP[p1, p7] + 2*MP[p3, p4] + 2*MP[p3, p5] + 
                 2*MP[p4, p5] - 2*MP[p5, p7] - (2*MP[p1, p7]*Spbb[p7, 
                    SumM[p1, p5], p6, p7])/Spbb[p7, p6, p1, p7])) + 
             2*MP[p3, p4]*Spbb[p3, p4, p2, p3]*((-Spbb[p3, p6, p2, p3] - 
                 (2*MP[p1, p7]*Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/
                  Spbb[p7, p6, p1, p7])*(2*MP[p5, p5] + 2*MP[p5, p6] + 
                 (2*MP[p1, p7]*Spbb[p7, p6, p5, p7])/Spbb[p7, p6, p1, p7]) - 
               Spbb[p3, SumM[p4, p5], p2, p3]*(2*MP[p1, p1] + 2*MP[p1, p2] + 
                 2*MP[p1, p7] + 2*MP[p2, p7] + 2*(2*MP[p2, p3] + 
                   2*MP[p2, p4] + 2*MP[p3, p4]) + (2*MP[p1, p7]*Spbb[p7, 
                    SumM[p1, p2], p6, p7])/Spbb[p7, p6, p1, p7]) + 
               (-Spbb[p3, SumM[p1, p7], p2, p3] + (2*MP[p1, p7]*Spbb[p3, p7]*
                   Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7])*
                (2*MP[p1, p5] + 2*MP[p1, p7] - 2*MP[p3, p4] - 2*MP[p3, p5] - 
                 2*MP[p4, p5] + 2*MP[p5, p7] + (2*MP[p1, p7]*Spbb[p7, 
                    SumM[p1, p5], p6, p7])/Spbb[p7, p6, p1, p7])))*
            Spbb[p3, p2, SumM[p3, p4], SumM[p2, p3, p4], p5, SumM[p3, p4], 
             p2, p3]))/(4*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4])*Spbb[p3, p2, p3, p4]*Spbb[p3, p4, p2, p3]^3*
          (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] - 
           (2*MP[p3, p4]*Spbb[p3, SumM[p4, p5], p2, p3])/Spbb[p3, p4, p2, 
             p3])*(2*MP[p5, p5] + 2*MP[p5, p6] + 
           (2*MP[p1, p7]*Spbb[p7, p6, p5, p7])/Spbb[p7, p6, p1, p7])*
          (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7] - 
           (2*MP[p3, p4]*(-Spbb[p3, SumM[p1, p7], p2, p3] + (2*MP[p1, p7]*
                Spbb[p3, p7]*Spbb[p7, p6, p2, p3])/Spbb[p7, p6, p1, p7]))/
            Spbb[p3, p4, p2, p3] + (2*MP[p1, p7]*Spbb[p7, SumM[p1, p2], p6, 
              p7])/Spbb[p7, p6, p1, p7]))))/(2*MP[p1, p7]*Spab[p7, p6, p7])
 
ATreeSAM["S", "+", "+", "+", "S", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p1, p1]*(-((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
           2*MP[p6, p7])*Spbb[p2, p1, p6, p2]) + 2*MP[p1, p6]*
         Spbb[p2, SumM[p3, p4, p5], p1, p2])*
       ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
         Spbb[p2, p4] - Spbb[p2, SumM[p3, p4, p5], p2, p4] - 
        Spbb[p2, SumM[p3, p4, p5], p3, p4]))/(2*MP[p4, p5]*
       (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*(2*MP[p6, p6] + 
        2*MP[p6, p7])*(2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
        2*MP[p6, p7])*Spaa[p3, p4]*((2*MP[p1, p6] + 2*MP[p1, p7] + 
          2*MP[p6, p6] + 2*MP[p6, p7])*Spab[p3, p1, p2] + 
        Spaa[p2, p3]*Spbb[p2, SumM[p3, p4, p5], p1, p2])) + 
     (Spbb[p2, p3]^3*((Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2]*
          Spbb[p2, p3, p1, p2]*(Spbb[p4, SumM[p6, p7], p5, p4] + 
           ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
             (Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2]*(
                -Spbb[p2, p1, SumM[p2, p3], p5, SumM[p6, p7], p4] + 
                Spbb[p2, p1, SumM[p2, p3], SumM[p6, p7], p5, p4]) + 
              (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spbb[p2, p1, 
                SumM[p2, p3], SumM[p6, p7], p5, SumM[p2, p3], p1, p2]))/
            Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2]^2)*
          (1 + (2*MP[p1, p7] + 2*MP[p2, p3] + 2*MP[p2, p7] + 2*MP[p3, p7] - 
             ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p4, p7, 
                SumM[p2, p3], p1, p2])/Spab[p4, SumM[p1, p2, p3], SumM[p2, 
                p3], p1, p2] - (2*MP[p2, p3]*Spbb[p2, p1, p7, p2])/
              Spbb[p2, p3, p1, p2] - (2*MP[p2, p3]*Spbb[p2, SumM[p3, p7], p1, 
                p2])/Spbb[p2, p3, p1, p2] + (((2*MP[p1, p2] + 2*MP[p1, p3] + 
                  2*MP[p2, p3])*Spab[p4, p5, SumM[p2, p3], p1, p2] + 
                2*MP[p4, p5]*Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, 
                  p2])*(Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2]^2*
                 Spbb[p4, SumM[p5, p6], p7, p4] + (2*MP[p1, p2] + 
                  2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p4, SumM[p1, p2, p3], 
                  SumM[p2, p3], p1, p2]*(-Spbb[p2, p1, SumM[p2, p3], p7, 
                    SumM[p5, p6], p4] + Spbb[p2, p1, SumM[p2, p3], 
                   SumM[p5, p6], p7, p4]) + (2*MP[p1, p2] + 2*MP[p1, p3] + 
                   2*MP[p2, p3])^2*Spbb[p2, p1, SumM[p2, p3], SumM[p5, p6], 
                  p7, SumM[p2, p3], p1, p2]))/(Spab[p4, SumM[p1, p2, p3], 
                SumM[p2, p3], p1, p2]*(Spab[p4, SumM[p1, p2, p3], SumM[p2, 
                    p3], p1, p2]^2*Spbb[p4, SumM[p6, p7], p5, p4] + 
                (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p4, 
                  SumM[p1, p2, p3], SumM[p2, p3], p1, p2]*
                 (-Spbb[p2, p1, SumM[p2, p3], p5, SumM[p6, p7], p4] + 
                  Spbb[p2, p1, SumM[p2, p3], SumM[p6, p7], p5, p4]) + 
                (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])^2*Spbb[p2, p1, 
                  SumM[p2, p3], SumM[p6, p7], p5, SumM[p2, p3], p1, p2])))/
            (2*MP[p6, p6] + 2*MP[p6, p7])))/((2*MP[p1, p2] + 2*MP[p1, p3] + 
           2*MP[p2, p3])*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 
           2*MP[p6, p7])*(2*MP[p4, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 
              2*MP[p2, p3])*Spab[p4, p5, SumM[p2, p3], p1, p2])/
            Spab[p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2])) + 
        ((-(Spab[p4, p5, SumM[p2, p3], p1, p2]*Spbb[p2, p3, p5, p4]) + 
           2*MP[p4, p5]*((2*MP[p2, p3] + 2*MP[p2, p5] + 2*MP[p3, p5])*
              Spbb[p2, p3, p1, p2] - 2*MP[p2, p3]*Spbb[p2, SumM[p3, p5], p1, 
               p2]))*(Spbb[p2, p3, p1, p2]*((2*MP[p1, p7] + 2*MP[p6, p6] + 2*
                MP[p6, p7])*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], SumM[p1, 
                p2, p3], SumM[p2, p3], p1, p2] - (2*MP[p1, p6] + 2*
                MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*Spbb[p2, p1, SumM[
                p2, p3], SumM[p1, p2, p3], p7, SumM[p2, p3], p1, p2]) + 
           2*MP[p2, p3]*(-(Spbb[p2, p1, p7, p2]*Spbb[p2, p1, SumM[p2, p3], 
                SumM[p4, p5], SumM[p1, p2, p3], SumM[p2, p3], p1, p2]) + 
             Spbb[p2, SumM[p3, p4, p5], p1, p2]*Spbb[p2, p1, SumM[p2, p3], 
               SumM[p1, p2, p3], p7, SumM[p2, p3], p1, p2])))/
         (2*MP[p4, p5]*(2*MP[p6, p6] + 2*MP[p6, p7])*
          ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spab[p4, p5, 
             SumM[p2, p3], p1, p2] + 2*MP[p4, p5]*Spab[p4, SumM[p1, p2, p3], 
             SumM[p2, p3], p1, p2])*(-((2*MP[p1, p6] + 2*MP[p1, p7] + 
              2*MP[p6, p6] + 2*MP[p6, p7])*Spbb[p2, p3, p1, p2]) + 
           2*MP[p2, p3]*Spbb[p2, SumM[p3, p4, p5], p1, p2]))))/
      (2*MP[p2, p3]*Spab[p4, p3, p2]*Spbb[p2, p1, p2, p3]*
       Spbb[p2, p3, p1, p2])
 
ATreeSAM["S", "+", "+", "S", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = 
    -(((2*MP[p1, p6] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7])*Spbb[p2, p1, p6, p2])/Spbb[p2, SumM[p3, p4, p5], 
           p1, p2])*(-1/2*(MP[p1, p1]*((2*MP[p1, p6] + 2*MP[p1, p7] + 2*
                MP[p6, p6] + 2*MP[p6, p7])*Spbb[p2, p5] - 
             Spbb[p2, SumM[p3, p4, p5], p2, p5] - Spbb[p2, SumM[p3, p4, p5], 
              p3, p5]))/(MP[p4, p5]*(2*MP[p3, p4] + 2*MP[p3, p5] + 
             2*MP[p4, p5])*Spaa[p3, p5]*(-Spaa[p2, p3] - 
             ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
               Spab[p3, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])) + 
         (MP[p1, p1]*((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
              2*MP[p6, p7])*Spbb[p2, p3] - Spbb[p2, SumM[p3, p4, p5], p2, 
             p3] - Spbb[p2, SumM[p3, p4, p5], p5, p3]))/
          (2*MP[p3, p4]*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
           Spaa[p3, p5]*(-Spaa[p2, p5] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*
                MP[p6, p6] + 2*MP[p6, p7])*Spab[p5, p1, p2])/
             Spbb[p2, SumM[p3, p4, p5], p1, p2])) - 
         (MP[p1, p1]*(((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*
                MP[p6, p7])*Spbb[p2, p5]*Spbb[p2, p1, p2, p3])/
             Spbb[p2, SumM[p3, p4, p5], p1, p2] + 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spbb[p2, p5]*Spbb[p2, p1, p5, p3])/Spbb[p2, SumM[p3, p4, p5], 
              p1, p2] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*
                MP[p6, p7])*Spbb[p2, p3]*Spbb[p2, p1, SumM[p2, p3, p4], p5])/
             Spbb[p2, SumM[p3, p4, p5], p1, p2] + Spbb[p3, p2, 
             SumM[p2, p3, p4], p5] + Spbb[p3, p5, SumM[p2, p3, p4], p5]))/
          (2*MP[p3, p4]*(-Spaa[p2, p3] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*
                MP[p6, p6] + 2*MP[p6, p7])*Spab[p3, p1, p2])/
             Spbb[p2, SumM[p3, p4, p5], p1, p2])*(-Spaa[p2, p5] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spab[p5, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])*
           (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, 
              p2]))))/((2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p6] + 
         2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7]))) + 
     (MP[p1, p1]*Spbb[p2, p3]^2*((2*MP[p1, p5] + 2*MP[p1, p7] + 
          2*MP[p5, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
            Spbb[p2, SumM[p3, p4, p6], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
            p2])*((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, SumM[p3, p4, p5], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2])*Spbb[p5, p6, p7, p5] + 2*MP[p5, p6]*
           ((MP[p1, p1]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, p5]^2)/Spbb[p2, SumM[p3, p4], p1, p2] + 
            ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p5]*
              Spbb[p2, p1, SumM[p2, p3, p4], p5])/Spbb[p2, SumM[p3, p4], p1, 
              p2] - Spbb[p5, SumM[p2, p3, p4], p1, p5])) - 
        (-2*MP[p5, p6]*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7]) + (2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7])*(2*MP[p1, p6] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4])*Spbb[p2, p1, p6, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2]) - 4*MP[p5, p6]*(2*MP[p1, p6] + 2*MP[p1, p7] + 
            2*MP[p6, p6] + 2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4])*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
             Spbb[p2, SumM[p3, p4], p1, p2]) - (2*MP[p6, p6] + 2*MP[p6, p7])*
           (2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7] - 
            ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, SumM[p3, 
                p4, p5], p1, p2])/Spbb[p2, SumM[p3, p4], p1, p2]) + 
          (2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p5] + 2*MP[p1, p6] + 
            2*MP[p5, p6] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, SumM[p3, p4, p7], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2]))*(-(((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
             Spbb[p2, p5]*Spbb[p2, p1, p6, p5])/Spbb[p2, SumM[p3, p4], p1, 
             p2]) + Spbb[p5, SumM[p2, p3, p4], p6, p5])))/
      (8*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*MP[p5, p6]*
       (2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p5, p6] + 2*MP[p5, p7] + 
        2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p2, p3] - 
        ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p3, p1, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])*(2*MP[p1, p6] + 2*MP[p1, p7] + 
        2*MP[p6, p6] + 2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4])*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])) + 
     (Spbb[p2, p3]^3*Spbb[p2, p3, p1, p2]*
       ((((MP[p1, p1]*Spbb[p5, p6, p7, p5]^2 + MP[p6, p6]*
              (-Spbb[p5, p1, p4, p5] - Spbb[p5, SumM[p2, p3], p4, p5] - 
                ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spbb[p2, p1, 
                   SumM[p2, p3], p5]*(-(MP[p1, p1]*Spbb[p2, p1, SumM[p2, p3], 
                      p5]) - Spbb[p2, p1, SumM[p2, p3], p4, p1, p5] - 
                   Spbb[p2, p1, SumM[p2, p3], p4, SumM[p2, p3], p5]))/
                 Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, 
                  SumM[p2, p3], p1, p2])^2)/((2*MP[p6, p6] + 2*MP[p6, p7])*
             (2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
             ((2*MP[p6, p6] + 2*MP[p6, p7])*(Spbb[p5, p4, p6, p5] + 
                ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spbb[p2, p1, 
                   SumM[p2, p3], p5]*Spbb[p2, p1, SumM[p2, p3], p4, p6, p5])/
                 Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, 
                  SumM[p2, p3], p1, p2]) - Spbb[p5, p6, p7, p5]*(
                2*MP[p4, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
                  Spbb[p2, p1, SumM[p2, p3], p5, p4, SumM[p2, p3], p1, p2])/
                 Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, 
                  SumM[p2, p3], p1, p2]))) - 
           ((1 + ((Spbb[p5, p4, p6, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 
                    2*MP[p2, p3])*Spbb[p2, p1, SumM[p2, p3], p5]*Spbb[p2, p1, 
                    SumM[p2, p3], p4, p6, p5])/Spbb[p2, p1, SumM[p2, p3], 
                   SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2])*
                (2*(MP[p4, p5] + MP[p4, p6] + MP[p5, p6]) + 
                 ((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*Spbb[p2, p1, 
                    SumM[p2, p3], SumM[p5, p6], p4, SumM[p2, p3], p1, p2])/
                  Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, 
                   SumM[p2, p3], p1, p2]))/((2*MP[p6, p6] + 2*MP[p6, p7])*
                 (Spbb[p5, p4, p6, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 
                     2*MP[p2, p3])*Spbb[p2, p1, SumM[p2, p3], p5]*Spbb[p2, 
                     p1, SumM[p2, p3], p4, p6, p5])/Spbb[p2, p1, SumM[p2, 
                     p3], SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2]) - 
                Spbb[p5, p6, p7, p5]*(2*MP[p4, p5] + ((2*MP[p1, p2] + 
                     2*MP[p1, p3] + 2*MP[p2, p3])*Spbb[p2, p1, SumM[p2, p3], 
                     p5, p4, SumM[p2, p3], p1, p2])/Spbb[p2, p1, SumM[p2, 
                     p3], SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2])))*
             (Spbb[p5, p4, p6, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 
                 2*MP[p2, p3])*Spbb[p2, p1, SumM[p2, p3], p5]*Spbb[p2, p1, 
                 SumM[p2, p3], p4, p6, p5])/Spbb[p2, p1, SumM[p2, p3], 
                SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2]))/
            (2*MP[p5, p6]*(2*MP[p4, p5] + ((2*MP[p1, p2] + 2*MP[p1, p3] + 
                 2*MP[p2, p3])*Spbb[p2, p1, SumM[p2, p3], p5, p4, 
                 SumM[p2, p3], p1, p2])/Spbb[p2, p1, SumM[p2, p3], 
                SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2])))*
          Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, SumM[p2, p3], p1, 
           p2])/((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          Spbb[p2, p3, p1, p2]^2*(2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4] - (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
            Spbb[p2, p3, p1, p2])) - (MP[p1, p1]*Spbb[p2, p1, SumM[p2, p3], 
           p5]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], SumM[p1, p2, p3], 
           SumM[p2, p3], p1, p2]*(1 + (2*MP[p4, p6] + 2*MP[p5, p6] + 
             (2*MP[p4, p5]*Spbb[p2, p1, SumM[p2, p3], p4, p6, SumM[p2, p3], 
                p1, p2])/Spbb[p2, p1, SumM[p2, p3], p5, p4, SumM[p2, p3], p1, 
               p2] + (Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], p6, 
                SumM[p2, p3], p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 
                2*MP[p2, p3] + (2*MP[p4, p5]*Spbb[p2, p1, SumM[p2, p3], 
                   SumM[p1, p2, p3], p4, SumM[p2, p3], p1, p2])/Spbb[p2, p1, 
                  SumM[p2, p3], p5, p4, SumM[p2, p3], p1, p2]))/
              Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], SumM[p1, p2, p3], SumM[
                p2, p3], p1, p2])/(2*MP[p6, p6] + 2*MP[p6, p7])))/
         (2*MP[p4, p5]*Spab[p5, p4, SumM[p2, p3], p1, p2]*
          Spbb[p2, p3, p1, p2]^2*(2*MP[p1, p6] + 2*MP[p1, p7] + 
           2*MP[p6, p6] + 2*MP[p6, p7] - (2*MP[p2, p3]*Spbb[p2, 
              SumM[p3, p4, p5], p1, p2])/Spbb[p2, p3, p1, p2])*
          (2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3] + 
           (2*MP[p4, p5]*Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, 
              SumM[p2, p3], p1, p2])/Spbb[p2, p1, SumM[p2, p3], p5, p4, 
             SumM[p2, p3], p1, p2]))))/(2*MP[p2, p3]*Spbb[p2, p1, p2, p3])
 
ATreeSAM["S", "+", "S", "+", "+", "s", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = 
    -(((2*MP[p1, p6] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7])*Spbb[p2, p1, p6, p2])/Spbb[p2, SumM[p3, p4, p5], 
           p1, p2])*(-1/2*(MP[p1, p1]*Spbb[p2, SumM[p4, p5], p3, p4])/
           (MP[p3, p4]*(2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*
            Spaa[p4, p5]*(-Spaa[p2, p5] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 
                2*MP[p6, p6] + 2*MP[p6, p7])*Spab[p5, p1, p2])/
              Spbb[p2, SumM[p3, p4, p5], p1, p2])) + 
         (MP[p1, p1]*(-Spbb[p2, p3, p2, p5] - Spbb[p2, p3, p4, p5] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spbb[p2, p5]*Spbb[p2, p1, p3, p2])/Spbb[p2, SumM[p3, p4, p5], 
              p1, p2]))/(Spaa[p4, p5]*(-Spaa[p2, p4] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spab[p4, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])*
           (2*MP[p2, p3] + ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*
                MP[p6, p7])*Spbb[p2, p1, p3, p2])/Spbb[p2, SumM[p3, p4, p5], 
              p1, p2])*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, 
              p2])) + (MP[p1, p1]*(-(((2*MP[p1, p6] + 2*MP[p1, p7] + 
                2*MP[p6, p6] + 2*MP[p6, p7])*Spbb[p2, p5]*Spbb[p2, p1, p3, 
                p4])/Spbb[p2, SumM[p3, p4, p5], p1, p2]) - Spbb[p4, p3, p2, 
             p5] - Spbb[p4, p3, p4, p5]))/(2*MP[p3, p4]*(-Spaa[p2, p4] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spab[p4, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, p2])*
           (-Spaa[p2, p5] - ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*
                MP[p6, p7])*Spab[p5, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, 
              p2])*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] - 
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7])*
              Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, 
              p2]))))/((2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p6] + 
         2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7]))) + 
     (MP[p1, p1]*Spbb[p2, p4]^2*((2*MP[p1, p5] + 2*MP[p1, p7] + 
          2*MP[p5, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
            Spbb[p2, SumM[p3, p4, p6], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
            p2])*((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, SumM[p3, p4, p5], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2])*Spbb[p5, p6, p7, p5] + 2*MP[p5, p6]*
           ((MP[p1, p1]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, p5]^2)/Spbb[p2, SumM[p3, p4], p1, p2] + 
            ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p5]*
              Spbb[p2, p1, SumM[p2, p3, p4], p5])/Spbb[p2, SumM[p3, p4], p1, 
              p2] - Spbb[p5, SumM[p2, p3, p4], p1, p5])) - 
        (-2*MP[p5, p6]*(2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7]) + (2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p6] + 
            2*MP[p6, p7])*(2*MP[p1, p6] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4])*Spbb[p2, p1, p6, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2]) - 4*MP[p5, p6]*(2*MP[p1, p6] + 2*MP[p1, p7] + 
            2*MP[p6, p6] + 2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4])*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
             Spbb[p2, SumM[p3, p4], p1, p2]) - (2*MP[p6, p6] + 2*MP[p6, p7])*
           (2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7] - 
            ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, SumM[p3, 
                p4, p5], p1, p2])/Spbb[p2, SumM[p3, p4], p1, p2]) + 
          (2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p5] + 2*MP[p1, p6] + 
            2*MP[p5, p6] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, SumM[p3, p4, p7], p1, p2])/Spbb[p2, SumM[p3, p4], p1, 
              p2]))*(-(((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
             Spbb[p2, p5]*Spbb[p2, p1, p6, p5])/Spbb[p2, SumM[p3, p4], p1, 
             p2]) + Spbb[p5, SumM[p2, p3, p4], p6, p5])))/
      (8*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*MP[p5, p6]*
       (2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p5, p6] + 2*MP[p5, p7] + 
        2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p2, p3] + 
        ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p1, p3, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])*(2*MP[p1, p6] + 2*MP[p1, p7] + 
        2*MP[p6, p6] + 2*MP[p6, p7] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4])*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])) - 
     (Spbb[p2, p1, p3, p2]*((MP[p6, p6]*Spbb[p4, p5]^2*(2*MP[p2, p3] + 
           2*MP[p2, p7] + 2*MP[p3, p7] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p7], 
              p1, p2])/Spbb[p2, p1, p3, p2] - 
           ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
             ((-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/Spbb[p2, p1, p3, 
                p2] - (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p2, p3], 
                 p4])/Spbb[p2, p1, p3, p2] - (2*MP[p2, p3]*Spbb[p2, p4]*
                Spbb[p2, p1, SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - 
              Spbb[p4, SumM[p2, p3], p1, p4] - Spbb[p4, SumM[p2, p3], SumM[
                p5, p6], p4]))/((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, 
                SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - 
             Spbb[p4, SumM[p2, p3], SumM[p5, p6], p4])))/
         (2*MP[p5, p6]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
          (2*MP[p4, p5] - ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
             ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/Spbb[p2, p1, 
                p3, p2] + Spbb[p4, p5, SumM[p2, p3], p4]))/
            ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6], p4])/
              Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], SumM[p5, p6], 
              p4]))*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 
           2*MP[p2, p3] - ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
             ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/Spbb[p2, p1, p3, p2] + 
              (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[
                p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], p1, p4]))/
            ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6], p4])/
              Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], SumM[p5, p6], 
              p4]))) - (MP[p6, p6]*((-2*MP[p2, p3]*Spbb[p2, p4]*
             Spbb[p2, p1, SumM[p5, p6, p7], p4])/Spbb[p2, p1, p3, p2] - 
           Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], p4])*
          ((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6, p7], p4])/
            Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p5, p6, p7], p1, p4])*
          ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p7, p4, p5])/
            Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], 
              p4]*Spbb[p2, p7, p4, p5])/Spbb[p2, p1, p3, p2] - 
           (2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p7, 
              SumM[p4, p6, p7], p5])/Spbb[p2, p1, p3, p2] - 
           (2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], p4]*Spbb[p2, p7, 
              SumM[p4, p6, p7], p5])/Spbb[p2, p1, p3, p2] + 
           MP[p1, p1]*Spbb[p4, p7, p4, p5] - MP[p1, p1]*Spbb[p4, p7, 
             SumM[p4, p6, p7], p5] + Spbb[p4, SumM[p2, p3], p1, p7, p4, p5] - 
           Spbb[p4, SumM[p2, p3], p1, p7, SumM[p4, p6, p7], p5]))/
         (2*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          MP[p5, p6]*(-2*MP[p1, p1] - 2*MP[p1, p2] - 2*MP[p1, p3] - 
           2*MP[p2, p3] + Spab[p4, SumM[p5, p6, p7], p4])*
          ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
            ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6], p4])/
              Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], SumM[p5, p6], 
              p4]) - (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6])*
            ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6, p7], p4])/
              Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], 
              p4]))*((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 
             2*MP[p2, p3])*(Spab[p5, SumM[p2, p3], p4] + 
             (2*MP[p2, p3]*Spab[p5, p1, p2]*Spbb[p2, p4])/Spbb[p2, p1, p3, 
               p2]) + Spaa[p4, p5]*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, 
                SumM[p5, p6, p7], p4])/Spbb[p2, p1, p3, p2] - 
             Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], p4]))) - 
        (MP[p1, p1]*((-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/
            Spbb[p2, p1, p3, p2] - (2*MP[p2, p3]*Spbb[p2, p4]*
             Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] - 
           Spbb[p4, SumM[p2, p3], p1, p4])*((MP[p6, p6]*Spbb[p4, p5]^3)/
            (2*MP[p6, p6] + 2*MP[p6, p7]) + 
           ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*(
                (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p6, p5])/Spbb[p2, p1, 
                  p3, p2] + Spbb[p4, SumM[p2, p3], p6, p5])*Spbb[p4, 
                SumM[p5, p6], p7, p4] + ((-2*MP[p2, p3]*Spbb[p2, p4]*
                  Spbb[p2, p1, SumM[p5, p6, p7], p4])/Spbb[p2, p1, p3, p2] - 
                Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], p4])*(Spbb[p4, 
                 SumM[p5, p6, p7], p7, p4, p6, p5] + Spbb[p4, SumM[p5, p6, 
                  p7], p7, p5, p6, p5]))^2/(2*MP[p5, p6]*(-2*MP[p1, p1] - 
              2*MP[p1, p2] - 2*MP[p1, p3] - 2*MP[p2, p3] + Spab[p4, SumM[p5, 
                p6, p7], p4])*((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 
                2*MP[p2, p3])*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, 
                   SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, 
                 SumM[p2, p3], SumM[p5, p6], p4]) - (2*MP[p4, p5] + 
                2*MP[p4, p6] + 2*MP[p5, p6])*((-2*MP[p2, p3]*Spbb[p2, p4]*
                  Spbb[p2, p1, SumM[p5, p6, p7], p4])/Spbb[p2, p1, p3, p2] - 
                Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], p4]))*
             ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*(
                Spab[p5, SumM[p2, p3], p4] + (2*MP[p2, p3]*Spab[p5, p1, p2]*
                  Spbb[p2, p4])/Spbb[p2, p1, p3, p2]) + Spaa[p4, p5]*(
                (-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6, p7], 
                   p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], 
                 SumM[p5, p6, p7], p4])))))/((2*MP[p1, p1] + 2*MP[p1, p2] + 
           2*MP[p1, p3] + 2*MP[p2, p3])*((-2*MP[p2, p3]*Spbb[p2, p4]*
             Spbb[p2, p1, SumM[p5, p6, p7], p4])/Spbb[p2, p1, p3, p2] - 
           Spbb[p4, SumM[p2, p3], SumM[p5, p6, p7], p4])*
          ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p4, p6, p5])/
            Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], 
              p4]*Spbb[p2, p4, p6, p5])/Spbb[p2, p1, p3, p2] + 
           (2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p5, p6, p5])/
            Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], 
              p4]*Spbb[p2, p5, p6, p5])/Spbb[p2, p1, p3, p2] + 
           MP[p1, p1]*Spbb[p4, p5, p6, p5] + (2*MP[p1, p1] + 2*MP[p1, p2] + 
             2*MP[p1, p3] + 2*MP[p2, p3])*((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[
                p2, p1, p6, p5])/Spbb[p2, p1, p3, p2] + 
             Spbb[p4, SumM[p2, p3], p6, p5]) + Spbb[p4, SumM[p2, p3], p1, p4, 
            p6, p5] + Spbb[p4, SumM[p2, p3], p1, p5, p6, p5])) + 
        (Spbb[p4, p5]^3*(((2*MP[p2, p3] + 2*MP[p2, p6] + 2*MP[p3, p6] + 
               (2*MP[p2, p3]*Spbb[p2, SumM[p3, p6], p1, p2])/Spbb[p2, p1, p3, 
                 p2])*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/
                Spbb[p2, p1, p3, p2] + Spbb[p4, p5, SumM[p2, p3], p4]) - 
             2*MP[p4, p5]*((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p6, p4])/
                Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], p6, p4]))*
            ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7] + 
               (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
                Spbb[p2, p1, p3, p2])*((-2*MP[p2, p3]*Spbb[p2, p4]*
                 Spbb[p2, p1, p5, p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, p5, 
                SumM[p2, p3], p4])*((2*MP[p2, p3]*Spbb[p2, p4]*
                 ((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p6, 
                     p7, SumM[p4, p5], p1, p2])/Spbb[p2, p1, p3, p2] + 
                  Spbb[p2, p1, SumM[p4, p5], p6, p7, SumM[p4, p5], 
                   SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] - (2*MP[p2, p3]*
                 Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p7, p6, SumM[p4, 
                   p5], SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, 
                SumM[p2, p3], SumM[p4, p5], p6, p7, SumM[p4, p5], 
                SumM[p2, p3], p4]) + (-((2*MP[p4, p5] + 2*MP[p4, p6] + 
                  2*MP[p5, p6])*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, 
                     p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, p5, SumM[p2, p3], 
                   p4])) + 2*MP[p4, p5]*((-2*MP[p2, p3]*Spbb[p2, p4]*
                   Spbb[p2, p1, SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - 
                 Spbb[p4, SumM[p2, p3], SumM[p5, p6], p4]))*
              ((2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], SumM[p2, p3], p4]*
                 ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, SumM[p4, 
                      p5], p1, p2])/Spbb[p2, p1, p3, p2] + MP[p1, p1]*
                   Spbb[p2, SumM[p4, p5], SumM[p2, p3], p4]))/Spbb[p2, p1, 
                 p3, p2] + (2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], SumM[p2, p3], 
                  p4]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], SumM[p2, p3], 
                  p4])/Spbb[p2, p1, p3, p2] + 2*MP[p4, p5]*
                ((2*MP[p2, p3]*Spbb[p2, p1, SumM[p2, p3], p4]*Spbb[p2, 
                    SumM[p4, p5], SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + 
                 (2*MP[p2, p3]*Spbb[p2, p4]*((2*MP[p1, p1]*MP[p2, p3]*
                      Spbb[p2, p4]*Spbb[p2, SumM[p4, p5], p1, p2])/Spbb[p2, 
                      p1, p3, p2] + MP[p1, p1]*Spbb[p2, SumM[p4, p5], 
                      SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] + 
                 (2*MP[p2, p3]*Spbb[p2, p4]*((2*MP[p2, p3]*Spbb[p2, p1, 
                       SumM[p2, p3], p4]*Spbb[p2, SumM[p4, p5], p1, p2])/
                     Spbb[p2, p1, p3, p2] - Spbb[p2, p1, SumM[p4, p5], p1, 
                     SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] - Spbb[p4, 
                  SumM[p2, p3], SumM[p4, p5], p1, SumM[p2, p3], p4]) + 
               (2*MP[p2, p3]*Spbb[p2, p4]*((2*MP[p2, p3]*Spbb[p2, SumM[p4, 
                      p5], p1, p2]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], 
                     SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] - Spbb[p2, p1, 
                   SumM[p4, p5], p1, SumM[p2, p3], SumM[p4, p5], SumM[p2, 
                    p3], p4]))/Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*
                 Spbb[p2, p4]*((2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], p1, p2]*
                    ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, SumM[p4, 
                         p5], p1, p2])/Spbb[p2, p1, p3, p2] + MP[p1, p1]*
                      Spbb[p2, SumM[p4, p5], SumM[p2, p3], p4]))/Spbb[p2, p1, 
                    p3, p2] + (2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], SumM[p2, 
                      p3], p4]*Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], p1, 
                     p2])/Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p4]*
                    ((2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], p1, p2]*Spbb[p2, 
                        p1, SumM[p2, p3], SumM[p4, p5], p1, p2])/Spbb[p2, p1, 
                       p3, p2] + Spbb[p2, p1, SumM[p4, p5], SumM[p2, p3], p1, 
                      SumM[p4, p5], p1, p2]))/Spbb[p2, p1, p3, p2] + 
                  Spbb[p2, p1, SumM[p4, p5], SumM[p2, p3], p1, SumM[p4, p5], 
                   SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] + Spbb[p4, 
                SumM[p2, p3], SumM[p4, p5], SumM[p2, p3], p1, SumM[p4, p5], 
                SumM[p2, p3], p4])) + (((2*MP[p1, p1] + 2*MP[p1, p2] + 
                 2*MP[p1, p3] + 2*MP[p2, p3])*(-2*MP[p1, p6] + 2*MP[p4, p5] + 
                 2*MP[p4, p6] + 2*MP[p5, p6] - (2*MP[p2, p3]*Spbb[p2, p1, p6, 
                    p2])/Spbb[p2, p1, p3, p2]) - (2*MP[p6, p6] + 
                 2*MP[p6, p7])*(2*MP[p2, p3] + 2*MP[p2, p7] + 2*MP[p3, p7] + 
                 (2*MP[p2, p3]*Spbb[p2, SumM[p3, p7], p1, p2])/Spbb[p2, p1, 
                   p3, p2]) + (2*(2*MP[p4, p5] + 2*MP[p4, p6] + 
                   2*MP[p5, p6]) + 2*MP[p6, p6] + 2*MP[p6, p7])*
                (2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 2*MP[p6, p7] + 
                 (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
                  Spbb[p2, p1, p3, p2]))*((-2*MP[p2, p3]*Spbb[p2, p4]*
                  Spbb[p2, p1, p5, p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, p5, 
                 SumM[p2, p3], p4])^2 + 2*MP[p4, p5]*((-2*MP[p2, p3]*
                 Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/Spbb[p2, p1, p3, p2] + 
               Spbb[p4, p5, SumM[p2, p3], p4])*((2*MP[p1, p6] - 
                 2*MP[p4, p5] - 2*MP[p4, p6] - 2*MP[p5, p6] + 
                 (2*MP[p2, p3]*Spbb[p2, p1, p6, p2])/Spbb[p2, p1, p3, p2])*
                ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/Spbb[p2, p1, p3, 
                   p2] + (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p2, 
                     p3], p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], 
                  p1, p4]) + (2*MP[p6, p6] + 2*MP[p6, p7])*
                ((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p7, p4])/
                  Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], p7, p4]) - 
               (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3] + 
                 2*(2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p6] + 
                   2*MP[p6, p7] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4, p5], 
                      p1, p2])/Spbb[p2, p1, p3, p2]))*((-2*MP[p2, p3]*
                   Spbb[p2, p4]*Spbb[p2, p1, SumM[p5, p6], p4])/Spbb[p2, p1, 
                   p3, p2] - Spbb[p4, SumM[p2, p3], SumM[p5, p6], p4])) + 
             4*MP[p4, p5]^2*((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/
                Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p4]*
                 Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + Spbb[
                p4, SumM[p2, p3], p1, p4])*((-2*MP[p2, p3]*Spbb[p2, p4]*
                 Spbb[p2, p1, SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - Spbb[
                p4, SumM[p2, p3], SumM[p5, p6], p4]))*
            ((2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], SumM[p2, p3], p4]*(
                Spbb[p2, p1, p6, SumM[p4, p5], SumM[p2, p3], p4] - 
                (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p6, p1, 
                   p2])/Spbb[p2, p1, p3, p2]))/Spbb[p2, p1, p3, p2] + 
             2*MP[p4, p5]*((2*MP[p2, p3]*Spbb[p2, p4]*(Spbb[p2, p1, p6, 
                   SumM[p4, p5], SumM[p2, p3], p4] - (2*MP[p2, p3]*Spbb[p2, 
                     p4]*Spbb[p2, p1, SumM[p4, p5], p6, p1, p2])/Spbb[p2, p1, 
                    p3, p2]))/Spbb[p2, p1, p3, p2] - (2*MP[p2, p3]*
                 Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p6, SumM[p2, p3], 
                  p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], 
                SumM[p4, p5], p6, SumM[p2, p3], p4]) - 
             (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p6, 
                SumM[p2, p3], SumM[p4, p5], SumM[p2, p3], p4])/
              Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p4]*(
                (2*MP[p2, p3]*Spbb[p2, SumM[p4, p5], p1, p2]*(Spbb[p2, p1, 
                    p6, SumM[p4, p5], SumM[p2, p3], p4] - (2*MP[p2, p3]*
                     Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], p6, p1, p2])/
                    Spbb[p2, p1, p3, p2]))/Spbb[p2, p1, p3, p2] + 
                (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, SumM[p4, p5], 
                   SumM[p2, p3], p6, SumM[p4, p5], p1, p2])/Spbb[p2, p1, p3, 
                  p2] + Spbb[p2, p1, SumM[p4, p5], SumM[p2, p3], p6, 
                 SumM[p4, p5], SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] + 
             Spbb[p4, SumM[p2, p3], SumM[p4, p5], SumM[p2, p3], p6, 
              SumM[p4, p5], SumM[p2, p3], p4])))/(4*MP[p4, p5]*
          (2*MP[p6, p6] + 2*MP[p6, p7])*(2*MP[p1, p6] + 2*MP[p1, p7] + 
           2*MP[p6, p6] + 2*MP[p6, p7] + (2*MP[p2, p3]*Spbb[p2, 
              SumM[p3, p4, p5], p1, p2])/Spbb[p2, p1, p3, p2])*
          ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/
             Spbb[p2, p1, p3, p2] + Spbb[p4, p5, SumM[p2, p3], p4])^3*
          (2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3] - 
           (2*MP[p4, p5]*((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/Spbb[p2, 
                p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, 
                 SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, 
                p3], p1, p4]))/((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, 
                p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, p5, SumM[p2, p3], p4]))*
          ((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p4, p5])/
            Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], p4, p5])*
          (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p6] - 
           (2*MP[p4, p5]*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, 
                 SumM[p5, p6], p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, 
                p3], SumM[p5, p6], p4]))/((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[
                p2, p1, p5, p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, p5, 
              SumM[p2, p3], p4])))))/(2*MP[p2, p3]*Spab[p2, p1, p2])
 
ATreeSAM["S", "+", "+", "S", "s", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/2*(Spbb[p6, p5, p7, p6]*
       (((-4*MP[p2, p4]*MP[p3, p4] + Spab[p2, p4, p3]*Spab[p3, p4, p2])*
          (-2*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*MP[p6, p7]*
            Spbb[p2, p6]*Spbb[p6, p5, p1, p2] + 
           ((2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p5, p5] + 2*MP[p5, p6] + 2*
                MP[p5, p7] + 4*MP[p6, p7])*Spbb[p2, SumM[p3, p4], p1, p2] + 
             (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, SumM[p6, 
                p7], p1, p2])*Spbb[p6, p5, p7, p6] + 2*MP[p6, p7]*
            Spbb[p2, SumM[p3, p4], p1, p2]*Spbb[p6, SumM[p1, p7], p5, p6]))/
         (2*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
          (2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*
          (2*MP[p3, p4]*Spab[p3, p1, p2] + 2*MP[p1, p2]*Spab[p3, p4, p2])*
          Spbb[p6, p5, p7, p6]) + (Spab[p3, p1, p2]*
          ((4*MP[p1, p2]^2*(Spbb[p2, p5, p4, p2] + Spbb[p2, SumM[p6, p7], p4, 
               p2]))/Spab[p3, p1, p2]^2 + (2*MP[p1, p2]*
             (-Spbb[p2, p4, p5, p3] - Spbb[p2, p4, SumM[p6, p7], p3] + 
              Spbb[p2, p5, p4, p3] + Spbb[p2, SumM[p6, p7], p4, p3]))/
            Spab[p3, p1, p2] + Spbb[p3, p5, p4, p3] + Spbb[p3, SumM[p6, p7], 
            p4, p3])*(1 + (2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p2, p6] + 
             2*MP[p2, p7] + 4*MP[p6, p7] - (2*MP[p1, p2]*(Spab[p3, 
                 SumM[p6, p7], p2] - (2*MP[p6, p7]*Spab[p3, p5, p6]*
                  Spbb[p2, p6])/Spbb[p6, p5, p7, p6]))/Spab[p3, p1, p2] - 
             ((2*MP[p3, p4]*Spab[p3, p1, p2] + 2*MP[p1, p2]*Spab[p3, p4, p2])*
               (4*MP[p1, p2]^2*((Spbb[p2, SumM[p6, p7], p4, p2] + Spbb[p2, 
                     SumM[p6, p7], p5, p2])*Spbb[p6, p5, p7, p6] - 
                  2*MP[p6, p7]*Spbb[p2, p6]*(-(MP[p5, p5]*Spbb[p2, p6]) + 
                    Spbb[p6, p5, p4, p2] + Spbb[p6, p5, SumM[p6, p7], p2])) + 
                Spab[p3, p1, p2]^2*((Spbb[p3, SumM[p6, p7], p4, p3] + 
                    Spbb[p3, SumM[p6, p7], p5, p3])*Spbb[p6, p5, p7, p6] - 
                  2*MP[p6, p7]*Spbb[p3, p6]*(-(MP[p5, p5]*Spbb[p3, p6]) + 
                    Spbb[p6, p5, p4, p3] + Spbb[p6, p5, SumM[p6, p7], p3])) + 
                2*MP[p1, p2]*Spab[p3, p1, p2]*((-Spbb[p2, p4, SumM[p6, p7], 
                      p3] - Spbb[p2, p5, SumM[p6, p7], p3] + Spbb[p2, 
                     SumM[p6, p7], p4, p3] + Spbb[p2, SumM[p6, p7], p5, p3])*
                   Spbb[p6, p5, p7, p6] + 2*MP[p6, p7]*(-(Spbb[p2, p6]*
                      Spbb[p6, p5, p4, p3]) - Spbb[p3, p6]*(-2*MP[p5, p5]*
                       Spbb[p2, p6] + Spbb[p6, p5, p4, p2] + Spbb[p6, p5, 
                       SumM[p6, p7], p2]) - Spbb[p2, p6]*Spbb[p6, p5, 
                      SumM[p6, p7], p3]))))/(Spab[p3, p1, p2]*(4*MP[p1, p2]^2*
                 (Spbb[p2, p5, p4, p2] + Spbb[p2, SumM[p6, p7], p4, p2]) + 
                2*MP[p1, p2]*Spab[p3, p1, p2]*(-Spbb[p2, p4, p5, p3] - 
                  Spbb[p2, p4, SumM[p6, p7], p3] + Spbb[p2, p5, p4, p3] + 
                  Spbb[p2, SumM[p6, p7], p4, p3]) + Spab[p3, p1, p2]^2*
                 (Spbb[p3, p5, p4, p3] + Spbb[p3, SumM[p6, p7], p4, p3]))*
               Spbb[p6, p5, p7, p6]) + (2*MP[p6, p7]*Spbb[p6, SumM[p1, p7], 
                p5, p6])/Spbb[p6, p5, p7, p6] + (2*MP[p6, p7]*Spbb[p6, 
                SumM[p2, p7], p5, p6])/Spbb[p6, p5, p7, p6])/
            (2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])))/
         (2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          (2*MP[p3, p4] + (2*MP[p1, p2]*Spab[p3, p4, p2])/
            Spab[p3, p1, p2]))))/(MP[p6, p7]*Spaa[p2, p3]*Spab[p6, p5, p6])
 
ATreeSAM["S", "+", "S", "+", "s", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/2*(Spbb[p2, p1, p3, p2]*
        (-1/4*(Spbb[p6, p5, p7, p6]*(((MP[p5, p5]*(2*MP[p2, p3]*Spbb[p2, p4]*
                    (MP[p1, p1]*Spbb[p2, p4] + Spbb[p2, p1, SumM[p2, p3], 
                      p4]) + Spbb[p2, p1, p3, p2]*Spbb[p4, SumM[p2, p3], p1, 
                     p4])^2)/Spbb[p2, p1, p3, p2]^2 + MP[p1, p1]*
                (Spbb[p4, SumM[p6, p7], p5, p4] + (2*MP[p6, p7]*Spbb[p4, p6]*
                    Spbb[p4, SumM[p5, p6, p7], p5, p6])/Spbb[p6, p5, p7, p6])^
                 2)/((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 
                2*MP[p2, p3])*(2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 
                2*MP[p6, p7])*((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 
                  2*MP[p2, p3])*((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, 
                     p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], p5, 
                   p4] - (2*MP[p6, p7]*Spbb[p4, p6]*((2*MP[p2, p3]*Spbb[p2, 
                        p4]*Spbb[p2, p1, p5, p6])/Spbb[p2, p1, p3, p2] + 
                     Spbb[p4, SumM[p2, p3], p5, p6]))/Spbb[p6, p5, p7, p6]) + 
                ((-2*MP[p2, p3]*Spbb[p2, p4]*(MP[p1, p1]*Spbb[p2, p4] + 
                     Spbb[p2, p1, SumM[p2, p3], p4]))/Spbb[p2, p1, p3, p2] - 
                  Spbb[p4, SumM[p2, p3], p1, p4])*(2*MP[p4, p5] - 
                  (2*MP[p6, p7]*Spbb[p6, p4, p5, p6])/Spbb[p6, p5, p7, 
                    p6]))) - (((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, 
                   p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], p5, 
                 p4] - (2*MP[p6, p7]*Spbb[p4, p6]*((2*MP[p2, p3]*Spbb[p2, p4]*
                     Spbb[p2, p1, p5, p6])/Spbb[p2, p1, p3, p2] + Spbb[p4, 
                    SumM[p2, p3], p5, p6]))/Spbb[p6, p5, p7, p6])*(1 + 
                (((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/
                    Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], p5, p4] - 
                   (2*MP[p6, p7]*Spbb[p4, p6]*((2*MP[p2, p3]*Spbb[p2, p4]*
                        Spbb[p2, p1, p5, p6])/Spbb[p2, p1, p3, p2] + 
                      Spbb[p4, SumM[p2, p3], p5, p6]))/Spbb[p6, p5, p7, p6])*
                  (4*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p2, p5] + 
                   2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5] + 
                   (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, p1, 
                     p3, p2] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p5], p1, p2])/
                    Spbb[p2, p1, p3, p2] - (2*MP[p6, p7]*((2*MP[p2, p3]*
                        Spbb[p2, p6]*Spbb[p2, p1, p5, p6])/Spbb[p2, p1, p3, 
                        p2] + Spbb[p6, p4, p5, p6] + Spbb[p6, SumM[p2, p3], 
                       p5, p6]))/Spbb[p6, p5, p7, p6]))/((2*MP[p1, p1] + 
                    2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
                   ((-2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p5, p4])/
                     Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], p5, p4] - 
                    (2*MP[p6, p7]*Spbb[p4, p6]*((2*MP[p2, p3]*Spbb[p2, p4]*
                         Spbb[p2, p1, p5, p6])/Spbb[p2, p1, p3, p2] + 
                       Spbb[p4, SumM[p2, p3], p5, p6]))/Spbb[p6, p5, p7, 
                      p6]) + ((-2*MP[p2, p3]*Spbb[p2, p4]*(MP[p1, p1]*
                        Spbb[p2, p4] + Spbb[p2, p1, SumM[p2, p3], p4]))/
                     Spbb[p2, p1, p3, p2] - Spbb[p4, SumM[p2, p3], p1, p4])*
                   (2*MP[p4, p5] - (2*MP[p6, p7]*Spbb[p6, p4, p5, p6])/
                     Spbb[p6, p5, p7, p6]))))/((2*MP[p2, p3] + 2*MP[p2, p4] + 
                2*MP[p3, p4] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
                 Spbb[p2, p1, p3, p2])*(2*MP[p4, p5] - (2*MP[p6, p7]*
                  Spbb[p6, p4, p5, p6])/Spbb[p6, p5, p7, p6]))))/
           (MP[p5, p6]*MP[p6, p7]) + (MP[p5, p5]*Spbb[p4, p6]*
           (1 + (4*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p2, p5] + 2*MP[p3, p4] + 
              2*MP[p3, p5] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
               Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p5], 
                 p1, p2])/Spbb[p2, p1, p3, p2] - (2*MP[p4, p5]*
                ((2*MP[p2, p3]*Spbb[p2, p6]*Spbb[p2, p1, p5, p6])/
                  Spbb[p2, p1, p3, p2] + Spbb[p6, SumM[p2, p3], p5, p6]))/
               Spbb[p6, p4, p5, p6] - ((2*MP[p6, p7] - (2*MP[p4, p5]*
                   Spbb[p6, p5, p7, p6])/Spbb[p6, p4, p5, p6])*
                ((2*MP[p2, p3]*Spbb[p2, p6]*Spbb[p2, p1, SumM[p4, p5], p6])/
                  Spbb[p2, p1, p3, p2] + Spbb[p6, SumM[p2, p3], SumM[p4, p5], 
                  p6]))/Spbb[p6, SumM[p4, p5], p7, p6])/(2*MP[p1, p1] + 
              2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3]))*
           Spbb[p6, SumM[p4, p5], p7, p6])/(2*MP[p4, p5]*(2*MP[p4, p5] + 
            2*MP[p4, p6] + 2*MP[p5, p6])*Spab[p4, p5, p6]*
           (2*MP[p6, p7] - (2*MP[p4, p5]*Spbb[p6, p5, p7, p6])/
             Spbb[p6, p4, p5, p6]))))/(MP[p2, p3]*Spab[p2, p1, p2]) + 
     (MP[p1, p1]*Spbb[p2, p4]^2*((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
         (2*MP[p6, p7]*Spbb[p2, p6]*Spbb[p2, p1, p5, p6] - 
          Spbb[p2, SumM[p3, p4, p5], p1, p2]*Spbb[p6, p5, p7, p6]) + 
        Spbb[p2, SumM[p3, p4], p1, p2]*((2*MP[p1, p6] + 2*MP[p1, p7] + 
            2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 4*MP[p6, p7])*
           Spbb[p6, p5, p7, p6] - 2*MP[p6, p7]*Spbb[p6, SumM[p2, p3, p4], p5, 
            p6])))/(8*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
       MP[p5, p6]*MP[p6, p7]*(2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 
        2*MP[p6, p7])*((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
         Spbb[p2, p1, p3, p2] + 2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2]))
 
ATreeSAM["S", "+", "S", "+", "s", "s", "+"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/2*(MP[p1, p1]*Spbb[p2, p7]*
        (-1/2*(MP[p1, p1]*Spbb[p2, p4]*Spbb[p2, SumM[p3, p4], SumM[p1, p7], 
             p2]*(1 + (2*MP[p3, p5] + 2*MP[p4, p5] + (2*MP[p3, p4]*
                 Spbb[p2, p3, p5, p2])/Spbb[p2, p4, p3, p2] + 
               (Spbb[p2, SumM[p3, p4], p5, p2]*(2*MP[p1, p2] + 2*MP[p1, p7] + 
                  2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p2, SumM[p4, p5, p6], p3, 
                     p2])/Spbb[p2, p4, p3, p2]))/Spbb[p2, SumM[p3, p4], 
                 SumM[p1, p7], p2])/(2*MP[p5, p5] + 2*MP[p5, p6])))/
           (MP[p3, p4]*Spab[p4, p3, p2]*(2*MP[p2, p3] + 2*MP[p2, p4] + 
             2*MP[p3, p4] + (2*MP[p1, p7]*Spbb[p2, SumM[p3, p4], p1, p2])/
              Spbb[p2, SumM[p1, p7], p1, p2])*(2*MP[p1, p2] + 2*MP[p1, p7] + 
             2*MP[p2, p7] - (2*MP[p3, p4]*Spbb[p2, SumM[p4, p5, p6], p3, p2])/
              Spbb[p2, p4, p3, p2])) - (Spbb[p2, SumM[p4, p5, p6], p3, p2]*
           (-1/2*((-(((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
                   Spbb[p2, p4]*Spbb[p2, p3, p5, p4])/Spbb[p2, SumM[p4, p5, 
                    p6], p3, p2]) + Spbb[p4, p3, p5, p4])*(1 + 
                ((2*(MP[p3, p4] + MP[p3, p5] + MP[p4, p5]) - 
                   ((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, 
                      SumM[p4, p5], p3, p2])/Spbb[p2, SumM[p4, p5, p6], p3, 
                     p2])*(-(((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*
                      Spbb[p2, p4]*Spbb[p2, p3, p5, p4])/Spbb[p2, SumM[p4, 
                       p5, p6], p3, p2]) + Spbb[p4, p3, p5, p4]))/
                 ((2*MP[p5, p5] + 2*MP[p5, p6])*(-(((2*MP[p1, p2] + 
                        2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4]*Spbb[p2, 
                        p3, p5, p4])/Spbb[p2, SumM[p4, p5, p6], p3, p2]) + 
                    Spbb[p4, p3, p5, p4]) - (2*MP[p3, p4] - ((2*MP[p1, p2] + 
                       2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4, p3, p2])/
                     Spbb[p2, SumM[p4, p5, p6], p3, p2])*Spbb[p4, p5, p6, 
                    p4])))/(MP[p4, p5]*(2*MP[p3, p4] - ((2*MP[p1, p2] + 
                   2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4, p3, p2])/
                 Spbb[p2, SumM[p4, p5, p6], p3, p2])) + 
            (MP[p1, p1]*Spbb[p4, p5, p6, p4]^2 + MP[p5, p5]*
               (((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4]*
                   (-(MP[p1, p1]*Spbb[p2, p4]) - Spbb[p2, p3, p2, p4] - 
                    Spbb[p2, p3, SumM[p1, p7], p4]))/Spbb[p2, SumM[p4, p5, 
                    p6], p3, p2] - Spbb[p4, p2, p3, p4] + Spbb[p4, 
                  SumM[p2, p5, p6], p3, p4])^2)/((2*MP[p5, p5] + 2*
                MP[p5, p6])*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*
                MP[p5, p6])*((2*MP[p5, p5] + 2*MP[p5, p6])*
                (-(((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4]*
                    Spbb[p2, p3, p5, p4])/Spbb[p2, SumM[p4, p5, p6], p3, 
                    p2]) + Spbb[p4, p3, p5, p4]) - (2*MP[p3, p4] - 
                 ((2*MP[p1, p2] + 2*MP[p1, p7] + 2*MP[p2, p7])*Spbb[p2, p4, 
                    p3, p2])/Spbb[p2, SumM[p4, p5, p6], p3, p2])*
                Spbb[p4, p5, p6, p4]))))/((2*MP[p1, p2] + 2*MP[p1, p7] + 
            2*MP[p2, p7])*(2*MP[p2, p3] - (2*MP[p1, p7]*Spbb[p2, p1, p3, p2])/
             Spbb[p2, SumM[p1, p7], p1, p2]))))/(MP[p1, p7]*
        Spab[p7, p1, p2]) + (MP[p1, p1]*Spbb[p2, p4]^2*
       (-(((-2*MP[p1, p5]*(2*MP[p5, p5] + 2*MP[p5, p6]) + 
             (-2*MP[p1, p5] - 2*MP[p1, p7] - 2*MP[p5, p7] + 2*MP[p6, p7])*
              (2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7]) + 
             2*MP[p1, p7]*(2*MP[p5, p5] + 2*MP[p5, p6] + 4*MP[p6, p7]))*
            Spbb[p2, SumM[p3, p4], p1, p2] + (2*MP[p2, p3] + 2*MP[p2, p4] + 
             2*MP[p3, p4])*((2*MP[p5, p5] + 2*MP[p5, p6])*Spbb[p2, p1, p5, 
               p2] + (2*MP[p5, p5] + 2*MP[p5, p6] + 4*MP[p6, p7])*
              Spbb[p2, p7, p1, p2] + (2*MP[p5, p5] + 2*MP[p5, p6] + 2*
                MP[p5, p7] + 2*MP[p6, p7])*Spbb[p2, SumM[p3, p4, p6], p1, 
               p2]))*((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
            Spbb[p2, p7]*Spbb[p2, p1, p6, p7] + Spbb[p2, SumM[p3, p4], p1, 
             p2]*Spbb[p7, p1, p6, p7])) + 
        (-((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p1, p6, 
             p2]) + 2*MP[p1, p6]*Spbb[p2, SumM[p3, p4], p1, p2])*
         ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
           (2*MP[p6, p7]*Spbb[p2, p7]*(MP[p1, p1]*Spbb[p2, p7] + 
              Spbb[p2, p1, SumM[p2, p3, p4], p7]) + Spbb[p2, p7, p1, p2]*
             Spbb[p7, p5, p6, p7]) + Spbb[p2, SumM[p3, p4], p1, p2]*
           (2*MP[p1, p7]*Spbb[p7, p5, p6, p7] - 2*MP[p6, p7]*
             Spbb[p7, SumM[p2, p3, p4], p1, p7]))))/
      (8*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
       (2*MP[p5, p5] + 2*MP[p5, p6])*MP[p6, p7]*(2*MP[p5, p5] + 
        2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*(2*MP[p2, p3] + 
        ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p1, p3, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])*(2*MP[p1, p7] + 
        ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p7, p1, p2])/
         Spbb[p2, SumM[p3, p4], p1, p2])*Spbb[p2, SumM[p3, p4], p1, p2]^2) - 
     (Spbb[p2, p1, p3, p2]*((MP[p1, p1]*(MP[p1, p1]*Spbb[p4, p7] + 
           (2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p7])/
            Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p7]*
             Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + 
           Spbb[p4, SumM[p2, p3], p1, p7])*(Spbb[p7, p1, p5, p7] + 
           Spbb[p7, p4, p5, p7] + Spbb[p7, SumM[p2, p3], p5, p7]))/
         ((2*MP[p5, p5] + 2*MP[p5, p6])*(Spab[p4, p1, p7] + 
           Spab[p4, SumM[p2, p3], p7])*(-2*MP[p5, p5] - 2*MP[p5, p6] + 
           Spab[p7, p1, p7] + Spab[p7, p4, p7] + Spab[p7, SumM[p2, p3], p7])*
          (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] + 
           (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, p1, p3, 
             p2])*(2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/
            Spbb[p2, p1, p3, p2] + ((2*MP[p5, p5] + 2*MP[p5, p6])*
             ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/Spbb[p2, p1, 
                p3, p2] - Spbb[p7, p1, p6, p7]))/(Spbb[p7, p1, p6, p7] + 
             Spbb[p7, p4, p6, p7] + Spbb[p7, SumM[p2, p3], p6, p7]))) + 
        (MP[p5, p5]*(-Spbb[p4, p1, p6, p7] - Spbb[p4, p7, p6, p7] - 
            Spbb[p4, SumM[p2, p3], p6, p7])^2*
          ((-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p7]^2)/Spbb[p2, p1, p3, p2] - 
            (2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, SumM[p2, p3], p7])/
             Spbb[p2, p1, p3, p2] - Spbb[p7, SumM[p2, p3], p1, p7])^2)/
         (2*(2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 
           2*MP[p5, p6])*((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 
             2*MP[p5, p6])*((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
              Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7]) + 
           (2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/
              Spbb[p2, p1, p3, p2])*(Spbb[p7, p1, p6, p7] + 
             Spbb[p7, SumM[p2, p3], p6, p7]))*
          (-((2*MP[p5, p5] + 2*MP[p5, p6])*(Spbb[p7, p1, p6, p7] + 
              Spbb[p7, SumM[p2, p3], p6, p7])) + 
           (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
            (Spbb[p7, p1, p6, p7] + Spbb[p7, p4, p6, p7] + 
             Spbb[p7, SumM[p2, p3], p6, p7]))) + 
        (MP[p1, p1]*(-Spbb[p4, p5, p1, p7] - Spbb[p4, p5, SumM[p2, p3], p7])*
          (-((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
             Spbb[p4, p7]) - Spbb[p4, p6, p1, p7] - Spbb[p4, p6, 
            SumM[p2, p3], p7]))/(2*(2*MP[p1, p1] + 2*MP[p1, p2] + 
           2*MP[p1, p3] + 2*MP[p2, p3])*MP[p4, p5]*(2*MP[p4, p5] + 
           2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*(2*MP[p1, p7] - 
           (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/Spbb[p2, p1, p3, p2] + 
           ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
             ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/Spbb[p2, p1, 
                p3, p2] - Spbb[p7, p1, p6, p7]))/(Spbb[p7, p1, p6, p7] + 
             Spbb[p7, SumM[p2, p3], p6, p7]))*(2*MP[p5, p5] + 2*MP[p5, p6] - 
           ((2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6])*
             (Spbb[p7, p1, p6, p7] + Spbb[p7, p4, p6, p7] + Spbb[p7, SumM[p2, 
                p3], p6, p7]))/(Spbb[p7, p1, p6, p7] + Spbb[p7, SumM[p2, p3], 
              p6, p7]))) - (((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
            Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7])*
          ((2*MP[p2, p3] + 2*MP[p2, p5] + 2*MP[p3, p5] + 
             (2*MP[p2, p3]*Spbb[p2, SumM[p3, p5], p1, p2])/Spbb[p2, p1, p3, 
               p2])*((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] + 
               (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, p1, p3, 
                 p2])*Spbb[p4, p5, p6, p4] + 2*MP[p4, p5]*
              ((-2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]^2)/Spbb[p2, p1, p3, 
                 p2] - (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p7, p4])/
                Spbb[p2, p1, p3, p2] - (2*MP[p2, p3]*Spbb[p2, p4]*
                 Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + Spbb[
                p4, p7, SumM[p2, p3], p4] - Spbb[p4, SumM[p2, p3], p1, p4]) - 
             (Spbb[p4, p7]*(2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, 
                   p2])/Spbb[p2, p1, p3, p2])*(-((2*MP[p2, p3] + 
                   2*MP[p2, p4] + 2*MP[p3, p4] + (2*MP[p2, p3]*Spbb[p2, 
                      SumM[p3, p4], p1, p2])/Spbb[p2, p1, p3, p2])*
                  Spbb[p4, p5, p6, p7]) + 2*MP[p4, p5]*((-2*MP[p2, p3]*
                    Spbb[p2, p4]*Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] - 
                  Spbb[p4, SumM[p2, p3], p6, p7])))/((2*MP[p2, p3]*
                 Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] - 
               Spbb[p7, p1, p6, p7])) - ((2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, 
                p1, p5, p4])/Spbb[p2, p1, p3, p2] + Spbb[p4, SumM[p2, p3], 
              p5, p4])*(-4*MP[p4, p5]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
                Spbb[p2, p1, p3, p2]) - (2*MP[p2, p3] + 2*MP[p2, p4] + 2*
                MP[p3, p4] + (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
                Spbb[p2, p1, p3, p2])*(2*MP[p5, p5] + 2*MP[p5, p6] - 
               ((2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/
                   Spbb[p2, p1, p3, p2])*Spbb[p7, p5, p6, p7])/
                ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
                  Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7])) + 
             (2*MP[p5, p5] + 2*MP[p5, p6] - ((2*MP[p1, p7] - (2*MP[p2, p3]*
                    Spbb[p2, p7, p1, p2])/Spbb[p2, p1, p3, p2])*Spbb[p7, p5, 
                  p6, p7])/((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
                  Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7]))*
              (2*MP[p2, p3] + 2*MP[p2, p6] + 2*MP[p3, p6] + (2*MP[p2, p3]*
                 Spbb[p2, SumM[p3, p6], p1, p2])/Spbb[p2, p1, p3, p2] + 
               ((2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/
                   Spbb[p2, p1, p3, p2])*((-2*MP[p2, p3]*Spbb[p2, p7]*
                    Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] - Spbb[p7, 
                   SumM[p2, p3], p6, p7]))/((2*MP[p2, p3]*Spbb[p2, p7]*
                   Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] - Spbb[p7, p1, 
                  p6, p7])) - 2*MP[p4, p5]*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*
                MP[p5, p5] + 2*MP[p5, p6] + ((2*MP[p1, p7] - (2*MP[p2, p3]*
                    Spbb[p2, p7, p1, p2])/Spbb[p2, p1, p3, p2])*
                 (Spbb[p7, p1, p6, p7] + Spbb[p7, SumM[p2, p3], p6, p7]))/
                ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
                  Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7])) + 
             (2*MP[p1, p5] + 2*MP[p1, p7] + 2*MP[p5, p7] - (2*MP[p2, p3]*
                 Spbb[p2, SumM[p5, p7], p1, p2])/Spbb[p2, p1, p3, p2] + 
               ((2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/
                   Spbb[p2, p1, p3, p2])*((-2*MP[p2, p3]*Spbb[p2, p7]*
                    Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] + Spbb[p7, p1, 
                   p6, p7] + Spbb[p7, p5, p6, p7]))/((2*MP[p2, p3]*
                   Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/Spbb[p2, p1, p3, p2] - 
                 Spbb[p7, p1, p6, p7]))*(2*MP[p4, p5] + 2*MP[p4, p6] + 2*
                MP[p5, p5] + 2*MP[p5, p6] + ((2*MP[p1, p7] - (2*MP[p2, p3]*
                    Spbb[p2, p7, p1, p2])/Spbb[p2, p1, p3, p2])*
                 (Spbb[p7, p1, p6, p7] + Spbb[p7, SumM[p2, p3], p6, p7]))/
                ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
                  Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7])))))/
         (4*MP[p4, p5]*Spab[p7, p6, p7]*(2*MP[p1, p7] - 
           (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/Spbb[p2, p1, p3, p2])*
          (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] + 
           (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, p1, p3, 
             p2])*(2*MP[p5, p5] + 2*MP[p5, p6] - 
           ((2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/Spbb[p2, p1, 
                p3, p2])*Spbb[p7, p5, p6, p7])/
            ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
              Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7]))*
          (2*MP[p4, p5] + 2*MP[p4, p6] + 2*MP[p5, p5] + 2*MP[p5, p6] + 
           ((2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/Spbb[p2, p1, 
                p3, p2])*(Spbb[p7, p1, p6, p7] + Spbb[p7, SumM[p2, p3], p6, 
               p7]))/((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
              Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7]))) - 
        ((-((MP[p1, p1]*Spbb[p4, p7]^4*(Spbb[p7, p1, p6, p7] + Spbb[p7, p4, 
                p6, p7] + Spbb[p7, SumM[p2, p3], p6, p7]))/
             ((2*MP[p1, p1] + 2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
              (-Spbb[p4, p1, p6, p7] - Spbb[p4, p7, p6, p7] - Spbb[p4, 
                SumM[p2, p3], p6, p7]))) - 
           ((-(MP[p1, p1]*Spbb[p4, p7]) - (2*MP[p1, p1]*MP[p2, p3]*
                 Spbb[p2, p4]*Spbb[p2, p7])/Spbb[p2, p1, p3, p2] - 
               (2*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p1, p4, p7])/Spbb[p2, p1, 
                 p3, p2] - (2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, 
                  SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] - Spbb[p4, 
                SumM[p2, p3], p1, p7] - Spbb[p4, SumM[p2, p3], p4, p7])^2*
             ((2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p7]^2)/Spbb[p2, p1, p3, p2] + 
               (2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p4, p7])/Spbb[p2, p1, 
                 p3, p2] + (2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, 
                  SumM[p2, p3], p7])/Spbb[p2, p1, p3, p2] + Spbb[p7, p4, p1, 
                p7] + Spbb[p7, SumM[p2, p3], p1, p7])^2)/
            ((Spab[p4, p1, p7] + Spab[p4, SumM[p2, p3], p7])*
             (-2*MP[p5, p5] - 2*MP[p5, p6] + Spab[p7, p1, p7] + 
              Spab[p7, p4, p7] + Spab[p7, SumM[p2, p3], p7])*
             (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] + (2*MP[p2, p3]*
                Spbb[p2, SumM[p3, p4], p1, p2])/Spbb[p2, p1, p3, p2])*
             (2*MP[p1, p7] - (2*MP[p2, p3]*Spbb[p2, p7, p1, p2])/Spbb[p2, p1, 
                p3, p2] + ((2*MP[p5, p5] + 2*MP[p5, p6])*
                ((2*MP[p2, p3]*Spbb[p2, p7]*Spbb[p2, p1, p6, p7])/
                  Spbb[p2, p1, p3, p2] - Spbb[p7, p1, p6, p7]))/(
                Spbb[p7, p1, p6, p7] + Spbb[p7, p4, p6, p7] + Spbb[p7, 
                 SumM[p2, p3], p6, p7]))))*Spbb[p7, p6, SumM[p5, p6], p5, p6, 
           p7])/((2*MP[p5, p5] + 2*MP[p5, p6])*(MP[p1, p1]*Spbb[p4, p7] + 
           (2*MP[p1, p1]*MP[p2, p3]*Spbb[p2, p4]*Spbb[p2, p7])/
            Spbb[p2, p1, p3, p2] + (2*MP[p2, p3]*Spbb[p2, p7]*
             Spbb[p2, p1, SumM[p2, p3], p4])/Spbb[p2, p1, p3, p2] + 
           Spbb[p4, SumM[p2, p3], p1, p7])*(Spbb[p7, p1, p6, p7] + 
            Spbb[p7, p4, p6, p7] + Spbb[p7, SumM[p2, p3], p6, p7])^2)))/
      (2*MP[p2, p3]*Spab[p2, p1, p2])
 
ATreeSAM["S", "+", "+", "+", "s", "S", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p5, p5]*((2*MP[p1, p1] + 2*MP[p1, p6] + 2*MP[p1, p7] + 
          2*MP[p6, p7])*Spbb[p2, p4] - Spbb[p2, SumM[p3, p4, p5], p2, p4] - 
        Spbb[p2, SumM[p3, p4, p5], p3, p4]))/(2*MP[p4, p5]*
       (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*(2*MP[p1, p1] + 
        2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*Spaa[p3, p4]*
       (Spaa[p2, p3] + ((2*MP[p1, p1] + 2*MP[p1, p6] + 2*MP[p1, p7] + 
           2*MP[p6, p7])*Spab[p3, p1, p2])/Spbb[p2, SumM[p3, p4, p5], p1, 
          p2])) + (Spbb[p2, p3]^3*
       ((MP[p5, p5]*Spbb[p2, p1, SumM[p2, p3], p4]^2)/(2*MP[p4, p5]*
          (2*MP[p1, p1] + 2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7] - 
           (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4, p5], p1, p2])/
            Spbb[p2, p3, p1, p2])*(2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4] - (2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])/
            Spbb[p2, p3, p1, p2] - ((2*MP[p1, p1] + 2*MP[p1, p6] + 
              2*MP[p1, p7] + 2*MP[p6, p7] - (2*MP[p2, p3]*Spbb[p2, 
                 SumM[p3, p4, p5], p1, p2])/Spbb[p2, p3, p1, p2])*
             Spbb[p2, p1, SumM[p2, p3], p4, SumM[p1, p2, p3], SumM[p2, p3], 
              p1, p2])/Spbb[p2, p1, SumM[p2, p3], SumM[p4, p5], 
             SumM[p1, p2, p3], SumM[p2, p3], p1, p2])) + 
        (Spbb[p2, p1, SumM[p2, p3], p4, SumM[p1, p2, p3], SumM[p2, p3], p1, 
           p2]*(4*MP[p2, p3]^2*Spbb[p2, SumM[p6, p7], p5, p2]*
            Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p1, p2]^2 - 
           2*MP[p2, p3]*Spbb[p2, p3, p1, p2]*Spbb[p2, p1, SumM[p2, p3], 
             SumM[p1, p2, p3], p1, p2]*(Spbb[p2, p1, SumM[p2, p3], 
              SumM[p1, p2, p3], p4, SumM[p6, p7], p5, p2] + 
             Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], 
              SumM[p6, p7], p5, p2] + Spbb[p2, SumM[p6, p7], p5, p4, 
              SumM[p1, p2, p3], SumM[p2, p3], p1, p2] + 
             Spbb[p2, SumM[p6, p7], p5, SumM[p2, p3], SumM[p1, p2, p3], 
              SumM[p2, p3], p1, p2]) + Spbb[p2, p3, p1, p2]^2*
            (-Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, p5, SumM[p6, 
                p7], SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], p1, p2] + 
             Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, SumM[p6, p7], 
              p5, p4, SumM[p1, p2, p3], SumM[p2, p3], p1, p2] + 
             Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p4, SumM[p6, p7], 
              p5, SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], p1, p2] + 
             Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], 
              SumM[p6, p7], p5, SumM[p2, p3], SumM[p1, p2, p3], SumM[p2, p3], 
              p1, p2])))/((2*MP[p1, p2] + 2*MP[p1, p3] + 2*MP[p2, p3])*
          (2*MP[p5, p5] + 2*MP[p5, p6] + 2*MP[p5, p7] + 2*MP[p6, p7])*
          (-2*MP[p2, p3]*Spab[p4, p1, p2]*(Spab[p4, p1, p2] + 
             Spab[p4, p3, p2]) + Spaa[p4, SumM[p2, p3], p1, p4]*
            Spbb[p2, p3, p1, p2])*(-2*MP[p2, p3]*Spbb[p2, SumM[p3, p5], p1, 
             p2]*Spbb[p2, p1, SumM[p2, p3], p4, SumM[p1, p2, p3], 
             SumM[p2, p3], p1, p2] - 2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, 
             p2]*Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p5, 
             SumM[p2, p3], p1, p2] + Spbb[p2, p3, p1, p2]*
            ((2*MP[p2, p3] + 2*MP[p2, p5] + 2*MP[p3, p5] + 2*MP[p4, p5])*
              Spbb[p2, p1, SumM[p2, p3], p4, SumM[p1, p2, p3], SumM[p2, p3], 
               p1, p2] + (2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
              Spbb[p2, p1, SumM[p2, p3], SumM[p1, p2, p3], p5, SumM[p2, p3], 
               p1, p2])))))/(2*MP[p2, p3]*Spbb[p2, p1, p2, p3]*
       Spbb[p2, p3, p1, p2])
 
ATreeSAM["S", "+", "+", "s", "+", "S", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = (MP[p4, p4]*Spbb[p2, p3]^2*Spbb[p5, SumM[p1, p7], p6, p5])/
      (2*MP[p2, p3]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
       (2*MP[p1, p1] + 2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*
       (2*MP[p3, p4] + ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
          Spbb[p5, p3, p4, p5])/Spbb[p5, SumM[p1, p6, p7], p4, p5])*
       (2*MP[p5, p6] + ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
          Spbb[p5, p4, p6, p5])/Spbb[p5, SumM[p1, p6, p7], p4, p5])) + 
     (MP[p4, p4]*Spbb[p3, p5]*(((-Spbb[p2, SumM[p5, p6, p7], p1, p2] + 
           (2*MP[p3, p4]*Spbb[p2, p5]*Spbb[p5, p4, p1, p2])/
            Spbb[p5, SumM[p3, p4], p4, p5] - 
           (Spbb[p2, p5]*Spbb[p5, SumM[p3, p4], p1, p2]*(2*MP[p5, p6] - 
              (2*MP[p3, p4]*Spbb[p5, p4, p6, p5])/Spbb[p5, SumM[p3, p4], p4, 
                p5]))/Spbb[p5, SumM[p1, p2, p7], p6, p5])*
          Spbb[p5, SumM[p1, p2, p7], p6, p5])/(2*MP[p1, p2]*
          (2*MP[p3, p4] + 2*MP[p3, p5] + 2*MP[p4, p5])*(2*MP[p5, p6] - 
           (2*MP[p3, p4]*Spbb[p5, p4, p6, p5])/Spbb[p5, SumM[p3, p4], p4, 
             p5])*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4] - 
           (Spbb[p5, p2, SumM[p1, p6, p7], p5]*(2*MP[p5, p6] - 
              (2*MP[p3, p4]*Spbb[p5, p4, p6, p5])/Spbb[p5, SumM[p3, p4], p4, 
                p5]))/Spbb[p5, SumM[p1, p2, p7], p6, p5] + 
           (2*MP[p3, p4]*Spbb[p5, SumM[p1, p6, p7], p4, p5])/
            Spbb[p5, SumM[p3, p4], p4, p5])) - 
        (MP[p4, p4]*Spbb[p2, p5]*(-Spbb[p5, p2, p6, p5] + 
           Spbb[p5, SumM[p1, p2, p7], p6, p5]))/
         ((2*MP[p1, p1] + 2*MP[p1, p6] + 2*MP[p1, p7] + 2*MP[p6, p7])*
          Spab[p2, SumM[p1, p6, p7], p5]*(2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4] + (2*MP[p3, p4]*Spbb[p5, SumM[p1, p6, p7], p4, p5])/
            Spbb[p5, SumM[p3, p4], p4, p5])*(2*MP[p5, p6] - 
           (2*MP[p3, p4]*Spbb[p5, p4, p6, p5])/Spbb[p5, SumM[p3, p4], p4, 
             p5] - (Spbb[p5, SumM[p1, p2, p7], p6, p5]*(2*MP[p2, p3] + 
              2*MP[p2, p4] + 2*MP[p3, p4] + (2*MP[p3, p4]*Spbb[p5, 
                 SumM[p1, p6, p7], p4, p5])/Spbb[p5, SumM[p3, p4], p4, p5]))/
            Spbb[p5, p2, SumM[p1, p6, p7], p5]))))/
      (2*MP[p3, p4]*Spab[p3, p4, p5]) - 
     (Spbb[p5, p4, p6, p5]*((MP[p4, p4]*Spbb[p2, p3]^2)/
         ((2*MP[p3, p4] - (2*MP[p5, p6]*Spbb[p5, p3, p4, p5])/
            Spbb[p5, p4, p6, p5])*(2*MP[p2, p3] + 2*MP[p2, p4] + 
           2*MP[p3, p4] - (2*MP[p5, p6]*Spbb[p5, SumM[p2, p3], p4, p5])/
            Spbb[p5, p4, p6, p5])*(2*MP[p2, p3] + (Spbb[p2, p3, p1, p2]*
             (-((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p5, p4, p6, 
                 p5]) + 2*MP[p5, p6]*Spbb[p5, SumM[p2, p3], p4, p5]))/
            (2*MP[p5, p6]*Spbb[p2, p5]*Spbb[p5, p4, p1, p2] + 
             (Spbb[p2, p3, p1, p2] + Spbb[p2, p4, p1, p2])*Spbb[p5, p4, p6, 
               p5]))) + (Spbb[p2, p3, p1, p2]*
          (2*MP[p5, p6]*Spbb[p5, SumM[p2, p3], p1, p2]*
            (MP[p4, p4]*Spbb[p5, SumM[p2, p3], p1, p2] + Spbb[p5, p4, p7, 
              SumM[p2, p3], p1, p2] + Spbb[p5, p4, SumM[p5, p6], 
              SumM[p2, p3], p1, p2]) + Spbb[p5, p4, p6, p5]*
            (-Spbb[p2, p1, SumM[p2, p3], p4, p7, SumM[p2, p3], p1, p2] + 
             Spbb[p2, p1, SumM[p2, p3], SumM[p5, p6], p4, SumM[p2, p3], p1, 
              p2])))/(2*MP[p1, p2]*(2*MP[p1, p2] + 2*MP[p1, p3] + 
           2*MP[p2, p3])*Spaa[p3, p2, p1, p3]*
          (-2*MP[p5, p6]*(Spbb[p2, p3, p1, p2]*Spbb[p5, p2, p4, p5] + 
             Spbb[p2, p3, p1, p2]*Spbb[p5, p3, p4, p5] + 2*MP[p2, p3]*
              Spbb[p2, p5]*Spbb[p5, p4, p1, p2]) + 
           ((2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p3, p1, p2] - 
             2*MP[p2, p3]*Spbb[p2, p4, p1, p2])*Spbb[p5, p4, p6, p5]))))/
      (2*MP[p5, p6]*Spab[p5, p4, p5])
 
ATreeSAM["S", "+", "+", "s", "S", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/2*(Spbb[p6, p5, p7, p6]*((MP[p7, p7]*Spbb[p2, p3]^2)/
         (2*MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
          (2*MP[p2, p3] - ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*
             Spbb[p2, p3, p1, p2])/Spbb[p2, SumM[p3, p4], p1, p2])) + 
        (Spbb[p2, p3, p1, p2]*(Spbb[p2, p1, SumM[p2, p3], p5, p4, 
            SumM[p2, p3], p1, p2] + Spbb[p2, p1, SumM[p2, p3], SumM[p6, p7], 
            p4, SumM[p2, p3], p1, p2]))/(2*MP[p1, p2]*(2*MP[p1, p2] + 
           2*MP[p1, p3] + 2*MP[p2, p3])*Spaa[p3, p2, p1, p3]*
          (2*MP[p2, p3]*Spbb[p2, p1, p4, p2] + (2*MP[p2, p4] + 2*MP[p3, p4])*
            Spbb[p2, p3, p1, p2]))))/(MP[p6, p7]*Spab[p6, p5, p6])
 
ATreeSAM["S", "+", "s", "+", "S", "+", "s"][p1_, p2_, p3_, p4_, p5_, p6_, 
     p7_] = -1/8*(MP[p3, p3]*Spbb[p2, p4]^2*Spbb[p6, p5, p7, p6])/
       (MP[p3, p4]*(2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*MP[p5, p6]*
        MP[p6, p7]*(2*MP[p2, p3] + ((2*MP[p2, p3] + 2*MP[p2, p4] + 
            2*MP[p3, p4])*Spbb[p2, p1, p3, p2])/Spbb[p2, SumM[p3, p4], p1, 
           p2])) - (Spbb[p2, p1, p3, p2]*
       (-1/4*(Spbb[p6, p5, p7, p6]*(2*MP[p2, p3]*Spbb[p2, p4]*
             (2*MP[p6, p7]*Spbb[p4, p6]*Spbb[p2, p1, p5, p6] + 
              Spbb[p2, p1, p5, p4]*Spbb[p6, p5, p7, p6]) + 
            Spbb[p2, p1, p3, p2]*(2*MP[p6, p7]*Spbb[p4, p6]*Spbb[p4, 
                SumM[p2, p3], p5, p6] + Spbb[p4, SumM[p2, p3], p5, p4]*Spbb[
                p6, p5, p7, p6])))/(MP[p5, p6]*MP[p6, p7]*
           ((2*MP[p2, p3] + 2*MP[p2, p4] + 2*MP[p3, p4])*Spbb[p2, p1, p3, 
              p2] + 2*MP[p2, p3]*Spbb[p2, SumM[p3, p4], p1, p2])*
           (2*MP[p6, p7]*Spbb[p6, p4, p5, p6] - 2*MP[p4, p5]*
             Spbb[p6, p5, p7, p6])) - (MP[p1, p1]*Spbb[p4, p6]*
          Spbb[p6, SumM[p4, p5], p7, p6])/(2*MP[p4, p5]*(2*MP[p4, p5] + 
           2*MP[p4, p6] + 2*MP[p5, p6])*Spab[p4, p5, p6]*
          (2*MP[p6, p7] - (2*MP[p4, p5]*Spbb[p6, p5, p7, p6])/
            Spbb[p6, p4, p5, p6]))))/(2*MP[p2, p3]*Spab[p2, p1, p2])
