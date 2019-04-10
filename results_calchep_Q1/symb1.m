(*
     ==============================
     *  CalcHEP  3.7.5 *
     ==============================
  process  e(p1)+E(p2)->t(p3)+T(p4)
*)

parameters={
 EE -> 3.13330000000*10^(-1)
,SW -> 4.74000000000*10^(-1)
,MW -> 8.03850000000*10^(1)
,Mtp -> 1.72500000000*10^(2)
,wZ -> 0.00000000000*10^(0)
           };

substitutions={
 MZ->MW/CW
,CW->Sqrt[1-SW^2]
              };

inParticles = {"e","E"}
outParticles = {"t","T"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =0^2;
p2/: SC[p2,p2] =0^2;
p3/: SC[p3,p3] =Mtp^2;
p2/: SC[p2,p3] = -1*(Mtp^2-0^2-0^2-Mtp^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((32*EE^4)/(3));
numerator =(2*SC[p1,p3]^2-2*SC[p1,p3]*SC[p1,p2]+SC[p1,p2]^2+SC[p1,p2]*Mtp^2
 );
denominator =(propDen[-p1-p2,0,0]^2);

addToSum;

(*
  Diagram  2 in subprocess 1
*)
totFactor = ((2*EE^4)/(3*CW^2*SW^2));
numerator =(64*SC[p1,p3]^2*SW^4-40*SC[p1,p3]^2*SW^2+6*SC[p1,p3]^2-64*
 SC[p1,p3]*SC[p1,p2]*SW^4+40*SC[p1,p3]*SC[p1,p2]*SW^2-12*SC[p1,p3]*SC[p1,p2]
 +32*SC[p1,p2]^2*SW^4-20*SC[p1,p2]^2*SW^2+6*SC[p1,p2]^2+32*SC[p1,p2]*Mtp^2*
 SW^4-20*SC[p1,p2]*Mtp^2*SW^2+3*SC[p1,p2]*Mtp^2);
denominator =(propDen[-p1-p2,0,0]*propDen[-p1-p2,MZ,wZ]);

addToSum;

(*
  Diagram  3 in subprocess 1
*)
totFactor = ((EE^4)/(12*CW^4*SW^4));
numerator =(256*SC[p1,p3]^2*SW^8-320*SC[p1,p3]^2*SW^6+200*SC[p1,p3]^2*SW^4-
 60*SC[p1,p3]^2*SW^2+9*SC[p1,p3]^2-256*SC[p1,p3]*SC[p1,p2]*SW^8+320*
 SC[p1,p3]*SC[p1,p2]*SW^6-296*SC[p1,p3]*SC[p1,p2]*SW^4+120*SC[p1,p3]*
 SC[p1,p2]*SW^2-18*SC[p1,p3]*SC[p1,p2]+128*SC[p1,p2]^2*SW^8-160*SC[p1,p2]^2*
 SW^6+148*SC[p1,p2]^2*SW^4-60*SC[p1,p2]^2*SW^2+9*SC[p1,p2]^2+128*SC[p1,p2]*
 Mtp^2*SW^8-160*SC[p1,p2]*Mtp^2*SW^6+64*SC[p1,p2]*Mtp^2*SW^4-12*SC[p1,p2]*
 Mtp^2*SW^2);
denominator =(propDen[-p1-p2,MZ,wZ]^2);

addToSum;

finishSum;
