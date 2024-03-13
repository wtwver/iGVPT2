%Chk=h2oB2PLYP
#P B2PLYP/DEF2TZVPP SCF(XQC,Tight) Opt(Tight)
#Units(Ang,Deg)

Comments

0   1
  H      0.594398   -0.430112   -0.202093
  H      0.130350    0.936362   -0.671422
  O      0.575252    0.493751    0.053280

--Link1--
%Chk=h2oB2PLYP
#B3LYP/chkbas SCF(XQC,Tight)  Freq iop(7/33=1)
#Test Guess(Read) Geom=AllCheck

