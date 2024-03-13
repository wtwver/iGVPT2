! TightOpt NumFreq 
! RIJCOSX B3LYP/G D3BJ def2-tzvp def2-tzvp/C TIGHTSCF  
%geom MaxIter 800 end
* xyz 0 1
C         0.000000    0.555291    0.000000
F         1.102715   -0.236609    0.000000
H         0.000000    1.154409   -0.907887
H         0.000000    1.154409    0.907887
F        -1.102715   -0.236609    0.000000
*
