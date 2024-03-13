memory,64,m
expec,dm
geomtyp=xyz
geometry={
nosym,
3 ! number of atoms 
GeomXYZ
H         0.748790   -0.462880    0.000000
H        -0.748790   -0.462880    0.000000
O         0.000000    0.115720    0.000000

}

basis=avtz
df-hf,basis=avtz
{df-lmp2,basis=avtz;polari,dm}
