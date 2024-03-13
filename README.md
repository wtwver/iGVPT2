# TODO
faster parallel water

understand header basis
Draw spectrum

# Requirements
downlaod ORCA from https://orcaforum.kofo.mpg.de/app.php/dlext/
```
sudo apt install -y openbabel
export LD_LIBRARY_PATH=/home/debian/orca 
```

## How to run example
Geometry in 2MR.ici can be `water.out` or `h2oB2PLYP.out`

```console
cd examples/Ex3-*
/home/debian/orca/orca water.inp
cat water.hess >> water.out
./xHybrid 2MR.ici && tail -n 30 2MR.out
```

## How to Compile from source :\
Requirements : cchemilib ([see github](https://github.com/allouchear/cchemi))

```console
cd src
./cleanigvpt2.sh
./compigvpt2.sh
```

## Methods

The description of the methods implemented in iGVPT2 is given in this document : https://arxiv.org/abs/1704.02144

## Citations

If you used Hybrid method in iGVPT2, please cite :

    Fast and Accurate Hybrid QM//MM Approach for Computing, Anharmonic Corrections to Vibrational Frequencies, Loïc Barnes, Baptiste Schindler, Isabelle Compagnon and Abdul-Rahman Allouche, Journal of Molecular Modeling 22, 285 (2016). https://doi.org/10.1007/s00894-016-3135-5

If you used Neural Networks potential implemented in iGVPT2, please cite :

    Combining quantum mechanics and machine-learning calculations for anharmonic corrections to vibrational frequencies. J Lam, S Abdul-Al, AR Allouche.  J. Chem. Theory Comput. 2020, 1681-1689.https://doi.org/10.1021/acs.jctc.9b00964.
    Parallel Multistream Training of High-Dimensional Neural Network Potentials. Singraber, A.; Morawietz, T.; Behler, J.; Dellago, C. J. Chem. Theory Comput. 2019, 15 (5), 3075–3092. https://doi.org/10.1021/acs.jctc.8b01092.

## Contributors

The code is written by Abdul-Rahman Allouche.\
The manual is writen by Loïc Barnes.\
Testing and Debugging : Loïc Barnes, Baptiste Schindler, Isabelle Compagnon, Abdul-Rahman Allouche
    
## User register

PLEASE register as a iGVPT2 User. This will help to keep up the support for iGVPT2. There is no commitment involved whatsoever.  [Click here to register.](https://docs.google.com/forms/d/e/1FAIpQLSc1wHwn9g3JN2rvYW0dNS-7Xf3dWE3WJ75jP3C6YP1aBgUEeQ/viewform)

 ## License

This software is licensed under the [GNU General Public License version 3 or any later version (GPL-3.0-or-later)](https://www.gnu.org/licenses/gpl.txt).


 