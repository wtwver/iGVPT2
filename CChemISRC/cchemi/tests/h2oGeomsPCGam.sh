mkdir /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomstmp1
cd /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomstmp1
cp /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomsFF_1.inp input
firefly -p -o /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomsFF_1.log
cd ..
mv PUNCH  /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomsFF_1.pun
/bin/rm -r  /home/allouche/CChemI/CChemI-100415/cchemi/tests/h2oGeomstmp1
