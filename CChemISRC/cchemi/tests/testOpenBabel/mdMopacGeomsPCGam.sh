mkdir /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp1
cd /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp1
cp /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_1.inp input
firefly -p -o /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_1.log
cd ..
mv PUNCH  /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_1.pun
/bin/rm -r  /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp1
mkdir /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp2
cd /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp2
cp /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_2.inp input
firefly -p -o /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_2.log
cd ..
mv PUNCH  /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomsFF_2.pun
/bin/rm -r  /home_nfs/ilm/allouche/tmp/CChemI-191215/cchemi/tests/testOpenBabel/mdMopacGeomstmp2
