#!/bin/sh

SLHC_version=6_2_0_SLHC25_patch2

source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh

cmsenv

cd /afs/cern.ch/work/a/amagnan/CMSSW_$SLHC_version/src/

eval `scramv1 runtime -sh`

cd -

cp /afs/cern.ch/work/a/amagnan/CMSSW_$SLHC_version/src/UserCode/HGCanalysis/test/hgcPhotonRecoConfig.py .
cp -r /afs/cern.ch/work/a/amagnan/CMSSW_$SLHC_version/src/UserCode/HGCanalysis/test/sample .

cmsRun hgcPhotonRecoConfig.py | tee hgcPhotonReco.out
cp Calib_singleGamma_0pu.root /afs/cern.ch/work/a/amagnan/CMSSW_$SLHC_version/src/UserCode/HGCanalysis/test/
