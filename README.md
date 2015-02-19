# HGCanalysis scripts

## For use of filter only for quickie PandoraPFA studies

cd ${CMSSW_BASE}/src
git clone https://github.com/sethzenz/HGCanalysis.git UserCode/HGCanalysis
git checkout origin/hacked-interactions-filter
cd Usercode ; scram b -j 9

## Installation 

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

## GEN-SIM-RECO production

Use the generateEventsFromCfi.sh to steer the generation.
It will call cmsDriver.py and customize the configuration file.
The options can be inspected by calling:

generateEventsFromCfi.sh -h

A full production can be ran locally or submitted to the batch using 
the submitLocalHGCalProduction.py wrapper script. Two examples are given below:

### Particle gun 

For muons it's enough a single energy

#default geometry
python scripts/submitLocalHGCalProduction.py -q 1nd -n 10 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single13_${CMSSW_VERSION} -p 13 -n 100 -e 100";
#change geometry scenario
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single13_v4_${CMSSW_VERSION} -p 13 -n 100 -e 100 -g Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco";

For regression use flat gun

python scripts/submitLocalHGCalProduction.py -q 2nw -n 2500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/FlatPtYSingle22_${CMSSW_VERSION}/RECO_a -c UserCode/HGCanalysis/python/particlePtYGun_cfi.py -n 200 -p 22 -f";
python scripts/submitLocalHGCalProduction.py -q 1nw -n 2500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/FlatPtYSingle22_${CMSSW_VERSION}/RECO_b -c UserCode/HGCanalysis/python/particlePtYGun_cfi.py -n 200 -p 22 -f";

# DECEMBER JAMBOREE PRODUCTION

## SIM

energies=(10 20 40 50 75 100 125 175 250 400 500)
pids=(22 211 11)
pu=(140 0 200)
for p in ${pu[@]}; do
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do
      python scripts/submitLocalHGCalProduction.py -q 1nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/ -p ${pid} -n 25 -e ${en} -s -a ${p}";
     done
done
done

pids=(111 15)
for p in ${pu[@]}; do
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do 
    	python scripts/submitLocalHGCalProduction.py -n 500 -q 2nd -s generateEventsFromCfi.sh -o "-c UserCode/HGCanalysis/python/jetGun_cfi.py -o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}/ -p ${pid} -n 25 -e ${en} -s -a ${p}"; 
    done
done
done

for p in ${pu[@]}; do
    python scripts/submitLocalHGCalProduction.py -q 2nw -n 1000 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/VBFtoH125toTauTau_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/VBFH125toTauTau_cfi.py -n 25 -p 25 -s -a ${p}";
    python scripts/submitLocalHGCalProduction.py -q 1nw -n 1000 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/QCDFlatPt15to3000_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/QCDForPF_14TeV_cfi.py -n 25 -p 1 -s -a ${p}";
done



#test alternative physics lists for pions
phys=("QGSP_FTFP_BERT_EML" "FTFP_BERT_EML" "FTFP_BERT_XS_EML" "QBBC")
pids=(211)
energies=(30)
for pid in ${pids[@]}; do
    for en in ${energies[@]}; do
    	for p in ${phys[@]}; do
        python scripts/submitLocalHGCalProduction.py -q 1nd -n 10 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${pid}_${CMSSW_VERSION}_${p} -p ${pid} -n 400 -e ${en} -l ${p}";
        done
    done
done


a=(lpchgcal/HGCAL_Samples/chgdPionFixedEAndEta_withPFRecHits_SLHC20_patch1_140PU lpchgcal/HGCAL_Samples/chgdPionFixedEAndEta_withPFRecHits_SLHC20_patch1_200PU lpchgcal/HGCAL_Samples/chgdPionFixedEAndEta_withPFRecHits_SLHC20_patch1_20PU lpchgcal/HGCAL_Samples/chgdPionFixedEAndEta_withPFRecHits_SLHC20_patch1_75PU lpchgcal/HGCAL_Samples/chgdPionFixedEAndEta_withPFRecHits_SLHC20_patch1_NoPU)
for i in ${a[@]}; do 
    cmsRun test/runHGCSimHitsAnalyzer_cfg.py ${i}; 
done

### Minimum bias (1000 events per file x 500 jobs, should be ok for later mixing with particle gun)

python scripts/submitLocalHGCalProduction.py -q 2nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/MinBias_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/minBias_cfi.py -n 500";

python scripts/submitLocalHGCalProduction.py -q 2nd -n 500 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/MinBias_v4_${CMSSW_VERSION} -c UserCode/HGCanalysis/python/minBias_cfi.py -n 500 -g Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco";

### Other processes

Can use the minimum bias example, just substitute the argument passed in the -c option to point to the new cfi snippet.

### Redigitization with pileup mixing (will run one job per file, randomizing the min.bias files at start)

tags=(Single211_${CMSSW_VERSION})
pu=(140 100 200)
for tag in ${tags[@]}; do
    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${tag}/RECO | awk '{print $5}'`)
    nFiles=${#inputFiles[@]};
    echo "Submitting $nFiles for ${tag}"
    for p in ${pu[@]}; do
    	python scripts/submitLocalHGCalProduction.py -n ${nFiles} -q 1nw -s digitizeAndMix.sh -o "-o /store/cmst3/group/hgcal/CMSSW/${tag}/ReRECO_PU${p} -m MinBias_${CMSSW_VERSION} -t ${tag}/RECO -p ${p}";
    done
done
    

## Producing analysis ntuples

The ntuples are produced by plugins/HGCSimHitsAnalyzer.cc.  Change the code according to your needs.
To submit the production of the ntuples you can use the following script (it will printout the options)

cmsRun runHGCSimHitsAnalyzer_cfg.py

Submit several jobs to the batch and store the output in EOS
tags=("Single211_CMSSW_6_2_0_SLHC21") # "Single2212_CMSSW_6_2_0_SLHC21")
#tags=("Single22_CMSSW_6_2_0_SLHC21")
for tag in ${tags[@]}; do
    inputFiles=(`cmsLs /store/cmst3/group/hgcal/CMSSW/${tag}/RECO-v4 | awk '{print $5}'`);
    nFiles=${#inputFiles[@]};	
    nJobs=$((nFiles/50));
    for i in `seq 0 ${nJobs}`; do
    	startFile=$((i*=50));
    	cmsRun test/runHGCSimHitsAnalyzer_cfg.py ${tag}/RECO-v4 ${startFile} 50 & 
    done
done


