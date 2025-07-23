#!/bin/bash



#SBATCH --job-name=DurhamDirTest                                                       # job name
#SBATCH -t 1-0:00                                                                       # time limit: (D-HH:MM) 
#SBATCH --output=outfiles/DurhamDirTest.out                                            # file to collect standard output
#SBATCH --error=outfiles/DurhamDirTest.err                                             # file to collect standard errors
#SBATCH --nodelist=p1cmp045                                                             # original = p1cmp045



#srun --pty -p largemem --time=0-12:00:00 bash                                           # only run first time!
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "root 6.22.00-python3-x86_64-centos7-gcc8-opt"

./VjetsPythia8 VjetsPythia8_kTDurham.cmnd -outroot DurhamDirTest > DurhamDirTest.txt