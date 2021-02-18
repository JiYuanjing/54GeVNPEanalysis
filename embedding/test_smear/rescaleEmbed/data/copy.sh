#! /bin/bash
array=("gamma2pi" "gamma2eta" "pi0" "dirpho" "eta")
# array=("gamma2pi" "gamma2eta")
cprcf() {
 scp -r jiyj@rftpexp.rhic.bnl.gov:$1 $2
}
for i in "${array[@]}"
do
  cprcf /star/u/jiyj/pwg/54GeV/embedding/PhotoE/test_smearpt/$i/embeddQa_${i}_smear.root .
done


