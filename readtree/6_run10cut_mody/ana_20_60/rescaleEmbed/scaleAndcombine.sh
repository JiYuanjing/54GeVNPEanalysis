#! /bin/bash
outDir="rescaleFile"
[ -d $outDir ] || mkdir $outDir

echo "scale pi0"
pi0DalitzBranch="0.01174"
# pi0Embeddfile="data/embeddQa_pi0_0918.root"
# pi0Embeddfile="data/embeddQa_pi0_0929.root"
# pi0Embeddfile="data/embeddQa_pi0_0930.root"
pi0Embeddfile="data/embeddQa_pi0_1016.root"
pi0Rescalefile="$outDir/rescale_pi0.root"
root -b -q generateRecaleHists.C\(\"$pi0Embeddfile\",\"$pi0Rescalefile\",$pi0DalitzBranch\)

echo "scale eta"
etaDalitzBranch="0.48*0.0069"
# etaEmbeddfile="data/embeddQa_eta_0918.root"
# etaEmbeddfile="data/embeddQa_eta_0929.root"
# etaEmbeddfile="data/embeddQa_eta_0930.root"
etaEmbeddfile="data/embeddQa_eta_1016.root"
etaRescalefile="$outDir/rescale_eta.root"
root -b -q generateRecaleHists.C\(\"$etaEmbeddfile\",\"$etaRescalefile\",$etaDalitzBranch\)

echo "scale gamma pi0"
gammaDalitzBranch="1"
# gammaDalitzBranch="0.88"
# gammaDalitzBranch="0.85"
# gammaDalitzBranch="0.5"
# gammaEmbeddfile="data/embeddQa_gamma_0918.root"
# gammaEmbeddfile1="data/embeddQa_gamma_pi0_0929.root"
# gammaEmbeddfile1="data/embeddQa_gamma_pi0_0930.root"
# gammaEmbeddfile1="data/embeddQa_gamma_pi0_1008.root"
# gammaEmbeddfile1="data/embeddQa_gamma_pi0_prod2.root"
# gammaEmbeddfile1="data/embeddQa_gamma_pi0_1014.root"
gammaEmbeddfile1="data/embeddQa_gamma_pi0_1016.root"
gammaRescalefile1="$outDir/rescale_gamma_pi0.root"
root -b -q generateRecaleHists.C\(\"$gammaEmbeddfile1\",\"$gammaRescalefile1\",$gammaDalitzBranch\)

echo "scale gamma eta"
# gammaDalitzBranch="1."
# gammaDalitzBranch="0.95"
# gammaEmbeddfile="data/embeddQa_gamma_0918.root"
# gammaEmbeddfile2="data/embeddQa_gamma_eta_0929.root"
# gammaEmbeddfile2="data/embeddQa_gamma_eta_0930.root"
# gammaEmbeddfile2="data/embeddQa_gamma_eta_1008.root"
# gammaEmbeddfile2="data/embeddQa_gamma_eta_prod2.root"
# gammaEmbeddfile2="data/embeddQa_gamma_eta_1014.root"
gammaEmbeddfile2="data/embeddQa_gamma_eta_1016.root"
gammaRescalefile2="$outDir/rescale_gamma_eta.root"
root -b -q generateRecaleHists.C\(\"$gammaEmbeddfile2\",\"$gammaRescalefile2\",$gammaDalitzBranch\)

echo "scale gamma dir pho"
# gammaDalitzBranch="1"
# gammaEmbeddfile="data/embeddQa_gamma_0918.root"
# gammaEmbeddfile3="data/embeddQa_gamma_dirpho_0929.root"
# gammaEmbeddfile3="data/embeddQa_gamma_dirpho_0930.root"
# gammaEmbeddfile3="data/embeddQa_gamma_dirpho_1002.root"
# gammaEmbeddfile3="data/embeddQa_gamma_dirpho_prod2.root"
gammaEmbeddfile3="data/embeddQa_gamma_dirpho_1016.root"
gammaRescalefile3="$outDir/rescale_gamma_dirpho.root"
root -b -q generateRecaleHists.C\(\"$gammaEmbeddfile3\",\"$gammaRescalefile3\",$gammaDalitzBranch\)

rm rescale_combine.root
hadd rescale_combine.root $pi0Rescalefile $etaRescalefile $gammaRescalefile1  $gammaRescalefile2 $gammaRescalefile3  
