cd ../../../../../TongueTwister2
make
cd ../TongueTwisterReader2
make
cd ../Experiments/Test/Tree1/GTR/Run1
rm Execute/out.config
rm Execute/*.json
rm Execute/out.tre
rm Execute/out.tsv
../../../../../TongueTwister2/tonguetwister -d config.json -o Execute/out
rm Execute/consensus.tre
rm Execute/alignments.nytril
rm Execute/alignments_full.nytril
../../../../../TongueTwisterReader2/ttread -b 0.1 -o Execute/out -f true -d Execute -n Execute
