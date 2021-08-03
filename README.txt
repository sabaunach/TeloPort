

# Generating input for wcd (subtelomeric sequences cut to length 60, with length at least 40)
# cut subtelomeric sequences
./Project/build/apps/sequenceQuality -i ./junctionFinder_out/B51/subTelSeq_r.fastq -o ./sequenceQuality_out/B51/subTelSeqCut60_r.fasta --ofmt fasta -c -60 -l 40
./Project/build/apps/sequenceQuality -i ./junctionFinder_out/B51/subTelSeq_f.fastq -o ./sequenceQuality_out/B51/subTelSeqCut60_f.fasta --ofmt fasta -c 60 -l 40
# concatenate into one file for wcd clustering
cat ./sequenceQuality_out/B51/subTelSeqCut60_f.fasta ./sequenceQuality_out/B51/subTelSeqCut60_r.fasta > ./sequenceQuality_out/B51/subTelSeqCut60.fasta
# remove temps
rm ./sequenceQuality_out/B51/subTelSeqCut*_f.fasta ./sequenceQuality_out/B51/subTelSeqCut*_r.fasta

# Generate wcd clusters:
wcd ./sequenceQuality_out/B51/subTelSeqCut60.fasta -l 40 -T 5 -H 0 --show_clusters --histogram > ./wcd_out/B51/cluster_c60_l40_t5.wcd

# Generate cluster multifastas for muscle input (and readable output info):
# Relevant options:
#  -w [ --wcd ] arg             specify wcd file.
#  -i [ --wcdin ] arg           specify fasta file used for wcd input.
#  -o [ --out ] arg             specify output file for displaying clusters
#  -r [ --outmfadir ] arg       specify output directory for cluster multifasta
#                               files (blank for no do)
#  -s [ --seqin ] arg           specify file with sequence info (blank for
#                               wcdin)
#  -f [ --seqfmt ] arg (=fasta) specify format of seqin
#  -c [ --cons ] arg            specify file with consensus sequence info and
#                               display consensus sequence info with clusters
#                               (not required)
#  --size arg (=2)              only use clusters of size >= size
#  --ids                        display IDs of sequences in clusters
#  --indices                    display indices of sequences in clusters
#  --aux                        display aux information from wcd output in
#                               stdout
#  --sort                       sort clusters by cluster size
# Controlling WHICH sequence is shown for sequences in cluster (subtelomeric, telomeric, full, cut, etc.) is through -s.

./Project/build/apps/wcdInterrogate -w ./wcd_out/B51/cluster_c60_l40_t5.wcd -i ./sequenceQuality_out/B51/subTelSeqCut60.fasta -s ./junctionFinder_out/B51/subTelSeq.fastq -f fastq -o ./wcdInterrogate_out/B51/cluster_c60_l40_t5.out -r ./wcdInterrogate_out/B51/mfacluster_c60_l40_t5/ --indices --sort --size=8

# Command for sequence alignment:
mkdir ./muscle_out/B51/mfacluster_c60_l40_t5/
muscle -in ./wcdInterrogate_out/B51/mfacluster_c60_l40_t5/cluster#.fasta -out ./muscle_out/B51/mfacluster_c60_l40_t5/cons#.fasta
# ISSUE: wcd takes care of reverse complement stuff but muscle alignment does not. Need to revc seqs beforehand!
# METHOD: split r and f reads before running wcd, or only do f, or revc r reads in junctionFinder
#	going to revc reads in junctionFinder
# New pipeline:
./Project/build/apps/junctionFinder -i ./telomericreadfinder_out/B51/telReads.fastq -o ./junctionFinder_out/B51_revc/ --revc=true --split=false
mkdir ./sequenceQuality_out/B51_revc/
./Project/build/apps/sequenceQuality -i ./junctionFinder_out/B51_revc/subTelSeq.fastq -o ./sequenceQuality_out/B51_revc/subTelSeq_c60_l40.fasta --ofmt fasta -c 60 -l 40
mkdir ./wcd_out/B51_revc/
wcd ./sequenceQuality_out/B51_revc/subTelSeq_c60_l40.fasta -l 40 -T 5 -H 0 --show_clusters --histogram > ./wcd_out/B51_revc/cluster_c60_l40_t5.wcd
mkdir ./wcdInterrogate_out/B51_revc/
./Project/build/apps/wcdInterrogate -w ./wcd_out/B51_revc/cluster_c60_l40_t5.wcd -i ./sequenceQuality_out/B51_revc/subTelSeq_c60_l40.fasta -s ./junctionFinder_out/B51_revc/subTelSeq.fastq -f fastq -o ./wcdInterrogate_out/B51_revc/cluster_c60_l40_t5.out -r ./wcdInterrogate_out/B51_revc/mfacluster_c60_l40_t5/ --indices --sort --size=8
mkdir ./muscle_out/B51_revc ./muscle_out/B51_revc/mfacluster_c60_l40_t5/
muscle -in ./wcdInterrogate_out/B51_revc/mfacluster_c60_l40_t5/cluster#.fasta -out ./muscle_out/B51_revc/mfacluster_c60_l40_t5/cons#.fasta
