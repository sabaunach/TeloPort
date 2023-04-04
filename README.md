# TeloPort

TeloPort / TeloReport is a collection of tools for extracting and clustering subtelomeric sequences from raw reads.

## About

Subtelomeres are often studied because they undergo frequent rearrangements. However, due to the repeating nature of telomeres, many reads containing valuable subtelomere information are difficult to place in the genome assembly. This project was created as a means of extracting and processing subtelomeric reads from raw read datasets for the purposes of analysis and for more accurate construction of subtelomeric regions in fully assembled genomes.

This project was developed in collaboration with Trey Stansfield (WKU) and the following mentors:

* UKY: [Dr. Mark Farman](https://plantpathology.ca.uky.edu/person/mark-farman) and [Dr. Jerzy Jaromczyk](https://www.engr.uky.edu/directory/jaromczyk-jerzy)
* EKU: [Dr. Patrick Calie](https://biology.eku.edu/people/calie)
* WKU: [Dr. Claire Rinehart](https://www.wku.edu/biology/staff/claire_rinehart)

## Getting Started

### Installation

#### Requirements

* [Boost](https://www.boost.org/) (required)
  * installing the `libboost-all-dev` package should do it
* [GCC](https://gcc.gnu.org/install/) and [Make](https://www.gnu.org/software/make/) (required; check versions using ```g++ -v``` and ```make -v```)
* [MUSCLE](https://www.drive5.com/muscle/downloads.htm) (optional; required for clustering and sequence alignment)
  * also available through [Bioconda](https://bioconda.github.io/)
* [wcdest](https://github.com/shaze/wcdest) (optional; required for sequence alignment and generating consensus sequence)

#### Building the Project

To compile the project, simply run the commmand ```make``` while in the main directory. This will create the **build/** directory. All of the included programs are within the **build/apps** directory.

# Using TeloPort

**IMPORTANT:**

It is recommended to add the **build/apps** to your ```$PATH```, and/or create a separate folder to store all of the output created by each program.

In a similar vein, remember that any files created in the **build/** directory will be removed upon running ```make```.

## telomereFinder

This program is the first step in the pipeline. It extracts all raw reads containing telomeric repeats. This is done by seeking a single telomere repeat sequence, then scoring the remaining hypothesized telomere repeat sequence using a Wagner-Fischer matrix.

Options:
```
  -h [ --help ]                         produce this message
  -s [ --separated ] arg                specify input from separated paired-end
                                        read fastq files (r1 then r2).
  -i [ --interleaved ] arg              specify input from an interleaved
                                        paired-end read fastq file.
  -o [ --out ] arg (=telomereFinder_out)
                                        specify output directory. default:
                                        telomereFinder_out/
  -n [ --nRatio ] arg (=0.05)           throw out reads with a ratio of N's
                                        greater than specified. default: 0.05
  -t [ --telRepeat] arg (=CCCTAA)       specify (forward) telomeric repeat
```

Examples:
```
mkdir ./genome1/
telomereFinder -i ./genome1_reads_interleaved.fastq -o ./genome1/tfout/
telomereFinder -s ./genome1_r1.fastq genome1_r2.fastq -o ./genome1/tfout/
```

## junctionFinder

This program is the second step in the pipeline. It finds the telomere junction, where telomeric sequence and subtelomeric sequence meet. The goal is to find the most precise boundary possible with little given information. This is done using a two-step sliding-window approach - the first window finds where the telomeric sequence begins to drop off, and the second window aims to find the exact boundary between telomeric and subtelomeric sequence.

Options:
```
  -h [ --help ]                         produce this message
  -i [ --in ] arg                       specify input from a fastq file.
  -o [ --out ] arg (=junctionFinder_out)
                                        specify output directory.
  -T [ --startTries ] arg (=10)         specify how many times to try a start
                                        window.
  -l [ --startSkipLen ] arg (=12)       specify how many bp to skip if the
                                        first start fails.
  -W [ --windowLen ] arg (=60)          specify window length on first pass.
  -L [ --stepLen ] arg (=1)             specify step length on first pass.
  -C [ --cutoff ] arg (=0.5)            specify cutoff penalty for first pass.
  -w [ --shortWindowLen ] arg (=15)     specify window size for second pass.
  -c [ --shortCutoff ] arg (=1.375)     specify cutoff penalty for second pass.
  -r [ --repeatsThrow ] arg (=4)        ignore reads with at least u telomere
                                        repeats remaining after cut. (=0 to
                                        turn off)
  --sort arg (=1)                       sort output based on telomere length
  --splitDir arg (=1)                   split output based on telomere
                                        orientation (forward/reverse)
  --splitJunc arg (=1)                  split output based on junction location
  --revc arg (=0)                       output reads in same orientation (by
                                        revc-ing revc reads)
  --fullMatches arg (=0)                include reads that are entirely
                                        telomeric sequence
  --pretty arg (=0)                     print each cut with color to stdout
  -t [ --telRepeat] arg (=CCCTAA)       specify (forward) telomeric repeat
```

Example:
```
junctionFinder -i ./genome1/tfout/telReads.fastq -o ./genome1/jfout/ --revc=true --splitDir=false
```

These settings will create a format more suitable for using with wcdest (after using sequenceQuality) later. The other settings can be used to adapt the output for your specific needs. A more detailed description of the algorithm and how to adjust window lengths and cutoffs can be found by running ```junctionFinder -h```.

## sequenceQuality

This program is used to trim sequences to a specified length for downstream applications.

Options:
```
  -h [ --help ]                 produce this message
  -i [ --in ] arg               specify input file.
  -o [ --out ] arg              specify output file. (blank for default)
  -f [ --ifmt ] arg (=fastq)    specify input format. (fastq/fasta)
  -g [ --ofmt ] arg (=fastq)    specify output format. (fastq/fasta)
  -c [ --cut ] arg (=0)         cut sequence to specified length. 0=no cut.
                                negative=cut from end.
  -l [ --lengththrow ] arg (=0) throw out sequences less than specified length.
                                0=no restriction. negative=throw out greater
  --log                         enable logging
```

Example:
```
mkdir ./genome1/sqout/
sequenceQuality -i ./genome1/jfout/subTelSeq.fastq -o ./genome1/sqout/subTelCut.fasta --ofmt fasta -c 60 -l 40
```

This will create a file containing trimmed subtelomeres ready to feed into wcdest. Enforcing a length of at least 40 and a cutoff of 60 helps wcd later on. Lost sequence information is restored after running wcdest through the wcdInterrogate program.

## wcdest

wcdest (pronounced "wicked") is a d2-clustering program created by [Scott Hazelhurst](https://github.com/shaze). We can use it to cluster the subtelomeres from our raw reads so that we can begin to set up a frequency table, analyze our subtelomeres, look for newly formed telomeres, and match subtelomeres to actual chromosome ends.

Example:
```
mkdir ./genome1/wcdout/
wcd ./genome1/sqout/subTelCut.fasta -l 40 -T 5 -H 0 --show_clusters --histogram > ./genome1/wcdout/genome1clusters.wcd
```

The options ```-T 5 -H 0``` are what worked best for our data; however, it is encouraged for adjustments to be made as needed to obtain clean clusters.

## wcdInterrogate

This program is used to link cluster info generated from wcdest to the original sequence data. These clusters can be analyzed by frequency or passed into the sequence alignment tool MUSCLE.

Options:
```
  -h [ --help ]                produce this message
  -w [ --wcd ] arg             specify wcd file.
  -i [ --wcdin ] arg           specify fasta file used for wcd input.
  -o [ --out ] arg             specify output file for displaying clusters
  -r [ --outmfadir ] arg       specify output directory for cluster multifasta
                               files (blank for no do)
  -s [ --seqin ] arg           specify file with sequence info (blank for
                               wcdin)
  -f [ --seqfmt ] arg (=fasta) specify format of seqin
  -c [ --cons ] arg            specify file with consensus sequence info and
                               display consensus sequence info with clusters
                               (not required)
  --size arg (=2)              only use clusters of size >= size
  --ids                        display IDs of sequences in clusters
  --indices                    display indices of sequences in clusters
  --aux                        display aux information from wcd output in
                               stdout
  --sort                       sort clusters by cluster size
  --log                        enable logging
```

Example:
```
mkdir ./genome1/wintout/
wcdInterrogate -w ./genome1/wcdout/genome1clusters.wcd -i ./genome1/sqout/subTelCut.fasta -s ./genome1/jfout/subTelSeq.fastq./genome1/jfout/subTelSeq.fastq -f fastq -o ./genome1/wintout/cluster.out -r ./genome1/wintout/mfacluster/ --indices --sort --size=8
```

This will create a multifasta file for each cluster containing each sequence as it appeared in the raw read dataset. These multifasta files can then be passed into MUSCLE.

## MUSCLE

[MUSCLE](https://www.drive5.com/muscle/) is a sequence alignment tool that is used to align sequences and generate a consensus sequence for each cluster.

Example:
```
muscle -in ./wcdInterrogate_out/B51_revc/mfacluster_c60_l40_t5/cluster#.fasta -out ./muscle_out/B51_revc/mfacluster_c60_l40_t5/cons#.fasta
```

## Presentations and Usage

[Kentucky Academy of Science 2020 Virtual Annual Meeting](https://mms.kyscience.org/members/publication/program_issue.php?iid=790313) (use ctrl-f and search 'TeloReport')

The program is currently being used by Trey Stansfield and the Farman Lab for continuing work in researching subtelomeric regions.

