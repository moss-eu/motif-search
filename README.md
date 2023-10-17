## Overview ##
Motif search tool for accurate data reconstruction for enzymatic synthesis-based, error-prone DNA storage. 

The details are introduced in the paper: Yan, Y., Pinnamaneni, N., Chalapati, S. et al. Scaling logical density of DNA storage with enzymatically-ligated composite motifs. Sci Rep 13, 15978 (2023). DOI https://doi.org/10.1038/s41598-023-43172-0

## Pre-requirement ##

### Intel TBB ###

- source: https://github.com/01org/tbb/releases/tag/2019_U5
- libtbb-dev package

### Installation ###

* clone
* Build it: `make`


## Command line ##

motif-search options:
```
  -h [ --help ]         Display help
  -t [ --ncpus ] arg    Number of threads to use
  -m [ --motifs ] arg   File containing motifs to be mapped
  -r [ --read ] arg     File(s) containing reads
  -l [ --kmerlen ] arg  Length of kmer (5)
  -o [ --output ] arg   Ouptut file(default stdout)
  --fp arg              Forward primer
  --rp arg              Reverse primer
  --s arg               Spacer
  --mode arg            Mode L(long) or S (short)
  --n arg               Number of motifs in each oligo
  --id arg              The motif source id list for each motif (0-index). No 
                        need to specify if all motifs come from sane motif 
                        source file.
  -i [ --itype ] arg    Input format (bonito or other), guppy by default
```

### Short mode example ###

1. Get all inferred oligos
```
path-to-motif-search/motif-search -t 12 -i other \
--s CTACAACGCAGATTACAACCTCAGT --mode S --n 2 --id 0 1 \
-m $path/motif/Pos1_Motifs_Exd$m1.fasta $path/motif/Pos3_Motifs_Exd$m2.fasta -r $path/read/BOA-cov$cov.string \
-o $path/res/BOA-$m1-$m2-cov$cov.res
```

2. Filter qualified oligos by consensus check (occur more than once)
```
short-consensus.py input  output
```

e.g. 
```
python path-to-motif-search/script/BOA/final/short-consensus.py \
$path/res/BOA-$m1-$m2-cov$cov.res $path/res/inferred-oligo-$m1-$m2-cov$cov.res
```

### Long mode example ###

1. Get all inferred oligos
```
path-to-motif-search/motif-search -t 12 \
--s CTACAACGCAGATTACAACCT --fp GATTACAACCT --rp CTACAACGCA --mode L --n 20 \
-m $path/motif/motifs-payload -r $path/read/read-$type-$cov.fastq \
-o $path/res/read-$type-$cov.res
```

2. Filter qualified oligos by consensus check (valid index, majority of each motif)
```
long-consensus.py nb_motif nb_oligo nb_oligo_in_binary_format input output
```

e.g. 
```
nb_motif=20
nb_oligo=2750
nb_oligo_bi=101010111110 #binary of 2750
decoded_motif_path=$path/res/read-$type-$cov.res
out_path=$path/res/read-$type-$cov.csv
python path-to-motif-search/script/EMP/long-consensus.py $nb_motif $nb_oligo $nb_oligo_bi $decoded_motif_path $out_path
```


## Simulation ##
To test the original oligos' recover rate, we simulate the Nanopore reads with `badread` (https://github.com/rrwick/Badread).
The reads are simulated in 3 types: good, medium, bad. The coverage is figured by `--quantity`.

### Good read ###
To simulate good read:
```
path-to-badread/Badread/badread-runner.py \
simulate --reference ref.fasta \
--quantity $cov --error_model random \
--qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 \
--chimeras 0 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" \
> $path/$cov/read-good-$cov.fastq
```

### Medium read ###
To simulate medium read we use the default configuration nanopore2020 error model:
```
path-to-badread/Badread/badread-runner.py \
simulate --reference /media/ssd/ngs-data-analysis/motif/mc_testbedenv_05/badread/ref.fasta \
--quantity $cov > $path/$cov/read-$cov.fastq
```

### Bad read ###
To simulate bad read:
```
path-to-badread/Badread/badread-runner.py \
simulate --reference /media/ssd/ngs-data-analysis/motif/mc_testbedenv_05/badread/ref.fasta \
--quantity $cov --glitches 1000,100,100 \
--junk_reads 5 --random_reads 5 --chimeras 10 --identity 75,90,8 \
> $path/$cov/read-bad-$cov.fastq
```

[comment]: <> (## Result ##)

[comment]: <> (### Result based on different quality ###)

[comment]: <> (These tests use the previous simulation command with `cov=80`.)

[comment]: <> (* Good read:)

[comment]: <> (    - Time: 4 min )

[comment]: <> (    - Recovered motifs: 100%)

[comment]: <> (    - Fully recovered oligos: 100%)

[comment]: <> (* Medium read:)

[comment]: <> (    - Time: 4 min)

[comment]: <> (    - Recovered motifs: 99.82%)

[comment]: <> (    - Fully recovered oligos: 96.62%)

[comment]: <> (* Bad read:)

[comment]: <> (    - Time: 2 min &#40;many reads are ignored because too long or too short&#41;)

[comment]: <> (    - Recovered motifs: 28.82%)

[comment]: <> (    - Fully recovered oligos: 0.4%)

[comment]: <> (### Result based on different coverage ###)

[comment]: <> (These tests are based on medium read with different coverage.)

[comment]: <> (* `cov=10`)

[comment]: <> (    - Time: 30 sec)

[comment]: <> (    - Recovered motifs: 67.85%)

[comment]: <> (    - Fully recovered oligos: 2.95%)
    
[comment]: <> (* `cov=20`)

[comment]: <> (    - Time: 1 min)

[comment]: <> (    - Recovered motifs: 87.32%)

[comment]: <> (    - Fully recovered oligos: 24.29%)
    
[comment]: <> (* `cov=40`)

[comment]: <> (    - Time: 2 min)

[comment]: <> (    - Recovered motifs: 98.19%)

[comment]: <> (    - Fully recovered oligos: 73.89%)
    
[comment]: <> (* `cov=80`)

[comment]: <> (    - Time: 4 min)

[comment]: <> (    - Recovered motifs: 99.82%)

[comment]: <> (    - Fully recovered oligos: 96.62%)

[comment]: <> (* `cov=200`)

[comment]: <> (    - Time: 10 min)

[comment]: <> (    - Recovered motifs: 99.99%)

[comment]: <> (    - Fully recovered oligos: 99.89%)
 

