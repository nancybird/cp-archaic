# cp-archaic
ChromoPainter-based method for inferring local ancestry, specifically suited for archaic ancestry.

Compile with:  gcc -Wall -o cp-archaic cp-archaic.c -lm -lz

# Example usage to infer Neanderthal segments in 1000 genomes

Method requires 'outgroup' assumed to have minimal amounts of introgression, and high-coverage archaic references (we use Altai, Chagyrskaya, Vindija and Denisovan, more is better!). These, along with target individuals must be phased and in ChromoPainter format. See https://github.com/nancybird/vcftochromopainter for code to convert to ChromoPainter format from a phased vcf.

**Step 1:** 

Paint both the outgroup and the archaics using the outgroup and archaics as donors. One of the archaics must be excluded as a donor (we use Altai Neanderthal, as simulations show it is better to have the archaic closer to the introgressing source, i.e. Vindija, as a donor). 

A popfile for step 1 could e.g. look like: 

Vindija D \
Chagyrskaya D \
Denisova D \
YRI D \
Vindija R \
Chagyrskaya R \
Denisova R \
AltaiNeandertal R \

**Another important note is when painting the archaics, you should leave one outgroup individual out, so that in the end all individuals are painted by the same number of outgroup individuals (because an outgroup individual cannot be painted by themselves)**

**Step 2:**

Paint the targets to get inferred archaic ancestry segments.

-g Chromoainter haps file (no default, required)

-r ChromoPainter recomrate file (no default, required)

-t ChromoPainter idfile (no default, required)

-f2 Popfile (see below for example)

-c chunklengths file from **Step 1**

-s 0 (don't make samples files)

-s0 print out start/end points of each inferred archaic segment

-l <double1> <double2> <double3>  inferred archaic segments must be >double1 before and >double2 after stitching across segments with gap <double3 between them (in kb); otherwise remove (default and recommended settings = 10 20 40)

-e  <double> in target individual, SNPs with Pr(archaic) >= this value will be called 'archaic' (initially, before remove segments/merging) if using '-f2' switch (recommended value= 0.1 from simulations)

-o output file


Example command for 1000 genomes:

./cp-archaic \
-g 1000G_masked_chr${chr}.chromopainter.inp.haps.gz \
-r 1000G_masked_chr${chr}.chromopainter.inp.recomrates \
-t 1000g_idfile.txt \
-f2 target.popfile.txt \
-c Painting/1000G_mask_all.chunklengths.out \
-s0 -l 10 20 40 -e 0.1 -s 0 \
-o 1000G_targets_CHB_masked_chr_${chr}


An example target.popfile.txt would be:

Vindija D NA \
Chagyrskaya D NA \
Denisova D NA \
YRI D NA \ 
YRI S 0.98 \
AltaiNeandertal A 0.02 \
CHB T NA \

where 0.02 gives the prior genome-wide probability of archaic ancestry. We find that varying the prior probabilty (up to 10%) does not significantly alter results. 

**Tip:** For speed, we recommend running cp-archaic on targets in batchs (e.g. 10 inds at a time) using specifications after the popfile e.g. -f2 target.popfile.txt 1 9 


If using cp-archaic please cite XX

