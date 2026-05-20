# cp-archaic
ChromoPainter-based method for inferring local ancestry, specifically suited for archaic ancestry.

Compile with:  gcc -Wall -o cp-archaic cp-archaic.c -lm -lz

# Example usage to infer Neanderthal segments in 1000 genomes

cp-archaic requires an 'outgroup' assumed to have minimal amounts of introgression (e.g. YRI in 1000G), and high-coverage archaic references (we use Altai, Chagyrskaya, Vindija and Denisovan, more is better!). These, along with target individuals must be merged, phased and in ChromoPainter format. See https://github.com/nancybird/vcftochromopainter for code to convert to ChromoPainter format from a phased vcf. 

For a note on phasing the archaics, please see our paper (link below), but essentially even with reference-free phasing of archaics in simulations, cp-archaic's performance was not significantly impacted. 

Additionally, we recommend applying an archaic mask file to the whole dataset to remove poorly called regions, e.g. ones found here (https://tigress-web.princeton.edu/AKEY/IBDmix_EXCLUDE_masks_hg19_Li_et_al_Science_2024/). 

**Step 1:** 

Paint both the outgroup and the archaics using the outgroup and archaics as donors. One of the archaics must be excluded as a donor (we use Altai Neanderthal, as simulations show it is better to have the archaic closer to the introgressing source, i.e. Vindija, as a donor). 

**Another important note is when painting the archaics, you should leave one outgroup individual out, so that in the end all individuals are painted by the same number of outgroup individuals (because an outgroup individual cannot be painted by themselves)** Therefore we recommend creating 2 separate ID files (one with one outgroup ind set to 0) and painting outgroup and archaics separately. 

An archaic popfile for step 1 could e.g. look like: 

Vindija D \
Chagyrskaya D \
Denisova D \
YRI D (WITH ONE IND LEFT OUT IN ID FILE)\
Vindija R \
Chagyrskaya R \
Denisova R \
AltaiNeandertal R 

An outgroup popfile for step 1 could e.g. look like: 

Vindija D \
Chagyrskaya D \
Denisova D \
YRI D (WITH ALL INDS IN ID FILE)\
YRI R

An example command for step 1 could be:

./cp-archaic \
-g 1000G_masked_chr${chr}.chromopainter.inp.haps.gz \
-r 1000G_masked_chr${chr}.chromopainter.inp.recomrates \
-t 1000g_idfile.txt \
-f2 outgroup.popfile.txt 0 0 (OR ARCHAIC POP FILE SEPARATELY) \
-c Painting/1000G_mask_all.chunklengths.out \
-s 0 \
-o 1000G_outgroup_chr${chr}


Once you have these outfiles, you want to cat archaic and outgroup chunklengths into one file containing paintings from both sets per chromosome, then you can run chromocombine (https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromocombine.html) to get the final chunklengths file for step 2. 



**Step 2:**

Paint the targets to get inferred archaic ancestry segments.

cp-archaic arguments:

-g Chromoainter haps file (no default, required) \
-r ChromoPainter recomrate file (no default, required) \
-t ChromoPainter idfile (no default, required) \
-f2 Popfile (see below for example). Numbers after pop file give indices of target individuals to paint, e.g. 0 0 means paint everyone, 1 9 means paint first 10 inds \
-c chunklengths file from **Step 1** \
-s 0 (don't make samples files) \
-s0 print out start/end points of each inferred archaic segment \
-l <double1> <double2> <double3>  inferred archaic segments must be >double1 before and >double2 after stitching across segments with gap <double3 between them (in kb); otherwise remove (default and recommended settings = 10 20 40) \
-e  <double> in target individual, SNPs with Pr(archaic) >= this value will be called 'archaic' (initially, before remove segments/merging) if using '-f2' switch (recommended value= 0.1 from simulations) \
-o output file \


Example command for 1000 genomes:

./cp-archaic \
-g 1000G_masked_chr${chr}.chromopainter.inp.haps.gz \
-r 1000G_masked_chr${chr}.chromopainter.inp.recomrates \
-t 1000g_idfile.txt \
-f2 target.popfile.txt 0 0 \
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


If using cp-archaic please cite Bird N, Walker E, Hellenthal G. The reliability of inferred archaic segments. Hum Popul Genet Genom. 2026;6(2):0007. https://doi.org/10.47248/hpgg2606020007

