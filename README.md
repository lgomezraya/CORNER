# CORNER
Code to perform the Corner ALgorithm for detection of situations where the four gamete test fails
[README.txt](https://github.com/lgomezraya/CORNER/files/8886391/README.txt)
README
This program is written in FORTRAN and can be compiled and executed with gfortran of the GNU.
This program identify situations where the four-gamete test fails (less that for haplotypes for a pair of SNPs) using the Corners algorithm and estimate linkage disequilibrium as described by Gomez-Raya and Rauw (GSE).
The programs is ready to run for 435 individuals. All individuals has to have a recorded genotype for every SNP (no missing genotype information!). The program can be adapted to another needs manipulating the variables (ni =number of individuals, nmar= number of markers, etc).
The program include program corner.f and routine em.f for the em algorithm. It can be executed with 
sh corner.total
the instruction in this file are:
rm ko
gfortran -o ko corner.f em.f
./ko
There are two input files freq.all and snpv11.all. The file  freq.all has input variables:
 nchro # number of chromosome
pos: position of the marker
k; group number not in use in this program
aa: genotype for one of the homozygous alleles
ab: genotype for heterozygotes
bb: genotype other homozygous allele
f11: frequency of aa
f12: frequency of ab
f22: frequency of bb
The input snp genotypic information is in file snpv11.all (DOI 10.5281/zenodo.6636369)
The input information on this snpv11.all is: 
ind: identification of the individual number
pos: base pair of the SNP marker
gr:  group number not in use in this program
mar: name of the SNP
nchro: number of chromosome (only 1 to 18)
alle1: allele 1 at the genotype
alle2: allele 2 at the genotype
This file should be sorted by marker and number of individual within marker. No animal or SNP can be missing.

One output file is  “corner.all” lists all chromosomal regions where any SNP pair combination fails within the region the four-gamete test. 
The output has the format
chromosme #,
initial number of SNP, 
ending number of SNP, 
total number of SNPs in the regions, 
initial position of SNP, 
ending position of SNP
Another output file is “consecutive” listing all consecutive SNP pairs that fails the four-gamete test. It has format:
chromosome #,
initial number of SNP;
ending number of SNP;
test for four gamete test: 0 not fail, 1 fails; 
initial position of SNP; 
ending position of SNP.
The program is currently being developed. Any use is at your own risk. No guarantees.

