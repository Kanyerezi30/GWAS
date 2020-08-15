#!/bin/bash
#Making a binary PED file
#The first thing we will do is to make a binary PED file. This more compact representation of the data saves space and speeds up\
#subsequent analysis. To make a binary PED file, use the following command.

plink --file hapmap1 --make-bed --out hapmap1

#If it runs correctly on your machine, you should see the following in your output:
     above as before
     ...
     Before frequency and genotyping pruning, there are 83534 SNPs
     Applying filters (SNP-major mode)
     89 founders and 0 non-founders found
     0 SNPs failed missingness test ( GENO > 1 )
     0 SNPs failed frequency test ( MAF < 0 )
     After frequency and genotyping pruning, there are 83534 SNPs
     Writing pedigree information to [ hapmap1.fam ]
     Writing map (extended format) information to [ hapmap1.bim ]
     Writing genotype bitfile to [ hapmap1.bed ]
     Using (default) SNP-major mode
     Analysis finished: Mon Jul 31 09:10:05 2006
#There are several things to note:
#When using the --make-bed option, the threshold filters for missing rates and allele frequency were automatically set to exclude nobody.\
#Although these filters can be specified manually (using --mind, --geno and --maf) to exclude people, this default tends to be wanted\
#when creating a new PED or binary PED file. The commands --extract / --exclude and --keep / --remove can also be applied at this stage.
#Three files are created with this command -- the binary file that contains the raw genotype data hapmap1.bed but also a\
#revsied map file hapmap1.bim which contains two extra columns that give the allele names for each SNP, and\
#hapmap1.fam which is just the first six columns of hapmap1.ped. You can view the .bim and .fam files -- but do not\
#try to view the .bed file. None of these three files should be manually editted.
#If, for example, you wanted to create a new file that only includes individuals with high genotyping (at least 95% complete),\
#you would run:

plink --file hapmap1 --make-bed --mind 0.05 --out highgeno

#which would create files

     highgeno.bed
     highgeno.bim
     highgeno.fam

#Working with the binary PED file
#To specify that the input data are in binary format, as opposed to the normal text PED/MAP format, just use the --bfile option instead\
#of --file. To repeat the first command we ran (which just loads the data and prints some basic summary statistics):

plink --bfile hapmap1
     Writing this text to log file [ plink.log ]
     Analysis started: Mon Jul 31 09:12:08 2006

     Options in effect:
             --bfile hapmap1

     Reading map (extended format) from [ hapmap1.bim ]
     83534 markers to be included from [ hapmap1.bim ]
     Reading pedigree information from [ hapmap1.fam ]
     89 individuals read from [ hapmap1.fam ]
     89 individuals with nonmissing phenotypes
     Reading genotype bitfile from [ hapmap1.bed ]
     Detected that binary PED file is v1.00 SNP-major mode
     Before frequency and genotyping pruning, there are 83534 SNPs
     Applying filters (SNP-major mode)
     89 founders and 0 non-founders found
     0 of 89 individuals removed for low genotyping ( MIND > 0.1 )
     859 SNPs failed missingness test ( GENO > 0.1 )
     16994 SNPs failed frequency test ( MAF < 0.01 )
     After frequency and genotyping pruning, there are 65803 SNPs
     Analysis finished: Mon Jul 31 09:12:10 2006
#The things to note here:
#That three files hapmap1.bim, hapmap1.fam and hapmap1.bed were loaded instead of the usual two files. That is, hapmap1.ped and\
#hapmap1.map are not used in this analysis, and could in fact be deleted now.
#The data are loaded in much more quickly -- based on the timestamp at the beginning and end of the log output, this took\
#2 seconds instead of 10.
#Summary statistics: missing rates
#Next, we shall generate some simple summary statistics on rates of missing data in the file, using the --missing option:

plink --bfile hapmap1 --missing --out miss_stat

#which should generate the following output:
     ...
     0 of 89 individuals removed for low genotyping ( MIND > 0.1 )
     Writing individual missingness information to [ miss_stat.imiss ]
     Writing locus missingness information to [ miss_stat.lmiss ]
     ...

#Here we see that no individuals were removed for low genotypes (MIND > 0.1 implies that we accept people with less than\
#10 percent missingness).
#The per individual and per SNP (after excluding individuals on the basis of low genotyping) rates are then output to the\
#files miss_stat.imiss and miss_stat.lmiss respectively. If we had not specified an --out option, the root output filename would have\
#defaulted to "plink".
#These output files are standard, plain text files that can be viewed in any text editor, pager, spreadsheet or statistics package\
#(albeit one that can handle large files). Taking a look at the file miss_stat.lmiss, for example using the more command which is present\
#on most systems:
#more miss_stat.lmiss
#we see

      CHR          SNP   N_MISS     F_MISS
        1    rs6681049        0          0
        1    rs4074137        0          0
        1    rs7540009        0          0
        1    rs1891905        0          0
        1    rs9729550        0          0
        1    rs3813196        0          0
        1    rs6704013        2  0.0224719
        1     rs307347       12   0.134831
        1    rs9439440        2  0.0224719
      ...

#That is, for each SNP, we see the number of missing individuals (N_MISS) and the proportion of individuals missing (F_MISS). Similarly:
#more miss_stat.imiss
#we see

         FID          IID MISS_PHENO     N_MISS     F_MISS
      HCB181            1          N        671 0.00803266
      HCB182            1          N       1156  0.0138387
      HCB183            1          N        498 0.00596164
      HCB184            1          N        412 0.00493212
      HCB185            1          N        329 0.00393852
      HCB186            1          N       1233  0.0147605
      HCB187            1          N        258 0.00308856
      ...

#The final column is the actual genotyping rate for that individual -- we see the genotyping rate is very high here.
