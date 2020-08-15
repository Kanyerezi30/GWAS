#!/bin/bash
plink --file hapmap1 #get the statistics about the ped and map files of hapmap1

#To re-run a previous job, use the --rerun option, which takes a PLINK LOG file as the parameter. This option will scan the LOG file,\
#extract the previous PLINK commands and re-execute them. If new commands are added to the command line, they will also be included; if the\
#command also appeared in the original file, any parameters will be taken from the newer version. For example, if the original command was

plink --file mydata --pheno pheno.raw --assoc --maf 0.05 --out run1

#then the command

plink --rerun run1.log --maf 0.1

#would repeat the analysis but with the new minor allele frequency threshold of 0.1, not 0.05. Note that commands in the old LOG file can be\
#overwritten but not removed with the rerun command.

#MS-DOS only allows command lines to be 127 characters in length -- sometimes, PLINK command lines can grow longer than this. In this case,\
#use the --script option, where the remaining options will be read from a text file. For example,

plink --script myscript1.txt

#where the file myscript1.txt is a plain text file containing

--ped ..\data\version1\50K\allsamples.ped 
--map ..\data\allmapfiles\finalversion\autosomal.map 
--out ..\results\working\sample-missingness-v1.22
--from rs66537222
--to rs8837323 
--geno 0.25
--maf 0.02
--missing

#As well as the --file command described above, PED and MAP files can be specified separately, if they have different names:

plink --ped mydata.ped --map autosomal.map

#Loading a large file (100K+ SNPs) can take a while (which is why we suggest converting to binary format). PLINK will give an error\
#message in most circumstances when something has gone wrong.

#The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype

#The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person. A PED file must have 1\
#and only 1 phenotype in the sixth column. The phenotype can be either a quantitative trait or an affection status column:\
#PLINK will automatically detect which type (i.e. based on whether a value other than 0, 1, 2 or the missing genotype code is observed).

#If an individual's sex is unknown, then any character other than 1 or 2 can be used. When new files are created\
#(PED, FAM, or other which contain sex) then the original coding will be preserved. However, these individuals will be dropped from\
#any analyses (i.e. phenotype set to missing also) and an error message will arise if an analysis that uses family information\
#is requested and an individual of 'unknown' sex is specified as a father or mother.

#To disable the automatic setting of the phenotype to missing if the individual has an ambiguous sex code, add the --allow-no-sex option.\
#When using a data generation command (e.g. --make-bed, --recode, etc) as opposed to an analysis command, then by default the phenotype\
#is not set to missing is sex is missing. This behaviour can be changed by adding the flag --must-have-sex.

#You can add a comment to a PED or MAP file by starting the line with a # character. The rest of that line will be ignored.\
#Do not start any family IDs with this character therefore.
#Affection status, by default, should be coded:

    -9 missing 
     0 missing
     1 unaffected
     2 affected

#If your file is coded 0/1 to represent unaffected/affected, then use the --1 flag:

plink --file mydata --1

#which will specify a disease phenotype coded:

     -9 missing
      0 unaffected
      1 affected

#The missing phenotype value for quantitative traits is, by default, -9 (this can also be used for disease traits as well as 0).\
#It can be reset by including the --missing-phenotype option:

plink --file mydata --missing-phenotype 99

#Other phenotypes can be swapped in by using the --pheno (and possibly --mpheno) option, which specify an alternate\
#phenotype is to be used, described below.

#Genotypes (column 7 onwards) should also be white-space delimited; they can be any character (e.g. 1,2,3,4 or A,C,G,T or anything else)\
#except 0 which is, by default, the missing genotype character. All markers should be biallelic.\
#All SNPs (whether haploid or not) must have two alleles specified. Either Both alleles should be missing (i.e. 0) or neither.\
#No header row should be given. For example, here are two individuals typed for 3 SNPs (one row = one person):

     FAM001  1  0 0  1  2  A A  G G  A C 
     FAM001  2  0 0  1  2  A A  A G  0 0 
     ...

#The default missing genotype character can be changed with the --missing-genotype option, for example:

plink --file mydata --missing-genotype N

#Different values to the missing phenotype or genotype code can be specified for output datasets created, with --output-missing-phenotype\
#and --output-missing-genotype.

#Different PED file formats: missing fields
#Sometimes data arrive in a number of different formats: for example, where the genotype information just has a single ID column\
#followed by all the SNP data, with the other family and phenotype information residing in a separate file.\
#Rather than have to recreate new files, it is sometimes possible to read in such files directly.\
#The standard behavior of PLINK when reading a PED file with --file or --ped can be modified to allow for the fact that\
#one or more of the normally obligatory 6 fields are missing:

--no-fid

#indicates there is no Family ID column: here the first field is taken to be individual ID, and the family ID is\
#automatically set to be the same as the individual ID (i.e. obviously, all individuals would be treated as unrelated).\
#In other files that require family and individual ID (e.g. alternate phenotype file and cluster files, for which this flag has no effect),\
#the individual ID would need to be entered also as the family ID therefore.

--no-parents

#indicates that there are no paternal and maternal ID codes; all individuals would be assumed to be founders in this case

--no-sex

#indicates that there is no sex field; all individuals set to have a missing sex code (which also sets that individual to missing unless\
#the allow-no-sex option is also used)

--no-pheno

#indicates that there is no phenotype filed; all individuals are set to missing unless an alternate phenotype file is specified.
#It is possible to use these flags together, so using all of them would specify the most simple kind of file mentioned above:\
#a single, unique ID code followed by all genotype data.

#If the genotype codes in a PED file are in the form AG rather than A G, for example, such that every genotype is exactly two characters\
#long, then then flag

./plink --file mydata --compound-genotypes

#can be added. Note that this only works for input for PED files (not TPED or LGEN files, and not for any output options, e.g. --recode, etc).

#To load the PED file from the standard input stream instead of a file, use the - symbol as the file name, e.g.

perl retrieve_data.pl | ./plink --ped - --map mymap.map --make-bed

#The MAP file still needs to be a normal file; this currently only works for --ped files.

#MAP files
#By default, each line of the MAP file describes a single marker and must contain exactly 4 columns:

     chromosome (1-22, X, Y or 0 if unplaced)
     rs# or snp identifier
     Genetic distance (morgans)
     Base-pair position (bp units)

#Genetic distance can be specified in centimorgans with the --cm flag. Alternatively, you can use a MAP file with the genetic distance\
#excluded by adding the flag --map3, i.e.

plink --file mydata --map3

#In this case, the three columns are expected to be

     chromosome (1-22, X, Y or 0 if unplaced)
     rs# or snp identifier
     Base-pair position (bp units)

#Base-pair positions are expected to correspond to positive integers within the range of typical human chromosome sizes

#Most analyses do not require a genetic map to be specified in any case; specifying a genetic (cM) map is most crucial for a set of\
#analyses that look for shared segments between individuals. For basic association testing, the genetic distance column can be set at 0.
#SNP identifers can contain any characters except spaces or tabs; also, you should avoid * symbols in names also.
#To exclude a SNP from analysis, set the 4th column (physical base-pair position) to any negative value (this will only work for MAP files,\
#not for binary BIM files).

     1  rs123456  0  1234555
     1  rs234567  0  1237793
     1  rs224534  0  -1237697        <-- exclude this SNP
     1  rs233556  0  1337456
     ...

#The MAP file must therefore contain as many markers as are in the PED file. The markers in the PED file do not need to be in genomic order:\
#(i.e. the order MAP file should align with the order of the PED file markers).

#Chromosome codes
#The autosomes should be coded 1 through 22. The following other codes can be used to specify other chromosome types:

     X    X chromosome                    -> 23
     Y    Y chromosome                    -> 24
     XY   Pseudo-autosomal region of X    -> 25
     MT   Mitochondrial                   -> 26

#The numbers on the right represent PLINK's internal numeric coding of these chromosomes: these will appear in all output rather than the\
#original chromosome codes.
#For haploid chromosomes, genotypes should be specified as homozygotes: for most analyses, PLINK will treat these appropriately.\
#For example, consider the following example PED file, containing two males (1 and 2) and two females (3 and 4):

     1 1 0 0 1   1   A A    A A    A A    A A    A A
     2 1 0 0 1   1   A C    A C    A C    A C    A C
     3 1 0 0 2   1   A A    A A    A A    A A    A A
     4 1 0 0 2   1   A C    A C    A C    A C    A C

#and MAP file

     1    snp1   0   1000
     X    snp2   0   1000
     Y    snp3   0   1000
     XY   snp4   0   1000
     MT   snp5   0   1000

#Generating frequencies for these SNPs,

plink --file test --freq

#we see plink.frq is

      CHR          SNP   A1   A2          MAF       NM
        1         snp1    C    A         0.25        8
       23         snp2    C    A          0.2        5
       24         snp3    C    A            0        1
       25         snp4    C    A         0.25        8
       26         snp5    C    A            0        2

#There are several things to note. First, the numeric chromosome codes are used in the output to represent X, Y, XY and MT. Second, haploid\
#chromosomes are only counted once (i.e. male X and Y chromosome SNPs and all MT SNPs).
#Third, several genotypes have been set to missing if they are not valid (female Y genotype, heterozygous haploid chromosome).
#The NM field represents the number of non-missing alleles for each SNP -- this is because invalid genotypes are automatically set to missing.
#We can see which genotypes have been set to missing by running the --recode command; however, usually PLINK preserves all genotypes
#when generating a new file (i.e. if one is just reformatting a file, say from text to binary format, it is not necessarily desirable to
#change any of the content; as above, summary statistic and analysis commands do set these genotypes missing automatically still).
#However, if we also add the --set-hh-missing flag, any invalid genotypes will be set to missing in the new file:

plink --file test --set-hh-missing

#which creates the new PED file plink.recode.ped

     1 1 0 0 1 1 A A A A A A A A A A
     2 1 0 0 1 1 C A 0 0 0 0 C A 0 0
     3 1 0 0 2 1 A A A A 0 0 A A A A
     4 1 0 0 2 1 C A C A 0 0 C A 0 0

#In other words, the actual alleles that PLINK pays attention to are shown in bold, all non-bold alleles are ignored.

     1 1 0 0 1   1   A A    A A    A A    A A    A A
     2 1 0 0 1   1   A C    A C    A C    A C    A C
     3 1 0 0 2   1   A A    A A    A A    A A    A A
     4 1 0 0 2   1   A C    A C    A C    A C    A C

#Allele codes
#By default, the minor allele is coded A1 and the major allele is coded A2 (this is used in many output files, e.g. from --freq or --assoc).
#By default this is based on all founders (unless --nonfounders is added) with sex-codes specified (unless --allow-no-sex is added).\
#This coding is applied after any other filters have been applied. It is sometimes desirable to prevent this automatic flipping of\
#A1 and A2 alleles, by use of the --keep-allele-order option. For example, if one wishes to dump the genotype counts by\
#use of the --model command, for two groups of individuals (using the --filter command), this ensures that the same minor allele\
#will always be used in grp1.model as grp2.model (which can facilitate downstream processing of these files, for instance).

plink --bfile --filter pop.dat POP1 --model --keep-allele-order --out pop-1-genotypes
plink --bfile --filter pop.dat POP2 --model --keep-allele-order --out pop-2-genotypes

#That is, for any SNP that happens to have a different minor allele in POP1 versus POP2, the output in the two .model files will\
#still line up in an easy manner.


#Transposed filesets
#Another possible file-format called a transposed fileset, containing two text files: one (TPED) containing SNP and genotype information\
#where one row is a SNP; one (TFAM) containing individual and family information, where one row is an individual.
#The first 4 columns of a TPED file are the same as a standard 4-column MAP file. Then all genotypes are listed for\
#all individuals for each particular SNP on each line.\
#The TFAM file is just the first six columns of a standard PED file. In otherwords, we have just taken the standard PED/MAP\
#file format, but swapped all the genotype information between files, after rotating it 90 degrees.\
#For each, the above example PED/MAP fileset

     <---- normal.ped ---->                  <--- normal.map --->
     1 1 0 0 1  1  A A  G T                  1  snp1   0  5000650
     2 1 0 0 1  1  A C  T G                  1  snp2   0  5000830
     3 1 0 0 1  1  C C  G G
     4 1 0 0 1  2  A C  T T
     5 1 0 0 1  2  C C  G T
     6 1 0 0 1  2  C C  T T

#would be represented as TPED/TFAM files:

     <------------- trans.tped ------------->      <- trans.tfam ->
     1 snp1 0 5000650 A A A C C C A C C C C C      1  1  0  0  1  1
     1 snp2 0 5000830 G T G T G G T T G T T T      2  1  0  0  1  1
                                                   3  1  0  0  1  1
                                                   4  1  0  0  1  2
                                                   5  1  0  0  1  2
                                                   6  1  0  0  1  2

#This kind of format can be convenient to work with when there are very many more SNPs than individuals (i.e. WGAS data).\
#In this case, the TPED file will be very long (as opposed to the PED file being very wide).
#To read a transposed fileset, use the command

plink --tfile mydata

#which implies mydata.tped and mydata.tfam exists; alternatively, if the files are differently named, they can be individually,\
#fully specified:

plink --tped mydata.tped --tfam pedinfo.tx

#You can generate transposed filesets with the --transpose option, described in the data management section


#Long-format filesets
#Another possible file-format called a long-format fileset, containing three text files:
#a LGEN file containing genotypes (5 columns, one row per genotype)
#a MAP file containing SNPs (4 columns, one row per SNP)
#a FAM file containing individuals (6 columns, one row per person)
#The MAP and FAM/PED files are described elsewhere this page. Consider the following example: A MAP file test.map

     1 snp2 0 2
     2 snp4 0 4
     1 snp1 0 1
     1 snp3 0 3
     5 snp5 0 1

#as described above. A FAM file test.fam

     1 1 0 0 1 2
     2 1 0 0 2 2
     2 2 0 0 1 1
     9 1 1 2 0 0

#as described below. Finally, an LGEN file, test.lgen

     1 1 snp1 A A
     1 1 snp2 A C
     1 1 snp3 0 0
     2 1 snp1 A A
     2 1 snp2 A C
     2 1 snp3 0 0
     2 1 snp4 A A
     2 2 snp1 A A
     2 2 snp2 A C
     2 2 snp3 0 0
     2 2 snp4 A A

#The columns in the LGEN file are
     family ID
     individual ID
     snp ID
     allele 1 of this genotype
     allele 2 of this genotype

#Not all entries need to be present in the LGEN file (e.g. snp5 or person 9/1) or snp4 for person 1/1.\
#These genotypes will be set to missing internally. The order also need not be the same in the LGEN file as for the MAP or FAM files.\
#If a genotype is listed more than once, the final version of it will be used.
#LGEN file can be reformatted as a standard PED file using the following command:

plink --lfile test --recode

#which creates these two files: a PED file, plink.recode.map

     1 1 0 0 1  2   A A  A C  0 0  0 0  0 0
     2 1 0 0 2  2   A A  A C  0 0  A A  0 0
     2 2 0 0 1  1   A A  A C  0 0  A A  0 0
     9 1 1 2 0  0   0 0  0 0  0 0  0 0  0 0

#and the MAP file, plink.recode.map (note: it has been put in genomic order)

     1       snp1    0       1
     1       snp2    0       2
     1       snp3    0       3
     2       snp4    0       4
     5       snp5    0       1

#All individuals must be uniquely identified by the combination of the family and individual IDs.
#To read a long-format fileset, use the command

plink --lfile mydata

#which implies mydata.lgen, mydata.map and mydata.map exist
#Currently, you cannot output a fileset in this format in PLINK.

#Additional options for long-format files
#If the LGEN file has specific allele codes, but as TG instead of T G (i.e. no spaces between the two alleles), add the flag

     --compound-genotypes

#It is possible to specify the reference allele with the --reference command when using long-format file input. This might be appropriate,\
#for example, if the data file contains calls for rare variants from a resequencing study. In this case, the majority of alleles will be\
#the reference, and so need not be repeated here. For example, consider this FAM file f1.fam

    1 1 0 0 1 1
    2 1 0 0 1 1
    3 1 0 0 1 1
    4 1 0 0 1 1
    5 1 0 0 1 1
    6 1 0 0 1 1

#and MAP file f1.map

    1       rs0001    0       1000001
    1       rs0002    0       1000002
    1       rs0003    0       1000003

#and LGEN file f1.lgen

    1 1 rs0001 C C
    2 1 rs0001 0 0
    6 1 rs0003 C C
    1 1 rs0002 G T
    4 1 rs0002 T T
    5 1 rs0002 G T

#then

plink --lfile f1 --recode

#would yield a file plink.ped that is as follows:

     1 1 0 0 1 1  C C  G T  0 0
     2 1 0 0 1 1  0 0  0 0  0 0
     3 1 0 0 1 1  0 0  0 0  0 0
     4 1 0 0 1 1  0 0  T T  0 0
     5 1 0 0 1 1  0 0  G T  0 0
     6 1 0 0 1 1  0 0  0 0  C C     

#If the reference all for each variant was set, e.g. with the following command

plink --lfile f1 --reference ref.txt --recode

#and the file ref.txt is

    rs0001 A
    rs0002 G
    rs0009 T

#then the output plink.ped will instead read:

     1 1 0 0 1 1  C C  T G  0 0
     2 1 0 0 1 1  0 0  G G  0 0
     3 1 0 0 1 1  A A  G G  0 0
     4 1 0 0 1 1  A A  T T  0 0
     5 1 0 0 1 1  A A  T G  0 0
     6 1 0 0 1 1  A A  G G  C C

#That is, the non-specified genotypes for the first two SNPs are now homozygous for the reference allele. Note: the word reference is used\
#in the context of the human genome reference allele, rather than for the calculation of an odds ratio.\
#The command to set the latter is --reference-allele {file}
#Also note in this example, that a) when an individual is set as explicitly missing in the LGEN file, they stay missing, b) that when a\
#reference allele is not set, then non-specified genotypes are missing (e.g. the third SNP, rs0003), c) that SNPs in the reference file\
#that are not present in the dataset (e.g. rs0009) are ignored.
#When reading a long-format file, the command

     --allele-count

#when specified along with --reference allows the data to be in the form of the number of non-reference alleles. For example, if input LGEN\
#file were

    1 1 rs0001 0
    2 1 rs0001 1 
    3 1 rs0001 2 
    4 1 rs0001 -1  
    5 1 rs0001 9 
    6 1 rs0001 X 

#this should translate into the first three individuals having the reference homozygote (0 non-reference alleles),\
#the heterozygote (1 non-reference allele) and the non-reference homozygote (2 non-reference alleles). The final three individuals\
#(FID 4 to 6) are all set to missing: this just indicates that any value other than a 0, 1 or 2 under this scheme is set to a missing\
#genotype. If the reference file only contains a single allele for that SNP, then the non-reference allele is coded as whatever is in the\
#reference allele plus a v character appended, e.g. just considering this one SNP:

      1 1 0 0 1 1   A  A
      2 1 0 0 1 1   A  Av
      3 1 0 0 1 1   Av Av
      4 1 0 0 1 1   0  0
      5 1 0 0 1 1   0  0
      6 1 0 0 1 1   0  0

#However, if the reference file contains two alleles, then the second is taken to be the non-reference allele, e.g. if ref.txt is

   rs0001 A  G

#then the output will read
     1 1 0 0 1 1 A A
     2 1 0 0 1 1 A G
     3 1 0 0 1 1 G G
     4 1 0 0 1 1 0 0
     5 1 0 0 1 1 0 0
     6 1 0 0 1 1 0 0

#Binary PED files
#To save space and time, you can make a binary ped file (*.bed). This will store the pedigree/phenotype information in\
#separate file (*.fam) and create an extended MAP file (*.bim)\
#(which contains information about the allele names, which would otherwise be lost in the BED file). To create these files use the command:

plink --file mydata --make-bed

#which creates (by default)

     plink.bed      ( binary file, genotype information )
     plink.fam      ( first six columns of mydata.ped ) 
     plink.bim      ( extended MAP file: two extra cols = allele names)

#The .fam and .bim files are still plain text files: these can be viewed with a standard text editor.\
#Do not try to view the .bed file however: it is a compressed file and you'll only see lots of strange characters on the screen...

#Do not make any changes any of these three files; e.g. setting the position to a negative value will not work to exclude a SNP for\
#binary files
#You can specify a different output root file name (i.e. different to "plink") by using the --out option:

plink --file mydata --out mydata --make-bed

#which will create

     mydata.bed
     mydata.fam
     mydata.bim

#To subsequently load a binary file, just use --bfile instead of --file

plink --bfile mydata

#When creating a binary ped file, the MAF and missingness filters are set to include everybody and all SNPs. If you want to change these,\
#use --maf, --geno, etc, to manually specify these options: for example,

plink --file mydata --make-bed --maf 0.02 --geno 0.1

#Alternate phenotype files
#To specify an alternate phenotype for analysis, i.e. other than the one in the *.ped file\
#(or, if using a binary fileset, the *.fam file), use the --pheno option:

plink --file mydata --pheno pheno.txt

#where pheno.txt is a file that contains 3 columns (one row per individual):

     Family ID
     Individual ID
     Phenotype

#The original PED file must still contain a phenotype in column 6 (even if this is a dummy phenotype, e.g. all missing),\
#unless the --no-pheno flag is given.
#If an individual is in the original file but not listed in the alternate phenotype file, that person's phenotype will be set to\
#missing. If a person is in the alternate phenotype file but not in the original file, that entry will be ignored. \
#The order of the alternate phenotype file need not be the same as for the original file.\
#If the phenotype file contains more than one phenotype, then use the --mpheno N option to specify the Nth phenotype is the\
#one to be used:

plink --file mydata --pheno pheno2.txt --mpheno 4

#where pheno2.txt contains 5 different phenotypes (i.e. 7 columns in total), this command will use the 4th for analysis (phenotype D):

     Family ID
     Individual ID
     Phenotype A
     Phenotype B
     Phenotype C
     Phenotype D
     Phenotype E

#Alternatively, your alternate phenotype file can have a header row, in which case you can use variable names to specify which\
#phenotype to use. If you have a header row, the first two variables must be labelled FID and IID.\
#All subsequent variable names cannot have any whitespace in them. For example,

     FID    IID      qt1   bmi    site  
     F1     1110     2.3   22.22  2     
     F2     2202     34.12 18.23  1     
     ...

#then

plink --file mydata --pheno pheno2.txt --pheno-name bmi --assoc

#will select the second phenotype labelled "bmi", for analysis
#Finally, if there is more than one phenotype, then for basic association tests, it is possible to specify that all phenotypes be\
#tested, sequentially, with the output sent to different files: e.g. if bigpheno.raw contains 10,000 phenotypes, then

plink --bfile mydata --assoc --pheno bigpheno.raw --all-pheno

#will loop over all of these, one at a time testing for association with SNP, generating a lot of output. You might want to use the\
#--pfilter command in this case, to only report results with a p-value less than a certain value, e.g. --pfilter 1e-3.

#Currently, all phenotypes must be numerically coded, including missing values, in the alternate phenotype file.\
#The default missing value is -9, change this with --missing-phenotype, but it must be a numeric value still\
#(in contrast to the main phenotype in the PED/FAM file).

#Creating a new binary phenotype automatically
#To automatically form a one-versus-others binary phenotype (note: binary meaning dichotomous here, rather than a BED/binary-PED file)\
#from a categorical covariate/phenotype file, use the command

plink --bfile mydata --make-pheno site.cov SITE3 --assoc

#which assumes the file

     site.cov

#contains exactly three fields

     Family ID
     Individual ID
     Code from which phenotype is created

#For example, if it were

     A1  1  SITE1
     B1  1  SITE1
     C1  1  SITE2
     D1  1  SITE3
     E1  1  SITE3
     F1  1  SITE4
     G2  1  SITE4

#then the above command would make individuals D1 and E1 as cases and everybody else as controls. However, if individuals present in\
#mydata were not specified in site.cov, then these people would be set to have a missing phenotype.
#An alternate specification is to use the * symbol instead of a value, e.g.

plink --bfile mydata --make-pheno p1.list * --assoc

#which assumes the file

     p1.list

#contains exactly two fields

     Family ID
     Individual ID

#In this case, anybody in the file p1.list would be made a case; all other individuals in mydata but not in p1.list would be set as a\
#control.
#"Loop association": automatically testing each group versus all others
#You may have a categorical factor that groups individuals (e.g. which plate they were genotyped on, or which sample they come from)\
#and want to test whether there are allele frequency differences between each group and all others. This can be accomplished with\
#the --loop-assoc command, e.g.

./plink --bfile mydata --loop-assoc plate.lst --assoc

#The file plate.lst should be in the same format as a cluster file, although it is only allowed to have a single variable\
#(i.e. 3 columns, FID, IID and the cluster variable). If this were

   10001  1   P1
   10002  1   P1
   10003  1   P2
   10004  1   P2
   10005  1   P3   
   10006  1   P3
   ...

#This command would test all P1 individuals against all others, then all P2 individuals against all others, etc.\
#Any of the main single SNP association tests for diseases can be supplied instead of --assoc (e.g. --fisher, --test-missing,\
#--logistic, etc). The output is written to different files for each group, e.g. in the format outputname.{label}.extension

     plink.P1.assoc
     plink.P2.assoc
     plink.P3.assoc

#Covariate files
#Certain PLINK commands support the inclusion of one or more covariates. Note that for stratified analyses,\
#namely using the CMH (--mh) options, the strata are specified using the --within option to define clusters, rather than --covar.
#To load a covariate use the option:

plink --file mydata --covar c.txt

#The covariate file should be formatted in a similar manner to the phenotype file. If an individual is not present in the covariate file,\
#or if the individual has a missing phenotype value (i.e. -9 by default) for the covariate, then that individual is set to\
#missing (i.e. will be excluded from association analysis).
#To select a particular subset of covariates, use one of the following commands, which either use numbers or names\
#(i.e. if a header row exists in the file),

plink --file mydata --covar c.txt --covar-number 2,4-6,8

#or

plink --file mydata --covar c.txt --covar-name AGE,BMI-SMOKE,ALC

#Note that ranges can be used in both cases, with the - hyphen symbol, e.g. if the first row were

     FID IID SITE AGE DOB BMI ETH SMOKE STATUS ALC 

#then both the above commands would have the same effect, i.e. selecting AGE, BMI, ETH, SMOKE, ALC.
#To output a new covariate file, possibly with categorical variables downcoded to binary dummy variables use\
#the --write-covar option as described here

#If the --gxe command is used, that selects only a single covariate, then use the command --mcovar, that works similarly to --mpheno\
#to select which single covariate to use: with the --gxe command, the --covar-name and --covar-number options will not work.

#Not all commands accept covariates, and PLINK will not always give you an error or warning.\
#The basic association (--assoc, --mh, --model, --tdt, --dfam, and --qfam) do not accept covariates, neither do the basic haplotype\
#association methods (--hap-assoc, --hap-tdt). Among the commands that do are --linear, --logistic, --chap and --proxy-glm.\
#Also --gxe accepts a single covariate only (the others listed here accept multiple covariates).
