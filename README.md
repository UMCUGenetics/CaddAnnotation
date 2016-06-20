# CaddAnnotation
The CADD annotation pipeline annotates the contents of a tab-delimted database into an input VCF (in the INFO field).
The chunkscript is used to run the annotation script in parallel, to annotate faster.

## Setup
```
git@github.com:CuppenResearch/CaddAnnotation.git
cd CaddAnnotation
virtualenv env
. env/bin/activate

pip install -r requirements.txt
```

## Usage
```
Usage: python CADDtoAnno_chunkscript.py -db <database_name> -v <input_vcf> -o <output_vcf> -q <queue> -t <maximum runtime> -m <e-mail> -td <temp_directory>

-db     Input database: default is CADD SNVs version 1.3
-v      Input VCF
-o      Output VCF
-q      Queue that job will use: default is all.q
-t      Maximum running time: default is 1 hour
-m      Email used for sending job-information (abortion/killing)
-td     Location of temp directory: default is current working directory
```

#### Example
```
python /path/to/chunkscript/CADDtoAnno_chunkscript.py -v /path/to/inputvcf/test.vcf -o /path/to/outputvcf/testoutput.vcf -t 1 -m example@mail.com
```

Contents of test.vcf:
```
1       22848972        rs61769198      A       C       43.17   PASS    AC=2;AF=1.00;AN=2;DB;DP=4;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=21.59;set=Intersection  GT:AD:DP:GQ:PL  1/1:0,2:2:6:69,6,0      ./.
```

Contents of testoutput.vcf:
```
1       22848972        rs61769198      A       C       43.17   PASS    AC=2;AF=1.00;AN=2;DB;DP=4;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=21.59;set=Intersection;CADDv1.3_GerpN=5.28;CADDv1.3_GerpS=-10.6;CADDv1.3_motifECount=NA;CADDv1.3_cHmmEnhG=0.008;CADDv1.3_dnaMGW=0.73;CADDv1.3_cHmmTssA=0.000;CADDv1.3_EncOCpolIIPVal=NA;CADDv1.3_mutIndex=-6;CADDv1.3_motifEName=NA;CADDv1.3_Segway=TF2;CADDv1.3_Length=0;CADDv1.3_GC=0.56;CADDv1.3_EncOCctcfPVal=NA;CADDv1.3_relcDNApos=NA;CADDv1.3_mirSVR-Score=NA;CADDv1.3_GerpRS=NA;CADDv1.3_motifDist=1.00;CADDv1.3_cHmmReprPCWk=0.000;CADDv1.3_ESP_AF=0.175;CADDv1.3_RawScore=-0.426843;CADDv1.3_cHmmTxFlnk=0.000;CADDv1.3_TFBS=NA;CADDv1.3_#Chrom=1;CADDv1.3_motifEScoreChng=NA;CADDv1.3_CCDS=NA|CCDS224.1;CADDv1.3_Ref=A;CADDv1.3_TG_ASN=NA;CADDv1.3_ESP_AFR=0.094;CADDv1.3_Grantham=NA;CADDv1.3_Exon=NA;CADDv1.3_SIFTcat=NA;CADDv1.3_targetScan=NA;CADDv1.3_mamPhCons=0.000;CADDv1.3_isDerived=TRUE;CADDv1.3_dnaProT=2.56;CADDv1.3_GeneName=ZBTB40-IT1|ZBTB40;CADDv1.3_TFBSPeaks=NA;CADDv1.3_cHmmBivFlnk=0.000;CADDv1.3_EncExp=NA;CADDv1.3_Dst2Splice=NA|-16;CADDv1.3_relProtPos=NA;CADDv1.3_Domain=NA;CADDv1.3_CpG=0.04;CADDv1.3_EncNucleo=2.70;CADDv1.3_EncOCFaireSig=NA;CADDv1.3_EncOCDNasePVal=NA;CADDv1.3_cHmmTssBiv=0.000;CADDv1.3_EncOCDNaseSig=NA;CADDv1.3_EncOCctcfSig=NA;CADDv1.3_TG_AMR=0.230;CADDv1.3_priPhyloP=-0.307;CADDv1.3_TFBSPeaksMax=NA;CADDv1.3_isTv=TRUE;CADDv1.3_Dst2SplType=NA|DONOR;CADDv1.3_protPos=NA;CADDv1.3_dnaHelT=-0.53;CADDv1.3_ConsScore=1|2;CADDv1.3_mapAbility20bp=1;CADDv1.3_nAA=NA;CADDv1.3_cHmmHet=0.000;CADDv1.3_oAA=NA;CADDv1.3_PolyPhenCat=NA;CADDv1.3_TG_EUR=0.200;CADDv1.3_cHmmTx=0.819;CADDv1.3_verPhCons=0.000;CADDv1.3_EncOCFairePVal=NA;CADDv1.3_cHmmEnhBiv=0.000;CADDv1.3_mamPhyloP=-4.713;CADDv1.3_ESP_EUR=0.216;CADDv1.3_scoreSegDup=NA;CADDv1.3_cDNApos=NA;CADDv1.3_bStatistic=497;CADDv1.3_relCDSpos=NA;CADDv1.3_motifEHIPos=NA;CADDv1.3_SIFTval=NA;CADDv1.3_EncOCmycSig=NA;CADDv1.3_EncH3K27Ac=3.76;CADDv1.3_TG_AF=0.120;CADDv1.3_GerpRSpval=NA;CADDv1.3_Anc=A;CADDv1.3_priPhCons=0.008;CADDv1.3_minDistTSE=4884;CADDv1.3_AnnoType=Intergenic|Transcript;CADDv1.3_minDistTSS=70500;CADDv1.3_PHRED=0.333;CADDv1.3_mirSVR-Aln=NA;CADDv1.3_mirSVR-E=NA;CADDv1.3_FeatureID=ENST00000438551|ENST00000404138;CADDv1.3_Intron=NA|17/18;CADDv1.3_dnaRoll=7.81;CADDv1.3_EncOCC=NA;CADDv1.3_cHmmTxWk=0.142;CADDv1.3_verPhyloP=-5.102;CADDv1.3_EncOCpolIISig=NA;CADDv1.3_TG_AFR=0.080;CADDv1.3_cHmmTssAFlnk=0.000;CADDv1.3_Type=SNV;CADDv1.3_tOverlapMotifs=1;CADDv1.3_fitCons=0.099652;CADDv1.3_cHmmReprPC=0.000;CADDv1.3_isKnownVariant=TRUE;CADDv1.3_mapAbility35bp=1;CADDv1.3_EncH3K4Me3=2.44;CADDv1.3_EncH3K4Me1=1.00;CADDv1.3_cHmmQuies=0.024;CADDv1.3_EncOCCombPVal=NA;CADDv1.3_Consequence=DOWNSTREAM|INTRONIC;CADDv1.3_EncOCmycPVal=NA;CADDv1.3_cHmmZnfRpts=0.008;CADDv1.3_ConsDetail=downstream|intron;CADDv1.3_Pos=22848972;CADDv1.3_CDSpos=NA;CADDv1.3_cHmmEnh=0.000;CADDv1.3_Alt=C;CADDv1.3_GeneID=ENSG00000237200|ENSG00000184677;CADDv1.3_PolyPhenVal=NA   GT:AD:DP:GQ:PL  1/1:0,2:2:6:69,6,0     ./.
```


#### Usage without the use of chunking
```
Usage: python CADDtoAnno_processingscript.py -db <database_name> -v <input_vcf> -o <output_vcf>

-db     Input database: default is CADD SNVs version 1.3
-v      Input VCF
-o      Output VCF
```

#### Example
```
python /path/to/processingscript/CADDtoAnno_processingscript.py -v /path/to/inputvcf/test.vcf -o /path/to/outputvcf/testoutput.vcf
```

# CADD report (runs locally)
Creates a CADD report that counts the total number of variants and the number of variants in functional categories.
Creates for each category a boxplot with the distribution of the CADD PHRED-score.

## Setup
```
cd CaddAnnotation
virtualenv env
. env/bin/activate
```

## Additional requirements for using the report script
```
reportlab==3.3.0
pymongo==3.2.2
matplotlib==1.5.1
pyPdf==1.13
```

Requires the use of a running [vcf-explorer](https://github.com/CuppenResearch/vcf-explorer).

## Usage
```
python CADD_report.py -sample <sample_name>
```

#### Example
```
python CADD_report.py -sample AB12345
```
