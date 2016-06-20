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
