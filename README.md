# CaddAnnotation

## Setup
```
cd vcf-explorer
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

## Example
```
python /path/to/chunkscript/CADDtoAnno_chunkscript.py -v /path/to/inputvcf/test.vcf -o /path/to/outputvcf/testoutput.vcf -t 1 -m example@mail.com
```
