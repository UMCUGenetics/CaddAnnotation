import gzip
import tabix
import subprocess
from itertools import islice
from helper import deep_update
import getopt, sys, argparse, os

"""
SHIRA VERDUIN
This program reads in a VCF and (CADD)database and annotates the VCF using the database.
"""

"""
This function creates the annotations for the VCF and writes this to file.
Uses the CADD_dictionary to combine keys and values as annotations (key=value;)
This CADD annotations are appended to INFO field and then written to output VCF.
"""
def createOutputFile(keys, vals, vcf_content, output_file):
    # Adds CADD annoations to the end of the INFO field.
    for pos in range(0, len(keys)):
	vcf_content[7] += ";CADDv1.3_{keys}={vals}".format(
	keys = keys[pos],
	vals = ",".join(vals[pos]))
    vcf_line = ""
    # Adds full INFO field to VCF.
    for header_length in range(0, len(vcf_content)):
	vcf_line += (vcf_content[header_length] + "\t")
    output_file.write(vcf_line.strip("\t") + "\n")


# CADD database is transferred to a dictionary containing CADD-info as key/value pairs.
"""
This function creates the header CADD content and writes this to output file.
Then, based on the variants in the chunk-VCF it will find the matching CADD annotations (if chr:pos:alt:ref is equal).
Output will be written to a database and this is returned to the createOutputFile function, where output is created.
"""
def CADDtoAnno(inputdb, inputvcf, outputfile):
    output_file = open(outputfile,  "w")
    caddtb = tabix.open(inputdb)

    # Write original header to output file.
    output_file.write(subprocess.check_output("cat {} | grep ^'##'".format(inputvcf), shell=True))

    # Define the keys of CADD header for VCF-header CADD info.
    dbkeys = ""
    with gzip.open(inputdb)  as db_header:
	for db_line in islice(db_header, 1, 2):
	    dbkeys += db_line.strip("\n")

    # Create the CADD-header information and write to output file.
    cadd_valuetype = ['String', 'Integer', 'String', 'String', 'String', 'String', 'Integer', 'String', 'String', 'String', 'String', 'Integer', 'String', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'String', 'Float', 'Float', 'Float', 'String', 'String', 'Float', 'Float', 'Float', 'Float', 'String', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'Integer','Integer', 'String', 'String', 'String', 'String', 'Float', 'Float', 'Float', 'Float', 'Float', 'Float', 'String', 'Float', 'String', 'String', 'String', 'String', 'String', 'Float', 'String', 'Float', 'String', 'Float', 'Float', 'Float']
    caddheader = []
    for valuetype in range (0, len(cadd_valuetype)):
	caddid = dbkeys.split("\t")
	output_file.write("##INFO=<ID=CADDv1.3_{0},Number=A,Type={1},Description={0}>\n".format(caddid[valuetype], cadd_valuetype[valuetype]))
	caddheader.append(caddid[valuetype])

    # Create main VCF header and write to output file.
    vcf_header = (subprocess.check_output("cat {} | grep ^'#CHROM'".format(inputvcf),shell=True).strip("\n"))
    comp_header = ""
    for header_val in vcf_header.split():
	comp_header += (header_val + "\t")
    output_file.write(comp_header[:-1] + "\n")

    # Open VCF.
    try:
	# For .vcf
	vcf_file = open(inputvcf, "r")
    except IOError:
	# For .vcf.gz
	vcf_file = gzip.open(inputvcf, "rb")
    except:
	print "VCF cannot be opened."

    with vcf_file:
	for vcf_line in vcf_file:
	    # Use only actual content (so no header content)
	    if vcf_line.startswith("#"):
		pass
	    else:
		vcf_line = vcf_line.strip("\n")
		vcf_content = vcf_line.split("\t")

		# Find CADD annotations using tabix query with variant information (chr:pos-pos).
		chrom = vcf_content[0]
		pos = vcf_content[1]
		tabix_pos = "{0}:{1}-{1}".format(chrom, pos)
		tabix_output = list(caddtb.querys(tabix_pos))

		# Find correct allelic CADD annotations (also for multiple alternative alleles).
		cadd_dict = {}
		alts = vcf_content[4].split(",")
		for alt in alts:
		    # Create CADD alt-dict
		    cadd_alt_dict = {}
		    for cadd_line in tabix_output:
			# If VCF ref:alt is equal to CADD ref:alt, continue
			if "{}:{}".format(vcf_content[3].strip(" "), alt) == "{}:{}".format(cadd_line[2], cadd_line[4]):
			    caddrange = range(0, len(caddheader))
			    for val in caddrange:
				# If key exists and the value is not equal to current value, append the value to line using a pipe;
				# in this way no identical (repeating) annotations occur.
				# If key does not exist, create new key/value pair.
				if caddheader[val] in cadd_alt_dict.keys():
				    if cadd_alt_dict[caddheader[val]][0] not in cadd_line[val]:
					cadd_alt_dict[caddheader[val]][0] += ("|" + cadd_line[val].replace(",","_"))
				else:
				    cadd_alt_dict[caddheader[val]] = [cadd_line[val].replace(",","_")]

		    # Update cadd_dictionary with alt_dictionary. Use keys and values for createOutputFile function.
		    deep_update(cadd_dict, cadd_alt_dict)
		keys = cadd_dict.keys()
		vals = cadd_dict.values()
		createOutputFile(keys, vals, vcf_content, output_file)
    output_file.close()
    vcf_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CADD annotation pipeline')
    # Arguments that are required and/or possible to define with running the script.
    parser.add_argument("-db", "--database_file", help="Input database used for annotations", default="hpc/cog_bioinf/common_dbs/CADD/whole_genome_SNVs_inclAnno.tsv.gz")
    parser.add_argument("-v", "--vcf_file", help="Input VCF-chunk to annotate", required=True)
    parser.add_argument("-o", "--output_file", help="Output VCF-file with annotations", required=True)
    args = parser.parse_args()
    
    CADDtoAnno(args.database_file, args.vcf_file, args.output_file)