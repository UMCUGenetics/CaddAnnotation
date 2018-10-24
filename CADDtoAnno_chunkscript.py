import argparse
import subprocess
import os
import sys

import helper


def chunkVCF(inputdb, inputvcf, outputfile, cadd_fields, queue, time, email, temp):
    """Read in the full unannotated VCF and chunks this in chunks of 25.000 variants.

    For every chunk the original header is added to the variants.
    The old chunks are then removed.
    """
    # Define temporary path were chunks will be stored and create directory.
    chunk_path = temp+"/chunks"
    os.system("mkdir {0}".format(temp))
    os.system("rm -rf {0}; mkdir -p {0}".format(chunk_path))
    # Determine length of header and actual header content.
    headerlines = subprocess.check_output("cat {} | egrep -c ^'#'".format(inputvcf), shell=True).split()
    header = subprocess.check_output("cat {} | egrep ^'##'".format(inputvcf), shell=True).strip("\n")
    # Chunk the VCF without its header.
    os.system("tail -n +{header} {vcf} | split -l 25000 - {path}/chunk_".format(
        header=(int(headerlines[0])+1), vcf=inputvcf, path=chunk_path)
    )

    # For every chunk a new file is made with the header content.
    for files in os.listdir(chunk_path):
        newfile = open(chunk_path + "/header{}.vcf".format(files), "w")
        chunks = open("{}/{}".format(chunk_path, files), "r")
        newfile.write(header + "\n")
        vcf_header = (subprocess.check_output("cat {} | grep ^'#CHROM'".format(inputvcf), shell=True).strip("\n"))

        # Write header content to file.
        comp_header = ""
        for header_val in vcf_header.split():
            comp_header += (header_val + "\t")
        newfile.write(comp_header[:-1])

        # Write variant content to file.
        for variant in chunks:
            newfile.write("\n" + variant.strip("\n"))

        # Close files.
        newfile.close()
        chunks.close()

    # Remove old chunks.
    os.system("rm -r {}/chunk*".format(chunk_path))
    createShells(inputdb, chunk_path, outputfile, cadd_fields, queue, email, time, temp)


def createShells(inputdb, chunk_path, outputfile, cadd_fields, queue, email, time, temp):
    """Create shell scripts that will run the processing part of CADDtoAnno for every chunk.

    Also creates a concat shell script that will merge the chunks after finishing annotating.
    """
    job_names = []
    concat_names = []
    script_path = os.path.dirname(os.path.realpath(__file__))

    # For every chunk a shell script is made that will run the processing script for this chunk.
    for files in os.listdir(chunk_path):
        if files.startswith("headerchunk") and files.endswith(".vcf"):
            chunks_shell = ("{0}/{1}.sh".format(chunk_path, files))
            sh_script = open(chunks_shell, "w")
            sh_script.write("#!/bin/bash\n. {}/env/bin/activate\n".format(script_path))
            sh_script.write("time python {path}/CADDtoAnno_processingscript.py -db {db} -v {chunkpath}/{chunkfile} -o {chunkpath}/{vcfname} -f {cadd_fields}".format(
                path=script_path,
                db=inputdb,
                chunkpath=chunk_path,
                chunkfile=files,
                vcfname=files.split(".")[0]+"_CADDv1.3.vcf",
                cadd_fields=' '.join(cadd_fields)
            ))
            sh_script.close()

            # Creates a list with the job names (files) and conat names (including path).
            chunknames = "chunk{}".format(files.split("_")[1])
            job_names.append(chunknames)
            concat_names.append("{0}/{1}".format(chunk_path, files.split(".")[0]+"_CADDv1.3.vcf"))

            # Submit the chunks to cluster.
            os.system("qsub -cwd -l h_rt={time}:0:0 -m a -M {email} -N {name} {script_shell}".format(
                time=time,
                email=email,
                name=chunknames,
                script_shell=chunks_shell)
            )

    # Create concat script.
    concat_script = open("{}/concat.sh".format(temp), "w")
    concat_script.write("#!/bin/bash\n. {}/env/bin/activate\n".format(script_path))
    concat_script.write("module load vcfbcf/vcftools/c7a7337")
    concat_script.write("\nvcf-concat {concat} > {output}".format(
        concat=" ".join(sorted(concat_names)),
        output=outputfile))
    concat_script.write("\nrm -r {}".format(temp))
    concat_script.close()

    # Submit concat script to cluster, but only when all chunk scripts are finished (hold).
    os.system("qsub -hold_jid {job} -cwd -l h_rt=1:0:0 -m a -M {email} {script_shell}".format(
        job=",".join(job_names),
        email=email,
        script_shell=temp+"/concat.sh"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CADD annotation pipeline')
    # Arguments that are required and/or possible to define with running the script
    parser.add_argument("-db", "--database_file", help="Input database used for annotations. Default is CADDv1.3_SNVs", default="/hpc/cog_bioinf/common_dbs/CADD/whole_genome_SNVs_inclAnno.tsv.gz")
    parser.add_argument("-v", "--vcf_file", help="Input name VCF-file to annotate", required=True)
    parser.add_argument("-o", "--output_file", help="Output name of VCF-file", required=True)
    parser.add_argument('-f', "--fields", help="CADD fields used to annotate vcf.", nargs='+', default=helper.cadd_columns.keys())
    parser.add_argument("-q", "--queue", help="Queue that job will use. Default is all.q", default="all.q")
    parser.add_argument("-t", "--time", help="Maximum running time (in hours) of job. Default is one hour", default="1")
    parser.add_argument("-m", "--email", help="E-mail used for sending job-info", required=True)
    parser.add_argument("-td", "--temp_directory", help="Temporary directory where chunks will be placed. Default is your current working directory", default="{}/Temp".format(os.getcwd()))
    args = parser.parse_args()

    for field in args.fields:
        if field not in helper.cadd_columns.keys():
            sys.exit('{0} is not a valid CADD column'.format(field))

    chunkVCF(args.database_file, args.vcf_file, args.output_file, args.fields, args.queue, args.time, args.email, args.temp_directory)
