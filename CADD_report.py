from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, PageBreak, ListFlowable, ListItem, Indenter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from bson.son import SON
from matplotlib.backends.backend_pdf import PdfPages
from pyPdf import PdfFileWriter, PdfFileReader
from reportlab.lib import colors

import pymongo
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import StringIO
import argparse
import re
import datetime, time
import os
import math

parser = argparse.ArgumentParser()
parser.add_argument("-sample", help="give the sample name of the patient")
args = parser.parse_args()

PAGE_HEIGHT=defaultPageSize[1]
PAGE_WIDTH=defaultPageSize[0]
styles = getSampleStyleSheet()

Title = "Report for CADD annotation pipeline"
pageinfo = "CADD report / {}".format(str(datetime.date.today()))

connection = pymongo.MongoClient("mongodb://localhost")
db = connection.vcf_explorer

"""
This function calculates basic statistics based on the VCF and returns each statistic.
"""
def stats():
    vcf = ('VCF name: {}'.format (db.vcfs.distinct('name')[0]))
    samples = ('Sample names: {}'.format (", ".join(db.variants.distinct("samples.sample"))))
    variants = ('Amount of variants for patient ({}): {:,}'.format (args.sample, db.variants.find( {"samples.sample": args.sample } ).count()))
    snvs = ('Amount of SNVs for patient: {:,}'.format (db.variants.find( {"samples.sample": args.sample, "$where": "this.alt.length == 1 && this.ref.length == 1" } ).count()))
    indels = ('Amount of InDels for patient: {:,}'.format (db.variants.find( {"samples.sample": args.sample, "$where": "this.alt.length != this.ref.length"} ).count()))

    return vcf, samples, variants, snvs, indels

"""
This function counts for each category the occurences of each group.
It returns the table title and table data.
"""
def countCategories(search_one):
    # First only variants found in the patient are selected
    # For the category the groups are selected and a UniqueSet is created (so group can only occur once for each variant)
    # Then these are grouped and counted and lastly sorted from high to low
    pipeline = [
            {"$match": {"samples.sample": args.sample}},
            {"$unwind": search_one},
            {"$group": {"_id": "$_id", "UniqueSet": {"$addToSet": search_one}}},
            {"$unwind": "$UniqueSet"},
            {"$group": {"_id": "$UniqueSet", "count": {"$sum": 1}}},
            {"$sort": SON([("_id", 1), ("count", -1)])}
            ]
    count = list(db.variants.aggregate(pipeline, allowDiskUse=True))
    segway_descr = {"C0" : "CTCF (strong)", "C1" : "CTCF (weak)", "D" : "dead zone", "E/GM" : "enhancer/gene middle", "F0" : "FAIRE only", "F1" : "FAIRE only", "GE0" : "gene body (end)", "GE1" : "gene end",
    "GE2" : "gene end", "GM0" : "gene body (middle)", "GM1" : "gene body (middle)", "GS" : "gene body (start)", "H3K9me1" : "H3K9me1 only", "L0" : "low zone", "L1" : "low zone", "NA" : "No segmentations found",
    "R0" : "repression", "R1" : "repression", "R2" : "repression", "R3" : "repression", "R4" : "repression", "R5": "repression", "TF0" : "transcription factor activity", "TF1" : "transcription factor activity", "TF2" : "transcription factor activity",
    "TSS" : "transcription start site"}

    # Output data (Category : Count) is added to the list 'data', which can be converted to a table
    data = []
    index = 0
    for element in count:
        if search_one == "$info.CADD.CADDv1-3_Segway":
            data.append([element["_id"], segway_descr[element["_id"]], str(element["count"])])
            index += 1
        else:
            data.append([element["_id"], str(element["count"])])

    cat_title = ("Amount of variants per " + search_one.split("_")[1] + " category:")

    return data, cat_title

"""
This function counts the amount of variants with a PHRED score above 25 and returns this score.
"""
def countHighPhred():
    # First only variants found in the patient are selected
    # Then all variants with a PHRED above 25 were selected and counted
    pipeline = [
        {"$match": {"samples.sample": args.sample}},
        {"$match": {"info.CADD.CADDv1-3_PHRED": {"$gt": 25}}},
        {"$group": {"_id": "PHRED", "count": {"$sum": 1}}}
        ]
    count = list(db.variants.aggregate(pipeline, allowDiskUse=True))

    # If there is at least 1 variant above 25 found the counted variants are used as vals
    # Else, vals is zero
    if len(count) != 0:
        for element in count:
            vals = (element["count"])
    else:
        vals = 0

    return ("Amount of variants with PHRED score above 25: %s") % (vals)

"""
This function counts the amount of variants with a proportion of at least 0.75 of cell types in a certain chromatin state.
It returns the table title and the table data.
"""
def countCHMM():
    chmm_cats = ["info.CADD.CADDv1-3_cHmmTssA", "info.CADD.CADDv1-3_cHmmTssAFlnk", "info.CADD.CADDv1-3_cHmmTxFlnk", "info.CADD.CADDv1-3_cHmmTx", "info.CADD.CADDv1-3_cHmmTxWk",          "info.CADD.CADDv1-3_cHmmEnhG", "info.CADD.CADDv1-3_cHmmEnh", "info.CADD.CADDv1-3_cHmmZnfRpts","info.CADD.CADDv1-3_cHmmHet", "info.CADD.CADDv1-3_cHmmTssBiv", "info.CADD.CADDv1-3_cHmmBivFlnk", "info.CADD.CADDv1-3_cHmmEnhBiv", "info.CADD.CADDv1-3_cHmmReprPC", "info.CADD.CADDv1-3_cHmmReprPCWk", "info.CADD.CADDv1-3_cHmmQuies"]

    chmm_descr = {"cHmmTssA" : "Active Transcription Start Site", "cHmmTssAFlnk" : "Flanking Active Transcription Start Site", "cHmmTxFlnk" : "Transcription at gene 5' and 3'",
                  "cHmmTx" : "Strong transcription", "cHmmTxWk" : "Weak transcription", "cHmmEnhG" : "Genic enhancers", "cHmmEnh" : "Enhancers", "cHmmZnfRpts" : "ZNF genes & repeats", "cHmmHet" : "Heterochromatin", "cHmmTssBiv" : "Bivalent/Poised Transcription Start Site", "cHmmBivFlnk" : "Flanking Bivalent Transcription Start Site/Enhancer",
                  "cHmmEnhBiv" : "Bivalent Enhancer", "cHmmReprPC" : "Repressed Polycomb", "cHmmReprPCWk" : "Weak Repressed Polycomb", "cHmmQuies" : "Quiescent/Low"}

    data = []
    # First only variants found in the patient were selected
    # Then the variant is selected when the chromatin state-score is above 0.75
    # Lastly, the variants are counted and grouped from high to low.
    index = 0
    for state in chmm_cats:
        pipeline = [
            {"$match": {"samples.sample": args.sample}},
            {"$match": {state: {"$gt": 0.75}}},
            {"$group": {"_id": state, "count": {"$sum": 1}}}
            ]
        count = list(db.variants.aggregate(pipeline, allowDiskUse=True))

        for element in count:
            data.append([element["_id"].split("_")[1], chmm_descr[str(element["_id"].split("_")[1])], str(element["count"])])
            index += 1

    # Output data (Category : Count) is added to the list 'data', which can be converted to a table
    chmm_title = ("Amount of variants for each chromatin state:")

    return data, chmm_title

"""
This function counts the amount of variant for the patient in each biotype.
"""
def customEnsembl():
    data = []
    for samples in (db.variants.find( {"samples.sample": args.sample } )):
        try:
            custom_line =  samples["info"]["ANN"]["custom"]
            for feature in custom_line:
                custom = re.sub('[A-Z]{4}[0-9]{11}_', '', feature)
                data.append(custom)
        except KeyError:
            pass

    count_dict = {}
    for item in data:
        if item in count_dict.keys():
            count_dict[item] += 1
        else:
            count_dict[item] = 1

    count_list = []
    for x in range (0, len(count_dict)):
        count_list.append([count_dict.keys()[x], str(count_dict.values()[x])])

    return count_list

"""
This function creates boxplots based on a .csv file created with mongodb.
This file contains all ConsDetail categories and PHRED scores.
The boxplots represent the distrubtion of the PHRED score for each category.
"""
def plots():
    plt.style.use('ggplot')
    # Transfer data to a float
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    val = {}
    bashcsv = open("./csv_data.sh", "w")
    bashcsv.write("cd ../mongodb-linux-x86_64-3.2.4/bin/")
    bashcsv.write("\n")
    bashcsv.write("./mongoexport --db vcf_explorer --collection variants --type=csv --fields info.CADD.CADDv1-3_ConsDetail,info.CADD.CADDv1-3_PHRED --out /home/sverduin/data/mongoexport.csv")
    bashcsv.close()

    os.system("bash ./csv_data.sh")
    time.sleep(10)
    with open("/home/sverduin/data/mongoexport.csv", 'r') as data:
        for line in data:
            if not line.startswith("i"):
                # Strip data from all characters except for the actual data
                values = filter(None, line.strip("\n").replace('"', '').replace('[', '').replace(']', '').split(","))
                # Create values for x-axis
                x_values = list(set(values))
                # Add values to val dictionary
                for x in range(0, len(x_values)):
                    if isfloat(x_values[x]):
                        pass
                    elif x_values[x] in val :
                        val[x_values[x]].append(float(values[-1]))
                    else:
                        val[x_values[x]] = [float(values[-1])]

    # Create output file for boxplot
    pp = PdfPages('Temp_boxplot.pdf')

    # lenkeys = float(len(val.keys()))
    num_plots = int(math.ceil(float(len(val.keys()))/5.))
    # Create plot for each category; 5 for each line
    fig, axes = plt.subplots(num_plots, ncols=1, figsize=(20, 20))

    # Creates rows for boxplot
    for x in range(0, num_plots):
        all_data = val.values()[x*5:(x+1)*5]

        # For every 5 categories create a boxplot per line
        bp = axes[x].boxplot(all_data, 0, 'gD', patch_artist=True)
        for box in bp['boxes']:
            # change outline color
            box.set( color='#000000', linewidth=2)
            # change fill color
            box.set( facecolor = '#5882FA' )

        # Set options for boxplot
        axes[x].yaxis.grid(True)
        axes[x].set_xticks([y+1 for y in range(len(all_data))])
        axes[x].set_ylabel('CADD PHRED Score')
        axes[x].set_ylim([0,50])

        # Create plots
        plt.setp(axes[x], xticks=[y+1 for y in range(len(all_data))],
                 xticklabels= val.keys()[x*5:(x+1)*5])

    # Set title and save plot to pdf
    fig.suptitle("Boxplots showing the distribution of the PHRED score for each consequence category in dataset", y=0.92, fontsize=17)
    pp.savefig()
    plt.close()
    pp.close()

"""
Function that creates the first page of final pdf.
"""
def myFirstPage(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Bold',16)
    canvas.drawCentredString(PAGE_WIDTH/2.0, PAGE_HEIGHT-108, Title)
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch,"Page 1 / %s" % pageinfo)
    canvas.restoreState()

"""
Function that creates any further pages for the final pdf.
"""
def myLaterPages(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch,"Page %d / %s" % (doc.page, pageinfo))
    canvas.restoreState()

"""
This function creates the final pdf.
In this function all the writing to pdf will be done.
"""
def go():
    customEnsembl()
    # Options for document
    doc = SimpleDocTemplate("Temp_report.pdf")
    style = styles["Normal"]
    elements = [Spacer(1,1*inch)]

    # Paragraphs are created of the vcf statistics
    vcf, samples, variants, snvs, indels = stats()
    vcf_story = Paragraph(vcf, style)
    samples_story = Paragraph(samples, style)
    variants_story = Paragraph(variants, style)
    snvs_story = Paragraph(snvs, style)
    indels_story = Paragraph(indels, style)
    phred_story = Paragraph(countHighPhred(), style)

    cons_data, cons_title = countCategories("$info.CADD.CADDv1-3_ConsDetail")
    custom_data = customEnsembl()
    coding = ["splice_synonymous", "synonymous", "missense_splice_NMD", "missense_splice", "stop_gained_NMD", "stop_gained_splice",
              "stop_retained", "missense", "incomplete_terminal_codon_coding_sequence", "stop_lost", "splice_synonymous_NMD",
              "missense_NMD", "coding_sequence", "stop_gained", "synonymous_NMD", "stop_lost_NMD", "initiator_codon", "initiator_codon_splice",
              "protein_coding", "antisense"]
    non_coding = ["splice_intron_nc", "intron_nc", "5_prime_UTR", "splice_5_prime_UTR_NMD", "3_prime_UTR_NMD", "non_coding_exon_nc",
                  "splice_intron", "nc", "5_prime_UTR_NMD", "splice_5_prime_UTR", "splice_donor_nc", "intron_NMD", "splice_donor",
                  "downstream", "splice_donor_non_coding_exon_nc", "upstream", "intergenic", "splice_non_coding_exon_nc", "3_prime_UTR",
                  "intron", "splice_donor_NMD", "splice_intron_NMD", "splice_acceptor_nc", "splice_acceptor_NMD", "splice_acceptor",
                  "splice_3_prime_UTR", "pseudogene", "polymorphic_pseudogene", "sense_intronic", "sense_overlapping"]
    regulatory = ["regulatory"]
    rna = ["mature_miRNA", "lincRNA", "miRNA", "snRNA", "rRNA", "misc_RNA", "snoRNA", "processed_transcript", "3prime_overlapping_ncrna"]

    coding_items = ""
    nc_items = ""
    regulatory_items = ""
    rna_items = ""
    others = ""

    for category in cons_data, custom_data:
        for item in category:
            if item[0] in coding:
                coding_items += str(": ".join(item))
                coding_items += "<br/>"
            elif item[0] in non_coding:
                nc_items += str(": ".join(item))
                nc_items += "<br/>"
            elif item[0] in regulatory:
                regulatory_items += str(": ".join(item))
                regulatory_items += "<br/>"
            elif item[0] in rna:
                rna_items += str(": ".join(item))
                rna_items += "<br/>"
            else:
                others += str(": ".join(item))
                others += "<br/>"

    coding_story = Paragraph(coding_items, style)
    nc_story = Paragraph(nc_items, style)
    regulatory_story = Paragraph(regulatory_items, style)
    rna_story = Paragraph(rna_items, style)

    if others != "":
        other_story = Paragraph(others, style)

    # Table title and table data is retrieved
    segway_data, segway_title = countCategories("$info.CADD.CADDv1-3_Segway")
    chmm_data, chmm_title = countCHMM()

    # Table title are converted to a paragraph
    cons_story = Paragraph('<b> {} </b>'.format(cons_title), style)
    segway_story = Paragraph('<b> {} </b>'.format(segway_title), style)
    chmm_story = Paragraph('<b> {} </b>'.format(chmm_title), style)

    # Options for table (makes sure text fits in boxes by wrapping)
    s = getSampleStyleSheet()
    s = s["BodyText"]
    s.wordWrap = 'CJK'

    # Creates for each element in the table data a paragraph
    cons_data2 = [[Paragraph(cell, s) for cell in row] for row in cons_data]
    segway_data2 = [[Paragraph(cell, s) for cell in row] for row in segway_data]
    chmm_data2 = [[Paragraph(cell, s) for cell in row] for row in chmm_data]

    # Creates the layout of the table
    table_style = [
        ('GRID', (0,0), (-1,-1), 1, colors.black),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('LEFTPADDING', (0,0), (-1,-1), 3),
        ('RIGHTPADDING', (0,0), (-1,-1), 3),
        ('FONTSIZE', (0,0), (-1,-1), 10),
        ('FONTNAME', (0,0), (-1,0), 'Times-Bold'),
        ('TEXTCOLOR', (0,0), (-1, 0), colors.red)
    ]

    # Creates for each table data an actual table with chosen table style
    cons_table=Table(cons_data2)
    cons_table.setStyle(table_style)
    segway_table=Table(segway_data2)
    segway_table.setStyle(table_style)
    chmm_table=Table(chmm_data2)
    chmm_table.setStyle(table_style)

    # Paragraphs are added to elements list
    elements.append(vcf_story)
    elements.append(samples_story)
    elements.append(variants_story)
    elements.append(snvs_story)
    elements.append(indels_story)
    elements.append(phred_story)

    elements.append(PageBreak())
    elements.append(Paragraph("<b> Hierarchy of the distribution of the variants </b>", style))
    elements.append(Paragraph("Coding:", style))
    elements.append(Indenter(36,0))
    elements.append(coding_story)
    elements.append(Indenter(-36,0))
    elements.append(Paragraph("Noncoding:", style))
    elements.append(Indenter(36,0))
    elements.append(nc_story)
    elements.append(Indenter(-36,0))
    elements.append(Paragraph("Regulatory:", style))
    elements.append(Indenter(36,0))
    elements.append(regulatory_story)
    elements.append(Indenter(36,0))

    for x in chmm_data:
        elements.append(Paragraph(('{}: {}'.format(x[1],x[2])), style))
    elements.append(Indenter(-72,0))

    elements.append(Paragraph("RNAs:", style))
    elements.append(Indenter(36,0))
    elements.append(rna_story)
    elements.append(Indenter(-36,0))

    if 'other_story' in locals():
        elements.append(Paragraph("Not defined by category:", style))
        elements.append(Indenter(36,0))
        elements.append(other_story)
        elements.append(Indenter(-36,0))

    # Table titles and tables itself are added to elements list
    elements.append(PageBreak())
    elements.append(cons_story)
    elements.append(cons_table)
    elements.append(PageBreak())
    elements.append(segway_story)
    elements.append(segway_table)
    elements.append(PageBreak())
    elements.append(chmm_story)
    elements.append(chmm_table)

    # Document is build based on the elements list and firstpage, laterpages functions
    doc.build(elements, onFirstPage=myFirstPage, onLaterPages=myLaterPages)

    """
    Merges the elements pdf with the violin pdf
    """
    # Violinplot document is created
    plots()
    packet = StringIO.StringIO()
    packet.seek(0)
    # Read pdf with elements list
    elements_pdf = PdfFileReader(file("Temp_report.pdf", "rb"))
    # Read pdf with boxplot
    violin_pdf = PdfFileReader(file("Temp_boxplot.pdf", "rb"))
    # Create new output file
    output = PdfFileWriter()
    # Add the "watermark" (which is the new pdf) on the existing page
    violin_page = violin_pdf.getPage(0)
    # For each page in elements list pdf, add to final pdf
    for num in range(0, elements_pdf.getNumPages()):
        output.addPage(elements_pdf.getPage(num))
    # Add boxplot pdf to new pdf
    output.addPage(violin_page)
    # Finally, write "output" to a real pdf file
    outputStream = file("CADD_Report.pdf", "wb")
    output.write(outputStream)
    outputStream.close()

if __name__ == "__main__":
    go()
