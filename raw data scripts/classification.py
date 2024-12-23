import sys

#Extract necessary info and save in csv file
with open ("3_assembly/1_stringtie_merged.gtf", "r") as r:
    with open ("5_classification/0_transcript_info_from_gtf.csv", "w") as w:
        w.write("Chromosome,type,end_1,end_2,strand,gene_id,transcript_id\n")
        r = r.readlines()
        for x in r[2:]:
            x = x.rstrip().split("\t")
            if len(x) < 9:
                print (x)
                break
            else:
                if len(x[8].split('"')) < 4:
                    print (x[8])
                    break
                else:
                    w.write(x[0] + "," + x[2] + "," + x[3] + "," + x[4] + "," + x[6] + "," + x[8].split('"')[1] + "," + x[8].split('"')[3] + "\n")


#Test the assumtion that same gene id should come from same locus in chromosome
test = {}
with open ("5_classification/0_transcript_info_from_gtf.csv", "r") as r:
    r = r.readlines()
    for x in r[1:]:
        x = x.rstrip().split(",")
        if x[1] == "transcript":
            if x[5] not in test.keys():
                test[x[5]] = [x[0],x[4],int(x[2]),int(x[3])]
            else:
                if test[x[5]][0] != x[0] or test[x[5]][1] != x[4]:
                    print (x[6])
                    sys.exit("Same gene but different chromosome or strand")
                else:
                    if ((min(int(x[3]),test[x[5]][3]))-(max(int(x[2]),test[x[5]][2]))) <= 0:
                        print (x[6], x[2], x[3], test[x[5]][2], test[x[5]][3])
                        sys.exit("Same gene but no overlap")
                    else:
                        test[x[5]][2] = min(int(x[2]),test[x[5]][2])
                        test[x[5]][3] = max(int(x[3]),test[x[5]][3])


#Classify transcrips into overlap, antisense, and intergenic groups

ref_genes_by_chrom = {}
novel_genes_by_chrom = {}

with open("5_classification/0_transcript_info_from_gtf.csv", "r") as r:
    r = r.readlines()
    for x in r[1:]:
        x = x.rstrip().split(",")
        if x[1] == "transcript":
            if 'PF3D7' in x[6]:
                if x[0] not in ref_genes_by_chrom:
                    ref_genes_by_chrom[x[0]] = []
                ref_genes_by_chrom[x[0]].append((int(x[2]), int(x[3]), x[6], x[4]))
            else:
                if x[0] not in novel_genes_by_chrom:
                    novel_genes_by_chrom[x[0]] = []
                novel_genes_by_chrom[x[0]].append((int(x[2]), int(x[3]), x[6], x[4]))

overlap_genes = []
ov = []
antisense_genes = []
an = []
intergenic_genes = []

for chrom in novel_genes_by_chrom:
    if chrom in ref_genes_by_chrom:
        for start, end, gene_name, strand in novel_genes_by_chrom[chrom]:
            overlap_found = False
            for ref_start, ref_end, ref_gene_name, ref_strand in ref_genes_by_chrom[chrom]:
                if start <= ref_end and end >= ref_start:
                    if strand == ref_strand:
                        overlap_genes.append(gene_name)
                        ov.append(ref_gene_name)
                        overlap_found = True
                        break
                    else:
                        antisense_genes.append(gene_name)
                        an.append(ref_gene_name)
                        overlap_found = True
                        break
            if not overlap_found:
                intergenic_genes.append(gene_name)
    else:
        intergenic_genes.extend(gene[2] for gene in novel_genes_by_chrom[chrom])

with open ("5_classification/1_transcript_classification.csv", "w") as w:
    # Write the header row
    w.write("Overlap,Ref_Overlap,Antisense,Ref_Antisense,Intergenic\n")
    
    # Write the data rows
    max_rows = max(len(overlap_genes), len(antisense_genes), len(intergenic_genes))
    for i in range(max_rows):
        row = ""
        if i < len(overlap_genes):
            row += overlap_genes[i] + "," + ov[i] + ","
        else:
            row += ",,"

        if i < len(antisense_genes):
            row += antisense_genes[i] + "," + an[i] + ","
        else:
            row += ",,"
        if i < len(intergenic_genes):
            row += intergenic_genes[i]
        w.write(row + "\n")

#Further classify overlap into outer_change or alternative splicing
exon_data = {}
with open("5_classification/0_transcript_info_from_gtf.csv", "r") as r:
    r = r.readlines()
    for x in r[1:]:
        x = x.rstrip().split(",")
        if x[1] == "exon":
            if x[6] not in exon_data.keys():
                exon_data[x[6]] = []
            exon_data[x[6]].extend([int(x[2]),int(x[3])])
with open ("5_classification/1_transcript_classification.csv", "r") as r:
    with open ("5_classification/2_overlap_classification.csv", "w") as w:
        w.write("Novel_id,Annotated_id,Novel_ss,Annotated_ss,Share_ss,5UTR_Diff,3UTR_Diff\n")
        r = r.readlines()
        for x in r[1:]:
            UTR_5 = "false"
            UTR_3 = "false"
            count = 0
            x = x.rstrip().split(",")
            novel = exon_data[x[0]]
            anno = exon_data[x[1]]
            total_novel = str(len(novel))
            total_anno = str(len(anno))
            if novel[0] != anno[0]:
                UTR_5 = "true"
            if novel[-1] != anno[-1]:
                UTR_3 = "true"
            for i in novel:
                if i in anno:
                    count += 1
            w.write(x[0] + "," + x[1] + "," + total_novel + "," + total_anno + "," + str(count) + "," + UTR_5 + "," + UTR_3 + "\n")
