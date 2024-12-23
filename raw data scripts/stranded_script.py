#This pipeline was designed for stranded Zymo-seq total RNA-seq by ZB lab (SBS, NTU)

#Necessary library
import os
import sys
import subprocess

#Input necessary files
print ("Start the pipeline")
core =str(input("Enter the number of core you want to use:"))
plasmodium = input("Enter the specie name (Option include: pf_3d7, pk_a1h1, pk_h):")
while plasmodium not in ["pf_3d7", "pk_a1h1", "pk_h"]:
    plasmodium = input("Enter the specie name. You have to input precisely the available option (pf_3d7, pk_a1h1, pk_h):")
raw_count = input("Do you need the raw read count beside the tpm/fpkm? Enter 'Yes'or 'No':")
classification = input("Do you need to classify the novel transcript to integenic,antisense,overlap? Enter 'Yes' or 'No':")
final_name = input("Enter the name you want for the folder containing all the output:")

#For re-run some specific steps
re_run_check = input("Are you re-runing the script? Enter 'Yes' or 'No':")
if re_run_check == "Yes":
    re_run = input("If you are re-runing the script, from which step you want to re-run? Enter 'Begining' or 'Trim' or 'Map' or 'Assembly':")
    os.system("mv ../"+final_name+"/* .")
    os.system("rm -r ../"+final_name)
    os.system("rm -r 4_count")
    os.system("rm -r 5_classification")
else:
    re_run = "No"
Begining = True
Trim = True
Map = True
Assembly = True
if re_run == "Trim":
    Begining = False
    os.system("rm -r 1_trim")
    os.system("rm -r 2_map")
    os.system("rm -r 3_assembly")
if re_run == "Map":
    Begining = False
    Trim = False
    os.system("rm -r 2_map")
    os.system("rm -r 3_assembly")
if re_run == "Assembly":
    Begining = False
    Trim = False
    Map = False
    os.system("rm -r 3_assembly")

#Unzip, merge, relocate, and rename raw file
if Begining == True:
    entry_start = input("Enter name of the raw zipped file or the unzipped folder with all raw file:")
    total_sample = []
    total_sample.append(entry_start)

    other_entry = "ask for more"
    while other_entry != "No":
        other_entry = input("Still have samples? Enter it here. If not, then enter exactly 'No': ")
        if other_entry != "No":
            total_sample.append(other_entry)

    os.system("mkdir all_raw_input")

    for entry_start in total_sample:
        if os.path.isfile(entry_start) is True:
            if ".zip" in entry_start:
                os.system("unzip " + entry_start)
                raw_dir = entry_start.split(".")[0]
            elif ".tar.gz" in entry_start:
                os.system("tar -xvzf " + entry_start)
                raw_dir = entry_start.split(".")[0]
            elif ".tar" in entry_start:
                os.system("tar -xvf " + entry_start)
                raw_dir = entry_start.split(".")[0]
            else:
                sys.exit("The zip files: " + entry_start + " are not supported. Please unzip the raw file and run again")
        elif os.path.isdir(entry_start) is True:
            raw_dir = entry_start
        else:
            sys.exit("Make sure your input file/folder: "+entry_start+" is kept in the same directory with this script!")
        os.system("mv "+raw_dir+" all_raw_input")

    os.system("mkdir 0_raw_data")
    file_list = []
    for root, dirs, files in os.walk("all_raw_input"):
        for file in files:
            if file.endswith(".fq.gz"):
                file_list.append(os.path.join(root))
    file_list = list(set(file_list))
    for entry in file_list: os.system("mv " + entry + " 0_raw_data")

    os.system("rm -r all_raw_input")
    os.system("mkdir 0_raw_reads")
    for entry in os.listdir("0_raw_data"):
        os.system("cat 0_raw_data/" + entry + "/*1.fq.gz > " + entry + "_1.fq.gz")
        os.system("cat 0_raw_data/" + entry + "/*2.fq.gz > " + entry + "_2.fq.gz")
    os.system("mv *fq.gz 0_raw_reads")
    os.system("rm -r 0_raw_data")

#Adapter and quality trimming
if Trim == True:
    print ("Start adapter and quality trimming")
    if os.path.isfile("adapter_UDI.csv") is True:
        UDI = {}
        with open ("UDI.csv", "r") as r:
            r = r.readlines()
            for x in r[1:]:
                x = x.rstrip().split(",")
                UDI[x[0]] = x[1]+","+x[2]
        with open ("adapter_UDI.csv", "r") as r:
            with open ("adapter.csv", "w") as w:
                w.write("Sample,i7_index,i5_index\n")
                r = r.readlines()
                for x in r[1:]:
                    x = x.rstrip().split(",")
                    w.write(x[0]+","+UDI[x[1]]+"\n")

    if os.path.isfile("adapter.csv") is False:
        sys.exit("Trimming requires the 'adapter.csv' file in the current directory before this script")

    os.system("mkdir 1_trim")
    adapter = {}
    with open ("adapter.csv", "r") as r:
        r = r.readlines()
        for x in r[1:]:
            x = x.rstrip().split(",")
            adapter[x[0]] = ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"+x[1]+"ATCTCGTATGCCGTCTTCTGCTTG","AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"+x[2]+"GTGTAGATCTCGGTGGTCGCCGTATCATT")


    with open ("1_trim/run.sh", "w") as w:
        for key, value in adapter.items():
            w.write("trim_galore -q 30 --phred33 --clip_R2 10 -a NNNNNNNNNN"+value[0]+" -a2 "+value[1]+" --stringency 15 --trim-n -e 0.2 -o ./1_trim --length 75 --paired 0_raw_reads/"+key+"_1.fq.gz 0_raw_reads/"+key+"_2.fq.gz\n")
            w.write("trim_galore -q 30 --phred33 --trim-n -o ./1_trim --length 35 --paired 1_trim/"+key+"_1_val_1.fq.gz 1_trim/"+key+"_2_val_2.fq.gz\n")
            w.write("rm 1_trim/"+key+"_1_val_1.fq.gz\n")
            w.write("rm 1_trim/"+key+"_2_val_2.fq.gz\n")
            w.write("mv 1_trim/"+key+"*1.fq.gz 1_trim/"+key+"_trimmed_1.fq.gz\n")
            w.write("mv 1_trim/"+key+"*2.fq.gz 1_trim/"+key+"_trimmed_2.fq.gz\n")
    
    command_per_file = (len(adapter)//int(core) + 1)*6
    script_count = 0
    processes = []
    with open ("1_trim/run.sh", "r") as r:
        r = r.readlines()
        while script_count <= int(core):
            script_file = "1_trim/run_"+str(script_count)+".sh"
            with open (script_file, "w") as w:
                for x in r[(script_count*command_per_file):((script_count + 1)*command_per_file)]:
                    w.write(x)
            process = subprocess.Popen(["bash",script_file])
            processes.append(process)
            script_count += 1

    for process in processes:
        process.wait()
    os.system("rm 1_trim/run_*")

#Trimming quality control
    print ("Start trimming QC")

    os.system("grep 'Total reads processed:' 1_trim/* | sort > 1_trim/reads_processed")
    os.system("grep 'Reads with adapters:' 1_trim/* | sort > 1_trim/reads_with_adapter")

    sample_qc = {}
    with open ("1_trim/reads_processed", "r") as r:
        r = r.readlines()
        for x in r:
            x = x.rstrip()
            if "_val" not in x:
                sample_qc[x.split(".fq.gz")[0][7:]] = x.split("  ")[-1].replace(",", "")
            else:
                sample_qc[x.split("_val")[0][7:]] += ("\t" + x.split("  ")[-1]).replace(",", "")
    with open ("1_trim/0_read_processed.txt", "w") as w:
        w.write("Sample\tSpecified_adapter_trim\tIllumina_adapter_trim\n")
        for key, value in sample_qc.items():
            w.write(key + "\t" + value + "\n")

    sample_qc = {}
    with open ("1_trim/reads_with_adapter", "r") as r:
        r = r.readlines()
        for x in r:
            x = x.rstrip()
            if "_val" not in x:
                sample_qc[x.split(".fq.gz")[0][7:]] = x.split("  ")[-1].split(" (")[0].replace(",", "")
            else:
                sample_qc[x.split("_val")[0][7:]] += ("\t" + x.split("  ")[-1]).split(" (")[0].replace(",", "")
    with open ("1_trim/0_read_with_adapter.txt", "w") as w:
        w.write("Sample\tSpecified_adapter_trim\tIllumina_adapter_trim\n")
        for key, value in sample_qc.items():
            w.write(key + "\t" + value + "\n")

    os.system("rm 1_trim/reads_processed")
    os.system("rm 1_trim/reads_with_adapter")

# Mapping to reference genome
if Map == True:
    print ("Start mapping and estimating gene coverage")

    os.system("mkdir 2_map")
    with open ("2_map/run.sh", "w") as w:
        for entry in os.listdir("1_trim"):
            if "_trimmed_1.fq.gz" in entry:
                sample = entry[:-16]
                w.write ("hisat2 --dta -S 2_map/"+sample+".sam --summary-file 2_map/"+sample+".summ --min-intronlen 20 --max-intronlen 3000 -x reference/"+plasmodium+"/"+plasmodium+" --rna-strandness RF -1 ./1_trim/"+entry+" -2 ./1_trim/"+entry.replace("1.fq.gz","2.fq.gz")+" -p "+core+"\n")
                w.write ("samtools view -@ "+core+" -b -o 2_map/"+sample+".bam 2_map/"+sample+".sam\n")
                w.write ("rm 2_map/"+sample+".sam\n")
                w.write ("samtools view -@ "+core+" -H 2_map/"+sample+".bam > 2_map/"+sample+".unq.sam\n")
                w.write ("samtools view -@ "+core+" -F 4 2_map/"+sample+".bam | grep 'NH:i:[123]' | grep -v NH:i:10 >> 2_map/"+sample+".unq.sam\n")
                w.write("samtools sort -@ "+core+" 2_map/"+sample+".unq.sam -o 2_map/"+sample+".unq.bam\n")
                w.write("rm 2_map/"+sample+".unq.sam\n")
    os.system("bash 2_map/run.sh")

# Mapping Quality Control
    print ("Start mapping QC")

    with open ("2_map/0_mapping_summary.csv", "w") as w:
        w.write("Sample,Total_read,Aligned_concordant_1,Aligned_concordant>1,Aligned_discordent_1,Aligned_1,Aligned>1\n")
        for entry in os.listdir("2_map"):
            if ".summ" in entry:
                with open ("2_map/"+entry, "r") as r:
                    r = r.readlines()
                    z = entry[:-5] +","
                    count = 0
                    for x in r:
                        if count in [0,3,4,7,12,13]:
                            y = x.split()
                            z += (str(y[0]) +",")
                            count += 1
                        else:
                            count += 1
                    w.write(z[:-1] +"\n")

#Assembly and quantification
if Assembly == True:
    print ("Start performing assembly with reference")
    os.system("mkdir 3_assembly")

    with open ("3_assembly/0_assembly.sh", "w") as w:
        for entry in os.listdir("2_map"):
            if ".unq.bam" in str(entry):
                w.write("stringtie -p "+core+" --rf -G reference/"+plasmodium+"/"+plasmodium+".gff -o 3_assembly/" + entry[:-8] + "_assembly.gtf 2_map/" + entry + "\n")
    os.system("bash 3_assembly/0_assembly.sh")
    #Merge
    with open ("3_assembly/1_mergelist.txt", "w") as w:
        for entry in os.listdir("3_assembly"):
            if "assembly.gtf" in entry:
                w.write("3_assembly/" + entry + "\n")
    os.system("stringtie --merge -p "+core+" -m 150 -f 0.1 -G reference/"+plasmodium+"/"+plasmodium+".gff -o 3_assembly/1_stringtie_merged.gtf 3_assembly/1_mergelist.txt")
    os.system ("rm 3_assembly/1_mergelist.txt")

    os.system("mkdir 3_assembly/expression")
    with open ("3_assembly/2_express.sh", "w") as w:
        for entry in os.listdir("2_map"):
            if ".unq.bam" in entry:
                w.write ("mkdir 3_assembly/expression/" + entry[:-8] + "\n")
                w.write ("stringtie -e -B -p "+core+" -G 3_assembly/1_stringtie_merged.gtf -o 3_assembly/expression/" + entry[:-8] + "/" + entry[:-8] +"_express.gtf 2_map/" + entry + "\n")
    os.system("bash 3_assembly/2_express.sh")

# Assembly Quality Control
    print ("Start assembly QC")
    with open ("3_assembly/3_assembly.summ", "w") as w:
        w.write("Sample,No of annotated transcript, No of novel transcript\n")
        for entry in os.listdir("3_assembly"):
            if ".gtf" in entry:
                with open ("3_assembly/"+entry, "r") as r:
                    r = r.readlines()
                    anno = 0
                    novel = 0
                    for x in r:
                        if "\ttranscript\t" in x:
                            if 'reference_id' in x or 'ref_gene_id' in x:
                                anno += 1
                            else:
                                novel += 1
                    w.write(entry + "," + str(anno) + "," + str(novel) + "\n")

# Summary tpm and fpkm value
    print ("Start summary the tpm and fpkm")
    tpm_list = {}
    fpkm_list = {}
    sample = []
    x = []
    for entry in os.listdir("3_assembly/expression"):
        if os.path.isdir("3_assembly/expression/" +entry) is True: x.append(entry)
    x.sort()

    for entry in x:
        sample.append(entry)
        with open ("3_assembly/expression/"+entry +"/"+ entry + "_express.gtf") as r:
            q = r.readlines()
            for x in q[2:]:
                y = x.split ("\t")
                if y[2] == "transcript":
                    z = y[8].split('"')
                    info = z[1] + "," + z[3]
                    if 'TPM' not in z[-3]:
                        print ("There is error with this: " + entry + " " + info)
                        if info not in tpm_list.keys():
                            tpm_list[info] = str("NA")
                            fpkm_list[info] = str("NA")
                        else:
                            tpm_list[info] += ("," + str("NA"))
                            fpkm_list[info] += ("," + str("NA"))
                    else:
                        if info not in tpm_list.keys():
                            tpm_list[info] = str(z[-2])
                            fpkm_list[info] = str(z[-4])
                        else:
                            tpm_list[info] += ("," + str(z[-2]))
                            fpkm_list[info] += ("," + str(z[-4]))

    with open ("3_assembly/expression/tpm.csv", "w") as w:
        with open ("3_assembly/expression/fpkm.csv", "w") as e:
            w.write("gene_id,transcript_id")
            e.write("gene_id,transcript_id")
            for x in sample:
                w.write("," + x)
                e.write("," + x)
            for key, value in tpm_list.items():
                w.write ("\n" + key + "," + value)
            for key, value in fpkm_list.items():
                e.write ("\n" + key + "," + value)

#Estimate read count for all gene
    if raw_count == "Yes" or raw_count == "yes":
        print("Start calculating read count")
        os.system("mkdir 4_count")
        os.system("python3 prepDE.py3 -i 3_assembly/expression -g 4_count/gene_count.csv -t 4_count/transcript_count.csv")

#Novel transcript classification
    if classification == "Yes" or classification == "yes":
        print("Start novel transcript classification")
        os.system("mkdir 5_classification")
        os.system("python3 classification.py")

# Move all output to a folder outside the package
    print ("Start moving all output to the directory outside this package")
    initial_script = ["script.py", "README.md", "reference", "UDI.csv", "prepDE.py3", "classification.py", "script_test.py", final_name]
    os.system("mkdir "+final_name)
    for entry in os.listdir("."):
        if entry not in initial_script:
            os.system("mv "+entry+" "+final_name)
    os.system("mv "+final_name+" ..")
    print ("All done!")
