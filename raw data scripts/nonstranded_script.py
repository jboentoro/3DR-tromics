#This pipeline was designed for Plasmodium falciparum 3D7 smartseq unstranded RNA-seq by ZB lab (SBS, NTU)

#Necessary library
import os
import sys
import subprocess

#Input necessary files
print ("Start the pipeline")
core = int(input("Enter the number of core you want to use:"))
plasmodium = input("Enter the specie name (Option include: pf_3d7, pk_a1h1, pk_h):")
while plasmodium not in ["pf_3d7", "pk_a1h1", "pk_h"]:
    plasmodium = input("Enter the specie name. You have to input precisely the available option (pf_3d7, pk_a1h1, pk_h):")
final_name = input("Enter the name you want for the folder containing all the output:")

#For re-run some specific steps
re_run_check = input("Are you re-runing the script? Enter 'Yes' or 'No':")
if re_run_check == "Yes":
    re_run = input("If you are re-runing the script, from which step you want to re-run? Enter 'Begining' or 'Trim' or 'Map' or 'Quantify':")
    os.system("mv ../"+final_name+"/* .")
    os.system("rm -r ../"+final_name)

Begining = True
Trim = True
Map = True
Quantify = True
if re_run == "Trim":
    Begining = False
    os.system("rm -r 1_trim")
    os.system("rm -r 2_map")
if re_run == "Map":
    Begining = False
    Trim = False
    os.system("rm -r 2_map")
if re_run == "Quantify":
    Begining = False
    Trim = False
    Map = False

#Unzip, merge, relocate, and rename raw file
if Begining == True:
    entry_start = input("Enter name of the raw file downloaded from sequencer or the unzipped folder:")
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

    if os.path.isfile("adapter_UDP.csv") is True:
        UDI = {}
        with open ("UDP.csv", "r") as r:
            r = r.readlines()
            for x in r[1:]:
                x = x.rstrip().split(",")
                UDI[x[0]] = x[1]+","+x[2]
        with open ("adapter_UDP.csv", "r") as r:
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
            adapter[x[0]] = ("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"+x[1]+"ATCTCGTATGCCGTCTTCTGCTTG","CTGTCTCTTATACACATCTGACGCTGCCGACGA"+x[2]+"GTGTAGATCTCGGTGGTCGCCGTATCATT")

    smartseq1 = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    smartseq2 = "CCCATGTACTCTGCGTTGATACCACTGCTT"

    with open ("1_trim/run.sh", "w") as w:
        for key, value in adapter.items():
            w.write("trim_galore -q 20 --phred33 -a "+value[0]+" -a2 "+value[1]+" --stringency 5 --trim-n -e 0.2 -o ./1_trim --length 35 --paired 0_raw_reads/"+key+"_1.fq.gz 0_raw_reads/"+key+"_2.fq.gz\n")
            for n in range(1,5):
                w.write("trim_galore -q 20 --phred33 -a "+smartseq1+" -a2 "+smartseq2+" --stringency 5 --trim-n -e 0.2 -o ./1_trim --length 35 --paired 1_trim/"+key+"_1_val"*n+"_1.fq.gz 1_trim/"+key+"_2_val"*n+"_2.fq.gz\n")
                w.write("rm 1_trim/"+key+"_1_val"*n+"_1.fq.gz\n")
                w.write("rm 1_trim/"+key+"_2_val"*n+"_2.fq.gz\n")
            w.write("mv 1_trim/"+key+"*1.fq.gz 1_trim/"+key+"_trimmed_1.fq.gz\n")
            w.write("mv 1_trim/"+key+"*2.fq.gz 1_trim/"+key+"_trimmed_2.fq.gz\n")

    command_per_file = (len(adapter)//core + 1)*15
    script_count = 0
    scripts = []

    with open ("1_trim/run.sh", "r") as r:
        r = r.readlines()
        while script_count <= core:
            script_file = "1_trim/run_"+str(script_count)+".sh"
            with open (script_file, "w") as w:
                for x in r[(script_count*command_per_file):((script_count + 1)*command_per_file)]:
                    w.write(x)
            scripts.append(script_file)
            script_count += 1
    with open ("script_output.log","w") as log_file:
        processes = [subprocess.Popen(["bash",script], stdout = log_file, stderr = log_file) for script in scripts]
    for process in processes:
        process.wait()

    print ("All running scripts have finished")
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
        w.write("Sample\tTrim1\tTrim2\tTrim3\tTrim4\tTrim5\n")
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
        w.write("Sample\tTrim1\tTrim2\tTrim3\tTrim4\tTrim5\n")
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
                w.write("hisat2 -S ./2_map/"+sample+".sam --summary-file ./2_map/"+sample+".summ --min-intronlen 30 --max-intronlen 3000 -x reference/"+plasmodium+"/"+plasmodium+" -1 ./1_trim/"+entry+" -2 ./1_trim/"+entry.replace("1.fq.gz","2.fq.gz")+" --fr -p "+core+"\n")
                w.write("samtools view -@ "+core+" -b -o ./2_map/"+sample+".bam ./2_map/"+sample+".sam\n")
                w.write("rm 2_map/"+sample+".sam\n")
                w.write("samtools view -@ "+core+" -H ./2_map/"+sample+".bam > ./2_map/"+sample+".unq.sam\n")
                w.write("samtools view -@ "+core+" -f 2 ./2_map/"+sample+".bam | fgrep -w NH:i:1 >> ./2_map/"+sample+".unq.sam\n")
                w.write("samtools sort -n -@ "+core+" ./2_map/"+sample+".unq.sam -o ./2_map/"+sample+".unq.bam\n")
                w.write("rm 2_map/"+sample+".unq.sam\n")
                w.write("bamToBed -bedpe -i ./2_map/"+sample+".unq.bam > ./2_map/"+sample+".unq.bed\n")
                w.write("rm 2_map/"+sample+".unq.bam\n")
                w.write("coverageBed -a reference/"+plasmodium+"/"+plasmodium+"_RNA.gff -b ./2_map/"+sample+".unq.bed > ./2_map/"+sample+".unq.cov\n")
                w.write("rm ./2_map/"+sample+".unq.bed\n\n")
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

# Calculate fpkm value
if Quantify == True:
    print ("Start calculating fpkm for each sample")
    os.system("perl calculate_fpkm.pl")

    name_sample = "Sample"
    transcript = {}
    for entry in os.listdir("2_map"):
        if ".fpkm.txt" in entry:
            name_sample += ("," + str(entry[:-9]))
            with open ("2_map/"+entry, "r") as r:
                r = r.readlines()
                for x in r[1:]:
                    transcript_id = x.rstrip().split("\t")[0]
                    transcript_fpkm = x.rstrip().split("\t")[1]
                    if transcript_id not in transcript.keys():
                        transcript[transcript_id] = ("," + str(transcript_fpkm))
                    else:
                        transcript[transcript_id] += ("," + str(transcript_fpkm))
    with open ("2_map/0_fpkm_summ.csv", "w") as w:
        w.write(name_sample + "\n")
        for key, value in transcript.items():
            w.write(key + value + "\n")

# Summary read count value
    print ("Start summary the read count")
    count_list = {}
    sample = []
    x = []
    for entry in os.listdir("2_map"):
        if ".unq.cov" in entry: x.append(entry.split(".")[0])
    x.sort()

    for entry in x:
        sample.append(entry)
        with open ("2_map/"+ entry + ".unq.cov") as r:
            q = r.readlines()
            for x in q:
                y = x.split ("\t")
                info = y[8].split(';')[1][7:] + "," + y[8].split(';')[0][3:]
                if info not in count_list.keys():
                    count_list[info] = str(y[9])
                else:
                    count_list[info] += ("," + str(y[9]))

    with open ("2_map/0_count_summ.csv", "w") as w:
        w.write("gene_id,transcript_id")
        for x in sample:
            w.write("," + x)
        for key, value in count_list.items():
            w.write ("\n" + key + "," + value)

# Move all output to a folder outside the package

    initial_script = ["script.py", "calculate_fpkm.pl", "README.md", "reference", "UDP.csv", "script_test.py", final_name]
    os.system("mkdir "+final_name)
    for entry in os.listdir("."):
        if entry not in initial_script:
            os.system("mv "+entry+" "+final_name)
    os.system("mv "+final_name+" ..")
    print ("All done!")
