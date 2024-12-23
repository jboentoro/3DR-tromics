1. Instruction
- This pipeline was designed for Plasmodium Zymoseq stranded RNA-seq by ZB lab (SBS, NTU).
- The library prep should be straned using zymoseq kit and the sequencer come from Novogene company. Any changes may lead to error when running script. 
- Any quetion or bugs, please contact me at tienquanghuy.duong@ntu.edu.sg

2. Manual

2.1 Preparation before running scripit
- Install Python3,Trim galore,HISAT2,Samtools, and Stringtie.
- Copy the whole folder (stranded_script containing: reference, pipeline.py, README.md, UDI.csv, prepDE.py3) to the working directory. Do not change any name or location of any file/folder.
- Move the raw data file downloaded from Novogene to this package_script folder. This script only handle one file at a time. 
   + The raw file always come in .zip, .tar, or .tar.gz file. This is ok to the script.
   + If the raw file not come in the above format, please unzip the file in advance and place the unzipped folder under this package_script folder.

- Prepare the "adapter.csv" file with exactly the same format as below. Header and the file name are important, please keep it unchanged.
Sample,i7_index,i5_index
Sample1,TGATTATACG,GTCGATTACA
Sample2,CAGCCGCGTA,ACTAGCCGTG
.............................

You can also provide the "adapter_UDI.csv" file instead of "adapter.csv". The format is exactly the same  as below. Header and the file name are important, please keep it unchanged.
Sample,UDI_index
Sample1,UDI_01
Sample2,UDI_34
.............................


2.2 Running script

- Command: python3 script.py
- This script will ask for raw zipped data. Please copy and enter the correct name of the raw file. If the file you enter is not in .zip, .tar, or .tar.gz, the script will ask the name of unzipped folder after a warning.
- Also, you can choose to get the raw read count beside the tpm and fpkm by input "Yes" when being asked
- Wait until the script finish its jobs.
- The final output will be saved in the same directory with this package_script folder.
