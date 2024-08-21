# Nucleotide-Polymorphism
The program identifies the consensus sequence and detects mutations for a given genome sequence. Depending on the mutations, Synonymous/non-synonymous mutations, ti/tv ratio, and dn/ds are found.

# What is a 'Nucleotide-Polymorphism' :
Nucleotide polymorphism is the change in a single base position in the genome, and it can be located in coding or non-coding regions. Otherwise known as Mutations. It can affect viral genomes in several ways, like the rate of multiplying (spreading), the severity of infection, etc.

# Installation of the Program

## Download Zip
   1. Go to the link https://github.com/MDash-NITAP/Nucleotide-Polymorphism or https://github.com/MDash-NITAP/Nucleotide-Polymorphism.git
   2. Select **Download Zip** option from the drop-down menu for Code<>      
   3. Extract the Zip folder in a desired location
   4. Prepare the data in the downloaded directories as below:
      
         ### Preparing directories for the program run:
         ### For NCBI file
         1. Inside the folder **' Supporting Files/ncbi** is a downloaded file **Covid19_NCBI_Seq.txt**
         2. Create the .csv file for **Covid19_NCBI_Seq.txt** by choosing the option 1, after running the program.
         ### For the Gene file
         1. Download two data folders from the Zenodo repository: https://doi.org/10.5281/zenodo.13355486
         2. Extract the folders **Aug2020-Sep2023** and **Dec2019-July2020**
         3. Move them into the downloaded repository on the local machine.
            
   5. Run the **_Polymorphism_CtoT.py_** file.

# Program Execution
The program can be run on a Command Prompt, Power Shell, or any GUI with an installed Python Interpreter.

# What the program offers:
It is a menu-driven program, which can be run on the Command Prompt/ PowerShell/ any GUI for Python Programming.
The outputs of the program are:
  1. It can create csv file from any text/fasta files (genome sequence/ Reference)
  2. Can align the genome sequence file by 'Fixed Length' Option
  3. Calculate Mutations and display their normalized value as described in the related journal paper
  4. Filters out Synonymous and non-synonymous mutations
  5. Calculates dN/dS (to know more about dN/dS: Please refer to the related journal paper)

# How it works:
         
Upon running the program, it will display a menu:

  1. Convert NCBI_Reference_Sequence.txt to .csv
  2. Create CSV From Fasta Files
  3. Fix Sequence Length Compared to NCBI
  4. Syn & Nsyn mutations 
  5. Normalized Frequencies 
  6. dN/dS
  7. Main Menu
  8. Exit
  Enter What You Want to Do (In Number):-

The user must enter the corresponding number to the required option and press Enter Key from the keyboard.
After processing data it will generate csv files and store them in the present working directory.

> **To execute any of the options from 4,5,6 in the menu, it's mandatory to execute options 1,2 and 3 sequentially.**

[ Specifically for the data provided in https://doi.org/10.5281/zenodo.13355486, options 2 and 3 are being executed, and the obtained .csv files are provided in the same repository, considering the volume of data. ]

## Input files and output files

### NCBI Input File
The genome sequence files have to be saved in the folder 'ncbi' before  executing the program. An example of a fasta/text file is given below:

      1 >266..805|ORF1ab|GU280_gp01|nsp1
      2 atggagagccttgtccctggtttca...
      3 >806..2719|ORF1ab|GU280_gp01|nsp2
      4 gcatacactcgctatgtcgataacaacttctgtggccctgatggctaccctctt...
      5 >2720..8554|ORF1ab|GU280_gp01|nsp3,former nsp1
      6 gcaccaacaaaggttacttttggtgatgacactgtgatagaag...
      .
      .
      .


### Genome sequence Input File
The input files have to be saved in the folder '_raw_fasta' before  executing the program. An example of a fasta/text file is given below:

    1 ># Query: nsp11||
    2 >
    3 >257437|13426|13464
    4 TCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTG
    5 >257436|13404|13442
    6 TCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTG

### Output Files
Upon execution, all the output will be available in the respective folders as csv files.





> For any questions or issues, please refer to the paper or contact madhusmita.dash81@gmail.com/ madhusmita.phd@nitap.ac.in.
