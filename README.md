
# Nucleotide-Polymorphism
The program identifies the consensus sequence and detects mutations for a given genome sequence. Depending on the mutations, Synonymous/non-synonymous mutations, ti/tv ratio, and dn/ds are found.

# What is a 'Nucleotide-Polymorphism' :
Nucleotide polymorphism is the change in a single base position in the genome, and it can be located in coding or non-coding regions. Otherwise known as Mutations. It can affect viral genomes in several ways, like the rate of multiplying (spreading), the severity of infection, etc.

# What the program offers:
It is a menu-driven program, which can be run on the Command Prompt/ PowerShell/ any GUI for Python Programming.
The outputs of the program are:
  1. It can create csv file from any text/fasta files (genome sequence/ Reference)
  2. Can align the genome sequence file by 'Fixed Length' Option
  3. Calculate Mutations and display their normalized value as described in the related journal paper
  4. Filters out Synonymous and non-synonymous mutations
  5. Calculates dN/dS (to know more about dN/dS: Please refer to the related journal paper)

# How it works:

## Preparing directories for the program run:
  ### For NCBI file
    1. Create a folder 'Supporting Files' in the current working directory (where the program file is saved).
    2. Create another folder 'ncbi' inside the folder 'Supporting Files'
  ### For the Gene file
    1. Create a folder with the name of the gene
    2. Inside the above folder, create another folder '_raw_fasta'
    3. Save the downloaded fasta/txt file(s) inside the folder '_raw_fasta'
         
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
Upon execution, all the output files will be available in the respective folders.
