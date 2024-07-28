import os
import re
import csv
import math
import glob
import numpy as np
import pandas as pd

def NCBI_length(file_list):
    for i in range(len(file_list)):
        file_name = file_list[i]
        # Reading NCBI Reference Sequence for all genes with gene-length
        ncbi_seq = pd.read_csv(file_name)
        # print(ncbi_seq.head())
        ncbi_seq = ncbi_seq.fillna(0)
        strain = ncbi_seq['strain name'].str.split('|', expand=True)
        for i in range(strain.shape[0]):
            if strain.iloc[i, 1] == 'ORF1ab':
                g_name = strain.iloc[i, 3][0:3]
                # gname = g_name
                strain.iloc[i, 3] = strain.iloc[i, 3].replace('nsp', '_')
                strain.iloc[i, 3] = strain.iloc[i, 3].split('_')[1]
                temp = re.findall(r'\d+', strain.iloc[i, 3])
                res = temp[0]
                # strain.iloc[i, 3] = strain.iloc[i, 3].extract('(\d+)')
                g_name = g_name + res
                strain.iloc[i, 1] = g_name
        # strain.rename(columns={1: 'gene_name'}, inplace=True)
        strain = strain.drop([0, 2, 3], axis=1)
        ncbi_seq = ncbi_seq.drop(['sl', 'strain name'], axis=1)
        ncbi_seq.insert(loc=0, column='gene_name', value=strain)
        # ncbi_seq = ncbi_seq.reset_index(drop=True)
        length = np.count_nonzero(ncbi_seq, axis=1) - 1
        ncbi_seq.insert(1, 'gene_length', length)
        ncbi_seq = ncbi_seq.set_index("gene_name")
        ncbi_seq.to_csv(file_name)
        # os.remove('Covid19_NCBI_Seq.csv')
        # ncbi_seq.loc['nsp12','gene_length'] = ncbi_seq.loc['nsp12','gene_length'] - 1 ####------------ TO BE REMOVED
        # print(ncbi_seq)


def create_csv(file_list, dir_ncbi_path): #  Only for files with ATGC
    # Organizing a directory for raw CSV files
    os.chdir(dir_ncbi_path)
    path_csv = '_raw_csv'
    if not os.path.exists(path_csv):
        # Create a new directory because it does not exist
        os.makedirs(path_csv)
    path_raw_csv = os.path.join(dir_ncbi_path,path_csv)
    os.chdir(path_raw_csv)

    for i in range(len(file_list)):
        file_name = file_list[i]
        print(file_name)

        # Create CSV from Fasta/Text files
        ln = 0
        CCount = 0
        nColumns = 0  # Counts the total Number of Coulumns to write the column numbers
        line_f = file_name
        line_f = line_f.strip()
        fil_ = os.path.basename(line_f)
        file = os.path.splitext(fil_)[0]
        extn = os.path.splitext(fil_)[1]
        sf = open(line_f, 'r')
        df = open('transit.txt', 'w')
        counter = 0
        for line in sf:
            print('\r Ececuting ' + str(counter) + ' row.......', end='')
            if line.startswith('>'):
                ln += 1 # Used for sl
                if (CCount > nColumns):
                    nColumns = CCount
                CCount = 0
                df.write("\n")
                df.write(str(ln) + ',')
                s = line.replace(',', '-')
                # t=s.replace(',','-')
                df.write(s[1:].rstrip('\n'))
                continue
            else:
                for li in line:
                    if (li == "\n"):
                        df.write(li.rstrip('\n'))

                    elif li not in ['A', 'T', 'G', 'C']:
                        lee = '-'
                        df.write(',' + lee.rstrip('\n'))
                        CCount += 1

                    else:
                        CCount += 1
                        lee = li.upper()
                        df.write(',' + lee.rstrip('\n'))

                    # print("Number of columns",CCount)

                    if (CCount > nColumns):
                        nColumns = CCount
            # print("Col", nColumns)
            # print("     .")
            counter += 1
        # print("Sequence alignment done  !!!.....")
        df.close()
        df = open('transit.txt', 'r')
        nf = open('transit1.txt', 'w')  # ####
        # print("Printing Column numbers ****")
        nf.write("sl,strain name")
        for j in range(1, nColumns + 1):
            nf.write(',' + str(j).rstrip('\n'))
            # print("     *")
        for line in df:
            nf.write(line)
        nf.close()
        df.close()
        sf.close()

        # To csv(using pandas module)
        df_to_csv = pd.read_csv('transit1.txt', low_memory= False)
        df_to_csv = df_to_csv.replace('NUN', None)
        df_to_csv.dropna(how="all", axis=1, inplace=True)
        #for col in range(1, df_to_csv.shape[1]):
            #df_to_csv.dropna(how='any', axis=0, inplace=True)
        csv_file = file + '_raw.csv'
        df_to_csv.to_csv(csv_file, index=False)
        print(f"\n{csv_file} File Created\n in {os.getcwd()}")
        #print(os.getcwd())
        os.remove('transit.txt')
        os.remove('transit1.txt')




# Finding the Consensus Sequence
def consensus(data_gene, file_name):
    # Reading From Input File
    data = data_gene.iloc[:, 1:]  ##  or data=df.loc[:,'1':]
    l = ['Total']
    lA = ['A']
    lT = ['T']
    lG = ['G']
    lC = ['C']

    col = data.shape[1]
    #print(col)
    for i in range(1, col+1):
        l.append(i)
        lA.append((data.iloc[:, i - 1]).str.count('A').sum())
        lT.append((data.iloc[:, i - 1]).str.count('T').sum())
        lG.append((data.iloc[:, i - 1]).str.count('G').sum())
        lC.append((data.iloc[:, i - 1]).str.count('C').sum())

    # Writting to csv WITH  csv Module
    atgc_file = 'Consensus_' + file_name + '.csv'
    with open(atgc_file, "w") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(l)
        wr.writerow(lA)
        wr.writerow(lT)
        wr.writerow(lG)
        wr.writerow(lC)
    #print(f"\nTotal number of A, T, G, C added to File:{atgc_file}")

    # Reference Sequence
    lstR = ['Consensus Sequence']
    df = pd.read_csv(atgc_file, index_col=0)
    df = df.reset_index()

    for i in range(1, df.shape[1]):
        # print(df.iloc[:,i].argmax())
        series = pd.Series(df.iloc[:, i])
        s = series.argmax()  # .argmax returns the row number(indices) of the largest value in a Series
        # t=s.nlargest(1)
        # print(s)
        t = df.iloc[s, 0]
        # print(t)
        lstR.append(t)
    # print(lstR)
    with open(atgc_file, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(lstR)
    #print("The Consensus Sequence Generated As:\n")

    df_consensus = pd.read_csv(atgc_file, index_col=0)
    #print("COnsensus Sequence File- \n", df_consensus)
    os.remove(atgc_file)
    return df_consensus                                        #-----------------End of Consensus Functions

# Function for calculationg nonZeroMutations
def mutation(df_consensus):
    #Calculating Separaed mutations in repeated Rows
    df_ref = df_consensus
    #print(df_ref)
    ref = df_ref.iloc[-1, 0:]  # extract the consensus sequence
    df_mutation = pd.DataFrame()
    df_mutation_row = 0
    for i in range(df_ref.shape[1]):
        # Mutation Calculation
        col = df_ref.loc['A':'C', str(i + 1)]
        muts = col[col != '0'].index
        if len(muts) != 0:
            for j in range(len(muts)):
                if (j != len(muts) and muts[j] != df_ref.loc['Consensus Sequence', str(i + 1)]):
                    df_mutation.loc[df_mutation_row, 'Position'] = i + 1
                    df_mutation.loc[df_mutation_row, 'Source'] = df_ref.loc['Consensus Sequence', str(i + 1)]
                    df_mutation.loc[df_mutation_row, 'Target'] = muts[j]
                    df_mutation_row += 1
        else:
            df_mutation.loc[df_mutation_row, 'Position'] = i + 1
            df_mutation.loc[df_mutation_row, 'Source'] = 0
            df_mutation.loc[df_mutation_row, 'Target'] = 0
            df_mutation_row += 1
    df_mutation = df_mutation.drop_duplicates()
    df_mutation = df_mutation.reset_index(drop=True)
    return(df_mutation)  #---------------------- Mutation Separated by adding rows with same informations
    #print(df_mutation)


def titv(df_mutation):
    data = df_mutation
    position = data.shape[1]  # ------------------------------getting number of columns-------------------

    # print(data.loc[])
    # print(strain.head())
    # ------------------------------------------------ Inserting new columns ------------------------------------------

    data.insert(position, "Transversion", 0)
    data.insert(position, "Transition", 0)
    data.insert(position, "T->G", 0)
    data.insert(position, "G->T", 0)
    data.insert(position, "C->A", 0)
    data.insert(position, "A->C", 0)
    data.insert(position, "C->G", 0)
    data.insert(position, "G->C", 0)
    data.insert(position, "T->A", 0)
    data.insert(position, "A->T", 0)
    data.insert(position, "T->C", 0)
    data.insert(position, "C->T", 0)
    data.insert(position, "G->A", 0)
    data.insert(position, "A->G", 0)

    # print(data.to_string())
    # ----------------------------------------------- Computing ti/tv -----------------------------------------------
    #lst_ref = []
    #lst_base = []

    for i in range(0, data.shape[0]):
        #print(i)
        ref = data.at[i, 'Source']
        base = data.at[i, 'Target']
        ag = ga = ct = tc = at = ta = gc = cg = ac = ca = gt = tg = ts = tv = syn = non_syn = 0
        #for j in range(0, len(lst_ref)):
        if (ref == 'A'):
            if (base == 'T'):
                at = at + 1
                data.at[i, 'A->T'] = at
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'C'):
                ac = ac + 1
                data.at[i, 'A->C'] = ac
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'G'):
                ag = ag + 1
                data.at[i, 'A->G'] = ag
                ts = ts + 1
                data.at[i, 'Transition'] = ts
        elif (ref == 'G'):
            if (base == 'T'):
                gt = gt + 1
                data.at[i, 'G->T'] = gt
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'C'):
                gc = gc + 1
                data.at[i, 'G->C'] = gc
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'A'):
                ga = ga + 1
                data.at[i, 'G->A'] = ga
                ts = ts + 1
                data.at[i, 'Transition'] = ts
        elif (ref == 'C'):
            if (base == 'A'):
                ca = ca + 1
                data.at[i, 'C->A'] = ca
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'G'):
                cg = cg + 1
                data.at[i, 'C->G'] = cg
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'T'):
                ct = ct + 1
                data.at[i, 'C->T'] = ct
                ts = ts + 1
                data.at[i, 'Transition'] = ts
        elif (ref == 'T'):
            if (base == 'A'):
                ta = ta + 1
                data.at[i, 'T->A'] = ta
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'G'):
                tg = tg + 1
                data.at[i, 'T->G'] = tg
                tv = tv + 1
                data.at[i, 'Transversion'] = tv
            elif (base == 'C'):
                tc = tc + 1
                data.at[i, 'T->C'] = tc
                ts = ts + 1
                data.at[i, 'Transition'] = ts
        else:
            # print("Its an ambiguous character--",lst_base[j])
            pass
    #print("------------------------------------ti and tv mutations calculated-------------------------")
    return data

def syn_nonSyn(df_consensus, df_titv):
    # Reading reference sequence
    df_ref = df_consensus
    df_ref = df_ref.reset_index()
    # print(df_ref)
    ref_seq = df_ref.iloc[4, 0:]
    #print("Reference Sequence: ",ref_seq)

    # Reading titv data
    #data = df_titv

    # New DataFrame for Syn Non-Syn Mutation
    #syn_nsyn = pd.DataFrame()
    syn_nsyn = df_titv
    # print(syn_nsyn.shape)
    position = syn_nsyn.shape[1]
    # print(position)
    # print("Adding Columns for Mutation Calculation")
    syn_nsyn.insert(position, "Non-Synonymous", 0)
    syn_nsyn.insert(position, "Synonymous", 0)
    syn_nsyn.insert(position, "target_AA", " ")
    syn_nsyn.insert(position, "source_AA", " ")
    syn_nsyn.insert(position, "target_Codon", " ")
    syn_nsyn.insert(position, "source_Codon", " ")
    syn_nsyn.insert(position, "Codon_Number", 0)

    # ------------------------------------------ Dictionary for Amino Acids--------------------
    amino_acid = {
        0: {'TTT': 'phenylalanine', 'TTC': 'phenylalanine'},
        1: {'TTA': 'leucine', 'TTG': 'leucine', 'CTT': 'leucine', 'CTC': 'leucine', 'CTA': 'leucine', 'CTG': 'leucine'},
        2: {'ATT': 'isoleucine', 'ATC': 'isoleucine', 'ATA': 'isoleucine'},
        3: {'ATG': 'methionine'},
        4: {'TGG': 'tryptophan'},
        5: {'GTT': 'valine', 'GTC': 'valine', 'GTA': 'valine', 'GTG': 'valine'},
        6: {'TCT': 'serine', 'TCC': 'serine', 'TCA': 'serine', 'TCG': 'serine', 'AGT': 'serine', 'AGC': 'serine'},
        7: {'CCT': 'proline', 'CCC': 'proline', 'CCA': 'proline', 'CCG': 'proline'},
        8: {'ACT': 'threonine', 'ACC': 'threonine', 'ACA': 'threonine', 'ACG': 'threonine'},
        9: {'GCT': 'alanine', 'GCC': 'alanine', 'GCA': 'alanine', 'GCG': 'alanine'},
        10: {'TAT': 'tyrosine', 'TAC': 'tyrosine'},
        11: {'CAT': 'histidine', 'CAC': 'histidine'},
        12: {'CAA': 'glutamine', 'CAG': 'glutamine'},
        13: {'AAT': 'asparagine', 'AAC': 'asparagine'},
        14: {'AAA': 'lysine', 'AAG': 'lysine'},
        15: {'GAT': 'aspartic_acid', 'GAC': 'aspartic_acid'},
        16: {'GAA': 'glutamic_acid', 'GAG': 'glutamic_acid'},
        17: {'TGT': 'cysteine', 'TGC': 'cysteine'},
        18: {'CGT': 'arginine', 'CGC': 'arginine', 'CGA': 'arginine', 'CGG': 'arginine', 'AGA': 'arginine',
             'AGG': 'arginine'},
        19: {'GGT': 'glycine', 'GGC': 'glycine', 'GGA': 'glycine', 'GGG': 'glycine'},
        20: {'TAA': 'stop_codon', 'TAG': 'stop_codon', 'TGA': 'stop_codon'}

    }

    # print("Dictionary for Amino Acids, Created.\n")
    # print(amino_acid.keys())
    # print("Adding Reference Base and Mutated Base from ti_tv file")
    for i in range(0, syn_nsyn.shape[0]):
        pos = syn_nsyn.loc[i, 'Position']
        codon_position = math.ceil(int(pos) / 3)
        #syn_nsyn.loc[i, 'Base_Position'] = (pos)
        syn_nsyn.loc[i, 'Codon_Number'] = codon_position
        #syn_nsyn.loc[i, 'source'] = data.loc[i, 'Source']
        #syn_nsyn.loc[i, 'target'] = data.loc[i, 'Target']

    # -------------------------------------------------Source Codon & Target Codon----------------------------------
    # print("Calculating Source and Target Codons")
    length = syn_nsyn.shape[0]

    for i in range(0, length):
        #print("Row : ", i)
        codon_no = syn_nsyn.at[i, 'Codon_Number']
        if codon_no >= 1:
            pos = ((syn_nsyn.at[i, 'Codon_Number'] - 1) * 3) + 1
            pos = int(pos)
            #print(pos)
            s_codon = ref_seq[pos] + ref_seq[pos + 1] + ref_seq[pos + 2]
            # print(s_codon)
            syn_nsyn.at[i, 'source_Codon'] = s_codon
            t_codon = ""
            for c_no in range(pos, pos + 3):
                if c_no == syn_nsyn.at[i, 'Position']:
                    t_codon = t_codon + syn_nsyn.at[i, 'Target']
                else:
                    t_codon = t_codon + ref_seq[c_no]

            syn_nsyn.at[i, 'target_Codon'] = t_codon

        else:
            s_codon = ""
            syn_nsyn.at[i, 'source_Codon'] = s_codon
            t_codon = ""
            syn_nsyn.at[i, 'target_Codon'] = t_codon
    # --------------------------------------Finding Changed amino acids -----------------------
    # print("Calculating Syn and NonSyn Mutations")

    for i in range(0, syn_nsyn.shape[0]):
        r_aa = syn_nsyn.loc[i, 'source_Codon']  # r_aa is reference/source amonoAcid
        m_aa = syn_nsyn.loc[i, 'target_Codon']  # m_aa is mutated/target amonoAcid
        for j in range(0, len(amino_acid)):
            if r_aa in amino_acid[j].keys():
                raa = amino_acid[j][r_aa]
                # print(raa)
                syn_nsyn.loc[i, 'source_AA'] = raa
                break
        for j in range(0, len(amino_acid)):
            if m_aa in amino_acid[j].keys():
                maa = amino_acid[j][m_aa]
                # print(raa)
                syn_nsyn.loc[i, 'target_AA'] = maa
                break
        if (syn_nsyn.loc[i, 'source_Codon'] != ""):
            if (syn_nsyn.loc[i, 'source_AA'] == syn_nsyn.loc[i, 'target_AA']):
                syn_nsyn.loc[i, 'Synonymous'] = 1
            elif (syn_nsyn.loc[i, 'source_AA'] != syn_nsyn.loc[i, 'target_AA']):
                syn_nsyn.loc[i, 'Non-Synonymous'] = 1
    syn_nsyn = syn_nsyn.fillna(0)
    return(syn_nsyn)

def normalized_mutations(data_synNsyn, data_consensus):
    # combination Source & Target
    df_mutation = data_synNsyn
    df_mutation['mutation'] = df_mutation['Source'] + '->' + df_mutation['Target']

    # Duplicate Position Count
    duplicate_position = df_mutation['Position'].value_counts()  # Series
    # print( duplicate_position.index)
    duplicate_position = duplicate_position.sort_index()
    duplicate_positions = pd.DataFrame(duplicate_position)
    duplicate_positions = duplicate_positions.reset_index()
    duplicate_positions.columns = ['Duplicate_Position', 'Counts']
    duplicate_positions = duplicate_positions[duplicate_positions['Counts'] > 1]
    duplicate_positions = duplicate_positions.reset_index(drop=True)
    # print(duplicate_positions)

    # A,T,G,C counts
    data_consensus = data_consensus.set_index('Total')
    data_consensus_T = data_consensus.transpose()
    atgc_count = data_consensus_T['Consensus Sequence'].value_counts()  # Series
    # print(atgc_count_loop)

    # A, T, G, C Count Adjustment --
    # 1. Search in duplivate_position, count>1  ->  Fetch the position number
    # 2. Match the position number in df_mutation_loop  ->  Fetch Source character
    # 3. Change the character value in atgc_count to 1 less than Duplicate_Position,counts
    for i in range(duplicate_positions.shape[0]):
        d_pos = duplicate_positions.loc[i, 'Duplicate_Position']
        count = duplicate_positions.loc[i, 'Counts']
        df_mutation_index = df_mutation[df_mutation['Position'] == d_pos].index
        df_mutation_index = df_mutation_index[0]
        source_base = df_mutation.at[df_mutation_index, 'Source']
        atgc_count[source_base] = atgc_count[source_base] + (count - 1)
        # print(atgc_count)

    # Count of 12 different mutations separately
    mutation_count = df_mutation['mutation'].value_counts()
    mut_lst = ['A->C', 'A->T', 'A->G', 'C->T', 'C->G', 'C->A', 'T->G', 'T->A', 'T->C', 'G->A', 'G->C', 'G->T']
    for i in mut_lst:
        if i not in mutation_count.index:
            mutation_count[i] = 0
    mutation_count = mutation_count.sort_index()
    mutation_counts = pd.DataFrame(mutation_count)
    mutation_counts = mutation_counts.reset_index()
    mutation_counts.columns = ['Mutations', 'Counts']
    mutation_counts.insert(mutation_counts.shape[1], "Normalized Frequency", 0)
    mutation_counts.insert(2, "Reference Neucleotide", 0)
    mutation_counts.insert(3, "Total Count of Ref Neucleotide", 0)

    for i in range(mutation_counts.shape[0]):
        if mutation_counts.loc[i, 'Mutations'][0] == 'A':
            mutation_counts.loc[i, 'Reference Neucleotide'] = 'A'
            mutation_counts.loc[i, 'Total Count of Ref Neucleotide'] = atgc_count['A']
            mutation_counts.loc[i, 'Normalized Frequency'] = round(mutation_counts.loc[i, 'Counts'] / atgc_count['A'],3)

        elif mutation_counts.loc[i, 'Mutations'][0] == 'T':
            mutation_counts.loc[i, 'Reference Neucleotide'] = 'T'
            mutation_counts.loc[i, 'Total Count of Ref Neucleotide'] = atgc_count['T']
            mutation_counts.loc[i, 'Normalized Frequency'] = round(mutation_counts.loc[i, 'Counts'] / atgc_count['T'],3)

        elif mutation_counts.loc[i, 'Mutations'][0] == 'G':
            mutation_counts.loc[i, 'Reference Neucleotide'] = 'G'
            mutation_counts.loc[i, 'Total Count of Ref Neucleotide'] = atgc_count['G']
            mutation_counts.loc[i, 'Normalized Frequency'] = round(mutation_counts.loc[i, 'Counts'] / atgc_count['G'],3)

        elif mutation_counts.loc[i, 'Mutations'][0] == 'C':
            mutation_counts.loc[i, 'Reference Neucleotide'] = 'C'
            mutation_counts.loc[i, 'Total Count of Ref Neucleotide'] = atgc_count['C']
            mutation_counts.loc[i, 'Normalized Frequency'] = round(mutation_counts.loc[i, 'Counts'] / atgc_count['C'],3)

    df_normFrequency = mutation_counts
    return(df_normFrequency)

###----------------------------------------------------------------------------------------- Main Body
current_path = os.getcwd() #gets current working directory
path_supportingFile = os.path.join(current_path,'Supporting Files')
path_ncbi = os.path.join(path_supportingFile,'ncbi')
#print(f"files inside {dir} - \n",os.listdir())


while 1:
    # Menu for File Formating
    menu = "1. Convert NCBI_Reference_Sequence.txt to .csv\n" \
           "2. Create CSV From Fasta Files\n" \
           "3. Fix Sequence Length Comparing to NCBI\n" \
           "4. Syn & Nsyn mutations \n" \
           "5. Normalized Frequencies \n" \
           "6. dN/dS\n" \
           "7. Main Menu\n" \
           "8. Exit\n"\
           "Enter What You Want to Do (In Number):- "
    choice = int(input(menu))
    if choice == 1:
        os.chdir(path_ncbi)
        path_fasta_file = path_ncbi + '\*.fasta'
        path_fas_file = path_ncbi + '\*.fas'
        path_text_file = path_ncbi + '\*.txt'

        if glob.glob(path_fasta_file):
            ncbi_file = glob.glob(path_fasta_file)
        elif glob.glob(path_fas_file):
            ncbi_file = glob.glob(path_fas_file)
        elif glob.glob(path_text_file):
            ncbi_file = glob.glob(path_text_file)

        ncbi_file = list(ncbi_file)
        create_csv(ncbi_file,path_ncbi)
        path_csv = '_raw_csv'
        path_raw_csv = os.path.join(path_ncbi, path_csv)
        os.chdir(path_raw_csv)
        csv_files = glob.glob("*_raw.csv")
        NCBI_length(csv_files)
        os.chdir(current_path)

    elif choice == 7:
        continue

    elif choice == 8:
        exit()

    else:
        os.chdir(path_ncbi)
        ncbi_seq = pd.read_csv("Covid19_NCBI_Seq.csv")
        ncbi_seq = ncbi_seq.set_index('gene_name')
        os.chdir(current_path)
        print(os.listdir())  # listdir() gives the list of files and folders
        dir = input(f"Enter the folder name :- ")
        path = os.path.join(current_path, dir)  # joins two strings as a path
        os.chdir(path)  # move into another directory
        path_raw_csv = os.path.join(path, '_raw_csv')
        dir_fixedLength = "_fixedLength"
        path_fixedLength = os.path.join(path, dir_fixedLength)
        dir_synNsyn = "_synNsyn"
        path_synNsyn = os.path.join(path, dir_synNsyn)
        path_normMutations = os.path.join(path,'_normMutations')

        #files = os.listdir()

        if choice == 2:  # Create CSV files
            dir_fasta = '_raw_fasta'
            path_raw_fasta = os.path.join(path,dir_fasta)
            os.chdir(path_raw_fasta)

            path_fasta_file = path_raw_fasta + '\*.fasta'
            path_fas_file = path_raw_fasta + '\*.fas'
            path_text_file = path_raw_fasta + '\*.txt'

            if glob.glob(path_fasta_file):
                files = glob.glob(path_fasta_file)
            elif glob.glob(path_fas_file):
                files = glob.glob(path_fas_file)
            elif glob.glob(path_text_file):
                files = glob.glob(path_text_file)

            create_csv(files, path)
            os.chdir(path)

        elif choice == 3: # Fix Sequence Length
            print("Reading Files...")
            if not os.path.exists(dir_fixedLength):
                # Create a new directory because it does not exist
                os.makedirs(dir_fixedLength)

            os.chdir(path_raw_csv)
            #print(glob.glob(path_raw_csv))
            if glob.glob('*_raw.csv'):  # glob.glob - collects all files with extension .csv
                #os.chdir(path_raw_csv)
                raw_csv_list = glob.glob('*_raw.csv')

                # Making CSV files Fixed Length as compared to Reference Sequence Files
                counter = 0
                for file in raw_csv_list:
                    data = pd.read_csv(file, low_memory=False)
                    '''
                    # Creating DataFrame by Reading chunks, to avoid low memory error
                    mylist = []
                    data = pd.DataFrame()
                    for chunk in pd.read_csv(file, sep=',', chunksize=20000):
                        mylist = chunk
                        data = pd.concat(mylist, axis=0)
                    del mylist
                    '''

                    data = data.drop(['sl'], axis=1)  # str(col),
                    print(f"Input file {file} is of shape {data.shape}") # file name format - <type of variant>_<protein>_raw.csv
                    data = data.dropna(subset=['1'])
                    data = data.reset_index(drop=True)
                    data['length'] = data.count(axis=1) - 1
                    gene_folder = os.path.splitext(os.path.basename(file))[0]  # Alpha_nsp1_raw
                    fileName_fixedColLength = gene_folder.split('_')[0]+ gene_folder.split('_')[1] + '_fixedLength.csv'
                    gene_folder = gene_folder.split('_')[1]  # nsp1
                    data = data[data.loc[:, 'length'] == ncbi_seq.loc[gene_folder, 'gene_length']]
                    data = data.drop("length", axis=1)
                    data = data.dropna(how = 'all',axis=1)
                    data = data.reset_index(drop=True)
                    filter = pd.DataFrame()
                    for col in range(30, data.shape[1]):     # Filtering out unwanted charasters
                        print(f'\rExecuting Column' + str(col) + ' .......', end='')
                        index_ = data[~(data[str(col)].isin(['A', 'T', 'G', 'C']))].index
                        df_index = pd.DataFrame(data.iloc[index_])
                        #print(df_index)
                        filter = filter.append(df_index)
                    #print(filter)
                    #filter.to_csv('filter.csv')
                    data = data[~data.index.isin(filter.index)] # filtering out rows in data, not present in filter
                    #print(f"After Fixing Length and Deleting Rows Having Unwanted Characters - {data.shape}")

                    os.chdir(path_fixedLength)
                    data.to_csv(fileName_fixedColLength, index=False)
                    os.chdir(path_raw_csv)
                    print(f"\nFile {fileName_fixedColLength} Created of shape - {data.shape} " )
                    counter+=1
            os.chdir(path)

        elif choice == 4: # From fixedLength to _synNsyn
            os.chdir(path_supportingFile)
            # Reading Synonymous values
            df_synCodon_Values = pd.read_csv("syn_codon_values.csv")
            os.chdir(path)
            if not os.path.exists('_synNsyn'):
                os.makedirs('_synNsyn')
            #path_synNsyn = os.path.join(path, '_synNsyn')
            if os.path.exists(dir_fixedLength):
                #path_fixedLength = path + dir_fixedLength
                os.chdir(path_fixedLength)
                #path_fixedLength_csv = path + '\*_fixedLength.csv'
                if glob.glob('*_fixedLength.csv'):
                    fixedLength_csv_files = glob.glob('*_fixedLength.csv')
                    for file in fixedLength_csv_files:
                        #b = len(fixedLength_csv_files)-i  # ----------------------Adding time counter
                        #c = '...'
                        #print('\rExecution stops in ' + str(b) + c, end='')
                        file_name = os.path.splitext(os.path.basename(file))[0]  # + os.path.splitext(os.path.basename(file))[1]
                        file_name = file_name.split('_')[0] + '_' + file_name.split('_')[1]
                        df_gene = pd.read_csv(file)
                        if df_gene.shape[0] > 0:
                            print(f"{file} of shape {df_gene.shape}\n")
                            print(f"Creating Consensus Sequence \n")
                            df_consensus = consensus(df_gene, file_name)
                            print(f"Calculating Mutations Spectra \n")
                            df_mutation = mutation(df_consensus)
                            # print("Mutations are: \n", df_mutation)
                            if df_mutation.shape[0] > 0:
                                print(f"Calculating Transitions and Transversions \n")
                                df_titv = titv(df_mutation)
                                # print(df_titv)
                                # df_titv.to_csv("titv1.csv", index = False)
                                print(f"Calculating Synonymous and non-Synonymous Mutations \n")
                                df_syn_nonSyn = syn_nonSyn(df_consensus, df_titv)
                                df_syn_nonSyn = df_syn_nonSyn.drop(
                                    ["T->G", "G->T", "C->A", "A->C", "C->G", "G->C", "T->A", "A->T", "T->C",
                                     "C->T", "G->A", "A->G"], axis=1)
                                # print(df_syn_nonSyn)

                                consensus_file = file_name + '_consensus.csv'
                                mutation_file = file_name + '_mutation.csv'
                                synNsyn_file = file_name + '_synNsyn.csv'
                                os.chdir(path_synNsyn)

                                df_consensus.to_csv(consensus_file)
                                df_mutation.to_csv(mutation_file, index=False)
                                df_syn_nonSyn.to_csv(synNsyn_file, index=False)
                            else:
                                consensus_file = file_name + '_consensus.csv'
                                mutation_file = file_name + '_mutation.csv'
                                os.chdir(path_synNsyn)
                                df_consensus.to_csv(consensus_file)
                                df_mutation.to_csv(mutation_file, index=False)

                            os.chdir(path_fixedLength)
                        else:
                            print(f"{file} Has No Data to Execute")

        elif choice == 5: #Normalized Frequencies
            os.chdir(path)
            if not os.path.exists('_normMutations'):
                # Create a new directory because it does not exist
                os.makedirs('_normMutations')
            if os.path.exists('_synNsyn'):
                os.chdir(path_synNsyn)
                synNsyn_files = glob.glob('*_synNsyn.csv')
                for file in synNsyn_files:
                    file_name = os.path.splitext(os.path.basename(file))[0]  # + os.path.splitext(os.path.basename(file))[1]
                    file_name = file_name.split('_')[0] + '_' + file_name.split('_')[1]

                    consensus_file = file_name + '_consensus.csv'
                    df_consensus = pd.read_csv(consensus_file)
                    df_syn_nonSyn = pd.read_csv(file)
                    print(f"Normalizing Mutations Spectra for {file} \n")
                    df_normalizedMutations = normalized_mutations(df_syn_nonSyn, df_consensus)
                    # print(df_normalizedMutations)
                    os.chdir(path_normMutations)
                    norm_mutation_file = file_name + '_normalized.csv'
                    df_normalizedMutations.to_csv(norm_mutation_file, index=False)
                    os.chdir(path_synNsyn)


        elif choice == 6: # Calculate dN/dS
            os.chdir(path_supportingFile)
            # Reading Synonymous values
            df_synCodon_Values = pd.read_csv("syn_codon_values.csv")
            df_dnds = pd.DataFrame(columns=['gene', 'NSyn_mutations', 'NSyn_sites', 'dN', 'Syn_mutations', 'Syn_sites', 'dS', 'dN/dS'])
            df_dnds_row = df_dnds.shape[0]
            os.chdir(path)
            if os.path.exists(dir_synNsyn):
                os.chdir(path_synNsyn)
                #path_fixedLength_csv = path + '\*_fixedLength.csv'
                if glob.glob('*_synNsyn.csv'):
                    synNsyn_files = glob.glob("*_synNsyn.csv")
                    for file in synNsyn_files:
                        file_name = os.path.splitext(os.path.basename(file))[0]  # + os.path.splitext(os.path.basename(file))[1]
                        file_name = file_name.split('_')[0] + '_' + file_name.split('_')[1]
                        df_dnds.loc[df_dnds_row, 'gene'] = file_name
                        df_syn_nonSyn = pd.read_csv(file)
                        print(f"{file} of shape {df_syn_nonSyn.shape}\n")
                        df_syn = df_syn_nonSyn.drop(df_syn_nonSyn[df_syn_nonSyn['Synonymous'] == 0].index)
                        df_syn = df_syn.reset_index(drop=True)
                        # codon_syn_site = df_syn_nonSyn.loc[:,'source_Codon']
                        # codon_syn_site['ds'] =
                        df_syn.insert(df_syn.shape[1], "dS", 0)
                        # print(df_syn.columns)
                        if df_syn.shape[0] == 0:
                            d_S = 0
                        else:
                            for i in range(df_syn.shape[0]):
                                print(f'\rExecuting Column' + str(i) + ' .......', end='')
                                codon = df_syn.loc[i, 'source_Codon']
                                index = df_synCodon_Values[df_synCodon_Values['Codon'] == codon].index
                                df_syn.loc[i, 'dS'] = df_synCodon_Values.loc[index[0], 'Syn_value']
                            d_S = df_syn['Synonymous'].sum() / df_syn['dS'].sum()
                        df_dnds.loc[df_dnds_row, 'Syn_mutations'] = df_syn['Synonymous'].sum()
                        df_dnds.loc[df_dnds_row, 'Syn_sites'] = df_syn['dS'].sum()
                        df_dnds.loc[df_dnds_row, 'dS'] = d_S

                        df_nsyn = df_syn_nonSyn.drop(df_syn_nonSyn[df_syn_nonSyn['Non-Synonymous'] == 0].index)
                        df_nsyn = df_nsyn.reset_index(drop=True)
                        df_nsyn.insert(df_nsyn.shape[1], "dN", 0)
                        if df_nsyn.shape[0] == 0:
                            d_N = 0
                        else:
                            for i in range(df_nsyn.shape[0]):
                                print(f'\rEcecuting Column' + str(i) + ' .......', end='')
                                codon = df_nsyn.loc[i, 'source_Codon']
                                index = df_synCodon_Values[df_synCodon_Values['Codon'] == codon].index
                                df_nsyn.loc[i, 'dN'] = 3 - df_synCodon_Values.loc[index[0], 'Syn_value']
                            d_N = df_nsyn['Non-Synonymous'].sum() / df_nsyn['dN'].sum()
                        df_dnds.loc[df_dnds_row, 'NSyn_mutations'] = df_nsyn['Non-Synonymous'].sum()
                        df_dnds.loc[df_dnds_row, 'NSyn_sites'] = df_nsyn['dN'].sum()
                        df_dnds.loc[df_dnds_row, 'dN'] = d_N
                        if d_S == 0:
                            dn_ds = 'NA'
                        else:
                            dn_ds = round(d_N / d_S, 5)
                        df_dnds.loc[df_dnds_row, 'dN/dS'] = dn_ds
                        df_dnds_row += 1
                        print(f"\ndN/dS Calculated for  {file_name}")
                        # Combining Syn and Nsyn dataframes
                        # df_synNsyn = df_syn.append(df_nsyn)

                    os.chdir(path)
                    dnds_file = dir + '_dnds.csv'
                    df_dnds.to_csv(dnds_file, index=False)


        else:
            print("Please Enter a Correct Option From  the MENU.")


























