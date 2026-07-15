from Bio import SeqIO
import re


def get_uniprot_sequence_dict(fasta_file_path, protein_list):

    """ Takes fasta file path, protein list and empty dictionary and returns a dictionary
        which contains the database sequences of all proteins in the list
    """
    seq_dict = {}
    string = "Reverse_"
    records = SeqIO.parse(fasta_file_path, 'fasta')
    for record in records:
        if not string in record.name:
            id = record.id.split('|')[1] if '|' in record.id else record.id
            if id in protein_list:
                seq_dict[id] = (record.seq, record.description)
    return seq_dict


def find_labeled_residue(row, seq_dict):

    """ Using the dictionary which contains the database sequences of the proteins of interest
        function returns the location of the labelled residue for the specific protein (uniprot) 
        and its sequence (peptide sequence)
    """
    uniprot = row[1]
    peptide_seq = row[3]

    if uniprot in seq_dict:
        database_sequence = seq_dict[uniprot]
        database_sequence = str(database_sequence)

        pepsplit = peptide_seq.split('.')[1]
        #get position of *
        pos_of_mod = pepsplit.find("*")
        #remove mod from peptide sequence
        pepsplit_clean = pepsplit.replace("*", "")
        #find pos of cleaned peptide in sequence
        pos_in_sequence = [i.start() for i in re.finditer(pepsplit_clean, database_sequence)]
        #make calculation pos 1 + pos 2 - 2
        pos_in_sequence[:] = [pos + pos_of_mod for pos in pos_in_sequence]
        #change list to string and and C at the beginning
        labelled_residues = (";").join(["C" + str(elem) for elem in pos_in_sequence])

    else:
        labelled_residues = "no_db_seq"

    return labelled_residues

def getAARange(row, seq_dict):
    sequence = row["SEQUENCE"]
    p_line_locus = row["uniprot_id"]

    uniprot_id = str(p_line_locus)
    peptide_seq = str(sequence)

    if uniprot_id in seq_dict:
        database_sequence = seq_dict[uniprot_id][0]
        database_sequence = str(database_sequence)
        #remove mod from sequence
        peptide_seq = re.sub("\((.*?)\)", "", peptide_seq)
        #split to get middle part
        pepsplit = peptide_seq.split('.')[1]
        #find pos of cleaned peptide in database sequence
        pos_in_sequence = [i.start() for i in re.finditer(pepsplit, database_sequence)]
        # #make calculation pos 1 + pos 2 - 2
        pos_in_sequence = pos_in_sequence[0]
        endpos_in_sequence = pos_in_sequence + len(pepsplit)
        #pos_in_sequence[:] = [pos + len(pepsplit) for pos in pos_in_sequence]
        # #change list to string and and C at the beginning
        aa_range = str(pos_in_sequence+1)+"-"+str(endpos_in_sequence)
    else:
        aa_range = "no_db_seq"
    return(aa_range)

def get_protein_description(row, seq_dict): 
    p_line_locus = row["uniprot_id"]
    uniprot_id = str(p_line_locus)

    if uniprot_id in seq_dict:
        descr = seq_dict[uniprot_id][1]
        descr = descr.split(" ", 1)[1]
        descr = descr.split(" ", 1)[0]
    else: 
        descr = "no_db_info"
    return descr