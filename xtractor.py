import json

# verbous FLAG
V = False

# Input PATH
GTF_PATH = "/home/grimidev/documents/dev-workspace/progetto-bioinformatica/src/input.gtf"

# Current absoulte path
CURRENT_PATH = None

# Path where create FASTA files
OUTPUT_PATH = None

# Big GTF data structure
gtf = {}

# FASTA GENOME data structure
fasta = {}

"""
gtf_sample = {
    "ENm006": {
        "$ID_GENE": {
            "$ID_TRASCRITTO_1" : {
                "exon" : [
                    {
                        "start_index" : 781851,
                        "end_index" : 781951,
                        "strand" : "+",
                        "frame" : None
                    }
                ]
            },
            "$ID_TRASCRITTO_2" : {
                #...
            },
            "$ID_TRASCRITTO_3" : {
                #...
            }
        },
        "$ID_GENE2" : {
            #...
        }
        # etc..
    }
}
"""


class xtractor:
    pass

def readGtf(gtf_path):
    tmp_gtf = {}
    with open(gtf_path) as raw_gtf:
        for record in raw_gtf:
            field = record.replace("\"", "").replace("\n", "").split("\t")
            genome = field[0]
            feature_name = field[2]
            start_index = field[3]
            end_index = field[4]
            strand = field[6]
            if field[7] == ".":
                frame = None
            else:
                frame = int(field[7])
            attr = field[8]
            gene_id = transcript_id = 0
            if attr[len(attr)-1] == ";":
                attr = attr[:len(attr)-1]
            attr = attr.split("; ")
            for i in range(0,len(attr)):
                key, value = attr[i].split(" ")
                if key == "transcript_id":
                    transcript_id = value
                elif key == "gene_id":
                    gene_id = value

            if genome not in tmp_gtf:
                tmp_gtf[genome] = {}

            if gene_id not in tmp_gtf[genome]:
                tmp_gtf[genome][gene_id] = {}

            if transcript_id not in tmp_gtf[genome][gene_id]:
                tmp_gtf[genome][gene_id][transcript_id] = {}

            if feature_name not in tmp_gtf[genome][gene_id][transcript_id]:
                tmp_gtf[genome][gene_id][transcript_id][feature_name] = []

            tmp_gtf[genome][gene_id][transcript_id][feature_name].append({
                "start_index" : int(start_index),
                "end_index" : int(end_index),
                "strand" : strand,
                "frame" : frame
            })

    return tmp_gtf

def readFastaGenome(fa_path):
    genome = {}
    genome_id = ""
    # Open and read fasta file located in fa_path
    with open(fa_path) as raw_fa:
        tmp_str = ""
        # read every single line till EOF
        for line in raw_fa:
            # if line is a fasta header
            if line[0] == ">":
                # extract genome name (id)
                genome_id = line[1:].replace("\n", "")
            else:
                # append nucleotides'sequence
                tmp_str = tmp_str + (line.replace("\n", ""))

        genome[genome_id] = tmp_str
    return genome


def build():
    global gtf

    for genome in gtf:
        for gene in genome:
            for transcript in gene:
                #
                pass


    pass

def reverseAndComplement (str):
    
    complementDictionary = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    
    newstr = []
    for i in range(0, len(str)):
        c = str[i]

        if c in complementDictionary:
            newstr.append(complementDictionary[c])
        else:
            print("fatal error in complementing transcript: found an unknown nucleobasis")
            return 

    # reverse
    return "".join(newstr)[::-1]

if __name__ == '__main__':
    pass
