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


def loadGTF():
    global gtf
    with open(GTF_PATH) as raw_gtf:
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

            if genome not in gtf:
                gtf[genome] = {}

            if gene_id not in gtf[genome]:
                gtf[genome][gene_id] = {}

            if transcript_id not in gtf[genome][gene_id]:
                gtf[genome][gene_id][transcript_id] = {}

            if feature_name not in gtf[genome][gene_id][transcript_id]:
                gtf[genome][gene_id][transcript_id][feature_name] = []

            gtf[genome][gene_id][transcript_id][feature_name].append({
                "start_index" : int(start_index),
                "end_index" : int(end_index),
                "strand" : strand,
                "frame" : frame
            })
    #print(json.dumps(gtf))

def build():
    global gtf

    for genome in gtf:
        for gene in genome:
            for transcript in gene:
                #
                pass


    pass

def reverseAndComplement (str):
    """
    complement = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    """
    newstr = []
    for i in range(0, len(str)):
        c = str[i]

        if c == "A":
            newstr.append("T")
        elif c == "C":
            newstr.append("G")
        elif c == "G":
            newstr.append("C")
        elif c == "T":
            newstr.append("A")
        else:
            print("fatal error in complementing transcript")

    # reverse
    return "".join(newstr)[::-1]

if __name__ == '__main__':
    loadGTF()
    build()
    pass
