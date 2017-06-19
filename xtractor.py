import argparse
import os

# verbose FLAG
V = False

# print on stdout instead of on files
STD_OUT = False

# wrap after XX characters in fasta files
FASTA_WRAP_AFTER_COLUMNS = 50

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
            if attr[len(attr) - 1] == ";":
                attr = attr[:len(attr) - 1]
            attr = attr.split("; ")
            for i in range(0, len(attr)):
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
                "start_index": int(start_index),
                "end_index": int(end_index),
                "strand": strand,
                "frame": frame
            })

    return tmp_gtf


def readFastaGenome(fa_path):
    """
    Reads genome file in fasta format.
    :param fa_path: Path to fasta file
    :return: Genome dictionary
    """
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


def grind(genome_path, gtf_path, output_dir=os.getcwd(), t_flag=False, e_flag=False, c_flag=False):
    gtf = readGtf(gtf_path)
    genome = readFastaGenome(genome_path)
    transcript = {}

    """
    transcripts = {
        $ID_GENOME : {
            $ID_GENE : {
                $ID_TRASCRITTO : $TRASCRITTO_STRING
            }
        }
    }
    """

    for genome_id in gtf:
        transcript[genome_id] = {}
        for gene_id in gtf[genome_id]:
            transcript[genome_id][gene_id] = {}
            for transcript_id in gtf[genome_id][gene_id]:
                # extract only exons features from current transcript
                exons = gtf[genome_id][gene_id][transcript_id]["exon"]
                tmp_dict = {}
                strand = '?'
                for exon in exons:
                    p0 = exon['start_index']
                    pf = exon['end_index']

                    if strand == '?':
                        strand = exon['strand']
                    elif strand != exon['strand']:
                        raise Exception(
                            "Fatal error: exons have different strands for the same transcript",
                            "(in genome_id: %s, gene_id: %s, transcript_id: %s" %
                            (genome_id, gene_id, transcript_id)
                        )

                    if p0 not in tmp_dict:
                        tmp_dict[p0] = genome[genome_id][p0:pf + 1]
                    else:
                        raise Exception(
                            "Fatal error: two different exons start at the same position",
                            "(in genome_id: %s, gene_id: %s, transcript_id: %s" %
                            (genome_id, gene_id, transcript_id)
                        )
                tmp_dict_keys = list(tmp_dict.keys())
                tmp_dict_keys.sort()
                curr_transcript = ""
                for key in tmp_dict_keys:
                    curr_transcript = curr_transcript + tmp_dict[key]
                if strand == '-':
                    curr_transcript = reverseAndComplement(curr_transcript)

                transcript[genome_id][gene_id][transcript_id] = curr_transcript

    # print(json.dumps(transcript))

    output = ""

    for genome_id in transcript:
        for gene_id in transcript[genome_id]:
            for transcript_id in transcript[genome_id][gene_id]:
                t = transcript[genome_id][gene_id][transcript_id]
                output = output + ">{0} gene_id={1} length={2}\n".format(transcript_id, gene_id, len(t))
                start_index = 0
                while start_index < len(t):
                    if start_index + FASTA_WRAP_AFTER_COLUMNS < len(t):
                        output = output + t[start_index:start_index + FASTA_WRAP_AFTER_COLUMNS] + "\n"
                    else:
                        output = output + t[start_index:] + "\n"
                        break
                    start_index = start_index + FASTA_WRAP_AFTER_COLUMNS
                output = output + "\n"

    if output[-2:] == "\n\n":
        output = output[0:-1]

    if not STD_OUT:
        #print to file
        pass
    else:
        print(output)


def reverseAndComplement(s):
    complementDictionary = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

    newstr = []
    for i in range(0, len(s)):
        c = s[i]

        if c in complementDictionary:
            newstr.append(complementDictionary[c])
        else:
            print("fatal error in complementing transcript: found an unknown nucleobasis")
            return

            # reverse
    return "".join(newstr)[::-1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='GTF to FASTA extractor (for transcripts, exome and cc)'
    )

    parser.add_argument(
        'genome_path',
        action='store',
        help='path to genome file'
    )

    parser.add_argument(
        'gtf_path',
        action='store',
        help='path to input file'
    )

    parser.add_argument(
        '-d', '--dir',
        action='store',
        default=os.getcwd(),
        required=False,
        help='Create files under specified directory',
        dest='output_dir'
    )

    parser.add_argument(
        '-t', '--transcripts',
        action='store_true',
        required=False,
        help='Create (only) transcripts file',
        dest='t_flag'
    )

    parser.add_argument(
        '-e', '--exome',
        action='store_true',
        required=False,
        help='Create (only) exome file',
        dest='e_flag'
    )

    parser.add_argument(
        '-c', '--cc',
        action='store_true',
        required=False,
        help='Create (only) cc file',
        dest='c_flag'
    )

    parser.add_argument(
        '-o',
        action='store_true',
        required=False,
        help='Print to standard output instead of creating files (allowed only whether exactly one flag between -t, -e or -c is used) nb: if -o is active none debug message could be printed',
        dest='std_out'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        required=False,
        help='Turn on debug printings (useless if -o flag is True)',
        dest='verbose'
    )

    args = parser.parse_args()

    STD_OUT = args.std_out
    V = args.verbose

    grind(
        args.genome_path,
        args.gtf_path,
        args.output_dir,
        args.t_flag,
        args.e_flag,
        args.c_flag
    )