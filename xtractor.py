import argparse
import os
import sys
import time
import json

# GLOBALS
# print on stdout instead of on files
STD_OUT = False

# wrap after XX characters in fasta files
FASTA_WRAP_AFTER_COLUMNS = 50


def readGtf(gtf_path):
    gtf = {}

    # generated gtf structure sample
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

    # open file passed as argument
    with open(gtf_path) as raw_gtf:
        # iterate over raw_gtf lines
        for record in raw_gtf:
            # split record on tabs (and remove quotes and newlines)
            # quotes are removed because they could be a cause for parsing problems
            # newlines are unnecessary
            field = record.replace("\"", "").replace("\n", "").split("\t")

            # extract record's fields
            genome = field[0]
            feature_name = field[2]
            start_index = field[3]
            end_index = field[4]
            strand = field[6]

            # convert "." to None for easier handling
            if field[7] == ".":
                frame = None
            else:
                frame = int(field[7])

            # extract attr field
            attr = field[8]

            gene_id = transcript_id = 0
            if attr[len(attr) - 1] == ";":
                # strip last ";" on attr field
                attr = attr[:len(attr) - 1]

            # spilt attr field on ";"
            attr = attr.split("; ")

            for i in range(0, len(attr)):
                # extract transcript_id and gene_id form attributes
                # others attributes don't matter for the project
                key, value = attr[i].split(" ")
                if key == "transcript_id":
                    transcript_id = value
                elif key == "gene_id":
                    gene_id = value

            # instantiate empty json-like structure
            if genome not in gtf:
                gtf[genome] = {}

            if gene_id not in gtf[genome]:
                gtf[genome][gene_id] = {}

            if transcript_id not in gtf[genome][gene_id]:
                gtf[genome][gene_id][transcript_id] = {}

            if feature_name not in gtf[genome][gene_id][transcript_id]:
                gtf[genome][gene_id][transcript_id][feature_name] = []

            # fill gtf field for current feature (actual gtf's record)
            gtf[genome][gene_id][transcript_id][feature_name].append({
                "start_index": int(start_index),
                "end_index": int(end_index),
                "strand": strand,
                "frame": frame
            })

    return gtf


def readFastaGenome(fa_path):
    """
    Reads genome file in fasta format.
    :param fa_path: Path to fasta file
    :return: Genome dictionary
    """
    genome = {}
    genome_id = ""

    """
    genome = {
        "$GENOME_ID" : $GENOME_STRING
    }    
    """

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


def grind(genome_path, gtf_path, output_dir=os.getcwd(), t_flag=False, e_flag=False, c_flag=False, s_flag=False):
    global STD_OUT, FASTA_WRAP_AFTER_COLUMNS

    flags_amount = 0

    if t_flag:
        flags_amount += 1
    if e_flag:
        flags_amount += 1
    if c_flag:
        flags_amount += 1
    if s_flag:
        flags_amount += 1

    # if none file flags is set, then set all of them to True to print all outputs
    if flags_amount == 0:
        t_flag = e_flag = c_flag = s_flag = True

    # if user wants the output to be streamed on the std out
    # but he/she chose more than one flag between -t, -e, -c and -s
    # then the script can no longer prints on stdout
    # so the user will be prompted to switch to file output
    if STD_OUT and flags_amount > 1:
        while True:
            # prompt the user to switch to file streaming
            i = input(
                "You chose to write onto the standard output (-o), but the stdout allows only one file (between -t, -e, -c and -s) at a time to be written on it. Write on files instead? (files will be created under the current folder) [Y/n]...\n").lower()
            if i == "y":
                # if choose YES then turn stdout off
                STD_OUT = False
                break
            elif i == "n":
                # if choose NO then quit the program
                sys.exit("User chose to abort the program")
            else:
                # else prompt the question again
                print("Wrong choice, please try again")

    # read gtf and genome from file
    gtf = readGtf(gtf_path)
    genome = readFastaGenome(genome_path)
    genome_id = list(genome.keys())[0]

    gene_strand_map = {}
    """
    gene_strand_map = {
        $GENE_1_ID : $GENE_1_STRAND, # here strand appears as +1 or -1
        $GENE_2_ID : $GENE_2_STRAND,
        #...
        $GENE_N_ID : $GENE_N_STRAND
    }
    """

    transcript_struct = {}
    """
    transcripts_struct= {
        $ID_GENOME : {
            $ID_GENE : {
                $ID_TRASCRITTO_1 : $TRASCRITTO_1_STRING,
                $ID_TRASCRITTO_2 : $TRASCRITTO_2_STRING,
                #...
                $ID_TRASCRITTO_N : $TRASCRITTO_N_STRING
            }
        }
    }
    """

    exon_struct = {}
    """
    exon_struct = {
        "$START_INDEX;$END_INDEX" : {
            "gene_id" : "$GENE_ID",
            "transcripts" : [
                "$TR_ID_1", "$TR_ID_2"
            ]
        }
        #...
    }
    """

    cc_struct = {}
    """
    cc_struct = {
        "$GENE_ID" : [
            # for each exon
            "$EXON_X1_START_INDEX;$EXON_Y1_END_INDEX",
            "$EXON_X2_START_INDEX;$EXON_Y2_END_INDEX",
            #...
            "$EXON_Xn_START_INDEX;$EXON_Yn_END_INDEX"
        ]
        #...
    }
    """

    cds_struct = {}
    """
    cds_struct = {
        $ID_GENE : {
            $ID_TRASCRITTO_1 : $CDS_1_SEQUENCE,
            $ID_TRASCRITTO_2 : $CDS_2_SEQUENCE,
            #...
            $ID_TRASCRITTO_N : $CDS_N_SEQUENCE
        }
    }
    """

    transcript_struct[genome_id] = {}
    cds_struct[genome_id] = {}

    # loop through every gene in the current genome
    for gene_id in gtf[genome_id]:
        transcript_struct[genome_id][gene_id] = {}
        cds_struct[genome_id][gene_id] = {}
        gene_strand_map[gene_id] = '?'

        # set unknown strand value for the current gene
        strand = '?'

        # loop through every transcript in the current gene
        for transcript_id in gtf[genome_id][gene_id]:
            # extract only exons features from current transcript
            exons = gtf[genome_id][gene_id][transcript_id]["exon"]

            # tmp dictionary for exons
            tmp_exons_dict = {}

            """
            tmp_exons_dict = {
                "$EXON_START_INDEX" : "$EXON_NUCLEOTIDES_SEQUENCE"
            }
            """

            # loop through every exon in the current transcript
            for exon in exons:

                # extract start and end index of current exon
                p0 = exon['start_index']
                pf = exon['end_index']

                if strand == '?':
                    # if strand is unknown set it as the strand
                    # of the first exon of the gene
                    strand = exon['strand']

                    # save gene strand in gene_strand_map
                    gene_strand_map[gene_id] = strand + str(1)

                elif strand != exon['strand']:
                    # we have to expect that all features of the same
                    # gene are read from the same genetic side, so if just one exon
                    # has a different strand abort the program with error message
                    raise Exception(
                        "Fatal error: exons have different strands for the same transcript",
                        "(in genome_id: %s, gene_id: %s, transcript_id: %s" %
                        (genome_id, gene_id, transcript_id)
                    )

                # if the p0 is not in the tmp_dict structure then the current exon
                # needs to be inserted into the tmp_dict structure
                # if p0 is already a key for tmp_dict it means that two differen
                # exons for the current transcript start at the same position
                if p0 not in tmp_exons_dict:
                    # save exon sequence in tmp_dict
                    tmp_exons_dict[p0] = genome[genome_id][p0 - 1:pf]
                else:
                    raise Exception(
                        "Fatal error: two different exons start at the same position",
                        "(in genome_id: %s, gene_id: %s, transcript_id: %s" %
                        (genome_id, gene_id, transcript_id)
                    )

                # "p0;pf" key for the current exon that will be used into exon_struct
                tmp_exon_struct_key = str(p0) + ";" + str(pf)
                if tmp_exon_struct_key not in exon_struct:
                    # if current exon is not in exon_struct then add it
                    exon_struct[tmp_exon_struct_key] = {
                        "gene_id": gene_id,
                        "transcripts": [
                            transcript_id
                        ]
                    }
                else:
                    # else add current transcript_id to its transcripts list
                    exon_struct[tmp_exon_struct_key]["transcripts"].append(transcript_id)

                # if not exists, instantiate array in cc_struct for the current gene
                if gene_id not in cc_struct:
                    cc_struct[gene_id] = []

                # append exon start_index and end_index as string "p0;pf" in cc_struct
                # array for the current gene
                if tmp_exon_struct_key not in cc_struct[gene_id]:
                    cc_struct[gene_id].append(tmp_exon_struct_key)

            # extract tmp_dict keys
            tmp_exons_dict_keys = list(tmp_exons_dict.keys())

            # sort tmp_exons_dict_keys list
            tmp_exons_dict_keys.sort()

            # current transcript nucleotides sequence
            curr_transcript = ""

            # read exons' sequence ordered by keys using tmp_exons_dict_keys list
            for key in tmp_exons_dict_keys:
                curr_transcript = curr_transcript + tmp_exons_dict[key]

            # reverse and complement transcript if it has negative strand
            if strand == '-':
                curr_transcript = reverseAndComplement(curr_transcript)

            # finally save transcript sequence in a transcript structure
            transcript_struct[genome_id][gene_id][transcript_id] = curr_transcript

            ################
            # start cds loop

            if "CDS" in gtf[genome_id][gene_id][transcript_id]:
                cdss = gtf[genome_id][gene_id][transcript_id]["CDS"]
                tmp_cdss_dict = {}

                for cds in cdss:
                    p0 = cds['start_index']
                    pf = cds['end_index']

                    # if the p0 is not in the tmp_dict structure then the current exon
                    # needs to be inserted into the tmp_dict structure
                    # if p0 is already a key for tmp_dict it means that two differen
                    # exons for the current transcript start at the same position
                    if p0 not in tmp_cdss_dict:
                        # save exon sequence in tmp_dict
                        tmp_cdss_dict[p0] = genome[genome_id][p0 - 1:pf]
                    else:
                        raise Exception(
                            "Fatal error: two different cds start at the same position",
                            "(in genome_id: %s, gene_id: %s, transcript_id: %s" %
                            (genome_id, gene_id, transcript_id)
                        )

                # extract tmp_dict keys
                tmp_cdss_dict_keys = list(tmp_cdss_dict.keys())

                # sort tmp_exons_dict_keys list
                tmp_cdss_dict_keys.sort()

                # current transcript nucleotides sequence
                curr_cds = ""

                # read exons' sequence ordered by keys using tmp_exons_dict_keys list
                for key in tmp_cdss_dict_keys:
                    curr_cds = curr_cds + tmp_cdss_dict[key]

                # reverse and complement transcript if it has negative strand
                if strand == '-':
                    curr_cds = reverseAndComplement(curr_cds)

                # finally save transcript sequence in a transcript structure
                cds_struct[genome_id][gene_id][transcript_id] = curr_cds


    # whole strings (file content) to output
    transcript_output = ""
    exome_output = ""
    cc_output = ""
    cds_output = ""

    #### format output for transcript file BEGIN
    if t_flag:
        # loop through all gene in current genome
        for gene_id in transcript_struct[genome_id]:
            # loop through all transcript found for the current gene
            for transcript_id in transcript_struct[genome_id][gene_id]:
                # extract trascript using current trascript_id
                t = transcript_struct[genome_id][gene_id][transcript_id]

                # append header for the current transcript
                # to transcript_output string
                # >/source=ENm006 /gene_id="ARHGAP4" /gene_strand=-1 /transcript_id=U52112.4-005 /length=642
                transcript_output += ">/source={0} /gene_id=\"{1}\" /gene_strand={2} /length={3} /transcript_id={4}\n".format(
                    genome_id,
                    gene_id,
                    gene_strand_map[gene_id],
                    len(t),
                    transcript_id
                )

                # variable to take into account how many characters
                # there are in the current line
                start_index = 0

                # loop through every character in current trascript sequence
                while start_index < len(t):
                    # if remaining cheracters to print are long enough to fill up
                    # an entire line FASTA_WRAP_AFTER_COLUMNS long
                    if start_index + FASTA_WRAP_AFTER_COLUMNS < len(t):
                        # then print a line FASTA_WRAP_AFTER_COLUMNS long
                        transcript_output += t[start_index:start_index + FASTA_WRAP_AFTER_COLUMNS] + "\n"
                    else:
                        # else print remaing characters in the current transcript sequence
                        transcript_output += t[start_index:] + "\n"
                        break
                    start_index += FASTA_WRAP_AFTER_COLUMNS

        # remove all trailing newlines at the end of the transcript output string
        while transcript_output[-1:] == "\n":
            transcript_output = transcript_output[0:-1]
    #### format output for transcript file END

    #### format output for exome file BEGIN
    if e_flag:
        # loop through every exon in exon_struct
        for exon_key, exon_val in exon_struct.items():
            # extract start_index and end_index placed in exon_key semicolon separated
            p0, pf = exon_key.split(";")
            p0 = int(p0)
            pf = int(pf)

            xseq = genome[genome_id][p0-1:pf]

            # append current exon header to output string
            exome_output += ">/source={0} gene_id=\"{1}\" /gene_strand={2} /length={3} /transcripts_id={4}\n".format(
                genome_id,
                exon_val["gene_id"],
                gene_strand_map[exon_val["gene_id"]],
                str(len(xseq)),
                "|".join(exon_val["transcripts"])
            )

            # variable to take into account how many characters
            # there are in the current line
            start_index = 0

            # loop through every character in current exon sequence
            while start_index < len(xseq):
                # if remaining cheracters to print are long enough to fill up
                # an entire line FASTA_WRAP_AFTER_COLUMNS long
                if start_index + FASTA_WRAP_AFTER_COLUMNS < len(xseq):
                    # then print a line FASTA_WRAP_AFTER_COLUMNS long
                    exome_output += xseq[start_index:start_index + FASTA_WRAP_AFTER_COLUMNS] + "\n"
                else:
                    # else print remaing characters in the current exon sequence
                    exome_output += xseq[start_index:] + "\n"
                    break
                start_index += FASTA_WRAP_AFTER_COLUMNS

        # remove all trailing newlines at the end of the exome output string
        while exome_output[-1:] == "\n":
            exome_output = exome_output[0:-1]
    #### format output for exome file END

    #### format output for cc file BEGIN
    if c_flag:
        # loop through every gene in cc_struct
        for gene_id in cc_struct:
            # temp sequence for cc of the current gene
            ccseq  = ""

            # extract current gene cc array
            gene_cc = cc_struct[gene_id]

            # sort it
            gene_cc.sort()

            # var to keep previuos start_index and end_index
            # in the following for cicle
            prev_p0 = -1
            prev_pf = -1

            # loop through every item in gene_cc
            for k in gene_cc:
                # split k with regards to ";"
                # reminder: k is format as "$start_index;$end_index"
                p0, pf = k.split(";")

                # convert them to int
                p0 = int(p0)
                pf = int(pf)

                # it's likely that following controls are redundant
                # or useless but I prefer redundant controls instead
                # of any risk of failure

                if p0 <= prev_pf:
                    if pf > prev_pf:
                        prev_pf = pf
                        continue
                else:
                    if abs(prev_pf - prev_p0) > 0:
                        ccseq += genome[genome_id][prev_p0-1:prev_pf]
                    prev_p0 = p0
                    if pf > prev_pf:
                        prev_pf = pf

            # cc_output header
            # >/source=ENm006 /gene_id="ARHGAP4" /gene_strand=-1 /transcript_id=U52112.4-005 /length=642
            cc_output += str.format(
                ">/source={0} /gene_id=\"{1}\" /gene_strand={2} /length={3}\n",
                genome_id,
                gene_id,
                gene_strand_map[gene_id],
                len(ccseq)
            )

            # variable to take into account how many characters
            # there are in the current line
            start_index = 0

            # loop through every character in current cc sequence
            while start_index < len(ccseq):
                # if remaining cheracters to print are long enough to fill up
                # an entire line FASTA_WRAP_AFTER_COLUMNS long
                if start_index + FASTA_WRAP_AFTER_COLUMNS < len(ccseq):
                    # then print a line FASTA_WRAP_AFTER_COLUMNS long
                    cc_output += ccseq[start_index:start_index + FASTA_WRAP_AFTER_COLUMNS] + "\n"
                else:
                    # else print remaing characters in the current cc sequence
                    cc_output += ccseq[start_index:] + "\n"
                    break
                start_index += FASTA_WRAP_AFTER_COLUMNS

        # remove all trailing newlines at the end of the cc output string
        while cc_output[-1:] == "\n":
            cc_output = cc_output[0:-1]
    #### format output for cc file END

    #### format output for cds file BEGIN
    if s_flag:
        # loop through every gene in cds_struct
        for gene_id in cds_struct:
            # temp sequence for cds of the current gene
            cdsseq  = ""

            # extract current gene cds array
            gene_cds = cds_struct[gene_id]

            # sort it
            gene_cds.sort()

            # var to keep previuos start_index and end_index
            # in the following for cicle
            prev_p0 = -1
            prev_pf = -1

            # loop through every item in gene_cds
            for k in gene_cds:
                # split k with regards to ";"
                # reminder: k is formatted as "$start_index;$end_index"
                p0, pf = k.split(";")

                # convert them to int
                p0 = int(p0)
                pf = int(pf)

                # it's likely that following controls are redundant
                # or useless but I prefer redundant controls instead
                # of any risk of failure

                if p0 <= prev_pf:
                    if pf > prev_pf:
                        prev_pf = pf
                        continue
                else:
                    if abs(prev_pf - prev_p0) > 0:
                        cdsseq += genome[genome_id][prev_p0-1:prev_pf]
                    prev_p0 = p0
                    if pf > prev_pf:
                        prev_pf = pf

            start_codon = cdsseq[0:3]
            stop_codon = cdsseq[-3:]

            # cds_output header
            cds_output += str.format(
                ">/source={0} /gene_id=\"{1}\" /gene_strand={2} /length={3} /start_codon={4} /stop_codon={5}\n",
                genome_id,
                gene_id,
                gene_strand_map[gene_id],
                len(cdsseq),
                start_codon,
                stop_codon
            )

            # variable to take into account how many characters
            # there are in the current line
            start_index = 0

            # loop through every character in current cc sequence
            while start_index < len(cdsseq):
                # if remaining cheracters to print are long enough to fill up
                # an entire line FASTA_WRAP_AFTER_COLUMNS long
                if start_index + FASTA_WRAP_AFTER_COLUMNS < len(cdsseq):
                    # then print a line FASTA_WRAP_AFTER_COLUMNS long
                    cds_output += cdsseq[start_index:start_index + FASTA_WRAP_AFTER_COLUMNS] + "\n"
                else:
                    # else print remaing characters in the current cc sequence
                    cds_output += cdsseq[start_index:] + "\n"
                    break
                start_index += FASTA_WRAP_AFTER_COLUMNS

        # remove all trailing newlines at the end of the cc output string
        while cds_output[-1:] == "\n":
            cds_output = cds_output[0:-1]
    #### format output for cc file END

    #### PRINT OUTPUTS
    if not STD_OUT:
        # print to file

        # add trailing slash to output_dir
        output_dir = os.path.join(output_dir, '')

        # check whether output_path exists
        if not os.path.exists(output_dir):
            raise Exception("You chose to create files under "
                + output_dir + " but that directory doesn't exist, "
                + "please retry with a new correct one.\n\nInstead, if you not had chose any directory, that is you not had "
                + "specified the -d flag, the program used the current directory by default. It sounds weird, "
                + "if the problem persists, avoid it trying specifying an existsing directory by using -d flag."
            )

        # file prefix is: "$GENOME_ID_$TIMESTAMP_"
        file_prefix = genome_id + "_" + str(int(time.time())) + "_"

        # if t_flag is on
        if t_flag:
            # create $GENOME_$TIMESTAMP_transcripts.fa file
            tmp_path = output_dir + file_prefix + "transcripts.fa"

            # create file only if it doesn't already exist
            if not os.path.exists(tmp_path):
                # create, write to and close file
                with open(tmp_path, 'w') as transcript_file:
                    transcript_file.write(transcript_output)
                    transcript_file.close()
                    print("Transcripts'file created at " + tmp_path)
            else:
                # current file already exists and this is crazy
                print("Error: trascripts'file " + tmp_path + " already exists and this behaviour is very weird")
                pass

        # if e_flag is on
        if e_flag:
            # create $GENOME_$TIMESTAMP_exome.fa file
            tmp_path = output_dir + file_prefix + "exome.fa"

            # create file only if it doesn't already exist
            if not os.path.exists(tmp_path):
                # create, write to and close file
                with open(tmp_path, 'w') as exome_file:
                    exome_file.write(exome_output)
                    exome_file.close()
                    print("Exome's file created at " + tmp_path)
            else:
                # current file already exists and this is crazy
                print("Error: exome's file " + tmp_path + " already exists and this behaviour is very weird")
                pass

        # if c_flag is on
        if c_flag:
            # create $GENOME_$TIMESTAMP_cc.fa file
            tmp_path = output_dir + file_prefix + "cc.fa"

            # create file only if it doesn't exist
            if not os.path.exists(tmp_path):
                # create, write to and close file
                with open(tmp_path, 'w') as cc_file:
                    cc_file.write(cc_output)
                    cc_file.close()
                    print("CC's file created at " + tmp_path)
            else:
                # current file already exists and this is crazy
                print("Error: CC's file " + tmp_path + " already exists and this behaviour is very weird")
                pass

        # if s_flag is on
        if s_flag:
            # create $GENOME_$TIMESTAMP_cds.fa file
            tmp_path = output_dir + file_prefix + "cds.fa"

            # create file only if it doesn't exist
            if not os.path.exists(tmp_path):
                # create, write to and close file
                with open(tmp_path, 'w') as cds_file:
                    cds_file.write(cds_output)
                    cds_file.close()
                    print("CDS' file created at " + tmp_path)
            else:
                # current file already exists and this is crazy
                print("Error: CDS' file " + tmp_path + " already exists and this behaviour is very weird")
                pass

        # if all flags are down
        if not (t_flag or e_flag or c_flag or s_flag):
            raise Exception(
                "Very weird error: the program was sure you chose to print something on files" +
                " but there is nothing to print :("
            )
    else:
        # print to std_out
        if t_flag:
            print(transcript_output)
        elif e_flag:
            print(exome_output)
        elif c_flag:
            print(cc_output)
        elif s_flag:
            print(cds_output)
        else:
            raise Exception(
                "Very weird error: the program was sure you chose to print on the Standard Output" +
                " with one flag (between -t, -e or -c) set, but all flags are down :("
            )


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
    # argparse setup
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
        '-s', '--cds',
        action='store_true',
        required=False,
        help='Create (only) cds file',
        dest='s_flag'
    )

    parser.add_argument(
        '-o',
        action='store_true',
        required=False,
        help='Print to standard output instead of creating files (allowed only whether exactly one flag between -t, -e or -c is used)',
        dest='std_out'
    )

    args = parser.parse_args()

    STD_OUT = args.std_out

    # run main procedure (grind(...))
    grind(
        args.genome_path,
        args.gtf_path,
        args.output_dir,
        args.t_flag,
        args.e_flag,
        args.c_flag,
        args.s_flag
    )
