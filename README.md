# Progetto Bioinformatica
Lo script si comporta secondo le specifiche in `docs/SPECIFICHE.md`.

Riceve in input un file GTF (`.gtf`), un file FASTA (`.fa`) produce in output 3 file FASTA (`.fa`).

## Dipendenze
Lo script utilizza i moduli `argparse`, `os` e `sys`che sono già inclusi di default in python. **Per PYTHON è richiesta una VERSIONE 3+** (una distribuzione della v3 o superiore).

# Output
Lo script può generare tre diversi file:
* `$GENOME_$TIMESTAMP_transcripts.fa` : file contenente i trascritti...
```
>$ID_TRASCRITTO_1 gene_id=$ID_GENE_1 length=$len
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

>$ID_TRASCRITTO_2 gene_id=$ID_GENE_1 length=$len
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXX

>$ID_TRASCRITTO_N gene_id=$ID_GENE_Y length=$len
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX

#etc...
```
* `$GENOME_$TIMESTAMP_exome.fa` : file contenente l'esoma del genoma passato in input..
* `$GENOME_$TIMESTAMP_cc.fa` : file contenente le coperture codificanti del genoma...

Con `$GENOME` che è l'ID del genoma estratto dall'header del genoma passato in input e `$TIMESTAMP` è il timestamp al momento della creazione del file (è uno stratagemma per essere sicuri di avere file con nomi univoci ad ogni esecuzione del programma).

# Utilizzo primario dello script

Lo script funziona da riga di comando nel modo seguente:
```
$ python xtractor.py path/to/GENOME.fa path/to/INPUT.gtf [ARGS]
``` 
Gli argomenti `path/to/GENOME.fa` e `path/to/INPUT.gtf` sono obbligatori e corrispondono ai percorsi (relativi o assoluti) rispettivamente del file del genoma (in formato FASTA) e del file della gtf (in formato GTF). L'**ordine dei parametri** `GENOME.fa` e `INTPU.gtf` **è importante**!

Lo script accetta anche una serie di parametri opzionali "**ARGS**" in base alle esigenze di gestione e creazione dei file in output.

Si digiti 
```
$ python xtractor.py -h
```
o, equivalentemente
```
$ python xtractor.py --help
```
per un elenco dettagliato delle flag disponibili.

**NB**: come in tutti gli script da riga di comando i parametri opzionali possono essere inseriti in qualsiasi ordine.

## Args
```
$ python xtractor.py -h
usage: xtractor.py [-h] [-d OUTPUT_DIR] [-t] [-e] [-c] [-o]
                   genome_path gtf_path

GTF to FASTA extractor (for transcripts, exome and cc)

positional arguments:
  genome_path           path to genome file
  gtf_path              path to input file

optional arguments:
  -h, --help            show this help message and exit
  -d OUTPUT_DIR, --dir OUTPUT_DIR
                        Create files under specified directory
  -t, --transcripts     Create (only) transcripts file
  -e, --exome           Create (only) exome file
  -c, --cc              Create (only) cc file
  -o                    Print to standard output instead of creating files
                        (allowed only whether exactly one flag between -t, -e
                        or -c is used)
```

# Utilizzare i metodi interni allo script
Importare da console python il modulo, ad esempio:
```
$ cd path/to/xtractor.py
$ python
Python 3.6.1 (default, Mar 27 2017, 00:27:06)
>>> import xtractor as x
```

## Metodi interni

### `readFastaGenome(str genome_path)`
Legge il genoma collocato nel file `genome_path` (formato `.fa`) e restituisce un `dict` come il seguente:
```python
genome_sample = {
    "$GENOME_ID" : "$GENOME_CONTENT_AS_STRING_WITHOUT_WRAPS" 
}
```

### `readGtf(str gtf_path)`
Legge la gtf collocata nel file `gtf_path` (formato `.gtf`) e restituisce una struttura dati json-like come la seguente:
```python
gtf_sample = {
    "$GENOME_ID": {
        "$GENE_ID_1" : {
            "$TRANSCRIPT_ID_1" : {
                "exon" : [
                    { # exon 1
                        "start_index" : 781851,
                        "end_index" : 781951,
                        "score" : None,
                        "strand" : "+",
                        "frame" : None,
                        "attributes" : {
                            "attr1" : "val1",
                            # ...
                            "attrN" : "valN"
                        }
                    },
                    {
                        # exon 2
                    },
                    # ...
                    {
                        # exon N
                    }
                ],
                "CDS" : [
                    # ...
                ],
                "5UTR" : [
                    # ...
                ]
                # and all other features'types
            },
            # ...
            "$TRANSCRIPT_ID_N" : {
                # ...
            }
        },
        #...
        "$GENE_ID_N" : {
            "$TRANSCRIPT_ID_xx" : {
                "exon" : [
                    # ...
                ]
                # etc...
            }
            # ...
        }
    }
}
```
Per visualizzare nel dettaglio tale struttura si provi a lanciare i seguenti comandi all'interno della console python:
```
import sys
import json
import xtractor as x
file = open("nome_qualsiasi.txt","w")
sys.stdout = file
print(json.dumps(x.readGtf("./src/input.gtf")))
file.close()
exit()
# uscire per resettare lo stdout
```
Poi si copi il contenuto del file `nome_qualsiasi.txt` in [JSON Formatter](https://jsonformatter.curiousconcept.com/) per visualizzare la gtf in modo dettagliato e _human readable_.

### `grind(str genome_path, str gtf_path, str output_dir=os.getcwd(), bool t_flag=False, bool e_flag=False, bool c_flag=False)`
Procedura principale per l'analisi dei file in input. Non restituisce nulla ma termina generando in output i 3 file descritti all'inizio.

_finire specifica di `grind()`_

### `reverseAndComplement(str sequence)`
Accetta una stringa di basi e ne il complemento in ordine opposto.