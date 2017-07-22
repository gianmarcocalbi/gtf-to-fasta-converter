# Progetto Bioinformatica
Lo script si comporta secondo le specifiche in `docs/SPECIFICHE.md`.

In generale, riceve in input un file GTF (`.gtf`) e un file FASTA (`.fa`) e produce in output 4 file FASTA (`.fa`).

## Dipendenze
Lo script utilizza i moduli `argparse`, `os`, `sys` e `time`. Essi **sono già inclusi** di default in python.

**È richiesta una VERSIONE 3+ di Python** (una distribuzione della v3 o superiore).

# Output
Lo script può generare i seguenti 4 diversi file.

_Con `$GENOME` che è l'ID del genoma estratto dall'header del genoma passato in input e `$TIMESTAMP` è il timestamp al momento della creazione del file (è uno stratagemma per essere sicuri di avere file con nomi univoci ad ogni esecuzione del programma)._

### `$GENOME_$TIMESTAMP_transcripts.fa`
File contenente i trascritti. Ad esempio:
```
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN /transcript_id=$TRASCRIPT_1_ID
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN /transcript_id=$TRASCRIPT_2_ID
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_Y_ID" /gene_strand=$GENE_Y_STRAND /length=$LEN /transcript_id=$TRASCRIPT_X_ID
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
#etc...
```

In cui, per ogni trascritto, vengono specificati:
- `$GENOME_ID` : id univoco del genoma di riferimento;
- `$GENE_N_ID` : id del gene a cui appartiene il trascritto (**NB**: **questo campo e solo questo campo** è circondato da doppi apici);
- `$GENE_N_STRAND` : strand del gene (`+1` o `-1`);
- `$TRASCRIPT_N_ID` : id del trascritto;
- `$LEN` : lunghezza della sequenza del trascritto.

### `$GENOME_$TIMESTAMP_exome.fa`
File contenente l'esoma del genoma passato in input. Ad esempio:
```
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN /transcripts=$id_tr1|$id_tr5|...
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_2_ID" /gene_strand=$GENE_2_STRAND /length=$LEN /transcripts=$id_tr11|$id_tr25|...
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
#etc...
```

In cui, per ogni esone, vengono specificati:
- `$GENOME_ID` : id univoco del genoma di riferimento;
- `$GENE_N_ID` : id del gene a cui appartiene l'esone (**NB**: **questo campo e solo questo campo** è circondato da doppi apici);
- `$GENE_N_STRAND` : strand del gene dell'esone (`+1` o `-1`);
- `$LEN` : lunghezza della sequenza dell'esone;
- `$id_trX|$id_trY|$id_trZ` : lista degli id dei trascritti in cui compare l'esone in questione, ogni id è separato da una barra verticale (`|`).


### `$GENOME_$TIMESTAMP_cc.fa`
File contenente le coperture codificanti del genoma. Ad esempio:
```
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_2_ID" /gene_strand=$GENE_2_STRAND /length=$LEN
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
#etc...
```

In cui, per ogni copertura codificante, vengono specificati:
- `$GENOME_ID` : id univoco del genoma di riferimento;
- `$GENE_N_ID` : id del gene a cui è riferita la copertura codificante (**NB**: **questo campo e solo questo campo** è circondato da doppi apici);
- `$GENE_N_STRAND` : strand del gene (`+1` o `-1`);
- `$LEN` : lunghezza della sequenza della copertura.


### `$GENOME_$TIMESTAMP_cds.fa`
File contenente le CDS di ciascun trascritto. Ad esempio:
```
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN /transcript_id=$TRASCRIPT_X_ID /start_codon=XXX /stop_codon=YYY
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_1_ID" /gene_strand=$GENE_1_STRAND /length=$LEN /transcript_id=$TRASCRIPT_Y_ID /start_codon=XXX /stop_codon=YYY
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXX
>/source=$GENOME_ID /gene_id="$GENE_Y_ID" /gene_strand=$GENE_Y_STRAND /length=$LEN /transcript_id=$TRASCRIPT_Z_ID /start_codon=XXX /stop_codon=YYY
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXX
#etc...
```

In cui, per ogni CDS, vengono specificati:
- `$GENOME_ID` : id univoco del genoma di riferimento;
- `$GENE_N_ID` : id del gene a cui appartiene la CDS (**NB**: **questo campo e solo questo campo** è circondato da doppi apici);
- `$GENE_N_STRAND` : strand del gene (`+1` o `-1`);
- `$LEN` : lunghezza della sequenza della CDS;
- `$TRASCRIPT_N_ID` : id del trascritto a cui appartiene la CDS;
- `start_codon` : codone di inizio della CDS;
- `stop_codon` : codone di fine della CDS.


_**NB**: il formato degli header dei file rispecchia lo standard utilizzato nelle banche dati._

# Utilizzo primario dello script

Lo script funziona da riga di comando nel modo seguente:
```
$ python xtractor.py path/to/GENOME.fa path/to/INPUT.gtf [ARGS]
``` 
Gli argomenti `path/to/GENOME.fa` e `path/to/INPUT.gtf` sono obbligatori e corrispondono ai percorsi (relativi o assoluti) rispettivamente del file del genoma (in formato FASTA) e del file della gtf (in formato GTF). L'**ordine dei parametri** `GENOME.fa` e `INPUT.gtf` **è importante**!

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
usage: xtractor.py [-h] [-d OUTPUT_DIR] [-t] [-e] [-c] [-s] [-o]
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
  -s, --cds             Create (only) cds file
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

#### `readFastaGenome(str genome_path)`
Legge il genoma collocato nel file `genome_path` (formato `.fa`) e restituisce un `dict` come il seguente:
```python
genome_sample = {
    "$GENOME_ID" : "$GENOME_CONTENT_AS_STRING_WITHOUT_WRAPS" 
}
```

#### `readGtf(str gtf_path)`
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

#### `grind(str genome_path, str gtf_path, str output_dir=os.getcwd(), bool t_flag=False, bool e_flag=False, bool c_flag=False, bool s_flag=False)`
Procedura principale per l'analisi dei file in input. Non restituisce nulla.

Riceve in input due stringhe obbligatorie:
- `genome_path` : percorso del genoma;
- `gtf_path` : percorso della gtf.

Accetta ulteriori 4 parametri opzionali:
- `output_dir` : percorso di una cartella ESISTENTE del sistema in cui verranno creati i file di output;
- `t_flag` : valore booleano, se True allora crea il file dei trascritti, altrimenti no;
- `e_flag` : valore booleano, se True allora crea il file dell'esoma, altrimenti no;
- `c_flag` : valore booleano, se True allora crea il file delle coperture codificanti, altrimenti no;
- `s_flag` : valore booleano, se True allora crea il file delle CDS, altrimenti no.


**NB**: i percorsi possono essere sia assoluti sia relativi (in tal caso devono far riferimento alla cartella in cui è contenuto il file `xtractor.py`). 

La procedura `grind(...)` produce gli output specificati ad inizio documento.

#### `reverseAndComplement(str sequence)`
Accetta una stringa di basi e ne restituisce il complemento in ordine opposto.

___
**EOF**