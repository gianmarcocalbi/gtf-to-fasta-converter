# Progetto Bioinformatica

Lo script agisce secondo le specifiche in `docs/SPECIFICHE.md`.

Riceve in input un file GTF (.gtf), un file FASTA (.fa) produce in output 3 file FASTA (.fa).

## Output
Lo script può generare tre diversi file:
* `$GENOME_$TIMESTAMP_transcripts.fa` : file contenente i trascritti...
* `$GENOME_$TIMESTAMP_exome.fa` : file contenente l'esoma del genoma passato in input..
* `$GENOME_$TIMESTAMP_cc.fa` : file contenente le coperture codificanti del genoma...

Con `$GENOME` che è l'ID del genoma estratto dall'header del genoma passato in input e `$TIMESTAMP` è il timestamp al momento della creazione del file (è uno stratagemma per essere sicuri di avere file con nomi univoci ad ogni esecuzione del programma).

## Utilizzo

Lo script funziona da riga di comando nel modo seguente:
```
$ python xtractor.py path/to/GENOME.fa path/to/INPUT.gtf [ARGS]
``` 
Gli argomenti `path/to/GENOME.fa` e `path/to/INPUT.gtf` sono obbligatori e corrispondono ai percorsi (relativi o assoluti) rispettivamente del file del genoma (in formato FASTA) e del file della gtf (in formato GTF). L'**ordine dei parametri** `GENOME.fa` e `INTPU.gtf` **è importante**!
```
$ python xtractor.py path/to/INPUT.gtf path/to/GENOME.fa [ARGS]
``` 

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

### Args
```
$ python xtractor.py -h
-h, --help  Mostra queste informazioni di aiuto

-d, --dir <path/to/dir>    Crea i file nella cartella situata al percorso path/to/dor 

-t, --transcripts   Genera in output solo il file dei trascritti
-e, --exome         Genera in output solo il file dell'esoma
-c, --cc            Genera in output solo il file delle coperture codificanti
```

# Struttura

La GTF è tenuta in una struttura dati allocata per intero in memoria principale. Essa, in linea teorica, sarebbe così strutturata:

```python
gtf_sample = {
    "ENm006": {
        "$ID_GENE" : {
            "$ID_TRASCRITTO_1" : {
                "exon" : [
                    {
                        "start_index" : 781851,
                        "end_index" : 781951,
                        "score" : None,
                        "strand" : "+",
                        "frame" : None,
                        "attributes" : {
                            "attr1" : "val1",
                            "attr2" : "val2"
                        }
                    }
                ]
            },
            "$ID_TRASCRITTO_2" : {
                #...
            },
            "$ID_TRASCRITTO_3" : {
                #...
            }
        }
    }
}
```

Rispetto ad essa ho eliminato alcuni campi superflui ai fini del progetto ottenendo una struttura dati più snella naturalmente:

```python
gtf_sample = {
    "ENm006": {
        "$ID_GENE" : {
            "$ID_TRASCRITTO_1" : {
                "$FEATURE_TYPE" : [
                    {
                        "start_index" : 781851,
                        "end_index" : 781951,
                        "strand" : "+",
                        "frame" : None
                    }
                ]
            },
            "$ID_TRASCRITTO_2" : {
                "exon" : [
                    #...
                ],
                "CDS" : [
                    #...
                ],
                "5UTR" : [
                
                ]
                # etc per tutte le altre features
            },
            "$ID_TRASCRITTO_3" : {
                #...
            }
        },
        # altri geni ...
    }
}
```

## Interfaccia "interna"

### readGtf(genome_path, gtf_path, output_dir, t_flag, e_flag, c_flag)


