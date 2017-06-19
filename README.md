# Progetto Bioinformatica

Lo script agisce secondo le specifiche in `docs/SPECIFICHE.md`.

Riceve in input un file GTF (.gtf), produce in output 3 file FASTA (.fa).

## Diversi tipi di utilizzo

Lo script funziona da riga di comando (`$`).

# MANCA LA GENOMICA!!!

### Tutti i file nella cartella corrente
```bash
$ pyhton xtractor.py FILE.gtf
```
Crea i tre file nella cartella corrente chiamati:
* transcripts_TIMESTAMP.fa
* exome_TIMESTAMP.fa
* cc_TIMESTAMP.fa

### Tutti i file in una cartella specifica
```bash
$ pyhton xtractor.py FILE.gtf PATH/
```
Crea i tre file nella cartella PATH/ chiamati:
* transcripts_TIMESTAMP.fa
* exome_TIMESTAMP.fa
* cc_TIMESTAMP.fa

## Aggiungere un po' di flags a caso

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

### readGtf(gtf_path)


