# Specifiche Progetto

## Input
Lo script Python deve prendere in input un **file GTF** contenente diversi geni strutturato come il file di esempio nella folder `src`.

### Formato GTF
File di puro testo con estensione .gtf; contiene record (righe separate da wraps). Ogni riga può far riferimento ad una delle seguenti feature:
* un esone;
* una  (coding sequence) o una sua parte;
* un 5' UTR (untranslated region) o una sua parte;
* un 3' UTR o una sua parte;
* uno start codon;
* uno stop codon.

Ogni record della GTF è composto sempre da 9 campi separati da tabulazioni, se un campo ha valore `.` (dot) allora è da considerarsi vuoto.

#### Campo #1 - ID Genomica
Il primo campo contiene l'ID della genomica di riferimento. Solitamente un file GTF contiene campi tutti riferiti alla stessa genomica ma questa caratteristica non p da prendersi come uno standard.

#### Campo #2 - Autore annotazione
Il secondo campo contiene la sorgente (autore) che ha prodotto l'annotazione. Può trattarsi di un software come di una persona.

#### Campo #3 - Nome feature
Il terzo campo contiene il nome della feature a cui fa riferimento il record attuale. Questo campo può assumere SOLO uno tra i seguenti 6 valori: 
```
[
    "exon",
    "CDS",
    "5UTR",
    "3UTR",
    "start_codon",
    "stop_codon"
]
```

#### Campo #4 - Posizione inizio feature
Il quarto campo contiene la posizione di inizio della feature corrente sulla genomica di riferimento. La genomica di riferimento è espressa nel campo #1.

#### Campo #5 - Posizione fine feature
Il quinto campo contiene la posizione di fine della feature corrente sulla genomica di riferimento.

_**Capire se l'indicizzazione delle posizioni inizia da 0 o 1 e se la posizione di fine è inclusa o no.**_

#### Campo #6 - Score feature
Indica lo score assegnato a tale feature (spesso è presente il dot per indicare che non vi è uno score assegnato alla feature corrente).

#### Campo #7 - Strand
Il settimo campo contiene lo strand. Può assumere valore:
* `+` : il frammento è stato preso dalla catena di trascrizione; 
* `-` : il frammento è stato preso dalla catena opposta a quella di trascrizione.

Per ottenere la sequenza di una feature con strand `-` è necessario:
1) estrarre la sottostringa della genomica di riferimento che corrisponde alla feature;
2) eseguire un'operazione di reverse&complement della sottostringa estratta.

#### Campo #8 - Frame
L'ottavo campo contiene il frame, esso è specificato solo per feature di tipo CDS, start e stop codon altrimenti è presente un dot. Il frame può assumere valori in {0,1,2}. Il frame indica se una tripla si ritrova spezzata sulla genomica di riferimento. Infatti è possibile che sequenze che sul gene sono "attaccate" risultino "spezzate" (dagli introni) sulla genomica di riferimento.

#### Campo #9 - Lista attributi
Il non campo contiene una lista di attributi (coppie chiave-valore) così formattata:
```
<attribute_name1> <value1>; <attribute_name2> <value2>; …
```
* uno spazio separa la chiave dal valore;
* il valore pare essere sempre di tipo stringa (doppi apici);
* due coppie chiave-valore sono separate sempre da punto e virgola + spazio (`; `);
* **c'è il punto e virgola dopo l'ultima coppia** chiave-valore.

Gli attributi obbligatori (sempre presenti) sono:
* `"gene_id"`;
* `"transcript_id"`.


## Output
Lo script produce in output **quattro file in formato FASTA** specificati in seguito.

### File .fa delle sequenze dei trascritti
File FASTA delle sequenze dei trascritti ricostruiti per tutti i geni contenuti nell'input. 

Nell’header FASTA di ogni trascritto devono comparire le seguenti informazioni: 
* lunghezza
* id del gene a cui il trascritto appartiene
* id del trascritto.

### File .fa dell'esoma
File FASTA dell’esoma, ovvero tutte le sequenze degli esoni contenuti nel GTF (ovviamente se un esone compare in due trascritti diversi deve comparire una volta sola nell’output).

Nell’header FASTA di ogni esone devono comparire le seguenti informazioni:
* lunghezza
* id del gene a cui l’esone appartiene
* lista degli id dei trascritti a cui l’esone appartiene

### File .fa delle coperture codificanti
File FASTA delle “coperture codificanti” di ciascun gene.

Nell’header FASTA di ogni copertura devono comparire le seguenti informazioni:
* lunghezza della copertura
* id del gene

### Copertura Codificante
La copertura codificante è la sequenza di basi tale per cui ogni base appartiene a una regione codificante del gene.

#### Esempio
Supponi di avere un gene con i seguenti tre trascritti:
- trascritto 1 —> concatenazione degli esoni A, B e C
- trascritto 2 —> concatenazione degli esoni A’, B e C dove A’ è un prefisso di A
- trascritto 3 —> concatenazione degli esoni A e C
- trascritto 4 —> concatenazione degli esoni A, D e C dove D non ha sovrapposizione con B e viene dopo B

La copertura sarà data dalla concatenazione di A, B, D e C

### File .fa delle sequenze delle CDS
File FASTA delle sequenze delle CDS ricostruite per ciascun trascritto nella GTF. 

Nell’header FASTA di ogni trascritto devono comparire le seguenti informazioni: 
* lunghezza
* id del gene a cui il trascritto appartiene
* id del trascritto
* start e stop codon.