# Specifiche

## Input
Lo script Python deve prendere in input un **file GTF** contenente diversi geni strutturato come il file di esempio nella folder `src`.

## Output
Lo script produce in output **tre file in formato FASTA** specificatin in seguito.

## File .fa delle sequenze dei trascritti
File FASTA delle sequenze dei trascritti ricostruiti per tutti i geni contenuti nell'input. 

Nell’header FASTA di ogni trascritto devono comparire le seguenti informazioni: 
* lunghezza
* id del gene a cui il trascritto appartiene
* id del trascritto.

## File .fa dell'esoma
File FASTA dell’esoma, ovvero tutte le sequenze degli esoni contenuti nel GTF (ovviamente se un esone compare in due trascritti diversi deve comparire una volta sola nell’output).

Nell’header FASTA di ogni esone devono comparire le seguenti informazioni:
* lunghezza
* id del gene a cui l’espone appartiene
* lista degli id dei trascritti a cui l’espone appartiene

## File .fa delle coperture codificanti
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