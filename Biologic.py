import logging
from re import A
from Bio import SeqIO, Entrez
import pathlib

import GUI

log = logging.getLogger(__name__)

def get_bonds():
    """ Gibt pro Paar die Anzahl der Wasserstoffbrücken zurück.
    Ab Aufgabe 6.
    """

class biologic:
    ALPHABETS = {
        "dna": "ACTG_",
        "rna": "AUGC_",
        "protein": None     # protein alphabet is added in __init__()
    }
    # transcription
    DNA_TO_RNA = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    # reverse transcription
    RNA_TO_DNA = {v: k for k,v in DNA_TO_RNA.items()}


    def __init__(self):
        """TODO: modify doc string! Konstruktor. Falls nötig, können hier Voreinstellungen gemacht werden. """
        # create translation table
        table_fp = GUI.SCRIPT_DIR/"data"/"genCode.txt"
        if not table_fp.exists():
            log.critical(f'File for codon table not found at {table_fp}')
            raise FileNotFoundError(f'No codon table file fount at {table_fp}')
        
        codon_table = {}
        with open(table_fp, "r") as ct_file:
            for line in ct_file.readlines():
                line.strip()
                split = line.split()
                aa = split[0]
                for codon in split[1:]:
                    codon_table[codon] = aa
        self.codon_table = codon_table
        log.info(f'Read codon table from {table_fp}')
        log.debug(f'Created codon table: {codon_table}')

        self.ALPHABETS['protein'] = list(codon_table.values())
        log.debug(f'Protein alphabet added to `ALPHABETS`')


    def set_dna(self,dna):
        """ Setze die interne DNA-Sequenz. Aufgabe 1"""
        self.dna = dna.upper()
        log.info(f"New DNA entered: {dna}")

    def validate_dna(self):
        """ Überprüfe die intern gespeicherte DNA. Aufgabe 1

        Checks if all letters of DNA sequence are valid according to 
        the alphabet.
        """
        if all(i in self.ALPHABETS['dna'] for i in self.dna):
            log.debug(f'DNA is valid. DNA: {self.dna}')
            return True
        else:
            log.warning(f"DNA is invalid. DNA: {self.dna}")
            return False


    def get_rna(self):
        """ Die intern gespeicherte RNA wird zurück gegeben. Aufgabe 1
        
        Returns:
            RNA sequence as string.
        """
        return self.rna
 

    def read_fasta_file(self, fp: str):
        """ Lese eine DNA-Sequenz aus einer FASTA-Datei und speichere sie intern ab. Aufgabe 2.
        
        Argument:
            fp: file path of fasta file
        
        Returns:
            str: DNA sequence from fasta file
        """
        with open(fp, "r") as file:
            # take first record from fasta file
            record = next(SeqIO.parse(file, "fasta"), None)

            if record is None:
                # no record was read
                log.warning(f'No record was read from fasta file at location {fp}')
            else:
                self.dna = str(record.seq).upper()
                log.info(f'DNA sequence was successfully read from fasta file: {self.dna}')
                return self.dna


    def read_from_ncbi(self, id: str, mail: str):
        """ Lese eine DNA-Sequenz von der NCBI-Nucleotid Datenbank. Aufgabe 2.
        
        Arguments:
            id: Accession Number der Sequenz
            mail: die zu verwendende email-Adresse für den Entrez-Zugriff
        
        Returns:
            str: DNA sequence from NCBI
        """
        Entrez.email = mail
        log.debug(f'{mail} set for Entrenz call.')
        handle = Entrez.efetch(
            db='nucleotide', 
            id=id, 
            rettype='fasta', 
            retmode='text'
        )
        record = SeqIO.read(handle, 'fasta')
        handle.close()

        log.debug(f'Fasta file for id {id} fetched from NCBI. {record}')
        self.dna = str(record.seq).upper()
        log.info(f'DNA sequence for {id} fetched from NCBI.')
        return self.dna


    def transcribe(self):  
        """ Transkribiere die intere DNA-Sequenz in RNA. Dabei sei die DNA der codogene Strang.
        Die RNA-Sequenz soll intern gespeichert und zurück gegeben werden. Aufgabe 3
        
        Returns:
            str - rna sequence
        """
        self.rna = "".join(self.DNA_TO_RNA[i] for i in self.dna)
        log.info(f'DNA transcribed to RNA: {self.rna}')
        return self.rna


    def translate(self,initPos:int = 0):    
        """ Translatiere die intern gespeicherte RNA in eine Aminosäurensequenz, ab der 
        vorgegebenen Position. Aufgabe 4

        Argument:
            initPos - starting position for translation

        Returns:
            :str: peptid sequence
        """
        peptid = []
        log.debug(f'Start translation.')
        for i in range(initPos, len(self.rna), 3):
            codon = self.rna[i:i+3]
            if len(codon) < 3:
                # No complete codon left
                continue
            
            try: 
                peptid.append(self.codon_table[codon])
            except KeyError:
                log.error(f'Codon {codon} not found in codon table.')
                peptid.append(f'?({codon})')
        
        self.protein = "".join(peptid)
        log.info(f'Protein sequence created. {self.protein}')
        return self.protein


    def get_proteins(self, rf_number):
       """ Erzeuge für einen vorgegebenen Leserahmen alle daraus ablesbaren Proteine. Ein Protein beginnt
       mit 'M' und endet mit einem Stop '_'. Die Leserahmen haben die Nummern 0 bis 5. Dabei sind 0-2 die Rahmen
       in der vorgegebenen Richtung und 3-5 die Rahmen der revertierten Sequenz. Aufgabe 5.
       Parameter: rf_number: Nummer des Leserahmens (reading frame)
       Ergebnis: die Liste der Proteine (Liste von Strings)
       """

    def eval_rna(self, struct):
        """ Berechne für eine vorgegebene Struktur die Anzahl der enthaltenen Wasserstoffbrücken.
        Aufgabe 6.
        Parameter: struct: eine vorgegebene Sekundärstruktur in Klammerschreibweise.
        Ergebnis: Anzahl der Wasserstoffbrücken
        """

    def foldRna(self) :
        """ Berechne die optimale Sekundärstruktur nach dem vorgegebenen Verfahren. Aufgabe 7.
        Ergebnis: die optimale Anzahl von Wasserstoffbrücken
        """

   
    def trace(self) :
        """ Berechne für die zuletzt durchgeführte Faltung die optimale Struktur. Aufgabe 8.
        Ergebnis: optimale Struktur in Klammerschreibweise.
        """
