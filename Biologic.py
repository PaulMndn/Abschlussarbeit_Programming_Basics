import logging
log = logging.getLogger(__name__)

def get_bonds():
    """ Gibt pro Paar die Anzahl der Wasserstoffbrücken zurück.
    Ab Aufgabe 6.
    """

class biologic:
    ALPHABETS = {
        "dna": "ACTG_",
        "rna": "AUGC_",
        "protein": "ACDEFGHIKLMNPQRSTVWY_"
    }

    def __init__(self):
        """ Konstruktor. Falls nötig, können hier Voreinstellungen gemacht werden. """

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
 

    def read_fasta_file(self, name):
        """" Lese eine DNA-Sequenz aus einer FASTA-Datei und speichere sie intern ab. Aufgabe 2.
        Parameter: name: der Name der FASTA-Datei (voller Pfadname) """


    def read_from_ncbi(self, id, mail):
        """ Lese eine DNA-Sequenz von der NCBI-Nucleotid Datenbank. Aufgabe 2.
        Parameter: id: Accession Number der Sequenz
                   mail: die zu verwendende email-Adresse für den Entrez-Zugriff
        """    
        

    def transcribe(self):  
        """ Transkribiere die intere DNA-Sequenz in RNA. Dabei sei die DNA der codogene Strang.
        Die RNA-Sequenz soll intern gespeichert und zurück gegeben werden. Aufgabe 3
        Ergebnis: die transkribierte RNA-Sequenz (String)
        """

    def translate(self,initPos):    
        """ Translatiere die intern gespeicherte RNA in eine Aminosäurensequenz, ab der 
        vorgegebenen Position. Aufgabe 4
        Parameter: initPos: die Startposition, ab der die Translation beginnt
        Ergebnis: die Aminosäurensequenz (String)
        """


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
