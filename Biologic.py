import logging
from Bio import SeqIO, Entrez
import pathlib
import re
import sys

import GUI

log = logging.getLogger(__name__)

def get_bonds():
    """Returns a dictionary with number of bonds of binding rna base pairs."""
    bonds = {
        "AU": 2,
        "GC": 3,
        "GU": 1
    }
    # add reversed pairs to dict
    keys = list(bonds.keys())
    for i in keys:
        bonds["".join(reversed(i))] = bonds[i]
    
    return bonds



class biologic:
    ## constants

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

    # translation
    START_CODON = "Met"
    STOP_CODON = "_"


    def __init__(self):
        """Create an instance of biologic. No parameters are required as they 
        are set later through methods.
        
        Raises:
            FileNotFoundError: genCode.txt not found at ./data/genCode.txt
        """
        # placeholders
        self.dna = None
        self.rna = None
        self.proteins = None
        
        # create translation table
        table_fp = GUI.SCRIPT_DIR/"data"/"genCode.txt"
        if not table_fp.exists():
            # genCode.txt not existent
            log.critical(f'File for codon table not found at {table_fp}')
            raise FileNotFoundError(f'No codon table file fount at {table_fp}')
        
        codon_table = {}
        with open(table_fp, "r") as ct_file:
            for line in ct_file.readlines():
                # line consists of amino acid and space seperated codons
                # stop codons code are represented in the peptide by an underscore
                line.strip()
                split = line.split()
                aa = split[0]
                for codon in split[1:]:
                    codon_table[codon] = aa
        self.codon_table = codon_table
        log.info(f'Codon table created from {table_fp}')
        log.debug(f'Created codon table: {codon_table}')

        # set protein alphabet
        self.ALPHABETS['protein'] = list(codon_table.values())
        log.debug(f'Protein alphabet added to `biologic.ALPHABETS`')


    def set_dna(self,dna):
        """ Setze die interne DNA-Sequenz. Aufgabe 1
        
        Arg:
            dna: DNA sequence as string
        """
        self.dna = dna.upper()
        log.info(f"New DNA entered: {dna}")

    def validate_dna(self):
        """ Überprüfe die intern gespeicherte DNA. Aufgabe 1

        Checks if all letters of DNA sequence are valid according to 
        the alphabet.

        Returns:
            bool: True if sequence is valid, False otherwise
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
        Reads DNA sequence from FASTA file and saves it.
        Sequence may be an empty string if no record was found in file.

        Arg:
            fp (str): file path of fasta file
        
        Returns:
            str: DNA sequence from fasta file
        """
        with open(fp, "r") as file:
            # take first record from fasta file
            record = next(SeqIO.parse(file, "fasta"), None)

            if record is None:
                # no record was read
                log.warning(f'No record was read from fasta file at location {fp}')
                self.set_dna("")
            else:
                self.set_dna(str(record.seq).upper())
                log.info(f'DNA sequence was successfully read from fasta file: {self.dna}')
            
            return self.dna


    def read_from_ncbi(self, id: str, mail: str):
        """ Lese eine DNA-Sequenz von der NCBI-Nucleotid Datenbank. Aufgabe 2.
        Read DNA sequence of specific accession number from NCBI nucleotide DB 
        and save it.
        
        Args:
            id: Accession Number der Sequenz
            mail: die zu verwendende email-Adresse für den Entrez-Zugriff
        
        Returns:
            str: DNA sequence from NCBI
        """
        Entrez.email = mail
        log.debug(f'{mail} set for Entrenz call.')
        
        # get db handle for accession number in fasta format
        handle = Entrez.efetch(
            db='nucleotide', 
            id=id, 
            rettype='fasta', 
            retmode='text'
        )
        # get SeqRecord from handle
        record = SeqIO.read(handle, 'fasta')
        handle.close()

        log.debug(f'Fasta file for id {id} fetched from NCBI. {record}')

        # save sequence
        self.set_dna(str(record.seq).upper())
        log.info(f'DNA sequence for {id} fetched from NCBI.')
        return self.dna


    def transcribe(self):  
        """ Transkribiere die intere DNA-Sequenz in RNA. Dabei sei die DNA der codogene Strang.
        Die RNA-Sequenz soll intern gespeichert und zurück gegeben werden. Aufgabe 3
        
        Returns:
            None if no DNA sequence is present
            str: rna sequence

        """
        if self.dna is None:
            # No DNA sequence present
            log.error('Could not transcribe DNA to RNA because DNA does not exists')
            return
        
        self.rna = "".join(self.DNA_TO_RNA[i] for i in self.dna)
        log.info(f'DNA transcribed to RNA: {self.rna}')
        return self.rna


    def translate(self,initPos: int = 0, reverte_rna: bool =False):    
        """ Translatiere die intern gespeicherte RNA in eine Aminosäurensequenz, ab der 
        vorgegebenen Position. Aufgabe 4

        Translate RNA sequence to peptide sequence. 
        If an unknown codon is incountered it is added into the peptide sequence
        as ``?(<codon>)``.

        Argument:
            initPos: Starting position for translation. Defaults to 0.
            reverte_rna: Whether RNA sequence is reversed for tranlsation.
                Defaults to False.

        Returns:
            None if no RNA sequence is present
            str: peptid sequence
        """
        if self.rna is None:
            # No RNA sequence present
            log.error("Could not translate DNA because RNA does not exist")
            return

        # reverte RNA if `reverte_rna` is True
        rna = self.rna if not reverte_rna else "".join(i for i in reversed(self.rna))

        peptid = []
        log.debug(f'Start translation.')

        # iterate over codons starting at `initPos`
        for i in range(initPos, len(rna), 3):
            codon = rna[i:i+3]
            if len(codon) < 3:
                # No complete codon left
                continue

            try: 
                # add corresponding amino acid to peptid sequence
                peptid.append(self.codon_table[codon])
            except KeyError:
                # unknown codon, should not be possible with tested codon table.
                log.warning(f'Codon {codon} not found in codon table.')
                peptid.append(f'?({codon})')
        
        protein = "".join(peptid)
        log.info(f'Protein sequence created. {protein}')
        return protein


    def get_proteins(self, rf_number) -> list[str]:
        """ Erzeuge für einen vorgegebenen Leserahmen alle daraus ablesbaren 
        Proteine. Ein Protein beginnt mit 'M' und endet mit einem Stop '_'. 
        Die Leserahmen haben die Nummern 0 bis 5. Dabei sind 0-2 die 
        Rahmen in der vorgegebenen Richtung und 3-5 die Rahmen der revertierten 
        Sequenz. Aufgabe 5. 
        
        Arg: 
            rf_number: Nummer des Leserahmens (reading frame).
        
        Returns: 
            None if rf_number is invalide.
            List of proteins.
        """
        # map rf_number to values for `reverte_rna` and `initPos` to pass into 
        # `self.translate()`
        rf_mapping = {
            0: (False, 0),
            1: (False, 1),
            2: (False, 2),
            3: (True, 0),
            4: (True, 1),
            5: (True, 2)
        }

        try:
            reversed, initPos = rf_mapping[rf_number]
        except KeyError:
            # Invalide rf_number, should not be possible as it is checked in GUI as well.
            log.error(f'Invalid reading frame: {rf_number}. Expected valule is 0-5.')
            return
        
        log.info(f"Trnaslating {'forward' if not reversed else 'reverse'} rna. "
            + f"Starting at position {initPos}.")
        
        # get full peptide
        translation = self.translate(initPos, reversed)
        if translation is None:
            # No RNA sequence was present
            log.error(f"RNA was not tranlsated to peptid. Can't match for proteins.")
            return
        
        log.debug(f"RNA tranlsated to peptid: {translation}")
        
        # find proteins in peptide
        prot_pattern = re.compile(r"Met[A-Za-z?\(\)]*_")
        self.proteins = re.findall(prot_pattern, translation)
        log.info(f"{len(self.proteins)} proteins found: {self.proteins}")
        
        return self.proteins


    def eval_rna(self, struct):
        """ Berechne für eine vorgegebene Struktur die Anzahl der enthaltenen Wasserstoffbrücken.
        Aufgabe 6.

        Arg: 
            struct: eine vorgegebene Sekundärstruktur in Klammerschreibweise.
        
        Returns:
            None: No RNA sequence is present
                or invalid base pairing in structure.
            int: Total number of hydrogen bonds in folded RNA.
        """
        if self.rna is None:
            # No RNA sequence is present
            log.error("No RNA sequence present to calculate bonds for.")
            return
        
        log.info(f"Start bond count calculation")
        log.debug(f"Bond count calculation with structure: {struct} on sequence: {self.rna}.")
        bond_count = 0

        # keep track of opening bracket indices
        opening_indices = []
        for i in range(len(struct)):
            char = struct[i]
            if char == ".":
                # no bond
                continue

            elif char == "(":
                # start base pair
                opening_indices.append(i)
                continue

            elif char == ")":
                # end base pair
                opening_index = opening_indices.pop()
                pair = self.rna[opening_index] + self.rna[i]
                try:
                    # add number of hydrogen bonds for base pair to total count
                    bond_count += get_bonds()[pair]
                except KeyError:
                    # Invalid base pairing
                    log.error(f"Invalid bond {pair} in struct {struct} for RNA sequence {self.rna}.")
                    return
        
        log.debug(f"Finished Calculation. Bond count: {bond_count}")
        return bond_count

    def foldRna(self) :
        """ Berechne die optimale Sekundärstruktur nach dem vorgegebenen Verfahren. Aufgabe 7.
        
        Calculates secundary structure using the provided algorythm. 
        The bond matrix is saved as an object variable.

        Returns: 
            int: Optimal number of hydrogen bonds.
        """
        log.debug("Start calculate bond-count-matrix")
        # init zero-matrix
        self.bond_mat = [[0 for i in self.rna] for i in self.rna]
        log.debug("Finshed initializing zero matrix")

        # iterate over bond_mat diagonally
        for i in range(4,len(self.rna)):
            for j in range(0,len(self.rna)-i):
                y,x = j, i+j
                # calc bond count from subsequence from pos y to x and save to bond_mat

                # get first and last bases of sequence and their bonds
                bases = self.rna[y] + self.rna[x]
                try:
                    bases_bonds = get_bonds()[bases]
                except KeyError:
                    # No base pair formation possible
                    bases_bonds = 0

                # number of bonds of first and last base + inner sequence (left down)
                first_last_bases_bonds = bases_bonds + self.bond_mat[y+1][x-1]
                
                if i == 4:
                    # if first diagonal just save to matrix and continue
                    # no sub-sequences possible (primary priximity contraint)
                    self.bond_mat[y][x] = first_last_bases_bonds
                    continue
                
                # calculate max bond counts from partial sub-sequences
                partial_sequences_bonds = max(
                    self.bond_mat[y][k] + self.bond_mat[k+1][x] \
                        for k in range(y,x)
                )

                # save max of first_and_last_bases_bind and 
                # partial_sub-sequences in bond_mat
                self.bond_mat[y][x] = max(
                    first_last_bases_bonds, 
                    partial_sequences_bonds
                )
                
                ## print-outs for debug
                # print(self.rna)
                # print(bases)
                # print("\n".join(str(i) for i in self.bond_mat))
                # print()
        
        log.debug(f"Finished creating bond matrix.")
        log.debug(str(self.bond_mat))

        log.info("Matrix of bonds of (sub)sequences calculated.")
        log.info(f"Max number of bonds for RNA sequence {self.rna} "
            + f"calculated to be {self.bond_mat[0][-1]}")
        
        # return max number of hydrogen bonds
        return self.bond_mat[0][-1]



    def struct_string(self,x,y):
        '''Recursive helper function to trace back through the bond count matrix
        and find the folding structure for subsequence from position y to x.

        Args:
            x: x-position in matrix, end of subsequence
            y: y-position in matrix, start of subsequence
        
        Returns:
            str: Folding structure of substring.
        
        Raises:
            Exception: if tracing algorithm failed
        '''
        bonds = self.bond_mat[y][x]
        # bond count for first and last base
        try:
            base_bonds = get_bonds()[self.rna[x]+self.rna[y]]
        except KeyError:
            base_bonds = None

        if bonds == 0:
            # no bonds in substring
            return "." * (abs(y-x)+1)
        
        if bonds == self.bond_mat[y+1][x]:
            # first base does not bind
            return "." + self.struct_string(x,y+1)
        
        if base_bonds is not None and bonds == (self.bond_mat[y+1][x-1] + base_bonds):
            # first and last bases bind
            return f"({self.struct_string(x-1, y+1)})"
        
        # none of the above
        # split into 2 subsequences that add to the bond count of current subsequence
        for k in range(y, x):
            if bonds == (
                self.bond_mat[y][k]
                + self.bond_mat[k+1][x]
            ):
                return self.struct_string(k,y) + self.struct_string(x, k+1)
        
        # tracing algorythm failed (fail safe)
        log.error("None of the possible tracing options for finding a "
            + "binding structure yielded a result.")
        
        raise Exception("None of the possible tracing options yielded a result.")
   


    def trace(self) :
        """ Berechne für die zuletzt durchgeführte Faltung die optimale Struktur. Aufgabe 8.
        Ergebnis: optimale Struktur in Klammerschreibweise.
        """
        log.debug(f"Start calculating structural string for RNA {self.rna}.")
        
        # start tracing bonds in bonding matrix starting with full RNA sequence
        binding_structure = self.struct_string(len(self.rna)-1, 0)

        return binding_structure


if __name__ == '__main__':
    #################################################
    ## ONLY FOR DEV PURPOSES, NO EXPLICIT FUNCTION ##
    #################################################

    # teststring ArgMetCysAlaLys_AsnGluMetPhe_Met
    # pattern "Met[A-Za-z?\(\)]*_"

    # bio = biologic()
    # bio.rna = "UCGACUCGGAG"
    # res = bio.foldRna()
    # print("\n".join(str(i) for i in bio.bond_mat))


    pass