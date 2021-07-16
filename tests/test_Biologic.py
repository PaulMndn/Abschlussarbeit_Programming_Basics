import pytest

import GUI  # import GUI before Biologic due to circular import conflict
import Biologic



def test_set_dna():
    bio = Biologic.biologic()
    bio.set_dna("AGTTCTAGGCGCCATACAAGATA")
    assert bio.dna == "AGTTCTAGGCGCCATACAAGATA"


def test_validate_dna():
    bio = Biologic.biologic()
    bio.set_dna("AGTTCTAGGCGCCATAC")
    assert bio.validate_dna() == True

    bio.set_dna("AGTTCTAGGCGNCCATAC")
    assert bio.validate_dna() == False

    bio.set_dna("ACTCCTUGTAGGTCA")
    assert bio.validate_dna() == False


def test_get_rna():
    bio = Biologic.biologic()
    bio.set_dna("AGTTCTAGGCGCCATACAAGATA")
    bio.translate()
    rna = bio.get_rna()
    # internally saved RNA is equal to returned RNA
    assert rna == bio.rna


def test_read_fasta_file():
    bio = Biologic.biologic()
    bio.read_fasta_file(str(GUI.SCRIPT_DIR/"data"/"blue_pigment.fasta"))

    assert bio.dna == "ATGGTCAACAACCGTAACCATGGGCTGGACTTACGGCTTGTCACCATTCCTTCATTCTTCTCCAAGAGTGCTTGCATCTACAATCCCATCATCTACTGCTTCATGAATGTAAAGCTCT"


def test_read_from_ncbi():
    bio = Biologic.biologic()
    bio.read_from_ncbi("M26172.1", "paul.menden@student.hswt.de")

    assert bio.dna == "ATGGTCAACAACCGTAACCATGGGCTGGACTTACGGCTTGTCACCATTCCTTCATTCTTCTCCAAGAGTGCTTGCATCTACAATCCCATCATCTACTGCTTCATGAATGTAAAGCTCT"


def test_transcription():
    bio = Biologic.biologic()
    bio.set_dna("AGTTCTAGGCGCCATACAAGATA")
    transcript = bio.transcribe()
    assert transcript == "UCAAGAUCCGCGGUAUGUUCUAU"


def test_tralslation():
    bio = Biologic.biologic()
    bio.rna = "AUGCGAAUAAAUUAAGUAG"
    prot = bio.translate(initPos=0)
    assert prot == "MetArgIleAsn_Val"

    prot = bio.translate(3)
    assert prot == "ArgIleAsn_Val"

    prot = bio.translate(15)
    assert prot == "Val"


def test_get_proteins():
    bio = Biologic.biologic()
    bio.rna = "AUGCGAAUAAAUUAAGUAG"
    proteins = bio.get_proteins(0)
    assert proteins == ["MetArgIleAsn_"]

    proteins = bio.get_proteins(4)
    assert proteins == ["MetAsn_"]

    bio.rna = "UAUGCGAUAAACGAUGUAAA"
    proteins = bio.get_proteins(1)
    assert proteins == ["MetArg_", "Met_"]


def test_eval_rna():
    bio = Biologic.biologic()
    bio.rna = "UUUCAGUAGCA"
    bonds = bio.eval_rna("..(..(...))")
    assert bonds == 5


def test_eval_rna_2():
    bio = Biologic.biologic()
    bio.rna = "UCGACUCGGAG"
    bonds = bio.eval_rna("(((...).)).")
    assert bonds == 8


def test_foldRna():
    bio = Biologic.biologic()
    bio.rna = "UUUCAGUAGCA"
    res_bonds = bio.foldRna()
    expect_mat = [
        [0,0,0,0,2,2,2,3,3,5,5],
        [0,0,0,0,0,1,1,2,3,3,5],
        [0,0,0,0,0,0,0,2,3,3,5],
        [0,0,0,0,0,0,0,0,3,3,3],
        [0,0,0,0,0,0,0,0,0,3,3],
        [0,0,0,0,0,0,0,0,0,3,3],
        [0,0,0,0,0,0,0,0,0,0,2],
        [0]*11,
        [0]*11,
        [0]*11,
        [0]*11,
    ]
    assert bio.bond_mat == expect_mat
    assert res_bonds == expect_mat[0][-1]

    assert bio.trace() == "..(..(...))"


def test_foldRna_2():
    bio = Biologic.biologic()
    bio.rna = "UCGACUCGGAG"
    res_bonds = bio.foldRna()
    expect_mat = [
        [0,0,0,0,0,0,3,6,7,8,8],
        [0,0,0,0,0,0,3,6,6,6,6],
        [0,0,0,0,0,0,3,3,3,3,5],
        [0,0,0,0,0,0,0,0,3,3,5],
        [0,0,0,0,0,0,0,0,3,3,5],
        [0,0,0,0,0,0,0,0,0,2,3],
        [0,0,0,0,0,0,0,0,0,0,3],
        [0]*11,
        [0]*11,
        [0]*11,
        [0]*11,
    ]
    print(i for i in bio.bond_mat)
    assert bio.bond_mat == expect_mat
    assert res_bonds == expect_mat[0][-1]

    assert bio.trace() == "(((...).))."

