import GUI  # import GUI before Biologic due to circular import conflict
import Biologic

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

