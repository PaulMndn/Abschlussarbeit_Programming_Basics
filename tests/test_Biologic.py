import GUI  # import GUI before Biologic due to circular import conflict
import Biologic

def test_foldRna():
    bio = Biologic.biologic()
    bio.rna = "UUUCAGUAGCA"
    res = bio.foldRna()
    expect = [
        [0,0,0,0,2,2,2,4,4,5,6],
        [0,0,0,0,0,2,2,2,4,4,5],
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
    assert res == expect

