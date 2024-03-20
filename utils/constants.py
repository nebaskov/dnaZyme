similar_compounds = {
    'Tris⋅HCl': ['Tris', 'Tris-HCl', 'Tris–HCl'],
    'HEPES': ['Na.Hepes', 'Hepes'],
    'Mn2+': ['MnCl2'],
    'NaCl': ['NaCI', 'NaNO3', 'NaCl'],
    'Mg2+': ['MgCl2'],
    'Zn2+': ['ZnCl2', 'Zn2+'],
    'Pb2+': ['PbOAc', 'Pb2+'],
    'Cd2+': ['CdCl2'],
    'Ni2+': ['NiCl2'],
    'alcohol': ['methanol'],
    'Ag+': ['AgNO3', 'Ag+'],
    'sodium cacodylate': ['cacodylate'],
    'Co2+': ['CoCl2'],
    'Cu2+': ['CuCl2']
}

properties = [
    'atomic_mass',
    'ionization_energy',
    'electron_affinity',
    'ionic_radii',
    'Z'
]

kmer_2_dict = {
    'AA': 0, 'AC': 0, 'AG': 0, 'AT': 1,
    'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
    'GA': 0, 'GC': 1, 'GG': 0, 'GT': 0,
    'TA': 0, 'TC': 0, 'TG': 1, 'TT': 0
}

sel_features = [
    'AC', 'AT', 'CC', 'CG', 'TA', 'TC', 'TT',
    'exactmw', 'amw', 'lipinskiHBD', 'NumRotatableBonds',
    'NumAtoms', 'FractionCSP3', 'CrippenMR', 'chi0v', 'kappa3',
    'NaCl', 'pH', 'KCl', 'cofactor concentration', 'temperature',
    'charge', 'electron_affinity', 'ionic_radii'
]
