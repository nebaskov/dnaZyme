import re
import pandas as pd
import numpy as np

import pymatgen.core as mg
from PyBioMed import Pydna

from autoencoder import (
    encoding,
    generate_latent_representations,
    filter_sequences,
    generate_rdkit_descriptors
)
from constants import (
    similar_compounds,
    properties,
    kmer_2_dict
)


def flatten(lst):
    return [item for sublist in lst for item in sublist]


def extract_conditions(df_ml, buffer_column_name):
    dict_conditions = {}
    for index, row in df_ml.iterrows():
        i = row[buffer_column_name]
        if pd.isna(i):
            i = ''
        i = i.replace('and', ',')
        lst = i.split(',')
        buffer_composition = {}
        for component in lst:
            component = component.strip()
            if 'pH' in component:
                ph = re.search(r"\b\d[.]\d+", component).group(0)
                buffer_composition['pH'] = ph
            compound = ' '.join(re.findall(r"\b[a-zA-Z]{2,}\S+", component))
            if compound == '':
                continue
            try:
                concentration = re.search(r"^.*?(?=['M'])", component).group(0) + 'M'
                if len(concentration) == 1:
                    continue
                if compound in concentration:
                    concentration = re.search(r"^.*?(?=['M'])", component[len(compound) + 1:]).group(0) + 'M'
            except AttributeError:
                try:
                    concentration = re.search(r"\b\d*%", component).group(0)
                except:
                    continue
            if compound in flatten(similar_compounds.values()):
                compound = [key for key, val in similar_compounds.items() if compound in val][0]

            buffer_composition[compound] = concentration
        dict_conditions[index] = buffer_composition
    df_conditions = pd.DataFrame(dict_conditions).T
    return df_conditions


def convert_values(old_value):
    converter_dict = {'M': 1, 'mM': 0.001, 'µM': 0.000001, 'nM': 0.000000001, 'μM': 0.000001}
    if type(old_value) == str:
        if '%' in old_value:
            new_value = float(re.findall(r"[+]?\d*\.?\d+|\d+", old_value)[0])
            return new_value
        else:
            if 'M' not in old_value:
                return float(old_value)
        old_value = old_value.replace('mM', ' mM')
        number_units = old_value.split()
        new_value = float(number_units[0]) * converter_dict[number_units[1]]
        return new_value
    else:
        return old_value


def merge_cofactors(row):
    row = row.dropna()
    cof = [i for i in row.keys() if '+' in i]
    concentration = row[cof].sum()
    return concentration


def get_pymatgen_descriptors(element, charge):
    try:
        mg.Element(element)
    except ValueError:
        return {'element': None, 'charge': None, 'atomic_mass': None, 'ionization_energy': None,
                'electron_affinity': None, 'ionic_radii': None, 'Z': None}
    charge = float(charge)
    descriptors_set = {'element': element, 'charge': charge}
    ion = mg.Species(element, charge)
    for method in properties:
        ans = getattr(ion.element, method)
        try:
            result = ans()
        except(TypeError, AttributeError):
            result = ans
        if method == 'ionic_radii':
            if element == 'Ni':
                result = 0.83
            elif element == 'Cr':
                result = 0.62
            else:
                result = result[int(charge)]
        descriptors_set[method] = float(result)
    return descriptors_set


def calculate_ion_descriptors(metal_ions):
    lst = metal_ions.strip("[]").replace("'", "").split(",")
    if ('metal ion dependency not reported' in metal_ions) | ('M2+-independent' in metal_ions) | (lst == ['']) | (
            'Mg2+-independent' in metal_ions):
        return {'element': 0, 'charge': 0, 'atomic_mass': 0, 'ionization_energy': 0, 'electron_affinity': 0,
                'ionic_radii': 0, 'Z': 0}
    elif len(lst) == 1:
        ion = lst[0]
        element = ion[0:2]
        charge = ion[2]
        charge = charge.replace('+', '1')
        return get_pymatgen_descriptors(element, charge)
    else:
        ds = []
        for ions in lst:
            ions = ions.strip()
            element = ions[0:2]
            charge = ions[2]
            charge = charge.replace('+', '1')
            d_n = get_pymatgen_descriptors(element, charge)
            ds.append(d_n)
        d = {}
        for k in d_n.keys():
            full = list(d[k] for d in ds)
            if k == 'element':
                d[k] = full
            else:
                d[k] = sum(full) / len(full)
        return d


def calculate_kmer(seq, k):
    dnaclass = Pydna.PyDNA(seq)
    kmer = list(dnaclass.GetKmer(k=k).values())
    return kmer


def calculate_autoencoder(df_ml, seq_column_name):
    filtered_sequences = filter_sequences(
        sequences=df_ml,
        max_length=96,
        sequences_column_name=seq_column_name,
        shuffle_seqs=False
    )
    descriptors_set = generate_rdkit_descriptors()
    encoded_sequences = encoding(
        sequences_list=filtered_sequences[seq_column_name],
        polymer_type='DNA',
        descriptors=descriptors_set,
        num=96
    )
    x_autoencoder = generate_latent_representations(
        encoded_sequences=encoded_sequences,
        path_to_model_folder=r'autoencoder/nucleic_acids'
    )
    return x_autoencoder


def get_full_descriptors_set(df, seq_column_name, buffer_column_name, cofactor_column_name, temperature_column_name):
    kmer_2_features = list(kmer_2_dict.keys())
    autoencoder_features = list(generate_rdkit_descriptors().columns)
    features_names = [*kmer_2_features, *autoencoder_features]
    df_conditions = extract_conditions(df_ml=df, buffer_column_name=buffer_column_name)
    for i in df_conditions.columns:
        df_conditions.loc[:, i] = df_conditions.loc[:, i].apply(convert_values)
    df_conditions['cofactor concentration'] = df_conditions.apply(lambda row: merge_cofactors(row), axis=1)
    df_conditions['temperature'] = df[temperature_column_name]
    df_conditions = df_conditions.dropna(thresh=30, axis=1).drop(['Mg2+', 'Zn2+', 'Mn2+', 'EDTA'], axis=1).fillna(0)
    df_cofactor = pd.DataFrame()
    df_cofactor[
        [
            'element',
            'charge',
            'atomic_mass',
            'ionization_energy',
            'electron_affinity',
            'ionic_radii',
            'Z'
        ]
    ] = df[cofactor_column_name].apply(
        calculate_ion_descriptors, axis=1, result_type='expand'
    ).fillna(0)
    df_sequence = pd.DataFrame(np.concatenate(
        (
            np.array([calculate_kmer(seq, 2) for seq in df[seq_column_name]]),
            calculate_autoencoder(df_ml=df, seq_column_name=seq_column_name)
        ), axis=1
    ))
    df_sequence.columns = [*features_names]
    df_full = pd.concat([df_sequence,
                         df_conditions,
                         df_cofactor], axis=1)
    return df_full
