import random
import numpy as np
import pandas as pd

from rdkit.Chem import rdMolDescriptors
from rdkit import Chem

from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
from keras.models import Model


aa_dict = {'dA': r'O=P(O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)C[C@@H]3O',
           'dT': r'CC1=CN(C(=O)NC1=O)C2CC(C(O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O',
           'dG': r'O=P(O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n1cnc2c1NC(=N/C2=O)\N)C[C@@H]3O',
           'dC': r'C1[C@@H]([C@H](O[C@H]1N2C=CC(=NC2=O)N)CO[P@@](=O)(O)O[P@@](=O)(O)OP(=O)(O)O)O'}


def generate_rdkit_descriptors(normalize: tuple = (-1, 1), monomer_dict_=aa_dict):
    descriptor_names = list(rdMolDescriptors.Properties.GetAvailableProperties())
    get_descriptors = rdMolDescriptors.Properties(descriptor_names)
    num_descriptors = len(descriptor_names)

    descriptors_set = np.empty((0, num_descriptors), float)

    for _, value in monomer_dict_.items():
        molecule = Chem.MolFromSmiles(value)
        descriptors = np.array(get_descriptors.ComputeProperties(molecule)).reshape((-1, num_descriptors))
        descriptors_set = np.append(descriptors_set, descriptors, axis=0)

    sc = MinMaxScaler(feature_range=normalize)
    scaled_array = sc.fit_transform(descriptors_set)
    return pd.DataFrame(scaled_array, columns=descriptor_names, index=list(monomer_dict_.keys()))


def filter_sequences(sequences: [list, pd.DataFrame],
                     max_length: int = 96,
                     sequences_column_name: str = None,
                     shuffle_seqs: bool = True,
                     aa_dict_=aa_dict):
    if type(sequences) is list:
        all_seqs = [seq.upper() for seq in sequences if len(seq) <= max_length]
        all_seqs = list(dict.fromkeys(all_seqs))
        filtered_seqs = [x for x in all_seqs if set(x).issubset(set(aa_dict_.keys()))]
        if shuffle_seqs:
            filtered_seqs = random.sample(filtered_seqs, len(filtered_seqs))
        return filtered_seqs

    elif type(sequences) is pd.DataFrame:
        sequences[sequences_column_name] = sequences[sequences_column_name].map(lambda x: x.replace(" ", ""))
        sequences[sequences_column_name] = sequences[sequences_column_name].str.upper()
        sequences = sequences[sequences[sequences_column_name].apply(lambda x: len(x) <= max_length)]
        peptide_subs = sequences[sequences[sequences_column_name].apply(lambda x:
                                                                        set(x).issubset(set(aa_dict_.keys())))]
        if shuffle_seqs:
            peptide_subs = peptide_subs.sample(frac=1)
        return peptide_subs


def seq_to_matrix_(sequence, polymer_type, descriptors, num):
    if polymer_type == 'DNA':
        prefix = 'd'
    else:
        print('Wrong polymer type')
        return

    rows = descriptors.shape[1]
    seq_matrix = tf.zeros(shape=[0, rows])
    for aa in sequence:
        aa_params = tf.constant(descriptors.loc[prefix+aa],
                                dtype=tf.float32)
        descriptors_array = tf.expand_dims(aa_params,
                                           axis=0)
        seq_matrix = tf.concat([seq_matrix, descriptors_array],
                               axis=0)
    seq_matrix = tf.transpose(seq_matrix)
    shape = seq_matrix.get_shape().as_list()[1]
    if shape < num:
        paddings = tf.constant([[0, 0], [0, num - shape]])
        add_matrix = tf.pad(
            seq_matrix,
            paddings=paddings,
            mode='CONSTANT',
            constant_values=-1
        )

        return add_matrix

    return seq_matrix


def encoding(sequences_list, polymer_type, descriptors, num):
    container = []
    for i, sequence in enumerate(sequences_list):
        if i % 3200 == 0:
            print(i * 100 / len(sequences_list), ' %')

        seq_matrix = tf.expand_dims(
            seq_to_matrix_(
                sequence=sequence,
                polymer_type=polymer_type,
                descriptors=descriptors,
                num=num
            ), axis=0
        )
        container.append(seq_matrix)

    encoded_seqs = tf.concat(
        container,
        axis=0
    )

    return encoded_seqs


def generate_latent_representations(
    sequant_encoded_sequences,
    path_to_model_folder=''
):
    trained_model = tf.keras.models.load_model(path_to_model_folder)

    layer_name = 'Latent'
    intermediate_layer_model = Model(inputs=trained_model.input,
                                     outputs=trained_model.get_layer(layer_name).output)
    latent_representation = intermediate_layer_model.predict(sequant_encoded_sequences)
    return latent_representation
