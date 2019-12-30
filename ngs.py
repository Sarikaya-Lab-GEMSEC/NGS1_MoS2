import os
import pickle

import pandas as pd
import numpy as np


def ngs_data(dataset='train', consensus=2, count_floor=8, y_type='class',
             directory='../data/ngs/',
             files=['Set1.csv', 'Set2.csv', 'Set3.csv'],
             output_file=None):
    # WARNING: consensus currently has no effect, all data pulls will use consensus 2
    # Check input values
    if not (dataset == 'train') | (dataset == 'test'):
        raise ValueError('dataset must equal train or test')
    if not (consensus == 2) | (consensus == 3):
        raise ValueError('consensus must be 2 or 3')
    if count_floor <= 0:
        raise ValueError('count_floor must be greater than 0')
    if not (y_type == 'class') | (y_type == 'counts') | (y_type == 'affinity')\
         | (y_type == 'center'):
        raise ValueError('y_type must be class, counts, or affinity')

    # Walk through sequences; populate a new dataframe with parsed data for
    # only sequences meeting the requirements
    query_file = os.path.join(directory, 
                  str(consensus) + str(count_floor) + dataset + '_set.pickle')

    if not os.path.isfile(query_file):  # Check if the import is precalculated
        print('Query does not exist: creating ' + query_file)
        # Get the full dataset:
        sequence_file = os.path.join(directory, 'sequence_dict.pickle')
        if not os.path.isfile(sequence_file):
            print('Dictionary does not exist: creating ' + sequence_file)
            sequences = {}
            for i, f in enumerate(files):
                df = pd.read_csv(os.path.join(directory, f))
                print('opened file ', f)
                for index, row in df.iterrows():
                    seq = row['AA_seq']
                    if not seq in sequences:
                        sequences[seq] = np.zeros((3, 5))
                    sequences[seq][i] = [row['CP1'], row['CP2'], row['CP3'], row['CE'],
                                      row['binding_affinity']]
            del df
            with open(sequence_file, 'wb') as f:
                pickle.dump(sequences, f, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(sequence_file, 'rb') as f:
                sequences = pickle.load(f)
            print('Loaded sequence dictionary')
    
        # Perform the test/train split using 10% of combined sequences for the test
        # set.
        seqs = np.array(list(sequences.keys()))
        seed = 42
        np.random.seed(seed)
        test_set_fraction = 0.1
        test_idx = np.random.random(size=len(seqs)) <= test_set_fraction
        train_idx = ~test_idx

        if dataset == 'train':
            print('Parsing training set')
            seqs = seqs[train_idx]
        elif dataset == 'test':
            print('Parsing test set')
            seqs = seqs[test_idx]

        # Walk through sequences; populate a new dataframe with parsed data
        df = []
        for i, seq in enumerate(seqs):  # Parse each sequence in the train or test set
            # Report progress?
            if not i % np.round(len(seqs)/100):
                print('Parsing sequences: ' +\
                      str(np.round(100 * i / len(seqs))) + '%')
            
            data = sequences[seq]  # Get the counts for sequence 'seq'

            # Enforce count_floor
            data = data[np.sum(data[:, :-1], axis=1) > count_floor]
            data = data[data[:, -1].argsort()]  # sort by affinity values

            # check consensus: binding_affinity must be within +-1
            rows = data[1:, -1] - data[:-1, -1]
            rows = np.less(rows, 1)
            if any(rows):  # if at least two rows satisfy consensus
                rows = np.insert(rows, 1, True)  # Make dims match
                entry = {}
                entry.update({'AA_seq': seq})
                entry.update({'binding_affinity': np.mean(data[rows, -1])})

                # Calculate count probability distribution
                P = np.sum(data[rows, :-1], axis=0, keepdims=False)
                P = P/P.sum()
                entry.update({'CP1': P[0], 'CP2': P[1], 'CP3': P[2], 'CE': P[3]})
                entry.update({'center_of_mass': (P * [0, 1, 2, 3]).sum()})
                
                # Assign classification label
                if entry['binding_affinity'] < 5:
                    entry.update({'class': 'weak'})
                elif entry['binding_affinity'] > 15:
                    entry.update({'class': 'strong'})
                else:
                    entry.update({'class': 'medium'})
                df.append(entry)  # add entry to the final dictionary

        # convert the dictionary to a dataframe and save it
        df = pd.DataFrame(df)
        with open(query_file, 'wb') as f:
            pickle.dump(df, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print('Loaded query from ' + query_file)
        with open(query_file, 'rb') as f:
            df = pickle.load(f)

    # Create excel file if requested
    if output_file:
        df.to_csv(output_file)

    # Return input/output lists
    X = df['AA_seq']
    if y_type == 'class':
        Y = df['class']
    elif y_type == 'counts':
        Y = df[['CP1', 'CP2', 'CP3', 'CE']]
    elif y_type == 'affinity':
        Y = df['binding_affinity']
    elif y_type == 'center':
        Y = df['center_of_mass']
    return X.to_numpy(), Y.to_numpy()


def sanger_data(data_file):
    """data file must be a csv with individual columns for amino acids
    in columns 'C:N' and response (y) values in column 'O'"""
    df = pd.read_csv(data_file)
    x = df[df.columns[2:2+12]]
    x = np.array([''.join(row) for row in x.itertuples(index=False)])
    y = df[df.columns[14]].to_numpy()
    return x, y