#!/usr/bin/env python3
import math as mt
import os
import itertools
import collections
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scipy.stats
import seaborn as sns
from Bio import SeqIO
import sys
import configparser
from data import Biodata
from analysis import Analysis
from visu import Visualization
import generic


class Kmer:

    """This class implements an alignment-free algorithm to correlate genetic sequences.

    The sequence A is cut into "words" of k nucleotides and for all the k reading frames;
    a vector then is filled with the occurrences of each possible words (4**k). The same 
    is made with the sequence B. Finally the two resulting vectors are correlated. If 
    a set of N sequences is given, the class generates N vector that are correlated 
    each other.

    The class includes a method to find bias-corrected and accelerated 
    confidence intervals for the correlation values.

    """

    def __init__(self, seq_dict=None, seq_dir=None):
        """It initializes the main attributes of the class.

        Attributes
        ----------
        seqs: 'list'
            A list containing the raw genetic sequences. Letters must be
            upper-case.
        length_seqs: 'list'
            A list containing the length of the genetic sequences in seqs.
            The length is measured in base pairs [bp].
        corr: 'str'
            Correlation function to use among the possible ones:

            - Pearson (P)
            - Spearman (S)
            - Kendall (T)
            - All (ALL)     <---- it correlates the sequences using the 
                                  three functions above. 

            It can be decided a priori. If None, the script will ask for
            one.
            WARNING: the script will either calculate with one correlation
            function or with all of them. The use of two correlations is 
            not possible (e.g. Pearson+Kendall).
        alphabet: 'str'
            The genetic alphabet. Only "ATCG" possible: reduced alphabet not
            supported.
        files: 'list'
            A list containing the sequences' names. If a sequence is
            downloaded from a database, the name usually corresponds to the
            file name.
        k: 'int'
            A parameter defining the length [bp] of an extracted "word".
        all_w: 'list'
            A list containing all the possible permutations given the 
            alphabet and the words' lengths. In this case, 4**k.
        corr_matrix: 'numpy 3-D array'
            A 3-D array that contains the correlation values. It resembles
            an array of matrices: each matrix refers to a correlation
            function.
        ordered_kmers: '2D list'
            A list containing a number of lists equal to the sequences' number.
            Each list contains the occurrences of all the possible words for a
            sequence. The occurrences of each list are sorted by all_w.

        """

        if seq_dir:
            seq_dict = Biodata(seq_dir)
            seq_dict = seq_dict.load_as_dict()

        self.ids = seq_dict["IDs"]
        self.sequences = seq_dict["sequences"]
        self.length_seqs = [len(x) for x in self.sequences]

        
        self.alphabet = "ATCG"
        self.k = 0
        self.corr_matrix = None
        self.ordered_kmers = None
    


    def optimal_k(self, max_k=None):
        """ Given a range of k values, the variety of the extracted
        words in a sequence changes. The method returns the (optimal)
        k(s) for which the variety (or richness) is maximum.
        The methodology is taken from:

        'Alignment-free genome comparison with feature frequency 
        profiles (FFP) and optimal resolutions', (Gregory E. Sims,
        Se-Ran Jun, Guohong A. Wu and Sung-Hou Kima).
        
        """
        min_k = 1
        if max_k is None:
            max_k = 8

        richness = np.zeros((len(self.sequences), max_k - 1))
        opt_k = {}
        k = 0
        for ind, seq in enumerate(self.sequences):
            for k in range(min_k, max_k):
                pos = 0
                end_pos = len(seq) - k + 1
                kmers = []
                for pos in range(0, end_pos):
                    sub = seq[pos:k+pos]
                    if not all(w in self.alphabet for w in sub):
                        continue
                    else:
                        kmers.append(str(sub))
                counting = collections.Counter(kmers)
                for values in counting.values():
                    if values >= 2:
                        richness[ind][k-1] += 1
            opt_k["{}".format(self.ids[ind])] = np.argmax(richness[ind]) + 1
        
        return opt_k

    def words_overlay(self, k=None):
        """The method extracts the words from each sequence given
        the parameter k. If k is None, the function will print
        the average k for the sequences based on the relation:

        k = log_4(sequence length),

        which theoretically finds the best k in the same
        fashion as in optimal_k (check the paper cited in the 
        latter for more informations). Then the user can choose
        the k to use.

        """
        if k is not None:
            self.k = k
        else:
            logs = 0
            average_k = 0
            for n in self.length_seqs:
                logs += mt.log(n, 4)
            average_k = logs / len(self.length_seqs)
            print("Average k: %.2f" % average_k)
            self.k = int(average_k)


        kmers_words = generic.generate_words(self.k, self.alphabet)

        print("Extracting words... ")
        self.ordered_kmers = [{} for x in range(len(self.sequences))]
        for index, sequence in enumerate(self.sequences):
            pos = 0
            end_pos = len(sequence) - self.k + 1
            kmers = []
            for pos in range(0, end_pos):
                sub = sequence[pos:self.k+pos]
                if not all(w in self.alphabet for w in sub):
                    continue
                else:
                    kmers.append(str(sub))
            unordered_dic_kmers = collections.Counter(kmers)
            for word in kmers_words:
                self.ordered_kmers[index][word] = unordered_dic_kmers[word]


        print("Words analysis completed.\n")



    def sKmer(self, binning=100):
        """This method cut sequences in subsequences: it is used when
        the user wants to look for local changes in a sequence. The
        method should be called before to perform any words extraction
        as well as correlation calculation. The subsequences get stored
        in the seqs attribute.

        Parameters
        ----------
        binning: 'int'
        It defines the length of the subsequences.
        
        """

        self.binning = binning
        subseqs = []
        cut = len(self.sequences[0]) // binning
        for ss in self.sequences:
            i = 0
            while i < len(ss) // binning:
                sub = ss[i*binning:(i+1)*binning]
                subseqs.append(sub)
                i += 1
        
        self.sequences = subseqs
        return cut



if __name__ == "__main__":

    quest = Kmer(seq_dir="test_seqs")
    #cut = quest.sKmer(binning=100)
    quest.words_overlay()
    ans = Analysis(quest.ordered_kmers)
    corr_matrix = ans.correlation_matrix(len(quest.sequences), correlation=["P"])
    vis = Visualization(k=quest.k)
    vis.heatmap(corr_matrix, quest.ids)
    #vis.histogram(quest.ordered_kmers)
    #vis.heatmap_sKmer(corr_matrix, cut)

