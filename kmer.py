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
import kmerutils


class Kmer:

    def __init__(self, seq_dict=None, seq_dir=None):
        """A class to extract kmer words from genetic sequences
        for alignment-free algorithms.

        Parameters
        ----------
        seq_dict : dict
            Dictionary with sequences' info and the sequences
            themselves.
        seq_dir : str
            Relative path to sequences' folder.

        Attributes
        ----------
        ids : list
            Contains IDs of the genetic sequences.
        sequences: list of Seq objects
            Each element is a Seq object representing a genetic
            sequence.
        length_seqs : list
            Lengths of the genetic sequences measured in
            base pairs [bp].
        alphabet : str
            Genetic alphabet. Only "ATCG" supported.
        ordered_kmers : list of dicts 
            Each dictionary in the list is structured as:
            `key`: kmer word
            `value`: # of occurrences for `key`
            Also, each dictionary corresponds to a unique sequence.
        """

        if seq_dir:
            seq_dict = Biodata(seq_dir)
            seq_dict = seq_dict.load_as_dict()

        self.ids = seq_dict["IDs"]
        self.sequences = seq_dict["sequences"]
        self.length_seqs = [len(x) for x in self.sequences]

        
        self.alphabet = "ATCG"
        self.ordered_kmers = None
    


    def optimal_k(self, max_k=None):
        """Given a range of k values, the variety of the extracted
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
        """It extracts words from genetic sequences given the parameter `k`.
        If `k` is None, it is computed for each sequence using:

        k = log_4(sequence length)

        Finally `k` is averaged.

        Parameters
        ----------
        k : int
        A parameter used to define the length of a kmer word.
        """

        if k is None:
            logs = 0
            average_k = 0
            for n in self.length_seqs:
                logs += mt.log(n, 4)
            average_k = logs / len(self.length_seqs)
            print("Average k: %.2f" % average_k)
            k = int(average_k)
        else:
            k = int(k)

        kmers_words = kmerutils.generate_words(k, self.alphabet)

        print("Extracting words... ")
        self.ordered_kmers = [{} for x in range(len(self.sequences))]
        for index, sequence in enumerate(self.sequences):
            pos = 0
            end_pos = len(sequence) - k + 1
            kmers = []
            for pos in range(0, end_pos):
                sub = sequence[pos:k+pos]
                if not all(w in self.alphabet for w in sub):
                    continue
                else:
                    kmers.append(str(sub))
            unordered_dic_kmers = collections.Counter(kmers)
            for word in kmers_words:
                self.ordered_kmers[index][word] = unordered_dic_kmers[word]


        print("Words analysis completed.\n")


class subKmer(Kmer):

    def __init__(self, seq_dict=None, seq_dir=None, binning=100):
        """A class to cut two genetic sequences into subsequences
        for further operations.

        Parameters
        ----------
        seq_dict : dict
            Dictionary with sequences' info and the sequences
            themselves.
        seq_dir : str
            Relative path to sequences' folder.
        binning : int
            Defines the length of a subsequence in base pairs [bp].
        
        Attributes
        ----------
        binning : int
            Defines the length of a subsequence in base pairs [bp].
        """
        self.binning = binning
        super().__init__(seq_dict, seq_dir)

    
    def sKmer(self):
        """It cuts the sequences into subsequences based on
        `binning`.

        Returns
        -------
        cut : int
            Integer necessary to plot the sKmer heatmap.
        """

        subseqs = []
        cut = len(self.sequences[0]) // self.binning
        for ss in self.sequences:
            i = 0
            while i < len(ss) // self.binning:
                sub = ss[i*self.binning:(i+1)*self.binning]
                subseqs.append(sub)
                i += 1
        
        self.sequences = subseqs
        return cut



if __name__ == "__main__":

    quest = Kmer(seq_dir="test_seqs")
    quest.words_overlay()
    ans = Analysis(quest.ordered_kmers)
    corr_matrix = ans.correlation_matrix(len(quest.sequences), correlation=["P"])
    vis = Visualization(k=quest.k)
    vis.heatmap(corr_matrix, quest.ids)
    vis.histogram(quest.ordered_kmers)
