#!/usr/bin/env python3
import math as mt
import os
import itertools
import collections
import numpy as np
from Bio import SeqIO
import sys
from data import Biodata
from analysis import Analysis
from visu import Visualization
import kmerutils as ku
from scipy import special

import time
import matplotlib.pyplot as plt

import fasta
import ffp

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

        start = time.time()
        if seq_dir:
            seq_dict = Biodata(seq_dir)
            seq_dict = seq_dict.load_as_dict()
        print("Time to load fasta genome is %s seconds" % (time.time() - start))
        self.ids = seq_dict["IDs"]
        self.alphabets = seq_dict["alphabets"]
        self.sequences = seq_dict["sequences"]
        self.length_seqs = [len(x) for x in self.sequences]

        ###to cancel after test
        
        self.ordered_counted_kmers = None
    

    def lower_k_res(self, mean=True):
        low_k = np.log(self.length_seqs) / np.log(4)
        if mean is True and low_k.shape[0] > 1:
            low_k = int(np.mean(low_k))
        
        return low_k
    
    def upper_k_res(self):
        k = 3
        #keys = ku.generate_words(k, "ATCG")
        #exp_freq = {word: None for word in keys}
        CRE = 1
        for index in range(len(self.sequences)):
            while CRE > 0.1:
                exp_l_kmer = {}
                keys = ku.generate_words(k, self.alphabets[index])
                start = time.time()
                l_kmer = self.FFP(k, norm=True)
                print("Time to FFP for k={} is {}".format(k, time.time() - start))
                l1_kmer = self.FFP(k - 1, norm=True)
                print("Time to FFP for k={} is {}".format(k - 1, time.time() - start))
                l2_kmer = self.FFP(k - 2, norm=True)
                print("Time to FFP for k={} is {}".format(k - 2, time.time() - start))
                sys.exit()
                for word in keys:
                    exp_freq = l1_kmer[index][word[1:]] * l1_kmer[index][word[:-1]] / l2_kmer[index][word[1:-1]]
                    if not np.isfinite(exp_freq) or np.isnan(exp_freq) or exp_freq == 0:
                        del l_kmer[index][word]
                    else:
                        exp_l_kmer[word] = l1_kmer[index][word[1:]] * l1_kmer[index][word[:-1]] / l2_kmer[index][word[1:-1]]
                a = np.array(list(exp_l_kmer.values()))
                b = np.array(list(l_kmer[index].values()))
                if k == 4:
                    pass
                    #print("Siamo a k 4",special.kl_div(a, b))
                    #print(list(l_kmer[index].keys()))
                kl_div = np.sum(special.kl_div(b, a))
                print(kl_div, "al giro", k)
                CRE += kl_div
                k += 1
        


    def FFP(self, k=None, norm=False):
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
            k = self.lower_k_res()
            print("k is: ", k)
        else:
            k = int(k)

        #print("Extracting words... ")
        self.ordered_counted_kmers = [{} for x in range(len(self.sequences))]
        for index, sequence in enumerate(self.sequences):
            end_pos = len(sequence) - k + 1
            kmers = [str(sequence[pos:k+pos]) for pos in range(end_pos)]
            unordered_counted_kmers = collections.Counter(kmers)
            kmers_words = ku.generate_words(k, self.alphabets[index])
            for word in kmers_words:
                self.ordered_counted_kmers[index][word] = unordered_counted_kmers[word]
            if norm is True:
                norm_factor = 1.0 / sum(self.ordered_counted_kmers[index].values())
                for key in self.ordered_counted_kmers[index].keys():
                    self.ordered_counted_kmers[index][key] = self.ordered_counted_kmers[index][key] * norm_factor
        
        return self.ordered_counted_kmers
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
    quest = Kmer(seq_dir="test_genome")
    seconds = []
    for k in range(1, 10):
        print("Starting k=%s..." % k)
        start = time.time()
        #quest.FFP(k)
        ffp.ffp(str(quest.sequences[0]), k)
        seconds.append(time.time() - start)
    print(seconds)
    plt.scatter(np.arange(1, 10), seconds)
    plt.title("FFP time profiling - NC_049958.1")
    plt.xlabel("Word length [bp]")
    plt.ylabel("Execution time [s]")
    plt.savefig("FFP_time_C.png", bbox_inches="tight")
    sys.exit()
    quest.upper_k_res()
    quest.FFP()
    ans = Analysis(quest.ordered_counted_kmers)
    corr_matrix = ans.correlation_matrix(len(quest.sequences), correlation=["P"])
    vis = Visualization(k=quest.k)
    vis.heatmap(corr_matrix, quest.ids)
    vis.histogram(quest.ordered_counted_kmers)
