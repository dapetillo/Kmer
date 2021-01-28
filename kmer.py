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
    
    def upper_k_res(self, k_init=3, k_max=6, mean=True):
        k_upper = np.array([])
        for index in range(len(self.sequences)):
            k = k_init
            kl_div = []
            CRE = []
            init_l_kmer = True
            l_kmer_chain = []
            while k <= k_max:
                exp_l_kmer = {}
                if init_l_kmer:
                    l_kmer_chain.append(ffp.ffp(str(self.sequences[index]), k))
                    l_kmer_chain.append(ffp.ffp(str(self.sequences[index]), k - 1))
                    l_kmer_chain.append(ffp.ffp(str(self.sequences[index]), k - 2))
                else:
                    l_kmer_chain.insert(0, ffp.ffp(str(self.sequences[index]), k))
                    l_kmer_chain.pop(-1)
                N = sum(l_kmer_chain[2].values()) / sum(l_kmer_chain[1].values()) / sum(l_kmer_chain[1].values())
                for word in l_kmer_chain[0].keys():
                    exp_freq = l_kmer_chain[1][word[1:]] * l_kmer_chain[1][word[:-1]] / l_kmer_chain[2][word[1:-1]] * N
                    if not np.isfinite(exp_freq) or exp_freq == 0:
                        del l_kmer_chain[0][word]
                    else:
                        exp_l_kmer[word] = exp_freq
                exp_list = np.array(list(exp_l_kmer.values()))
                l_kmer_list = np.array(list(l_kmer_chain[0].values())) / sum(l_kmer_chain[0].values())
                kl_div.append(np.sum(l_kmer_list * np.log2(l_kmer_list / exp_list)))
                k += 1
                init_l_kmer = False
            for i in range(len(kl_div)):
                CRE.append(np.sum(kl_div[i:]))
            np.append(k_upper, CRE.index(min(CRE)) + k_init)
        
        if mean is True and k_upper.shape[0] > 1:
            k_upper = int(np.mean(k_upper))

        return k_upper

        


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
    quest.upper_k_res()
    sys.exit()
    quest.FFP()
    ans = Analysis(quest.ordered_counted_kmers)
    corr_matrix = ans.correlation_matrix(len(quest.sequences), correlation=["P"])
    vis = Visualization(k=quest.k)
    vis.heatmap(corr_matrix, quest.ids)
    vis.histogram(quest.ordered_counted_kmers)
