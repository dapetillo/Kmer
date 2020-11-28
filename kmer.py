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

    def __init__(self, seq_dict=None, corr=None):
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

        if seq_dict is not None:
            self.ids = seq_dict["IDs"]
            self.sequences = seq_dict["sequences"]
            self.length_seqs = [len(x) for x in self.sequences]
       
        if corr not in ["P", "S", "T", "ALL"]:
            self.corr = input("Correlation functions: \n\n-Pearson (P) \
             \n-Spearman (S) \n-Kendall (T) \n-All (ALL) \n\nChoose one of them: ")
        else:
            self.corr = corr
        
        self.alphabet = "ATCG"
        self.k = 0
        self.all_w = None
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

        self.all_w = np.empty(4**self.k, dtype=object)
        for index, items in enumerate(itertools.product(self.alphabet, repeat=self.k)):
            self.all_w[index] = ''.join(items)

        print("Extracting words... ")
        self.ordered_kmers = [[] for lists in range(len(self.sequences))]

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
            for key in self.all_w:
                self.ordered_kmers[index].append(unordered_dic_kmers[key])

        print("Words analysis completed.\n")


    def bootstrapping_BCa(self, alpha=0.04549, tolerance=10, B=10, BCa=True):
        """The method calculates confidence intervals for specific 
        correlation values of a 'model' sequence against N - 1 sequences,
        using bootstrapping and bootstrapping BCa. The parameters must be
        set according to the experiment (type of sequences, statistical 
        significance, computational power, etc.).
        The confidence levels are saved in *.txt format.

        WARNING: the bootstrapping theory makes use of samples generated 
        from the original one. However, between two sequences the corr. 
        value is only one: to have a sample of corr. values, 
        one need a sample of sequences, and the sequences being in 
        such a way that do not differentiate too much from the original ones. 
        This is accomplished sampling with replacement the words 
        (together with the corresponding occurrences) of the original 
        sequences.

        Parameters
        ----------
        alpha: 'float'
        It represents the significance level.
        tolerance: 'int'
        Because of the warning above, it happens that the sequences in
        S do not have the same length of their original ones. The tolerance
        variable forces the system to provide sequences with length L equal
        to L +/- tolerance. Smaller the tolerance, more reliable the 
        bootstrapping, but the computational time rises (exponentially).
        B: 'int'
        Number of bootstraps. Greater the value, stronger the statistics.
        BCa: 'boolean'
        A switch to either perform BCa or not.

        """
        print("Number of bootstraps: ", B)
        print("New sample size`s tolerance: +/- ", tolerance, "occurences")
        CL = 1- alpha
        print("Confidence level: ", CL*100, "%")

        theta_low, theta_up, corr_func = ([[] for l in range(0, 3)] for i in range(3))
        print(theta_low)
        z = scipy. stats.norm.ppf(alpha)
        zop = scipy.stats.norm.ppf(CL)
        x = 0
        if self.corr == "ALL":
            stop = 3 
            #corr_values = [[] for l in range(0,3)]
        else:
            stop = 1      #corr_values = []
            theta_low[1], theta_up[1], corr_func[1] = ([mt.nan]*(len(self.sequences) - 1) for l in range(3))
            theta_low[2], theta_up[2], corr_func[2] = ([mt.nan]*(len(self.sequences) - 1) for l in range(3))
        for x in range(0, 1):
            for y in range(1, len(self.sequences)):
                corr_values = [[] for l in range(0, stop)]
                n = 0
                for n in range(0, B):
                    N = 0
                    M = 0
                    new_size_N = 0
                    new_size_M = 0

                    
                    low_tol_N = sum(self.ordered_kmers[x]) - tolerance
                    up_tol_N = sum(self.ordered_kmers[x]) + tolerance
                    low_tol_M = sum(self.ordered_kmers[y]) - tolerance
                    up_tol_M = sum(self.ordered_kmers[y]) + tolerance
        
                    #the while loop forces the size within the tolerance level
                    while (new_size_N <= low_tol_N  or new_size_N >= up_tol_N) and (
                            new_size_M <= low_tol_M  or new_size_M >= up_tol_M):

                        sample_keys = np.random.choice(self.all_w, len(self.all_w))
                        values_N = []
                        values_M = []

                        for key in sample_keys:
                            index = np.where(self.all_w == key)
                            values_N.append(self.ordered_kmers[x][index[0][0]])
                            values_M.append(self.ordered_kmers[y][index[0][0]])

                        new_size_N = sum(values_N)
                        new_size_M = sum(values_M)

                    if self.corr == "S":
                         corr_values[0].append(scipy.stats.spearmanr(values_N, values_M)[0])
                    elif self.corr == "T":
                        corr_values[0].append(scipy.stats.kendalltau(values_N, values_M)[0])
                    elif self.corr == "P":
                        corr_values[0].append(scipy.stats.pearsonr(values_N, values_M)[0])
                    else:
                        corr_values[0].append(scipy.stats.spearmanr(values_N, values_M)[0])
                        corr_values[1].append(scipy.stats.kendalltau(values_N, values_M)[0])
                        corr_values[2].append(scipy.stats.pearsonr(values_N, values_M)[0])


                corr_values = [sorted(corr_values[h]) for h in range(0, stop)] #len(corr_values))]

                if BCa is True:

                    jack_theta = [[] for l in range(0, stop)]#len(corr_values))]
                    jack_theta_average = [[] for l in range(0, stop)]#len(corr_values))]
                    lower_alpha = []
                    upper_alpha = []
                    k = 0
                    for k in range(0, stop):#len(corr_values)):
                        count = 0
                        j = 0
                        for j in range(0, len(corr_values[k])): #It runs over the elements of a correlation function
                            jack_theta[k].append(sum(l for l in corr_values[k] if l != corr_values[k][j]))
                            if corr_values[k][j] < self.corr_matrix[k][x][y]:
                                count += 1
                                
                        jack_theta[k] = [l / (len(jack_theta[k]) - 1) for l in jack_theta[k]]        
                        jack_theta_average[k].append(sum(jack_theta[k]) / len(jack_theta[k]))
                        accel = (sum(jack_theta_average[k] - l for l in jack_theta[k])**3)/6*((
                            sum(jack_theta_average[k] - l for l in jack_theta[k])**2)**(3/2))
                        z_zero = scipy.stats.norm.ppf(count/B)

                        lower_alpha.append(int(scipy.stats.norm.cdf(z_zero + (z_zero + z)/
                                                                    (1 - accel*(z_zero + z)))*B))
                        upper_alpha.append(scipy.stats.norm.cdf(z_zero + (z_zero + zop)/
                                                                (1 - accel*(z_zero + zop)))*B)
                        if not upper_alpha[k][0].is_integer():
                            upper_alpha[k] = int(upper_alpha[k] + 1)
                        else: upper_alpha[k] = int(upper_alpha[k])

                        theta_low[k].append(corr_values[k][lower_alpha[k] - 1])
                        theta_up[k].append(corr_values[k][upper_alpha[k] - 1])
                        corr_func[k].append(self.corr_matrix[k][x][y])

                else:
                    lower_alpha = int(B*(1 - CL))
                    upper_alpha = CL*B
                    if not upper_alpha.is_integer():
                        upper_alpha = int(upper_alpha + 1)
                    else: 
                        upper_alpha = int(upper_alpha)
                    k = 0    
                    for k in range(0, stop):#len(corr_values)):
                        theta_low[k].append(corr_values[k][lower_alpha - 1])
                        theta_up[k].append(corr_values[k][upper_alpha - 1])
                        corr_func[k].append(self.corr_matrix[k][x][y])

    

        with open("Conf_int.txt", "wb+") as datafile_id:
            data = np.array([theta_low[0], theta_up[0], corr_func[0],
                             theta_low[1], theta_up[1], corr_func[1],
                             theta_low[2], theta_up[2], corr_func[2]])
            data = data.T
            np.savetxt(datafile_id, data, fmt="%f", delimiter="    ", header="SpearCIlow,\
                       SpearCIup, Spear, KenCIlow, KenCIup, Ken, PearCIL, PearCIup, Pear")




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
        for ind, ss in enumerate(self.sequences):
            if ind == 0:
                self.limit = len(ss) // binning
            i = 0
            while i < len(ss) // binning:
                sub = ss[i*binning:(i+1)*binning]
                subseqs.append(sub)
                i += 1
        self.sequences = subseqs





    def histogram(self):
        """It saves/shows the words distribution for a sequence.
        The words extraction must be performed before to call the
        method.

        """
        words = np.arange(4**self.k)
        occurr = self.ordered_kmers

        for ind in range(0, len(occurr)):

            occurr[ind] = [x / sum(occurr[ind]) for x in occurr[ind]]
            plt.clf()
            plt.bar(words, occurr[ind], align="center")
            plt.xticks(words, self.all_w, rotation="vertical")
            plt.title("Set title")
            plt.xlabel("Words")
            plt.ylabel("Frequencies")
            plt.savefig("Namefile{}.png".format(ind), bbox_inches="tight")
            #plt.close()
            #plt.show()



    def heatmap(self, matrix=np.array([]), title="Heatmap", fout="out.png"):
        """It visualizes the matrix correlation values via heatmap.
        Each row represents a sequence as well as each column.
        The labels present the ids' names.
        """
        
        if not matrix.any():
            matrix = self.corr_matrix

        if self.corr == "ALL":
            stop = len(self.corr_matrix)
        else:
            stop = 1

        for ind in range(0, stop):
            plt.clf()
            plt.figure()
            sns.heatmap(self.corr_matrix[ind], square=True, vmin=-1, vmax=1,
                        xticklabels=self.ids, yticklabels=self.ids, cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks":[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]}, fmt=".2f",
                        annot=False, annot_kws={"size": 9})
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.title(title)
            plt.savefig(fout, bbox_inches="tight")



    def heatmap_sKmer(self):
        """ It visualizes the matrix correlation values via heatmap when
        sKmer is applied. Each row represents the subsequences of a sequence Y
        while each column represents the subsequences of a sequence X.

        """
        if self.corr == "ALL":
            stop = len(self.corr_matrix)
            name_corr = ["Spearman", "Kendall", "Pearson"]
        else:
            stop = 1
            name_corr = [self.corr]
        for ind in range(0, stop):
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sns.heatmap(self.corr_matrix[ind][0:self.limit, self.limit:],
                        square=True, vmin=-1, vmax=1, cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks":[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]},
                        fmt=".2f", annot=False, xticklabels=10, yticklabels=10,
                        annot_kws={"size": 9})
            ax.plot([0, ax.get_ylim()[1]], [ax.get_ylim()[1], 0], ls="--", color=".3",
                    linewidth=1.)

            plt.title("sKmer - {} corr. for k = {} and bin = {} bp".format(
                      name_corr[ind], self.k, self.binning))
            plt.xlabel("Subsequences X [(x+1)*{} bp]".format(self.binning))
            plt.ylabel("Subsequences Y [(y+1)*{} bp]".format(self.binning))
            plt.savefig("Namefile{}.png".format(ind), bbox_inches="tight")
            #plt.pause(0.001)
            #plt.close()



if __name__ == "__main__":

    bdata = Biodata(seq_dir="test_seqs")
    bdata.load_as_dict()
    quest = Kmer(corr="P", seq_dict=bdata.biodata)
    quest.words_overlay()
    ans = Analysis(quest.ordered_kmers)
    corr = ans.correlation_matrix(len(quest.sequences))
    print(corr)
    sys.exit()
    quest.heatmap()


                
