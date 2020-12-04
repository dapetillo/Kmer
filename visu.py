import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import kmerutils

class Visualization:

    def __init__(self, k, save=True):
        """A class to visualize results from alignment-free
        algorithms.

        Parameters
        ----------
        k : int
            A parameter used to define the length of a kmer word.
        save : bool
            Whether save a result to file.
        
        Attributes
        ----------
        k : int
            A parameter used to define the length of a kmer word.
        """
        self.k = k


    def histogram(self, kmers, alphabet="ATCG"):
        """Builds the histogram of kmer words for a sequence.

        Parameters
        ----------
        kmers : list of dicts
            Each dictionary in the list is structured as:
            `key`: kmer word
            `value`: # of occurrences for `key`
            Also, each dictionary corresponds to a unique sequence.
        alphabet : str
            Genetic alphabet. Only "ATCG" supported.
        """

        kmers = kmerutils.unpack_kmers(kmers)
        kmers_words = kmerutils.generate_words(self.k, alphabet)
        
        number_of_words = np.arange(4**self.k)
        for ind in range(0, len(kmers)):

            kmers[ind] = [x / sum(kmers[ind]) for x in kmers[ind]]
            plt.clf()
            plt.bar(number_of_words, kmers[ind], align="center")
            plt.xticks(number_of_words, kmers_words, rotation="vertical")
            plt.title("Set title")
            plt.xlabel("Words")
            plt.ylabel("Frequencies")
            plt.savefig("Namefile{}.png".format(ind), bbox_inches="tight")


    def heatmap(self, corr_matrix, ids, title="Heatmap", fout="out.png"):
        """Builds the heatmap of a correlation matrix.
        
        Parameters
        ----------
        corr_matrix : list of (2,) ndarray
            Each element in list is a correlation matrix among sequences.
        ids : list
            List of IDs' sequences.
        title : str
            Heatmap's title.
        fout : str
            Output namefile.
        """
    
        for ind in range(0, len(corr_matrix)):
            mask = np.zeros_like(corr_matrix[0], dtype=np.bool)
            mask[np.triu_indices_from(mask)] = True
            plt.clf()
            plt.figure()
            sns.heatmap(corr_matrix[ind], mask=mask, square=True, vmin=-1, vmax=1,
                        xticklabels=ids, yticklabels=ids,
                        cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks": np.arange(-1, 1.25, 0.25)},
                        fmt=".2f", annot=False, annot_kws={"size": 9})
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.title(title)
            plt.savefig(fout, bbox_inches="tight")


    def heatmap_sKmer(self, corr_matrix, cut, binning=100):
        """ Builds the heatmap of a correlation matrix.
        Each row represents the subsequences of a sequence Y
        while each column represents the subsequences of a sequence X.

        Parameters
        ----------
        corr_matrix : list of (2,) ndarray
            Each element in list is a correlation matrix
            between subsequences.
        cut : int
            Limit to the correlation matrix to plot.
        binning : int
            Binning used to cut a sequence.
        """
    
        for ind in range(0, len(corr_matrix)):
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sns.heatmap(corr_matrix[ind][0:cut, cut:],
                        square=True, vmin=-1, vmax=1, cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks": np.arange(-1, 1.25, 0.25)},
                        fmt=".2f", annot=False, xticklabels=10, yticklabels=10,
                        annot_kws={"size": 9})
            ax.plot([0, ax.get_ylim()[0]], [0, ax.get_ylim()[0]], ls="--", color=".3",
                     linewidth=1.)

            plt.title("sKmer - {} corr. for k = {} and bin = {} bp".format(
                      "Pearson", self.k, binning))
            plt.xlabel("Subsequences X [(x+1)*{} bp]".format(binning))
            plt.ylabel("Subsequences Y [(y+1)*{} bp]".format(binning))
            plt.savefig("Namefile{}.png".format(ind), bbox_inches="tight")
            #plt.pause(0.001)
            #plt.close()
