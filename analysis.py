import scipy.stats
import numpy as np
import kmerutils

class Analysis:

    def __init__(self, kmers):
        """A class to perform statistical analyses on sequences
        manipulated with alignment-free algorithms.

        Parameters
        ----------
        kmers : list or list of dicts
            List of kmer words occurrences of list of dictionaries
            containing the occurrences.
        """

        self.kmers = kmerutils.unpack_kmers(kmers)


    def correlation_matrix(self, matrix_size, correlation=["P", "S", "T"]):
        """It correlates N sequences among each other using kmer words
        occurrences.

        Parameters
        ----------
        matrix_size : int
            Size of the correlation matrix.
        correlation : list of str
            List of correlations to calculate. Current correlations
            supported:
            - "P": Pearson
            - "S": Spearman
            - "T": Kendall's tau
        
        Returns
        -------
        corr_matrix : list of (2,) ndarray
            Each element in list is a correlation matrix among sequences.
        """

        corr_matrix = [np.zeros((matrix_size, matrix_size)) 
                       for l in range(0, len(correlation))]
        print("Calculating correlations...")
        for x in range(0, matrix_size):
            y = 0
            while x >= y:
                index_step = 0
                if "S" in correlation:
                    value = scipy.stats.spearmanr(self.kmers[x], self.kmers[y])[0]
                    corr_matrix[index_step][x][y] = value
                    index_step += 1
                if "T" in correlation:
                    value = scipy.stats.kendalltau(self.kmers[x], self.kmers[y])[0]
                    corr_matrix[index_step][x][y] = value
                    index_step += 1
                if "P" in correlation:
                    value = scipy.stats.pearsonr(self.kmers[x], self.kmers[y])[0]
                    corr_matrix[index_step][x][y] = value
                y += 1

        # fills the other half of the matrix (symmetry)
        for ind in range(0, len(correlation)):
            corr_matrix[ind] = (corr_matrix[ind] + corr_matrix[ind].T
                                - np.diag(corr_matrix[ind].diagonal()))

        return corr_matrix


