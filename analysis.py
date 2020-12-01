import scipy.stats
import numpy as np
import generic

class Analysis:

    def __init__(self, kmers):
        self.kmers = generic.unpack_kmers(kmers)


    def correlation_matrix(self, matrix_size=None, correlation=["P", "S", "T"]):
        """It correlates N sequences among each other using the words
        occurrences. Given the symmetric nature of the corr. functions, only
        N((N-1)/2 + 1) values are calculated.

        """
        if matrix_size is None:
            raise TypeError("Matrix size not specified!")
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


