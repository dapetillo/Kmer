import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Visualization:

    def __init__(self, k, save=True, fout="out.png"):
        self.k = k


    def heatmap(self, corr_matrix, ids, title="Heatmap", fout="out.png"):
        """It visualizes the matrix correlation values via heatmap.
        Each row represents a sequence as well as each column.
        The labels present the ids' names.
        """

        for ind in range(0, len(corr_matrix)):
            plt.clf()
            plt.figure()
            sns.heatmap(corr_matrix[ind], square=True, vmin=-1, vmax=1,
                        xticklabels=ids, yticklabels=ids,
                        cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks": [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]},
                        fmt=".2f", annot=False, annot_kws={"size": 9})
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.title(title)
            plt.savefig(fout, bbox_inches="tight")


    def heatmap_sKmer(self, corr_matrix, cut, binning=100):
        """ It visualizes the matrix correlation values via heatmap when
        sKmer is applied. Each row represents the subsequences of a sequence Y
        while each column represents the subsequences of a sequence X.

        """
    
        for ind in range(0, len(corr_matrix)):
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sns.heatmap(corr_matrix[ind][0:cut, cut:],
                        square=True, vmin=-1, vmax=1, cmap="RdBu_r", linewidths=.1,
                        cbar_kws={"ticks":[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]},
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
