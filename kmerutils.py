import itertools


def generate_words(k, alphabet="ATCG"):
    """Computes all possible combinations (with repetition)
    of kmer words from a given genetic alphabet.

    Parameters
    ----------
    k : int
        Length of the kmer word.
    alphabet : str
        Genetic alphabet.

    Returns
    -------
    kmers_words : list of str
        List of all possible combinations of kmer words.
    """

    kmers_words = []
    for items in itertools.product(alphabet, repeat=k):
        kmers_words.append(''.join(items))

    return kmers_words


def unpack_kmers(kmers):
    """It unpacks kmers words from complex data structure.

    Parameters
    ----------
    kmers : list or list of dicts
        List of kmer words occurrences of list of dictionaries
        containing the occurrences.

    Returns
    -------
    kmers : list of lists
        Each nested list contains occurrences of kmer words
        for a unique sequence.
    """

    if all(isinstance(d, dict) for d in kmers):
        kmers = [list(d_kmers.values()) for d_kmers in kmers]

    return kmers
