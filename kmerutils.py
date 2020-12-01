import itertools


def generate_words(k, alphabet="ATCG"):
    kmers_words = []
    for items in itertools.product(alphabet, repeat=k):
        kmers_words.append(''.join(items))

    return kmers_words


def unpack_kmers(kmers):
    if all(isinstance(d, dict) for d in kmers):
        kmers = [list(d_kmers.values()) for d_kmers in kmers]

    return kmers
