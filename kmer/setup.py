from distutils.core import setup, Extension


def main():
    fasta_module = Extension("fasta", ["read/fasta.c"])
    ffp_module = Extension("ffp", ["ffp/ffp.c", "ffp/hash.c"])
    setup(name="multi", version="0.0.1",
          description="Python interface for fasta C library function",
          author="Daniele Petillo",
          ext_modules=[fasta_module, ffp_module])


if __name__ == "__main__":
    main()