from distutils.core import setup, Extension


def main():
    fasta_module = Extension("fasta", ["fasta.c"])
    setup(name="fasta", version="0.0.1",
          description="Python interface for fasta C library function",
          author="Daniele Petillo",
          author_email="daniele.petillo@gmail.com",
          ext_modules=[fasta_module])


if __name__ == "__main__":
    main()