from Bio import SeqIO
import os
import configparser


class Biodata:

    def __init__(self, seq_dir=None, gb_features=["source", "CDS"]):
        """A class to read sequences from different file formats and
        corresponding metadata.

        Parameters
        ----------
        seq_dir : str
            Relative path to folder's sequences
        gb_features : list of str
            List of features to read from GenBank files. A list of
            supported features can be found in **features.ini**
        
        Attributes
        ----------
        gb_features : list of str
            List of features to read from GenBank files. A list of
            supported features can be found in **features.ini**
        full_fname : list of str
            Absolute path to files' sequences
        biodata : dict
            It stores sequences and metadata.
        """

        feat_parser = configparser.ConfigParser()
        feat_parser.read("features.ini")

        for element in gb_features:
            value = feat_parser["genbank"].get(element)
            if value is None:
                raise ValueError(
                    "{} is an invalid feature for GenBank".format(element))
        
        self.gb_features = gb_features
        
        if seq_dir is not None:
            self.full_fname = [os.path.abspath(os.path.join(seq_dir, x)) for x in os.listdir(seq_dir)]
        
        self.biodata = {"IDs": [], "features": [], "sequences": []}
    

    def parse_genbank(self, gb):
        """It parses GenBank files based on features.

        Parameters
        ----------
        gb : str
            Absolute path to Genbank file.
        

        Returns
        -------
        ids : list of str
            IDs of the sequences.
        feats : list of str
            Features' names of the sequences.
        seqs : list of Seq objects
            The sequences themselves.
        """
    
        feat_parser = configparser.ConfigParser()
        feat_parser.read("features.ini")
        ids = []
        feats = []
        seqs = []

        for record in SeqIO.parse(gb, "genbank"):
            for feat in record.features:
                if feat.type in self.gb_features:
                    sequence = feat.location.extract(record).seq
                else:
                    continue

                ids.append(record.id)
                feats.append(feat.type)
                seqs.append(sequence)

        return (ids, feats, seqs)

    def parse_fasta(self, fa, acc_version="NC_"):
        """It parses GenBank files based on features.

        Parameters
        ----------
        fa : str
            Absolute path to FASTA file.
        acc_version : str
            Identifier to parse ID from FASTA file.

        Returns
        -------
        ids : list of str
            IDs of the sequences.
        feats : list of str
            Features' names of the sequences.
        seqs : list of Seq objects
            The sequences themselves.
        """
        ids = []
        feats = []
        seqs = []
        for record in SeqIO.parse(fa, "fasta"):
            seqs.append(record.seq.upper())
            split_id = record.id.split('|')
            acc_version_id = [x for x in split_id if x.startswith(acc_version)]
            if acc_version_id:
                    ids.append(acc_version_id[0])
            elif record.id.startswith("gi"):
                ids.append(split_id[1])  # extract GI number
            feats.append("source")
        
        return (ids, feats, seqs)

    def load_as_dict(self):
        """Calls parsing functions and collects output
        into dictionary.

        Returns
        -------
        biodata : dict
            It stores sequences and metadata.

        Raises
        ------
        ValueError
            If file format not supported.
        """

        for name in self.full_fname:
            if name.endswith((".gb", ".genbank")):
                ids, feats, seqs = self.parse_genbank(name)
            elif name.endswith((".fa", ".fasta")):
                ids, feats, seqs = self.parse_fasta(name)
            else:
                raise ValueError("File format not supported.")
            
            self.biodata["IDs"] += ids
            self.biodata["features"] += feats
            self.biodata["sequences"] += seqs
        
        return self.biodata


if __name__ == "__main__":
    bdata = Biodata(seq_dir="test_seqs")
    bdata = bdata.load_as_dict()
    print(bdata)
