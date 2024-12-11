## read environment variables
import os
from dotenv import load_dotenv
from dataclasses import dataclass


def fetch_sequence(accid, db="nucleotide", rettype="fasta", retmode="text"):
    """Fetch a sequence from NCBI using Entrez."""
    from Bio import Entrez

    Entrez.email = os.getenv("EMAIL")
    handle = Entrez.efetch(db=db, id=accid, rettype=rettype, retmode=retmode)
    return handle.read()


@dataclass
class Reference:

    accid: str
    filepath: str

    def __post_init__(self):
        self._fetch_sequence()

    def _fetch_sequence(self):
        if not os.path.exists(self.filepath):
            with open(self.filepath, "w", encoding="utf-8") as f:
                f.write(fetch_sequence(self.accid))

    def __str__(self):
        return f"Reference({self.accid})"

    @property
    def filename(self):
        """
        Return the filename of the reference sequence."""
        return os.path.basename(self.filepath)


class RefManager:
    """
    Manage reference sequences.
    Takes a list of accession IDs and a directory to store the sequences.
    """

    def __init__(self, accids, refdir):
        self.accids = accids
        self.refdir = refdir
        self.references = self._load_references()

    def _load_references(self):
        refs = []
        for accid in self.accids:
            filepath = os.path.join(self.refdir, f"{accid}.fasta")
            refs.append(Reference(accid, filepath))
        return refs


def samtools_flagstat_to_excel(flagstat, outfile):
    """
    Convert samtools flagstat output to an Excel file.
    """

    import pandas as pd

    lines = flagstat.split("\n")
    data = {}
    for line in lines:
        if line:
            key, value = line.split("(")
            data[key.strip()] = int(value.split()[0])

    df = pd.DataFrame(data, index=[0])
    df.to_excel(outfile, index=False)
