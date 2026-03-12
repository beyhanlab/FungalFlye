#!/usr/bin/env python3
"""Tools for parsing a FASTA file or stream to a dict-like collection
of strings.  This module is intended for lightweight string-oriented
manipulation of sequences.  Use Sequence.py for more"biological"
representations of protein and nucleotide sequences.
"""

import re

comp = {"a":"t","t":"a","c":"g","g":"c","n":"n",
        "A":"T","T":"A","C":"G","G":"C","N":"N"}
def complement(sequence):
    """Return the reverse complement of sequence, stripping
    any non-nucleotide (atgcn) characters."""
    retval = ""
    for i in range(len(sequence) - 1, -1, -1):
        try:
            retval += comp[sequence[i]]
        except KeyError:
            pass
    return retval

class FastaFile:
    """A set of named sequences with optional annotations, read
    from a FASTA formatted text file or stream.  Sequences are
    Python strings rather than higher-level sequence objects."""
    def __init__(self, file):
        """Initialize from filename or file stream in FASTA format.
        """

        # Make a compiled regular expression to parse the file
        # (Note: This is __way__ faster than parsing line by line)

        fasta_re = re.compile(
            # Header lines are marked by the '>' sign
            # We allow headers to begin in the middle of a line
            # We remove whitespace at the ends of the header
            ">[\s]*(?P<header>.*?)[\s]*$"+
            # Sequence is anything that is not a header line
            # We currently count any text in comments (';.*$') as
            # sequence (could parse this out at the same time as
            # whitespace)
            "(?P<seq>[^>]*)",
            # Use multiline mode to parse an entire FASTA file in one go
            re.M)
        name_re = re.compile(
            "^(?P<name>[\S]+)")

        seq_re = re.compile("[^A-Za-z]+")

        # Parse file
        if(isinstance(file, str)):
            text = file
        else:
            text = file.read()
            if(isinstance(text, bytes)):
                text = text.decode()
        parsed = fasta_re.findall(text)
        
        # Strip garbage and index sequences/annotations by name
        self.seqs = {}
        self.headers = {}
        self.key_order = []
        for i in parsed:
            k = name_re.search(i[0]).group("name")
            self.key_order.append(k)
            if(k in self.seqs):
                print("WARNING: overwriting", k)
            # Note: could use a user-supplied re for garbage removal
            self.seqs[k] = seq_re.sub("",i[1])
            self.headers[k] = i[0]

    def __repr__(self):
        return "<FASTA file object "+self.seqs.__repr__()+">"
    def __str__(self):
        return "FASTA file with %d sequences" % (len(self.seqs))
    def __len__(self):
        return len(self.seqs)
    def __getitem__(self, i):
        return self.seqs[i]
    def __iter__(self):
        for key in self.key_order:
            yield (key, self.seqs[key])

# Input/output stream versions.
# TODO: fold these into the original FastaFile class
#       to avoid redundant implementation.

def fasta_stream(fp):
    header = None
    seq = None
    for line in fp:
        if(line.startswith(">")):
            if(seq is not None):
                yield header, seq
            header = line[1:].rstrip("\r\n")
            seq = ""
        else:
            assert(header is not None)
            seq += re.sub("[\s]+","",line)
    if(header is not None):
        assert(seq is not None)
        yield header, seq

class FastaFormatter:
    def __init__(self, out):
        self.out = out
    def write(self, stanza):
        (header, seq) = stanza
        self.out.write(">%s\n%s\n" % (
            header,
            re.sub(r"(.{80})", r"\1\n", seq)))

if(__name__ == "__main__"):
    print("Hello, world")
