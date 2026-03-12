#!/usr/bin/env python3
"""Reader/writer for GFFv3 files as defined by
http://www.sequenceontology.org/gff3.shtml
"""

import re
from collections import defaultdict
from urllib.parse import unquote
from .Locus import Locus

class Gff3record:
    seqid_str = "[a-zA-Z0-9.:^*$@!+_?\-|]+"
    seqid_re = re.compile("^"+seqid_str+"$")
    field_str = "[^\t\n\r;=%&,]*"
    field_re = re.compile("^"+field_str+"$")
    float_str = (
        # sign
        "[+-]?(?:"+
        # mantissa >= 1
        "(?:[\d]+\.?[\d]*)|"+
        # mantissa < 1
        "(?:\.[\d]+))"+
        # exponent
        "(?:[eE][+-]?[\d]+)?")
    gff3_re = re.compile("^"+"(?P<seqid>"+seqid_str+")\t"+
                         "(?P<source>"+field_str+")\t"+
                         "(?P<_type>"+field_str+")\t"+
                         "(?P<start>[\d]+)\t"+
                         "(?P<end>[\d]+)\t"+
                         "(?P<score>(?:\.|(?:"+float_str+")))\t"+
                         "(?P<strand>[+-.?])\t"+
                         "(?P<phase>[012.])"+
                         "(?:\t(?P<attributes>[^\t]+))?$")
    def __init__(self,
                 seqid,
                 source,
                 _type,
                 start,
                 end,
                 score,
                 strand,
                 phase,
                 attributes):
        # Landmark for local coordinate system
        self.seqid = seqid
        # Feature generator (algorithm or method).
        # In practice, a qualifier to the type
        self.source = source
        # (Previously, "method").  Data type from a Sequence Ontology
        self.type = _type
        # First feature position relative to landmark, counting from 1
        self.start = int(start)
        # Last feature position relative to landmark, counting from 1
        #   end >= start.  (For circular landmarks, end wraps past the
        #   last position of the landmark; length zero features imply
        #   an interstitial feature on the positive side of start)
        self.end = int(end)
        # Optional floating point score
        try:
            self.score = float(score)
        except:
            self.score = "."
        # Strand: +,-,. (unstranded), ? (unknown)
        self.strand = str(strand)
        # Phase: 0,1,2 (for CDS features only)
        try:
            self.phase = int(phase)
        except:
            self.phase = "."
        # dict of one to many key=[value,...] lists
        self.attributes = attributes
        
    @classmethod
    def fromString(cls, s):
        """Parse a GFF3 record based on the standard defined at:
        http://www.sequenceontology.org/gff3.shtml."""
        p = cls.gff3_re.search(s)
        try:
            assert(p is not None)
        except AssertionError:
            print(s)
            raise
        args = p.groupdict()
        if(args["attributes"] is not None):
            t = defaultdict(lambda : [])
            for i in p.group("attributes").split(";"):
                (key,val) = [j.strip() for j in i.split("=")]
                t[key].append(unquote(val))
            args["attributes"] = t
        else:
            args["attributes"] = {}

        return cls(**args)

    def isValid(self):
        """Return True if the state of this record conforms to the GFF3
        standard, False otherwise.

        N.B. The current implementation of this function is incomplete.
        N.B. A current copy of the Sequence Ontology must be available
             for validation of the type field."""

        return all((
            Gff3record.seqid_re.search(self.seqid) is not None,
            Gff3record.field_re.search(self.source) is not None,
            Gff3record.field_re.search(self.type) is not None,
            self.end >= self.start,
            (self.score == ".") or (type(self.score) == float),
            self.strand in "+-.?",
            # TODO: Expand _type check for SO version(s) of CDS features
            #       (Actually, this may be okay as is)
            ((self.phase == ".") and (self.type != "CDS"))
            or (type(self.phase) == int),
            isinstance(self.attributes, dict),
            all(isinstance(i, str) for i in self.attributes.keys()),
            all(isinstance(i, list) for i in self.attributes.values()),
            # TODO: Enforce additional assertions for reserved keys
            ))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "\t".join(map(str,(
            self.seqid,
            self.source,
            self.type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            self.attribute_string())))

    def attribute_string(self):
        if(self.attributes is None):
            return ""
        retval = []
        for (key, vals) in self.attributes.items():
            for val in vals:
                retval.append("%s=%s" % (key,val))

        return ";".join(retval)

    # For compatibility with GffFile.GffRecord
    rest = property(attribute_string)

    def Locus(self, genome = None):
        return Locus(ref = self.seqid,
                     start = self.start,
                     stop = self.end,
                     strand = self.strand,
                     genome = genome)

    def Write(self, fp):
        fp.write("%s\n" % self)

class Gff3file:
    """Collection of GFF3 records."""
    def __init__(self, features = None):
        if(features is not None):
            self.features = features
        else:
            self.features = []

        # Generate name, type, and child indices
        #   (taking care not to insert new keys into the attributes
        #    defaultdicts)

        self.name_to_feature = {}
        self.type_to_features = defaultdict(lambda : [])
        for feature in self.features:
            self.type_to_features[feature.type].append(feature)
            if(feature.attributes is not None):
                if("ID" in feature.attributes):
                    for name in feature.attributes["ID"]:
                        self.name_to_feature[name] = feature
                elif("Name" in feature.attributes):
                    for name in feature.attributes["Name"]:
                        self.name_to_feature[name] = feature

        self.children = defaultdict(lambda : [])
        for child in self.features:
            if((child.attributes is not None) and
               ("Parent" in child.attributes)):
                for parent in child.attributes["Parent"]:
                    self.children[self.name_to_feature[parent]].append(child)

    @classmethod
    def fromFile(cls, fp):
        """Create a Gff3file from an open filestream in GFF3 format."""
        features = []
        for i in fp:
            if(i[0] != "#"):
                if(len(i.strip()) < 1):
                    # Case 0: blank line
                    continue
                else:
                    # Case 1: feature
                    features.append(Gff3record.fromString(i.rstrip("\r\n")))
            elif(i[1] == "#"):
                # Case 2: directive
                
                # Currently recognized directives from the GFF3 spec:
                # ##FASTA
                #      This notation indicates that the annotation portion
                #      of the file is at an end and that the remainder of
                #      the file contains one or more sequences (nucleotide
                #      or protein) in FASTA format. This allows features
                #      and sequences to be bundled together.

                if(i.startswith("##FASTA")):
                    break

            else:
                # Case 3: comment
                continue

        return cls(features)

    @classmethod
    def fromLoci(cls, loci,
                 source = "source",
                 _type = "type",
                 score = ".",
                 phase = ".",
                 attributes = None):
        """Generate a Gff3file from an iterable of loci."""

        if(type(source) == str):
            source = lambda i, locus, s = source: s
        if(type(_type) == str):
            _type = lambda i, locus, s = _type: s
        if(type(score) == str):
            score = lambda i, locus, s = score: s
        if(type(phase) == str):
            phase = lambda i, locus, s = phase: s
        if(attributes is None):
            attributes = lambda i, locus : {"Name":("Feature%05d" % i,)}
        elif(not hasattr(attributes, "__call__")):
            attributes = lambda i, locus, s = attributes: s

        return cls([Gff3record(seqid = locus.Locus().ref,
                               source = source(i, locus),
                               _type = _type(i, locus),
                               start = locus.Locus().start,
                               end = locus.Locus().stop,
                               score = score(i, locus),
                               strand = locus.Locus().strand,
                               phase = phase(i, locus),
                               attributes = attributes(i, locus))
                    for (i, locus) in enumerate(loci)])

    def __len__(self):
        return len(self.features)

    def __getitem__(self, i):
        return self.features[i]

    def Write(self, fp, header = True):
        if(header):
            fp.write("##gff-version 3\n")
        for i in self.features:
            i.Write(fp)

    def attribute_dot(self):
        """Return attribute dependency DAG in GraphViz dot format."""
        typemap = defaultdict(lambda : set())

        for (parent, children) in self.children.items():
            for child in children:
                typemap[parent.type].add(child.type)
        
        retval = """digraph G\n{\nsize="8,10";\ncenter=true;\nrankdir=LR;\n"""
        for (parent, children) in typemap.items():
            for child in children:
                retval += '"%s"->"%s";\n' % (parent,child)
        retval += """}\n"""

        return retval

if(__name__ == "__main__"):
    print("Hello, world")
    # TODO: add read/write test
