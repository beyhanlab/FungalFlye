#!/usr/bin/env python3
"""
Tools for manipulating expression profiling data, based on the
JavaTreeView CDT file format.
"""

from MsvUtil import median, prototypeable, transpose
from SafeMath import safelog, safesub, safemean
from TempWrapper import TempDir

import csv
import os.path
import shutil
import sys, csv

try:
    from numpy import array, arange
except ImportError:
    sys.stderr.write("WARNING: numpy not available\n")

try:
    # Prefer local build of Michael de Hoon's Cluster3 python bindings
    import Pycluster
except ImportError:
    try:
        # Fall back on Biopython's packaging of these bindings
        import Bio.Cluster as Pycluster
    except:
        sys.stderr.write("WARNING: Pycluster not available\n")

# Function dropped from Bio.Cluster
from cluster_methods import _savetree
        
def toFloat(x):
    if(x == None):
        return 0.0
    assert(type(x) == float)
    return x

def toMask(x):
    if(x == None):
        return 0
    return 1

@prototypeable
class CdtRow:
    def __init__(self, gid, uniqid, name = None, gweight = 1.0, ratios = [],
                 extra = [], leaf_offset = None, **kw):
        self.gid = gid
        self.uniqid = uniqid
        self.name = name
        self.gweight = gweight
        self.ratios = [None]*len(ratios)
        for (i, r) in enumerate(ratios):
            try:
                self.ratios[i] = float(r)
            except:
                pass
            
        self.extra = extra
        if(leaf_offset is not None):
            self.length = float(self.extra[leaf_offset])
        else:
            self.length = None
        
    def Gid(self):
        return self.gid
    
    def Uniqid(self):        
        return self.uniqid
    
    def Name(self):
        return self.name

    def Gweight(self):
        return self.gweight

    def Ratios(self):
        return self.ratios

    def Extra(self):
        return self.extra

    def Completeness(self):
        if(len(self.ratios) == 0):
            return None
        return sum(1 for i in self.ratios if(i is not None))/float(len(self))

    def __len__(self):
        return len(self.ratios)

    def __getitem__(self, i):
        return self.ratios[i]

    # hack for GtrNode
    def Depth(self):
        return 0

    #def __str__(self):
    #    return str(self.Uniqid())

@prototypeable
class CdtFile:
    def __init__(self,
                 probes, # Data as CdtRows (including uid and gweight)
                 fieldnames = None, # Column names
                 aids = None, # IDs for mapping columns to atr files
                 eweights = None, # Column weights for column clustering
                 extranames = None, # Extra annotatation column names
                 **kw): # swallow any additional arguments (for prototype)
        
        self.probes = probes

        # N.B.: can't make empty CdtFile
        if(len(self.probes) == 0):
            raise ValueError("No probes!")
        cols = len(self.probes[0])

        if(fieldnames is not None):
            self.fieldnames = list(fieldnames)
        else:
            self.fieldnames = [("V%d" % i) for i in range(cols)]
        assert(len(self.fieldnames) == cols)

        if(extranames is not None):
            self.extranames = list(extranames)
        else:
            self.extranames = []
            
        self.aids = aids
        if(self.aids is not None):
            assert(len(self.aids) == cols)

        if(eweights is not None):
            self.eweights = eweights
        else:
            self.eweights = [1.0]*cols
        assert(len(self.eweights) == cols)

        self.uids = dict(map(lambda x: (x.Uniqid(), x),
                             self.probes))
               
    @classmethod
    def fromCdt(cls, cdtfile):
        """Parse generalized CDT file, based on the description in
        chapter 2 of the JavaTreeView user manual."""
        
        if(isinstance(cdtfile, str)):
            fp = csv.reader(open(cdtfile), dialect = csv.excel_tab)
        else:
            fp = csv.reader(cdtfile, dialect = csv.excel_tab)

        header = next(fp)

        # If the first column of the header is "GID", then it is assigned
        # to the gid vector for mapping to row based dendrograms.  Otherwise,
        # it is assumed to be missing.  In this case, we choose to share
        # the UID field with GID for the purpose of linking rows.

        # The names of the UNIQID and NAME columns are not specified, so
        # they are purely determined by location.
        
        if(header[0] == "GID"):
            gid = 0
            uniqid = 1
            name = 2
        else:
            sys.stderr.write("Warning, didn't find a GID column!  Sharing UID with GID.\n")
            if(header[1] == "NAME"):
                gid = 0
                uniqid = 0
                name = 1
            else:
                sys.stderr.write("Warning, didn't find a NAME column!  Sharing UID with NAME.\n")
                gid = 0
                uniqid = 0
                name = 0

        # If GWEIGHT is given, columns between NAME and GWEIGHT may contain
        # arbitrary annotation data.  Otherwise, GWEIGHT is assumed to be
        # uniform, and ratio data is assumed to start after the NAME column.

        try:
            gweight = header.index("GWEIGHT")
            ratiocols = range(gweight+1,len(header))
            specialcols = range(name+1,gweight)
        except ValueError:
            gweight = None
            ratiocols = range(name+1,len(header))
            specialcols = []

        fieldnames = [header[i] for i in ratiocols]
        extranames = [header[i] for i in specialcols]

        try:
            leaf_offset = extranames.index("LEAF")
        except ValueError:
            leaf_offset = None
        
        probes = []

        aids = None
        eweights = None

        # "Bugs" relative to the JavaTreeView spec:
        #   Not currently supporting special rows other than AID
        #   Currently allowing random mix-in of EWEIGHT and AID

        for fields in fp:
            if(len(fields) < 1):
                continue

            assert(len(fields) == len(header))

            field0 = fields[0]

            args = {"gid": fields[gid], "uniqid": fields[uniqid]}
            if(name):
                args["name"] = fields[name]
            if(gweight):
                args["gweight"] = fields[gweight]

            extra = [fields[i] for i in specialcols]
            fields = [fields[i] for i in ratiocols]
            
            if(field0 == "AID"):
                aids = fields
            elif(field0 == "EWEIGHT"):
                eweights = [float(i) for i in fields]
            else:
                args["ratios"] = fields
                args["extra"] = extra
                args["leaf_offset"] = leaf_offset
                probes.append(CdtRow(**args))

        return cls(probes, fieldnames, aids, eweights,
                   extranames)

    @classmethod
    def fromBagel(cls, bagel):
        """Initialize from BAGEL bar format output or a BagelData object.

        Note that the load step in BagelData drops oligos with insufficient
        data."""
        
        from BagelData import BagelData
        if(not isinstance(bagel, BagelData)):
            if(isinstance(bagel, str)):
                bagel = BagelData(open(bagel))
            else:
                bagel = BagelData(bagel)
        # Find dividing line between expression values and stats
        c = 2
        while(bagel.data.header[c].find("97.5") < 0):
            c += 1

        return cls(probes = [CdtRow(gid = "GENE%dX" % (i+1),
                                    uniqid = row["Unique ID"],
                                    # TODO: Replace w/ annotation
                                    name = row["Unique ID"],
                                    gweight = 1.0,
                                    ratios = [safelog(j) for j in row[2:c]],
                                    extra = row[c:])
                             for (i,row) in enumerate(bagel.data)],
                   fieldnames = bagel.data.header[2:c],
                   extranames = bagel.data.header[c:])

    @classmethod
    def fromMINiML(cls, xml, path, col_names = "Channel",
                   log2 = False):
        """Init CdtFile from 2 channel MINiML data from GEO.
        xml = file handle for GSEX_family.xml
        path = location of external tables
        col_names: where to get column names
        log2: whether to apply log2 transform to the VALUE column
           "Channel": ratio of channel sources
           "Title":   channel titles"""
        # TODO: read directly from tarball rather than directory tree
        # TODO: add support for embedded table data (we currently assume
        #       External-Data for GPL and GPR equivalents)
        import lxml.etree as etree
        NSS = {u"geo": u"http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
        miniml = etree.parse(xml).getroot()
        platforms = miniml.xpath("//geo:Platform", namespaces = NSS)
        if(len(platforms) < 1):
            # Hack for old (pre 6/25/2012) GEO namespace
            NSS = {u"geo": u"http://www.ncbi.nlm.nih.gov/projects/geo/info/MINiML"}
            platforms = miniml.xpath("//geo:Platform", namespaces = NSS)
        
        # Array-level data
        fieldnames = []
        extranames = None

        # Probe-level data
        id2index = {}
        id2name = {}
        id2extra = {}

        for platform in platforms:
            curnames = []
            for(i, col) in enumerate(platform.xpath("geo:Data-Table/geo:Column",
                                                    namespaces = NSS)):
                assert(int(col.attrib["position"]) == i+1)
                name = col.xpath("geo:Name", namespaces = NSS)[0].text.strip()
                if(i == 0):
                    assert(name == "ID")
                else:
                    assert(name != "ID")
                    curnames.append(name)
            if(extranames is None):
                extranames = curnames
            else:
                assert(extranames == curnames)

            # Parse GPL equivalent
            for (i,line) in enumerate(
                open(path+"/"+platform.xpath("geo:Data-Table/geo:External-Data",
                                             namespaces = NSS)[0].text.strip())):
                fields = line.rstrip("\r\n").split("\t")
                try:
                    assert(id2name[fields[0]] == fields[1])
                except KeyError:
                    id2name[fields[0]] = fields[1]

                assert(len(fields) == len(extranames)+1)
                try:
                    assert(id2extra[fields[0]] == fields[1:])
                except KeyError:
                    id2extra[fields[0]] = fields[1:]

        uids = sorted(id2name.keys())
        ratios = []
        uid2index = dict((i,n) for (n,i) in enumerate(uids))

        for sample in miniml.xpath("//geo:Sample", namespaces = NSS):
            try:
                channel = "%s/%s" % (
                    sample.xpath("geo:Channel[@position=1]/geo:Source",
                                 namespaces = NSS)[0].text.strip(),
                    sample.xpath("geo:Channel[@position=2]/geo:Source",
                                 namespaces = NSS)[0].text.strip())
            except:
                channel = "No channel"
                
            try:
                title = sample.xpath("geo:Title",
                                     namespaces = NSS)[0].text.strip()
            except:
                title = "No title"
                
            iid = sample.attrib["iid"]
            
            if(col_names == "Channel"):
                fieldnames.append(channel)
            elif(col_names == "Title"):
                fieldnames.append(title)
            elif(col_names == "ID"):
                fieldnames.append(iid)
            elif(col_names.find("{") > -1):
                fieldnames.append(col_names.format(channel = channel,
                                                   title = title,
                                                   iid = iid))
            else:
                raise ValueError("Bad col_names value: %s" % col_names)
            ratios.append([None]*len(uids))
            table = sample.xpath("geo:Data-Table", namespaces = NSS)[0]
            assert(table.xpath("geo:Column[@position=1]/geo:Name",
                               namespaces = NSS)[0].text.strip() == "ID_REF")
            assert(table.xpath("geo:Column[@position=2]/geo:Name",
                               namespaces = NSS)[0].text.strip() == "VALUE")
            for line in open(
                path+"/"+table.xpath("geo:External-Data",
                                     namespaces = NSS)[0].text.strip()):
                fields = line.rstrip("\r\n").split("\t")
                ratio = fields[1]
                if(log2):
                    ratio = safelog(ratio)
                ratios[-1][uid2index[fields[0]]] = ratio

        ratios = transpose(ratios)

        assert(len(id2name) == len(uids))
        assert(len(ratios) == len(uids))
        assert(len(id2extra) == len(uids))

        return cls(
            probes = [CdtRow(gid = "GENE%dX" % (i+1),
                             uniqid = uid,
                             name = id2name[uid],
                             ratios = r,
                             extra = id2extra[uid])
                      for (i, (uid, r)
                           ) in enumerate(zip(uids, ratios))],
            eweights = [1.0]*len(fieldnames),
            fieldnames = fieldnames,
            extranames = extranames)

    @classmethod
    def fromSingleChannelMINiML(cls, xml, path, sample2fieldname = None,
                                quality_filter = None):
        """Construct a CdtFile from a single-channel, single-platform GEO
           series.  
        
           xml  = Path to MINiML XML file describing the series
           path = Path to external platform and sample tables
           sample2fieldname = function taking sample as lxml element
                              and returning fieldname as a string
                              (default is to use channel 1 "source")
        """
        if(sample2fieldname is None):
            sample2fieldname = lambda x: ("%s" % (
                x.xpath("geo:Channel[@position=1]/geo:Source",
                        namespaces = NSS)[0].text.strip(),))

        import lxml.etree as etree
        NSS = {u"geo": u"http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
        miniml = etree.parse(xml).getroot()
        platforms = miniml.xpath("//geo:Platform", namespaces = NSS)
        if(len(platforms) < 1):
            # Hack for old (pre 6/25/2012) GEO namespace
            NSS = {u"geo": u"http://www.ncbi.nlm.nih.gov/projects/geo/info/MINiML"}
            platforms = miniml.xpath("//geo:Platform", namespaces = NSS)
        
        assert(len(platforms) == 1)
        platform = platforms[0]

        # Array-level data
        fieldnames = []
        extranames = []

        # Probe-level data
        id2index = {}
        uids = []
        names = []
        ratios = []
        extra = []

        for(i, col) in enumerate(platform.xpath("geo:Data-Table/geo:Column",
                                                namespaces = NSS)):
            assert(int(col.attrib["position"]) == i+1)
            name = col.xpath("geo:Name", namespaces = NSS)[0].text.strip()
            if(i == 0):
                assert(name == "ID")
            else:
                assert(name != "ID")
                extranames.append(name)
                extra.append([])

        # Parse GPL equivalent
        for (i,line) in enumerate(
            open(path+"/"+platform.xpath("geo:Data-Table/geo:External-Data",
                                         namespaces = NSS)[0].text.strip())):
            fields = line.rstrip("\r\n").split("\t")
            id2index[fields[0]] = i
            uids.append(fields[0])
            names.append(fields[1])
            assert(len(fields) == len(extra)+1)
            for (e,f) in zip(extra, fields[1:]):
                e.append(f)

        assert(len(id2index) == max(id2index.values())+1)

        for sample in miniml.xpath("//geo:Sample", namespaces = NSS):
            fieldnames.append(sample2fieldname(sample))
            ratios.append([None]*len(id2index))
            table = sample.xpath("geo:Data-Table", namespaces = NSS)[0]
            assert(table.xpath("geo:Column[@position=1]/geo:Name",
                               namespaces = NSS)[0].text.strip() == "ID_REF")
            assert(table.xpath("geo:Column[@position=2]/geo:Name",
                               namespaces = NSS)[0].text.strip() == "VALUE")
            for line in open(
                path+"/"+table.xpath("geo:External-Data",
                                     namespaces = NSS)[0].text.strip()):
                fields = line.rstrip("\r\n").split("\t")
                if((quality_filter is None) or quality_filter(fields)):
                    ratios[-1][id2index[fields[0]]] = fields[1]

        ratios = transpose(ratios)
        extra = transpose(extra)

        assert(len(names) == len(uids))
        assert(len(ratios) == len(uids))
        assert(len(extra) == len(uids))

        return cls(
            probes = [CdtRow(gid = "GENE%dX" % (i+1),
                             uniqid = uid,
                             name = name,
                             ratios = r,
                             extra = e)
                      for (i, (uid, name, r, e)
                           ) in enumerate(zip(uids, names, ratios, extra))],
            eweights = [1.0]*len(fieldnames),
            fieldnames = fieldnames,
            extranames = extranames)

    def __len__(self):
        return len(self.probes)

    def __getitem__(self, i):
        return self.probes[i]

    def GetUid(self, i):
        return self.uids[i]

    def Rows(self):
        return len(self)

    def Cols(self):
        return len(self.fieldnames)        

    def transpose(self):
        """Return a CdtFile with the same heatmap as this one, but
        transposed.  Currently, annotations are discarded (they will
        be restored once we support column annotations."""

        return CdtFile(probes = [
            CdtRow(gid = str(n), uniqid = field, ratios = i[:],
                   gweight = w)
            for (n,(field, i, w)) in enumerate(zip(self.fieldnames,
                                                   transpose(self),
                                                   self.eweights))],
                       fieldnames = [i.Uniqid() for i in self],
                       aids = [i.Gid() for i in self],
                       eweights = [i.Gweight() for i in self])

    def masked_pair(self, uids = None, cols = None,
                    return_rows = False, return_cols = False):
        """Return this data matrix factored into
        two numpy arrays (floats and NA masks) as in PyCluster.
        """

        if(cols == None):
            cols = range(self.Cols())
        if(uids == None):
            rows = [i for i in self]
            uids = [i.Uniqid() for i in rows]
        else:
            rows = [self.GetUid(i) for i in uids]
        
        a = array([[toFloat(i.ratios[j]) for j in cols]
                   for i in rows], dtype = "float")
        m = array([[toMask(i.ratios[j]) for j in cols]
                   for i in rows], dtype = "bool")

        retval = [a,m]
        
        if(return_rows):
            retval.append(rows)
        if(return_cols):
            retval.append(cols)

        return tuple(retval)
    
    def write(self, fp = sys.stdout, eweight = True, gid = True):
        self.writeCdt(fp, eweight, gid)

    def writeCdt(self, fp = sys.stdout, eweight = True, gid = True):
        if(isinstance(fp, str)):
            fp = open(fp,"w")
            
        out = csv.writer(fp, dialect = csv.excel_tab)
        header = []
        if(gid):
            header.append("GID")
        header += (["UNIQID","NAME"]+
                   self.extranames+
                   ["GWEIGHT"]+
                   self.fieldnames)
        out.writerow(header)
        if(self.aids != None):
            out.writerow(["AID",None,None]+
                         [None]*len(self.extranames)+
                         [None]+
                         self.aids)

        if(eweight):
            if(self.eweights is not None):
                eweights = self.eweights
            else:
                eweights = [1.0]*len(self.fieldnames)
            out.writerow(["EWEIGHT",None,None]+
                         [None]*len(self.extranames)+
                         [1.0]+
                         eweights)

        for row in self.probes:
            cur = []
            if(gid):
                cur.append(row.Gid())
            cur += ([row.Uniqid(),row.Name()]+
                    row.Extra()+
                    [row.Gweight()]+
                    row.Ratios())
            out.writerow(cur)

    def sample(self, k, seed = None):
        """Return a CdtFile with a random sample of k rows from
        this CdtFile.  If seed is given, it is used to initialize
        the random number generator.  CdtRow objects are shared
        between this CdtFile and the returned CdtFile."""

        # Instantiate our own generator to avoid re-seeding
        # any generator in the client.
        from random import Random
        rand = Random(x = seed)

        return CdtFile.fromPrototype(
            self,
            probes = rand.sample(self.probes, k))

    def get_sig(self, lfc = 1., fdr = .05, ri = None, direction = "either"):
        """Return a copy of this CdtFile retaining only rows with at least one
        absolute ratio value >= lfc with the corresponding p-value < fdr.
        p-values are taken from the extra columns named as p(x) for each ratio
        column x.  If ri is given, it is a list of strings specifying which
        ratio columns to use.  Otherwise, all ratio columns with a corresponding
        p-value are used.  By default, genes that are up or down in one of the
        contrasts are returned.  Set direction to "up" or "down" to specify
        only up or down regulated genes.
        """
        def either(x, lfc):
            return (x is not None) and (abs(x) >= lfc)
        def up(x, lfc):
            return (x is not None) and (x >= lfc)
        def down(x, lfc):
            return (x is not None) and (x <= -lfc)
        
        isdiff = {"either": either, "up":up, "down":down}[direction]

        if(ri is None):
            pi = [n for (n,i) in enumerate(self.extranames)
                  if(i.startswith("p(") and i.endswith(")"))]
            ri = [self.fieldnames.index(self.extranames[i][2:-1])
                  for i in pi]
        else:
            pi = [self.extranames.index("p(%s)" % self.fieldnames[i])
                  for i in ri]
        return CdtFile.fromPrototype(
            self,
            probes = [i for i in self
                      if(any([(isdiff(i[j],lfc) and (float(i.extra[k]) < fdr))
                              for (j,k) in zip(ri,pi)]))])
    
    def bicluster(self, outprefix, dist = "e", method = "m"):
        """Perform both row and column clustering of this heatmap,
        writing the results to outprefix.{cdt,gtr,atr}.
        """

        (a,m) = self.masked_pair()

        # Column cluster

        d = Pycluster.distancematrix(a, mask = m, dist = dist,
                                     transpose = True)
        tree = Pycluster.treecluster(distancematrix = d, method = method,
                                     transpose = True, data = None)
        tree.scale()
        record = Pycluster.Record()
        # Restore "None" values
        record.data = array([list(i) for i in self])
        record.geneid = [i.Uniqid() for i in self]
        record.genename = [i.Name() for i in self]
        record.gweight = None
        record.gorder = None
        record.expid = self.fieldnames[:]
        record.eweight = self.eweights[:]
        record.eorder = None
        record.uniqid = "UNIQID"

        with TempDir() as tmp_dir:
            
            record.save(os.path.join(tmp_dir,"col"), expclusters = tree)

            # Row cluster

            tree = self.cluster(dist=dist,method=method)
            self.writeCdtGtr(os.path.join(tmp_dir,"row"),tree)

            # Copy trees to permanent locations

            shutil.copy(os.path.join(tmp_dir,"col.atr"),
                        outprefix+".atr")
            shutil.copy(os.path.join(tmp_dir,"row.gtr"),
                        outprefix+".gtr")
                        
            # Read back partial CDTs
            #  (Note: In principal, we only require the trees)

            rc = [i.rstrip("\r\n").split("\t")
                  for i in open(os.path.join(tmp_dir,"row.cdt"))]
            cc = [i.rstrip("\r\n").split("\t")
                  for i in open(os.path.join(tmp_dir,"col.cdt"))]

        # Merge and write final CDT
            
        di = rc[0].index("GWEIGHT")+1
        # Map column names to their original offsets
        name2offset = dict((i,n) for (n,i) in enumerate(rc[0][di:]))

        out = open(outprefix+".cdt","w")
        # column header
        out.write("\t".join(rc[0][:di]+cc[0][3:])+"\n")
        # column IDs
        out.write("\t".join(["AID"]+[""]*(di-1)+cc[1][3:])+"\n")
        # eweights
        out.write("\t".join(["EWEIGHT"]+[""]*(di-1)+cc[2][3:])+"\n")
        for i in rc[2:]:
            out.write("\t".join(i[:di]+[i[name2offset[j]+di]
                                        for j in cc[0][3:]])+"\n")
        out.close()
            
    def cluster(self, outprefix = None, cols = None, uids = None,
                dist = "e", method = "m", distfile = None):
        """Return a Pycluster.Tree object corresponding to a hierarchical
        clustering of the CDT data with default parameters.

        If outprefix is given, CDT and GTR output files will be generated
        using Pycluster's output functions (so extended annotation columns
        will be dropped).

        If cols is given (as an iterable of integers) it specifies a subset
        of columns to consider for the clustering.

        If uids is given (as an iterable of strings) it specifies a subset
        of rows to cluster.

        The indexing of the tree (as generated by Pycluster._savetree)
        corresponds to the ordering of the uids list, if given, or else the
        traversal order of the CDT rows.

        distfile, if given, is a filename or filestream to which the lower
        triangle of of the distance matrix is written in CSV format.
        """
        # Note: could add filter here
        sys.stderr.write("Building array...\n")

        (a,m,rows,cols) = self.masked_pair(uids = uids, cols = cols,
                                           return_rows = True,
                                           return_cols = True)
        
        sys.stderr.write("Building distance matrix...\n")
        if(self.eweights is None):
            d = Pycluster.distancematrix(a, mask = m, dist = dist)
        else:
            d = Pycluster.distancematrix(
                a, mask = m, dist = dist,
                weight = array([self.eweights[j] for j in cols]))

        if(distfile is not None):
            if(isinstance(distfile, str)):
                fp = open(distfile,"w")
            else:
                fp = distfile
            # Distance matrix is a list of 1D arrays giving the lower
            #   triangle of the matrix, not including the diagonal.
            #   (with first element length zero)
            out = csv.writer(fp)
            for row in d:
                out.writerow(row)
            del out
            if(isinstance(distfile, str)):
                fp.close()

        sys.stderr.write("Clustering...\n")
        tree = Pycluster.treecluster(distancematrix = d, method = method,
                                     data = None)
        tree.scale()

        if(outprefix is not None):        
            sys.stderr.write("Writing...\n")
            record = Pycluster.Record()
            # Restore "None" values
            record.data = array([[i.ratios[j] for j in cols]
                                       for i in rows])
            record.geneid = [i.Uniqid() for i in rows]
            record.genename = [i.Name() for i in rows]
            record.gweight = None
            record.gorder = None
            record.expid = [self.fieldnames[j] for j in cols]
            record.eweight = [self.eweights[j] for j in cols]
            record.eorder = None
            record.uniqid = "UNIQID"
            record.save(outprefix, geneclusters = tree)

        return tree

    def writeCdtGtr(self, prefix, tree, uids = None):
        """Given a gene tree from Pycluster, write the corresponding
        prefix.cdt and prefix.gtr files with matching GIDs."""

        if(uids == None):
            uids = [i.Uniqid() for i in self]

        # Write GTR
        geneindex = _savetree(prefix, tree, arange(len(self)), 0)

        # Write CDT
        out = csv.writer(open(prefix+".cdt","w"), dialect = csv.excel_tab)
        out.writerow(["GID","UNIQID","NAME"]+
                     self.extranames+
                     ["GWEIGHT"]+
                     self.fieldnames)
        if(self.aids != None):
            pass
        out.writerow(["EWEIGHT",None,None]+
                     [None]*len(self.extranames)+
                     [1.0]+
                     # Note that we are killing any existing eweight values
                     [1.0]*len(self.fieldnames))

        for i in geneindex:
            row = self.GetUid(uids[i])
            out.writerow(["GENE%dX" % i,row.Uniqid(),row.Name()]+
                         row.Extra()+
                         [row.Gweight()]+
                         row.Ratios())
        del out

    def add_extra_column(self, name, annotations):
        """Return a CdtFile with new annotation column "name"
        right appended.  annotations is either a list of annotation
        strings with the same ordering as the rows of the CdtFile
        or a dictionary of {uid:string} with uids matching a subset
        of the rows of the CdtFile (missing uids will have empty
        annotations).
        """
        if(hasattr(annotations, "keys")):
            x = []
            for i in self:
                try:
                    val = annotations[i.Uniqid()]
                except KeyError:
                    val = ""
                x.append(val)
            annotations = x
            
        return CdtFile.fromPrototype(
            self,
            probes = [CdtRow.fromPrototype(i, extra = i.extra+[j])
                      for (i,j) in zip(self, annotations)],
            extranames = self.extranames + [name])

    def add_ratio_column(self, name, ratios, eweight = 1.):
        """Return a CdtFile with new heatmap column "name"
        right appended.  ratios is either a list of float|None
        with the same ordering as the rows of the CdtFile
        or a dictionary of {uid:float|None} with uids matching a subset
        of the rows of the CdtFile (missing uids will be set to None).
        """
        if(hasattr(ratios, "keys")):
            x = []
            for i in self:
                try:
                    val = ratios[i.Uniqid()]
                except KeyError:
                    val = None
                x.append(val)
            ratios = x
            
        return CdtFile.fromPrototype(
            self,
            probes = [CdtRow.fromPrototype(i, ratios = i.ratios+[j])
                      for (i,j) in zip(self, ratios)],
            fieldnames = self.fieldnames + [name],
            eweights = self.eweights + [eweight])

    def median_normalize_cols(self):
        """Return a (median) column-centered version of this heatmap."""

        medians = [median(i[j] for i in self) 
                   for j in range(len(self.fieldnames))]

        return CdtFile.fromPrototype(
            self,
            probes = [CdtRow.fromPrototype(
                i, ratios = [safesub(j,k) for (j,k) in zip(i,medians)])
                      for i in self])

    def mean_normalize_rows(self):
        """Return a (mean) row-centered version of this heatmap.

        This transform is particularly useful for factoring out
        mean expression levels from (log(FPKM)) RnaSeq data.
        """
        def rownorm(ratios):
            m = safemean(ratios)
            return [safesub(i,m) for i in ratios]

        return CdtFile.fromPrototype(
            self, 
            probes = [CdtRow.fromPrototype(i, ratios = rownorm(i.ratios))
                      for i in self])

    def median_normalize_rows(self):
        """Return a (median) row-centered version of this heatmap.
        """
        def rownorm(ratios):
            m = median(ratios)
            return [safesub(i,m) for i in ratios]

        return CdtFile.fromPrototype(
            self, 
            probes = [CdtRow.fromPrototype(i, ratios = rownorm(i.ratios))
                      for i in self])

class GtrNode:
    def __init__(self, gid, left = None, right = None,
                 length = 0.0, pos = None, uid = None,
                 name = None, annotation = None):
        self.gid = gid
        self.left = left
        self.right = right
        self.length = float(length)
        self.uid = uid
        self.name = name
        self.annotation = annotation

        if(self.uid == ""):
            self.uid = None
        if(self.name == ""):
            self.name = None
        if(self.annotation == ""):
            self.annotation = None

        if((self.left == None) and
           (self.right == None)):
            self.depth = 0
            self.pos = pos

        else:
            self.depth = max(self.left.Depth(),
                             self.right.Depth())+1
            if(pos != None):
                self.pos = pos
            else:
                self.pos = float(self.left.Pos() + self.right.Pos())/2.0

    def Gid(self):
        return self.gid

    def Uid(self):
        return self.uid

    def Depth(self):
        return self.depth

    def Pos(self):
        return self.pos

    def Name(self):
        return self.name

    def Annotation(self):
        return self.annotation

    def __str__(self):
        return " ".join(filter(None, (self.Gid(), self.Uid(), self.Name())))

    def DfsIterator(self, chain = False, onPop = None):
        """Generator that yields this node and its children as a
        depth-first-search traversal."""
        node = self
        stack = []
        visited = set()

        while(1):
            if(node not in visited):
                visited.add(node)
                if(chain):
                    yield (node, stack)
                else:
                    yield node

            elif((node.left != None) and
                 (node.left not in visited)):
                stack.append(node)
                node = node.left

            elif((node.right != None) and
                 (node.right not in visited)):
                stack.append(node)
                node = node.right

            else:
                if(onPop != None):
                    if(chain):
                        onPop(node, stack)
                    else:
                        onPop(node)
                        
                if(len(stack) > 0):
                    node = stack[-1]
                    stack = stack[:-1]
                else:
                    break

class GtrFile:
    def __init__(self, gtrfile, cdt):
        # TODO: implement tree in terms of nx

        if(isinstance(gtrfile, str)):
            fp = csv.reader(open(gtrfile), dialect = csv.excel_tab)
        else:
            fp = csv.reader(gtrfile, dialect = csv.excel_tab)

        self.cdt = cdt

        self.nodes = {}
        self.levels = {0:[]}

        for (i, leaf) in enumerate(self.cdt):
            self.nodes[leaf.Gid()] = node = GtrNode(leaf.Gid(), pos = i,
                                                    uid = leaf.Uniqid())
            self.levels[0].append(node)

        # TODO: Add support for extended tree file format by sniffing
        #       for header line.

        # GTR mode: 0 --> TreeFile format
        #           1 --> JTV extended TreeFile format
        mode = 0

        for (i, fields) in enumerate(fp):
            # Note: this is the magic for recognizing an extended
            #       GTR file that is described on page 21 of the
            #       JavaTreeView user's manual.
            if((i == 0) and (fields[0] == "NODEID")):
                columns = dict((k.strip().upper(),i)
                               for (i,k) in enumerate(fields))
                mode = 1
                continue
            try:
                if(mode == 1):
                    gid = fields[columns["NODEID"]]
                    left = fields[columns["LEFT"]]
                    right = fields[columns["RIGHT"]]
                    try:
                        length = fields[columns["CORRELATION"]]
                    except KeyError:
                        length = fields[columns["TIME"]]
                    try:
                        name = fields[columns["NAME"]]
                    except KeyError:
                        name = None
                    try:
                        annotation = fields[columns["ANNOTATION"]]
                    except KeyError:
                        annotation = None
                else:
                    (gid, left, right, length) = fields[:4]
            except ValueError:
                print("Bad GTR line:")
                print(i,fields)
                continue

            if(mode == 1):
                self.nodes[gid] = node = GtrNode(
                    gid, self.nodes[left], self.nodes[right], length,
                    name = name, annotation = annotation)
            else:
                self.nodes[gid] = node = GtrNode(
                    gid, self.nodes[left], self.nodes[right], length)
            try:
                self.levels[node.Depth()].append(node)
            except KeyError:
                self.levels[node.Depth()] = [node]

        self.depth = max(self.levels.keys())+1

        #print "Parsed %d nodes on %d levels" % (len(self.nodes),
        #                                        self.depth)

        
    def Depth(self):
        return self.depth

    def Root(self):
        return self.levels[self.depth - 1][0]

class CdtTransform:
    """A linear transform on a CdtFile."""
    def __init__(self, domain_names, range_names, matrix):
        """Initialize from a linear transform.

           domain_names: field names of input vector (length M)
           range_names:  field names of output vector (length N)
           matrix:       matrix of N rows and M columns, giving the
                         linear transform.
        """
        self.domain_names = tuple(domain_names)
        self.range_names = list(range_names)
        self.matrix = matrix
        assert(len(matrix) == len(self.range_names))
        assert(len(matrix[0]) == len(self.domain_names))

    def __call__(self, cdt):
        """Return a transformed copy of cdt.

        The results of the transform are appended as additional ratio
        columns, using range_names from __init__.  This is a non-destructive
        operation (i.e., the "domain" columns are retained).
        """
        indices = [cdt.fieldnames.index(i) for i in self.domain_names]
        return CdtFile.fromPrototype(
            cdt,
            fieldnames = cdt.fieldnames + self.range_names,
            eweights = cdt.eweights + [1.0]*len(self.range_names),
            probes = [
            CdtRow.fromPrototype(
            probe,
            ratios = probe.ratios + [sum(probe[i]*term
                                         for (i,term) in zip(indices,vector))
                                     for vector in self.matrix])
                      for probe in cdt])

if(__name__ == "__main__"):
    print("Hello, world")
