import re
from collections import defaultdict

from .utils import defline_parser


class SamParser(object):
    def __init__(self, file_name, out= None, **opts):
        self.file_name = file_name
        self.opts = opts
        if out != None:
            self.opts['out'] = out
        self.reference_headers = {"*": NullHeader()}
        self.meta_headers = []
        self.reads = defaultdict(list)
        self.line = None
        self.parse_file()

    def parse_file(self):
        with open(self.file_name, 'r') as handle:
            for line in handle:
                line = line.replace("\n", '')
                self.line = line
                if "@SQ" in line:
                    result = SamHeader.parse(line)
                    try:
                        self.reference_headers[result.definition.split(":")[1]] = result
                    except IndexError, e:
                        # stub to fix missing key
                        self.reference_headers[result.definition] = result
                elif "@" in line[0]:
                    result = SamHeader.parse(line)
                    self.meta_headers.append(result)
                else: 
                    result = SamRead.parse(line)
                    self.reads[result.rname].append(result)

    def alter_reference(self, target, replacement):
        new_headers = {}
        for name, header in self.reference_headers.items():
            if re.search(target, name):
                header.definition = re.sub(target, replacement, name)
                new_headers[re.sub(target, replacement, name)] = header
            else:
                new_headers[name] = header
        new_reads = {}
        for name, reads in self.reads.items():
            if re.search(target, name):
                for read in reads:
                    read.rname = re.sub(target, replacement, name)
                new_reads[re.sub(target, replacement, name)] = reads
            else:
                new_reads[name] = reads

        self.meta_headers[0].rest = ["SO:unsorted"]
        self.reference_headers = new_headers
        self.reads = new_reads

    def remove_reference(self, target, keep = False):
        new_headers = None
        new_reads = None
        if not keep:
            new_headers = {name:header for name, header in self.reference_headers.items() if not re.search(target, name)}
            new_reads = {name:reads for name, reads in self.reads.items() if not re.search(target,name)}
        else: 
            new_headers = {name:header for name, header in self.reference_headers.items() if re.search(target, name)}
            new_reads = {name:reads for name, reads in self.reads.items() if re.search(target,name)}

        self.meta_headers[0].rest = ["SO:unsorted"]
        self.reference_headers = new_headers
        self.reads = new_reads

    def write_output(self):
        outfile = None
        if  "out" in self.opts:
            outfile = self.opts['out']
        else:
            outfile = self.file_name[:-3] + 'filt.sam'
        with open(outfile, 'w') as handle:
            handle.write(str(self.meta_headers[0]) + "\n")
            for head_line in self.reference_headers.values():
                if ("is_null" in head_line.fields): 
                    continue
                handle.write(str(head_line) + "\n")
            for head_line in self.meta_headers[1:]:
                handle.write(str(head_line) + "\n")
            for rname in self.reads:
                for read in self.reads[rname]:
                    handle.write(str(read) + "\n")
        return outfile

class SamHeader(object):
    @classmethod
    def parse(self, line):
        fields = line.split('\t')
        result = SamHeader(fields)
        return result

    def __init__(self, fields):
        self.tag = fields[0]
        self.definition = fields[1]
        self.fields = {}
        if len(fields) > 2:
            self.rest = fields[2:]
            # Capture fields without caring what they are.
            for field in fields[2:]:
                fields_data = field.split(":")
                fname = fields_data[0]
                fval = fields_data[1:]
                self.fields[fname] = "".join(fval).replace('\n','')
                # Whitelisted Integer Fields
                if fname in ['LN']:
                    self.fields[fname] = int(self.fields[fname])
        else: 
            self.rest = []

    def __repr__(self):
        rest = '\t'.join(self.rest)
        rep = "{tag}\t{definition}\t{rest}".format(**{"tag": self.tag, "definition": self.definition, "rest": rest})
        return rep

class NullHeader(SamHeader):
    def __init__(self):
        self.tag = "@SQ"
        self.definition = "*"
        self.fields = {"LN": float("inf"), "is_null": True}
        self.rest = [""]



class SamRead(object):
    @classmethod
    def parse(self, line):
        fields = line.split("\t")
        aln =  (fields[11:])
        fields = fields[:11]
        fields.append(aln)
        result = SamRead(*fields)
        return result

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, aln):
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = int(pos)
        self.mapq = float(mapq)
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = int(tlen)
        self.seq = seq
        self.qual = qual
        self.aln = '\t'.join(aln)

    def __repr__(self):
        rep = "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{aln}".format( **self.__dict__ )
        return rep
