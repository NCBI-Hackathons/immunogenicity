import os
import re
import json
from collections import defaultdict

# 3rd Party Imports
VCFFilterBase = None
_Filter = None
class PyVCFStub(object):
    def __init__(self, *args, **kwargs):
        pass
    @classmethod
    def __call__(self, *args, **kwargs):
        raise NotImplementedError("These features are not available. Please install PyVCF: https://github.com/jamescasbon/PyVCF")
try:
    import vcf
    from vcf.parser import _Filter
    from vcf.filters import Base as VCFFilterBase
except ImportError, e:
    # Stubbing so if the user does not install vcf they can still see the help menu from pathoscope
    vcf = PyVCFStub()
    vcf.model = PyVCFStub()
    vcf.model._Record = PyVCFStub()
    VCFFilterBase = PyVCFStub
    _Filter = PyVCFStub

from pathovar.utils import defline_parser

class AnnotatedVariant(vcf.model._Record):
    def __init__(self, _record, annotations = None):
        if annotations == None:
            annotations = []
        vcf.model._Record.__init__(self, _record.CHROM, _record.POS, _record.ID, _record.REF,
         _record.ALT, _record.QUAL, _record.FILTER, _record.INFO, _record.FORMAT, 
         _record._sample_indexes, _record.samples)
        self.INFO['GENE'] = [x.to_info_field() for x in annotations]
        self.annotations = annotations

    def __repr__(self):
        rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
        return rep

    def __str__(self):
        rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
        return rep

def purify_vcf(file_path, output_file = None):
    if output_file is None:
        output_file = file_path + '.pure.vcf'
    inp = vcf.Reader(open(file_path))
    out = vcf.Writer(open(output_file, 'w'), inp)
    for var in inp:
        if len(var.FILTER) == 0:
            out.write_record(var)
    return output_file


## filter_vcf
# Based on vcf_filter.py from PyVCF. Compatible with instantiated
# PyVCF _Filter objects.
# @param keep If True, keep sequences that pass a filter, else exclude sequences that pass the filter
def filter_vcf(file_path, filters, short_circuit = False, drop_filtered = False, invert = False, output_file = None):
    if output_file is None:
        output_file = file_path + '.filt.vcf'
    inp = vcf.Reader(open(file_path, 'r'))

    # build filter chain
    chain = []
    for filter_obj in filters:
        chain.append(filter_obj)
        short_doc = filter_obj.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        # add a filter record to the output
        inp.filters[filter_obj.filter_name()] = _Filter(filter_obj.filter_name(), short_doc)

    # output must be created after all the filter records have been added
    output = vcf.Writer(open(output_file, 'w'), inp)

    # apply filters
    for record in inp:
        output_record = True
        for filt in chain:
            result = filt(record)
            if (result is None) or (result is not None and invert): continue

            # save some work by skipping the rest of the code
            if drop_filtered:
                output_record = False
                break

            record.add_filter(filt.filter_name())
            if short_circuit: break
        
        # If the record is to be kept (not dropping filtered, or record passed all filters)
        if output_record:
            output.write_record(record)
    return output_file

# Only want filters to return a non-None value when the record fails to pass their filtering constraints
def filter_vcf_in_memory(variants, filters, short_circuit = False, drop_filtered = False, invert = False, **kwargs):

    chain = []
    for filter_obj in filters:
        chain.append(filter_obj)
        short_doc = filter_obj.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        # add a filter record to the output
        try:
            variants.filters[filter_obj.filter_name()] = _Filter(filter_obj.filter_name(), short_doc)
        except:
            pass
    filtered_records = []
    for record in variants:
        output_record = True
        for filt in chain:
            result = filt(record)
            if (result is None) or (result is not None and invert): continue

            # save some work by skipping the rest of the code
            if drop_filtered:
                output_record = False
                break

            record.add_filter(filt.filter_name())
            if short_circuit: break
        
        # If the record is to be kept (not dropping filtered, or record passed all filters)
        if output_record:
            filtered_records.append(record)

    return filtered_records

def get_variant_genes(vcf_path):
    genes = {}
    variant_by_gene = defaultdict(list)
    intergenic_variants = defaultdict(list)
    reader = vcf.Reader(open(vcf_path,'r'))
    for var in reader:
        gene_info = var.INFO.get('GENE', '')
        if re.search(r"Intergenic", gene_info ):
            idents = defline_parser(var.CHROM)
            org_name = idents.get('org_name', var.CHROM).replace('_', ' ')
            gid = idents['gene_id']
            intergenic_variants[org_name + ":" + gid].append(var)
        else:
            gene_info_groups = gene_info.split('||')
            general_info = gene_info_groups[0].split('|')
            
            # Fields are key:value pairs. Certain fields have special handling
            # The first position will have a leading paren to exclude
            general_info[0] = general_info[0][1:]
            general_info = [pair.split(':') for pair in general_info]
            gene_dict = {val[0]: val[1] for val in general_info}
            for key in gene_dict:
                if key not in ('acc', 'gi'):
                    gene_dict[key] = gene_dict[key].replace('_', ' ')
            genes[gene_dict['gi']] = gene_dict
            variant_by_gene[gene_dict['gi']].append(var)

    return(genes, variant_by_gene, intergenic_variants)

def vcf_to_gene_report(vcf_path):
    genes, variant_by_gene = get_variant_genes(vcf_path)

    report = open(os.path.splitext(vcf_path)[0] + '.variant_report.tsv', 'w')
    report.write('Gene\tVariant Positions\n')
    for gene in genes.values():
        line = "gi|{gi}|" 
        if 'acc' in gene:
            line +="ref|{acc}|" 
        line += " {title}"
        line = line.format(**gene)
        for var in variant_by_gene[gene['gi']]:
            line += '\t({POS}:{REF}->{ALT})'.format(**var.__dict__)
        report.write(line + '\n')
    report.close()

class VCFFilterException(Exception):
    pass

## VCF Filter Classes
class FilterByComparisonVCF(VCFFilterBase):
    '''Filter a VCF File by comparing its sites to another VCF File and operating on the intersection/difference'''
    name = ''
    default_ref_vcfs = None
    default_intersection = False

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--ref-vcfs', type=str, nargs="+", help="A VCF file against which to compare (may specify more than one)")
        parser.add_argument('--intersection', action = "store_true", default=None, help="Instead of excluding intersecting sites, keep them and drop sites not found in both files.")

    def __init__(self, args):
        self.reference_vcfs = args.ref_vcfs 
        if self.reference_vcfs is None:
            raise VCFFilterException("No reference VCFs provided to FilterByComparisonVCF")
            self.reference_vcfs = []
        self.reference_variants =[]
        for ref_vcf in self.reference_vcfs:
            self.reference_variants += [var for var in vcf.Reader(open(ref_vcf))]
        self.intersection = args.intersection if args.intersection else FilterByComparisonVCF.default_intersection
        self.reference_dict = defaultdict(lambda : defaultdict(bool))
        for var in self.reference_variants:
            self.reference_dict[var.CHROM][(var.start, var.end, tuple(map(str, var.ALT)))] = True
        self.count = 0

    def __call__(self, record):
        res = self.reference_dict[record.CHROM][(record.start,record.end, tuple(map(str, record.ALT)))]
        if self.intersection and not res:
            return record
        elif res:
            return record

    def filter_name(self):
        mode = "conc-" if self.intersection else "dis-"
        return self.name + mode + '-'.join(self.reference_vcfs)

class FilterByChromMatch(VCFFilterBase):
    '''Filter a VCF File by regular expresion match over its CHROM column'''
    name = "regmatch"
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--pattern', type=str, required = False,
                help='Regular Expression to match the CHROM column')

    def __init__(self, args):
        if args.pattern is None:
            raise VCFFilterException("Missing init value")
        self.pattern = args.pattern

    def __call__(self, record):
        if not re.search(self.pattern, record.CHROM):
            return record

    def filter_name(self):
        return "%s-/%s/" % (self.name, self.pattern)

    def __repr__(self):
        return self.filter_name()

class FilterByAltCallDepth(VCFFilterBase):
    '''Filter a VCF File by limiting variants to only those with at least X% Alt calls'''
    
    name = 'alt-call-depth'
    default_alt_depth = 0.4


    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--alt-depth', type=float, default=None, 
            help="The minimum percentage of all calls for a locus that must be an alternative allele (MAF)")
    
    def __init__(self, args):
        if args.alt_depth is None:
            raise VCFFilterException("Missing init value")
        self.alt_depth = args.alt_depth if args.alt_depth else FilterByAltCallDepth.default_alt_depth

    def filter_name(self):
        return "%s-%.2f" % (self.name, self.alt_depth)

    def __call__(self, record):
        total_reads = sum(record.INFO["DP4"])
        ref_reads = sum(record.INFO["DP4"][:2])
        alt_reads = sum(record.INFO["DP4"][2:])

        if not (alt_reads / float(total_reads) >= self.alt_depth):
            return record

class FilterByReadDepth(VCFFilterBase):
    name = 'call-depth'
    default_min_depth = 5
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-depth', type=int, default = None, 
            help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")


    def __init__(self, args):
        if args.min_depth is None:
            raise VCFFilterException("Missing init value")
        self.min_depth = args.min_depth if args.min_depth else FilterByReadDepth.default_min_depth

    def filter_name(self):
        return "%s-%d" % (self.name, self.min_depth)

    def __call__(self, record):
        if not (sum(record.INFO['DP4']) >= self.min_depth):
            return record

class FilterByCallQuality(VCFFilterBase):
    name = 'call-qual'
    default_min_qual = 20
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-qual', type=float, default=None, 
            help="The minimum Phred scale Quality score of a variant call")

    def __init__(self, args):
        if args.min_qual is None:
            raise VCFFilterException("Missing init value")
        self.min_qual = args.min_qual

    def filter_name(self):
        return "%s-%f" % (self.name, self.min_qual)

    def __call__(self, record):
        if not record.QUAL >= self.min_qual:
            return record

class FilterByMappingQuality(VCFFilterBase):
    name = 'map-qual'
    default_min_mq = 20
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-mq', type=float, default=None, 
            help="The minimum Phred scale Quality score of a read mapping for a variant call")

    def __init__(self, args):
        if args.min_mq is None:
            raise VCFFilterException("Missing init value")
        self.min_qual = args.min_mq

    def filter_name(self):
        return "%s-%f" % (self.name, self.min_qual)

    def __call__(self, record):
        if not record.INFO['MQ'] >= self.min_qual:
            return record

EXPOSED_FILTERS = [FilterByComparisonVCF,FilterByAltCallDepth,FilterByReadDepth,FilterByCallQuality,FilterByMappingQuality,]
def main():
    import argparse
    arg_parser = argparse.ArgumentParser(prog='filter-vcf')
    for filter_type in EXPOSED_FILTERS:
        filter_type.customize_parser(arg_parser)
    arg_parser.add_argument('vcf_file', help="Target VCF to filter")
    arg_parser.add_argument('-d', '--drop', action = 'store_true', default = False, help = "Drop variants that fail any filters [default=False]")
    arg_parser.add_argument('-p', '--purify', action = 'store_true', default = False, help = "Do not apply any filters to the input VCF, but drop all\
     rows that have not passed previous filters")
    arg_parser.add_argument('-o', dest='output_file', default=None, help="The name of the output file. Defaults to the input file name + '.filt.vcf'")
    #arg_parser.add_argument('-c', '--config', default = None, help = "Read parameters from a pathovar configuration file?")
    args = arg_parser.parse_args()

    out = None
    if not args.purify:
        filters = []
        for filter_type in EXPOSED_FILTERS:
            try:
                filt = filter_type(args)
                filters.append(filt)
            except VCFFilterException, e:
                pass
                # print("Filter Failed", filter_type, e)
                # Ignore failing to build a filter, since it means that the filter 
                # was not initialized properly, and it should not be used in that case.
        out = filter_vcf(args.vcf_file, filters, drop_filtered = args.drop, output_file = args.output_file)
    else:
        out = purify_vcf(args.vcf_file, output_file = args.output_file)
    print out

if __name__ == '__main__':
    main()

#END            