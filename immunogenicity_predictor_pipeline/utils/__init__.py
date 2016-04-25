import re

__all__ = ['fasta_utils', 'vcf_utils', 'sam_utils']


# Namespace
# For compatibility with classes expecting
# input from argparse


class Namespace:
    def __init__(self):
        pass

    def __repr__(self):
        return str(self.__dict__)


##
gene_id = re.compile(r'gi\|([^\|]+)\|')
##
genbank_id = re.compile(r'gb\|([^\|]+)\|')
ref_id = re.compile(r'ref\|([^\|]+)\|')
##
tax_id = re.compile(r'ti\|([^\|]+)\|')
##
org_name = re.compile(r'org[^\|]*\|([^\|]+)\|')
##
star_pat = re.compile(r"^\*$")

##
#


def defline_parser(line):
    result = {}
    org_name_matches = org_name.findall(line)
    if len(org_name_matches) > 0:
        result['org_name'] = org_name_matches[0]
    tax_id_matches = tax_id.findall(line)
    if len(tax_id_matches) > 0:
        result['tax_id'] = tax_id_matches[0]
    gene_id_matches = gene_id.findall(line)
    if len(gene_id_matches) > 0:
        result['gene_id'] = gene_id_matches[0]
    ref_id_matches = ref_id.findall(line)
    if len(ref_id_matches) > 0:
        result['ref_id'] = ref_id_matches[0]
    genbank_id_matches = genbank_id.findall(line)
    if len(genbank_id_matches) > 0:
        result['genbank_id'] = genbank_id_matches[0]
    # * represents unaligned
    if line == "*":
        result['org_name'] = "*"
        result['gene_id'] = "*"
        result['tax_id'] = "*"
    if len(result.keys()) == 0:
        raise DeflineUnparsedException("Could Not Parse Defline: %s" % line)
    return result


class DeflineUnparsedException(Exception):
    pass
