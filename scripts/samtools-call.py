#!/usr/bin/env python

import os
import argparse
import subprocess

GRCH38_PATH = os.environ["GRCH38_PATH"]


def mpileup(bam_file, output_path, max_read_depth=100):
    opts = {}
    opts['bam_file'] = bam_file
    opts['final_vcf'] = bam_file[:-11] + '.vcf'
    opts['intermediary'] = bam_file + '.raw.bcf'
    opts["ref_genome"] = GRCH38_PATH
    opts["max_read_depth"] = max_read_depth
    proc = subprocess.Popen(
        ("samtools mpileup -ugf {ref_genome} {bam_file}"
         " | bcftools call -vm -o {final_vcf}").format(**opts), shell=True)
    proc.wait()
    #TODO Add Error Check
    return opts['final_vcf']

def guess_outpath(bam_file):
    dir = os.path.dirname(bam_file)
    base = os.path.basename(bam_file)
    name = os.path.splitext(base)[0]
    if name.endswith("sorted"):
        name = os.path.splitext(name)[0]
    return "%s%s%s.vcf" % (dir, os.sep, name)

app = argparse.ArgumentParser("sam-call")
app.add_argument("bam_path")
app.add_argument("-o", "--out-path", required=False, default=None)



if __name__ == "__main__":
    import sys
    args = app.parse_args()
    if args.out_path is None:
        args.out_path = guess_outpath(args.bam_path)
    print args
    print mpileup(args.bam_path, args.out_path)
