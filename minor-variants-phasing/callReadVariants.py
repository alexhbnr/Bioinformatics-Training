#!/usr/bin/env python
##############################################################################
# Returns for every PacBio read the alleles for all variant positions provided
# by a list. The reads are processed by PacBio's SMRT portal and exported as
# h5 files.
##############################################################################

from __future__ import division, print_function

import argparse, sys
import itertools
import re
import csv
import warnings

import numpy as np
import pandas as pd
import h5py

from operator import itemgetter
from pbcore.io import CmpH5Reader, FastaReader

from Bio import pairwise2

baseDict = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W",
    "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B",
    "ACGT": "N", "": ""}
gapBase = "-"

def main(args=sys.argv[1:]):
    ''' Analyses each read from a PacBio sequencing run provided in a hp5 container
        and returns the alleles for all list of provided variants in a tabular format
    '''
    options = parseArguments(args)
    cmpH5 = loadCmpH5(options.cmp)
    variantsCsv = pd.read_csv(options.variants)
    contigVariants = contigVariantPositions(variantsCsv)

    result = callVariants(cmpH5, contigVariants, options.reference, options.radius)
    with open(options.output, "wb") as csvfile:
        output = csv.writer(csvfile, delimiter=',')

        for row in result:
            for i in range(0, len(result[row])):
                a = result[row][i]
                a.insert(0, row)
                output.writerow(a)

def loadCmpH5(filename):
    ''' Returns a instance of the CMP H5 Reader from PacBio's python library
    '''
    filename = os.path.abspath(os.path.expanduser(filename))
    handle = h5py.File(filename, "r")
    return CmpH5Reader(handle)

def contigVariantPositions(variantsCsv):
    ''' Extracts the contigs, positions, reference and alternative alleles from
        the variants file and builds a dictionary of the form:
            {contig: [(position, refBase, altBases)]}
        The positions are stored in 0-indexed format.
    '''
    contigRefPositions = list(variantsCsv[['CONTIG', 'POSITION', 'REF', 'ALT']].values.tolist())
    byContig = itertools.groupby(contigRefPositions, key=itemgetter(0))
    contigVariants = dict(
        (contig, sorted(set(
            (p - 1, "".join(uniq(cpra[2] for cpra in posData)), "".join(uniq(cpra[3] for cpra in posData)))
            for p, posData in (
                (a, list(b))
                for a, b in itertools.groupby(contigData, key=itemgetter(1)))), key=itemgetter(0)))
        for contig, contigData in byContig)
    return contigVariants

def uniq(iterable):
    ''' Removes duplicates from a sorted list
    '''
    return sorted(set(iterable))

def callVariants(cmpH5, contigVariants, referenceFasta, radius):
    ''' Iterates over all contigs in the reference FastA file and extracts all
        reads that aligned against these contigs. For each aligned read the 
        function to determine the variant bases is invoked.
    '''
    readVariants = {}
    for contigName, variantSites in contigVariants.iteritems():
        contigName = contigName.strip("\"")
        with FastaReader(referenceFasta) as fr:
            for contig in fr:
                if contig.name == contigName:
                    refSeq = contig.sequence
                    try:
                        refInfo = cmpH5.referenceInfo(contigName)
                        refId = refInfo.ID
                        refLen = refInfo.Length
                        readVariants[contigName] = [["READ"] + 
                                                    ["{0:s}{1:d}{2:s}".format(ref, pos + 1, baseDict.get(alt, "X")) for pos, ref, alt in variantSites]]
                        readVariants[contigName].extend(
                            callReadVariants(refSeq, cmpH5[readRow], variantSites, radius, nRead)
                            for nRead, readRow
                            in enumerate(cmpH5.readsInRange(refId, 0, refLen, justIndices=True), start=1))
                    except KeyError:
                        warnings.warn("Cannot find contig: \"{0:s}\" in the cmpH5: \"{1:s}\"".format(contig, cmpH5.filename))
                else:
                    warnings.warn("Cannot find contig: \"{}\" ".format(contigName) +
                                  "in the FASTA: \"{}".format(referenceFasta))
    return readVariants

def callReadVariants(reference, read, variantSites, radius, nRead):
    ''' Returns for each read the variant bases by extracting short region around
        the defined variant sites (length determined by the radius), doing a 
        pairwise global alignment and identifying the alignment with the highest
        score and the shortest length.
    '''
    readName = ccsName(read) if read.cmpH5.readType == "CCS" else read.readName
    resultRow = [readName]

    readString = read.read(orientation="genomic")
    refPositions = read.referencePositions(orientation="genomic")
    alnLen = len(readString)

    def baseInWindow(start, end):
        return lambda i: (
            refPositions[i] >= windowStart and
            refPositions[i] < windowEnd and
            readString[i] != gapBase)

    assert alnLen == len(refPositions)

    for site, refBase, altBases in variantSites:
        windowStart = max(site - radius, 0)
        windowEnd = min(site + radius + 1, len(reference))

        if not read.spansReferenceRange(windowStart, windowEnd):
            resultRow.append("")
            continue

        isValidBase = baseInWindow(windowStart, windowEnd)

        refWindow = reference[windowStart:windowEnd]
        readWindow = "".join(readString[i] for i in xrange(alnLen) if isValidBase(i))

        if len(readWindow) > 0:
            # Adaption of original script by using Biopython for pairwise alignment
            aln = pairwise2.align.globalxx(refWindow, readWindow)
            # Identify highest scoring alignment
            scores = np.asarray([alignment[2] for alignment in aln])
            lengths = np.asarray([alignment[4] for alignment in aln])
            best_aln = np.where(lengths == lengths[scores == np.max(scores)].min())
            align_ref, align_qry, align_score, align_begin, align_end = aln[best_aln[0][0]]
            bases = returnBase(align_ref, align_qry, site, windowStart, radius)
            ambigBase = baseDict.get(bases, "X")
            resultRow.append(ambigBase)

    return resultRow

def ccsName(read):
    ''' Returns a read identifier for each name
    '''
    return "{0:s}/{1:d}/ccs".format(read.movieInfo.Name, read.HoleNumber)

def returnBase(align_ref, align_qry, site, windowStart, radius):
    ''' Returns the base of the read at the variant position
    '''
    gap_pos_ref = [m.start() for m in re.finditer("-", align_ref)]
    correction = 0
    for match in gap_pos_ref:
        if match <= radius+1:
            correction += 1
    return align_qry[site - windowStart + correction]

def parseArguments(args):
    ''' Argument parser based on argparse
    '''
    parser = argparse.ArgumentParser(description="Call which read support which variants " +
                                     "provided in the variant list")
    parser.add_argument("--radius", type=int, default=10, help="local realignment radius [10]")
    parser.add_argument("--cmp", required=True, help="input cmp.h5 file")
    parser.add_argument("--variants", required=True, help="input minor_variants.csv file")
    parser.add_argument("--reference", required=True, help="input reference FastA file")
    parser.add_argument("--output", required=True, help="csv file")
    return parser.parse_args(args)

if __name__ == "__main__":
    main()
