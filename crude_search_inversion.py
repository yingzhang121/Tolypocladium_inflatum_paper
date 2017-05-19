#!/usr/bin/evn python

import sys, argparse
from itertools import chain

def build_dict(infile):

    coords = {}

    for line in infile:
        fields = line.split()
        ref_st, ref_end, qry_st, qry_end = map(int, fields[0:4])
        ref_chr, qry_chr, qry_size, qry_strand = fields[13], fields[14], int(fields[8]), int(fields[12])
        if ref_chr not in coords: coords[ref_chr] = {}
        if qry_chr not in coords[ref_chr]: coords[ref_chr][qry_chr] = {"coords":{0:[],1:[]}, "size":qry_size, "strands":[] }
        coords[ref_chr][qry_chr]["coords"][0].append([ref_st, ref_end])
        coords[ref_chr][qry_chr]["coords"][1].append(sorted([qry_st, qry_end]))
        coords[ref_chr][qry_chr]["strands"].append(qry_strand)

    return coords

def get_parser():
    """ Parse input """
    desc = "Very aggressive clustering of homologous regions based on mummer coords file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('coords', nargs="?", type=argparse.FileType('r'), default=sys.stdin, help="mummer show-coords output" )
    parser.add_argument('output', nargs="?", type=argparse.FileType('wt'), default=sys.stdout, help="output file, default stdout")
    parser.add_argument('-d', '--dist', type=int, default=5000, help="maximal distance for clustering")
    return parser

def get_gaps( rows ):

    n = len(rows) - 1
    gaps = [ rows[i+1][0]-rows[i][1] for i in range(n) ]
    return gaps

def get_block_size( coords ):
    return [ x[1]-x[0] for x in coords ]

def _scan_reverse(gaps, center, dist):

    for i in range( 0, center ):
        idx_gap = center - 1 - i
        gap = gaps[idx_gap]
        if gap >= dist: return idx_gap+1
    return 0

def _scan_forward( gaps, center, dist ):

    n = len(gaps)
    for i in range( center, n ):
        idx_gap = i
        gap = gaps[idx_gap]
        if gap >= dist: return idx_gap+1
    return n+1

def scan( gaps, center, dist ):
    return _scan_reverse( gaps, center, dist ), _scan_forward( gaps, center, dist )

def cluster_regions( coords, dist ):

    gaps = get_gaps( coords )
    blocks = get_block_size( coords )

    max_block = max( blocks )
    center = blocks.index(max_block)

    stblock, endblock = scan( gaps, center, dist )
    return get_range( coords[stblock:endblock] )

def get_size(l):
    st, end = get_range(l)
    return end-st

def clustering( rcoords, qcoords, dist ):

    newqsize = 0
    length = get_size(qcoords)
    while True:
        refst, refend = cluster_regions( rcoords, dist )
        qryst, qryend = cluster_regions( qcoords, dist )
        qsize = qryend - qryst
        #print(dist, length, qsize)
        if qsize/length >= 0.95 or qsize == newqsize:
            return refst, refend, qryst, qryend
        else:
            dist = 2*dist
            newqsize = qsize

def determine_orientation( lcoords, strands ):
    """Majority vote to identify the orientation"""

    check1 = sum(strands)
    size_strand = [ (lcoords[i][1]-lcoords[i][0])*strands[i] for i in range(len(strands)) ]
    check2 = sum(size_strand)
    if check1 >= 0 or check2 >= 0: return 1
    else: return -1

def find_idexes( strands, match ):

    idxst = []
    idxend = []
    direction = match
    i = 0
    while i < len(strands):
        try: idx = strands[i:].index( match )
        except: break
        if match == direction: idxst.append( idx+i )
        else: idxend.append( idx+i )
        i += idx+1
        match = 0-match
    if len( idxst ) != len( idxend ):
        idxend.append( len(strands) )

    return idxst, idxend

def _processing_qchr( coords, rchr, dist, outf ):

    qry_chrs = list( coords.keys() )
    for qchr in qry_chrs:
        homo_coords = coords[qchr]["coords"]
        size = coords[qchr]["size"]
        strands = coords[qchr]["strands"]

        # determine orientation using majority vote
        direction = determine_orientation( homo_coords[1], strands )
        match = 0-direction
        if match not in strands:
            #print("no inversions between %s and %s" % (rchr, qchr), file=sys.stderr)
            continue
        idxst, idxend = find_idexes( strands, match )
        print("processing:", rchr, qchr, file=sys.stderr)

        for st, end in zip( idxst, idxend ):
            refcoords = homo_coords[0][st:end]
            qrycoords = homo_coords[1][st:end]
            refst, refend, qryst, qryend = \
                clustering( refcoords, sorted(qrycoords), dist )
            if qryend - qryst <= 1000: continue
            print("%s\t%d\t%d\t%s\t%d\t%d\t%d" % \
                (rchr, refst, refend, qchr, qryst, qryend, size), file=outf)

def processing_data( infile, dist, outf ):

    coords = build_dict( infile )
    ref_chrs = list( coords.keys() )
    for rchr in ref_chrs:
        _processing_qchr( coords[rchr], rchr, dist, outf )

def main():

    parser = get_parser()
    args = parser.parse_args()
    infile = args.coords
    dist = args.dist
    outf = args.output

    processing_data( infile, dist, outf )

if __name__=="__main__":
    main()
