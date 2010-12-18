#!/usr/bin/env python
#Author: lidaof@gmail.com
#Date: 12-18-2010
'''%prog -i <input nucleotide sequence file> 
input a nucleotide (DNA or RNA) sequence file, output a protein sequence file containing ORFs found'''

from Bio import SeqIO
from ORF import *

__version__ = '0.1'

def main():
    import sys
    from optparse import OptionParser
    
    parser = OptionParser(usage=__doc__, version=__version__)
    parser.add_option("-i","--input",dest="input",metavar="INPUT",help="input a nucleotide sequence (DNA/RNA) file for ORFs finding [required]")
    parser.add_option("-f","--input-format",dest="informat",help="format of input sequence file,[default: %default]")
    #parser.add_option("-intype","--input-type",dest="intype",help="type of input sequence file: DNA or RNA") -> auto detect
    parser.add_option("-o","--output",dest="output",metavar="OUTPUT",help="output the result protein sequence file,[default: STDOUT]")
    parser.add_option("-F","--output-format",dest="outformat",help="format of output sequence file,[default: %default]")
    parser.add_option("-c","--codon-table",type="int",dest="codontable",help="codon table used when translating,[default: %default]")
    parser.add_option("-s","--strand",dest="strand",help="translate on plus/minus strand,[default: %default]")
    parser.set_defaults(informat='fasta')
    #parser.set_defaults(intype='DNA')
    #parser.set_defaults(output=sys.stdout)
    parser.set_defaults(outformat='fasta')
    parser.set_defaults(codontable=1)
    parser.set_defaults(strand='both')
    (options, args) = parser.parse_args()
    if options.input == None:
       parser.error("must specify an input file,use -h to see parameters")
    informat = options.informat
    outformat = options.outformat
    table = options.codontable
    if options.output == None:
        out = sys.stdout
    else:
        out = open(options.output,'w')
    fas = SeqIO.parse(options.input, informat)
    for fa in fas:
        orfobj = ORF(fa)
        if options.strand == 'both':
            orfobj.getORFs(table)
            for rec in orfobj.plus_orfs:
                SeqIO.write(rec,out,outformat)
            for rec in orfobj.minus_orfs:
                SeqIO.write(rec,out,outformat)
        elif options.strand == 'plus':
            plus_orfs = orfobj.oneStrandORF(table,False,0)
            for rec in plus_orfs:
                SeqIO.write(rec,out,outformat)
        elif options.strand == 'minus':
            minus_orfs = orfobj.oneStrandORF(table,True,0)
            for rec in minus_orfs:
                SeqIO.write(rec,out,outformat)
    out.close()

if __name__ == "__main__":
    main()
