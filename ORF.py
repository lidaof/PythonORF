#/usr/bin/env python
#Author: lidaof@gmail.com
#Date: 12-18-2010
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ORF:
    '''ORF class'''
    def __init__(self, fa, type='DNA'):
        self.fa = fa
        self.type = type
        if self.type == 'DNA':
            self.startcodon = 'ATG'
        else:
            self.startcodon = 'AUG'
        self.plus_orfs = []
        self.minus_orfs = []
        self.orf_number = 0

    def oneStrandORF(self,table,reverse_flag,index,length):
        """input a SeqRecord object
        reverse_flag: true or false, true means reverse it,
        index: start from 0 when plus, add plus when reverse"""
        recs = []
        loclis = []
        startpos = 0
        table = int(table)
        index = int(index)
        seq = self.fa.seq
        if reverse_flag:
            seq = Seq.reverse_complement(seq)
        while seq.count(self.startcodon,startpos) > 0:
            pos = seq.find(self.startcodon,startpos)
            loclis.append(pos)
            startpos = pos + len(self.startcodon)
        i = 1
        pend_d = {} #check children orfs belong to another long orf -- aviod this
        mend_d = {}
        for loc in loclis:
            alreadyexsit = True
            proseq = Seq.translate(seq[loc:], table = table, to_stop=True) 
            if len(proseq) >= int(length):
                proid = self.fa.id + '_'  + str(i+index)
                end = loc + (len(proseq) * 3)
                if reverse_flag:
                    desc = '[' + str(len(seq)-loc) + ' - ' + str(len(seq)-end+1) + ']' + ' (REVERSE SENSE) ' + ' '.join(self.fa.description.split()[1:])
                    mend = len(seq)-end+1
                    if mend not in mend_d:
                        mend_d[mend] = 1
                        alreadyexsit = False
                else:
                    desc = '[' + str(loc+1) + ' - ' + str(end) + '] ' + ' '.join(self.fa.description.split()[1:])
                    if end not in pend_d:
                        pend_d[end] = 1
                        alreadyexsit = False
                prorec = SeqRecord(proseq, id=proid, description=desc)
                #SeqIO.write(prorec,out,format)
                if not alreadyexsit:
                    recs.append(prorec)
                    i += 1
        return recs
    def getORFs(self,table,length):
        self.plus_orfs = self.oneStrandORF(table,False,0,length)
        self.minus_orfs = self.oneStrandORF(table,True,len(self.plus_orfs),length)
        self.orf_number = len(self.plus_orfs) + len(self.minus_orfs)
