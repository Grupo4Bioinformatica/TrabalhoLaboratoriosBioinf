# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 00:59:34 2016

@author: Pedro
"""
#blast_record.alignments[0].accession

import glob, os
from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO

class trata_fastas:
    
    def __init__(self):
         self.caminho="fastas\\"
         self.E_VALUE_THRESH = 0.05

    def blasta(self):
        print("A tratar de fastas individuais....")
        for file in glob.glob(self.caminho+"*.txt"):
            if not os.path.isfile(file+".xml"):            
                record=SeqIO.read(file,"fasta")
                cds=record.seq
                result_handle = NCBIWWW.qblast("blastp", "swissprot", cds, record.format("fasta"), hitlist_size=20,entrez_query="homo sapiens[organism]")
                save_file = open(""+file[:-4]+ ".xml", "w")
                save_file.write(result_handle.read()) 
                save_file.close() 
                result_handle.close()     
        return None            
    
    
    def blasta_result(self):        
        print("\n\n\n")  
        text=[]
        for file in glob.glob(self.caminho+"*.xml"):
            #print("caminho:"+file)
                results = open(file) 
                blast_record = NCBIXML.read(results)
                #print(blast_record)    
                #print(blast_record.database)
                #print(blast_record.gap_penalties)
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < self.E_VALUE_THRESH:
                            text.append('****Alignment****</br>')
                            text.append('sequence:'+alignment.title+"</br>")
                            text.append('length:'+ str(alignment.length)+"</br>")
                            text.append('e value:'+ str(hsp.expect)+"</br>")
                            text.append( hsp.query[0:75] + '...</br>' )
                            text.append( hsp.match[0:75] + '...</br>' )
                            text.append( hsp.sbjct[0:75] + '...</br>' )
        return text