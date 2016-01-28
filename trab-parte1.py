# -*- coding: utf-8 -*-
"""
Created on Fri Jan 1 15:42:16 2016

@author: Dany
"""
#from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO
from Web import html
from ProteinValidationTable import ProteinValidationTable
import os.path
import re
import os
#from Bio.Align.Applications import ClustalwCommandline
from fastaBlast import trata_fastas
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo
from bioservices import UniProt
from Bio.Align.Applications import ClustalwCommandline
'''
treponema palidum - Sífilis
Genes: locus tag : TP_RS01795 - TP_RS02340)
Zona do Genoma: 381901 - 501800
'''


pt = ProteinValidationTable()
site = html()
fastaBl=trata_fastas()
protein_id=""
listaGenes=[]
record=[]
listaProteinas=[]
alinhamentos=[]
ec_lista=[]
u=UniProt()

# Acede a UniProt pelo nome do Gene    
def uniprot_get2(gene):
    res=u.search("geneid:"+gene,frmt='txt')
    return res


# funçao criada pelo grupo do Sequeira
def make_kegg(ec_list):
    from Bio.KEGG import REST
    for ec in ec_list:
        print("Tratar de kegg...")
        try:
            keggname = "Kegg\\"+ec+".txt"     
            if not os.path.isfile(keggname):
                request = REST.kegg_get(ec)
                open(keggname, 'wb').write(request.read())
                print("Kegg SUCCESS!!!")
        except:
                print('kegg request failed or file already exists')
            




def parseUniprot(protein_id):
    dic = {}
    
    protname = "UniProt\\" + protein_id + ".txt" 
    with open(protname) as f:
        lines = f.read().splitlines()
    j=0
    geneName=""
    subcellarLocation = ""
    function = ""
    naoFazerSuce = False
    aFazerSubce = False
    aFazerFun = False
    geneOnt = ""
    ecNumber = ""
    nomeProteina = ""
    grauRevisao = ""
    nAA = ""
    uniprotID = ""
    keyWords = ""
    mycomments = ""
    while j < len(lines):
        lines[j] = re.sub( '\s+', ' ', lines[j])
        
        if lines[j].startswith("ID"):
            aux = lines[j].split(";")
            #nome proteina
            nomeProteina = aux[0].split(' ')[1]
            #Grau de Revisão
            grauRevisao = aux[0].split(' ')[2]
            nAA = aux[1][1:].replace(".","")
        #uniprotID 
            
        if "EC=" in lines[j]:
            ecNumber=lines[j].split("EC=")[1].split("{")[0]
            while(ecNumber.endswith(";") or ecNumber.endswith("-") or ecNumber.endswith(" ") or ecNumber.endswith(".")):
            #if ecNumber.endswith(";") or ecNumber.endswith("-") or ecNumber.endswith(" "):
                ecNumber=ecNumber[:-1]
    
        if lines[j].startswith("AC"):
            uniprotID = lines[j].split(";")[0][3:]
        
        if lines[j].startswith("DE"):
            geneName = lines[j].split("DE")[1][4:]
        #localizacao celular
        if lines[j].startswith("CC"):
            if "----" in lines[j]:
                naoFazerSuce = True
            if not naoFazerSuce:
                if "-!-" in lines[j]:
                    lines[j] = lines[j][4:]
                    if "SUBCELLULAR LOCATION" in lines[j]:
                        aFazerSubce = True
                        aFazerFun = False
                        lines[j] = lines[j][22:]
                    elif "FUNCTION" in lines[j]:
                        aFazerFun = True
                        aFazerSubce = False
                        lines[j] = lines[j][10:]
                    else:
                        aFazerSubce = False
                        aFazerFun = False
                        
            if aFazerSubce:
                subcellarLocation += lines[j][3:] + "<br>"
            if aFazerFun:
                function += lines[j][3:] + "\n"
        #GeneOntology
        if lines[j].startswith("DR GO"):
            geneOnt += lines[j][7:] + "<br>"
                
        #keywords  
        if lines[j].startswith("KW"):
            keyWords += lines[j][3:] + "\n"
                
        
                
        #print(str(lines[j]))
        if j>=0 and lines[j].startswith("//"):
            break
        j+=1
        
    if nomeProteina == "":
        nomeProteina = "NA"        
    dic['nomeProteina'] = nomeProteina
    
    if geneName == "":
        nomeProteina = "NA"        
    dic['geneName'] = geneName
    
    if grauRevisao == "":
        grauRevisao = "NA"  
    dic['grauRevisao'] = grauRevisao
    
    if nAA == "":
        nAA = "NA"  
    dic['nAA'] = nAA
    
    if uniprotID == "":
        uniprotID = "NA"  
    dic['uniprotID'] = uniprotID
    
    dic['accessionN'] = protein_id
    
    if subcellarLocation == "":
        subcellarLocation = "NA"
    else:
        subcellarLocation = subcellarLocation[:-1]
        
    dic['subcellarLocation'] = subcellarLocation
    
    if function == "":
        function = "NA"
    else:
        function = function[:-1]
        
    dic['function'] = function
    
    
    if keyWords == "":
        keyWords = "NA"
    else:
        keyWords = keyWords[:-1]
        
    dic['keyWords'] = keyWords
    
    if geneOnt == "":
        geneOnt = "NA"
    else:
        geneOnt = geneOnt[:-1]
        
    dic['geneOnt'] = geneOnt
    
    
    global ec_lista
    if ecNumber == "":
        ecNumber = "NA"
    else:
        ec_lista.append(ecNumber)
        dic['accessionN'] = protein_id+"KEGG"
   
    dic['ecNumber'] = ecNumber
    
    
    
    desc = "Comments\\" + protein_id + ".txt" 
    with open(desc) as f:
        for line in f:
            mycomments += line + "<br>"
            
    if mycomments == "":
        mycomments = "NA"
    else:
        mycomments = mycomments[:-1]
    dic['mycomments'] = mycomments
    
    return dic
    

  

def filtra_locus(tag):
        global listaProteinas        
        prot=tag.split('S')[1]
        listaProteinas.append(prot)
        num=int(prot)
        #print("numero de loc:"+str(num))
        return num


def createFasta(name,seq):
    print("Adicionar fasta ao ficheiro fastas.faa...") 
    single_Fasta= ">%s\n%s\n" % (name,seq)      
    return single_Fasta




def create_muscle():
    print("Tratar de muscle aln...")
    
    if not os.path.isfile("C:\\Users\Pedro\Desktop\laboratorios\TrabalhoLBFinal\myfasta.aln"):
        muscle_exe = r"C:\Users\Pedro\Desktop\laboratorios\TrabalhoLBFinal\muscle.exe"
        infile1 = r"C:\Users\Pedro\Desktop\laboratorios\TrabalhoLBFinal\fasta_output.faa"
        outfile1 = r"C:\Users\Pedro\Desktop\laboratorios\TrabalhoLBFinal\myfasta.aln"
        cline1 = MuscleCommandline(muscle_exe, input=infile1, out = outfile1, clwstrict=True )
        cline1()

 



def Arvore():
    print("Tratar da arvore...")
    align = AlignIO.read("myfasta.aln", "clustal")
    tree=Phylo.read("backs\clustalw.dnd", "newick")
    tree.rooted = True
    Phylo.draw(tree)
    
 
    
def run():   
    print("A arrancar...")
    fastaFile_name="fasta_output.faa"
    #global f
    f=open(fastaFile_name,"w")
    global record
    record= SeqIO.read("sifilis.gb", "genbank")
    global count
    count=0
   

    
    
    for i in range(len(record.features)):
       # print("-> " + str(i))
       
       #ja esta cortado o ficheiro gb       
       locus=record.features[i].qualifiers.get('locus_tag')
       if locus!=None:
            numlocus=filtra_locus(''.join(locus))    
        
        #em comentario temos a possibilidade de filtrar ao locus_tag
       if record.features[i].type=="CDS": #and numlocus>=1795 and numlocus<=2340:      
    
            lista = record.features[i].qualifiers
            
           
            loc=record.features[i].location
            luco=str(loc).split('(')[1][:-1]
            lista.update({'strand':luco})
                
            if 'product' in lista:
                
                    if 'exception' in lista:
                        #print("Don't have sequence")
                        #TRATAR DISTO DEPOIS
                        None
                    else:
                    
                        cds = ''.join(lista['translation'])
                        protein_id = ''.join(record.features[i].qualifiers['protein_id'])
                        fasta_out=createFasta(protein_id[:-2],cds)
                                              
                        f.write(fasta_out)
                        site.dicionario[protein_id]=lista
                        geneID=lista['db_xref'][1].split(':')[1]
                        listaGenes.append(geneID+"\n")
                        
                        protname = "UniProt\\" + protein_id + ".txt"     
                        
                        print("id da proteina:"+ protein_id)                        
                        count+=1                        
                        
                        if not os.path.isfile(protname):
                            save = open(protname, "w")
                            pesquisaUniprot=uniprot_get2(geneID)
                            save.write(str(pesquisaUniprot))
                            save.close()
                        else:
                            dic = parseUniprot(protein_id)
                            lista.update(dic)
                                                        
                        if ''.join(lista['product']) == "hypothetical protein":
                    
                            fname = "Blasts\\Hypothetical\\" + protein_id + ".xml"
                                                                         
                            if not os.path.isfile(fname):
                                print("Tratar de Blast em falta...") 
                                result_handle = NCBIWWW.qblast("blastp", "swissprot", cds, record.format("fasta"), hitlist_size=20)
                                save_file = open("Blasts\\Hypothetical\\" + protein_id + ".xml", "w")
                                save_file.write(result_handle.read()) 
                                save_file.close()
                                result_handle.close()
                                
                            else:
                                #ficheiro existe tratar o blast
                                None
                        else:
                           
                        
                            fname = "Blasts\\Not_Hypothetical\\" + protein_id + ".xml"
                        
                            if not os.path.isfile(fname):   
                                
                                print("Tratar de Blast em falta...") 
                                result_handle = NCBIWWW.qblast("blastp", "swissprot", cds, record.format("fasta"), hitlist_size=20)
                                save_file = open("Blasts\\Not_Hypothetical\\" + protein_id + ".xml", "w")
                                save_file.write(result_handle.read()) 
                                save_file.close() 
                                result_handle.close()
                            else:
                                #ficheiro existe tratar o blast
                                None
                        
            else:
                print("Not Found 'product' in qualifiers")
                        
       elif record.features[i].type=="gene":
            None
       else: 
            #TRATAR DISTO DEPOIS
            #print("Sequencia não curada")
            None
    f.close()
    make_kegg(ec_lista)
    create_muscle()
    Arvore()
    #fastaBl.blasta()
    
################   MAIN   ######################################
 
   
run()

# CRIAÇAO DO SITE
site.fazTudo()






    
