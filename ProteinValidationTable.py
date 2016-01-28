# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:38:45 2016

@author: DavidMiguel
"""
import csv

class ProteinValidationTable:
    fname = 'ProteinTable.txt'
    lista = {}
    def __init__(self):
        with open(self.fname) as f:
            reader = csv.reader(f, delimiter="\t")
            self.lista = list(reader)
        
    def existeProteina(self,nome):
        for i in range(len(self.lista)):
            if nome in self.lista[i]:
                return i
        return -1
        
    def getProteina(self,i):
        if i > len(self.lista):
            return {}
        return self.lista[i]
        