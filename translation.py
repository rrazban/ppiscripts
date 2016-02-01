#!/usr/bin/python

#RNA to AA

import ppisettings

rnataamap = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",	#Flatworm Mitochondrial
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"W", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"N", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"S", "AGG":"S",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
aminoacids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

aathydromap = {"C":0.255409, "M":0.281532, 
	"F":0.332123, "I":0.309393, "L":0.333128, "V":0.275519,
	"W":0.275334, "Y":0.249067, "A":0.196047, "G":0.170884,
	"T":0.177372, "S":0.157779, "N":0.152413, "Q":0.161195,
	"D":0.143849, "E":0.144312, "H":0.197540, "R":0.163297,
	"K":0.125262, "P":0.169253}

aatnummap = {"C":0, "M":1, "F":2, "I":3, "L":4, "V":5, 
		"W":6, "Y":7, "A":8, "G":9, "T":10, "S":11,
		"N":12, "Q":13, "D":14, "E":15, "H":16, 
		"R":17, "K":18, "P":19}

#assume sequence multiple of three and all base pairs are codons
def rnatoaa(rna):
	aas=''
	for i in range(ppisettings.aalen):
		codon=rna[3*i:3*i+3]
		aas+=rnataamap[codon]
	return aas
