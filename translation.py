#!/usr/bin/python

#RNA to AA

import ppisettings

mapaa = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",	#Flatworm Mitochondrial
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
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
aminoacids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

#assume sequence multiple of three and all base pairs are codons
def rnatoaa(rna):
	aas=''
	for i in range(ppisettings.aalen):
		codon=rna[3*i:3*i+3]
		aas+=mapaa[codon]
	return aas
