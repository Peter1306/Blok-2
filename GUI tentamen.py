import tkinter
from tkinter import *

import re
import matplotlib.pyplot as plt
import numpy as np
import tkinter
from tkinter import *

def editfile():
    """ Deze functie scheid de headers van de sequenties """
    plant = open("TAIR10_pep_20101214.fa", "r")
    headers = []
    seqs = []
    seq = ""
    for line in plant:
        if ">" in line:
            if seq != "":
                seqs.append(seq)
                seq = ""
            headers.append(line.replace("\n", ""))
        else:
            seq += line.replace("\n", "")
    seqs.append(seq)
    return headers, seqs

def regex(headers,seqs):
    """Deze functie zoekt naar een match in de sequentie en koppelt deze aan een header en zet die in de sequentie"""
    findaccess = []
    for i in range(len(seqs)):
        match = re.search("[LIVMFYC].[HY].D[LIVMFY]K..N[LIVMFYCT]{3}", str(seqs[i]))
        if match:
            findaccess.append(headers[i])
    #print(findaccess)
    return findaccess

def getaccess(findaccess):
    """Deze functie split de gevonden headers em zorgt ervoor dat alleen de accessiecode gelezen wordt"""
    accesscodes = []
    for a in findaccess:
        access = a.split(" ")
        accesscodes.append(access[0][1:])
    print(accesscodes)
    return accesscodes

def opengfffile():
    """Deze functie opent het GFF3 bestand en split deze in kolommen en zet deze in een lijst"""
    with open("TAIR10_GFF3_genes.gff") as f:
        gff3genes = []
        for line in f:
            chrinfo = line.split("\t")
            gff3genes.append(gff(chrinfo))
    return gff3genes


class gff:

    def __init__(self, chrinfo):

        self.seqid = chrinfo[0]
        self.source = chrinfo[1]
        self.type = chrinfo[2]
        self.start = int(chrinfo[3])
        self.end = int(chrinfo[4])
        self.score = chrinfo[5]
        self.strand = chrinfo[6]
        self.phase = chrinfo[7]
        self.attributes = chrinfo[8]
        self.id = None
        self.parent = None

        splitattributes = chrinfo[8].rstrip().split(";")

        for i in splitattributes:
            if "Parent=" in i:
                self.parent = i[7:].split(",")
            if "ID=" in i:
                self.id = i[3:]

    def getseqid(self):
        return self.seqid

    def getsource(self):
        return self.source

    def gettype(self):
        return self.type

    def getstart(self):
        return self.start

    def getend(self):
        return self.end

    def getscore(self):
        return self.score

    def getstrand(self):
        return self.strand

    def getphase(self):
        return self.phase

    def getattributes(self):
        return self.attributes

    def getid(self):
        return self.id

    def getparent(self):
        return self.parent

def GUI(accesscodes):
    plant = Tk()
    plant.title("Information on Genes with Serine/Threonine protein kinase active site")
    plant.geometry("600x600")

    access = StringVar()
    access.set("Choose you Access-Code")

    dropdown = OptionMenu(plant, access, accesscodes)
    dropdown.pack()

    def show():
        myplant = Label(plant, text="Amount of Exons:""\n""Total length:""\n""Chromosome nr:").pack()

    infobutton = Button(plant, text="Information", command=show).pack()

    plant.mainloop()

def main():
    headers,seqs = editfile()
    findaccess = regex(headers,seqs)
    accesscodes = getaccess(findaccess)
    gff3genes = opengfffile()
    GUI(accesscodes)


main()