import re
import matplotlib.pyplot as plt
import numpy as np

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
    with open("TAIR10_GFF3_genes.gff") as f:
        gff3genes = []
        for line in f:
            chrinfo = line.split("\t")
        gff3genes.append(chrinfo)
        return chrinfo, gff3genes

class gff:

    def __init__(self,chrinfo):

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

        for i in self.attributes:
            splitattributes = chrinfo[8].split(";")
            if "Parent=" in splitattributes:
                self.parent = i[3:]
            if "ID=" in splitattributes:
                self.id = i[8:].split(",")

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

    def getmatch(self,gff3genes,accesscodes):
        matches = []
        for ac in accesscodes:
            for obj in gff3genes:
                if obj.getid == ac or obj.getparent[0] == ac:
                    matches.append(obj.getseqid)
        print(matches)
        return matches

#def graph(chromosomes):
 #   nr_of_chromosomes = []
  #  values = []
   # for entry in chromosomes:
    #    if entry not in nr_of_chromosomes:
     #       nr_of_chromosomes.append(entry)
      #      values.append(1)
       # elif entry in nr_of_chromosomes:
        #    for index in range(len(nr_of_chromosomes)):
         #       if entry == nr_of_chromosomes[index]:
          #          values[index] += 1
    #return values, nr_of_chromosomes

def main():
    headers,seqs = editfile()
    findaccess = regex(headers,seqs)
    accesscodes = getaccess(findaccess)
    chrinfo, gff3genes = opengfffile()
    matches = getmatch(gff3genes, accesscodes)
    #chr, matchcodes = gfffile(accesscodes)
    #chromosomes = match(matchcodes, chr, accesscodes)
    values, nr_of_chromosomes = graph(chr)
    y_pos = np.arange(len(nr_of_chromosomes))
    plt.bar(y_pos, values)
    plt.xticks(y_pos, nr_of_chromosomes)
    plt.xlabel("Chromosome nr.")
    plt.ylabel("Amount of genes with Serine/Threonine kinase active site")
    plt.title('Amount of genes with Serine/Threonine kinase active site per chromosome')
    plt.show()
main()
