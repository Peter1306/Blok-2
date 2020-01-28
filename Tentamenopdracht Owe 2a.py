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
    try:
        findaccess = []
        for i in range(len(seqs)):
            match.lower = re.search("[LIVMFYC].[HY].D[LIVMFY]K..N[LIVMFYCT]{3}", str(seqs[i]))
            if match:
                findaccess.append(headers[i])
    #print(findaccess)
        return findaccess
    except: 

def getaccess(findaccess):
    """Deze functie split de gevonden headers em zorgt ervoor dat alleen de accessiecode gelezen wordt"""
    accesscodes = []
    for a in findaccess:
        access = a.split(" ")
        accesscodes.append(access[0][1:])
    print(accesscodes)
    return accesscodes

def gfffile(accesscodes):
    """Deze functie zoekt in dit bestand alle chromosomen en accessiecodes en zet deze apart in een lijst"""
    matchcodes = []
    chr = []
    with open("TAIR10_GFF3_genes.gff") as f:
        for line in f:
            chrinfo = line.split("\t")
            chr.append(chrinfo[0])
            matchcodes.append(chrinfo[8])
    print(matchcodes)
    return matchcodes, chr

def match(matchcodes, chr, accesscodes):
    chromosomes = []
    for ac in accesscodes:
        accessmatch = re.search(ac, matchcodes)
        if accessmatch:
            chromosomes.append(chr[ac])
    print(chromosomes)
    return chromosomes

def graph(chromosomes):
    nr_of_chromosomes = []
    values = []
    for entry in chromosomes:
        if entry not in nr_of_chromosomes:
            nr_of_chromosomes.append(entry)
            values.append(1)
        elif entry in nr_of_chromosomes:
            for index in range(len(nr_of_chromosomes)):
                if entry == nr_of_chromosomes[index]:
                    values[index] += 1
    return values, nr_of_chromosomes

def main():
    headers,seqs = editfile()
    findaccess = regex(headers,seqs)
    accesscodes = getaccess(findaccess)
    chr, matchcodes = gfffile(accesscodes)
    chromosomes = match(matchcodes, chr, accesscodes)
    values, nr_of_chromosomes = graph(chr)
    y_pos = np.arange(len(nr_of_chromosomes))
    plt.bar(y_pos, values)
    plt.xticks(y_pos, nr_of_chromosomes)
    plt.xlabel("Chromosome nr.")
    plt.ylabel("Amount of genes with Serine/Threonine kinase active site")
    plt.title('Amount of genes with Serine/Threonine kinase active site per chromosome')
    plt.show()
main()
