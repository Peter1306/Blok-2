# Auteur: Peter Bos
# Datum: 28-01-2020 - 31-01-2020
# Thematoets Python II - Course 2a 2019/2020

import matplotlib.pyplot as plt
from tkinter import *
from tkinter import Canvas


class Gff:
    """Deze functie neemt de data van het GFF3-bestand en maakt hier
    variabelen van in het object.
    De namen van de kolommen (bijv. self.seq_id) komt van:
    https://www.ensembl.org/info/website/upload/gff3.html#fields"""
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

        split_attributes = chrinfo[8].rstrip().split(";")

        for i in split_attributes:
            if "Parent=" in i:
                self.parent = i[7:].split(",")
            if "ID=" in i:
                self.id = i[3:]

    def get_seqid(self):
        return self.seqid

    def get_source(self):
        return self.source

    def get_type(self):
        return self.type

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_score(self):
        return self.score

    def get_strand(self):
        return self.strand

    def get_phase(self):
        return self.phase

    def get_attributes(self):
        return self.attributes

    def get_id(self):
        return self.id

    def get_parent(self):
        return self.parent

    def get_length(self):
        self.length = self.end - self.start
        return self.length


def edit_file():
    """ Deze functie scheid de headers van de sequenties en zet beide
    in een aparte lijst """
    try:
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
    except FileNotFoundError:
        print("Cannot find the file")
        exit()


def regex(headers, seqs):
    """Deze functie zoekt naar een match in de sequentie en koppelt
    deze aan een header en zet die in de sequentie"""
    find_access = []
    for i in range(len(seqs)):
        match = re.search("[LIVMFYC].[HY].D[LIVMFY]K..N[LIVMFYCT]{3}",
                          str(seqs[i]))
        if match:
            find_access.append(headers[i])
    return find_access


def get_access(find_access):
    """Deze functie split de gevonden headers en zorgt ervoor dat
    alleen de accessiecode gelezen wordt"""
    access_codes = []
    for a in find_access:
        access = a.split(" ")
        access_codes.append(access[0][1:])
    return access_codes


def open_gff_file():
    """Deze functie opent het GFF3 bestand en split deze in kolommen
    en zet deze in een lijst"""
    try:
        with open("TAIR10_GFF3_genes.gff") as f:
            gff3_genes = []
            for line in f:
                chr_info = line.split("\t")
                gff3_genes.append(Gff(chr_info))
        return gff3_genes
    except FileNotFoundError:
        print("Cannot find the file")
        exit()


def make_graph_data(gff3_genes, access_codes):
    """Deze functie kijkt in hoeverre de accessiecodes uit de
    fasta-file overeenkomen met de eiwitten in het GFF3 bestand.
    Het chromosoomnummer (Chr1/2/3/4/5) wordt vervolgens in een lijst
    gezet en bij elkaar op geteld."""
    chromosome_proteins = {}
    for gen in gff3_genes:
        if gen.get_type() == "protein":
            for access_code in access_codes:
                if access_code in gen.get_attributes():
                    if gen.get_seqid() not in chromosome_proteins:
                        chromosome_proteins[gen.get_seqid()] = 1
                    else:
                        chromosome_proteins[gen.get_seqid()] += 1
    return chromosome_proteins

def graph(chromosomeproteins):
    """Deze functie zet data uit de lijst chromosome_proteins in een
    grafiek en geeft ook naam aan(de labels van) de grafiek"""
    chrnr = []
    proteins = []
    for chromosome, protein in sorted(chromosomeproteins.items()):
        chrnr.append(chromosome)
        proteins.append(protein)
    plt.bar(chrnr, proteins)
    plt.xlabel("Chromosome nr.")
    plt.ylabel("Amount of genes with Serine/Threonine kinase active "
               "site")
    plt.title('Amount of genes with Serine/Threonine kinase active '
              'site per chromosome')
    plt.show()


def gui(gff3genes, accesscodes):
    """Deze functie maakt een GUI met een dropdown-menu(met
    accessiecodes) met een button die, als je er op klikt informatie
    weergeeft over de gekozen accessiecode"""
    plant = Tk()
    plant.title("Information on Genes with Serine/Threonine protein "
                "kinase active site")
    plant.geometry("600x600")

    access = StringVar()
    access.set("Choose you Access-Code")

    dropdown = OptionMenu(plant, access, *sorted(accesscodes))
    dropdown.pack()

    info = StringVar()
    myplant = Label(plant, textvariable=info).pack()

    for entry in gff3genes:
        if chromosome_number == entry.get_seqid():
            if entry.get_type() == "Chromosome" and entry.get_start()\
                    = 1:
                start = entry.get_start()
                end = entry.get_end()

    # start == 1
    # end == 100

    #Deze in de forloop van show

    for i in gff3genes:
        if access.get() in i.get_attributes():
            if i.get_type() == "gene":
                location = i.get_end()
                location_on_gene = (start/location)/100





    def show(info, gff3genes, access):
        """Deze functie genereert de informatie over de
        accessiecode. Dus de lengte, het aantal Exonen en het
        Chromosoom nummer"""
        exons = 0
        for gene in gff3genes:
            if access.get() in gene.get_attributes():
                length = gene.get_length()
                chromosome_number = gene.get_seqid()
                if "exon" in gene.get_type():
                    exons += 1
        info.set("Amount of Exons: " + str(exons) + "\n""Total length: "
                 + str(length) + "\n""Chromosome nr: " +
                 str(chromosome_number))

    infobutton = Button(plant, text="Information", command=lambda: \
    show(info, gff3genes, access)).pack()

    plant.mainloop()


def main():
    headers, seqs = edit_file()
    findaccess = regex(headers, seqs)
    accesscodes = get_access(findaccess)
    gff3_genes = open_gff_file()
    chromosomeproteins = make_graph_data(gff3_genes, accesscodes)
    graph(chromosomeproteins)
    gui(gff3_genes, accesscodes)


main()
