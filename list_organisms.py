#!/usr/bin/env python
import sys
import csv
from operator import itemgetter
import warnings

FAST_SEARCH = True

def download_org_genome(organism):
    from Bio import SeqIO
    # Entrez.email = ""
    # first get the genbank id for the organism genome
    handle = Entrez.esearch(db="nucleotide", term="%s[Organism]" % organism,
                            rettype="gb", retmode="xml")
    r = Entrez.read(handle)
    handle.close()
    genome_id = r['IdList'][0]
    # now download the genome in fasta format
    handle = Entrez.efetch(db="nucleotide", id=genome_id,
                           rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    genome_file = organism.replace(' ', '_') + '.fasta'
    SeqIO.write(seq_record, genome_file, 'fasta')

    return genome_file


def get_organism(gb_id):
    from Bio import Entrez
    Entrez.email = "zagordi.osvaldo@virology.uzh.ch"
    handle = Entrez.efetch(db="nucleotide", id=str(gb_id),
                           rettype="gb", retmode="xml")
    r = Entrez.read(handle)
    handle.close()
    return r[0]['GBSeq_organism']


class DBhit:
    '''Parses the single hits and exposes them'''
    def __init__(self, full_line):
        self.name = full_line[0]
        self.hit = full_line[8].split('_')[1]
        self.hit_name = full_line[8]
        self.cov = float(full_line[3])
        self.iden = float(full_line[7][:-1])
        try:
            self.strand = full_line[6]
        except:
            self.strand = 'unknown'
    # self.start, self.stop = full_line[9:11]

    def add_organism(self, org):
        self.organism = org

class DBmatch:
    '''Parses all hits to the database'''
    def __init__(self, results_file):
        import csv
        self.matches = {}
        csvfile = open(results_file, 'rb')
        for line in csv.reader(csvfile, delimiter='\t'):
            self.matches[line[0]] = DBhit(line)
    

    def out(self):
        kh = self.matches.keys()[0]
        for k in vars(self.matches[kh]):
            print k, self.matches[kh].__dict__[k]


    def get_organisms(self):
        # group the genbank ids so to do less interrogations
        count_ids = {}
        for m, v in self.matches.items():
            count_ids[v.hit] = count_ids.get(v.hit, 0) + 1
        ms = sorted(count_ids.iteritems(), key=itemgetter(1), reverse=True)

        print >> sys.stderr, "Doing %d interrogations to Entrez" % len(ms)
        if len(ms) > 50:
            warnings.warn("This might take a while")
        orgs_count = {}  # only used to count
        orgs_dict = {}  # used to match genbank id with organism name
        for s in ms:
            org = get_organism(s[0])
            print >> sys.stderr, s[0], org
            orgs_dict[s[0]] = org
            orgs_count[org] = orgs_count.get(org, 0) + s[1]

        self.orgs_list = sorted(orgs_count.iteritems(), key=itemgetter(1),
                                reverse=True)

        for m, v in self.matches.items():
            try:
                v.add_organism(orgs_dict[v.hit])
            except KeyError:
                v.add_organism('single_read')


    def list_organisms(self):
        tv = sum([p[1] for p in self.orgs_list])
        print 'organism,reads,percent'
        for org_hit in self.orgs_list:
            print '%s,%d,%4.1f' % (org_hit[0], org_hit[1],
                                   100*(float(org_hit[1]))/tv)


db = DBmatch(sys.argv[1])
#db.out()
db.get_organisms()
db.list_organisms()
