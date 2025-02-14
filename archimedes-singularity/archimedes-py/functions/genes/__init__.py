import os.path
import csv

with open(os.path.join( os.path.dirname(__file__), 'GRCm38.ensembl93.alias.txt')) as mouse_tsv:
    GRCm38_ensembl93 = [ row for row in csv.reader(mouse_tsv, delimiter='\t') ]
