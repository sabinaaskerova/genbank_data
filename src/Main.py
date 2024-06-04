from Fouille import Fouille
from Tree import Root
from Log import Log

cache_name = "logs/aborted-parsing.txt"

# Initialise l'arborescence
arbre = Root()

# Retrouve les séquences de tous les organismes de l'arborescence

regions = ['CDS', 'intron', 'telomere', 'centromere','mobile_element','ncRNA','rRNA','tRNA','3\'UTR','5\'UTR']
fouille = Fouille(regions)

print("Peuplement de l'arborescence...")

abort_log = Log(cache_name)
fouille.peuplement(arbre, abort_log, abort_log.gen_list())
abort_log.delete()

