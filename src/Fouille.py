from Bio import Entrez, SeqIO
from TestsExtraction import RegionExtractor
from TestsExtraction import TestSequence
from concurrent.futures import ThreadPoolExecutor
import threading
from os_compatibility import format_path
import http.client  # For handling http.client.IncompleteRead
from urllib.error import URLError

class Fouille:
    def __init__(self, regions):
        self.regions = regions
        self.total_organisms = None
        self.processed_organisms = 0
        self.total_ncs = 0
        self.processed_ncs = 0
        self.lock = threading.Lock()  # Create a threading lock

        filename = format_path("Results/overview.txt")
        with open(filename, 'r') as file:
            self.total_organisms = sum(1 for _ in file) - 1

        print("Calculating total number of NCS...")
        files = ["Eukaryota.ids", "Bacteria.ids", "Archaea.ids", "Viruses.ids"]
        for file in files:
            filename = format_path("Results/"+file)
            with open(filename, 'r') as file:
                self.total_ncs += sum(1 for _ in file) - 1


    
    def search(self, organisme, db = "nucleotide"):

        Entrez.email = "thomas.fischer5@etu.unistra.fr"

        if organisme == "":
            return []
        
        requete = "(" + organisme + "[Organism] AND NC_000001:NC_999999[ACCN]"
        if self.regions != []:
            requete += " AND (" + "[FKEY] OR ".join(self.regions) + "[FKEY]" # filtrer par les regions qui nous interessent
        requete += ")"
        
        handle = Entrez.esearch(
            db = db, 
            term=requete,
            usehistory='y', 
            idtype="acc"
        )
        resultat = Entrez.read(handle)
        handle.close()

        return resultat["IdList"]

    def fetchFromID(self, ids, db="nuccore", rettype="gbwithparts"):
        Entrez.email = "thomas.fischer5@etu.unistra.fr"

        handle = None

        finBool = False
        while not finBool:
            try:
                handle = Entrez.efetch(
                    db=db, 
                    id=ids, 
                    rettype=rettype, 
                    retmode="text", 
                )
            except:
                if(handle!=None):
                    try:
                        handle.close()
                    except:
                        pass
                return None, 0
            else:
                finBool = True

        return handle
    
    def fetchAll(self, organisme):
        ids = ','.join(self.search(organisme))
        lg = 0 if ids == "" else len(ids.split(','))
        # print(f"{organisme}: {lg} found.")
        if (len(ids) > 0):
            return self.fetchFromID(ids), 0
        else:
            return [], 1
        
    def process_record(self, extractor, record, path):
        extractor.process_record(record, path)
        with self.lock:
            self.processed_ncs += 1
        
        print(f"Processed {self.processed_ncs} records.")
        # print(f"Percentage: {self.processed_ncs / self.total_ncs * 100:.2f}%")
    
    def fetch_and_process(self, organisme, path, abort_log):
        try:
            handle, ret = self.fetchAll(organisme)
            if not handle:
                return

            extractor = RegionExtractor(self.regions)
            if not TestSequence.test_origin(handle):
                print("No sequence found.")
                abort_log.add()
                return

            # print("Fetching and processing records.")
            with ThreadPoolExecutor() as record_executor:
                for record in SeqIO.parse(handle, "genbank"):
                    record_executor.submit(
                        self.process_record, 
                        extractor, 
                        record, 
                        path)

        except (http.client.IncompleteRead, IOError, ValueError, Entrez.ParserError, URLError, ConnectionError, TimeoutError) as e:
            print(f"An error occurred: {e}")
            abort_log.add()
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            abort_log.add()
        
        return

    def peuplement(self, tree, abort_log, ignore=[], progress_callback=None):
        with ThreadPoolExecutor() as executor:
            futures = []
            for child in tree.children :
                if child.children == []:
                    if len(ignore) > 0:
                        print(f"{child.data} already parsed in aborted run.")
                        del ignore[0]
                    else:
                        future = executor.submit(
                            self.fetch_and_process, 
                            child.data, 
                            child.path, 
                            abort_log)
                        futures.append(future)
                else:
                    self.peuplement(child, abort_log, ignore)
                self.processed_organisms += 1
                
                if progress_callback: # send the progress to the PyQT5 interface 
                    progress_callback(
                        self.total_organisms, 
                        self.processed_organisms, 
                        self.total_ncs, 
                        self.processed_ncs)

            # Attendre l'achèvement de toutes les tâches
            for future in futures:
                future.result()
                
    def get_total_organisms(self):
        return self.total_organisms

    def get_processed_organisms(self):
        return self.processed_organisms

    def get_total_ncs(self):
        return self.total_ncs
    
    def get_processed_ncs(self):
        return self.processed_ncs