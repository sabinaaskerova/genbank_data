from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, BeforePosition, AfterPosition, ExactPosition
from Bio import Entrez, SeqIO
Entrez.email = "thomas.fischer5@etu.unistra.fr"
import http.client
import datetime
import os
import threading

class SequenceHandler:
    def __init__(self, record):
        self.record = record
        self.organism_name = self.get_organism_name_and_NC()[0]
        self.NC_number = self.get_organism_name_and_NC()[1]
    
    def is_sequence_all_N(self):
        return all(nucleotide == 'N' for nucleotide in self.record.seq) # time-consuming but needed to reject

    def get_organism_name_and_NC(self):
        try:
            organism_name = self.record.annotations['organism'].replace(' ', '_')
            NC_number = self.record.id
            return organism_name, NC_number
        except KeyError as e:
            print(f"Error getting organism name and NC number: {e}")
            return None, None

    def get_sequence(self, location):
        try:
            start = int(location.start)
            end = int(location.end)
            if start >= end:
                return None
            return self.record.seq[start:end]
        except AttributeError as e:
            print(f"Error getting sequence: {e}")
            return None

    def get_sequence_complement(self, feature):
        try:
            sequence = ""
            for part in feature.location.parts:
                sequence += self.record.seq[part.start:part.end].reverse_complement()
            return sequence
        except AttributeError as e:
            print(f"Error getting sequence complement: {e}")
            return None
    
class TestSequence:
    def __init__(self, record):
        self.record = record
    
    def check_codon(sequence):
        valid_codons = ['ttg', 'atg', 'ctg', 'gtg', 'tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg', 'att', 'atc', 'ata']
        if str(sequence[:3]).lower() not in valid_codons: # Check if the sequence starts with a valid codon
            return False
        return True

    def test_origin(handle):
        try:
            for record in SeqIO.parse(handle, "genbank"):
                if record.seq:
                    return True
                else:
                    return False
        except (ValueError, IOError, http.client.IncompleteRead): # file parsing or network-related errors
            return False
        
            
    def is_complement_join(feature):
        try:
            if isinstance(feature.location, CompoundLocation) and feature.location.operator == 'join':
                for part in feature.location.parts:
                    if isinstance(part, FeatureLocation) and part.strand == -1:
                        return True
            return False
        except Exception as e:
            print(f"An error occurred in is_complement_join: {e}")
            return False

    def test_feature_bounds(feature):
        try:
            if isinstance(feature.location, CompoundLocation):
                parts = feature.location.parts

                for i in range(len(parts) - 1):
                    last_superior_bound = int(parts[i].end)
                    next_inferior_bound = int(parts[i+1].start)
                    if last_superior_bound >= next_inferior_bound:
                        # print(f"Last superior bound {last_superior_bound} is not lower than next inferior bound {next_inferior_bound}.")
                        return False
                
                for part in parts:  # inferiority for each pair of bounds
                    inferior_bound = int(part.start)
                    superior_bound = int(part.end)
                    if inferior_bound >= superior_bound: 
                        # print(f"Inferior bound {inferior_bound} is not lower than superior bound {superior_bound}.")
                        return False
                return True
                    
            elif isinstance(feature.location.start, ExactPosition) and isinstance(feature.location.end, ExactPosition):
                lower_bound = int(feature.location.start)
                upper_bound = int(feature.location.end)

                if lower_bound < upper_bound: # valid bounds
                    return True
                else:
                    # print(f"Bounds are invalid: {lower_bound} > {upper_bound}.")
                    return False

            else: # syntax error or no bounds or not ingerers
                return False
        except Exception as e:
            print(f"An error occurred in test_feature_bounds: {e}")
            return False


class FileHandler:
    def __init__(self, feature_type, organism_name, NC_number):
        self.feature_type = feature_type
        self.organism_name = organism_name
        self.NC_number = NC_number
        self.index_lock = threading.Lock()  # Create a threading lock
    
    # pour chaque join, il y a une suite d'introns/exons qui sont compte par index
    def write_introns_exons(self, feature_type, choice,record, feature, path):
        last_end = None
        index = 1
        print(f"{path}/{self.feature_type}_{self.organism_name}_{self.NC_number}.txt")
        try:
            for part in feature.location.parts:
                if choice == 'exon' or (last_end is not None and choice == 'intron'):
                    sequence = record.seq[part.start:part.end] if choice == 'exon' else record.seq[last_end:part.start]
                    if feature.location.operator == 'complement':
                        sequence = sequence.reverse_complement()
                    if TestSequence.check_codon(sequence):
                        with open(f"{path}/{feature_type}_{self.organism_name}_{self.NC_number}.txt", "a") as out_file:
                            location_str = ','.join([f"{p.start+1}..{p.end}" for p in feature.location.parts])
                            if isinstance(feature.location, CompoundLocation) and feature.location.strand == -1 and feature.location.operator == 'join':
                                out_file.write(f"{feature_type} {self.organism_name} {self.NC_number}: complement(join({location_str})) {choice.capitalize()} {index}\n{sequence}\n")
                            else:
                                out_file.write(f"{feature_type} {self.organism_name} {self.NC_number}: join({location_str}) {choice.capitalize()} {index}\n{sequence}\n")
                        
                        # index += 1
                        with self.index_lock:
                            index += 1
                last_end = part.end
        except Exception as e:
            print(f"An error occurred while writing introns/exons: {e}")


    def write_to_file(self, sequence, location, path):
        print(f"{path}/{self.feature_type}_{self.organism_name}_{self.NC_number}.txt")
        try:
            with open(f"{path}/{self.feature_type}_{self.organism_name}_{self.NC_number}.txt", "a") as out_file:
                if isinstance(location, CompoundLocation) and location.strand == -1:
                    location_str = f"{location.start+1}..{location.end}"
                    out_file.write(f"{self.feature_type} {self.organism_name} {self.NC_number}: complement({location_str})\n")
                else:
                    location_str = f"{location.start+1}..{location.end}"
                    out_file.write(f"{self.feature_type} {self.organism_name} {self.NC_number}: {location_str}\n")

                out_file.write(str(sequence))
                out_file.write("\n")
        except Exception as e:
            print(f"An error occurred while writing to file: {e}")
            
    def needUpdate(self, record, path):
        # Vérification de l'existance d'un résultat pour cette région
        filename = \
            "_".join([self.feature_type, self.organism_name, self.NC_number])
        path = path + "/" + filename + ".txt"
        if not os.path.exists(path): # Pas de fichier pour cette région, on fait le parsing
            return True

        date = record.annotations["date"]
        try:
            '''a day (%d), a three-letter month (%b), 
            and a four-digit year (%Y), 
            separated by hyphens.'''
            date_obj = datetime.datetime.strptime(date, "%d-%b-%Y")
            file_last_mod = \
                datetime.datetime.fromtimestamp(os.path.getmtime(path))
            if date_obj > file_last_mod: # Notre séquence est plus récente que notre fichier => màj
                return True
            else: # Fichier local à jour
                return False
        except ValueError: # Date de dernière mod inaccessible, màj dans le doute
            return True


class RegionExtractor:
    def __init__(self, regions):
        self.regions = regions

            
    def process_record(self, record, path):
        sequence_handler = SequenceHandler(record)
        organism_name, NC_number = sequence_handler.get_organism_name_and_NC()
        try:
            if sequence_handler.is_sequence_all_N(): # on rejette si la sequence est entierement composee de N
                print("The entire sequence is 'N'")
            else:
                for feature in record.features: # on parcourt tous les types de region (e.g. CDS et/ou intron)
                    if feature.type in self.regions: # tester si le type de region correspond
                        # print(feature.type)
                        file_handler = FileHandler(feature.type, organism_name,NC_number)
                        if not file_handler.needUpdate(record, path):
                            # Le fichier est déjà à jour, on passe à la suite
                            print("File doesn't need update")
                            continue
                        if feature is None or feature.location is None: # on rejette si les features sont invalides
                            continue 
                        if not TestSequence.test_feature_bounds(feature): # Test syntaxe des paires de bornes
                            continue
                            
                        if isinstance(feature.location, CompoundLocation) and feature.location.operator == 'complement':
                            sequence = sequence_handler.get_sequence_complement(record, feature)
                            if TestSequence.check_codon(sequence):
                                file_handler.write_to_file(sequence, feature.location, path=path)

                        elif isinstance(feature.location, CompoundLocation) and feature.location.operator == 'join':
                            if feature.type == "CDS": #exons
                                file_handler.write_introns_exons(feature.type, "exon", record, feature,path=path)
                                if "intron" in self.regions:
                                    file_handler.write_introns_exons("intron", "intron", record, feature,path=path)

                        elif isinstance(feature.location, CompoundLocation) and TestSequence.is_complement_join(feature):
                            if feature.type == "CDS":                   
                                file_handler.write_introns_exons(feature.type, "exon", record, feature,path=path)
                                if "intron" in self.regions:
                                    file_handler.write_introns_exons("intron", "intron", record, feature,path=path)
                        
                        else:
                            
                            for location in feature.location.parts:
                                sequence = sequence_handler.get_sequence(location)
                                
                                if sequence is not None and TestSequence.check_codon(sequence):
                                    file_handler.write_to_file(sequence,location,path=path)     
                    else: # on rejette si le type de region ne correspond pas
                        continue # accelere le traitement
        except Exception as e:
            print(f"An error occurred: {e}")