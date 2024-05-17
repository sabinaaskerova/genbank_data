import requests
import csv
import os
import datetime
import shutil
from os_compatibility import format_filename, format_path
import dateparser
import pytz
class Node():
    def __init__(self, data, path, parent=None):
        self.data = data
        self.children = []
        self.parent = parent
        self.path = path #le chemin relatif du nœud dans le système de fichier, depuis le repertoire racine du projet
        self.nc_records_count = 0

    def addChild(self, obj):
        obj.parent = self 
        self.children.append(obj)
        return obj
    
class Root(Node):# Contient les méthodes de construction de l'arbre des organismes et des répertoire.
    # filename: le chemin relatif du fichier overview.txt  
    def download_file(self, url, filename): # for downloadig NC ids files
        response = requests.get(url)
        if response.status_code == 200:
            if not os.path.isdir('Results'):
                os.mkdir('Results')
            with open(os.path.join('Results', filename), 'wb') as f:
                f.write(response.content)
        else:
            print(f'Error: Unable to get {filename}')


    def __init__(self) -> None: #Initialise le nœud racine et appelle toutes les méthodes de création/mise à jour de l'arbre.
        # super().__init__('.', '/Results')
        super().__init__('.', '.')
        base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/"
        files = ["Archaea.ids", "Bacteria.ids", "Eukaryota.ids", "Viruses.ids"]

        for file in files:
            url = base_url + file
            self.download_file(url, file)

        baseurl = 'https://ftp.ncbi.nlm.nih.gov/' # Ouverture de l'overview
        path = 'genomes/GENOME_REPORTS/'
        filename = 'overview.txt'
        self.filename = "Results/" + filename
        
        if os.path.exists(self.filename): # Si overview.txt existe déjà : mise à jour 
            print("Overview already exists, looking for update...")
            self.__parseOverview(bar=True)
            self.__updateOverview(baseurl + path + filename)
        else: # Sinon initialisation
            self.__createOverview(baseurl + path + filename)
            self.__parseOverview(bar=True)

    def __updateOverview(self, url): #Met à jour l'aborescence locale à partir du fichier overview situé  à l'url en paramètre.
        r = requests.head(url)  # Obtention des dates de dernières modification 
        date_header = r.headers["last-modified"]
        print(date_header)

        url_last_mod = dateparser.parse(date_header)
        update_date_known = True
        file_last_mod = os.path.getmtime(self.filename)
        if file_last_mod is not None:
            formatted_time = datetime.datetime.fromtimestamp(file_last_mod).strftime("%a, %d %b %Y %H:%M:%S")
            file_last_mod = pytz.UTC.localize(datetime.datetime.strptime(formatted_time, "%a, %d %b %Y %H:%M:%S"))
        else:  # Handle the case where modification time couldn't be retrieved
            print("Error: Couldn't retrieve modification time for", self.filename)
            file_last_mod = datetime.datetime.now(pytz.UTC)  # Use current time as a fallback
            update_date_known = False

        if url_last_mod > file_last_mod or not update_date_known: # Si la version locale n'est plus à jour
            print("Updating overview file")
            r = requests.get(url)
            if (r.status_code == 200):
                old_path = 'Results/temp.csv' # Calcule du delta entre l'ancien et le nouvel overview
                os.rename(self.filename, old_path)
                with open(self.filename, "w") as file:
                    file.writelines(r.text)
                    file.close()
                with open(self.filename, "r") as new, \
                      open(old_path, "r") as old:
                    added = []
                    removed = []
                    new_lines = new.readlines()
                    old_lines = old.readlines()
                    for line in new_lines:
                        if line not in old_lines:
                            added.append(line.split('\t'))
                    for line in old_lines:
                        if line not in new_lines:
                            removed.append(line.split('\t'))
                    new.close()
                    old.close()
                
                for row in added: # Insertions et déletions éventuelles
                    if row != [] and len(row) > 3:
                        self.__addBranch(row)
                for row in removed:
                    if row != [] and len(row) > 3:
                        self.__rmBranch(row)
                os.remove(old_path)

            else:
                print('Error: Unable to get overview.txt')
        else:
            print("Overview file already up to date")

    def __createOverview(self, url): # Télécharge le fichier overview à partir de l'url.
        r = requests.get(url)
        if (r.status_code == 200):
            if not os.path.isdir('Results'):
                os.mkdir('Results')
            with open(self.filename, "w") as file:
                file.writelines(r.text)
                file.close()
        else:
            print('Error: Unable to get overview.txt')

    def __parseOverview(self, bar=False):
        """Crée l'arbre et les répertoires correspondants à partir de 
        l'attribut filename."""

        print("Parsing overview file")
        if bar:
            with open(self.filename, "r") as f:
                size = sum(1 for _ in f)
        with open(self.filename) as csv_file:
            reader = csv.reader(csv_file, delimiter='\t')
            next(reader, None)
            if bar:
                i = 0
                for row in reader:
                    if row != [] and len(row) > 3:
                        i += 1
                        self.__addBranch(row)
            else:
                for row in reader:
                    if row != [] and len(row) > 3:
                        self.__addBranch(row)
            csv_file.close()

    def __addBranch(self, row): # row : ["Kingdom", "Group", "Subgroup", "Organism"]
        current_node = self
        path = 'Results'
        for i in [1,2,3,0]:
            name = row[i]
            path = path + '/' + format_filename(name)
            path = format_path(path)
            if not os.path.isdir(path):
                os.mkdir(path)
            current_node = current_node.addChild(Node(name, path))

    def __rmBranch(self, row): # Supprime une feuille de l'arbre à partir d'une ligne.
        print("Removing branch: " + str(row[:4]))
        current_node = self
        for i in [1,2,3,0]:
            nb_children = len(current_node.children)
            k = 0
            found = False
            while k < nb_children and not found:
                if current_node.children[k].data == row[i]:
                    if i == 0:
                        print("delete...")
                        child = current_node.children[k]
                        shutil.rmtree(child.path, ignore_errors=True)
                        del current_node.children[k]
                    else:
                        current_node = current_node.children[k]
                    found = True
                k += 1
        
    def print_tree(self, node, sep = ''):
        """Affiche l'arbre."""
        if(node == None):
            print("No tree to print")
            return
        print(sep + node.data)
        for child in node.children:
            self.print_tree(child, sep + '  ')
            
    def printProgressBar (
            self, iteration, total, decimals = 1, length = 50, fill = '█'):
        """Affiche une barre de progression, en fonction de l'état actuel
        et du nombre total d'itérations."""

        printEnd = "\r" if iteration + 1 < total else "\n"
        percent = ("{0:." + str(decimals) + "f}").format(
            100 * (iteration / float(total)))
        filledLength = int(length * iteration // total) + 1
        bar = fill * filledLength + ' ' * (length - filledLength)
        print(f'\r |{bar}| {percent}% ', end = printEnd)
