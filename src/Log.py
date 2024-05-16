import os
from os_compatibility import format_path
from pathlib import Path

class Log():
    def __init__(self, fname):
        self.create(fname)
        self.fname = format_path(fname)
    
    def add(self):
        with open(self.fname, "+a") as file:
            file.write("0\n")
    
    def create(self, fname):
        reps = fname.split('/')[:-1]
        for rep in reps:
            try:
                os.mkdir(rep)
            except FileExistsError:
                pass
        Path(fname).touch()

    def delete(self):
        os.remove(self.fname)
    
    def gen_list(self):
        with open(self.fname) as file:
            return file.readlines()