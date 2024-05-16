import os

def is_windows():
    return os.name == "nt"

def format_filename(filename):
    if is_windows():
        for symb in "<>:\\\"/|?*":
            filename = filename.replace(symb, '-')
    else:
        filename = filename.replace('/', '-')
    return filename

def format_path(path):
    if is_windows():
        return path.replace('/', '\\')
    return path