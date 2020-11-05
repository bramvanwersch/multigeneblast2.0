#!/usr/bin/env python3

"""
Constants that are used by various functions of multigeneblast.

Original creator: Marnix Medena
Recent contributor: Bram van Wersch

Copyright (c) 2012 Marnix H. Medema
License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""


# imports
import os
import sys
import platform

# shared constants between files
HITS_PER_PAGE = 50
SCREENWIDTH = 1024
FASTA_EXTENSIONS = (".fasta", ".fas", ".fa", ".fna", ".faa")
EMBL_EXTENSIONS = (".embl", ".emb")
GENBANK_EXTENSIONS = (".gbk", ".gb", ".genbank", ".gbff")
SVG_CORE_EXTENSION = 5000
PROT_DATABASE_EXTENSIONS = ["_database_index.pickle", "_contigs.tar.gz", ".dmnd"]
NUC_DATABASE_EXTENSIONS = ["_database_index.pickle", "_contigs.tar.gz", ".nsq", ".nin", ".nhr"]

CHUNK = 256 * 1024

# smoothing paramater that allows slight overlap between genes when filtering redundant blast hits
ALLOWED_OVERLAP = 500

# name of the folder containing the final results
OUT_FOLDER_NAME = "multigeneblast_result"

# last part of the name of the text output
TEXT_OUT_NAME = "_cluster_text.mgb"


# path constants
def get_mgb_path():
    """
    Get the path to the file where mgb is installed. If multigeneblast is not
    added to the path environmental variable then calling this method outside of
    the folder with mgb in it will result in an exception

    :return: a string path or exception
    """
    # Find path to mgb files if run from another directory
    pathfolders = os.getenv('PATH').split(os.pathsep)
    pathfolders.reverse()
    pathfolders.append(os.getcwd())
    pathfolders.reverse()
    mgb_path = ""
    for folder in pathfolders:
        try:
            if "guilib.py" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and \
                    "multigeneblast_gui.py" in os.listdir(folder):
                mgb_path = folder
                break
        except Exception:
            pass
    try:
        if mgb_path == "" and os.sep in sys.argv[0] and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
            mgb_path = os.path.join(os.getcwd(), sys.argv[0].rpartition(os.sep)[0])
    except Exception:
        pass
    if mgb_path == "" or not os.path.exists(mgb_path):
        raise Exception("Error: Please add the MultiGeneBlast installation directory to"
                        " your $PATH environment variable before running the executable from another folder.")
    return mgb_path


def get_appdata_path():
    """
    Find a path to the appdata of the operating system

    :return: a string path or exception
    """
    # Find path to Application Data
    if platform.system() == "Windows":
        # roaming appdata folder of windows
        appdata = os.getenv("APPDATA")
    elif platform.system() == "Darwin":
        appdata = os.path.expanduser("~") + "/Library/Application Support"
    elif platform.system() == "Linux":
        try:
            if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
                appdata = get_mgb_path() + os.sep + "multigeneblast_data"
            else:
                os.mkdir(get_mgb_path() + os.sep + "multigeneblast_data")
                appdata = get_mgb_path() + os.sep + "multigeneblast_data"
        except Exception:
            try:
                if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
                    appdata = os.getcwd() + os.sep + "multigeneblast_data"
                else:
                    os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
                    appdata = os.environ['HOME'] + os.sep + "multigeneblast_data"
            except Exception:
                raise Exception("No permission to write to installation folder. "
                                "Please change user or save somewhere else.")
    else:
        raise Exception("MultiGeneBlast does not support {}".format(platform.system()))
    if platform.system() == "Darwin" or platform.system() == "Windows":
        try:
            os.mkdir(appdata + os.sep + 'MultiGeneBlast')
            appdata = appdata + os.sep + 'MultiGeneBlast'
        except Exception:
            if os.path.exists(appdata + os.sep + 'MultiGeneBlast'):
                appdata = appdata + os.sep + 'MultiGeneBlast'
    return appdata


def get_temp_data():
    """
    Return the path to a mgb_temp folder in the temp folder of the operating
    system

    :return: a file path to that folder.

    It is not guaranteed the folder exists. It is common practice to remove the
    folder to clean up files.
    """
    # Find path to temporary files
    if platform.system() == "Windows":
        temp = os.getenv('TEMP')
    elif platform.system() == "Darwin":
        temp = os.getenv('TMPDIR')
    else:
        try:
            os.mkdir(os.getenv('HOME') + os.sep + ".mgbtemp")
            temp = os.getenv('HOME') + os.sep + ".mgbtemp"
        except Exception:
            temp = APPDATA
    # return a unique folder name to make sure that no unwanted files are deleted when cleaning the folder
    return temp + os.sep + "mgb_temp"


def get_exec_folder():
    """
    Locate the exec folder for the current operating system

    :return: string name that is an absolute path to the folder with all
    the external program
    """
    base = get_mgb_path()
    if platform.system() == "Windows":
        return "{}{}exec{}windows_exec".format(base, os.sep, os.sep)
    if platform.system() == "Darwin":
        return "{}{}exec{}mac_exec".format(base, os.sep, os.sep)
    if platform.system() == "Linux":
        return "{}{}exec{}linux_exec".format(base, os.sep, os.sep)


def get_muscle_prog_name():
    """
    Get the correct name of the muscle program based on the operating system
    and the bitmode python is running in.

    :return: a string that is the correct name of the muscle instalation
    """
    if platform.system() == "Windows":
        return "muscle3.8.31_i86win32.exe"
    bit64_mode = sys.maxsize > 2**32
    if bit64_mode:
        if platform.system() == "Darwin":
            return "muscle3.8.31_i86darwin64"
        else:
            return "muscle3.8.31_i86linux64"
    else:
        if platform.system() == "Darwin":
            return "muscle3.8.31_i86darwin32"
        else:
            return "muscle3.8.31_i86linux32"


# path constants
CURRENTDIR = os.getcwd()
APPDATA = get_appdata_path()
TEMP = get_temp_data()
MUSCLE_PROG_NAME = get_muscle_prog_name()
