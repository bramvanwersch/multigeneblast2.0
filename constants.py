
#imports
import os, sys


#shared constants between files
HITS_PER_PAGE = 50
FASTA_EXTENSIONS = (".fasta",".fas",".fa",".fna")
EMBL_EXTENSIONS = (".embl",".emb")
GENBANK_EXTENSIONS = (".gbk",".gb",".genbank")
SVG_CORE_EXTENSION = 5000
DATABASE_EXTENSIONS = ["_database_index.pickle","_contigs.tar.gz",".psq",".psi",".psd",".pin",".phr"]

#path constants
def get_mgb_path():
    #Find path to mgb files if run from another directory
    pathfolders = os.environ['PATH'].split(os.pathsep)
    pathfolders.reverse()
    pathfolders.append(os.getcwd())
    pathfolders.reverse()
    mgb_path = ""
    for folder in pathfolders:
        try:
            if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
                mgb_path = folder
                break
        except:
            pass
    try:
        if  mgb_path == "" and os.sep in sys.argv[0] and "read_input_gui.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
            mgb_path = sys.argv[0].rpartition(os.sep)[0]
            os.chdir(mgb_path)
    except:
        pass
    if mgb_path == "":
        raise Exception("Error: Please add the MultiGeneBlast installation directory to"
                        " your $PATH environment variable before running the executable from another folder.")
    return mgb_path


def get_appdata_path():
    #Find path to Application Data
    if sys.platform == ('win32'):
        appdata = os.environ['ALLUSERSPROFILE'] + os.sep + 'Application Data'
    elif sys.platform == ('darwin'):
        appdata = os.path.expanduser("~") + "/Library/Application Support"
    else:
        try:
            if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
                appdata = os.getcwd() + os.sep + "multigeneblast_data"
            else:
                os.mkdir(os.getcwd() + os.sep + "multigeneblast_data")
                appdata = os.getcwd() + os.sep + "multigeneblast_data"
        except:
            try:
                if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
                    appdata = os.getcwd() + os.sep + "multigeneblast_data"
                else:
                    os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
                    appdata = os.environ['HOME'] + os.sep + "multigeneblast_data"
            except:
                raise Exception("No permission to write to installation folder. Please change user or save somewhere else.")

    if sys.platform == ('darwin') or sys.platform == ('win32'):
        try:
            os.mkdir(appdata + os.sep + 'MultiGeneBlast')
            appdata = appdata + os.sep + 'MultiGeneBlast'
        except:
            if os.path.exists(appdata + os.sep + 'MultiGeneBlast'):
                appdata = appdata + os.sep + 'MultiGeneBlast'
    return appdata

def get_temp_data():
    #Find path to temporary files
    if sys.platform == ('win32'):
        temp = os.environ['TEMP']
    elif sys.platform == ('darwin'):
        temp = os.environ['TMPDIR']
    else:
        try:
            os.mkdir(os.environ['HOME'] + os.sep + ".mgbtemp")
            temp = os.environ['HOME'] + os.sep + ".mgbtemp"
        except:
            temp = APPDATA
    return temp

#path constants
MGBPATH = get_mgb_path()
CURRENTDIR = os.getcwd()
APPDATA = get_appdata_path()
TEMP = get_temp_data()