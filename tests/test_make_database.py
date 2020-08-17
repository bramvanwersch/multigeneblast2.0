#Basic test for the databse creation. Passing these tests means that all files
#are correctly created. This does not guarantee correct creation

import os
import shutil
from constants import DATABASE_EXTENSIONS

def create_database(command):
    os.system(command)

    #make sure that even when chaning directories the directory can be found
    my_path = os.path.abspath(os.path.dirname(__file__))
    db_folder = os.path.join(my_path, "test_data_base")

    #check if all files are present
    db_folder_list = os.listdir(db_folder)
    expected_files = ["test" + ext for ext in DATABASE_EXTENSIONS]
    for file in expected_files:
        assert file in db_folder_list

    #make sure to remove the database
    shutil.rmtree(db_folder)


if __name__ == "__main__":
    #basic command using gbk
    command1 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.gbk"
    create_database(command1)
    print("Finished running gbk test...")
    print()

    #command for embl
    command2 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.embl"
    create_database(command2)
    print("Finished running embl test...")
    print()

    #command for combined
    command3 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.gbk tests\\GCA_000204155.1_ASM20415v1_genomic.embl"
    create_database(command3)
    print("Finished running combined test...")
    print()
    print("3/3 tests succesfully finished.")