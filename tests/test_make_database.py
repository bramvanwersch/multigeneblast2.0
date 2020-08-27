#Basic test for the databse creation. Passing these tests means that all files
#are correctly created. This does not guarantee correct creation

import os
import shutil
from constants import PROT_DATABASE_EXTENSIONS, APPDATA, NUC_DATABASE_EXTENSIONS

def test_protein_database(command):
    assert os.system(command) == 0

    #make sure that even when chaning directories the directory can be found
    my_path = os.path.abspath(os.path.dirname(__file__))
    db_folder = os.path.join(my_path, "test_data_base")

    #check if all files are present
    db_folder_list = os.listdir(db_folder)
    expected_files = ["test" + ext for ext in PROT_DATABASE_EXTENSIONS]
    for file in expected_files:
        assert file in db_folder_list

    #make sure to remove the __database_file_label
    shutil.rmtree(db_folder)

def test_nucleotide_database(command):
    assert os.system(command) == 0

    #make sure that even when chaning directories the directory can be found
    my_path = os.path.abspath(os.path.dirname(__file__))
    db_folder = os.path.join(my_path, "test_data_base")

    #check if all files are present
    db_folder_list = os.listdir(db_folder)
    expected_files = ["test" + ext for ext in NUC_DATABASE_EXTENSIONS]
    for file in expected_files:
        assert file in db_folder_list

    #make sure to remove the __database_file_label
    shutil.rmtree(db_folder)


if __name__ == "__main__":
    #basic command using gbk
    print("Starting Makedatabase with gbk file:")
    command1 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.gbk -inf all"
    test_protein_database(command1)
    print("Finished running gbk test...")
    print()

    #command for embl
    print("Starting Makedatabase with embl file:")
    command2 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.embl -inf all"
    test_protein_database(command2)
    print("Finished running embl test...")
    print()

    #command for combined
    print("Starting Makedatabase with gbk and embl file file:")
    command3 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\GCA_000204155.1_ASM20415v1_genomic.gbk tests\\GCA_000204155.1_ASM20415v1_genomic.embl -inf all"
    test_protein_database(command3)
    print("Finished running combined test...")
    print()

    #command for testing wgs master record
    print("Starting Makedatabase with WGS record file:")
    command4 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\WGS_master_record_test.gb -inf all"
    test_protein_database(command4)
    print("Finished running WGS test...")
    print()

    #command for testing rare cases (missing sequences incomplete entries etc.)
    print("Starting Makedatabase with rare cases file:")
    command5 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\rare_case_test.embl -inf all"
    test_protein_database(command5)
    print("Finished running rare case test...")
    print()

    #command for testign nucleotide database creation
    print("Starting Makedatabase with raw nucleotide data file:")
    command6 = "..\\make_database.py -o tests\\test_data_base -n test -i tests\\rare_case_test.embl tests\\nucleotide_db_test.fasta tests\\GCA_000204155.1_ASM20415v1_genomic.gbk -inf all -t nucl"
    test_nucleotide_database(command6)
    print("Finished running nucleotide database test...")
    print()

    print("6/6 tests succesfully finished.")