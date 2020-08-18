

import os
import shutil

def test_multigeneblast(command, remove=True):
    os.system(command)

    expected_files = ['clusterblast_output.mgb', 'run.log', 'jquery-1.4.2.min.js', 'jquery.svg.js', 'jquery.svgdom.js', 'style.css']
    my_path = os.path.abspath(os.path.dirname(__file__))
    db_folder = os.path.join(my_path, "test_run")

    #check if all files are present
    for path, subdirs, files in os.walk(db_folder):
        #ignore checking the image directory
        if os.path.basename(path) == "images":
            continue
        for file in files:
            if os.path.basename(path) == "visual":
                assert file in expected_files or (file.startswith("displaypage") and file.endswith(".xhtml"))
            else:
                assert file in expected_files

    #remove the folder
    if remove:
        shutil.rmtree(db_folder)




if __name__ == "__main__":

    #command using gb query
    print("Starting run with genbank file:")
    command1 = "..\\multigeneblast.py -in test_querys\\V_maris_query.gb -from 0 -to 59320 -out tests\\test_run_gbk -db test_data_base\\v_maris.pal -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 -sw 0.5 -m y -op 3"
    test_multigeneblast(command1, remove=False)
    print("Finished running gb query test...")
    print()

    #command using embl query
    print("Starting run with embl file:")
    command2 = "..\\multigeneblast.py -in test_querys\\V_maris_query.embl -from 0 -to 59320 -out tests\\test_run_embl -db test_data_base\\v_maris.pal -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 -sw 0.5 -m y -op 3"
    test_multigeneblast(command2, remove=False)
    print("Finished running embl query test...")
    print()

    #command using architecture search fasta
    print("Starting run with architecture mode(fasta file):")
    command3 = "..\\multigeneblast.py -in test_querys\\V_maris_query_partial.fasta -out tests\\test_run_fasta -db test_data_base\\v_maris.pal -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 -sw 0.5 -m y -op 3"
    test_multigeneblast(command3, remove=False)
    print("Finished running architecture mode test...")
    print()