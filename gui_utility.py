import subprocess
import os
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showerror

def run_extrenal_command(command, outbox, base_frame):
    """
    Run an external command on the command line using subprocess. Report progress in the outbox MessageBox

    :param command: a string that is a command to be deployed on the command line
    :param outbox: a MessageBox object for reporting progress
    :param base_frame: the frame where the MessageBox is on top of.
    :return: the exit code of the program (0 is succes) and a boolean telling if the exit of the program was expected
    """
    # create a place to push messages to
    base_frame.update()
    popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    lines_iterator = iter(popen.stdout.readline, b"")
    while popen.poll() is None:
        for line in lines_iterator:
            outbox.text_insert(line)
            base_frame.update()
    # check exit message
    output, error = popen.communicate()
    if error != b'':
        expected = False
        string_error = error.decode('utf-8')
        error_lines = string_error.split("\n")

        outbox.text_insert(string_error + "\n", "Error")
        # expected error. Notify user of it.
        if "raise MultiGeneBlastException" in error_lines[-1] or "raise MultiGeneBlastException" in error_lines[-2] \
                or "raise MultiGeneBlastException" in error_lines[-3]:
            outbox.text_insert("\nMultigeneblast exitied because an input variable was not complete. Make sure to fix it "
                "and run Multigeneblast again OR if you think this is not a problem with the input "
                "please click the button below to send an error report.")
            expected = True
        outbox.change_errormessage(error)
        outbox.add_ok_button()
        outbox.add_error_button()
        base_frame.update()
    else:
        expected = True
        outbox.add_ok_button()
        base_frame.update()
    return popen.returncode, expected

def select_out_directory():
    """
    Select the output directory.
    """
    selected = askdirectory(mustexist=False)
    if selected == "":
        return
    #if not files are present test if you can write in the folder
    #easier to ask forgiveniss then permission
    try:
        with open(selected + os.sep + "test.txt", "w") as f:
            f.write("test")
        os.remove(selected + os.sep + "test.txt")
    except PermissionError:
        showerror("Error", "No permission to write to this folder. "
            "Please choose a directory in which you have writing permissions.")
        select_out_directory()
    return selected
