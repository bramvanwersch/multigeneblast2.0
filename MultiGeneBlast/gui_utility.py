#!/usr/bin/env python3

"""
Utility functions used by several classes from the gui.

Original creator: Marnix Medena
Recent contributor: Bram van Wersch

Copyright (c) 2012 Marnix H. Medema
License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""


import subprocess
import os
import urllib.parse
import webbrowser
from tkinter import StringVar, Entry, IntVar, Checkbutton, Spinbox, Toplevel, TOP, Scrollbar, RIGHT, Y, Listbox, \
    EXTENDED, LEFT, END, BOTTOM, Button, Text, DISABLED, NORMAL, INSERT
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showerror
from tkinter.ttk import Frame, Scale, Label


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
            outbox.text_insert("\nMultigeneblast exited because an input variable was not complete. Make sure to fix "
                               "it and run Multigeneblast again OR if you think this is not a problem with the input "
                               "please click the button below to send an error report.")
            expected = True
        outbox.change_errormessage(error)
        outbox.add_ok_button()
        outbox.add_report_button()
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
    # if not files are present test if you can write in the folder easier to ask forgiveniss then permission
    try:
        with open(selected + os.sep + "test.txt", "w") as f:
            f.write("test")
        os.remove(selected + os.sep + "test.txt")
    except PermissionError:
        showerror("Error", "No permission to write to this folder. Please choose a directory in which you have"
                           " writing permissions.")
        select_out_directory()
    return selected


def about():
    """
    Open a webpage that displayes information about how to use MultiGeneBlast
    """
    webbrowser.open('http://multigeneblast.sourceforge.net/')


class ScaleBar(Frame):
    """
    Custom widget that is a Scale and Entry widget combined
    """
    def __init__(self, parent, minimum, maximum, default="0", scale_command=None, input_type="int"):
        """
        :param parent: The parent widget holding this one
        :param minimum: The minimum value of the scale
        :param maximum: The maximum value of the scale
        :param default: the startvalue of the scale
        :param scale_command: Optional function object that can be calles when the scale is issued a command
        :param input_type: the input type of the widget
        """
        Frame.__init__(self, parent)

        self.parent = parent
        self.minimum = minimum
        self.maximum = maximum
        self.default = default
        self.command = scale_command
        self.input_type = input_type

        # place holder value to signify that this attribute is expected to exist
        self.scale = None
        self.entry = None

        self.var = StringVar()
        self.lastvar = StringVar()

        self.init_widgets()

    def init_widgets(self):
        """
        Innitialize widgets
        """
        scale = Scale(self, from_=self.minimum, to=self.maximum, command=self.on_scale, length=200)
        scale.grid(row=0, column=1)
        self.scale = scale
        self.var.set(self.default)
        self.lastvar.set(self.var.get())
        scale.set(self.default)
        self.entry = Entry(self, textvariable=self.var, width=8)
        self.entry.bind("<FocusOut>", self.on_validate)
        self.entry.bind("<Return>", self.on_validate)
        self.entry.grid(row=0, column=0)

    def set_scale(self, start, end, value):
        """
        Set the values of the scale widget

        :param start: start value
        :param end: end value
        :param value: current value
        """
        self.scale.configure(from_=start)
        self.scale.configure(to=end)
        self.var.set(str(value))
        self.scale.set(value)

    def on_validate(self):
        """
        Called when a user finishes enterign something in the entry widget
        """
        if "-" in self.var.get():
            self.var.set(self.minimum)
        if not str(self.var.get()).isdigit():
            self.var.set(self.lastvar.get())
        if int(self.var.get()) <= int(self.minimum):
            self.var.set(self.minimum)
        if int(self.maximum) - int(self.var.get()) < 0:
            self.var.set(self.maximum)
        self.on_scale(self.var.get())
        self.scale.set(self.var.get())
        self.lastvar.set(str(self.var.get()))

    def on_scale(self, val):
        """
        Called when the entry widget requires the scale bar to change
        """
        if self.input_type == "int":
            v = int(float(val))
        else:
            v = float("".join(str(val).partition(".")[0:2]) + str(val).partition(".")[2][:2])
        self.var.set(str(v))
        if self.command is not None:
            self.command()

    def setval(self, val):
        """
        Set the value of the scale bar

        :param val: The value
        """
        self.var.set(val)
        self.scale.set(val)

    def getval(self):
        """
        Get the value of the scale bar

        :return: the variable of the scale bar
        """
        return self.var.get()


class CheckBox(Frame):
    """
    Custom implementation of a checkbox
    """
    def __init__(self, parent, description):
        """
        :param parent: The parent widget holding this one
        :param description: an optional text in front of the button
        """
        Frame.__init__(self, parent)
        self.description = description
        self.var = IntVar()
        self.init_widgets()
        self.var.set(0)

    def init_widgets(self):
        """
        innitialize the widgets
        """
        cb = Checkbutton(self, text=self.description, variable=self.var)
        cb.select()
        cb.grid(row=0, column=0)

    def getval(self):
        """
        Get the value of of the button 0 or 1
        :return: an integer
        """
        return self.var.get()


class CustomSpinBox(Frame):
    """
    Custom SpinBox widget
    """
    def __init__(self, parent, minimum, maximum, incrval, default=0):
        """
        :param parent: The parent widget holding this one
        :param minimum: The minimum value of the spinbox
        :param maximum: The maximum value of the spinbox
        :param incrval: The value with which the spinbox should increment
        :param default: The default value for the spinbox
        """
        Frame.__init__(self, parent)

        self.parent = parent
        self.minimum = minimum
        self.maximum = maximum
        self.incrval = incrval
        self.default = default
        self.var = IntVar()
        self.init_widgets()

    def init_widgets(self):
        """
        innitialize widgets
        """
        self.var.set(self.default)

        cb = Spinbox(self, from_=self.minimum, to=self.maximum, increment=self.incrval, textvariable=self.var)
        cb.grid(row=0, column=0)

    def getval(self):
        """
        Get the value of the spinbox

        :return: an Integer
        """
        return self.var.get()


class ListBoxChoice(Toplevel):
    """
    Custom ListBox
    """
    def __init__(self, master, title=None, message=None, list_=None):
        """
        :param master: the parent widget
        :param title: an optional title
        :param message: an optional message
        :param list_: a list of entries to start with
        """
        super().__init__(master)
        self.value = None
        if list_ is None:
            self.list_ = []
        else:
            self.list_ = list_
        self.list_box = None
        self.init_widgets(title, message)

    def init_widgets(self, title, message):
        """
        innitialize widgets

        :param title: optinal title
        :param message: optional message
        """
        self.transient()
        self.grab_set()

        self.bind("<Return>", self._choose)
        self.bind("<Escape>", self._cancel)

        if title:
            self.title(title)

        if message:
            Label(self, text=message).pack(padx=5, pady=5)

        list_frame = Frame(self)
        list_frame.pack(side=TOP, padx=5, pady=5)

        scroll_bar = Scrollbar(list_frame)
        scroll_bar.pack(side=RIGHT, fill=Y)
        self.list_box = Listbox(list_frame, selectmode=EXTENDED)
        self.list_box.pack(side=LEFT, fill=Y)
        scroll_bar.config(command=self.list_box.yview)
        self.list_box.config(yscrollcommand=scroll_bar.set)
        self.list_.sort()
        for item in self.list_:
            self.list_box.insert(END, item)

        button_frame = Frame(self)
        button_frame.pack(side=BOTTOM)

        choose_button = Button(button_frame, text="Choose", command=self._choose)
        choose_button.pack()

        cancel_button = Button(button_frame, text="Cancel", command=self._cancel)
        cancel_button.pack(side=RIGHT)

    def _choose(self):
        """
        Save all selected entries in a the value property
        """
        try:
            if len(self.list_box.curselection()) == 1:
                selected = self.list_box.curselection()[0]
                self.value = self.list_[int(selected)]
            else:
                selected = self.list_box.curselection()
                self.value = ";".join([self.list_[int(idx)] for idx in selected])
        except IndexError:
            self.value = None
        self.destroy()

    def _cancel(self):
        """
        Destroy this TopLevel
        """
        self.destroy()

    def return_value(self):
        """
        Value returned before the widget is destroyed

        :return: string of the self.value property
        """
        self.master.wait_window(self)
        return self.value


class MessageBox(Toplevel):
    """
    Customm Text widget for displaying
    """
    def __init__(self, master, title=None, list_=None):
        """
        :param master: the parent widget
        :param title: an optional title
        :param list_: a list of entries to start with
        """
        super().__init__(master)
        self.master = master
        self.value = None
        self.errormessage = ""
        if list_ is None:
            self.list_ = []
        else:
            self.list_ = list_

        self.messageFrame = Frame(self)
        self.buttonFrame = Frame(self)
        self.init_widgets(title)

    def init_widgets(self, title):
        """
        Innitialize the widgets

        :param title: optional title
        """
        self.transient()
        self.grab_set()

        self.bind("<Escape>", self._cancel)
        self.protocol('WM_DELETE_WINDOW', self._cancel)

        if title:
            self.title(title)

        # define a new frame and put a text area in it
        self.messageFrame.pack(side=TOP)
        self.messageFrame.text = Text(self.messageFrame, height=40, width=100, background='white', state=DISABLED)
        self.messageFrame.text.tag_config('Error', foreground="red")
        self.messageFrame.text.tag_config('Warning', foreground="blue")

        # put a scroll bar in the frame
        self.messageFrame.scroll = Scrollbar(self.messageFrame)
        self.messageFrame.text.configure(yscrollcommand=self.messageFrame.scroll.set)

        # pack everything
        self.messageFrame.text.pack(side=LEFT)
        self.messageFrame.scroll.pack(side=RIGHT, fill=Y)

    def text_insert(self, text, tag=None):
        """
        Insert a message into the MessageBox. Allow a potential tag to color certain messages.

        :param text: a string to add as message
        :param tag: two options avaialble 'Error' and 'Warning', these tags color the text displayed
        """
        if isinstance(text, bytes):
            text = text.decode("utf8")
        self.messageFrame.text.see(END)
        self.messageFrame.text.config(state=NORMAL)
        if tag is not None:
            self.messageFrame.text.insert(INSERT, text, tag)
        elif "WARNING" in text:
            self.messageFrame.text.insert(INSERT, text, "Warning")
        elif "ERROR" in text:
            self.messageFrame.text.insert(INSERT, text, "Error")
        else:
            self.messageFrame.text.insert(INSERT, text)
        self.messageFrame.text.config(state=DISABLED)

    def add_ok_button(self):
        """
        Add an Ok button to close the window
        """

        self.buttonFrame.pack(side=BOTTOM)

        ok_button = Button(self.buttonFrame, text="OK", command=self._cancel, width=10)
        ok_button.pack()

    def add_report_button(self):
        """
        Add a button to allow the user to report a bug by providing
        a button that redirects them to an email
        """
        self.buttonFrame = Frame(self)
        self.buttonFrame.pack(side=BOTTOM)

        error_button = Button(self.buttonFrame, text="Send Error Report", command=self.send_report, width=20)
        error_button.pack()

    def change_errormessage(self, errormessage):
        """
        The message send with the error report

        :param errormessage: a string with the error
        :return:
        """
        self.errormessage = errormessage

    def send_report(self):
        """
        Function triggered when the send error report button is clicked
        """
        # TODO make this more proper by sending automatically and add a log file. Logging module offers options to
        #  automatically send on error
        try:
            webbrowser.open("mailto:multigeneblast@gmail.com?SUBJECT=ErrorReport&BODY=" +
                            urllib.parse.quote(self.errormessage.encode("utf8")))
        except Exception:
            webbrowser.open("sourceforge.net/tracker/?func=add&group_id=565495&atid=2293721")
        else:
            pass

    def _cancel(self, event=None):
        """
        Destroy this TopLevel

        :param event: optional argument
        """
        self.destroy()
