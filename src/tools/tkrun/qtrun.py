#! /usr/bin/env python
#
import sys
import os
from PyQt5 import QtCore, QtWidgets
import argparse
import re
import subprocess

_version = "4-oct-2023"
_debug = False

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, parameters, input_file, filetype):
        super(MainWindow, self).__init__()

        self.groups = parameters
        self.input_file = input_file
        self.input_file_type = filetype
        self.radio_groups = []
        self.slider_multiplier = []
        self.sliders = []

        self.initUI()
    
    def initUI(self):
        self.pagelayout = QtWidgets.QVBoxLayout()       #page layout

        # run, save, load, quit, help buttons -> located in a toolbar
        toolbar = self.addToolBar("ToolBar")

        run_action = QtWidgets.QAction('Run', self)
        dry_action = QtWidgets.QAction('Dryrun', self)
        save_action = QtWidgets.QAction('Save', self)
        load_action = QtWidgets.QAction('Load', self)
        quit_action = QtWidgets.QAction('Quit', self)
        help_action = QtWidgets.QAction('Help', self)

        toolbar.addAction(run_action)
        toolbar.addSeparator()
        toolbar.addAction(dry_action)
        toolbar.addSeparator()
        toolbar.addAction(save_action)
        toolbar.addSeparator()
        toolbar.addAction(load_action)
        toolbar.addSeparator()
        toolbar.addAction(quit_action)
        toolbar.addSeparator()
        toolbar.addAction(help_action)

        run_action.triggered.connect(self.run)
        dry_action.triggered.connect(self.dry)
        save_action.triggered.connect(self.save)
        load_action.triggered.connect(self.load)
        quit_action.triggered.connect(self.quit)
        help_action.triggered.connect(self.help)
        
        # set the main page layout
        widget = QtWidgets.QWidget()
        widget.setLayout(self.pagelayout) 

        # add scrollbar
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(widget)
        self.setCentralWidget(scroll)

        self.createWidgetsFromGroups()
    
    # runs the file. takes the input file and runs it, piping the set paraments into the file as args
    # it assumes the script can self-execute with the appropriate #! directive on the first line
    def run(self):
        print("running: " + self.input_file)
        contents = self.gather_data()
        param = self.input_file
        for line in contents:
            for key, value in line.items():
                param += f" {key}={value}"
        subprocess.run(param.split())

    def dry(self):
        cmd = self.input_file
        contents = self.gather_data()
        param = self.input_file
        for line in contents:
            for key, value in line.items():
                cmd += f" {key}={value}"
        print(cmd)
    
    # saves the options into a separate file named inputfilename.key. This will be saved in the format
    # key=value and formatted according to the input file type.
    def save(self):
        contents = self.gather_data()
        default_file_path = self.input_file + ".key"
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", default_file_path, "All Files (*)")
        if file_path:
            with open(file_path, "w") as file:
                file.write("# written by qtrun.py - do not edit\n")
                for line in contents:
                    for key, value in line.items():
                        file.write(f"{key}={value}")
                    file.write("\n")
            print("saved to " + file_path)
    
    # loads the saved options from a previous session into the current gui. Will default to the file
    # named inputfilename.key
    def load(self):
        options = QtWidgets.QFileDialog.Options()
        load_file = self.input_file +".key" if os.path.exists(self.input_file +".key") else ""
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a File", load_file, "All Files (*)", options=options)
        
        #parse the loaded file
        if file:
            #turn the parameters in the file into a dictionary
            default_values = {}
            with open(file, "r") as f:
                for line in f:
                    if line[0] == '#':
                        continue
                    if self.input_file_type == 'csh':
                        line = re.sub("set","",line,count=1)
                    label, value = line.strip().split("=")
                    if value.startswith('"') and value.endswith('"') or value.startswith("'") and value.endswith("'"):
                        value = value[1:-1]
                    value = value.split(",") if "," in value else [value]
                    default_values[label] = value
            
            #go through the elements in the widget and alter it to the values specified
            for elements in range(self.pagelayout.count()):
                element = self.pagelayout.itemAt(elements).layout()
                if element is not None:
                    for widget_index in range(element.count()):
                        widget = element.itemAt(widget_index).widget()
                        
                        if isinstance(widget, QtWidgets.QLineEdit):
                            try:
                                widget.setText(''.join(default_values[widget.objectName()]))
                            except:
                                print("Skipping resetting new entry",widget.objectName())                                

                        elif isinstance(widget, QtWidgets.QRadioButton) or isinstance(widget, QtWidgets.QCheckBox):
                            try:
                                if widget.text() in default_values[widget.objectName()]:
                                    widget.setChecked(True)
                            except:
                                print("Skipping resetting new button",widget.objectName())

                        elif isinstance(widget, QtWidgets.QSlider):
                            try:
                                multiplier = self.slider_multiplier.pop(0)
                                widget.setValue(int(float(''.join(default_values[widget.objectName()]))*multiplier))
                                self.slider_multiplier.append(multiplier)
                            except:
                                print("Skipping resetting new slider",widget.objectName())

    def quit(self):
        self.close()

    def help(self):
        print("qtrun version %s" % _version)
        print("")
        print("Run:     run the script withe the selected `key=val` arguments")
        print("Dryrun:  show the command that would be run")
        print("Save:    save the key=val settings in a key file for later retrieval")
        print("Load:    load a keyfile previously saved")
        print("Quit:    quit the program")
        print("Help:    this help")        
        print("")
        print("Hover with the mouse over the keyword names to get some help on the keywords")
        print("The current script is: ",args.input_file[0])
        # print("The current keyfile is: ",args.input_file[1])

    # function to create each widget from the input file
    def createWidgetsFromGroups(self):
        for group in self.groups:
            group_type, group_name, options, default_option, help = group

            if group_type == "RADIO":
                new_group = QtWidgets.QButtonGroup()
                self.radio_groups.append(new_group)
                group_layout = QtWidgets.QHBoxLayout()
                label = QtWidgets.QLabel(group_name+":")
                label.setToolTip(help)
                group_layout.addWidget(label)
                for option in options:
                    option = option.strip()
                    radio_button = QtWidgets.QRadioButton(option)
                    radio_button.setObjectName(group_name)
                    new_group.addButton(radio_button)
                    group_layout.addWidget(radio_button)
                    if option in default_option:
                        radio_button.setChecked(True)
                self.pagelayout.addLayout(group_layout)

            elif group_type == "IFILE" or group_type == "OFILE" or group_type == "IDIR" or group_type == "ODIR":
                group_layout = QtWidgets.QHBoxLayout()
                label = QtWidgets.QLabel(group_name+":")
                label.setToolTip(help)
                group_layout.addWidget(label)
                btn = QtWidgets.QPushButton(self)
                btn.setText("browse...")
            
                txt = QtWidgets.QLineEdit(self)
                txt.setText(default_option)
                txt.setObjectName(group_name)
                if group_type == "OFILE" or group_type == "IFILE":
                    self.ofile = group_name
                    btn.clicked.connect(lambda _, edit=txt: self.browse("FILE", edit))
                else:
                    btn.clicked.connect(lambda _, edit=txt: self.browse("DIR", edit))
                group_layout.addWidget(btn)
                group_layout.addWidget(txt)
                self.pagelayout.addLayout(group_layout)

            elif group_type == "CHECK":
                group_layout = QtWidgets.QHBoxLayout()
                label = QtWidgets.QLabel(group_name+":")
                label.setToolTip(help)
                group_layout.addWidget(label)
                for option in options:
                    option = option.strip()
                    checkbox = QtWidgets.QCheckBox(option, self)
                    checkbox.setObjectName(group_name)
                    group_layout.addWidget(checkbox)
                    if option in default_option:
                        checkbox.setChecked(True)
                self.pagelayout.addLayout(group_layout)

            elif group_type == "ENTRY":
                group_layout = QtWidgets.QHBoxLayout()
                label = QtWidgets.QLabel(group_name+":")
                label.setToolTip(help)
                group_layout.addWidget(label)
                txt = QtWidgets.QLineEdit(self)
                txt.setObjectName(group_name)
                txt.setText(default_option)
                group_layout.addWidget(txt)
                self.pagelayout.addLayout(group_layout)
            
            elif group_type == "SCALE":
                group_layout = QtWidgets.QHBoxLayout()
                label = QtWidgets.QLabel(group_name+":")
                label.setToolTip(help)
                group_layout.addWidget(label)
                options = ''.join(options)
                options = options.split(':')

                #creates a horizontal decimal slider
                decimals = len(str(options[2]).split('.')[1]) if '.' in str(options[2]) else 0
                multiplier = 10**decimals
                slider = QtWidgets.QSlider(self)
                slider.setOrientation(QtCore.Qt.Horizontal)
                slider.setSingleStep(int(float(options[2])*multiplier))
                slider.setPageStep(int(float(options[2])*multiplier))       #moves the slider when clicking or up/down
                slider.setRange(int(options[0])*multiplier, int(options[1])*multiplier)
                slider.setValue(int(float(default_option)*multiplier))

                slider_label = QtWidgets.QLabel(f"{slider.value()/multiplier}", self)

                slider.valueChanged.connect(lambda: self.update_label()) # updates to current value
                slider.setObjectName(group_name)
                self.slider_multiplier.append(multiplier)
                self.sliders.append((slider, slider_label, multiplier))
                group_layout.addWidget(slider_label)
                group_layout.addWidget(slider)
                self.pagelayout.addLayout(group_layout)

            # Create a visual separator (horizontal line)
            separator = QtWidgets.QFrame()
            separator.setFrameShape(QtWidgets.QFrame.HLine)
            separator.setFrameShadow(QtWidgets.QFrame.Sunken)
            separator.setFixedHeight(1)  # Set a fixed height for the separator
            self.pagelayout.addWidget(separator)

            if args.debug:
                print(f"{group_name} created as {group_type}")

    def update_label(self):
        for slider, label, multiplier in self.sliders:
            label.setText(f"{slider.value()/multiplier}")

    def browse(self, gtype, txt):
        options = QtWidgets.QFileDialog.Options()
        file = None
        dir = None
        if gtype == 'FILE':
            file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a File", "", "All Files (*)", options=options)
        if gtype == 'DIR':
            dir = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory")
        if file:
            txt.setText(file)
            print(file + " selected")
        if dir:
            txt.setText(dir)
            print(dir + " selected")
        
    def gather_data(self):
        layout_data = []

        for hbox_layout_index in range(self.pagelayout.count()):
            hbox_layout = self.pagelayout.itemAt(hbox_layout_index).layout()
            
            if hbox_layout is not None:
                defaults = {}
                key = None

                for widget_index in range(hbox_layout.count()):
                    widget = hbox_layout.itemAt(widget_index).widget()
                    
                    if widget_index == 0 and isinstance(widget, QtWidgets.QLabel):
                        key = widget.text().split(':')[0]
                        if self.input_file_type == "csh":
                            key = "set " + key
                        defaults[key] = []

                    elif isinstance(widget, QtWidgets.QLineEdit):
                        value = widget.text()
                        defaults[key].append(value)

                    elif isinstance(widget, QtWidgets.QRadioButton) or isinstance(widget, QtWidgets.QCheckBox):
                        if widget.isChecked():
                            value = widget.text()
                            defaults[key].append(value)

                    elif isinstance(widget, QtWidgets.QSlider):
                        multiplier = self.slider_multiplier.pop(0)
                        value = str(widget.value()/multiplier)
                        defaults[key].append(value)
                        self.slider_multiplier.append(multiplier)

                values = defaults[key]
                values = ','.join(values)
                
                if self.input_file_type == "python":
                    defaults[key] = "'" + values + "'"
                else:
                    defaults[key] = values
                layout_data.append(defaults)
        return layout_data

def parsefile(file):
    filetype = "sh"
    with open(file, "r") as f:
        lines = f.readlines()

    groups = []
    # group 1 = set or None
    # group 2 = name of widget
    # group 3 = default values, may have "" around
    # group 4 = # help or None
    # group 5 = widget type
    # group 6 (unused) = name=value if old format, otherwise None
    # group 7 = widget parameters or None
    pattern = '^\s*(set\s+)?([^#]+)\s*=([^#]+)(#.*)?#>\s+([^\s]+)(.*=[^\s]*)?(.+)?$'
    for line in lines:
        match = re.match(pattern, line)
        if match:
            if match.group(1):
                filetype = "csh"
            group_type = match.group(5).strip()
            group_name = match.group(2).strip()
            default_option = match.group(3).strip()
            #check for quotations
            if (default_option[0] == '"' and default_option[-1] == '"') or (default_option[0] == "'" and default_option[-1] == "'"):
                default_option = default_option[1:-1]
                filetype = "python"
            options = match.group(7).split(',') if match.group(7) else ""
            help = match.group(4).split('#')[1].strip() if match.group(4) else ""
            groups.append((group_type, group_name, options, default_option, help))
    
    return groups, filetype

if __name__ == '__main__':
    global args # global to access args outside
    parser = argparse.ArgumentParser(description="Dynamic GUI Builder - version %s" % _version,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("input_file", help="Script containing `key=val` parameters and #> GUI directives\n" +
                        "optionally add a keyfile [not implemented yet]",
                        default=None, nargs='*')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('-v', '--version', action='store_true', help='Show version')

    args = parser.parse_args()
    if args.version:
        print(_version)
        sys.exit(0)

    _debug = args.debug

    input_file = args.input_file[0]
    if len(args.input_file) > 1:
        print("multiple args not yet supported")
    if _debug:
        print(input_file)
    groups, filetype = parsefile(input_file)        

    app = QtWidgets.QApplication(sys.argv)
    
    w = MainWindow(groups, os.path.abspath(input_file), filetype)
    # w.inputFile = args.input_file
    w.adjustSize()  #adjust to fit elements accordingly

    #sets a minimum window size
    w.setMinimumWidth(500) 
    w.setMinimumHeight(200)
    w.show()

    try:
        if _debug:
            print('opening window')
        sys.exit(app.exec())
    except SystemExit:
        if _debug:
            print('closing window')
