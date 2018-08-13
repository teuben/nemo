#! /bin/env python
# Last Modified on 08-18-08 by Brian J. Prager
   
class Test:   
   class PyRun:                                                                     # Begin Class to build GUI
      """class to run functions under different environments"""
      def __init__(self, fun, infile, lengthflag, programflag):
   		self._fun = fun                                                            # Store the Function
        
   		if lengthflag == 1:
   			print 'Please Provide Function containing GUI information.'
   			self.Quit()
        
   		self._argcnt = fun.func_code.co_argcount                                   # Store Number of Arguments in Function passed to PyRun
   		self._radcount = 0                                                         # Variable to work with Radio Widgets
   		self._checkcount = 0                                                       # Variable to work with Check Widgets
   		self._statecount = 0                                                       # Variable to work with Check Widgets
   		self._place = 0                                                            # Variable to help with window geometry
   		self._HELPList = []                                                        # List to Store Help Info from Doc String
   		self.Ifieldentry =0                                                        # Variable to help with Infile Entries
   		self.Ifieldentries =0                                                      # Variable to help with Infiles Entries
   		self.Ofieldentry =0                                                        # Variable to help with Outfile Entries
   		self.DirFieldentry=0                                                       # Variable to help with Directory Entries
   		self.StringNumber =[]                                                      # A list to keep track of InFiles Entries in order to split string
   		self.Program = str(infile)                                                 # Parent Program
   		self.ProgramFlag = programflag                                             # Flag to run a function or entire program
   		self.Infile = infile
        
   		tmp = self.Program.split("'")
   		self.blah = tmp[1]
   		self.ProgramName = "%s.py" %tmp[1]
        
      def GetValues(self):
         """Function to Grab Entries in the GUI"""
         runlist  = []                                                              # Make List to Store Dynamic Variables
         DynPassList = []                                                           # Create Dynamic List to Store Values to be passed
         self.RadioTypeList = []
         PassList = []                                                              # Static List of Values to be passed
         radcount = 0                                                               # Iterative Value to find Radio Values
         checkcount = 0                                                             # Iterative Value to find Check Values
         statecount = 0                                                             # Iterative Value of States that Store Check Values
         Passcount = 0                                                              # Store Length of Variable List to be Passed
        
         for a in range(self._argcnt):                                              # Loop over number of Arguments needed to form Dynamic Variables
            
            FieldClassi = "FieldClass = self.Field%d.__class__" %a                  # Create String to store current argument's Widget Type
            exec(FieldClassi)                                                       # Create Dynamic Variable to determine Tkinter Class
            FieldType = str(FieldClass).split('.')                                  # Split Variable into Class alone
            
            if FieldType[1] == "Entry" or FieldType[1] == "Scale":
               Passcount = Passcount+1
               runlist.append("runtmp%d = self.Field%d.get()" %(a,a))
            elif FieldType[1] == "Radiobutton":
               Passcount = Passcount+1
               radcount = radcount + 1
               self.RadioTypeList.append(a)
               runlist.append("runtmp%d = self.RadioVar%d.get()" %(a,radcount))
            elif FieldType[1] == "Checkbutton":
               tmpChkList = []
               checkcount = checkcount + 1
               statecount = statecount + 1
               DynChkList = "TmpChkList = map((lambda CheckVar%d: CheckVar%d.get()), self._states%d)"\
               %(checkcount,checkcount,statecount)                                  # Make Dynamic List holding boolean values representing selections
               DynChkList = DynChkList + "\nChkList%d = []" %(checkcount)
               exec(DynChkList)
               for b in range(len(TmpChkList)):
                  if TmpChkList[b] == '1':
                     runlist.append("ChkList%d.append(self.CheckChange%d[%d])"\
                     %(checkcount,a,b))                                             # Reassign Boolean Value 1 to proper variable
                  else:
                     runlist.append("ChkList%d.append(None)" %(checkcount))         # Reassign Boolean Value 0 to None
               runlist.append("runtmp%d = ChkList%d" %(a,checkcount))
               Passcount = Passcount+1
         for b in range(len(runlist)):
            exec(runlist[b])                                                        # Store Values User entered in static variable
         for c in range(Passcount):
            DynPassList.append("PassList.append(runtmp%d)" %c)
            exec(DynPassList[c])
         for d in range(len(PassList)):                                             # Change Entry Types for Check Button Lists
            if type(PassList[d]) == list:
               tmpCheckList = []
               for e in range(len(PassList[d])):
                  if type(PassList[d][e]) != NoneType:
                     if self.EntryTypeOrder[d] == 'BLANK':
                        tmpCheckList.append(PassList[d][e])
                     elif self.EntryTypeOrder[d] == 'f':
                        tmpCheckList.append(float(PassList[d][e]))
                     elif self.EntryTypeOrder[d] == 'i':
                        tmpCheckList.append(int(PassList[d][e]))
                     elif self.EntryTypeOrder[d] == 's':
                        tmpCheckList.append(str(PassList[d][e]))
                     elif self.EntryTypeOrder[d] == 'm':
                        MixOrderList = self.MixEntryTypeOrder[d].split(':')
                        if len(PassList[d]) != len(MixOrderList):
                           print "Incorrect Number of Arguments for Mixed Data Types!"
                           self.Quit()
                        if MixOrderList[e] == 'f':
                           tmpCheckList.append(float(PassList[d][e]))
                        elif MixOrderList[e] == 'i':
                           tmpCheckList.append(int(PassList[d][e]))
                        elif MixOrderList[e] == 's':
                           tmpCheckList.append(str(PassList[d][e]))
                        else:
                           print 'Unsupported Format Used in Mix Data Types!'
                           self.Quit()
               PassList[d] = tmpCheckList
         for f in range(self._argcnt):
            if self.StringNumber[f] == 'Multiple':                                  # Find Infiles entry to split into a proper list
               tmpFilesList = PassList[f].split('||')
               PassList[f] = tmpFilesList[0:len(tmpFilesList)-1]
            if PassList[f] == 'BLANK':                                              # Remove 'BLANK' so it is not passed to program
               PassList[f] = ''
         for g in range(self._argcnt):                                              # Change type for Entry Widgets if user indicates them
            FieldClassi = "FieldClass = self.Field%d.__class__" %g
            exec(FieldClassi)
            FieldType = str(FieldClass).split('.')
            
            if FieldType[1] == "Entry":
               if self.EntryTypeOrder[g] == 'BLANK':
                  pass
               elif self.EntryTypeOrder[g] == 'i':
                  PassList[g] = int(PassList[g])
               elif self.EntryTypeOrder[g] == 'f':
                  PassList[g] = float(PassList[g])
               elif self.EntryTypeOrder[g] == 's':
                  pass
               else:
                  print "Wrong variable format put into Entry Format!"
                  self.Quit()
         if len(self.RadioTypeList) > 0:                                            # Change type for Radio entries
            for z in range(len(self.RadioTypeList)):
               idx = self.RadioTypeList[z]
               if PassList[idx] != '':
                  if self.EntryTypeOrder[idx] == 'm':
                     ParamList = self.ParamlistOrder[idx].split(',')
                     MixList = self.MixEntryTypeOrder[idx].split(':')
                     if len(ParamList) != len(MixList):
                        print "Incorrect Number of Arguments for Mixed List!"
                        self.Quit()
                     for b in range(len(ParamList)):
                        if PassList[idx] == ParamList[b]:
                           if MixList[b] == 's':
                              TypeChangeRadio = PassList[idx]
                           elif MixList[b] == 'f':
                              TypeChangeRadio = float(PassList[idx])
                           elif MixList[b] == 'i':
                              TypeChangeRadio = int(PassList[idx])
                           else:
                              print "Wrong variable format put into mixed formating!"
                              self.Quit()
                  else:
                     TypeChangeRadio = PassList[idx]
                  PassList[idx] = TypeChangeRadio
         return PassList                                                            # Pass List of GUI Entries
     
     
      def Reload(self):                                                             # Reload current function
         """Reload the current function"""
         print "reload: %s" % self.Infile
         global ProgCnt
         ProgCnt = 1
         self.root.destroy()
   
      def Runbutton(self):                                                          # Create Function that returns values from GUI to needed program
         """Run Button Function"""
         FncList = self.GetValues()                                                 # Store List of entered GUI values
         FncRunStr = "self._fun("
         ProgramRunStr = "%s " %(self.ProgramName)
         if self.ProgramFlag == 1:
            if RunTypeFlag == 'System':
               for a in range(self._argcnt):
                  ProgramRunStr = ProgramRunStr + "%s " %FncList[a]
               print ProgramRunStr
               os.system(ProgramRunStr)
            elif RunTypeFlag == 'Miriad':
               for b in range(self._argcnt):
                  if FncList[b] != '':
                     ProgramRunStr = ProgramRunStr + "%s=%s " %(self.VarlistOrder[b],FncList[b])
                  else:
                     pass
               os.system(ProgramRunStr) 
         else:
            for b in range(self._argcnt):
               if b < (self._argcnt - 1):
                  FncRunStr = FncRunStr + "FncList[%d]," %b
               else:
                  FncRunStr = FncRunStr + "FncList[%d])" %b
            exec(FncRunStr)                                                         # Execute Function Called by User with GUI values
   
     
      def Quit(self):                                                               # Quit Button Function
         self.root.quit()
     
           
      def helper(self):
         """Displays a Help Window compiled from the Doc String"""
         helpwindowtop = Toplevel()                                                 # Create a new widget
         helpwindowtop.geometry("595x330")
         helpwindowtop.title("Variable Help")                                       # Name New Window
         helpwindow = Frame(helpwindowtop)
         helpwindow.grid_rowconfigure(0, weight=0)
         helpwindow.grid_columnconfigure(0, weight=0)
         text=Text(helpwindow)
         for b in range(len(self._HELPList)):
            FormatHELP = self._HELPList[b].split()                                  # Split Current Line in List of Help Strings For Formatting Reasons
            HelpString = "%s " %str(FormatHELP[0])
            for c in range(len(FormatHELP) - 1):
               HelpString = HelpString + "%s " %str(FormatHELP[c+1])
            text.insert(END,"%s \n" %HelpString)                                    # Print Current List Item in Left Justified Format
        
         scroll = Scrollbar(helpwindow)
         text.config(state=DISABLED,yscrollcommand=scroll.set)
         scroll.config(command=text.yview)
         text.pack(side=LEFT)
         scroll.pack(side=LEFT,fill=Y)
         helpwindow.pack(anchor=N+W,fill=BOTH)
   
      def Check(self):
         """Create a Window with some Useful Checks"""
         ArgFlag=0                                                                  # A Flag to alert user about defining variables in a call to function
         DocFlag=0                                                                  # A flag to alert user that Document string does not give information for each variable
         CheckWindowtop = Toplevel()                                                # Make a new widget
         CheckWindow = Frame(CheckWindowtop)
         CheckWindow.grid_rowconfigure(0, weight=0)
         CheckWindow.grid_columnconfigure(0, weight=0)
         CheckWindowtop.title("Check of Docstring and Arguments")                   # Name the new Window
         CheckWindowtop.geometry("595x330")
         text=Text(CheckWindow)
              
         if self._fun.func_defaults != None:
            text.insert(END,"Default Arguments Passed to Function:\n")
            if len(self._fun.func_code.co_varnames) == len(self._fun.func_defaults):   # Check if Each Argument has been defined in call to function
               for a in range(len(self._fun.func_code.co_varnames)):
                  text.insert(END,"%2s= %s \n" %(self._fun.func_code.co_varnames[a],\
                  self._fun.func_defaults[a]))
            else:                                                                      # If Each Argument is not Defined, then Display Values found by Docstring
               text.insert(END,"Note: Not all arguments defined in pass to "+\
               "function.Missing values provided by docstring.\n")
               ArgFlag=1
               for b in range(len(self._fun.func_code.co_varnames)):
                  text.insert(END,"%2s= %s \n" %(self._fun.func_code.co_varnames[b],\
                  self.ValuelistOrder[b]))
         else:
            ArgFlag=1 
         
         text.insert(END,"\nArguments Quoted by Docstring:\n")
         for c in range(len(self.VarlistOrder)):
            if self.ValuelistOrder[c] != 'BLANK':
               text.insert(END,"%2s= %s \n" %(self.VarlistOrder[c],self.ValuelistOrder[c]))
            else:
               text.insert(END,"%2s= (Not Defined in Doc String) \n" %(self._fun.\
               func_code.co_varnames[c]))
               DocFlag=1
         
         text.insert(END,"\nArguments To be passed to function by GUI:\n")
         ArgCheckList = self.GetValues()                                            # Get Current Values in the GUI
         for d in range(len(self._fun.func_code.co_varnames)):
            text.insert(END,"%2s= %s \n" %(self._fun.func_code.co_varnames[d],ArgCheckList[d]))
         if ArgFlag == 1:
            text.insert(END,"\nFor best results, please pass defaults values to each argument.")
         if DocFlag == 1:
            text.insert(END,"\nFor best results, please write a complete document string.")
         
         scroll = Scrollbar(CheckWindow)
         text.config(state=DISABLED,yscrollcommand=scroll.set)
         scroll.config(command=text.yview)
         scroll.grid(row=0, column=0, columnspan=2)
         text.pack(side=LEFT)
         scroll.pack(side=LEFT,fill=Y)
         CheckWindow.pack(anchor=N+W,fill=BOTH)
        
      def tkrun(self):                                                              # Function to Create GUI
         global ProgCnt
         ProgCnt = 0
   
         DISP        = []                                                           # Create List of Dynamic Tkinter Classes depending on input function
         typelist    = []                                                           # List of relevant lines from the help document
         Fieldlist   = []                                                           # List of Tkinter Class for each variable
         Varlist     = []                                                           # Store which variable each line corresponds to
         Paramlist   = []                                                           # Store list of Parameters used for certain Tkinter Classes
         Valuelist   = []                                                           # Store value in document string for each variable
         EntryType   = []                                                           # Store User defined Variable Type
         MixEntryType = []                                                          # Store User defined Mixed Variable Types
         self.FieldlistOrder   = []                                                 # Reordered Fieldlist to represent how variables are called in function
         self.VarlistOrder     = []                                                 # Reordered varlist to represent how variables are called in function
         self.ParamlistOrder   = []                                                 # Reordered Paramlist to represent how variables are called in function
         self.ValuelistOrder   = []                                                 # Reordered Valuelist to represent how variables are called in function
         self.EntryTypeOrder   = []                                                 # Reordered EntryType to represent how variables are called in function
         self.MixEntryTypeOrder   = []                                              # Reordered MixEntryType to represent how variables are called in function
         
         if self._fun.__doc__ == None:
            print "ERROR: Missing Document String"
            self.Quit()
         document = self._fun.__doc__                                               # Read in Document String for Function
         if document[0][0] == '!':
            inputflag = 1
         else:
            inputflag = 0
         doclist = document.split('\n')
         for a in range(len(doclist)):
            SpltDoci = doclist[a].split('#> ')
            if len(SpltDoci) > 1:
               typelist.append(SpltDoci[1])                                         # Store Relevant lines from Document String
            else:
               self._HELPList.append(SpltDoci[0])
         for b in range(len(typelist)):
            typelisti = typelist[b].split(' ')
            typelistj = typelisti[1].split('=')
            Varlist.append(typelistj[0])
           
            if len(typelisti[0].split('.')) == 1:
               Fieldlist.append(typelisti[0])
               EntryType.append('')
            else:
               Fieldlist.append(typelisti[0].split('.')[0])
               EntryType.append(typelisti[0].split('.')[1])
            
            if len(typelistj) > 1:                                                  # Store a value for each variable depending on Tkinter class
               Valuelist.append(typelistj[1])
            else:
               Valuelist.append('')
            
            if len(typelisti) > 2:                                                  # Store Parameters for certain Tkinter Classes
               Paramlist.append(typelisti[2])
            else:
               Paramlist.append('')
   
            if EntryType[b] == 'm':
               MixEntryType.append(typelisti[3])
            else:
               MixEntryType.append('')
       
         for c in range(self._fun.func_code.co_argcount):                           # Begin loop to sort document information into order it appears in function
            flag = 0
            self.FieldlistOrder.append('')
            self.VarlistOrder.append('')
            self.ParamlistOrder.append('')
            self.ValuelistOrder.append('')
            self.EntryTypeOrder.append('')
            self.MixEntryTypeOrder.append('')
            
            for d in range(len(Varlist)):
               if self._fun.func_code.co_varnames[c] == Varlist[d]:
                  self.VarlistOrder[c]   = Varlist[d]
                  self.FieldlistOrder[c] = Fieldlist[d]
                  self.ParamlistOrder[c] = Paramlist[d]
                  
                  if Valuelist[d] != '':
                     self.ValuelistOrder[c] = Valuelist[d]
                  else:
                     self.ValuelistOrder[c] = 'BLANK'
                 
                  if EntryType[d] != '':
                     self.EntryTypeOrder[c] = EntryType[d]
                  else:
                     self.EntryTypeOrder[c] = 'BLANK'
   
                  if MixEntryType[d] != '':
                     self.MixEntryTypeOrder[c] = MixEntryType[d]
                  else:
                     self.MixEntryTypeOrder[c] = 'BLANK'
                    
                  flag = 1
               elif self._fun.func_code.co_varnames[c] != Varlist[d] and flag == 0:
                  self.VarlistOrder[c]      = 'BLANK'
                  self.FieldlistOrder[c]    = 'BLANK'
                  self.ParamlistOrder[c]    = 'BLANK'
                  self.ValuelistOrder[c]    = 'BLANK'
                  self.EntryTypeOrder[c]    = 'BLANK'
                  self.MixEntryTypeOrder[c] = 'BLANK'              
               else:
                  pass
         for e in range(self._fun.func_code.co_argcount):
            if self.VarlistOrder[e] == 'BLANK':
               self.VarlistOrder[e] = self._fun.func_code.co_varnames[e]
            else:
               pass
   
         paramlength = 0
         extrarow    = 0
         for e in range(len(self.ParamlistOrder)):
            paramtmp = self.ParamlistOrder[e].split(',')
            if len(paramtmp) > paramlength:
               paramlength = len(paramtmp)                                       # Store the number maximum number of parameters for one variable
            else:
               pass
            if len(paramtmp) > 13:
               extrarow = extrarow + long(len(paramtmp)/13)                      # Find how many extra lines on the window are needed for extra rows of parameters
         
         if (self._fun.func_code.co_argcount+extrarow) >= 17:                    # Set default window height if over sMaximum number of lines for a window.
            WINheight = 1000
         else:
            WINheight = (50+40*self._fun.func_code.co_argcount+40*extrarow)      # Adjust Height if under maximum.
         
         if paramlength >= 13:                                                   # Set default window width if over sMaximum number of lines for a window.
            WINwidth = 1250
         else:
            if paramlength <= 5:                                                 # Set minimum window width
               WINwidth = 555
            else:
               WINwidth = (515+75*(paramlength-5))                               # Set intermediate window widths
         
         root = Tk()                                                             # Create Parent Window
         self.root = root
         if self.ProgramFlag == 0:
            root.title("%s" %self._fun.__name__)                                 # Name Window with Function Name
         else:
            root.title("%s" %self.ProgramName)
         root.geometry("%dx%d" %(WINwidth,WINheight))                            # Change Size of GUI to represent amount of variables
   
         menu = Menu(root,relief=RIDGE)                                          # Create a menu
         root.config(menu=menu)
         filemenu = Menu(menu)
         helpmenu = Menu(menu)
         menu.add_cascade(label="File", menu=filemenu)
         filemenu.add_command(label="Check", command=self.Check)
         filemenu.add_command(label="Reload", command=self.Reload)
         filemenu.add_command(label="Quit", command=self.Quit)
         menu.add_command(label="Run", command=self.Runbutton)
         menu.add_command(label="Help", command=self.helper)
        
         frame = Frame(root, bd=0)                                               # Create Frame to Pack Canvas and Scrollbar that will hold widgets
         frame.grid_rowconfigure(0, weight=0)                                    # Maintain Proper Dimensions with Frame
         frame.grid_columnconfigure(0, weight=0)                                 # Maintain Proper Dimensions with Frame
         canvas = Canvas(frame)
         scrollbar = Scrollbar(frame)
         canvas.config(yscrollcommand=scrollbar.set)
         canvas['background'] = 'gray80'
         canvas['width'] = (WINwidth-20)
         canvas['height'] =(50+40*self._fun.func_code.co_argcount+40*extrarow)
         canvas['scrollregion'] = (0, 0, 0, (60+40*self._fun.func_code.co_argcount+40*extrarow))
         canvas.pack(side=LEFT)
         scrollbar.config(command=canvas.yview)
         scrollbar.pack(side=LEFT,fill=Y)
         frame.pack(anchor=N+W,fill=BOTH)
        
        
         if self._fun.func_defaults != None:
            if len(self._fun.func_defaults) == (self._fun.func_code.co_argcount) and inputflag == 0: #Determine where to grab default entries from
               DefaultSource = 1
            else:
               DefaultSource=0
         else:
            DefaultSource=0
   
         for a in range(self._fun.func_code.co_argcount):
            textvar = "tmp%d = StringVar()\n" %a
            if DefaultSource == 1:
               textvar = textvar + "tmp%d.set(%s.%s.func_defaults[%d])" %(a,self._fun.__module__,self._fun.__name__,a)
            else:
               textvar = textvar + "tmp%d.set(self.ValuelistOrder[%d])" %(a,a)
                    
            exec(textvar)                                                                 # Make Static Version of each Dynamic Variable
            
            DISP.append("self.Label%d = Label(frame,text='%s:',relief=RIDGE)"\
            %(a,self._fun.func_code.co_varnames[a]))   # Create Label for each Variable
            DISP.append("canvas.create_window(28, (15+%d), width=51, height=20,window=self.Label%d,anchor=CENTER)" %((a*40+40*self._place),a))
           
            if self.FieldlistOrder[a] == "ENTRY":                                         # Make Entry Class if quoted in Docstring for current Variable
               DISP.append("self.Field%d = Entry(root,textvariable=\
               tmp%d,width=50,bd=1,background='White',selectbackground='Black')" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == "SCALE":                                       # Make Scale Class if quoted in Docstring for current Variable
               Scalelist = self.ParamlistOrder[a].split(':')
               
               DISP.append("self.Field%d = Scale(root,from_=%d,to=%d,resolution=%f,\
               orient=HORIZONTAL,length=350,width=9,font=8,variable=tmp%d)"\
               %(a,float(Scalelist[0]),float(Scalelist[1]),float(Scalelist[2]),a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=40,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == "CHECK":                                       # Make Check Class if quoted in Docstring for current Variable
               Checklist = self.ParamlistOrder[a].split(',')
               self._checkcount = self._checkcount+1
               self._statecount = self._statecount+1
               
               DynCheck = "self.CheckChange%d = Checklist" %a
               DynStates = "self._states%d = []" %self._statecount
               exec(DynCheck)
               exec(DynStates)
               
               for b in range(len(Checklist)):
                  DISP.append("self.CheckVar%d = StringVar()" %(self._checkcount))
                  
                  DISP.append("self.Field%d = Checkbutton(root, text='%s',variable=self.CheckVar%d)" %(a,str(Checklist[b]),self._checkcount))
                  DISP.append("canvas.create_window(%d, (15+%d), width=60, height=20,window=self.Field%d,anchor=CENTER)"\
                  %((150+long((b-13*long(b/13))*75)),(a*40+long(b/13)*40+40*self._place),a))
                  DISP.append("self._states%d.append(self.CheckVar%d)" %(self._statecount,self._checkcount))
               
               if len(Checklist) >= 13:
                  self._place=self._place+long(len(Checklist)/13)
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == "RADIO":                                       # Make Radio Class if quoted in Docstring for current Variable
               Radiolist = self.ParamlistOrder[a].split(',')
               self._radcount = self._radcount+1
               
               for c in range(len(Radiolist)):
                  if self.EntryTypeOrder[a] == 'BLANK':
                     DynRadioVar = "self.RadioVar%d = StringVar()" %(self._radcount)
                  elif self.EntryTypeOrder[a] == 'f':
                     DynRadioVar = "self.RadioVar%d = DoubleVar()" %(self._radcount)
                  elif self.EntryTypeOrder[a] == 'i':
                     DynRadioVar = "self.RadioVar%d = IntVar()" %(self._radcount)
                  elif self.EntryTypeOrder[a] == 's':
                     DynRadioVar = "self.RadioVar%d = StringVar()" %(self._radcount)
                  elif self.EntryTypeOrder[a] == 'm':
                     DynRadioVar = "self.RadioVar%d = StringVar()" %(self._radcount)
                  else:
                     print 'Unsupported Format Used in Mix Data Types!'
                     self.Quit()
                  exec(DynRadioVar)
   
                  DISP.append("self.Field%d = Radiobutton(root, text='%s', value='%s',\
                  variable=self.RadioVar%d)" %(a,str(Radiolist[c]),Radiolist[c],\
                  self._radcount))
                  DISP.append("canvas.create_window(%d, (15+%d), width=60, height=20,window=self.Field%d,anchor=CENTER)"\
                  %((150+long((c-13*long(c/13))*75)),(a*40+long(c/13)*40+40*self._place),a))
               
               if len(Radiolist) >= 13:
                  self._place=self._place+long(len(Radiolist)/13)
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == "INFILE":                                      # Make Infile if quoted in Docstring for current Variable
               DynInFile = "self.InFile%d = StringVar()" %(self.Ifieldentry)
               exec(DynInFile)
               InFileSTR = ""
               InFileSTR = InFileSTR + "def GetFile%d(self,event):\n" %a
               InFileSTR = InFileSTR + "   currentdir = os.getcwd()\n"
               InFileSTR = InFileSTR + "   filename = askopenfilename(initialdir=currentdir)\n"
               InFileSTR = InFileSTR + "   if filename:\n"
               InFileSTR = InFileSTR + "      self.InFile%d.set(filename)\n" %(self.Ifieldentry)
               InFileSTR = InFileSTR + "Test.PyRun.GetFile%d = GetFile%d\n" %(a,a)
               InFileSTR = InFileSTR + "self.InFile%d = tmp%d" %(self.Ifieldentry,a)
               exec(InFileSTR)
               
               DISP.append("self.Field%d = Entry(root,textvariable=self.InFile%d,width=50,bd=1,background='White',selectbackground='Black')"%(a,self.Ifieldentry))
               DISP.append("self.Field%d.bind('<Button-1>', self.GetFile%d)" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               DISP.append("self.Label%d_%d = Label(root,text='In',relief=RIDGE)" %(a,a))
               DISP.append("canvas.create_window(470, (15+%d), width=25, height=20,window=self.Label%d_%d,anchor=CENTER)" %((a*40+40*self._place),a,a))
               self.Ifieldentry = self.Ifieldentry+1
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == "INFILES":                                      # Make Infile if quoted in Docstring for current Variable
               DynInFiles = "self.InFiles%d = StringVar()" %(self.Ifieldentries)
               exec(DynInFiles)
               InFilesSTR = ""
               InFilesSTR = InFilesSTR + "def GetFiles%d(self,event):\n" %a
               InFilesSTR = InFilesSTR + "   currentdir = os.getcwd()\n"
               InFilesSTR = InFilesSTR + "   filenames = askopenfilenames(initialdir=currentdir)\n"
               InFilesSTR = InFilesSTR + "   if filenames:\n"
               InFilesSTR = InFilesSTR + "      filenamesstr = ''\n"
               InFilesSTR = InFilesSTR + "      for a in range(len(filenames)):\n"
               InFilesSTR = InFilesSTR + "         filenamesstr = filenamesstr + filenames[a] + '||'\n"
               InFilesSTR = InFilesSTR + "      self.InFiles%d.set(filenamesstr)\n" %(self.Ifieldentries)
               InFilesSTR = InFilesSTR + "Test.PyRun.GetFiles%d = GetFiles%d\n" %(a,a)
               InFilesSTR = InFilesSTR + "self.InFiles%d = tmp%d" %(self.Ifieldentries,a)
               exec(InFilesSTR)
              
               DISP.append("self.Field%d = Entry(root,textvariable=self.InFiles%d,width=50,bd=1,background='White',selectbackground='Black')"%(a,self.Ifieldentries))
               DISP.append("self.Field%d.bind('<Button-1>', self.GetFiles%d)" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               DISP.append("self.Label%d_%d = Label(root,text='InFiles',relief=RIDGE)" %(a,a))
               DISP.append("canvas.create_window(482, (15+%d), width=40, height=20,window=self.Label%d_%d,anchor=CENTER)" %((a*40+40*self._place),a,a))
               self.Ifieldentries = self.Ifieldentries+1
               self.StringNumber.append('Multiple')
            elif self.FieldlistOrder[a] == "OUTFILE":                                     # Make Outfile if quoted in Docstring for current Variable
               DynOutFile = "self.OutFile%d = StringVar()" %(self.Ofieldentry)
               exec(DynOutFile)
               OutFileStr = ""
               OutFileStr = OutFileStr + "def SaveFile%d(self,event):\n" %a
               OutFileStr = OutFileStr + "   currentdir = os.getcwd()\n"
               OutFileStr = OutFileStr + "   filename = asksaveasfilename(initialdir=currentdir)\n"
               OutFileStr = OutFileStr + "   if filename:\n"
               OutFileStr = OutFileStr + "      self.OutFile%d.set(filename)\n" %(self.Ofieldentry)
               OutFileStr = OutFileStr + "Test.PyRun.SaveFile%d = SaveFile%d\n" %(a,a)
               OutFileStr = OutFileStr + "self.OutFile%d = tmp%d" %(self.Ofieldentry,a)
               exec(OutFileStr)
               
               DISP.append("self.Field%d = Entry(root,textvariable=self.OutFile%d,width=50,bd=1,background='White',selectbackground='Black')" %(a,self.Ofieldentry))
               DISP.append("self.Field%d.bind('<Button-1>', self.SaveFile%d)" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               DISP.append("self.Label%d_%d = Label(root,text='Out',relief=RIDGE)" %(a,a))
               DISP.append("canvas.create_window(472, (15+%d), width=25, height=20,window=self.Label%d_%d,anchor=CENTER)" %((a*40+40*self._place),a,a))
               self.Ofieldentry = self.Ofieldentry+1
               self.StringNumber.append('')
            elif self.FieldlistOrder[a] == 'OPENDIR':
               DynDirFile = "self.DirFile%d = StringVar()" %(self.DirFieldentry)
               exec(DynDirFile)
               DirFileSTR = ""
               DirFileSTR = DirFileSTR + "def GetDir%d(self,event):\n" %a
               DirFileSTR = DirFileSTR + "   currentdir = os.getcwd()\n"
               DirFileSTR = DirFileSTR + "   filename = askdirectory(initialdir=currentdir)\n"
               DirFileSTR = DirFileSTR + "   if filename:\n"
               DirFileSTR = DirFileSTR + "      self.DirFile%d.set(filename)\n" %(self.DirFieldentry)
               DirFileSTR = DirFileSTR + "Test.PyRun.GetDir%d = GetDir%d\n" %(a,a)
               DirFileSTR = DirFileSTR + "self.DirFile%d = tmp%d" %(self.DirFieldentry,a)
               exec(DirFileSTR)
               
               DISP.append("self.Field%d = Entry(root,textvariable=self.DirFile%d,width=50,bd=1,background='White',selectbackground='Black')"%(a,self.DirFieldentry))
               DISP.append("self.Field%d.bind('<Button-1>', self.GetDir%d)" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               DISP.append("self.Label%d_%d = Label(root,text='Dir',relief=RIDGE)" %(a,a))
               DISP.append("canvas.create_window(472, (15+%d), width=25, height=20,window=self.Label%d_%d,anchor=CENTER)" %((a*40+40*self._place),a,a))
               self.DirFieldentry = self.DirFieldentry+1
               self.StringNumber.append('')
            else:                                                                         # Default to Entry Class if not in Docstring
               DISP.append("self.Field%d = Entry(root,textvariable=tmp%d,width=50,bd=1,background='White',selectbackground='Black')" %(a,a))
               DISP.append("canvas.create_window(275, (15+%d), width=350, height=20,window=self.Field%d,anchor=CENTER)" %((a*40+40*self._place),a))
               self.StringNumber.append('')
         for d in range(len(DISP)):
            exec(DISP[d])
   
   
         root.mainloop()
   
   class Begin:
      def __init__(self,arg):
         self.Sys1 = arg
      def ProgramHelp(self):
         print "Code Designed to Automatically Create a GUI for functions with the proper doc string."
         print "Example Usage: GUIBuild.py pyrunpull.foo [Where pyrunpull is a program, and foo is a function.]"
         print "Example Document String Format:\n\n"
         print "#> ENTRY a=1"
         print "a is an example Entry Widget which may default to 1 if no value is passed through the function. It can also take '.i' and '.f' formating"+\
               "to change the type to integer or float."
         print "#> SCALE b=2 0:10:.1"
         print "b is an example Scale Widget. It may default to 2 if no value is passed through the function, and it spans 0 to 10 in increments of .1."
         print "#> CHECK.m c= String,1,3.14 s:i:f"
         print "c is an example of the Checkbutton Widget. It has no default, but will display a list of options (String,1, or 3.14). CHECK may be called"+\
               "'.m','.s','.f','.i', or nothing at all. If nothing is provided, the entries will default to string format. '.s' will also provide string"+\
               "format, whereas '.f' and '.i' correspond to float and integer. '.m' indicates a mixed type check list. If '.m' is used, all entries must"+\
               "be designated a type, and this is to be placed into the format as shown above. One space away from the list of entry values, and put into"+\
               "a list of datatypes separated by a colon. '.m' is the only format that strictly requires the list of data types to be present"
         print "\n#> Tells the program where to look for information needed to build the GUI. Any line without this will be parsed into the help menu.\n"
         print "Additional Widgets:"
         print "   INFILE : Allows user to right click an entry form to bring up file browser to choose single file to read."
         print "   OUTFILE : Allows user to right click an entry form to bring up file browser to choose single file to write to."
         print "   OPENDIR : Allows user to right click and bring up a directory browser to choose a directory."
         print "   INFILES : Allows user to right click and bring up a file browser where they may choose multiple files. Files will display as a string"+\
               "with files separated by a '||' to be parsed by the GUI before passing it to the called function as a list."
         print "   RADIO : A radio checklist. It takes the similar formating options as CHECK."
   
      
      def Grab(self):
     
         inFileFull  = self.Sys1
         global ProgCnt
         ProgCnt = 1
      
         if inFileFull == "--help":
            self.ProgramHelp()
         else:
            ProgramList = inFileFull.split(':')
   
            for ProgIdx in range(len(ProgramList)):
               while ProgCnt == 1:
                  inFileSplt  = ProgramList[ProgIdx].split('.')
                  inFile      = inFileSplt[0]
      
                  importstr = "import %s" %inFile
                  exec(importstr)
         
                  if len(inFileSplt) == 1:
                     LengthFlag = 1
                  else:
                     LengthFlag = 0
   
                  if len(sys.argv) == 3 and sys.argv[2] == 'MirRUN' or len(sys.argv) == 3 and sys.argv[2] == 'SysRUN':
                     ProgramFlag = 1
                     if sys.argv[2] == 'MirRUN':
                        RunTypeFlag = 'Miriad'
                     else:
                        RunTypeFlag = 'System'
                  else:
                     ProgramFlag = 0
   
                  command = "p = Test.PyRun(%s,%s,%s,%s)\n" %(ProgramList[ProgIdx],inFile,LengthFlag,ProgramFlag)
                  command = command + "p.tkrun()"
                  #print "DEBUG: exec(%s)" % command
                  exec(command)
                  if ProgCnt == 1:
                     reloadstr = "reload(%s)" %inFile
                     exec(reloadstr)
   
            print "All done with %s" % ProgramList

from Tkinter import *
from tkFileDialog import askopenfilename
from tkFileDialog import askopenfilenames
from tkFileDialog import asksaveasfilename
from tkFileDialog import askdirectory
import os,sys

if __name__ == '__main__':
   passer = sys.argv[1]
   a = Test.Begin(passer)
   a.Grab()
