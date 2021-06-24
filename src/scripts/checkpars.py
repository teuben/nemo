#! /usr/bin/env python
# Author: Parker Tewell
# ptewell@terpmail.umd.edu
# 
# apr 15 2021: First Version
# apr 30 2021: Stopped flagging files w/o help=h
#              Added -v verbose flag
# may 6  2021: -h flag added
# may 29 2021: -f flag added

import re, os, subprocess, getopt, sys

# Global flags
VERBOSE = False # If True prints man & help outputs for bad files
HELP = False
TASKLIST = "src/scripts/tasklist"

def get_man_matches(file):
    try:
        man_doc = open("man/man1/"+file+".1", "r") # Read man file
    except:
        return None

    # Get man matches
    scan_flag = False
    man_matches = []

    for line in man_doc.readlines():
        if re.search(r'\.TP',line): # If we encounter a .TP, scan next line for a command
            scan_flag = True
        elif scan_flag: #\\fB(\w|#)*= \\f[a-zA-Z][a-zA-Z]+[0-9]*
            match = re.findall(r'\\fB[\w|#|/]*=',line)
            if not match: # If the .TP isn't followed by a command, flag file as bad
                return 'Non-conformant: ' + line
            else:
                man_matches.append(match[0][:-1])
            scan_flag = False
        
    # Transform man data for comparison
    man_matches = [match[3:] for match in man_matches]

    man_doc.close()
    return man_matches

def get_help_matches(file):
    help_out = os.popen(file+" help=h").read() # Read help output

    # Grab keywords
    help_matches = []
    for line in help_out.splitlines():
        if(line[0] != ' '): # Make sure the line doesn't start with a space
            keyword = line.split()[0]
            if keyword != "VERSION":
                help_matches.append(keyword) 
            else: 
                break
    
    return help_matches

def help_exists(file):
    return True if "VERSION" in subprocess.getoutput(file+" help=h") else False

def checkMan():
    files_read, files_read_names = 0, []
    bad_files, bad_file_names = 0,[]
    tasklist = open(TASKLIST)

    # Iterate over tasklist
    for file in tasklist:
        file = file.split()[0]

        # Ignoring files with no help=h output
        if help_exists(file):
            files_read+=1
            files_read_names.append(file)

            # Check if the commands and the order they appear in match
            # Using a try catch in case some files have no matching man
            man_out = get_man_matches(file)
            help_out = get_help_matches(file)

            # Flag a file as bad if their keywords don't match or if one of then unsuccesfully ran
            if man_out != help_out or not(man_out and help_out):
                bad_files+=1
                bad_file_names.append(file)

                if VERBOSE:
                    print(file)
                    print("man: "+str(man_out))
                    print("bin: "+str(help_out))
                    print() 
    
    print("Files with help=h & man mismatches")
    print("Files read: " + str(files_read))
    print("Bad files found: " + str(bad_files))
    print("Bad files: " + str(bad_file_names))

    tasklist.close()
    return files_read_names

def checkBin(tasklist):
    files_read, bad_files, bad_file_names = 0,0,[]
    binlist = open("src/scripts/bins.list")

    for file in binlist:
        file = file.rstrip()
        if not file in tasklist:
            bad_files+=1
            bad_file_names.append(file)
        files_read+=1
        
    print("Bin files with no matching tasklist file:")
    print("Files read: " + str(files_read))
    print("Bad files found: " + str(bad_files))
    print("Bad files: " + str(bad_file_names))

def readFlags():
    global VERBOSE, HELP, TASKLIST
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv,"vhf:")
        for opt, arg in opts:
            if opt in ["-v"]:
                VERBOSE = True
            elif opt in ["-h"]:
                HELP = True
            elif opt in ["-f"]:
                TASKLIST = arg
    except:
        print("getopt error")

def help():
    print("checkpars V1.3")
    print("-v           prints help and man keywords for mismatched files")
    print("-f <file>    allows you to specify tasklist to scan")
    print("-h           prints help page")

def main():
    readFlags()

    if HELP: # Prints help and exits
        help()
        return

    # Change working directory
    os.chdir(os.environ['NEMO'])

    # Run tests
    checkMan()

if __name__ == "__main__":
    main()
