#! /usr/bin/env python
# Author: Parker Tewell
# ptewell@terpmail.umd.edu
# 
# apr 15 2021: First Version

import re, os
from sys import breakpointhook

def get_man_matches(file):
    try:
        man_doc = open("man/man1/"+file+".1", "r") # Read man file
    except:
        return None

    # Get man matches
    scan_flag = False
    man_matches = []

    for line in man_doc.readlines():
        if re.search(r'.TP',line): # If we encounter a .TP, scan next line for a command
            scan_flag = True
        elif scan_flag:
            match = re.findall(r'\\f[a-zA-Z][a-zA-Z]+[0-9]*',line)
            if not match: # If the .TP isn't followed by a command, flag file as bad
                return None
            else:
                man_matches.append(match[0])
            scan_flag = False
        
    # Transform man data for comparison
    man_matches = [match[3:] for match in man_matches]

    man_doc.close()
    return man_matches

def get_help_matches(file):
    try:
        help_out = os.popen(file+" help=h").read() # Read help output
    except:
        return None

    # Grab keywords
    help_matches = []
    for line in help_out.splitlines():
        keyword = line.split()[0]
        help_matches.append(keyword) if keyword != "VERSION" else breakpointhook
    
    return help_matches

def checkMan():
    files_read, files_read_names = 0, []
    bad_files, bad_file_names = 0,[]
    tasklist = open("src/scripts/tasklist")

    # Iterate over tasklist
    for file in tasklist:
        file = file.split()[0]

        # Check if the commands and the order they appear in match
        # Using a try catch in case some files have no matching man
        man_out = get_man_matches(file)
        help_out = get_help_matches(file)

        # Flag a file as bad if their keywords don't match or if one of then unsuccesfully ran
        if man_out != help_out or not(man_out and help_out):
            bad_files+=1
            bad_file_names.append(file)

        files_read+=1
        files_read_names.append(file)
    
    print("Files with help=h & man mismatches")
    print("Files read: " + str(files_read))
    print("Bad files found: " + str(bad_files))
    print("Bad files: " + str(bad_file_names))

    tasklist.close()
    return files_read_names

def checkBin(tasklist):
    files_read = 0
    bad_files, bad_file_names = 0, []
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

def main():
    # Change working directory
    os.chdir("..")
    os.chdir("..")

    # Run tests
    tasklist = checkMan()
    print()
    checkBin(tasklist)

if __name__ == "__main__":
    main()
