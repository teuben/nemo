
CRCK(1)             UNIX Programmer's Manual              CRCK(1)

NAME
     crck (version 2) - CRC value of files

SYNOPSIS
     crck [-t|u|c|i] [filename ... ]

OPTIONS
     -t   Perform CRCK on a CP/M text file that is stored as a Unix file.
          CP/M record delimiters (CR/LF) are substituted for Unix
          newlines and the file is padded to a 128-byte sector boundary
          with ^Z (0x1A).

     -u   Perform CRCK on a regular Unix file.  Takes whatever is there.
          This is the default option.

     -c   Perform CRCK on a CP/M non-text file (such as .COM or .LBR).
          Handled the same as the -u flag except 128-byte sectors are
          reported rather than bytes.

     -i   Perform CRCK on a CP/M non-text file stored in ITS binary
          format. Ignores the first four bytes of the file which contain
          an ITS code.

     Wildcards may, of course, be used and the program may be
     used in pipelines.

DESCRIPTION
     Crck calculates and prints a 16-bit CRC value and either a
     CP/M 128-byte sector count or the byte count for the named
     files.  It is typically used to validate files transferred
     between Unix and CP/M type operating systems.  Crck uses a
     CRC compatible with the CP/M and CP/M-86 versions of crck
     prior to 5.0, as well as the crck command in 8 bit versions
     of YAM (Yet Another Modem program).

     When comparing copies of a text file on Unix and CP/M or
     similar systems, account must be taken of the format of the
     file on the Unix system to determine which option to use to
     check the file. The file may be in CP/M format, with car-
     riage return/line feed line terminators and padding with
     control-Z to make the length divisible by 128 (-u option);
     or it may be in Unix format, with newline line delimiters
     and no padding (-t option); or in ITS-binary format (-i
     option).

     Since CP/M binary files are not 'translated' to Unix, these
     may be checked directly. Note, though, that any binary file
     to be sent to a CP/M system must have a length of a multiple
     of 128 bytes to be checked correctly on both ends. If it
     does not, the CP/M system will 'fill' the rest of the last
     128-byte record, and the CRC checks will be different, even
     though the data transferred correctly.

DIAGNOSTICS
     Prints proper usage if an incorrect flag is given.

BUGS
     Has not been rigorously tested.  May not produce correct
     results with all processor types.  This CRC detects all sin-
     gle bit errors and most, but not all, multiple bit errors.
     This CRC does not detect all two bit errors, unlike the CRC
     used with standard synchronous communications protocols.

WARNING
     Files with different data may still produce the same crck
     value.
