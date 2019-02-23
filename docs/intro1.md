# A birds eye view of NEMO

## Installing

To use it (this assumes somebody has installed it for you)

   source /somewhere/nemo/nemo_start.sh           (or .csh)

To install it in a bash shell:

   wget https://teuben.github.io/nemo/install_nemo
   chmod +x install_nemo
   ./install_nemo 
   source nemo/nemo_start.sh

## Using

The source command will have modified your shell and added a number
of new commands to your unix shell. Much like other unix commands
these NEMO commands will:

      * have a unix man page, e.g.          man mkplummer
      * have a --help option, e.g.          mkplummer --help
      * have a help= system keyword, e.g.   mkplummer help=h
        to describe the program keywords
      * use unix pipes                      mkplummer - 10 | tsf -


## Example 1: 
