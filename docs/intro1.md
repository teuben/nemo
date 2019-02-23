# A birds eye view of NEMO

## Installing

To use it (*this assumes somebody has installed it for you*) in csh shell:

       source /somewhere/nemo/nemo_start.csh

To install it in a bash shell:

       wget https://teuben.github.io/nemo/install_nemo
       chmod +x install_nemo
       ./install_nemo 
       source nemo/nemo_start.sh

The source command can then be added to your **.cshrc** (csh shell) or
**.bashrc** (linux bash) or **.bash_profile** (mac bash)

## Using

Your unix shell will now have been modified and a large number
of new commands are available. Much like other unix commands
these NEMO commands will:
    

      * live in $NEMOBIN, e.g.                    ls $NEMOBIN
      * have a unix man page, e.g.                man mkplummer
      * have a --help option, e.g.                mkplummer --help
      * have a help= system keyword, e.g.         mkplummer help=\?
        and describe the program keywords, e.g.   mkplummer help=h
      * use unix pipes                            mkplummer - 10 | tsf -


## Example 1: 
