#! /usr/bin/env python
#
# docopt


"""Example11

Usage:
  ex11.py [options]

Options:
  -h --help       Show this help
  --version       Show version
  -a --aaa=A      Set the A value [Default: 1]
  -b --bbb=B      Set the B value [Default: 2]
  -c              Enable C

There is a lot more we could say about this example, for example,
we could even some real examples of usage.

   ex11.py -a 10
   ex11.py --aaa 10
   ex11.py --aaa=10
"""
from docopt import docopt

def main():
    args = docopt(__doc__, version='example11 1.0')
    print(args)
    print('a=',args['--aaa'])
    print('b=',args['--bbb'])
    print('c=',args['-c'])

if __name__ == '__main__':
    main()

