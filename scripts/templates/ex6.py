#! /usr/bin/env python
#
import click


@click.command()    
@click.option('--aaa', default='Hello' ,help='The AAA par')
@click.option('--bbb',                  help='The BBB par')

def main(aaa,bbb):
    """
    This is the main command,
    aaa and bbb are parameters
    """
    click.echo(f"aaa={aaa}  bbb={bbb}")

if __name__ == '__main__':
    main()
