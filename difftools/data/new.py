#!/bin/python3
""" Original input file for testing diff tools.

This file is for testing the output and rendering of tools for generating and
showing diffs. This is the modified file after changes have been applied against

`old.py`. It's a python file so you can test rendering of diffs with syntax
highlighting.

testing what happens if you delete a single blank line out of many in a
docstring below.

"""
import argparse

CONST1 = 1.0
CONST2 = 2.0

CONST3 = 3

class Foo(object):
  """A test class.

  This class includes most syntactic features of a class.
  """
  CLS_ATTR1 = 'string'
  CLS_ATTR2 = 1.5e-10

  def __init__(self, a, b, c=None):
    """Initialize a Foo instance."""
    if c is None: c = []
    self.a, self.b, self.c = a, b, c

  attra = "another attr for moving to try and trigger zebra mode"
  def baz(self, a, b):
    self.a, self.b = a, b
    return a*b

  @property
  def bar(self):
    return 'bar'

  def addc(self, c):
    # adds c to the end of self.c
    self.c.append(c)

  # This is another multi-line comment for testing moved "zebra" lines, so it
  # seems to need to span multiple lines. I'm not entirely sure why, but maybe
  # single-line moves are not sufficient to trigger the zebra coloring?


# This is a single-line comment that is long for testing moved lines.
if __name__ == '__main__':
  cmdline = argparse.ArgumentParser(description='Test Foo objects.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  cmdline.add_argument('-a', type=int, default=0,
      help='The integer number to use for a.')
  cmdline.add_argument('-b', type=float, default=None,
      help='The float number to use for b.')

  args=cmdline.parse_args()

  foo=Foo(a=args.a, b=args.b)
  print(f'{foo=} {foo.a=} {foo.b=} {foo.c=}')
  print(f'{foo.bar=}')
  print(f'{foo.baz(111,111)=}')
  foo.addc(1)
  foo.addc(2.0)
  print(f'{foo=} {foo.a=} {foo.b=} {foo.c=}')
