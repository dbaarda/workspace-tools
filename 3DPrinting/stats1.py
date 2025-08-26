#!/usr/bin/python3
# Data structures and methods for analysis of data.
from collections import defaultdict

try:
  # Use numpy if available to support vector versions.
  from numpy import inf, sqrt, log, exp, fmin, fmax, array, printoptions
except ImportError:
  # If it's not available use our own minimal array vector.
  from collections.abc import Iterable
  from math import inf, sqrt as _sqrt, exp as _exp, log as _log

  class array(list):
    """A lightweight pypy3 compatible numpy.array alternative."""

    def __format__(self, format_spec=''):
      return f'[{", ".join(f"{v:{format_spec}}" for v in self)}]'

    def __mul__(self, other):
      if isinstance(other, Iterable):
        return self.__class__(x*y for x,y in zip(self,other,strict=True))
      return self.__class__(x*other for x in self)

    def __truediv__(self, other):
      if isinstance(other, Iterable):
        return self.__class__(x/y for x,y in zip(self,other,strict=True))
      return self.__class__(x/other for x in self)

    def __add__(self, other):
      if isinstance(other, Iterable):
        return self.__class__(x+y for x,y in zip(self,other,strict=True))
      return self.__class__(x+other for x in self)

    def __sub__(self, other):
      if isinstance(other, Iterable):
        return self.__class__(x-y for x,y in zip(self,other,strict=True))
      return self.__class__(x-other for x in self)

    def __rtruediv__(self,other):
      if isinstance(other, Iterable):
        return self.__class__(x/y for x,y in zip(other,self,strict=True))
      return self.__class__(other/x for x in self)

    def __rsub__(self, other):
      if isinstance(other, Iterable):
        return self.__class__(x-y for x,y in zip(other,self,strict=True))
      return self.__class__(other-x for x in self)

    __rmul__ = __mul__
    __radd__ = __add__
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __itruediv__ = __truediv__

  ## define numpy-like vectorizing versions of math functions.

  def sqrt(x):
    if isinstance(x,array):
      return x.__class__(_sqrt(v) for v in x)
    return _sqrt(x)

  def log(x):
    """log to the base e function."""
    if isinstance(x,array):
      return x.__class__(_log(v) for v in x)
    # Like math.log() except returns -inf for 0.0 instead of ValueError.
    return _log(x) if x else -inf

  def exp(x):
    if isinstance(x,array):
      return x.__class__(_exp(v) for v in x)
    return _exp(x)

  def fmin(a,b):
    if isinstance(a,array):
      if isinstance(b,array):
        return a.__class__(min(x,y) for x,y in zip(a,b,strict=True))
      return a.__class__(min(x,b) for x in a)
    if isinstance(b,array):
      return b.__class__(min(a,x) for x in b)
    return min(a,b)

  def fmax(a,b):
    if isinstance(a,array):
      if isinstance(b,array):
        return a.__class__(max(x,y) for x,y in zip(a,b,strict=True))
      return a.__class__(max(x,b) for x in a)
    if isinstance(b,array):
      return b.__class__(max(a,x) for x in b)
    return max(a,b)

  def clip(v,vmin,vmax):
    return fmax(vmin,fmin(vmax,v))


class P(object):
  """A wrapper for float 'p' format-spec support for numbers and arrays.

  The 'p' format specifier formats like 'f' but also strips off any trailing
  zeros or decimal point. The '#' alternate form will also strip any leading
  zero. It can wrap numbers, tuples, lists, or arrays. The default format is
  'g'. This can be used like this;

  f'{P(1.1001):.3p}' -> '1.1'
  f'{P(0.1234):#.3p}' -> '.123'
  f'{P([1,2,3,4]):.1f}' -> '[1.0, 2.0, 3.0, 4.0]'
  """
  def __init__(self, v):
    self.v = v

  def __repr__(self):
    return f'{self.__class__.__name__}({repr(self.v)})'
  
  def __str__(self):
    return format(self)
  
  def __format__(self, format_spec=''):
    self.spec=format_spec if format_spec else 'g'
    return self._f(self.v)

  def _fp(self, v):
    """format a float using self.spec as the format spec."""
    # TODO: support 'P' like 'F'.
    if self.spec.endswith('p'):
      s = f'{v:{self.spec[:-1] + "f"}}'
      if '.' in s: s = s.rstrip('0').rstrip('.')
      if '#' in self.spec and any(s.startswith(p) for p in ('0.', '-0.', '+0.', ' 0.')):
        s = s.replace('0', '', 1)
      return s
    return f'{v:{self.spec}}'

  def _f(self, v):
    """format anything using self.spec as the format spec."""
    if isinstance(v, (float,int)):
      return self._fp(v)
    if isinstance(v, tuple):
      return f"({' '.join(self._f(e) for e in v)})"
    if isinstance(v, list):
      return f"[{' '.join(self._f(e) for e in v)}]"
    if hasattr(v, 'ndim'):
      # assume it's a numpy.array, use printoptions.
      with printoptions(edgeitems=2, threshold=5, formatter=dict(
          float_kind=self._fp, int_kind=self._fp)):
        return format(v)
    # fallback to str().
    return str(v)


class Sample(object):
  """A simple sample stats collector.

  This can be used to collect stats about values by adding them with add() or
  in batches using update(). The values can be multiple value arrays (lists
  and tuples will be converted to arrays), and all the corresponding
  attributes will also be arrays. Note arrays will be numpy.array if numpy is
  available, otherwise they will be a numpy-array-like list subclass.

  Attributes:
    num: the number of values added.
    sum: the sum of values added.
    sum2: the sum of the square of the values.
    min: the minimum value added.
    max: the maximum value added.
    avg: the average of the values added.
    var: the variance of the values added.
    dev: the stddev of the values added.
  """

  def __init__(self, data=None):
    self.num = 0
    self.sum = 0
    self.sum2 = 0
    self.min = inf
    self.max = -inf
    if data:
      self.update(data)

  def add(self, v, n=1):
    """Add a value."""
    # Convert tuples and lists to np.array.
    if isinstance(v,(tuple,list)): v=array(v)
    self.num += n
    self.sum += n*v
    self.sum2 += n*v*v
    self.min = fmin(self.min, v)
    self.max = fmax(self.max, v)

  def update(self, data):
    """Add all the values in data."""
    for v in data:
      self.add(v)

  @property
  def avg(self):
    """Get the average of the values."""
    return self.sum / self.num

  @property
  def var(self):
    """Get the variance of the values."""
    # clamp >=0 to correct for rounding errors when var=0.
    return fmax(0.0,(self.sum2 - self.sum*self.sum/self.num) / self.num)

  @property
  def dev(self):
    """Get the standard deviation of the values."""
    return sqrt(self.var)

  def atdev(self, dev):
    """Get the value at +devs."""
    return self.avg + dev*self.dev

  def lognorm(self):
    """Get the mu, sigma parameters for a lognormal distribution."""
    r = self.dev/self.avg
    s2 = log(r*r + 1.0)
    return log(self.avg) - 0.5*s2, sqrt(s2)

  def atlogdev(self, dev):
    """Get the value at +devs for a lognormal distribution."""
    m, s = self.lognorm()
    return exp(m + dev*s)

  def __repr__(self):
    return f'{self.__class__.__name__}()'

  def __str__(self):
    """Get a simple string rendering of the stats."""
    if self.num:
      return f'num={self.num:g} sum={P(self.sum)} min/avg/max/dev={P(self.min)}/{P(self.avg)}/{P(self.max)}/{P(self.dev)}'
    else:
      return "num=0 sum=0"

  def getstate(self):
    """Get a snapshot of the Sample state."""
    return (self.num, self.sum, self.sum2, self.min, self.max)

  def setstate(self, state):
    """Set the state of the Sample from a previos snapshot."""
    self.num, self.sum, self.sum2, self.min, self.max = state


class Histogram(Sample):
  """A Sample that also builds a histogram.

  This has all the general stats attributes collected by Sample, but also has
  a histogram breakdown of the values. The histogram has buckets that can be
  fixed or exponentially sized, with the bucket lower bounds defined by the
  formula;

    bmin[i] = width*i + scale*base^(i-1)   # for i>0
    bmin[i] = 0                            # for i=0
    bmin[i] = width*i + scale*base^(-i-1)  # for i<0

  By default width, scale, and base are all 0. Setting width>0 gives fixed
  size buckets. Setting scale > 0, and base>1 gives exponential buckets that
  are mirrored for negative values.

  The buckets can also be limited between hmin (default -inf) and hmax
  (default inf) values. If these are set the histogram is limited in range,
  with the lowest bucket index being the one that would include hmin, and the
  largest bucket index the one that would include hmax.
  """

  def __init__(self, data=None, width=0, scale=0, base=0, hmin=-inf, hmax=inf):
    """Initialize a Histogram instance."""
    if not width and not scale: width=1
    self.width, self.scale, self.base, self.hmin, self.hmax = width, scale, base, hmin, hmax
    self.data = defaultdict(int)
    super().__init__(data=data)

  def __str__(self):
    """Get a simple string rendering of the stats and histogram."""
    scale = int(self.num/len(self.data)) if self.data else 1
    return super().__str__() + '\n' + ''.join(f'[{bmin:8.1f},{bmax:8.1f}): {c:8d} {"*" * (10*c//scale)}\n' for bmin,bmax,c in self)

  def _getindex(self, v):
    """Get the index i of the bucket that v belongs in."""
    if v < self.hmin:
      return self._getindex(self.hmin)
    elif self.hmax < v:
      return self._getindex(self.hmax)
    elif v < 0:
      return next(i for i in range(-1,-10000,-1) if self._getbound(i) <= v)
    else:
      return next(i for i in range(1,10000) if self._getbound(i) > v) - 1

  def _getbound(self, i):
    """Get the lower bound of bucket i."""
    if i == 0:
      v = 0.0
    elif i > 0:
      v = self.width*i + self.scale*self.base**(i-1)
    else:
      v = self.width*i - self.scale*self.base**(-i-1)
    return -inf if v < self.hmin else inf if v > self.hmax else v

  def __iter__(self):
    """Iterate through the (bmin, bmax, count) histogram buckets."""
    imin,imax = (min(self.data),max(self.data)) if self.data else (0,0)
    bmax = self._getbound(imin)
    for i in range(imin,imax+1):
      bmin, bmax, c = bmax, self._getbound(i+1), self.data[i]
      yield bmin, bmax, c

  def add(self, v, n=1):
    """Add a value."""
    super().add(v, n)
    self.data[self._getindex(v)] += n


if __name__ == '__main__':
  s = Sample()
  s.add(1.0)
  s.add(2.0)
  print(s)
