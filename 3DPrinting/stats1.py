#!/usr/bin/python3
# Data structures and methods for analysis of data.

from math import inf,sqrt
from collections import defaultdict

class Sample(object):
  """A simple sample stats collector.

  This can be used to collect stats about values by adding them with add() or in batches using update().

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
    self.num += n
    self.sum += v*n
    self.sum2 += v*v*n
    self.min = min(self.min, v)
    self.max = max(self.max, v)

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
    return (self.sum2 - self.sum*self.sum/self.num) / self.num

  @property
  def dev(self):
    """Get the standard deviation of the values."""
    return sqrt(self.var)

  def __repr__(self):
    return f'{self.__class__.__name__}()'

  def __str__(self):
    """Get a simple string rendering of the stats."""
    if self.num:
      return f'num={self.num:g} sum={self.sum:g} min/avg/max/dev={self.min:g}/{self.avg:g}/{self.max:g}/{self.dev:g}'
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
