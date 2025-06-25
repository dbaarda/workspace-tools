#!/bin/pypy3
"""
gcodegen.py gcode generator.

This attempts to define a pressure to flow rate model that is better than the
Linear Advance model at accounting for all the non-linear behaviour people
have encountered. It is not intended to be perfect, just better enough to
justify the additional complexity.

The Linear Advance model is...

https://klipper.discourse.group/t/modification-of-pressure-advance-for-high-speed-bowden-printers/13053?u=dbaarda

The new model is an extension of the linear advance model that adds
back-pressure from the print surface. The nozzle extrudes into a bead of
cooling filament that is pressed against the print surface. This puts back
pressure against the nozzle. The back pressure will be a function of be bead
size; the more pressure, the bigger the bead.

At zero-to-low speeds the pressure does not define the flow rate, just the
bead size. What does define the flow rate is the nozzle travel velocity. As
the nozzle moves it smears the bead, leaving a bead-diameter-wide track of
filament behind, which is the flow rate.

We will assume the back-pressure is a function of only the bead diameter,
which is the track width, and is independent of the nozzle velocity. Note this
could be very wrong, with faster speeds giving the bead less time to cool and
harden, lowering the back-pressure, but we assume this for now. We also assume
that the back-pressure for a given material and layer height is a linear
function with a constant offset of the bead diameter. I could, and have,
invent all sorts of rationalizations for this assumption but I'm basicly
making stuff up so I won't include it here. We'll assume this for simplicity
for now.

The flow through the nozzle also requires pressure. I strongly suspect that
this is Bernoulli’s principle which has pressure proportional to velocity
squared, but the existing PA model assumes Poiseuille’s law with pressure
proportional to velocity. For now we assume Poiseuille’s law, but see the
following discussion for details;

https://klipper.discourse.group/t/modification-of-pressure-advance-for-high-speed-bowden-printers/13053/18?u=dbaarda

This gives us;

  Pe = Pn + Pb      # (1)
  Pn = Kf*vn        # (2)
  Pb = Kb*Db        # (3)
  vn = v*Db*h/Af    # (4)
  Db = r*w          # (5)
  Af = pi*(Df/2)^2  # (6)

Where:

  Pe is the extruder pressure in mm of filament advance.
  Pn is the nozzle flow pressure in mm of filament advance.
  Pb is the bead back-pressure in mm of filament advance.
  Kf is the nozzle pressure factor.
  Kb is the bead backpressure factor.
  vn is the nozzle outflow velocity in mm/s of filament.
  ve is the extruder velocity in mm/s of filament.
  Df is the filament diameter in mm.
  Db is the bead diameter in mm.
  Af is the filament cross-section area in mm^2.
  r is the extrusion ratio.
  w is the line width in mm.
  h is the layer height in mm.
  v is the nozzle velocity in mm/s

For a constant track width at steady state where extruder velocity `ve` equals
nozzle flow velocity `nv` this can be simplified to;

  Pe = Kf*ve + Cf
  Cf = Kb*Db

Where

  Cf is the pressure linear-advance constant.

We can calculate outflow vn as a function of extruder pressure and nozzle
velocity like this;

  Db = vn/v*Af/h                   # (7) from re-arranging (5)
  Pe = Kf*vn + Kb*Db               # from (1) substituting in (2), (3), and (4)
     = Kf*vn + Kb*vn/v*Af/h        # substituting in (7)
     = (Kf + Kb/v*Af/h)*vn
  (Kf + Kb/v*Af/h)*vn = Pe
     vn = Pe /(Kf + Kb/v*Af/h)
     vn = Pe * v/(Kf*v + Kb*Af/h) # (8)

This means we only get nozzle outflow when `Pe > 0` and `v > 0`. As nozzle
velocity v approaches infinity the outflow approaches `Pe/Kf`. The outflow
is half that max rate when `v = Kb/Kf*Af/h`.

Note the bead diameter as a function of extruder pressure Pe and nozzle velocity v is;

  Db = vn/v*Af/h                                # from (4) substituting (7)
     = (Pe * v/(Kf*v + Kb*Af/h))/v * Af/h  # substituting in (8)
     = Pe/(Kf*v + Kb*Af/h) * Af/h
     = Pe/(Kf/Af*h*v + Kb))

This means at nozzle velocity `v=0` the bead diameter is `Pe/Kb`, and as
the nozzle velocity approaches infinity it tends towards zero. The bead
diameter is half the `v=0` max when `v = Kb/Kf*Af/h`.

The rate of change in extruder pressure is;

dPe/dt = ve - vn
       = de/dt - vn
   dPe = de - dt * vn
       = de - dt * Pe * v/(Kf*v + Kb*Af/h)
       = de - Pe * dt/(Kf + Kb*Af/(v*h))

This is an exponentially decaying process with time-constant `T = Kf +
Kb*Af/(v*h)`. This means with no additional extrusion `ve=0`, extruder
pressure `Pe` will decay towards 0 over a time roughly of `Kf + Kb*Af/(h*v)`.
Note that for `v=0` that time is infinity, which means it will not decay at
all.

However, that assumes there is already a bead with backflow pressure Pb to
match Pe. When you first start extruding after a restore, there is no bead and
no backpressure, so there will be outflow till the bead is formed.

If there is no bead yet, there is no backpressure, and there will be outflow
to form a bead even if nozzle velocity v is zero.

Some Adv3 numbers from earlier testing with h=0.3 w=0.6 (db=0.6,eb=0.035mm);

de=+2mm seems to be the minimum for a normal starting dot.
de=+1mm barely seems to register a dot at all.
de=+0.5mm doesn't seem to show any difference from minimum drool.

This suggests roughly Re>0.8 to account for backlash, and Kb=2.0.

de=+2.3mm? for vl=10mm/s, ve=0.75mm/s
de=+3mm for vl=30mm/s, ve=2.25mm/s
de=+4mm for vl=60mm/s, ve=4.5mm/s
de=+5mm seems about right advance for vl=100mm/s ve=7.5mm/s

Suggests roughly Kf=0.4 with Cf=Pb=Kb*0.6+Re=2.0.

Testing suggests you need extra retraction on top of this to avoid drool and
stringing. I suspect this drool is "heat creep expansion", where heat travels
up the filament feed causing expansion. For 0.5x0.2x1.0 at about 60mm/s print
speeds giving `ve=2.5mm/s`, retraction needs to be at least 5mm to avoid this
drool. This suggests settings should be about;

Kf=0.4
Kb=2.0
Re=3.2

However, other testing suggests Kb=4.0 might be more accurate, giving;

Kf=0.4
Kb=4.0
Re=2.4

Adding in a bit of extra retraction headroom would give;

Kf=0.4
Kb=4.0
Re=3.0
"""
import argparse
import re
import sys
from math import e, pi, inf, sqrt
import vtext
from typing import NamedTuple
from functools import cached_property
from pprint import *
from collections import deque

nl = '\n'  # for use inside fstrings.


def attrs(o):
  """ Get all attributes off an object as a dict."""
  return {n:getattr(o,n) for n in dir(o) if not n.startswith('__')}


def ftime(t):
  """Format a time in seconds as hh:mm:ss.sss."""
  h,s = divmod(t,3600)
  m,s = divmod(s,60)
  return f'{int(h):02d}:{int(m):02d}:{s:06.3f}'


def fbool(b):
  """ Format a boolean as 'Y', 'N', or '?' for True, False, None (unknown). """
  return '?' if b is None else 'Y' if b else 'N'


def fstr(s, globals=None, locals=None):
  """ Format a string as an f-string. """
  return eval(f'rf"""{s}"""',globals,locals)


def fesc(s):
  """ Escape a string to prevent f-string expansion problems. """
  return str(s).replace('{', '{{').replace('}','}}')


def isneareq(a, b, n=6):
  """ Is a nearly equal to b?"""
  return abs(a-b) < 10**-n


def isnearle(a, b, n=6):
  """ Is a nearly less than or equal to b? """
  return a <= b + 10**-n


def getnear(a, n=6):
  """ Get the rounded nearest value with n decimal digits."""
  # This will pass through None and also turns -0.0 into 0.0.
  return None if a is None else round(a,n) or 0.0


def acircle(d):
  """ Get the area of a circle from it's diameter. """
  return pi*(d/2)**2


def solveqe(a,b,c):
  """ Return the two roots of a quadratic eqn. """
  d = 0.5/a
  v1 = -b*d
  v2 = sqrt(b**2 - 4*a*c)*d
  return v1 - v2, v1 + v2


class MoveBase(NamedTuple):
  """ Move instructions base type.

  This includes the immutable base attributes of a move instruction from which
  all other attributes can be derived.
  """
  dx: float
  dy: float
  dz: float
  de: float
  v: float  # target velocity in mm/s.
  h: float = None  # The line or end of move height above the current layer.
  w: float = None  # The line width.


class Move(MoveBase):
  """Move instructions.

  These represent gcode moves enhanced to include the start v0, end v1, and
  middle vm velocities. They can partition a basic gcode "intent" move command
  into "physical" acceleration, constant velocity, deceleration components
  based on the printers acceleration, velocity, and cornering specs. They
  assume the printer interprets gcode into movements as described at;

  https://mmone.github.io/klipper/Kinematics.html

  Note for printers with different specs this can be subclassed with the `Ap`,
  `Vp`, and `Gp` overridden by setting them as arguments when subclassing.
  """
  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Ap = 500  # printer acceleration in mm/s^2.
  Vp = 100  # printer max velocity in mm/s.
  Gp = 0.1  # printer cornering deviation distance in mm.
  Fa = acircle(Fd)  # Fillament area in mm^2.
  Na = acircle(Nd)  # Nozzle area in mm^2.

  # The cached properties that depend on velocity attributes.
  __vprops = frozenset('vm dt dt0 dt1 dtm dl0 dl1 dlm'.split())

  def __init_subclass__(cls, Fd=Fd, Nd=Nd, Ap=Ap, Vp=Vp, Gp=Gp):
    super().__init_subclass__()
    cls.Fd, cls.Nd, cls.Ap, cls.Vp, cls.Gp = Fd, Nd, Ap, Vp, Gp
    cls.Fa = acircle(Fd)
    cls.Na = acircle(Nd)
    cls.a0, cls.a1, cls.am = Ap, -Ap, 0.0

  def __new__(cls, *args, v0=None, v1=None, **kwargs):
    """ Create a new Move. """
    # Round arguments to clean out fp errors.
    args = tuple(getnear(a) for a in args)
    kwargs = {k:getnear(v) for k,v in kwargs.items()}
    self = super().__new__(cls, *args, **kwargs)
    assert 0.0 <= self.v <= cls.Vp
    self.v0 = self.v1 = 0.0
    self.setv(v0, v1)
    return self

  def _replace(self, *args, v0=None, v1=None, **kwargs):
    """ Create a new modified Move. """
    # Round arguments to clean out fp errors.
    args = tuple(getnear(a) for a in args)
    kwargs = {k:getnear(v) for k,v in kwargs.items()}
    m = super()._replace(*args,**kwargs)
    m.v0, m.v1 = self.v0, self.v1
    m.setv(v0, v1)
    return m

  def __repr__(self):
    return super().__repr__()[:-1] + f', v0={self.v0}, vm={self.vm}, v1={self.v1})'

  def __str__(self):
    dz, dl, de, h, w, r = self.dz, self.dl, self.de, self.h, self.w, self.r
    v, v0, v1, vm = round(self.v), round(self.v0), round(self.v1), round(self.vm)
    ve, ve0, ve1, vem = round(self.ve, 2), round(self.ve0, 2), round(self.ve1, 2), round(self.vem, 2)
    gostr=f'{dl=:.1f}@{v}({v0}<{vm}>{v1}) {de=:.3f}@{ve}({ve0}<{vem}>{ve1}) l={h:.2f}x{w:.2f}x{r:.2f}'
    if self.isdraw:
      return f'draw {gostr}'
    elif self.ishopup:
      return f'hopup {dz=:.2f}@{v} {h=}'
    elif self.ishopdn:
      return f'hopdn {dz=:.2f}@{v} {h=}'
    elif self.isretract:
      return f'retract {de=:.4f}@{ve}'
    elif self.isrestore:
      return f'restore {de=:.4f}@{ve}'
    else:
      return f'move {gostr}'
    #return f'{self.__class__.__name__}({dx=:.2f}, {dy=:.2f}, {dz=:.2f}, {de=:.4f}, {v=:.2f}, {v0=:.2f}, {v1=:.2f}, {h=:.2f}, {w=:.2f})'

  def setv(self, v0=None, v1=None):
    """ Set the v0, v1 start and end velocities. """
    v0 = self.v0 if v0 is None else getnear(v0)
    v1 = self.v1 if v1 is None else getnear(v1)
    assert self.isgo or v0 == v1 == 0.0
    assert 0.0 <= v0 <= self.v, f'0.0 <= {v0} <= {self.v} for {self!r}'
    assert 0.0 <= v1 <= self.v, f'0.0 <= {v1} <= {self.v} for {self!r}'
    if (v0, v1) != (self.v0, self.v1):
      # Delete cached properties that depend on v0 or v1.
      for a in self.__vprops:
        if a in self.__dict__: delattr(self,a)
      self.v0, self.v1 = v0, v1

  def change(self, *args, s=None, **kwargs):
    """ Change or scale a move. """
    if s is not None:
      # add scaled fields to kwargs if not in args or kwargs.
      fields = self._fields[len(args):4]
      kwargs.update({n:s*getattr(self,n) for n in fields if n not in kwargs})
    return self._replace(*args,**kwargs)

  @classmethod
  def _eb(cls, h, w, r=1.0):
    """ Volume in filament length of bead."""
    return h * acircle(w*r) / cls.Fa

  @property
  def f(self):
    """ The 'f' mm/min speed of a move. """
    return round(self.v * 60)

  @cached_property
  def vm(self):
    """ The speed of the middle phase of a move. """
    if not self.isgo:
      return 0.0
    v, v0, v1 = self.v, self.v0, self.v1
    # This is the limit velocity assuming constant acceleration at a/2.
    # This implements smoothed look-ahead to reduce spikey velocity.
    vm = sqrt((self.dl*self.Ap + v0**2 + v1**2)/2)
    # If vm is less than v0 or v1, we don't need acceleration at both ends,
    # don't set vm higher than the move's speed, and round it.
    vm = getnear(min(v, max(v0, v1, vm)))
    assert v0 <= vm <= v
    assert v1 <= vm <= v
    return vm

  @cached_property
  def r(self):
    """ The extrusion ratio of a move. """
    if self.isdraw:
      return getnear(self.de/self.el)
    return 0.0

  @cached_property
  def dl(self):
    """ The horizontal length of a move. """
    return sqrt(self.dx*self.dx + self.dy*self.dy)

  @cached_property
  def dl0(self):
    """Horizontal length of acceleration phase."""
    dl0 = 0.5*(self.v0+self.vm)*self.dt0
    assert dl0 <= self.dl, f'{dl0=} <= {self.dl=}, {self.dt0=} {self!r}'
    return dl0

  @cached_property
  def dl1(self):
    """Horizontal length of deceleration phase."""
    dl1 = 0.5*(self.vm + self.v1)*self.dt1
    assert dl1 <= self.dl, f'{dl1=} <= {self.dl=}, {self.dt1=} {self!r}'
    return dl1

  @cached_property
  def dlm(self):
    """Horizontal length of middle phase."""
    dlm = self.dl - self.dl0 - self.dl1
    assert 0.0 <= dlm <= self.dl, f'0.0 <= {dlm=} <= {self.dl=}, {self.dl0=} {self.dl1=} {self!r}'
    return dlm

  @cached_property
  def dt0(self):
    """Duration of acceleration phase."""
    return self.dv0 / self.a0

  @cached_property
  def dt1(self):
    """Duration of deceleration phase."""
    return self.dv1 / self.a1

  @cached_property
  def dtm(self):
    """Duration of middle phase."""
    return self.dlm/self.vm if self.vm else abs(self.dz if self.dz else self.de)/self.v

  @cached_property
  def dt(self):
    """ Get the execution time for a move. """
    return self.dt0 + self.dtm + self.dt1

  @property
  def dv0(self):
    """ The acceleration phase change in velocity. """
    return self.vm - self.v0

  @property
  def dv1(self):
    """ The decceleration phase change in velocity. """
    return self.v1 - self.vm

  @property
  def dvm(self):
    """ The mid phase change in velocity. """
    return 0.0

  @property
  def de0(self):
    """ The acceleration phase change in e. """
    return self.de * self.dl0/self.dl if self.dl else 0.0

  @property
  def de1(self):
    """ The decceleration phase change in e. """
    return self.de * self.dl1/self.dl if self.dl else 0.0

  @property
  def dem(self):
    """ The mid phase change in e. """
    return self.de * self.dlm/self.dl if self.dl else self.de

  @property
  def ve(self):
    """ The average extrusion rate of a move. """
    dt = self.dt
    return self.de/dt if dt else 0.0

  @property
  def ve0(self):
    """ The extrusion rate at the start of a move. """
    dl = self.dl
    return self.de/dl * self.v0 if dl else self.ve

  @property
  def vem(self):
    """ The extrusion rate at the middle of a move. """
    dl = self.dl
    return self.de/dl * self.vm if dl else self.ve

  @property
  def ve1(self):
    """ The extrusion rate at the end of a move. """
    dl = self.dl
    return self.de/dl * self.v1 if dl else self.ve

  @property
  def dh0(self):
    """ Change in height of acceleration phase."""
    return self.dz * self.dl0/self.dl if self.dl else 0.0

  @property
  def dh1(self):
    """ Change in height of decceleration phase."""
    return self.dz * self.dl1/self.dl if self.dl else 0.0

  @property
  def dhm(self):
    """ Change in height of mid phase."""
    return self.dz * self.dlm/self.dl if self.dl else self.dz

  @property
  def vl(self):
    """ The average line speed of a move. """
    dt = self.dt
    return self.dl/dt if dt else 0.0

  @property
  def db(self):
    """ The extrusion bead diameter. """
    return self.w*self.r

  @cached_property
  def eb(self):
    """ The extrusion bead volume. """
    return self._eb(self.h, self.w, self.r)

  @cached_property
  def el(self):
    """ The line volume not including r under/over extrusion. """
    return self.h*self.w*self.dl/self.Fa

  @cached_property
  def isgo(self):
    """ Is this a horizontal move? """
    return self.dx or self.dy

  @cached_property
  def isdraw(self):
    # Note we treat moving retractions as pressure adjusted draws.
    return self.isgo and self.de != 0.0

  @property
  def ismove(self):
    return self.isgo and self.de == 0.0

  @property
  def isadjust(self):
    return self.de and not self.isgo

  @property
  def isretract(self):
    return self.de < 0 and not self.isgo

  @property
  def isrestore(self):
    return self.de > 0 and not self.isgo

  @property
  def ishop(self):
    return self.dz and not self.isgo

  @property
  def ishopup(self):
    return self.dz > 0 and not self.isgo

  @property
  def ishopdn(self):
    return self.dz < 0 and not self.isgo

  @property
  def isspeed(self):
    """ Is this just a speed change? """
    return not (self.isgo or self.dz or self.de)

  @property
  def isaccel(self):
    """ Is this a partitioned acceleration or deceleration phase move? """
    dt = self.dt
    return self.dt0/dt > 0.9 or self.dt1/dt > 0.9

  def costh(self, m):
    """ Get cos(a) of the angle a between two moves.

    Note this is not the change in direction angle, but the angle between the
    two moves. So 180deg means no change in direction, and 0deg is a full
    reverse.
    """
    # cos(th) = v0.v1/(|v0|*|v1|)
    dx0, dy0, dl0 = self.dx, self.dy, self.dl
    dx1, dy1, dl1 = m.dx, m.dy, m.dl
    # if a move has no horizontal movement, treat as 0deg.
    return -(dx0*dx1 + dy0*dy1)/(dl0 * dl1) if dl0 and dl1 else 1.0

  def joind(self, m):
    """ Calculate the max path deviation when joining two moves.

    This calculates the max deviation distance between two moves and the
    direct path if you join them into a single move.

    The largest deviation of the line C from the pair of lines A then B, is
    the smallest distance between the point at the end of line A and the line
    C. If we make the start point of A the origin, then A, B, and C are
    vectors of length a,b, and c.

    Any point on the line C can be represented as t*C, which is just scaling
    the vector C. Points for t<0 are before the start of C, and for t>1 are
    after the end of line C. The shortest distance between point A and line C
    is when t is;

    t = A.C/C.C = A.C/c^2

    For 0<t<1 the shortest distance is between the points t.C and A, which is
    the length of the normal vector from C to A;

    d = |A.unitnormal(C)| which is the same as d = |CxA|/|C|
    """
    # The x and y values for triangle sides a, b, and c.
    ax,ay = self.dx, self.dy
    bx,by = m.dx, m.dy
    cx,cy = ax+bx, ay+by
    c2 = cx*cx + cy*cy      # C.C AKA c^2
    t = (ax*cx + ay*cy)/c2  # t = A.C/C.C
    if 0<t:
      # point A is "before" the start of C, shortest dist is from start of C to A, or just a.
      return self.dl
    elif 1<t:
      # point A is "after" the end of C, shortest dist is from end of C to A, or just b.
      return m.dl
    else:
      # point A is "between" ends of C, shortest dist is normal from C to A.
      return abs(ax*cy - ay*cx)/sqrt(c2)

  def cornerv(self, m):
    """ Calculate the max corner velocity between two moves.

    This calculates the max velocity given acceleration `a` for a circular
    path between two moves that deviates `g` distance from the corner. We
    ignore dz as insignificant and use horizontal velocities only. See;

    https://onehossshay.wordpress.com/2011/09/24/improving_grbl_cornering_algorithm/
    """
    # Limit sinth2 to just below 1 to avoid divide by zero.
    # sin(th/2) = sqrt((1-cos(th))/2)
    sinth2 = min(sqrt((1 - self.costh(m))/2), 0.9999)
    # R = g * sin(th/2)/(1-sin(th/2))
    R = self.Gp*sinth2/(1-sinth2)
    # vmax = sqrt(a*R)
    vmax = sqrt(self.Ap*R)
    return min(vmax, self.v, m.v)

  def join(self, m):
    """ Add two moves together. """
    return self.change(dx=self.dx+m.dx, dy=self.dy+m.dy, dz=self.dz+m.dz,
        de=self.de+m.de, v=max(self.v, m.v), v0=self.v0, v1=m.v1)

  def split(self, dlmin=1.0):
    """ Partition a move into a list of moves for the acceleration phases. """
    v0,v1,vm,dl = self.v0, self.v1, self.vm, self.dl
    if v0 == vm == v1 or not self.isdraw:
      # A non-moving, constant speed, or non-extruding move, don't partition it.
      return [self]
    # Get distances for move phases, with a little bit of extra room for acceleration.
    dl0, dlm, dl1 = self.dl0 + 1e-6, self.dlm - 2e-6, self.dl1 + 1e-6
    # initialize middle phase v0 and v1 to just vm.
    mv0 = mv1 = vm
    # if phase distances are small don't bother partitioning them.
    if dl0 < dlmin or dlm < dlmin:
      dl0, dlm, mv0 = 0.0, dl0+dlm, v0
    if dlm < dlmin or dl1 < dlmin:
      dl1, dlm, mv1 = 0.0, dlm+dl1, v1
    # Make a list of the s,v0,v1 values for each phase.
    phases = [(dl0/dl, v0, vm), (dlm/dl, mv0, mv1), (dl1/dl, vm, v1)]
    # Return a list moves for the non-zero phases.
    return [self.change(v0=v0, v1=v1, s=s) for s,v0,v1 in phases if s]

  def canjoin(self, m, lh=None):
    """ Can this move be reasonably joined with m?

    This needs the layer height to figure out if these are hopped moves which
    can always be joined. We can't arbitrarily join moves at layer height
    because they could be avoiding outlines.
    """
    if lh is None: lh = self.Nd
    return (
        # Join moves together if ...
        (self.ismove and m.ismove and
        # ... they are hopped or straight enough.
        ((self.h > lh and m.h > lh) or self.joind(m) < self.Gp)) or
        # Join draws together if ...
        (self.isdraw and m.isdraw and
        # ... they have the same line settings, and ...
        (self.h,self.w,self.r) == (self.h,self.w,self.r) and
        # The path deviation from joining them is small.
        self.joind(m) < self.Gp) or
        # Always join hopup/hopdn together.
        (self.ishop and m.ishop) or
        # Always join retract/restore together.
        (self.isadjust and m.isadjust))


  def fixv(self, v0=0.0, v1=0.0, m0=None, m1=None):
    """ Update move velocities for moves before/after.

    This updates a move's velocities to what the printer will do based on on
    the moves before and after. The v0 or v1 values will be set less than the
    provided velocities or moves if necessary to fit within the Ap
    acceleration and deceleration limits.

    if m0 is set the v0 start velocity is set to the corner velocity.
    If m1 is set the v1 end velocity is set to the corner velocity.
    If v0, or v1 are still not set they default to 0.0.
    """
    if m0 is not None: v0 = m0.cornerv(self)
    if m1 is not None: v1 = self.cornerv(m1)
    dl = self.dl
    if not dl:
      # A stop move, v0, v1, and vm are zero.
      self.v0 = self.v1 = self.vm = 0.0
    else:
      # limit v0 and v1 for acceleration with a tiny bit of headroom.
      if v0 < v1:
        v1 = min(v1, sqrt(v0**2 + 2*dl*self.Ap) - 1e-6)
      elif v0 > v1:
        v0 = min(v0, sqrt(v1**2 + 2*dl*self.Ap) - 1e-6)
      self.setv(v0, v1)


def ldt(l):
  """ Get the execution time for a list of moves. """
  return sum(m.dt for m in l if isinstance(m,Move))


class Layer(NamedTuple):
  """ Layer definition object. """
  # Note we always start layers from a high "hop" position.
  n : int   # layer number (0 for pre-extrude).
  z : float # height of bottom of layer.
  Vp: int   # default layer print speed.
  Vt: int   # etc.
  Vz: int
  Ve: int
  Vb: int
  h : float # layer height.
  w : float # default line width.
  r : float # default line extrusion ratio.


class GCodeGen(object):
  """ GCodeGen gcode generator.

  The general operation of this is;

  1. An instance is created with the desired arguments.
  2. Doing resetfile() initializes the file state between "runs". The nextcmds
     argument can be used to set self.nextcmds with a list of pending cmds.
     Setting finalize=True will mean any cmds added will be finalized and
     rendered before they are added.
  3. Commands are added with add(code), which executes them to update state
     and records them in self.prevcmds. Commands can be Layers, Moves, tuples
     with gcode instructions, strings of gcode lines, binary blobs or other
     arbitrary things added/supported by subclasses.
  5. After the commands are all added, getGcode() can be called to get the
     encoded binary gcode blob. This applys several post-processing
     transformations before finally transforming them into byte encoded gcode.
     The transformations applied depends on initialization arguments.

  The string commands can include Python f-strings that are evaluated in the
  final transformation to gcode in a namespace containing all the execution
  state at that point. This means they can include complex expressions
  evaluated from the current execution state for the gcode they generate.

  There are many helper methods that generate and add() commands to do various
  different things. These include;

  * startFile()/endFile(): add gcode initialization/finalization cmds.
  * startLayer()/endLayer(): add layer initialization/finalization cmds.
  * nextLayer(): shorthand for ending the previous layer and starting a new one.
  * cmd(): add an arbitrary gcode command.
  * cmt(): add a gcode comment line.
  * log(): add a gcode comment log line, depending on self.verb setting.
  * draw(), move(), retract(), restore(), hopup(), hopdn(): add move cmds
    for basic actions.
  * preExt(), fbox(), dbox(), dot(), line(), text(): add multiple move
    commands for complex output.

  The add(cmd) method is the underlying method to execute and record a cmd in
  self.prevcmds. It filters out empty cmds that evaluate to False. It always
  just adds encoded bytes. For other commands it normally executes them with
  inc(cmd) and appends them to prevcmds, but if resetfile() was called with
  finalise=True, then it adds them by calling fadd(cmd).

  The inc(cmd) method is used to execute a command and increment the state. It
  dispatches to different inc<cmd>(code) methods for executing the different
  types of commands. The default behaviour is to not change the state.

  The fadd(cmd) method will finalize a command and render it into encoded
  binary that it appends to prevcmds with add(cmd). For Moves it dispatches to
  the addmove() method. For anything else it executes fmt(cmd) to format the
  command, inc(cmd) to increment the state, and finally recursively calls
  add(fcmd). Note fmt() is called before the state is incremented, and the
  result is recursively formatted and encoded after the state is updated. This
  means fmt() can return strings with content evaluated based on the state
  before the cmd is executed, and include f-strings that evaluate based on the
  state after it's executed. After recursion fmt() will eventually render it
  into encoded binary that add() appends to prevcmds.

  The faddmove(m) method will apply final adjustments to the move based on
  cmdline arguments, render it with fmove(), execute it with inc(), and then
  add() the rendered format string. It will also add() additional comments
  and logging if necessary.

  The fmt(cmd) method takes a command and renderes it into a text format
  string, or for strings into formatted and encoded binary. It dispatches to
  f<cmd>(cmd) methods for handling the different types of cmds.

  Post-processing is done by mod<feature>() post-processor methods which
  process and update or replace self.prevcmds. These are selectively executed
  by a top-level mod() post-processor depending on initialization arguments.
  mod() is called by getCode() before appending all prevcmds into the final
  encoded result.

  Post-processors can directly modify or replace self.prevcmds, or they can do
  `resetfile(nextcmds=self.prevcmds)`, and then pop the commands from the
  front of self.nextcmds, modify them, and replay them with add() back into
  self.prevcmds. This allows them to use the state generated by the modified
  cmds, lookbacks into self.prevcmds, and lookaheads into self.nextcmds, to
  drive the modifications of each cmd.

  The final transformation into encoded gcode is implemented as a modfmt()
  post-processor. It does resetfile(nextcmds=self.prevcmds, initialize=True)
  and replays the commands with add(), which executes and transforms them into
  their encoded binary with fadd(). This means faddmove() can use
  self.nextcmds for lookaheads when finalizing the move commands.

  Attributes:
    Te: Extruder temp (degC).
    Tp: Platform temp (degC).
    Fe: Extruder fan speed (0.0->1.0).
    Fc: Case fan speed (0.0->1.0).
    fKf: Firmware Linear Advance factor (0.0->4.0 mm/mm/s).
    Kf: Linear Advance factor (0.0->4.0 mm/mm/s).
    Kb: Bead backpressure factor (??? mm/mm).
    Re: Retraction distance (0.0->10.0 mm).
    Vp: Speed when printing (mm/s).
    Vt: Speed when traveling (mm/s).
    Vz: Speed when raising/lowering (mm/s).
    Ve: Speed when retracting (mm/s).
    Vb: Speed when restoring (mm/s).
    Lr: Default extrusion ratio (0.1 -> 10.0).
    Lh: Default line height (0.1->0.4mm).
    Lw: Default line width (0.2->0.8mm).
    relext: relative extrusion enabled.
    dynret: linear advance dynamic retraction enabled.
    dynext: linear advance dynamic extrusion enabled.
    optmov: path optimizing enabled.
    f_t: current file execution time.
    f_l: current file printed line length.
    f_e: current file extruded filament length.
    l_t: current layer execution time.
    l_l: current layer printed line length.
    l_e: current layer extruded filament length.
    pe: current extruder pressure or retraction.
    eb: current extruded bead volume.
    db: current extruded bead diameter (property).
    vn: nozzle flow velocity (property).
    ve: current extruder velocity.
    vl: current line velocity.
    dt: the execution time of the last command.
    dl: the line length of the last command.
    de: the extruder change of the last command.
    en: the extruded volume of the last command.
    c: the last command.
    m: the last move command.
    d: the last move command that was a draw.
    x,y,z,e,f: Current x,y,z,e,f position values.
    h: current height above the layer base.
    w: current default line width (property).
    r: current default extrusion ratio (property).
  """

  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Fa = acircle(Fd)  # Fillament area.
  Na = acircle(Nd)  # Nozzle area.
  Zh = Nd/2  # z-hop height.

  startcode = """\
;gcode_flavor: flashforge
;machine_type: Adventurer 3 Series
;filament_diameter0: {Fd}
;right_extruder_temperature: {Te}
;platform_temperature: {Tp}
;extruder_fan: {round(Fe*100)}%
;case_fan: {round(Fc*100)}%
;base_print_speed: {round(Vp)}
;travel_speed: {round(Vt)}
;vertical_speed: {round(Vz)}
;retract_speed: {round(Ve)}
;restore_speed: {round(Vb)}
;layer_height: {Lh:.2f}
;line_width: {Lw:.2f}
;extrusion_ratio: {Lr:.2f}
;command_line: {_getCmdline()}
;start gcode
M118 X{75.0:.2f} Y{75.0:.2f} Z{150.0:.2f} T0
M140 S{round(Tp)} T0
M104 S{round(Te)} T0
M104 S0 T1
M107
M900 K{fKf:.3f} T0
G90
G28
M132 X Y Z A B
G1 Z50.000 F420
G161 X Y F3300
M7 T0
M6 T0
{_setFc(Fc)}
M108 T0
{_setFe(Fe)+nl if Fe else ""}\
{"M83"+nl if relext else ""}\
"""
  endcode = """\
;end gcode
M107
M104 S0 T0
M140 S0 T0
G162 Z F1800
G28 X Y
M132 X Y A B
M652
G91
M18
"""

  def _getCmdline(self):
    return ' '.join(sys.argv)

  def _setFc(self, fc):
    return f'M652 S{round(fc*255)}' if fc else 'M651'

  def _setFe(self, fe):
    return f'M106' if fe == 1.0 else f'M106 S{round(fe*255)}' if fe else 'M107'

  def __init__(self,
      Te=210, Tp=50, Fe=1.0, Fc=1.0, fKf=0.0,
      Kf=0.0, Kb=0.0, Re=5,
      Vp=60, Vt=100, Vz=7, Ve=40, Vb=30,
      Lh=0.2, Lw=Nd, Lr=1.0,
      relext=False, optmov=False, segvel=False, fixvel=0,
      diff=False, verb=False,
      dynret=False, dynext=False):
    # Temp and Fan settings.
    self.Te, self.Tp, self.Fe, self.Fc, self.fKf = Te, Tp, Fe, Fc, fKf
    # Linear advance and retraction settings.
    self.Kf, self.Kb, self.Re = Kf, Kb, Re
    # Default velocities.
    self.Vp, self.Vt, self.Vz, self.Ve, self.Vb = Vp, Vt, Vz, Ve, Vb
    # Default line height, width, and extrusion ratio.
    self.Lh, self.Lw, self.Lr = Lh, Lw, Lr
    # Processing mode options.
    self.relext, self.optmov, self.segvel, self.fixvel = relext, optmov, segvel, fixvel
    # Logging mode options.
    self.diff, self.verb = diff, verb
    # Pressure advance mode options.
    self.dynret, self.dynext = dynret, dynext
    class GMove(Move, Fd=self.Fd, Nd=self.Nd): pass
    self.GMove = GMove
    self.resetfile()

  @property
  def db(self):
    """ Get the current bead diameter from the bead volume."""
    eb, h = self.eb, self.h
    # h<=0 can happen for initializing moves or moves at the start of a new
    # layer before a hop.
    return 2*sqrt(eb*self.Fa/(h*pi)) if h > 0.0 else 0.0

  @property
  def vn(self):
    """Get nozzle flow velocity from the nozzle pressure and back pressure."""
    pe, db = self.pe, self.db
    pb = max(0.0, self.Kb*db)            # bead backpressure.
    pn = max(0.0, pe - pb)               # nozzle pressure.
    return pn/self.Kf if self.Kf else self.ve if pe>=pb else 0.0

  @property
  def w(self):
    """ Get the current default width."""
    # This is the width of the previous draw or the layer default.
    return self.d.w if self.d else self.layer.w

  @property
  def r(self):
    """ Get the current default extrusion ratio."""
    # This is the ratio of the previous draw or the layer default.
    return self.d.r if self.d else self.layer.r

  def islayer(self, cmd):
    return isinstance(cmd, Layer)

  def ismove(self, cmd):
    return isinstance(cmd, Move)

  def iscmd(self, cmd):
    return type(cmd) == tuple

  def isstr(self, cmd):
    return isinstance(cmd, str)

  def isbin(self, cmd):
    return isinstance(cmd, bytes)

  def resetfile(self, nextcmds=None, finalize=False):
    """ Reset all file state variables.

    This gets ready to process or finalize a queue of nextcmds or new commands.
    """
    # Zero all positions so the first move's deltas are actually absolute.
    self.x, self.y, self.z, self.e, self.f = 0.0, 0.0, 0.0, 0.0, 0
    self.f_t = self.f_l = self.f_e  = 0.0
    self.pe = self.eb = self.vl = self.ve = self.h = 0.0
    self.prevcmds, self.nextcmds, self.finalize = deque(), nextcmds, finalize
    self.resetlayer()

  def resetlayer(self):
    """ Reset all layer state variables. """
    # Set layer to a default base layer for any layer-less moves.
    self.layer = Layer(
        -1, 0.0, self.Vp, self.Vt, self.Vz, self.Ve, self.Vb, self.Lh, self.Lw, self.Lr)
    self.l_t = self.l_l = self.l_e = 0.0
    self.c = self.m = self.d = None
    self.resetstate()

  def resetstate(self):
    """ Reset all move state variables. """
    self.dt = self.de = self.dl = self.en = 0.0

  def iterstate(self, dt, de, dl, dv, dh):
    """ Increment the state for a small dt interval. """
    # Note filament doesn't get sucked back into the nozzle if the pressure
    # advance is negative, and the bead cannot give negative backpressure.
    # Get filament extruded into the nozzle and height.
    pe, eb, db, h = self.pe, self.eb, self.db, self.h
    pe += de
    pb = max(0, self.Kb * db)                 # bead backpressure.
    pn = max(0, pe - pb)                      # nozzle pressure.
    en = pn * dt/(self.Kf + dt)               # filament extruded out the nozzle.
    el = min(db * h * dl / self.Fa, eb + en)  # fillament removed from bead to line,
    self.pe = pe - en
    self.eb = max(0.0, eb + en - el)
    self.dt += dt
    self.de += de
    self.en += en
    self.dl += dl
    self.vl += dv
    self.h += dh

  def incstate(self, dt, de=0.0, dl=0.0, dv=0.0, dh=0.0):
    """ Increment the state for a large dt interval. """
    # Keep the expected final state to check against and use later.
    dt1, de1, dl1, vl1, h1 = self.dt+dt, self.de+de, self.dl+dl, self.vl+dv, self.h+dh
    t, idt = dt, 0.001
    if dl:
      # integrate a horizontal move.
      al, de_dl, dh_dl = dv/dt, de/dl, dh/dl
      while t > 0.0:
        idt = min(idt, t)
        idv = al * idt
        idl = (self.vl + idv/2) * idt  # Do trapezoidal integration for dl.
        ide = idl * de_dl
        idh = idl * dh_dl
        self.iterstate(idt, ide, idl, idv, idh)
        t -= idt
    else:
      # integrate a non-horizontal change.
      ve, vz = de/dt, dh/dt
      while t > 0.0:
        idt = min(idt, t)
        ide = ve * idt
        idh = vz * idt
        self.iterstate(idt, ide, 0.0, 0.0, idh)
        t -= idt
    assert isneareq(self.dt, dt1)
    assert isneareq(self.de, de1)
    assert isneareq(self.dl, dl1)
    assert isneareq(self.vl, vl1)
    assert isneareq(self.h, h1)
    # Use the expected final states to remove small integration errors.
    self.dt, self.de, self.dl, self.vl, self.h = dt1, de1, dl1, vl1, h1

  def inc(self, code):
    if self.islayer(code):
      self.inclayer(code)
    elif self.ismove(code):
      self.incmove(code)
    elif self.iscmd(code):
      self.inccmd(code)
    self.c = code

  def inct(self, dt):
    """ Increment time. """
    self.f_t += dt
    self.l_t += dt
    self.dt = dt

  def inclayer(self, l):
    """ Increment state for starting a Layer. """
    self.resetlayer()
    self.layer = l
    self.h = getnear(self.z - l.z)

  def incmove(self, m):
    """ Increment state for executing a Move. """
    assert isneareq(self.h, self.z - self.layer.z), f'{self.h=} != {self.z - self.layer.z}'
    assert isneareq(self.vl, m.v0), f'{self.vl=} != {m.v0=} for {m} at {ftime(self.f_t)}'
    self.resetstate()
    # Save the expected final h to assert against later.
    h1 = self.h + m.dz
    # Update the state for the three move phases.
    if m.dt0: self.incstate(m.dt0, m.de0, m.dl0, m.dv0, m.dh0)
    if m.dtm: self.incstate(m.dtm, m.dem, m.dlm, m.dvm, m.dhm)
    if m.dt1: self.incstate(m.dt1, m.de1, m.dl1, m.dv1, m.dh1)
    # Check that the end states are as expected.
    assert isneareq(self.dt, m.dt)
    assert isneareq(self.de, m.de)
    assert isneareq(self.dl, m.dl)
    assert isneareq(self.vl, m.v1)
    assert isneareq(self.h, h1)
    # Update and round positions to remove floating point errors.
    self.x = getnear(self.x+m.dx)
    self.y = getnear(self.y+m.dy)
    self.z = getnear(self.z+m.dz)
    self.e = getnear(self.e+m.de)
    self.f = m.f
    self.f_l += m.dl if self.en else 0.0
    self.f_e += self.en
    self.l_l += m.dl if self.en else 0.0
    self.l_e += self.en
    self.vl = m.v1
    self.ve = m.ve1
    self.h = getnear(self.z - self.layer.z)
    self.m = m
    if m.isdraw: self.d = m
    self.inct(m.dt)

  def inccmd(self, cmd):
    """ Increment state for special gcode commands. """
    self.resetstate()
    cmd, kwargs = cmd
    if cmd == 'G92':
      # This is a "set XYZE to ..." cmd, often used to keep E < 1000.
      if 'X' in kwargs: self.x = kwargs['X']
      if 'Y' in kwargs: self.y = kwargs['Y']
      if 'Z' in kwargs: self.z = kwargs['Z']
      if 'E' in kwargs: self.e = kwargs['E']
    elif cmd == 'G4':
      # This is a pause command.
      dt = kwargs['P']/1000
      self.incstate(dt=dt)
      self.inct(dt)

  def fmt(self, code):
    if self.islayer(code):
      return self.flayer(code)
    elif self.ismove(code):
      return self.fmove(code)
    elif self.iscmd(code):
      cmd, kwargs = code
      return self.fcmd(cmd, **kwargs)
    elif isinstance(code, str):
      # Format and encode strings.
      return self.fstr(code.strip()).encode()
    else:
      # It's probably Flashprint binary header stuff.
      return code

  def flayer(self, l):
    """ Format a layer as a comment string. """
    # layer 0 is the preExtrude
    return self.fcmt('{"layer" if layer.n else "preExtrude"}:{layer.h:.2f}')

  def fmove(self, m):
    """ Format a Move as a gcode command. """
    dx,dy,dz,de = m[:4]
    # Don't include args in cmd if they are unchanged.
    # Include both x and y if either are changed.
    x = self.x + dx if dx or dy else None
    y = self.y + dy if dy or dx else None
    z = self.z + dz if dz else None
    # If relext or diff mode is on, use relative extrusion.
    e = (de if self.relext or self.diff else self.e + de) if de else None
    f = m.f if m.f != self.f else None
    cmd = self.fcmd('G1', X=x, Y=y, Z=z, E=e, F=f)
    if self.verb:
      # Add additional command comment suffix.
      cmd = f'{cmd:40s}; {{m}} r={{en/m.el if m.el else db/m.w:.2f}}'
    return cmd

  def fcmd(self, cmd, **kwargs):
    """ Format a gcode command. """
    fmts = dict(X='.2f',Y='.2f',Z='.3f',E='.4f',F='d')
    args = ' '.join(f'{k}{v:{fmts.get(k,"")}}' for k,v in kwargs.items() if v is not None)
    return f'{cmd} {args}' if args else cmd

  def fstr(self, str):
    """ Format a string as an f-string with attributes of self."""
    try:
      return fstr(str, locals=attrs(self))
    except Exception as e:
      raise Exception(f'Error formatting {str!r}.') from e

  def fcmt(self, cmt):
    """ Format a comment. """
    return f';{cmt}'

  def flog(self, txt):
    """ Conditionally format a log comment. """
    #print(self.fstr(txt))
    if self.verb:
      if self.diff:
        return self.fcmt('{dt:.3f}: ' + txt)
      else:
        return self.fcmt('{ftime(f_t)}: ' + txt)

  def faddmove(self, m):
    """ Do a final format and add of a Move. """
    if m.isretract and self.dynret:
      # Adjust retraction to Re and also relieve advance pressure.
      #self.log(f'adjusting {m} for {pe=:.4f}')
      m = m.change(de=-self.Re - self.pe)
      # If de is zero, skip adding this.
      if not m.de:
        self.log('skipping empty retract.')
        return
    elif m.isrestore and self.dynret:
      # Get the pe and eb for the next draws.
      pe, eb = self._getNextPeEb()
      # Adjust existing restore and add the starting bead.
      # self.log(f'adjusting {c} for {pe=:.4f} {eb=:.4f}')
      # Lets also try P control with Kp=0.8 for restores.
      m = m.change(de=0.8*pe - self.pe + eb)
      # If de is zero, skip adding this.
      if not m.de:
        self.log('skipping empty restore.')
        return
    elif m.isdraw and self.dynext:
      # For calculating the pressure pe needed, use the ending ve1.
      pe = self._calc_pe(m.ve1, m.db)
      # Adjust de to include required change in pe over the move.
      # self.log(f'adjusting {c} {c.ve=:.4f}({c.ve0:.4f}<{c.vem:.4f}>{c.ve1:.4f}) {c.isaccel=}')
      # TODO: This tends to oscillate, change to a PID or PD controller.
      # Lets try the P control first with Kp=0.8 to try and reduce the overshoot.
      m = m.change(de=m.de + 0.8*(pe - self.pe))
    if self.e + m.de > 1000.0:
      # Reset the E offset before it goes over 1000.0
      self.cmd('G92', E=0)
    fmove = self.fmove(m)
    self.inc(m)
    self.add(fmove)
    self.log('{pe=:.4f} {eb=:.4f} {ve=:.4f} {vn=:.4f} {vl=:.3f} {db=:.2f}')

  def fadd(self, code):
    """ Do a final format and add. """
    if self.ismove(code):
      # For moves use faddmove().
      self.faddmove(code)
    else:
      # Format it, increment state, and recursively add() the formatted code.
      fcode = self.fmt(code)
      self.inc(code)
      self.add(fcode)

  def add(self, code):
    """ Add code Layer, Move, Wait, or format-string instances. """
    if not code:
      # Skip empty commands.
      return
    elif isinstance(code, bytes):
      # Just append encoded binary.
      self.prevcmds.append(code)
    elif self.finalize:
      # In finalize mode add it with fadd().
      self.fadd(code)
    else:
      # Increment state and append it.
      self.inc(code)
      self.prevcmds.append(code)

  def cmd(self, cmd, **kwargs):
    """ Add a gcode command as a tuple. """
    # Strip out None arguments.
    kwargs = {k:v for k,v in kwargs.items() if v is not None}
    self.add((cmd, kwargs))

  def cmt(self, cmt):
    """ Add a comment. """
    self.add(self.fcmt(cmt))

  def log(self, txt):
    """ Conditionally add a log comment. """
    #print(self.fstr(txt))
    self.add(self.flog(txt))

  def startFile(self):
    self.add(self.fstr(self.startcode))
    # Note Flashprint always generates and assumes pre-extrude has h=0.2mm.
    self.startLayer(0, h=0.2)

  def endFile(self):
    # Retract to -2.5mm for the next print.
    self.hopup(pe=-2.5)
    self.endLayer()
    self.filestats()
    self.add(self.fstr(self.endcode))

  def startLayer(self, n=None, z=None,
      Vp=None, Vt=None, Vz=None, Ve=None, Vb=None,
      h=None, w=None, r=1.0):
    if n is None: n = self.layer.n + 1
    if z is None: z = self.layer.z + self.layer.h if n > 1 else 0.0
    l = Layer(
        n, z,
        Vp or self.Vp,
        Vt or self.Vt,
        Vz or self.Vz,
        Ve or self.Ve,
        Vb or self.Vb,
        h or self.Lh,
        w or self.Lw,
        r*self.Lr)
    self.add(l)
    self.log('layer={layer.n} h={layer.h:.2f} w={layer.w:.2f} r={layer.r:.2f}')
    if not isneareq(l.r, self.Lr, 2):
      self.cmt(f'extrude_ratio:{round(layer.r, 2):.3g}')

  def endLayer(self):
    self.layerstats()

  def nextLayer(self, **largs):
    self.endLayer()
    self.startLayer(**largs)

  def move(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      v=None, s=1.0, h=None, w=None, r=0.0):
    """ Move the extruder.

    Movement can be absolute using `x,y,z,e` arguments, or relative using
    `dx,dy,dz,de` arguments. The z axis can also be specified with `h` height
    above the base of the current layer. The extrusion can also be specified
    using an `r` extrusion ratio muliplier for a `h*w` line, where `h` and `w`
    default to the current layer's if not specified, and r is multiplied by
    the current layer's extrusion ratio.  Speed can be specified directly
    using the `v` velocity argument, or the `s` speed multiplier of the
    default `Vp,Vt,Vz,Ve,Vb` speeds for print, travel, vertical, retract, and
    restore movement.

    Note this is a simple raw movement. Any pressure advance compensation
    will be added in post-processing.
    """
    if z is None and dz is None and h is not None: z = self.layer.z + h
    dx = x - self.x if x is not None else dx if dx is not None else 0.0
    dy = y - self.y if y is not None else dy if dy is not None else 0.0
    dz = z - self.z if z is not None else dz if dz is not None else 0.0
    dl = sqrt(dx**2 + dy**2)
    # Note: setting h and z or dz means explicit h overrides z height for de calcs.
    if h is None: h = self.z + dz - self.layer.z
    if w is None: w = self.w
    r = self.layer.r * r
    el = dl * w * h / self.Fa  # line filament volume not including extrusion ratio.
    de = e - self.e if e is not None else de if de is not None else r * el
    if v is None:
      if de:
        v = s * (self.layer.Vp if dl else self.layer.Ve if de < 0 else self.layer.Vb)
      else:
        v = s * (self.layer.Vt if dl else self.layer.Vz)
    else:
      # Always keep explicit f values as they can be used for future moves.
      self.f = round(v*60)
    m = self.GMove(dx,dy,dz,de,v,h=h,w=w)
    # Only add it if it does anything.
    if any(m[:4]):
      #if m.isdraw and not isneareq(m.r, self.r, 2):
      #  self.cmt(f'extrude_ratio:{round(m.r, 2):.3g} {m.r=} {self.r=}')
      self.add(m)
    else:
      self.log(f'skipping empty move {m}')

  def draw(self, x=None, y=None, z=None, r=None, **kwargs):
    """ Draw a line.

    This is the same as move() but with the default previous draw or layer r
    value so it draws instead of travels by default.
    """
    # We need r before scaling by layer.r
    if r is None: r = self.r / self.layer.r
    self.move(x, y, z, r=r, **kwargs)

  def retract(self, e=None, de=None, ve=None, s=1.0, pe=None):
    """Do a retract.

    This is a retraction that by default relieves any pressure and retracts to
    pe which defaults to -self.Re. Setting de or e will do a fixed retraction
    ignoring the current pressure. Note compensation for advance pressure is
    adjusted in post processing if enabled.
    """
    if pe is None: pe = -self.Re
    if de is None: de = -self.pe + pe
    #if de is None:
    #  de = -self.Re if pe is None else -self.pe + pe
    # If dynamic retraction is enabled and de is zero or more, set de to a
    # tiny retraction so that it doesn't get optimized away or seen as a
    # restore and can be dynamically adjusted later.
    if self.dynret and de >= 0: de = -0.00001
    self.move(e=e, de=de, v=ve, s=s)

  def restore(self, e=None, de=None, vb=None, s=1.0, pe=None):
    """ Do a restore.

    This is a restore that by default reverts any retraction and restores
    pressure to pe which defaults to 0.0. Setting e or de will do a fixed
    restore ignoring current retraction. Note adding starting dots and
    compensating for advance pressure is added in post-processing if enabled.
    """
    if pe is None: pe = 0.0
    if de is None: de = -self.pe + pe
    #if de is None:
    #  de = self.Re if pe is None else -self.pe + pe
    # If dynamic retraction is enabled and de is zero or less, set de to a
    # tiny restore so that it doesn't get optimized away or seen as a
    # retraction and can be dynamicly adjusted later.
    if self.dynret and de <= 0: de = 0.00001
    self.move(e=e, de=de, v=vb, s=s)

  def hopup(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      vt=None, vz=None, ve=None, vb=None,
      s=1.0, h=None, pe=None):
    """ Do a retract and raise. """
    # default hopup is to Zh above layer height.
    if h is None: h = self.layer.h + self.Zh
    self.retract(e=e, de=de, ve=ve, s=s, pe=pe)
    self.move(z=z, dz=dz, v=vz, s=s, h=h)

  def hopdn(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      vt=None, vz=None, ve=None, vb=None,
      s=1.0, h=None, pe=None):
    """ Do a move, drop, and restore. """
    # default hop down is to layer height.
    if h is None: h = self.layer.h
    if (x, y, dx, dy) != (None, None, None, None):
      self.move(x=x, y=y, dx=dx, dy=dy, v=vt, s=s)
    self.move(z=z, dz=dz, v=vz, h=h, s=s)
    self.restore(e=e, de=de, vb=vb, s=s, pe=pe)

  def wait(self, t):
    """ Do a pause. """
    self.cmd('G4', P=round(t*1000))

  def _calc_pe(self, ve, db):
    """ Get steady-state pressure advance pe needed for target ve and db. """
    pn = max(0, self.Kf*ve)
    pb = max(0, self.Kb*db)
    return pn+pb

  def _getNextPeEb(self):
    "Get pe and eb for the next draws."
    dt = dl = vedl = dbdl = eb = 0.0
    # Get the line-average ve and db for the next 10mm and 1s of moves.
    for m1 in (m for m in self.nextcmds if self.ismove(m)):
      if not m1.isdraw or (dl > 10 and dt > 1.0): break
      dt+=m1.dt
      dl+=m1.dl
      vedl+=m1.ve*m1.dl
      dbdl+=m1.db*m1.dl
      # Get eb from the first move.
      if not eb:
        eb = m1.eb
        # For dynext, only use the first move.
        if self.dynext: break
    if dl:
      # Calculate pe using the line-average ve and db.
      return self._calc_pe(vedl/dl, dbdl/dl), eb
    return 0.0, 0.0

  def modjoin(self):
    """Postprocess gcode to join moves."""
    # Note this resets the file but doesn't re-process the cmds, it just appends them.
    nextcmds = self.prevcmds
    self.resetfile()
    l = self.layer
    i = m = None
    for c in nextcmds:
      # Keep the layer for checking moves against.
      if self.islayer(c):
        l = c
      elif self.ismove(c):
        # Skip joining moves before layer 1 to not optimize pre-extrudes.
        if l.n > 0 and m and m.canjoin(c, l.h):
          # replace this entry with the new joined m.
          self.prevcmds[i] = m = m.join(c)
          # Replace the joined move with a log entry.
          self.log(f'Joined {c} into previous Move.')
          continue
        else:
          # Try to join this with later moves.
          i, m = len(self.prevcmds), c
      # Never join across other gcode commands.
      elif self.iscmd(c):
        m = None
      self.prevcmds.append(c)

  def modfixv(self):
    """ Fix v0,v1 velocities of moves in gcode. """
    # This adjusts moves in self.prevcmds in-place.

    def _rfixv(prevs, m, v0, m1):
      """ Fix m and previous moves for acceleration limits. """
      # Note if m1 is None this will use v1=0.0 instead.
      m.fixv(v0=v0, v1=0.0, m1=m1)
      # fix previous moves if their v1 doesn't match the next v0.
      for m0 in reversed(prevs):
        # if v1 is already set to the target, stop.
        if m0.v1 == m.v0: break
        m0.fixv(v0=m0.v0, v1=m.v0)
        m = m0

    prevs = []
    v0, m = 0.0, None
    # Only go through non-comment commands using None for non-moves.
    cmds = (m for m in self.prevcmds if not isinstance(m, str))
    for m1 in (m if self.ismove(m) else None for m in cmds):
      if m:
        # fixv m and possibly previous moves for next move m1.
        # This will use v1=0.0 to stop if m1 is not a move.
        _rfixv(prevs, m, v0=v0, m1=m1)
        # Append m to the list of previous moves for _rfix().
        prevs.append(m)
        # v0 is the last m.v1
        v0 = m.v1
      # m1 is the next m.
      m = m1
    if m:
      # Fix the last move for final v1.
      _rfixv(prevs, m, v0=v0, m1=None)

  def modsegvel(self):
    """ Partition moves into accel/cruise/decel phases. """
    splits = (c.split() if self.ismove(c) else [c] for c in self.prevcmds)
    self.prevcmds = deque(p for s in splits for p in s)

  def modfixvel(self):
    """ Fix v velocities for acceleration limits according to fixvel level. """
    # This replaces self.prevcmds with the adjusted commands.

    def _fixvel(m):
      """Fix v depending on fixvel level of velocity fix required."""
      if self.ismove(m) and m.isgo and (m.vm < m.v) and (
          # Any slow move.
          (self.fixvel == 3) or
          # Accelerating move.
          (self.fixvel >= 2 and m.v0<m.vm==m.v1) or
          # Decelerating move.
          (self.fixvel >= 1 and m.v0==m.vm>m.v1)):
        # Set v to vm rounded down to an int.
        v = int(m.vm)
        return m.change(v=v, v0=min(m.v0, v), v1=min(m.v1, v))
      return m

    # set velocity to vm depending on self.fixvel.
    self.prevcmds = deque(_fixvel(c) for c in self.prevcmds)
    # Run modfixv() again to update v0,v1.
    self.modfixv()

  def modfmt(self):
    """Finalize and format the commands."""
    self.resetfile(nextcmds=self.prevcmds, finalize=True)
    while self.nextcmds:
      self.add(self.nextcmds.popleft())

  def mod(self):
    if self.optmov:
      # Join together any moves that can be joined.
      self.modjoin()
    # Fix all the velocities.
    self.modfixv()
    if self.fixvel:
      # set velocity to vm for acceleration limits.
      self.modfixvel()
    if self.segvel:
      # Partition all the moves into phases
      self.modsegvel()
    # finalize and return the resulting gcode.
    self.modfmt()

  def getCode(self):
    self.mod()
    return b'\n'.join(self.prevcmds) + b'\n'

  def layerstats(self):
    h='{layer.h:.2f}'
    w='{layer.w:.2f}'
    r='{(Fa * l_e) / (layer.h * layer.w * l_l) if l_l else layer.r:.2f}'
    v='{l_l/l_t if l_t else 0.0:.1f}'
    ve='{l_e/l_t if l_t else 0.0:.3f}'
    self.log('layer={layer.n} finished.')
    self.log(f'  avg:{h}x{w}x{r}@{v} ve={ve}')
    self.log('  t={ftime(l_t)} l={l_l:.1f} e={l_e:.1f}')

  def filestats(self):
    h='{(layer.z+layer.h)/layer.n:.2f}'
    w='{Lw:.2f}'
    r='{(Fa * f_e) / ((layer.z+layer.h)/(layer.n or 1) * Lw * f_l) if f_l else Lr:.2f}'
    v='{f_l/f_t if f_t else 0.0:.1f}'
    ve='{f_e/f_t if f_t else 0.0:.3f}'
    lt='{ftime(f_t/(layer.n or 1))}'
    self.log('file finished with {layer.n} layers.')
    self.log(f'  avg:{h}x{w}x{r}@{v} ve={ve} lt={lt}')
    self.log('  t={ftime(f_t)} l={f_l:.1f} e={f_e:.1f}')

  def preExt(self,x0,y0,x1,y1,le=120.0, lw=60.0, m=10.0, ve=20.0, vw=20, h=0.2, r=4):
    """Preextrude around a box with margin m.

    This extrudes for `le` mm at speed `ve` mm/s, then wipes for another `lw`
    mm at speed `vw` mm/s without extruding, and finally if we did a wipe
    another 10mm at 1mm/s to stick the last bit of wipe drool to the plate.
    This should leave the nozzle primed with no residual advance pressure, and
    prevent the wipe from becoming a big piece of stringing that messes with
    the print.

    Set `lw=0` to turn off the wipe and stick phases. Note this includes a
    hopdn before and hopup at the end.
    """
    self.cmt('structure:pre-extrude')
    x0,x1 = sorted([x0,x1])
    y0,y1 = sorted([y0,y1])
    x0-=m
    x1+=m
    y0-=m
    y1+=m
    x, y = x0, y0
    l, v = le, ve
    self.hopdn(x=x,y=y,h=h)
    while l:
      if self.x == x0 and self.y<y1:
        dy = min(y1 - self.y, l)
        self.draw(dy=dy,v=v,r=r)
        l-=dy
      elif self.y == y1 and self.x < x1:
        dx = min(x1 - self.x, l)
        self.draw(dx=dx,v=v,r=r)
        l-=dx
      elif self.x == x1 and self.y > y0:
        dy = min(self.y - y0, l)
        self.draw(dy=-dy, v=v, r=r)
        l-=dy
      else:
        dx = min(self.x - x0, l)
        self.draw(dx=-dx, v=v, r=r)
        l-=dx
      if (self.x, self.y) == (x1,y1):
        # Reached top-right, widen y0,y1
        y0-=2
        y1+=2
      elif (self.x, self.y) == (x1,y0):
        # reached bottom-right, widen x0,x1
        x0-=2
        x1+=2
      if l == 0 and r > 0:
        # We have finished extruding, go another lw distance without extruding
        # to relieve pressure and wipe ooze.
        l, v, r = lw, vw, 0
      elif l == 0 and r == 0 and v>1 and lw:
        # We have finished the wipe, go another 10mm at 1mm/s to stick the
        # final drool to the plate.
        l, v, r = 10,1,0
    self.hopup()

  def fbox(self,x0,y0,x1,y1,shells=None,**kwargs):
    """ Fill a box with frame lines. """
    if shells is None: shells = inf
    w = kwargs.get('w', self.w)
    # make sure x0<x1 and y0<y1
    x0,x1 = sorted((x0,x1))
    y0,y1 = sorted((y0,y1))
    x0+=w/2
    x1-=w/2
    y0+=w/2
    y1-=w/2
    s=0
    self.hopdn(x0,y0)
    while x1-x0>2*w and y1-y0>2*w and s<shells:
      if shells==0:
        self.cmt('structure:shell-outer')
      elif s==1:
        self.cmt('structure:shell-inner')
      self.draw(x0,y1,**kwargs)
      self.draw(x1,y1,**kwargs)
      self.draw(x1,y0,**kwargs)
      self.draw(x0+w,y0,**kwargs)
      x0+=w
      x1-=w
      y0+=w
      y1-=w
      s+=1
    if y1-y0 > 2*w and s<shells:
      fw=x1-x0
      self.move(x=(x0+x1)/2,y=y0+(fw-w)/2)
      self.draw(y=y1-(fw-w)/2, w=fw)
    elif s<shells:
      fw=y1-y0
      self.move(x=x0+(fw-w)/2,y=(y1+y0)/2)
      self.draw(x=x1-(fw-w)/2, w=fw)
    self.hopup()

  def dbox(self,x0,y0,x1,y1,shells=3,**kwargs):
    """ fill a box with diagonal lines. """
    w = kwargs.get('w', self.w)
    x0,x1 = sorted((x0,x1))
    y0,y1 = sorted((y0,y1))
    if shells:
      self.fbox(x0,y0,x1,y1,shells,**kwargs)
    b=(shells+0.5)*w
    x0,y0,x1,y1 = getnear(x0+b), getnear(y0+b), getnear(x1-b), getnear(y1-b)
    self.cmt('structure:infill-solid')
    dx=w*2**0.5
    # Diagonal direction depends on odd vs even layers.
    if self.layer.n % 2:
      # initialize left and right diagonal ends at top left corner.
      lx,ly=rx,ry=x0,y1
      dy = -dx
    else:
      # initialize left and right diagonal ends at bottom left corner.
      lx,ly=rx,ry=x0,y0
      dy = dx
    self.hopdn(rx,ry)
    while lx<x1:
      if rx<x1:
        # right diagonal end is following top or bottom edge right.
        rx+=dx
        if rx>x1:
          # hit the right corner, start going along right side.
          if ry == y0:
            ry+=(rx-x1)
          else:
            ry-=(rx-x1)
          rx=x1
      else:
        # right diagonal end is following right side.
        ry+=dy
        ry=min(y1,max(y0,ry))
      if lx==x0:
        # left diagonal end is following left side.
        ly+=dy
        if ly<y0:
          # hit the bottom corner, start following bottom.
          lx+=y0-ly
          ly=y0
        elif y1<ly:
          # hit the top corner, start following top.
          lx+=ly-y1
          ly=y1
      else:
        # left diagonal end is following the top.
        lx+=dx
        if lx>x1: lx=x1
      #self.log(f'pos={{(x,y)}} l={(lx,ly)} r={(rx,ry)}')
      if self.y == ry or self.x == rx or (rx - self.x) < (self.x - lx):
        # starting diagonal from right point.
        self.draw(x=rx,y=ry,**kwargs)
        self.draw(x=lx,y=ly,**kwargs)
      else:
        # starting diagonal from left point.
        self.draw(x=lx,y=ly,**kwargs)
        self.draw(x=rx,y=ry,**kwargs)
    self.hopup()

  def dot(self, x, y, r=1.0, **kwargs):
    # Note this doesn't seem to render anything in FlashPrint.
    self.hopdn(x,y)
    # Draw a tiny line so FlashPrint renders something.
    self.draw(dy=-0.2,r=r, **kwargs)
    self.hopup()

  def line(self, l, r=1.0, **kwargs):
    """Draw a line from a sequence of points."""
    x0,y0 = l[0]
    # If there is only one point, draw a dot instead.
    if len(l) == 1 or (len(l) == 2 and l[0] == l[1]):
      return self.dot(x0, y0, r, **kwargs)
    self.hopdn(x0,y0)
    for x,y in l[1:]:
      self.draw(x,y,**kwargs)
    self.hopup()

  def text(self, t, x0, y0, x1=None, y1=None, fsize=5, angle=0, **kwargs):
    self.cmt('structure:shell-outer')
    self.log(f'text {[(x0,y0),(x1,y1)]} {fsize=} {angle=} {t!r}')
    w = kwargs.get('w', self.w)
    v = vtext.ptext(t,x0,y0,x1,y1,fsize,angle,w)
    for l in v:
      self.line(l,**kwargs)


def RangeType(min=-inf, max=inf):
  def init(s):
    v=type(min)(s)
    if v<min or max<v:
      raise argparse.ArgumentTypeError(f'must be between {min} and {max}')
    return v
  return init


def GCodeGenArgs(cmdline):
  """Add argparse cmdline arguments for GCodeGen params."""
  cmdline.add_argument('-Te', type=RangeType(0, 265), default=245,
      help='Extruder temperature.')
  cmdline.add_argument('-Tp', type=RangeType(0, 100), default=100,
      help='Platform temperature.')
  cmdline.add_argument('-Fe', type=RangeType(0.0,1.0), default=1.0,
      help='Extruder fan speed between 0.0 to 1.0.')
  cmdline.add_argument('-Fc', type=RangeType(0.0,1.0), default=0.0,
      help='Case fan speed between 0.0 to 1.0.')
  cmdline.add_argument('-fKf', type=RangeType(0.0,4.0), default=0.0,
      help='Firmware Linear Advance factor between 0.0 to 4.0 in mm/mm/s.')
  cmdline.add_argument('-Kf', type=RangeType(0.0,4.0), default=0.0,
      help='Linear Advance factor between 0.0 to 4.0 in mm/mm/s.')
  cmdline.add_argument('-Kb', type=RangeType(0.0,10.0), default=0.0,
      help='Bead backpressure factor between 0.0 to 10.0 in mm/mm.')
  cmdline.add_argument('-Re', type=RangeType(0.0,10.0), default=1.0,
      help='Retraction distance between 0.0 to 10.0 in mm.')
  cmdline.add_argument('-Vp', type=RangeType(5,100), default=60,
      help='Base printing speed in mm/s.')
  cmdline.add_argument('-Vt', type=RangeType(5,100), default=100,
      help='Base travel speed in mm/s.')
  cmdline.add_argument('-Vz', type=RangeType(1,10), default=7,
      help='Base raise/lower speed in mm/s.')
  cmdline.add_argument('-Ve', type=RangeType(1,50), default=40,
      help='Base retract speed in mm/s.')
  cmdline.add_argument('-Vb', type=RangeType(1,50), default=30,
      help='Base restore speed in mm/s.')
  cmdline.add_argument('-Lh', type=RangeType(0.1,0.4), default=0.3,
      help='Layer height in mm.')
  cmdline.add_argument('-Lw', type=RangeType(0.2,0.8), default=0.5,
      help='Line width in mm.')
  cmdline.add_argument('-Lr', type=RangeType(0.0,10.0), default=1.0,
      help='Line extrusion ratio between 0.0 to 10.0.')
  cmdline.add_argument('-E', action='store_true',
      help='Enable relative extrusion.')
  cmdline.add_argument('-O', action='store_true',
      help='Enable optimize, merging draw/move cmds where possible.')
  cmdline.add_argument('-S', action='store_true',
      help='Enable segmenting moves into accel/cruise/decel phases.')
  cmdline.add_argument('-V', type=RangeType(0,3), default=0,
      help='Enable fix velocities for acceleration limits (0=off,1=decel,2=decel+accel,3=all).')
  cmdline.add_argument('-d', action='store_true',
      help='Enable diff mode output.')
  cmdline.add_argument('-v', action='store_true',
      help='Enable verbose logging output.')
  cmdline.add_argument('-R', action='store_true',
      help='Enable Linear Advance dynamic retract/restore.')
  cmdline.add_argument('-P', action='store_true',
      help='Enable Linear Advance dynamic extrusion rates.')


def GCodeGetArgs(args):
  """ Get the dict of standard GCode arguments. """
  return dict (
      Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc, fKf=args.fKf,
      Kf=args.Kf, Kb=args.Kb, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve, Vb=args.Vb,
      Lh=args.Lh, Lw=args.Lw, Lr=args.Lr,
      relext=args.E, optmov=args.O, segvel=args.S, fixvel=args.V,
      diff=args.d, verb=args.v,
      dynret=args.R, dynext=args.P)


if __name__ == '__main__':
  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=GCodeGen(**GCodeGetArgs(args))

  gen.startFile()
  gen.preExt(-10,-10,10,10)
  gen.nextLayer()
  gen.dbox(-10,-10,10,10)
  gen.nextLayer()
  gen.dbox(-10,-10,10,10)
  gen.nextLayer()
  gen.text("Hello\nWorld!",-8,8,8,-8)
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
