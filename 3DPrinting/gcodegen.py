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
       = de - Pe * dl/(Kf*v + Kb*Af/h)
       = de - Pe * dl*dt/(Kf*dl + Kb*Af/h*dt)
       = de - Pe * dt/(Kf + Kb*Af/(v*h))

This is an exponentially decaying process over time and line-distance. Looking
at it over time it has a decay time-constant `T = Kf + Kb*Af/(v*h)`. Looking
at it over distance it has a decay distance-constant of `L = Kf*v + Kb*Af/h`.
This means with no additional extrusion `ve=0`, extruder pressure `Pe` will
decay towards 0 over a time roughly of `Kf + Kb*Af/(h*v)`. Similarly it will
decay towards zero over a distance of `Kf*v + Kb*Af/h`. Note that for `v=0`
the time to decay is infinity, which means it will not decay at all.

However, that assumes there is already a bead with backflow pressure Pb to
match Pe. When you first start extruding after a restore, there is no bead and
no backpressure, so there will be outflow till the bead is formed.

If there is no bead yet, there is no backpressure, and there will be outflow
to form a bead even if nozzle velocity v is zero.

Some Adv3 numbers from earlier testing with h=0.3 w=0.6 (db=0.6,eb=0.035mm);

de=+2mm seems to be the minimum for a normal starting dot.
de=+1mm barely seems to register a dot at all.
de=+0.5mm doesn't seem to show any difference from minimum ooze.

This suggests roughly Re>0.8 to account for backlash, and Kb=2.0.

de=+2.3mm? for vl=10mm/s, ve=0.75mm/s
de=+3mm for vl=30mm/s, ve=2.25mm/s
de=+4mm for vl=60mm/s, ve=4.5mm/s
de=+5mm seems about right advance for vl=100mm/s ve=7.5mm/s

Suggests roughly Kf=0.4 with Cf=Pb=Kb*0.6+Re=2.0.

Testing suggests you need extra retraction on top of this to avoid ooze and
stringing. I suspect this ooze is "heat creep expansion", where heat travels
up the filament feed causing expansion. For 0.5x0.2x1.0 at about 60mm/s print
speeds giving `ve=2.5mm/s`, retraction needs to be at least 5mm to avoid this
ooze. This suggests settings should be about;

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
import sys
from math import e, pi, inf, sqrt, copysign
from typing import NamedTuple
from functools import cached_property
from pprint import *
from collections import deque
import vtext
import stats1

nl = '\n'  # for use inside fstrings.

# These are short formatting functions for use in fstrings.
from stats1 import P


# Short fstring functions for formatting numbers and numpy.array values.
def f0(v, o=''):
  return f'{P(v):{o}.0f}'

def f1(v, o=''):
  return f'{P(v):{o}.1p}'

def f2(v, o=''):
  return f'{P(v):{o}.2p}'

def f3(v, o=''):
  return f'{P(v):{o}.3p}'

def f6(v, o=''):
  return f'{P(v):{o}.6p}'

def fzp(v, p):
  """ Format a float with max p decimal places and strip the leading zero."""
  return f'{P(v):#.{p}p}'


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
  return eval(f'rf"""{s}"""', globals, locals)


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


def clamp(v, vmin, vmax):
  """Clamp a value between vmin and vmax limits."""
  return vmin if v < vmin else vmax if v > vmax else v


def acircle(d):
  """ Get the area of a circle from it's diameter. """
  return pi/4*(d)**2


def dcircle(a):
  """ Get the diameter of a circle from it's area. """
  return sqrt(a*4/pi)


def solveqe(a,b,c):
  """ Return the two roots of a quadratic eqn. """
  # (-b +- sqrt(b^2 - 4*a*c))/(2*a)
  d = 0.5/a   # 1/(2*a)
  v1 = -b*d   # -b/(2*a)
  v2 = sqrt(b**2 - 4*a*c)*d  # sqrt(b^2 - 4*a*c)/(2*a)
  return v1 - v2, v1 + v2


class AttrDictMixin(object):
  """A mixin to add dict-like access to attributes for fstr()."""

  def __getitem__(self, key):
    try:
      return getattr(self, key)
    except AttributeError:
      raise KeyError(key)


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

  These represent gcode moves enhanced to include the target extrusion ratio
  `r` and the start/middle/end velocities `v0`/`vm`/`v1`. The default `r` is
  calculated from `de` but can be explicitly set to indicate the intent when
  `de` includes pressure adjustments. The velocities partition a basic gcode
  "intent" move command into "physical" acceleration, constant velocity,
  deceleration components based on the printers acceleration, velocity, and
  cornering specs. They assume the printer interprets gcode into movements as
  described at;

  https://mmone.github.io/klipper/Kinematics.html

  Note for printers with different specs this can be subclassed with the `Ap`,
  `Vp`, `Ae`, `Ve`, `Gp`, and `Je` overridden by setting them as arguments
  when subclassing. The values used here are taken from the OrcaSlicer
  defaults for the Flashforge Adventurer 3 with changes `Vp = 150-> 100 mm/s`,
  `Ve = 30 -> 40 mm/s` to match published and/or apparent specs.

  Restore's are assumed to draw a circular dot with a target diameter w*r and
  Retractions, Hops and Moves don't draw anything. Draws draw a line of target
  width w*r and length dl with a half-dot "byte" out of the start from the
  previous restore or draw, and corresponding "bulge" on the end. The default
  r for draws and moves is calculated from the ratio of extrusion to
  line-volume, otherwise it is the retract/restore distance.

  The `de` extrusion for any of these can be adjusted to negative or positive
  values for pressure adjustments and will not necessarily correspond with the
  target `r` extrusion ratios.
  """
  Fd = 1.75 # Fillament diameter.
  Nd = 0.4  # Nozzle diameter.
  Ap = 500  # printer acceleration in mm/s^2.
  Vp = 100  # printer max velocity in mm/s.
  Ae = 500  # extruder acceleration in mm/s^2.
  Ve = 40   # extruder max velocity in mm/s.
  Gp = 0.1  # printer cornering deviation distance in mm.
  Je = 2.5  # extruder "jerk" speed deviation in mm/s.
  Fa = acircle(Fd)  # Fillament area in mm^2.
  Na = acircle(Nd)  # Nozzle area in mm^2.
  La = 0.5*Vp**2/Ap # Acceleration to Vp distance in mm.
  # These are constant middle phase acceleration and velocity change attributes.
  am = aem = dvm = dvem = 0.0

  # The cached properties that depend on velocity attributes.
  __vprops = frozenset('vm dt dt0 dt1 dtm dl0 dl1 dlm ve0 vem ve1 de0 dem de1'.split())

  def __init_subclass__(cls, Fd=Fd, Nd=Nd, Ap=Ap, Vp=Vp, Gp=Gp, Ae=Ae, Ve=Ve, Je=Je):
    super().__init_subclass__()
    cls.Fd, cls.Nd = Fd, Nd
    cls.Ap, cls.Vp, cls.Gp = Ap, Vp, Gp
    cls.Ae, cls.Ve, cls.Je = Ae, Ve, Je
    cls.Fa = acircle(Fd)
    cls.Na = acircle(Nd)
    cls.La = 0.5*cls.Vp**2/cls.Ap

  def __new__(cls, *args, r=None, v0=None, v1=None, **kwargs):
    """ Create a new Move. """
    # Round arguments to clean out fp errors.
    args = tuple(getnear(a) for a in args)
    kwargs = {k:getnear(v) for k,v in kwargs.items()}
    self = super().__new__(cls, *args, **kwargs)
    if r is None:
      r = 0.0 if self.de <=0 else self.de/self.el if self.isgo else 1.0
    self.r = getnear(r)
    self.v0 = self.v1 = 0.0
    self.setv(v0, v1)
    assert 0.0 <= self.r, f'{self.r=} for {self!r}'
    assert 0.0 <= self.v <= cls.Vp, f'{self.v=} {cls.Vp=} for {self!r}'
    return self

  def _replace(self, *args, r=None, v0=None, v1=None, **kwargs):
    """ Create a new modified Move. """
    # Round arguments to clean out fp errors.
    args = tuple(getnear(a) for a in args)
    kwargs = {k:getnear(v) for k,v in kwargs.items()}
    m = super()._replace(*args,**kwargs)
    m.r = self.r if r is None else getnear(r)
    m.v0, m.v1 = self.v0, self.v1
    m.setv(v0, v1)
    return m

  def __repr__(self):
    return super().__repr__()[:-1] + f', r={self.r}, v0={self.v0}, v1={self.v1})'

  def __str__(self):
    lstr = f'({f2(self.h)}x{f2(self.w)}x{f2(self.r)})'
    pstr = f' dl={self.dl:.1f}@{f0(self.v)}({f0(self.v0)}<{f0(self.vm)}>{f0(self.v1)})' if self.dl else ''
    hstr = f' dz={self.dz:.2f}@{f0(self.v)}' if self.dz else ''
    estr = f' de={self.de:.4f}@{f2(self.ve)}({f2(self.ve0)}<{f2(self.vem)}>{f2(self.ve1)})' if self.de else ''
    # TODO: add support for agumenting/reporting actual nozzle output.
    rstr = f' re={f2(self.re)}'
    if self.isdraw:
      mstr = f'draw'
    elif self.ishopup:
      mstr = f'hopup'
    elif self.ishopdn:
      mstr = f'hopdn'
    elif self.isretract:
      mstr = f'retract'
    elif self.isrestore:
      mstr = f'restore'
    else:
      mstr = f'move'
    return ''.join((mstr,lstr,pstr,hstr,estr,rstr))

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
  def _db(cls, eb, h):
    """ Get the bead diameter from the bead volume."""
    # h<=0 can happen for initializing moves or moves at the start of a new
    # layer before a hop, and h>Nd*2 can happen at the start of a print when
    # moving to the first layer, so always clamp h. For h>Nd*2 it spagetti's
    # instead of beading, so that's the effective upper limit on bead height.
    return dcircle(eb*cls.Fa/clamp(h, 0.01, cls.Nd*2))

  @classmethod
  def _eb(cls, h, db):
    """ Volume in filament length of bead."""
    return clamp(h, 0.01, cls.Nd*2) * acircle(db) / cls.Fa

  @classmethod
  def _el(cls, l, h, db):
    """ Volume in filament length of line."""
    return l * clamp(h, 0.01, cls.Nd*2) * db / cls.Fa

  @property
  def f(self):
    """ The 'f' mm/min speed of a move. """
    return round(self.v * 60)

  @cached_property
  def re(self):
    """Extruder ratio (not target ratio)."""
    return self.de/self.el

  @cached_property
  def dt0(self):
    """Duration of acceleration phase."""
    dt0 = self.dv0 / self.a0 if self.dl else self.dve0 / self.ae0 if self.de else 0.0
    assert dt0 >= 0.0, f'{self.dv0=} {self.a0=} {self.dl=} {self.dve0=} {self.ae0=} {self.de=} {self!r}'
    return dt0

  @cached_property
  def dt1(self):
    """Duration of deceleration phase."""
    dt1 = self.dv1 / self.a1 if self.dl else self.dve1 / self.ae1 if self.de else 0.0
    assert dt1 >= 0.0, f'{self.dv1=} {self.a1=} {self.dl=} {self.dve1=} {self.ae1=} {self.de=} {self!r}'
    return dt1

  @cached_property
  def dtm(self):
    """Duration of middle phase."""
    dtm = self.dlm/self.vm if self.vm else self.dem/self.vem if self.vem else abs(self.dz/self.v)
    assert dtm >= 0.0, f'{self.dlm=} {self.vm=} {self.dem=} {self.vem=} {self!r}'
    return dtm

  @cached_property
  def dt(self):
    """ Get the execution time for a move. """
    return self.dt0 + self.dtm + self.dt1

  @property
  def a0(self):
    """ The acceleration phase change in velocity. """
    # Line acceleration should never be limited by extruder acceleration, even with PA compensating moves.
    return self.Ap if self.dl else 0.0

  @property
  def a1(self):
    """ The acceleration phase change in velocity. """
    return -self.a0

  @property
  def vl(self):
    """ The target line speed of a move. """
    return self.v if self.dl else 0.0

  @cached_property
  def vm(self):
    """ The speed of the middle phase of a move. """
    vl, v0, v1 = self.vl, self.v0, self.v1
    # This is the limit velocity assuming constant acceleration at a/2.
    # This implements smoothed look-ahead to reduce spikey velocity.
    vm = sqrt((self.dl*self.Ap + v0**2 + v1**2)/2)
    # If vm is less than v0 or v1, we don't need acceleration at both ends,
    # don't set vm higher than the move's speed, and round it.
    vm = getnear(min(vl, max(v0, v1, vm)))
    assert v0 <= vm <= vl
    assert v1 <= vm <= vl
    return vm

  @property
  def dv0(self):
    """ The acceleration phase change in velocity. """
    return self.vm - self.v0

  @property
  def dv1(self):
    """ The decceleration phase change in velocity. """
    return self.v1 - self.vm

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

  @property
  def vz(self):
    """ The target vertical velocity of a move."""
    return self.dz/self.dl * self.v if self.dl else copysign(self.v, self.dz) if self.dz else 0.0

  @property
  def dz0(self):
    """ Change in z of acceleration phase."""
    return self.dz * self.dl0/self.dl if self.dl else 0.0

  @property
  def dz1(self):
    """ Change in z of decceleration phase."""
    return self.dz * self.dl1/self.dl if self.dl else 0.0

  @property
  def dzm(self):
    """ Change in z of mid phase."""
    return self.dz * self.dlm/self.dl if self.dl else self.dz

  @property
  def ae0(self):
    """ The extrusion acceleration rate at the start of a move. """
    dl, de = self.dl, self.de
    return de/dl * self.a0 if dl else copysign(self.Ae, de) if de else 0.0

  @property
  def ae1(self):
    """ The extrusion acceleration rate at the end of a move. """
    dl, de = self.dl, self.de
    return de/dl * self.a1 if dl else -copysign(self.Ae, de) if de else 0.0

  @property
  def ve(self):
    """ The target extrusion rate of a move. """
    return self.de/self.dl * self.v if self.dl else copysign(self.v, self.de) if self.de else 0.0

  @cached_property
  def ve0(self):
    """ The extrusion rate at the start of a move. """
    dl = self.dl
    # We assume pure retract/restore moves have to start with ve=0
    return self.de/dl * self.v0 if dl else 0.0

  @cached_property
  def ve1(self):
    """ The extrusion rate at the end of a move. """
    dl = self.dl
    # We assume pure retract/restore moves have to end with ve=0
    return self.de/dl * self.v1 if dl else 0.0

  @cached_property
  def vem(self):
    """ The extrusion rate at the middle of a move. """
    dl, de = self.dl, self.de
    # We assume there are no vertical extruding moves.
    if dl:
      return de/dl * self.vm
    ve, ve0, ve1 = self.ve, self.ve0, self.ve1
    # This is the limit velocity assuming constant acceleration at a/2.
    # This implements smoothed look-ahead to reduce spikey velocity.
    vem = copysign(sqrt((self.de*self.ae0 + ve0**2 + ve1**2)/2), self.de)
    # If vm is less than v0 or v1, we don't need acceleration at both ends,
    # don't set vm higher than the move's speed, and round it.
    if ve >= 0:
      vem = getnear(min(ve, max(ve0, ve1, vem)))
      assert 0 <= ve0 <= vem <= ve
      assert 0 <= ve1 <= vem <= ve
    else:
      vem = getnear(max(ve, min(ve0, ve1, vem)))
      assert ve <= vem <= ve0 <= 0
      assert ve <= vem <= ve1 <= 0
    return vem

  @property
  def dve0(self):
    return self.vem - self.ve0

  @property
  def dve1(self):
    return self.ve1 - self.vem

  @cached_property
  def de0(self):
    """ The acceleration phase change in e. """
    de0 = self.de * self.dl0/self.dl if self.dl else 0.5*(self.ve0+self.vem)*self.dt0
    assert 0 <= de0 <= self.de or self.de <= de0 <= 0, f'0.0 <= {de0=} <= {self.de=}, {self.dt0=} {self!r}'
    return de0

  @cached_property
  def de1(self):
    """ The decceleration phase change in e. """
    de1 = self.de * self.dl1/self.dl if self.dl else 0.5*(self.vem + self.ve1)*self.dt1
    assert 0 <= de1 <= self.de or self.de <= de1 <= 0, f'0 <= {de1=} <= {self.de=}, {self.dt1=} {self!r}'
    return de1

  @cached_property
  def dem(self):
    """ The mid phase change in e. """
    dem = self.de - self.de0 - self.de1
    assert 0 <= dem <= self.de or self.de <= dem <= 0, f'0 <= {dem=} <= {self.de=}, {self.de0=} {self.de1=} {self!r}'
    return dem

  @property
  def dd(self):
    """ The total distance of a move.

    This is a combined total travel distance that moves can be averaged over.
    It is the combined horizontal line distance and vertical distance.
    Retracts and restores are considered to have a line distance matching the
    volume of their bead. Hops can have dz and dl for slope/spiral hops.
    """
    dl, dz = self.dl if self.isgo else 0.0 if self.dz else pi/4*self.w, self.dz
    return sqrt(dl*dl + dz*dz)

  @property
  def vd(self):
    """ The target total distance velocity."""
    dt = m.dl/m.v if m.isgo else m.dz/m.v if m.dz else m.de/m.v
    return self.dd/dt

  @property
  def db(self):
    """ The target extrusion bead diameter. """
    return self.w * self.r

  @cached_property
  def eb(self):
    """ The extrusion bead volume not including r under/over extrusion. """
    eb = self._eb(self.h, self.w)
    assert eb >=0, f'{eb=} {self!r}'
    return eb

  @cached_property
  def el(self):
    """ The line volume not including r under/over extrusion. """
    # Note this is the bead volume for adjusts and the extra dz bead for non moving hops.
    el = self._el(self.dl, self.h, self.w) if self.isgo else self._eb(abs(self.dz), self.w) if self.dz else self.eb
    assert el >=0, f'{el=} {self!r}'
    return el

  @property
  def el0(self):
    """ The start line volume not including r under/over extrusion. """
    el0 = self.el*(self.dl0/self.dl if self.isgo else self.dz0/self.dz if self.dz else self.de0/self.de)
    assert el0 >=0, f'{el0=} {self!r}'
    return el0

  @property
  def el1(self):
    """ The end line volume not including r under/over extrusion. """
    el1 = self.el*(self.dl1/self.dl if self.isgo else self.dz1/self.dz if self.dz else self.de1/self.de)
    assert el1 >=0, f'{el1=} {self!r}'
    return el1

  @property
  def elm(self):
    """ The mid line volume not including r under/over extrusion. """
    elm = self.el*(self.dlm/self.dl if self.isgo else self.dzm/self.dz if self.dz else self.dem/self.de)
    assert elm >=0, f'{elm=} {self!r}'
    return elm

  @property
  def vt(self):
    """Target nozzle extrusion velocity."""
    return self.et/self.dl * self.v if self.dl else self.et/self.de * self.v if self.de else 0.0

  @cached_property
  def et(self):
    """Target nozzle extrusion."""
    et = self.el * self.r
    assert et >=0, f'{et=} {self!r}'
    return et

  @property
  def et0(self):
    """Target nozzle extrusion during acceleration phase."""
    et0 = self.el0 * self.r
    assert et0 >=0, f'{et0=} {self!r}'
    return et0

  @property
  def et1(self):
    """Target nozzle extrusion during deceleration phase."""
    et1 = self.el1 * self.r
    assert et1 >=0, f'{et1=} {self!r}'
    return et1

  @property
  def etm(self):
    """Target nozzle extrusion during middle phase."""
    etm = self.elm * self.r
    assert etm >=0, f'{etm=} {self!r}'
    return etm

  @cached_property
  def isgo(self):
    """ Is this a horizontal move? """
    return self.dx or self.dy

  @cached_property
  def isdraw(self):
    # Note there can be sloped draws.
    return self.isgo and self.r > 0.0

  @property
  def ismove(self):
    # Any dz means it's a slope or spiral hop.
    return self.isgo and self.r <= 0 and not self.dz

  @property
  def isadjust(self):
    return self.de and not (self.isgo or self.dz)

  @property
  def isretract(self):
    return self.de < 0 and not (self.isgo or self.dz)

  @property
  def isrestore(self):
    return self.de > 0 and not (self.isgo or self.dz)

  @property
  def ishop(self):
    # Note there can be slope or spiral hops.
    return self.dz and self.r <= 0

  @property
  def ishopup(self):
    return self.dz > 0 and self.r <= 0

  @property
  def ishopdn(self):
    return self.dz < 0 and self.r <= 0

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

  def extruderv(self, m):
    """ Calculate the max velocity for extruder 'jerk' between two moves.

    This calculates the max velocity between two moves so extrusion velocity
    ve changes instantly by at most `Je`. We ignore dz as insignificant and
    use horizontal velocities only. This matches the marlin 'jerk' mechanism
    for handling sudden extruder velocity changes.
    """
    # vmax = Je/abs(self.de/self.dl - m.de/m.dl)
    dedl = abs(self.de*m.dl - m.de*self.dl)
    vmax = self.Je*self.dl*m.dl/dedl if dedl else inf
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
        (self.h,self.w,self.r) == (m.h,m.w,m.r) and
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

    if m0 is set the v0 start velocity is set to the max corner/extruder velocity.
    If m1 is set the v1 end velocity is set to the max corner/extruder velocity.
    If v0, or v1 are still not set they default to 0.0.
    """
    if m0 is not None: v0 = min(m0.cornerv(self), m0.extruderv(self))
    if m1 is not None: v1 = min(self.cornerv(m1), self.extruderv(m1))
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


class CmdStats(object):
  """ Stats for a sample of moves."""
  def __init__(self, name):
    self.name = name
    self.c = stats1.Sample() # cmd stats (dt,dd,dl,dz,de,el,et,en,v,v0,vm,v1 per cmd).
    self.l = stats1.Sample() # cmd line stats (vl,vz,ve,vn,h,w,r per dd).
    self.o = stats1.Sample() # out line stats (vl,vz,ve,vn,re,rt,rn,er per dd).
    self.r = stats1.Sample() # extruder stats (re,rt,rn,er per el).

  def _atdev(self, s, d):
    """Get a Sample's values at d stddevs clipped between min and max."""
    return stats1.clip(s.atdev(d), s.min, s.max)

  def __str__(self):
    num = self.c.num
    dt, dd, dl, dz, de, el, et, en, *_  = self.c.sum
    er = (en-et)/el
    assert isneareq(self.l.num, dd)
    assert isneareq(self.o.num, dd)
    assert isneareq(self.r.num, el)
    cmdsums = f'totals {num=} dt={ftime(dt)} {dd=:.1f} {dl=:.1f} {dz=:.2f} {de=:.3f} {el=:.3f} {et=:.3f} {en=:.3f} {er=:+.3f}'
    *_, v, v0, vm, v1 = self.c.avg
    vl, vz, ve, vn, h, w, r = self.l.avg
    vl0, vz0, ve0, vn0, h0, w0, r0 = self._atdev(self.l,-2)
    vl1, vz1, ve1, vn1, h1, w1, r1 = self._atdev(self.l,0)
    vl2, vz2, ve2, vn2, h2, w2, r2 = self._atdev(self.l,2)
    cmdlines = f'move args: {v=:.0f}({v0:.0f}<{vm:.0f}>{v1:.0f}) {h=:.2f}({h0:.2f}<{h1:.2f}<{h2:.2f}) {w=:.2f}({w0:.2f}<{w1:.2f}<{w2:.2f}) {r=:.2f}({r0:.2f}<{r1:.2f}<{r2:.2f})'
    vlstr=f' {vl=:.0f}({vl0:.0f}<{vl1:.0f}<{vl2:.0f})' if vl else ''
    vzstr=f' {vz=:.1f}({vz0:.1f}<{vz1:.1f}<{vz2:.1f})' if vz else ''
    vestr=f' {ve=:.3f}({ve0:.3f}<{ve1:.3f}<{ve2:.3f})' if ve else ''
    vnstr=f' {vn=:.3f}({vn0:.3f}<{vn1:.3f}<{vn2:.3f})' if vn else ''
    cmdvels = f'target speeds per dd:{vlstr}{vzstr}{vestr}{vnstr}'
    vl, vz, ve, vn, re, rt, rn, er = self.o.avg
    vl0, vz0, ve0, vn0, re0, rt0, rn0, er0 = self._atdev(self.o,-2)
    vl1, vz1, ve1, vn1, re1, rt1, rn1, er1 = self._atdev(self.o,0)
    vl2, vz2, ve2, vn2, re2, rt2, rn2, er2 = self._atdev(self.o,2)
    vlstr=f' {vl=:.0f}({vl0:.0f}<{vl1:.0f}<{vl2:.0f})' if vl else ''
    vzstr=f' {vz=:.1f}({vz0:.1f}<{vz1:.1f}<{vz2:.1f})' if vz else ''
    vestr=f' {ve=:.3f}({ve0:.3f}<{ve1:.3f}<{ve2:.3f})' if ve else ''
    vnstr=f' {vn=:.3f}({vn0:.3f}<{vn1:.3f}<{vn2:.3f})' if vn else ''
    actvels = f'actual speeds per dd:{vlstr}{vzstr}{vestr}{vnstr}'
    restr=f' {re=:.2f}({re0:.2f}<{re1:.2f}<{re2:.2f})' if re else ''
    rtstr=f' {rt=:.2f}({rt0:.2f}<{rt1:.2f}<{rt2:.2f})' if rt else ''
    rnstr=f' {rn=:.2f}({rn0:.2f}<{rn1:.2f}<{rn2:.2f})' if rn else ''
    erstr=f' {er=:+.2f}({er0:+.2f}<{er1:+.2f}<{er2:+.2f})'
    lratios = f'volume ratios per dd:{restr}{rtstr}{rnstr}{erstr}'
    re, rt, rn, er = self.r.avg
    re0, rt0, rn0, er0 = self._atdev(self.r,-2)
    re1, rt1, rn1, er1 = self._atdev(self.r,0)
    re2, rt2, rn2, er2 = self._atdev(self.r,2)
    restr=f' {re=:.2f}({re0:.2f}<{re1:.2f}<{re2:.2f})' if re else ''
    rtstr=f' {rt=:.2f}({rt0:.2f}<{rt1:.2f}<{rt2:.2f})' if rt else ''
    rnstr=f' {rn=:.2f}({rn0:.2f}<{rn1:.2f}<{rn2:.2f})' if rn else ''
    erstr=f' {er=:+.2f}({er0:+.2f}<{er1:+.2f}<{er2:+.2f})'
    eratios = f'volume ratios per el:{restr}{rtstr}{rnstr}{erstr}'
    return f"""{self.name} stats:
  {cmdsums}
  {cmdlines}
  {cmdvels}
  {actvels}
  {lratios}
  {eratios}"""

  def addinc(self, m, dt, dl, dz, de, en):
    """ Update the stats for part of a move."""
    # get the fraction of a move from the distance ratios.
    s = (dl/m.dl if dl else dz/m.dz if dz else de/m.de)
    # scale dist and volumes by the fraction of the move completed.
    dd = s * m.dd; el = s * m.el; et = s * m.et
    assert 0<dd<500 and 0<el<500 and 0<=et<500 and 0<=en<500, f'{dd=} {el=} {et=} {en=} for {m!r}'
    vl, vz, ve, vn = dl/dt, dz/dt, de/dt, en/dt
    re, rt, rn = de/el, et/el, en/el
    er = rn - rt
    self.o.add((vl, vz, ve, vn, re, rt, rn, er), dd)
    self.r.add((re, rt, rn, er), el)

  def addcmd(self, m, en):
    """ Update the stats for a whole move."""
    self.c.add((m.dt,m.dd,m.dl,m.dz,m.de,m.el,m.et,en,m.v,m.v0,m.vm,m.v1))
    self.l.add((m.vl,m.vz,m.ve,m.vt,m.h,m.w,m.r), m.dd)


class PrintStats(object):
  """ Print stats collection and rendering."""

  def __init__(self):
    # Stats on speeds.
    self.total = CmdStats('total')
    self.draw = CmdStats('draw')
    self.move = CmdStats('move')
    self.hopup = CmdStats('hopup')
    self.hopdn = CmdStats('hopdn')
    self.retract = CmdStats('retract')
    self.restore = CmdStats('restore')

  def __str__(self):
    return f"""{self.total}
{self.draw}
{self.move}
{self.hopup}
{self.hopdn}
{self.retract}
{self.restore}"""

  def addinc(self, m, dt, dl, dz, de, en):
    """ Update the line stats for part of a move."""
    self.total.addinc(m, dt, dl, dz, de, en)
    if m.isdraw: self.draw.addinc(m, dt, dl, dz, de, en)
    elif m.ismove: self.move.addinc(m, dt, dl, dz, de, en)
    elif m.ishopup: self.hopup.addinc(m, dt, dl, dz, de, en)
    elif m.ishopdn: self.hopdn.addinc(m, dt, dl, dz, de, en)
    elif m.isretract: self.retract.addinc(m, dt, dl, dz, de, en)
    elif m.isrestore: self.restore.addinc(m, dt, dl, dz, de, en)

  def addcmd(self, m, en):
    """ Update the stats for a whole move."""
    self.total.addcmd(m, en)
    if m.isdraw: self.draw.addcmd(m, en)
    elif m.ismove: self.move.addcmd(m, en)
    elif m.ishopup: self.hopup.addcmd(m, en)
    elif m.ishopdn: self.hopdn.addcmd(m, en)
    elif m.isretract: self.retract.addcmd(m, en)
    elif m.isrestore: self.restore.addcmd(m, en)


class GCodeGen(AttrDictMixin):
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
    pn: current nozzle pressure (property).
    pb: current bead backpressure (property).
    eb: current extruded bead volume.
    db: current extruded bead diameter (property).
    vn: current nozzle flow velocity (property).
    ve: current extruder velocity.
    vl: current line velocity.
    dt: the execution time of the last command.
    dl: the line length of the last command.
    de: the extruder change of the last command.
    en: the extruded volume of the last command.
    db0: the bead diameter at the start of the last move.
    eb0: the bead volume at the start of the last move.
    r0,rm,r1: the extrusion ratios of the last move (property)
    c: the last non-comment command.
    m: the last move command.
    d: the last move command that was a draw or restore.
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

  # Slicer Types for formating gcode.
  SlicerTypes = ('FlashPrint', 'OrcaSlicer')
  FlashPrint, OrcaSlicer = SlicerTypes

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
    return f'M651 S{round(fc*255)}' if fc else 'M652'

  def _setFe(self, fe):
    return f'M106' if fe == 1.0 else f'M106 S{round(fe*255)}' if fe else 'M107'

  def __init__(self,
      Slicer=FlashPrint,
      Te=210, Tp=50, Fe=1.0, Fc=1.0, fKf=0.0,
      Kf=0.0, Kb=0.0, Re=5,
      Vp=60, Vt=100, Vz=7, Ve=40, Vb=30,
      Lh=0.2, Lw=Nd, Lr=1.0,
      relext=False, optmov=False, segvel=False, fixvel=0,
      diff=False, verb=False,
      dynret=False, dynext=False):
    assert Slicer in self.SlicerTypes
    # Slicer gcode flavour to use.
    self.Slicer = Slicer
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

  def Pb(self, db):
    """Get the bead backpressure from db."""
    return max(0.0, self.Kb*min(db,2*self.Nd))

  def Pn(self, vn):
    """Get the nozzle pressure from nozzle extrude velocity."""
    return max(0.0, self.Kf*vn)

  def Pe(self, vn, db):
    """ Get pressure advance pe for target vn and db. """
    return self.Pn(vn) + self.Pb(db)

  @property
  def pn(self):
    """ The current nozzle pressure."""
    # Get nozzle pressure from pe minus pb.
    return max(0.0, self.pe - self.pb)

  @property
  def pb(self):
    """ The current bead backpressure."""
    return self.Pb(self.db)

  @property
  def db(self):
    """ Get the current bead diameter from the bead volume."""
    return self.GMove._db(self.eb, self.hl)

  @property
  def vn(self):
    """Get nozzle flow velocity from the nozzle pressure and back pressure."""
    return self.pn/self.Kf if self.Kf else self.ve if self.pe>=self.pb else 0.0

  @property
  def rn(self):
    """ The nozzle extrusion ratio for the last move."""
    assert isneareq(self.en, self.en0+self.enm+self.en1), f'{self.en=} {self.en0=} {self.enm=} {self.en1=}'
    assert isneareq(self.m.el, self.m.el0+self.m.elm+self.m.el1), f'{self.m.el=} {self.m.el0=} {self.m.elm=} {self.m.el1=}'
    return self.en/self.m.el if self.m else 0.0

  @property
  def rn0(self):
    """ The nozzle extrusion ratio for the accel phase of the last move."""
    return self.en0/self.m.el0 if self.m and self.m.el0 else self.rnm

  @property
  def rnm(self):
    """ The nozzle extrusion ratio for the mid phase of the last move."""
    return self.enm/self.m.elm if self.m and self.m.elm else 0.0

  @property
  def rn1(self):
    """ The nozzle extrusion ratio for the decel phase of the last move."""
    return self.en1/self.m.el1 if self.m and self.m.el1 else self.rnm

  @property
  def hl(self):
    """ Get the current line height."""
    return self.m.h if self.m else self.layer.h

  @property
  def wl(self):
    """ Get the current line width."""
    return self.m.w if self.m else self.layer.w

  @property
  def rl(self):
    """ Get the current line extrusion ratio."""
    return self.m.r if self.m else self.layer.r

  @property
  def hd(self):
    """ Get the last draw line height."""
    return self.d.h if self.d else self.layer.h

  @property
  def wd(self):
    """ Get the last draw line width."""
    return self.d.w if self.d else self.layer.w

  @property
  def rd(self):
    """ Get the last draw line extrusion ratio."""
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

  def iswait(self, cmd):
    return self.iscmd(cmd) and cmd[0] == 'G4'

  def resetfile(self, nextcmds=None, finalize=False):
    """ Reset all file state variables.

    This gets ready to process or finalize a queue of nextcmds or new commands.
    """
    # Zero all positions so the first move's deltas are actually absolute.
    self.x, self.y, self.z, self.e, self.f = 0.0, 0.0, 0.0, 0.0, 0
    self.f_t = self.f_l = self.f_e  = 0.0
    # We assume prints start retracted.
    self.pe = -self.Re
    self.eb = self.vl = self.ve = self.h = 0.0
    self.prevcmds, self.nextcmds, self.finalize = deque(), nextcmds, finalize
    self.c = self.m = self.d = None
    self.resetlayer()
    self.resetstats()

  def resetstats(self):
    self.stats = PrintStats()

  def resetlayer(self):
    """ Reset all layer state variables. """
    # Set layer to a default base layer for any layer-less moves.
    self.layer = Layer(
        0, 0.0, self.Vp, self.Vt, self.Vz, self.Ve, self.Vb, self.Lh, self.Lw, self.Lr)
    self.l_t = self.l_l = self.l_e = 0.0
    self.resetstate()

  def resetstate(self):
    """ Reset all move state variables. """
    self.dt = self.dl = self.de = self.en = 0.0


  def iterstate(self, dt, dl, de, dh, dvl, dve):
    """ Increment the state for a small dt interval. """
    # Note filament doesn't get sucked back into the nozzle if the pressure
    # advance is negative, and the bead cannot give negative backpressure.
    # Get filament extruded into the nozzle and height.
    self.pe += de
    en = self.pn * dt/(self.Kf + dt)          # filament extruded out the nozzle.
    el = min(self.GMove._el(dl, self.hl, self.db), self.eb + en)  # fillament removed from bead to line,
    self.pe -= en
    self.eb = max(0.0, self.eb + en - el)  # Need to make sure float errors don't make this negative.
    self.dt += dt
    self.dl += dl
    self.de += de
    self.en += en
    self.vl += dvl
    self.ve += dve
    self.h += dh
    self.stats.addinc(self.m, dt, dl, dh, de, en)

  def incstate(self, dt, dl=0.0, de=0.0, dh=0.0, dvl=0.0, dve=0.0):
    """ Increment the state for a large dt interval and return en. """
    # Keep the expected final state to check against and use later.
    dt1, dl1, de1, h1, vl1, ve1 = self.dt+dt, self.dl+dl, self.de+de, self.h+dh, self.vl+dvl, self.ve+dve
    # Only do a fine-grained simulation when finalizing.
    en0, t, idt = self.en, dt, 0.001
    if dl or de:
      # integrate a horizontal or extruding move.
      al, ae, dh_dle = dvl/dt, dve/dt, dh/dl if dl else dh/de
      while t > 0.0:
        idt = min(idt, t)
        idvl = al * idt
        idve = ae * idt
        idl = (self.vl + idvl/2) * idt  # Do trapezoidal integration for dl.
        ide = (self.ve + idve/2) * idt  # Do trapezoidal integration for de.
        idh = idl * dh_dle if dl else ide * dh_dle
        self.iterstate(idt, idl, ide, idh, idvl, idve)
        t -= idt
    else:
      # integrate a vertical move.
      vz = dh/dt
      while t > 0.0:
        idt = min(idt, t)
        idh = vz * idt
        self.iterstate(idt, 0.0, 0.0, idh, 0.0, 0.0)
        t -= idt
    assert isneareq(self.dt, dt1), f'{self.dt=} {dt1=}'
    assert isneareq(self.dl, dl1), f'{self.dl=} {dl1=}'
    assert isneareq(self.de, de1), f'{self.de=} {de1=}'
    assert isneareq(self.h, h1), f'{self.h=} {h1=}'
    assert isneareq(self.vl, vl1), f'{self.vl=} {vl1=}'
    assert isneareq(self.ve, ve1), f'{self.ve=} {ve1=}'
    # Use the expected final states to remove small integration errors.
    self.dt, self.dl, self.de, self.h, self.vl, self.ve = dt1, dl1, de1, h1, vl1, ve1
    return self.en - en0

  def inc(self, code):
    if self.islayer(code):
      self.inclayer(code)
    elif self.ismove(code):
      self.incmove(code)
    elif self.iscmd(code):
      self.inccmd(code)
    if not (self.isstr(code) or self.isbin(code)):
      self.c = code

  def inct(self, dt):
    """ Increment time. """
    self.f_t += dt
    self.l_t += dt
    self.dt = dt

  def inclayer(self, l):
    """ Increment state for starting a Layer. """
    # reset stats before layer 1 to exlude pre-extrudes.
    if self.layer.n < 1:
      self.resetstats()
    self.resetlayer()
    self.layer = l
    self.h = getnear(self.z - l.z)
    assert self.h >= 0.0, f'{self.z=} {l}'

  def incmove(self, m):
    """ Increment state for executing a Move. """
    # Check that the start states are as expected.
    assert isneareq(self.h, self.z - self.layer.z), f'{self.h=} != {self.z - self.layer.z}'
    assert isneareq(self.vl, m.v0), f'{self.vl=} != {m.v0=} for {m} at {ftime(self.f_t)}'
    assert isnearle(abs(self.ve - m.ve0), self.GMove.Je), f'{self.ve=} {m.ve0=} {self.m!s} -> {m!s}'
    self.resetstate()
    # Update starting state vars and the current move.
    self.db0, self.eb0, self.m = self.db, self.eb, m
    # Only run the detailed simulation when finalizing the output.
    if self.finalize:
      # Save the expected final h to assert against later.
      h1 = self.h + m.dz
      # "jerk" ve from the previous ve1 to ve0.
      self.ve = m.ve0
      # Update the state for the three move phases.
      self.en0 = self.incstate(m.dt0, m.dl0, m.de0, m.dz0, m.dv0, m.dve0) if m.dt0 else 0.0
      self.enm = self.incstate(m.dtm, m.dlm, m.dem, m.dzm, m.dvm, m.dvem) if m.dtm else 0.0
      self.en1 = self.incstate(m.dt1, m.dl1, m.de1, m.dz1, m.dv1, m.dve1) if m.dt1 else 0.0
      self.stats.addcmd(m, self.en)
      # Check that the end states are as expected.
      assert isneareq(self.dt, m.dt)
      assert isneareq(self.dl, m.dl)
      assert isneareq(self.de, m.de)
      assert isneareq(self.vl, m.v1)
      assert isneareq(self.ve, m.ve1)
      assert isneareq(self.h, h1)
    # Update and round positions to remove floating point errors.
    self.x = getnear(self.x+m.dx)
    self.y = getnear(self.y+m.dy)
    self.z = getnear(self.z+m.dz)
    self.e = getnear(self.e+m.de)
    self.f = m.f
    self.f_l += m.dl
    self.f_e += self.en
    self.l_l += m.dl
    self.l_e += self.en
    self.vl = m.v1
    self.ve = m.ve1
    self.h = getnear(self.z - self.layer.z)
    if m.r > 0 : self.d = m
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
    if self.Slicer == self.FlashPrint:
      # layer 0 is the preExtrude
      return self.fcmt('{"layer" if layer.n else "preExtrude"}:{layer.h:.2f}')
    else:
      return '\n'.join([
      self.fcmt('LAYER_CHANGE'),
      self.fcmt('Z:{f3(layer.z+layer.h)}'),
      self.fcmt('HEIGHT:{f3(layer.h)}')])

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
      cmd = f'{cmd:40s}; {{dt=:.3f}} {{m}} {{rn=:.2f}}({{rn0:.2f}}<{{rnm:.2f}}>{{rn1:.2f}})'
    return cmd

  def _farg_FlashPrint(self, k, v):
    # This formats always with the number of digits for the desired precision.
    # Number of decimals for cmd arguments.
    argfmt = dict(X='.2f',Y='.2f',Z='.3f',E='.4f',F='d').get(k,'')
    return k + (v if isinstance(v, str) else f'{v:{argfmt}}')

  def _farg_OrcaSlicer(self, k, v):
    # This formats using the least digits necessary for the desired precision.
    argres = dict(X=3,Y=3,Z=3,E=5,F=0).get(k,6)
    return k + (v if isinstance(v, str) else fzp(v, argres))

  def farg(self, k, v):
    # This replaces itself with the right formatter after first use.
    if self.Slicer == 'FlashPrint':
      self.farg=self._farg_FlashPrint
    else:
      self.farg=self._farg_OrcaSlicer
    return self.farg(k,v)

  def fcmd(self, cmd, **kwargs):
    """ Format a gcode command.

    Args set to None like X=None are omitted. Args set to a string value like
    X='100' use the string value without any re-formating like "X100". Args
    set to an empty string like A='' will give the arg without a value like
    "A". Numeric value arguments will generally be formatted with decimal
    places appropriate for their required precision.
    """
    args = ' '.join(self.farg(k,v) for k,v in kwargs.items() if v is not None)
    return f'{cmd} {args}' if args else cmd

  def fstr(self, str):
    """ Format a string as an f-string with attributes of self."""
    try:
      return fstr(str, locals=self)
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
        return self.fcmt('{ftime(f_t)} {f_l:.1f}mm: ' + txt)

  def faddmove(self, m):
    """ Do a final format and add of a Move. """
    if self.dynret and m.isretract:
      # Adjust retraction to Re and also relieve advance pressure.
      m = m.change(de=-self.Re - self.pe)
      # self.log(f'adjusted {m0} -> {m}')
      # If de is zero, skip adding this.
      if not m.de:
        self.log('skipping empty retract.')
        return
    elif self.dynret and m.isrestore:
      # Get the pe and eb for the next draws.
      pe, eb = self._getNextPeEb()
      # Adjust existing restore and add the starting bead.
      m = m.change(de=pe - self.pe + eb)
      # self.log(f'adjusted {m0} -> {m}')
      # If de is zero, skip adding this.
      if not m.de:
        self.log('skipping empty restore.')
        return
    elif self.dynext and self.layer.n and (m.isdraw or (self.pe > 0.0 and m.isgo)):
      # Skip adjusting pressures for preextrudes before layer 1.
      # Adjust extrusion for any draw or non-retracted move.
      # Get the ve and db before and after this move.
      ve0 = self.ve
      ve1 = next((m1.ve0 for m1 in self.nextcmds if self.ismove(m1)), 0.0)
      # Get the pe and eb for the next moves.
      pe, eb = self._getNextPeEb(m)
      # Get the pressure adjustment required, scaled by P=0.8.
      # TODO: Maybe use a PID or PD controller instead.
      dpe = 0.8*(pe + eb - (self.pe + self.eb))
      # Only adjust the move if the change required is large enough.
      if abs(dpe) > 0.01:
        # Calculate the adjusted m.de from the average ve to accumulate the
        # required pressure change over the next max(m.La/2, m.dl).
        de = m.de + dpe*m.dl/max(m.La/2, m.dl)
        Je = 0.8*m.Je
        if m.v0 > 0.0:
          # Constrain de within jerk limits of previous move's final ve.
          mve0 = de/m.dl*m.v0
          mve0 = clamp(mve0, ve0-Je, ve0+Je)
          de = mve0*m.dl/m.v0
        if m.v1 > 0.0:
          # Constrain de within jerk limits to next move's initial ve.
          mve1 = de/m.dl*m.v1
          mve1 = clamp(mve1, ve1 - Je, ve1 + Je)
          de = mve1*m.dl/m.v1
        m0,m = m, m.change(de=de)
        assert abs(m.ve0) <= m.Ve
        assert abs(m.vem) <= m.Ve
        assert abs(m.ve1) <= m.Ve
        assert abs(m.ve0 - ve0) <= m.Je, f'{m.ve0=} {ve0=} {m.ve0-ve0=} {m}'
        assert abs(m.ve1 - ve1) <= m.Je, f'{m.ve1=} {ve1=} {m.ve1-ve1=} {m}'
        #self.log(f'adjusted for target {pe=:.4f} {m0} -> {m}')
    if self.e + m.de > 1000.0:
      # Reset the E offset before it goes over 1000.0
      self.cmd('G92', E=0)
    fmove = self.fmove(m)
    self.incmove(m)
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
    """ Add any kind of command code. """
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
    # convert G0 and G1 cmds into moves.
    if cmd in ('G0', 'G1'):
      # Use the previous f value and machine limits when inserting raw move cmds.
      kwargs['v'] = min(kwargs.pop('F', self.f)/60, self.GMove.Vp)
      kwargs = {k.lower():v for k,v in kwargs.items()}
      self.move(**kwargs)
    else:
      self.add((cmd, kwargs))

  def cmt(self, cmt):
    """ Add a comment. """
    self.add(self.fcmt(cmt))

  def log(self, txt):
    """ Conditionally add a log comment. """
    #print(self.fstr(txt))
    self.add(self.flog(txt))

  def startFile(self):
    self.add(self.startcode)
    # Note Flashprint always generates and assumes pre-extrude has h=0.2mm.
    self.startLayer(0, h=0.2)

  def endFile(self):
    self.endLayer()
    self.filestats()
    self.add(self.endcode)

  def startLayer(self, n=None, z=None,
      Vp=None, Vt=None, Vz=None, Ve=None, Vb=None,
      h=None, w=None, r=1.0):
    if n is None: n = self.layer.n + 1
    if z is None: z = getnear(self.layer.z + self.layer.h) if n > 1 else 0.0
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
    self.log('layer={layer.n} z={layer.z:.2f} h={layer.h:.2f} w={layer.w:.2f} r={layer.r:.2f}')

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
    if h is None: h = getnear(self.h + dz, 3)
    if w is None: w = self.wl
    # r is scaled by layer.r so changing layer.r scales all draws.
    r = self.layer.r * r
    de = e - self.e if e is not None else de if de is not None else self.GMove._el(dl,w,h*r)
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
      if self.Slicer == self.FlashPrint and m.r > 0 and not isneareq(m.r, self.rd, 2):
        self.cmt(f'extrude_ratio:{f2(m.r)}')
      if self.Slicer != self.FlashPrint and m.r > 0 and not isneareq(m.w, self.wd):
        self.cmt(f'WIDTH:{getnear(m.w+self.layer.h*(1-pi/4))}')
      self.add(m)
    else:
      self.log(f'skipping empty move {m!r}')

  def draw(self, x=None, y=None, z=None, h=None, w=None, r=None, **kwargs):
    """ Draw a line.

    This is the same as move() but with the default previous draw or layer r
    value so it draws instead of travels by default.
    """
    # For defaults use last move's h,w,r if it was a draw, else use layer's.
    if h is None: h = self.hl if self.rl else self.layer.h
    if w is None: w = self.wl if self.rl else self.layer.w
    if r is None: r = self.rl if self.rl else self.layer.r
    # Note we need r before scaling by layer.r.
    self.move(x, y, z, h=h, w=w, r=r/self.layer.r, **kwargs)

  def retract(self, e=None, de=None, ve=None, s=1.0, pe=None):
    """Do a retract.

    This is a retraction that by default relieves any pressure and retracts to
    pe which defaults to -self.Re. Setting de or e will do a fixed retraction
    ignoring the current pressure. Note compensation for advance pressure is
    adjusted in post processing if enabled.
    """
    #if pe is None: pe = -self.Re
    #if de is None: de = -self.pe + pe
    if de is None:
      de = -self.Re if pe is None else -self.pe + pe
    # If dynamic retraction is enabled and de is zero or more, set de to a
    # tiny retraction so that it doesn't get optimized away or seen as a
    # restore and can be dynamically adjusted later.
    if self.dynret and de >= 0: de = -0.00001
    self.move(e=e, de=de, v=ve, s=s)

  def restore(self, e=None, de=None, vb=None, s=1.0, pe=None, h=None, w=None, r=None):
    """ Do a restore.

    This is a restore that by default reverts any retraction and restores
    pressure to pe which defaults to 0.0. Setting e or de will do a fixed
    restore ignoring current retraction. Note adding starting dots and
    compensating for advance pressure is added in post-processing if enabled.
    """
    #if pe is None: pe = 0.0
    #if de is None: de = -self.pe + pe
    if de is None:
      de = self.Re if pe is None else -self.pe + pe
    # If dynamic retraction is enabled and de is zero or less, set de to a
    # tiny restore so that it doesn't get optimized away or seen as a
    # retraction and can be dynamicly adjusted later.
    if self.dynret and de <= 0: de = 0.00001
    self.draw(e=e, de=de, v=vb, s=s, h=h, w=w, r=r)

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

  def _getNextPeEb(self,m=None):
    "Get pe and eb for the next draws including m."
    # Get the exponentially decaying average pe and eb for the future moves up
    # to the next retract/restore. This needs to calculate the increasing
    # decay as we go further into the future. Decay can be over time,
    # distance, or both, however, we just use decay over distance because
    # that's what matters for the print output quality.
    Tl = self.GMove.La  # the decay time-constant acceleration distance.
    pe = eb = 0.0     # the averaged pe and eb.
    K = 1.0           # the accumulated decay multiplier.

    def iterdl(ve, db, meb, dl):
      nonlocal K, pe, eb
      # accumulate decayed averages and update decay multiplier.
      a = K * dl/(Tl + dl)
      pe += a*self.Pe(ve, db)
      eb += a*meb
      K -= a

    c,l,t,m1 = 0, 0, 0, m
    def iterm(m):
      nonlocal c,l,t,m1
      # accumulate acceleration/deceleration phases using average ve.
      if m.dt0: iterdl((m.ve0 + m.vem)/2, m.db, m.eb, m.dl0)
      if m.dtm: iterdl(m.vem, m.db, m.eb, m.dlm)
      if m.dt1: iterdl((m.vem + m.ve1)/2, m.db, m.eb, m.dl1)
      # Note we use the last draw's middle phase for the last values because
      # Orca sometimes puts wipe-moves before retracts, and we want the
      # pressure of the last move, not the zero it tries and fails to
      # decelerate to before the retract.
      if m.isdraw:
        m1 = m
      c+=1
      l+=m.dl
      t+=m.dt

    if m:
      iterm(m)
    for m in (mi for mi in self.nextcmds if self.ismove(mi)):
      if m.isadjust or K < 0.05: break
      iterm(m)
    if m1:
      # treat last draw's middle pe and eb as lasting forever.
      pe += K*self.Pe(m1.vem, m1.db)
      eb += K*m1.eb
    #self.log(f'target {pe=:.4f} {eb=:.4f} for next {Tl:.1f}mm from {l:.1f}mm and {t:.3f}s of {c} moves to {m1}')
    return pe, eb

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
    # Only go through movement related commands using None for waits.
    cmds = (m for m in self.prevcmds if self.ismove(m) or self.iswait(m))
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
    self.log('file finished with {layer.n} layers.')
    self.log('{str(stats).replace(nl,nl+";")}')

  def preExt(self,x0,y0,x1,y1,le=120.0, lw=60.0, m=10.0, ve=20.0, vw=20, h=0.2, r=4):
    """Preextrude around a box with margin m.

    This extrudes for `le` mm at speed `ve` mm/s, then wipes for another `lw`
    mm at speed `vw` mm/s without extruding, and finally if we did a wipe
    another 10mm at 1mm/s to stick the last bit of wipe ooze to the plate.
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
        # final ooze to the plate.
        l, v, r = 10,1,0
    self.hopup()

  def fbox(self,x0,y0,x1,y1,shells=None,**kwargs):
    """ Fill a box with frame lines. """
    if shells is None: shells = inf
    w = kwargs.get('w', self.layer.w)
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
    w = kwargs.get('w', self.layer.w)
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
        ry=clamp(ry, y0, y1)
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
    w = kwargs.get('w', self.layer.w)
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
  cmdline.add_argument('-Slicer', choices=GCodeGen.SlicerTypes, default=GCodeGen.FlashPrint,
      help='Slicer GCode format for comments and hints.')
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
      Slicer=args.Slicer,
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
