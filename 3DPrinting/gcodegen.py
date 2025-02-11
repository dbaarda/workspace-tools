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
from math import e, pi, inf, sqrt
import vtext
from typing import NamedTuple
from functools import cached_property
from pprint import *

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


def isneareq(a,b,e=1e-6):
  """ Is a nearly equal to b?"""
  return abs(a-b) < e


def isnearle(a,b,e=1e-6):
  """ Is a nearly less than or equal to b? """
  return a <= b+e


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
    assert 0.0 <= v0 <= self.v
    assert 0.0 <= v1 <= self.v
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
    return (self.vm - self.v0) / self.a0

  @cached_property
  def dt1(self):
    """Duration of deceleration phase."""
    return (self.v1 - self.vm) / self.a1

  @cached_property
  def dtm(self):
    """Duration of middle phase."""
    return self.dlm/self.vm if self.vm else abs(self.dz if self.dz else self.de)/self.v

  @cached_property
  def dt(self):
    """ Get the execution time for a move. """
    return self.dt0 + self.dtm + self.dt1

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
    return self.isgo and self.de > 0

  @property
  def ismove(self):
    return self.isgo and self.de <= 0

  @property
  def isretract(self):
    return self.de < 0 and not self.isgo

  @property
  def isrestore(self):
    return self.de > 0 and not self.isgo

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
    # v = sqrt(a*R)
    return min(sqrt(self.Ap*R), self.v1, m.v0)

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

  def canjoin(self,m):
    # Always join moves together.
    return ((self.ismove and m.ismove) or
      # Join draws together if ...
      (self.isdraw and m.isdraw and
      # ... they have the same line settings, and ...
      (self.h,self.w,self.r) == (self.h,self.w,self.r) and
      # The path deviation from joining them is less than Gp.
      self.joind(m) < self.Gp))


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
    v = self.v
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


def ljoin(l):
  # TODO: make this and lfixv methods of GCodeGen and make ljoin not optimze
  # pre-extrudes or non-hopped moves (moves not crossing edges).
  i = m = None
  i1 = 0
  while i1<len(l):
    m1 = l[i1]
    if isinstance(m1, Move):
      if m is None:
        # This is the first m to join with.
        i, m = i1, m1
      else:
        # We have the next move m1, join m and m1 if they can.
        if m.canjoin(m1):
          # replace the i1 entry with the new joined m1.
          l[i1] = m1 = m.join(m1)
          # delete the old m entry, and adjust i1 to still point at m1
          del l[i]
          i1-=1
        # m1 is the next m to try and join.
        i, m = i1, m1
    i1+=1
  return l


def _rfixv(l, v1):
  """ Reverse fix velocities for deceleration limits. """
  for m in reversed(l):
    if isinstance(m, Move):
      # if v1 is already set to the target, stop.
      if m.v1 == v1:
        break
      m.fixv(v0=m.v0, v1=v1)
      v1 = m.v0


def lfixv(l, v0=0.0, v1=0.0):
  """ Fix velocities in a list of moves/comments/commands. """
  i = m = None
  for i1,m1 in enumerate(l):
    if isinstance(m1, Move):
      if m is None:
        # This is the first m to fix.
        i, m = i1, m1
      else:
        # We have the next move m1, fixv m.
        m.fixv(v0=v0, m1=m1)
        # if v0 changed, reverse fix velocities for acceleration limits.
        if m.v0 != v0:
          _rfixv(l[:i], m.v0)
        # v0 is the last m.v1, m1 is the next m to fix,
        v0, i, m = m.v1, i1, m1
  if m is not None:
    # Fix the last move.
    m.fixv(v0=v0, v1=v1)
  return l


def lsplit(l):
  """ Partition moves in a list of moves/comments/commands. """
  return [p for m in l for p in (m.split() if isinstance(m, Move) else [m])]


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

  Attributes:
    Te: Extruder temp (degC).
    Tp: Platform temp (degC).
    Fe: Extruder fan speed (0.0->1.0).
    Fc: Case fan speed (0.0->1.0).
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
    pe: current advance/retract length for extruder pressure.
    eb: current extruded bead volume for nozzle backpressure.
    f_t: current file execution time.
    f_l: current file printed line length.
    f_e: current file extruded filament length.
    l_t: current layer execution time.
    l_l: current layer printed line length.
    l_e: current layer extruded filament length.
    vl: move/draw velocity at the end of the last move.
    ve: extruder velocity at the end of the last move.
    vn: nozzle flow velocity at the end of the last move.
    db: extruded bead diameter at the end of the last move.
    x,y,z,e,f: Current x,y,z,e,f position values.
    h: current height above layer base (property).
    _log_t: time of last logging output for diff mode output.
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
;start gcode
M118 X{75.0:.2f} Y{75.0:.2f} Z{150.0:.2f} T0
M140 S{round(Tp)} T0
M104 S{round(Te)} T0
M104 S0 T1
M107
M900 K{Kf:.3f} T0
G90
{"M82"+nl if relext else ""}\
G28
M132 X Y Z A B
G1 Z50.000 F420
G161 X Y F3300
M7 T0
M6 T0
{_setFc(Fc)}
M108 T0
{_setFe(Fe)+nl if Fe else ""}\
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

  def _setFc(self, fc):
    return f'M652 S{round(fc*255)}' if fc else 'M651'

  def _setFe(self, fe):
    return f'M106' if fe == 1.0 else f'M106 S{round(fe*255)}' if fe else 'M107'

  def __init__(self,
      Te=210, Tp=50, Fe=1.0, Fc=1.0, Kf=0.0, Kb=0.0, Re=5,
      Vp=60, Vt=100, Vz=7, Ve=40, Vb=30,
      Lh=0.2, Lw=Nd, Lr=1.0,
      relext=False, dynret=False, dynext=False, optmov=False, fixvel=False, diff=False, verb=False):
    # Temp and Fan settings.
    self.Te, self.Tp, self.Fe, self.Fc = Te, Tp, Fe, Fc
    # Linear advance and retraction settings.
    self.Kf, self.Kb, self.Re = Kf, Kb, Re
    # Default velocities.
    self.Vp, self.Vt, self.Vz, self.Ve, self.Vb = Vp, Vt, Vz, Ve, Vb
    # Default line height, width, and extrusion ratio.
    self.Lh, self.Lw, self.Lr = Lh, Lw, Lr
    # Processing mode options.
    self.relext = relext
    self.dynret, self.dynext, self.optmov, self.fixvel, self.diff, self.verb = dynret, dynext, optmov, fixvel, diff, verb
    class GMove(Move, Fd=self.Fd, Nd=self.Nd): pass
    self.GMove = GMove
    self.resetfile()

  @property
  def h(self):
    """ Get the current height above the layer base."""
    return self.z - (self.layer.z if self.layer else 0.0)

  @property
  def w(self):
    """ Get the current default width."""
    return (self.layer.w if self.layer else self.Lw)

  @property
  def r(self):
    """ Get the current default extrusion ratio."""
    return (self.layer.r if self.layer else self.Lr)

  def resetfile(self):
    """ Reset all state variables. """
    # Zero all positions so the first move's deltas are actually absolute.
    self.x, self.y, self.z, self.e, self.f = 0.0, 0.0, 0.0, 0.0, 0
    self.f_t = self.f_l = self.f_e = self._log_t = 0.0
    self.pe = 0.0  # current extruder linear advance pressure or retracted length mm.
    self.eb = 0.0  # current extruded bead volume in mm of fillament.
    self.gcode = []
    self.resetlayer()

  def resetlayer(self):
    """ Reset all layer state variables. """
    self.layer = None
    self.l_t = self.l_l = self.l_e = 0.0
    self.lastdraw = None

  def inclayer(self, l):
    """ Increment state for starting a Layer. """
    self.resetlayer()
    self.layer = l

  def inct(self, dt):
    """ Increment time. """
    self.f_t += dt
    self.l_t += dt

  def incmove(self, m):
    """ Increment state for executing a Move. """
    dx,dy,dz,de,v = m[:5]
    v0, v1 = m.v0, m.v1
    dt, dl = m.dt, m.dl
    # Round positions to remove floating point errors.
    self.x = getnear(self.x+dx)
    self.y = getnear(self.y+dy)
    self.z = getnear(self.z+dz)
    self.e = getnear(self.e+de)
    self.f = m.f
    self.pe, self.eb, en = self._calc_pe_eb_en(m)
    self.f_l += dl if en else 0.0
    self.f_e += en
    self.l_l += dl if en else 0.0
    self.l_e += en
    self.vl = v1
    self.ve = de/dl*v1 if dl else de/dt if dt else 0.0
    self.vn = self._vn()
    self.db = self._db()
    self.en = en
    self.inct(dt)
    if m.isdraw: self.last_draw = m

  def inccmd(self, cmd):
    """ Handle special string gcode commands. """
    if m := re.match(r'G92(?: X([-+0-9.]+))?(?: Y([-+0-9.]+))?(?: Z([-+0-9.]+))?(?: E([-+0-9.]+))?', cmd):
      # This is a "set XYZE to ..." cmd, often used to keep E < 1000.
      x,y,z,e = m.groups()
      if x is not None: self.x = float(x)
      if y is not None: self.y = float(y)
      if z is not None: self.z = float(z)
      if e is not None: self.e = float(e)

  def fstr(self, str):
    """ Format a string as an f-string with attributes of self."""
    try:
      return fstr(str, locals=attrs(self))
    except Exception as e:
      raise Exception(f'Error formatting {str!r}.') from e

  def fcmd(self, cmd, **kwargs):
    """ Format a gcode command. """
    fmts = dict(X='.2f',Y='.2f',Z='.3f',E='.4f',F='d')
    args = ' '.join(f'{k}{v:{fmts.get(k,"")}}' for k,v in kwargs.items() if v is not None)
    return f'{cmd} {args}' if args else cmd

  def fmove(self, m):
    """ Format a Move as a gcode command. """
    dx,dy,dz,de,v = m[:5]
    # Don't include args in cmd if they are unchanged.
    # Include both x and y if either are changed.
    x = self.x + dx if dx or dy else None
    y = self.y + dy if dy or dx else None
    z = self.z + dz if dz else None
    # If relext or diff mode is on, use relative extrusion.
    if self.relext or self.diff:
      e = de if de else None
    else:
      e = self.e + de if de else None
    f = m.f
    if f == self.f: f = None
    cmd = self.fcmd('G1', X=x, Y=y, Z=z, E=e, F=f)
    if self.verb:
      cmd = f'{cmd:40s}; {m} r={{en/{m.el} if {m.el} else db/{m.w}:.2f}}'
    return cmd

  def fadd(self, code):
    """ Do a final format to string and add. """
    if isinstance(code, Layer):
      self.inclayer(code)
    elif isinstance(code, Move):
      out = self.fmove(code)
      # Flashprint interprets extrude_ratio comments as layer height changes.
      if code.isdraw and (not self.last_draw or code.h != self.last_draw.h):
        self.fadd(f';extrude_ratio:{round(code.h/self.layer.h, 4)}')
      self.incmove(code)
      self.fadd(out)
      self.flog('{pe=:.4f} {eb=:.4f} {ve=:.4f} {vn=:.4f} {vl=:.3f} {db=:.2f}')
    elif isinstance(code, int):
      # ints are waits.
      self.inct(code/1000)
      self.fadd(self.fcmd('G4', P=code))
    elif isinstance(code, str):
      self.inccmd(code)
      self.add(self.fstr(code.strip()).encode())
    else:
      # It's probably Flashprint binary header stuff.
      self.add(code)

  def add(self, code):
    """ Add code Layer, Move, Wait, or format-string instances. """
    if isinstance(code, Layer):
      self.inclayer(code)
    elif isinstance(code, Move):
      self.incmove(code)
    elif isinstance(code, int):
      # ints are waits.
      self.inct(code/1000)
    elif isinstance(code, str):
      self.inccmd(code)
    self.gcode.append(code)

  def cmd(self, cmd, **kwargs):
    self.add(self.fcmd(cmd,**kwargs))

  def cmt(self, cmt):
    self.add(f';{cmt}')

  def log(self, txt):
    #print(self.fstr(txt))
    if self.verb:
      if self.diff:
        self.cmt('{f_t - _log_t:.3f}: ' + txt)
      else:
        self.cmt('{ftime(f_t)}: ' + txt)

  def flog(self, txt):
    """Add a formatted log entry."""
    if self.verb:
      self.log(txt)
      # We need to fstr() reformat that last log entry...
      self.fadd(self.gcode.pop())
      # Update the log time.
      self._log_t = self.f_t

  def startFile(self):
    self.add(self.fstr(self.startcode))
    # Note Flashprint always generates and assumes pre-extrude has h=0.2mm.
    self.startLayer(0, h=0.2)

  def endFile(self):
    # Retract to -2.5mm for the next print.
    self.retract(de=-2.5)
    self.filestats()
    self.add(self.fstr(self.endcode))

  def startLayer(self, n=None, z=None,
      Vp=None, Vt=None, Vz=None, Ve=None, Vb=None,
      h=None, w=None, r=1.0):
    if n is None: n = self.layer.n + 1 if self.layer else 0
    if z is None: z = self.layer.z + self.layer.h if self.layer and self.layer.n else 0.0
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
    if l.n:
      self.cmt('layer:{layer.h:.2f}')
      self.log('layer={layer.n} h={layer.h:.2f} w={layer.w:.2f} r={layer.r:.2f}')

  def endLayer(self, **upargs):
    self.hopup(**upargs)
    self.layerstats()

  def _calc_pe(self, db, ve):
    """ Get steady-state pressure advance pe needed for target db and ve. """
    pb = max(0, self.Kb*db)
    pn = max(0, self.Kf*ve)
    #pn = self.Kf*ve * e**(-dt/self.Kf)
    return pn+pb

  def _db(self, eb=None, h=None):
    """ Get bead diameter from bead volume."""
    h = getnear(h) # remove accumulating floating point errors.
    if eb is None: eb = self.eb
    if h is None: h = self.h
    assert eb >= 0.0
    # h<=0 can happen for initializing moves or moves at the start of a new
    # layer before a hop.
    return 2*sqrt(eb*self.Fa/(h*pi)) if h > 0.0 else 0.0

  def _vn(self,pe=None,eb=None,h=None,ve=None):
    """Get nozzle flow velocity from pe, eb, h, and ve."""
    if pe is None: pe = self.pe
    if ve is None: ve = self.ve
    db = self._db(eb,h)                # diameter of bead
    pb = max(0, self.Kb*db)            # bead backpressure.
    pn = max(0, pe - pb)               # nozzle pressure.
    return pn/self.Kf if self.Kf else ve if pe>=pb else 0.0

  def _calc_pe_eb_en(self, m):
    """ Calculate updated pe, eb, and en values after a move."""
    if (self.Kf or self.Kb) and (self.pe>0 or self.eb>0 or m.de>=0):
      # Calculate change in pressure advance 'pe' and bead volume `eb` by
      # iterating over the time interval with a dt that is small relative to
      # the Kf timeconstant.
      dt, dl = m.dt, m.dl
      assert dt >= 0.0
      if not dl:
        # A hop, drop, retract, restore, or set speed.
        assert m.v0 == m.v1 == m.vm == 0.0
        v, vz, ve = 0.0, m.dz/dt if dt else 0.0, m.de/dt if dt else 0.0
        # There is only one phase.
        phases =((dt,0),)
      else:
        v, dz_dl, de_dl = m.v0, m.dz/dl, m.de/dl
        # Calculate accelerate, move, decelerate phases.
        phases = (m.dt0, m.a0),(m.dtm, m.am),(m.dt1, m.a1)
      pe, eb, el, l, h = self.pe, self.eb, 0.0, 0.0, self.h
      for t,a in phases:
        while t>0:
          dt = min(0.001, t)  # Use this iteration time interval.
          dv = a*dt
          dl = (v + dv/2)*dt  # Do trapezoidal integration for dl.
          assert dl >= 0.0
          v += dv
          l += dl
          assert isnearle(0.0, v)
          # Note filament doesn't get sucked back into the nozzle if the pressure
          # advance is negative, and the bead cannot give negative backpressure.
          # Get filament extruded into the nozzle and height.
          if dl:
            pe += de_dl * dl
            h += dz_dl * dl
          else:
            pe += ve * dt
            h += vz * dt
          db = self._db(eb, h)                # diameter of bead.
          pb = max(0, self.Kb*db)             # bead backpressure.
          pn = max(0, pe - pb)                # nozzle pressure.
          nde = pn*dt/(self.Kf+dt)            # Get filament extruded out the nozzle.
          lde = min(db*h*dl/self.Fa, eb+nde)  # Get volume removed from bead to line,
          pe -= nde
          eb = max(0.0, eb + nde - lde)
          el += lde
          t -= dt
      assert isneareq(v, m.v1), f'{m!r} {v=} != {m.v1=}'
      assert isneareq(l, m.dl), f'{m!r} {l=} != {m.dl=}'
    else:
      # With Kf and Kb zero, there is no nozzle or bead pressure, so all
      # filament gets extruded and pe<=0 (can be negative when retracted). The
      # bead size is max the target move bead size, but can be less if
      # insuficient filament is extruded. The filament added to the line is
      # all the filament extruded minus the bead.
      ex = self.pe + m.de
      pe = min(0.0, ex)
      en = max(0.0, ex)
      eb = min(m.eb, self.eb + en)
      el = en + self.eb - eb
    #print(f'{self.Kf=} {self.Kb=} {self.Re=} {self.pe=} {self.eb=} {m} -> {pe=} {eb=} {el=}')
    assert isneareq(pe + eb + el, self.pe + self.eb + m.de), f'{pe=}+{eb=}+{el=} != {self.pe=}+{self.eb=}+{m.de=}'
    # Note we return pe,eb,en not pe,eb,el, because "what came out the nozzle"
    # is more useful than "what came out the nozzle not including the bead".
    # We also recalculate en from the raw inputs to avoid accumulated integration errors.
    return pe, eb, self.pe + m.de - pe

  def draw(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      v=None, s=1.0, h=None, w=None, r=1.0):
    """ Draw a line.

    Movement can be absolute using `x,y,z,e` arguments, or relative using
    `dx,dy,dz,de` arguments. The z axis can also be specified with `h` height
    above the base of the current layer. The extrusion can also be specified
    using an `r` extrusion ratio muliplier for a `h*w` line, where `h` and `w`
    default to the current layer's if not specified, and r is multiplied by
    the current layer's extrusion ratio.  Speed can be specified directly
    using the `v` velocity argument, or the `s` speed multiplier of the
    default `Vp,Vt,Vz,Ve,Vb` speeds for print, travel, vertical, retract, and
    restore movement.

    Note this is the simple raw movement. Any pressure advance compensation
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
      self.add(m)
    else:
      self.log(f'skipping empty move {m}')

  def move(self, x=None, y=None, z=None, r=0.0, **kwargs):
    """ Move the extruder.

    This is the same as draw() but with the default r=0.0 so it travels instead
    of prints by default.
    """
    self.draw(x, y, z, r=r, **kwargs)

  def retract(self, e=None, de=None, ve=None, s=1.0):
    """Do a retract.

    This is a simple retraction with default self.Re. Note compensation for
    advance pressure is added in post processing if enabled.
    """
    if de is None: de = -self.Re
    # If dynamic retraction is enabled and de is zero, set de to a tiny
    # retraction so that it doesn't get optimized away and can be dynamically
    # adjusted later.
    if self.dynret and not de: de = -0.00001
    self.move(e=e, de=de, v=ve, s=s)

  def restore(self, e=None, de=None, vb=None, s=1.0):
    """ Do a restore.

    This is a simple restore with default self.Re. Note adding starting dots
    and compensating for advance pressure is added in post-processing if
    enabled.
    """
    if de is None: de = self.Re
    # If dynamic retraction is enabled and de is zero, set de to a tiny
    # restore so that it doesn't get optimized away and can be dynamically
    # adjusted later.
    if self.dynret and not de: de = 0.00001
    self.move(e=e, de=de, v=vb, s=s)

  def hopup(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      vt=None, vz=None, ve=None, vb=None,
      s=1.0, h=None):
    """ Do a retract, raise, and move. """
    # default hopup is to Zh above layer height.
    if h is None: h = self.layer.h + self.Zh
    self.retract(e=e, de=de, ve=ve, s=s)
    self.move(z=z, dz=dz, v=vz, s=s, h=h)

  def hopdn(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      vt=None, vz=None, ve=None, vb=None,
      s=1.0, h=None):
    """ Do a move, drop, and restore. """
    # default hop down is to layer height.
    if h is None: h = self.layer.h
    if (x, y, dx, dy) != (None, None, None, None):
      self.move(x=x, y=y, dx=dx, dy=dy, v=vt, s=s)
    self.move(z=z, dz=dz, v=vz, h=h, s=s)
    self.restore(e=e, de=de, vb=vb, s=s)

  def wait(self, t):
    self.add(round(t*1000))

  def _getNextPeEb(self, gcode):
    "Get pe and eb for the next draws."
    dt = dl = vedl = dbdl = eb = 0.0
    # Get the line-average ve and db for the next 10mm and 1s of moves.
    for m1 in (m for m in gcode if isinstance(m, Move)):
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
      # Calculate pe using the line-average db and ve.
      return self._calc_pe(dbdl/dl, vedl/dl), eb
    return 0.0, 0.0

  def getCode(self):
    if self.optmov:
      # Join together any moves that can be joined.
      ljoin(self.gcode)
    # Fix all the velocities.
    lfixv(self.gcode)
    if self.fixvel:
      # set velocity to vm for decelerating moves.
      self.gcode = [m.change(v=m.vm) if isinstance(m,Move) and m.v>m.v0==m.vm>m.v1 else m for m in self.gcode]
    if self.dynext:
      # Partition all the moves into phases
      self.gcode = lsplit(self.gcode)
    # finalize and return the resulting gcode.
    gcode = self.gcode
    self.resetfile()
    for i,c in enumerate(gcode):
      # Adjust Move's to include dynamic retraction and extrusion.
      if isinstance(c, Move):
        if c.isretract and self.dynret:
          # Adjust retraction to Re and also relieve advance pressure.
          c = c.change(de=-self.Re - self.pe)
          # If de is zero, skip adding this.
          if not c.de:
            self.flog('skipping empty retract.')
            continue
        elif c.isrestore and self.dynret:
          # Get the pe and eb for the next draws.
          pe, eb = self._getNextPeEb(gcode[i+1:])
          # Adjust existing restore and add the starting bead.
          #self.flog(f'adjusting {c} for {pe=:.4f} {eb=:.4f}
          c = c.change(de=pe - self.pe + eb)
          # If de is zero, skip adding this.
          if not c.de:
            self.flog('skipping empty restore.')
            continue
        elif c.isdraw and self.dynext:
          # For calculating the pressure pe needed, use the ending ve1.
          pe = self._calc_pe(c.db, c.ve1)
          # Adjust de to include required change in pe over the move.
          #self.flog(f'adjusting {c} {c.ve=:.4f}({c.ve0:.4f}<{c.vem:.4f}>{c.ve1:.4f}) {c.isaccel=}')
          # TODO: This tends to oscillate, change to a PID or PD controller.
          c = c.change(de=c.de + pe - self.pe)
      self.fadd(c)
    return b'\n'.join(self.gcode) + b'\n'

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
    self.cmt('preExtrude:{h:.2f}')
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
          # hit the bottom corner, start following top or bottom.
          lx+=y0-ly
          ly=y0
        elif y1<ly:
          # hit the bottom corner, start following top or bottom.
          lx+=ly-y1
          ly=y1
      else:
        # left diagonal end is following the top.
        lx+=dx
        if lx>x1: lx=x1
      #self.log(f'pos=({{x}},{{y}}) l=({lx},{ly}) r=({rx},{ry})')
      if self.y==ry or self.x==rx :
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
  cmdline.add_argument('-R', action='store_true',
      help='Enable Linear Advance dynamic retract/restore.')
  cmdline.add_argument('-P', action='store_true',
      help='Enable Linear Advance dynamic extrusion rates.')
  cmdline.add_argument('-O', action='store_true',
      help='Enable optimize, merging draw/move cmds where possible.')
  cmdline.add_argument('-V', action='store_true',
      help='Enable fix velocities for acceleration limits.')
  cmdline.add_argument('-d', action='store_true',
      help='Enable diff mode output.')
  cmdline.add_argument('-v', action='store_true',
      help='Enable verbose logging output.')


def GCodeGetArgs(args):
  """ Get the dict of standard GCode arguments. """
  return dict (
      Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc,
      Kf=args.Kf, Kb=args.Kb, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve, Vb=args.Vb,
      Lh=args.Lh, Lw=args.Lw, Lr=args.Lr,
      dynret=args.R, dynext=args.P, optmov=args.O, fixvel=args.V,
      relext=args.E, diff=args.d, verb=args.v)


if __name__ == '__main__':
  import sys

  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=GCodeGen(**GCodeGetArgs(args))

  gen.startFile()
  gen.preExt(-10,-10,10,10)
  gen.startLayer()
  gen.dbox(-10,-10,10,10)
  gen.endLayer()
  gen.startLayer()
  gen.dbox(-10,-10,10,10)
  gen.endLayer()
  gen.startLayer()
  gen.text("Hello\nWorld!",-8,8,8,-8)
  gen.endLayer()
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
