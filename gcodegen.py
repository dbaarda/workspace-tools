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
  """ Move instructions.

  These represent gcode moves enhanced to include the start and end
  velocities. They can partition a basic gcode "intent" move command into
  "physical" acceleration, constant velocity, deceleration components based on
  the printers acceleration, velocity, and cornering specs. They assume the
  printer interprets gcode into movements as described at;

  https://mmone.github.io/klipper/Kinematics.html

  Note for printers with different specs this can be subclassed with the `Ap`,
  `Vp`, and `Gp` overridden.
  """
  dx: float
  dy: float
  dz: float
  de: float
  v: float  # target velocity in mm/s.
  v0: float = None # start velocity in mm/s.
  v1: float = None # end velocity in mm/s.
  vm: float = None # mid velocity in mm/s.
  h: float = None  # The line or end of move height above the current layer.
  w: float = None  # The line width.
  # Note the extrusion ratio r is calculated as a property.

  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Ap = 500  # printer acceleration in mm/s^2.
  Vp = 100  # printer max velocity in mm/s.
  Gp = 0.1  # printer cornering deviation distance in mm.
  Fa = acircle(Fd)  # Fillament area in mm^2.
  Na = acircle(Nd)  # Nozzle area in mm^2.
  am = 0.0 # acceleration for mid phase of movements.
  log = None  # optional log(msg) function to use.

  def __str__(self):
    dx,dy,dz,de,v,v0,v1,vm,h,w = self
    dl, r = self.dl, self.r
    v,v0,v1,vm = round(v),round(v0),round(v1),round(vm)
    if self.isdraw:
      return f'draw {dl=:.2f}@{v}({v0}<{vm}>{v1}) l={h:.2f}x{w:.2f}x{r:.2f}'
    elif self.isup:
      return f'hopup {dz=:.2f}@{v} {h=}'
    elif self.isdn:
      return f'hopdn {dz=:.2f}@{v} {h=}'
    elif self.isretract:
      return f'retract {de=:.4f}@{v}'
    elif self.isrestore:
      return f'restore {de=:.4f}@{v}'
    else:
      return f'move {dl=:.2f}@{v}({v0}<{vm}>{v1}) {de=:.4f} {h=:.2f}'
    #return f'{self.__class__.__name__}({dx=:.2f}, {dy=:.2f}, {dz=:.2f}, {de=:.4f}, {v=:.2f}, {v0=:.2f}, {v1=:.2f}, {h=:.2f}, {w=:.2f})'

  def set(self, *args, s=None, **kwargs):
    """ Change or scale a move. """
    if s is not None:
      # add scaled fields to kwargs if not in args or kwargs.
      fields = self._fields[len(args):4]
      kwargs.update({n:s*getattr(self,n) for n in fields if n not in kwargs})
    return self._replace(*args,**kwargs)

  @classmethod
  def _dvdl(cls, v0, v1, a=None):
    """ The distance travelled accelerating between two velocities. """
    if a is None: a = cls.Ap if v1>=v0 else -cls.Ap
    return (v1**2 - v0**2)/(2*a)

  @classmethod
  def _dvdt(cls, v0, v1, a=None):
    """ The time taken accelerating between two velocities. """
    if a is None: a = cls.Ap if v1>=v0 else - cls.Ap
    return (v1-v0)/a

  @classmethod
  def _eb(cls, h, w, r=1.0):
    """ Volume in filament length of bead."""
    return h * acircle(w*r) / cls.Fa

  @property
  def r(self):
    """ The extrusion ratio of a move. """
    if self.isdraw:
      return (self.de * self.Fa) / (self.h * self.w * self.dl)
    return 0.0

  @property
  def dl(self):
    """ The horizontal length of a move. """
    return sqrt(self.dx*self.dx + self.dy*self.dy)

  @property
  def dl0(self):
    """Horizontal length of acceleration phase."""
    return min(self._dvdl(self.v0, self.vm), self.dl) if self.isgo else 0.0

  @property
  def dl1(self):
    """Horizontal length of deceleration phase."""
    return min(self._dvdl(self.vm, self.v1), self.dl) if self.isgo else 0.0

  @property
  def dlm(self):
    """Horizontal length of middle phase."""
    dlm = self.dl - self.dl0 - self.dl1
    assert 0 <= dlm <= self.dl
    return dlm

  @property
  def dt0(self):
    """Duration of acceleration phase."""
    return 2*self.dl0/(self.v0 + self.vm) if self.isgo else 0.0

  @property
  def dt1(self):
    """Duration of deceleration phase."""
    return 2*self.dl1/(self.vm + self.v1) if self.isgo else 0.0

  @property
  def dtm(self):
    """Duration of middle phase."""
    return abs(self.dlm/self.vm if self.isgo else (self.dz if self.dz else self.de)/self.v)

  @property
  def dt(self):
    """ Get the execution time for a move. """
    return self.dt0 + self.dtm + self.dt1

  @property
  def a0(self):
    """ The acceleration at the start of a move. """
    dt0 = self.dt0
    return (self.vm - self.v0)/dt0 if dt0 else 0.0

  @property
  def a1(self):
    """ The acceleration at the start of a move. """
    dt1 = self.dt1
    return (self.v1 - self.vm)/dt1 if dt1 else 0.0

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
    return self.de/dl * self.v if dl else self.ve

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

  @property
  def eb(self):
    """ The extrusion bead volume. """
    return self._eb(self.h, self.w, self.r)

  @property
  def el(self):
    """ The line volume not including r under/over extrusion. """
    return self.h*self.w*self.dl/self.Fa

  @property
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
  def isup(self):
    return self.dz > 0 and not self.isgo

  @property
  def isdn(self):
    return self.dz < 0 and not self.isgo

  @property
  def isspeed(self):
    return not (self.isgo or self.dz or self.de)

  @property
  def isgo(self):
    return self.dx or self.dy

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
    return min(sqrt(self.Ap*R), self.v, m.v)

  def fixv(self, v=None, v0=0.0, v1=0.0, m0=None, m1=None):
    """ Update move velocities for moves before/after.

    This updates a move's velocities to what the printer will do based on
    on the moves before and after.

    if v is not set it uses m.v.
    if m0 is set the v0 start velocity is set to the m0.v1 end velocity.
    If m1 is set the v1 end velocity is set to the corner velocity.
    If v0, or v1 are still not set they default to 0.0.
    """
    if v is None: v = self.v
    if m0 is not None: v0 = m0.v1
    if m1 is not None: v1 = self.cornerv(m1)
    dl = self.dl
    if not dl:
      # A stop move, v0 and v1 are zero.
      m = self.set(v=v, v0=0, v1=0, vm=0)
    else:
      if self.log and self._dvdl(v0, v1) > dl:
        self.log(f'WARNING:cannot accelerate from {v0}mm/s to {v1}mm/s over {dl}mm with a={self.Ap}mm/s^2.')
      #assert self._dvdl(v0, v1) <= dl, f'Cannot accelerate {v0=}mm/s to {v1=}mm/s over {dl=}mm.'
      # This is the limit velocity assuming constant acceleration at a/2.
      # This implements smoothed look-ahead to reduce spikey velocity.
      vm = sqrt((dl*self.Ap + v0**2 + v1**2)/2)
      # If vlim is less than v0 or v1, we don't need acceleration at both ends.
      vm = max(v0, v1, vm)
      vm = min(v, vm)  # don't set vm higher than requested.
      assert 0 <= v0 <= vm <= self.Vp
      assert 0 <= v1 <= vm <= self.Vp
      m = self.set(v=v,v0=v0,v1=v1,vm=vm)
    return m

  def join(self, m):
    """ Add two moves together. """
    return self._replace(dx=self.dx+m.dx, dy=self.dy+m.dy, dz=self.dz+m.dz,
        de=self.de+m.de, v=min(self.v, m.v),v0=min(self.v0, m.v0), v1=min(m.v1, m.v1), vm=min(self.vm, m.vm))

  def split(self):
    """ Partition a move into a list of moves for the acceleration phases. """
    v0,v1,vm,dl = self.v0, self.v1, self.vm, self.dl
    if not dl or v0 == vm == v1 or self.de <= 0:
      # A non-moving, constant speed, or non-extruding move, don't partition it.
      return [self]
    # Get distances for acceleration phases.
    dl0, dl1 = self.dl0, self.dl1
    # initialize middle phase v0 and v1 to just vm.
    vm0, vm1 = vm, vm
    # if acceleration distances are small don't bother partitioning them.
    if dl0 < 1.0: dl0, vm0 = 0.0, v0
    if dl1 < 1.0: dl1, vm1 = 0.0, v1
    # Get the middle phase distance.
    dlm = dl - dl0 - dl1
    assert 0 <= dlm <= dl
    # Make a list of the s,v0,v1 values for each phase.
    phases = [(dl0/dl, v0, vm), (dlm/dl, vm0, vm1), (dl1/dl, vm, v1)]
    # Return a list moves for the non-zero phases.
    l=[self.set(v0=v0, v1=v1, s=s) for s,v0,v1 in phases if s]
    return l

  def canjoin(self,m):
    # Always join moves together.
    return ((self.ismove and m.ismove) or
      # Join draws together if ...
      (self.isdraw and m.isdraw and
      # ... they have the same line settings, and ...
      (self.h,self.w,self.r) == (self.h,self.w,self.r) and
      # The path deviation from joining them is less than Gp.
      self.joind(m) < self.Gp))


class Move(MoveBase):
  __slots__ = ()

  def __init_subclass__(cls, Fd=1.75, Nd=0.4, Ap=500, Vp=100, Gp=0.1, log=None):
    super().__init_subclass__()
    cls.Fd, cls.Nd, cls.Ap, cls.Vp, cls.Gp, cls.log = Fd, Nd, Ap, Vp, Gp, log
    cls.Fa = acircle(Fd)
    cls.Na = acircle(Nd)

  def __new__(cls, *args, **kwargs):
    """ Make v0 and v1 default to v for moves and 0.0 for non-moves. """
    # Round all arguments to 6 decimal places to clean out fp errors.
    args = tuple(round(a,6) for a in args)
    kwargs = {k:round(v,6) for k,v in kwargs.items()}
    m = super().__new__(cls, *args, **kwargs)
    assert 0 <= m.v <= cls.Vp
    assert m.v0 is None or 0 <= m.v0 <= m.v
    assert m.v1 is None or 0 <= m.v1 <= m.v
    assert m.vm is None or 0 <= m.vm <= m.v
    if m.v0 is None or m.v1 is None or m.vm is None:
      v0 = m.v if m.v0 is None else m.v0
      v1 = m.v if m.v1 is None else m.v1
      vm = m.v if m.vm is None else m.vm
      if m.isgo:
        return m._replace(v0=v0, v1=v1, vm=vm)
      else:
        return m._replace(v0=0.0, v1=0.0, vm=0.0)
    return m


def ljoin(l):
  # TODO: make this and lfixv methods of GCodeGen and make ljoin not optimze
  # pre-extrudes.
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


def lfixv(l, v=None, v0=0.0, v1=0.0):
  """ Fix velocities in a list of moves/comments/commands. """
  i = m = None
  for i1,m1 in enumerate(l):
    if isinstance(m1, Move):
      if m is None:
        # This is the first m to fix.
        i, m = i1, m1
      else:
        # We have the next move m1, fix m in the list.
        l[i] = m = m.fixv(v=v, v0=v0, m1=m1)
        # v0 is the last m.v1, m1 is the next m to fix,
        v0, i, m = m.v1, i1, m1
  if m is not None:
    # Fix the last move.
    l[i] = m.fixv(v=v, v0=v0, v1=v1)
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
    en_relext: relative extrusion enabled.
    en_dynret: linear advance dynamic retraction enabled.
    en_dynext: linear advance dynamic extrusion enabled.
    en_optmov: path optimizing enabled.
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
  """

  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Fa = acircle(Fd)  # Fillament area.
  Na = acircle(Nd)  # Nozzle area.
  Zh = Nd    # z-hop height.

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
      en_relext=False, en_dynret=False, en_dynext=False, en_optmov=False, en_verb=False):
    # Temp and Fan settings.
    self.Te, self.Tp, self.Fe, self.Fc = Te, Tp, Fe, Fc
    # Linear advance and retraction settings.
    self.Kf, self.Kb, self.Re = Kf, Kb, Re
    # Default velocities.
    self.Vp, self.Vt, self.Vz, self.Ve, self.Vb = Vp, Vt, Vz, Ve, Vb
    # Default line height, width, and extrusion ratio.
    self.Lh, self.Lw, self.Lr = Lh, Lw, Lr
    # Processing mode options.
    self.en_relext = en_relext
    self.en_dynret, self.en_dynext, self.en_optmov, self.en_verb = en_dynret, en_dynext, en_optmov, en_verb
    class GMove(Move, Fd=self.Fd, Nd=self.Nd, log=self.log): pass
    self.GMove = GMove
    self.resetfile()

  @property
  def h(self):
    """ Get the current height above the layer base."""
    return self.z - (self.layer.z if self.layer else 0.0)

  def resetfile(self):
    """ Reset all state variables. """
    # Zero all positions so the first move's deltas are actually absolute.
    self.x, self.y, self.z, self.e, self.f = 0.0, 0.0, 0.0, 0.0, 0
    self.f_t = self.f_l = self.f_e = 0.0
    self.pe = 0.0  # current extruder linear advance pressure or retracted length mm.
    self.eb = 0.0  # current extruded bead volume in mm of fillament.
    self.gcode = []
    self.resetlayer()

  def resetlayer(self):
    """ Reset all layer state variables. """
    self.layer = None
    self.l_t = self.l_l = self.l_e = 0.0

  def inclayer(self, l):
    """ Increment state for starting a Layer. """
    self.resetlayer()
    self.layer = l

  def incmove(self, m):
    """ Increment state for executing a Move. """
    dx,dy,dz,de,v,v0,v1,*_ = m
    dt, dl = m.dt, m.dl
    self.x += dx
    self.y += dy
    self.z += dz
    self.e += de
    self.f = round(v*60)
    self.pe, self.eb, en = self._calc_pe_eb_en(m)
    self.f_t += dt
    self.f_l += dl if en else 0.0
    self.f_e += en
    self.l_t += dt
    self.l_l += dl if en else 0.0
    self.l_e += en
    self.vl = v1
    self.ve = de/dl*v1 if dl else de/dt if dt else 0.0
    self.vn = self._vn()
    self.db = self._db()
    self.dn = en

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
    if self.en_relext:
      e = de if de else None
    else:
      e = self.e + de if de else None
    f = round(60 * v)
    if f == self.f: f = None
    cmd = self.fcmd('G1', X=x, Y=y, Z=z, E=e, F=f)
    if self.en_verb:
      cmd = f'{cmd:40s}; {m} r={{dn/{m.el} if {m.dl} else db/{m.w}:.2f}}'
    return cmd

  def fadd(self, code):
    """ Do a final format to string and add. """
    if isinstance(code, Layer):
      self.inclayer(code)
    elif isinstance(code, Move):
      out = self.fmove(code)
      self.incmove(code)
      self.fadd(out)
      self.flog('{pe=:.4f} {eb=:.4f} {ve=:.4f} {vn=:.4f} {vl=:.3f} {db=:.2f}')
    elif isinstance(code, int):
      # ints are waits.
      self.f_t += code/1000
      self.l_t += code/1000
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
      self.f_t += code/1000
      self.l_t += code/1000
    elif isinstance(code, str):
      self.inccmd(code)
    self.gcode.append(code)

  def cmd(self, cmd, **kwargs):
    self.add(self.fcmd(cmd,**kwargs))

  def cmt(self, cmt):
    self.add(f';{cmt}')

  def log(self, txt):
    #print(self.fstr(txt))
    if self.en_verb:
      self.cmt('{ftime(f_t)}: ' + txt)

  def flog(self, txt):
    """Add a formatted log entry."""
    if self.en_verb:
      self.log(txt)
      # We need to fstr() reformat that last log entry...
      self.fadd(self.gcode.pop())

  def startFile(self):
    self.add(self.fstr(self.startcode))
    self.startLayer(0)

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
    self.up(**upargs)
    self.layerstats()

  def _calc_pe(self, db, ve):
    """ Get steady-state pressure advance pe needed for target db and ve. """
    pb = max(0, self.Kb*db)
    pn = max(0, self.Kf*ve)
    #pn = self.Kf*ve * e**(-dt/self.Kf)
    return pn+pb

  def _db(self, eb=None, h=None):
    """ Get bead diameter from bead volume."""
    if eb is None: eb = self.eb
    if h is None: h = self.h
    # h=0 can happen for moves at the start of a new layer before a hop.
    # When this happens, assume h is the current layer height.
    if h < 0.01: h = self.layer.h
    assert eb >= 0.0
    return 2*sqrt(eb*self.Fa/(h*pi)) if eb else 0.0

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
        assert m.v0 == m.v1 == m.vm
        v, vz, ve = 0.0, m.dz/dt if dt else 0.0, m.de/dt if dt else 0.0
        # There is only one phase.
        phases =((dt,0),)
      else:
        v, dz_dl, de_dl = m.v0, m.dz/dl, m.de/dl
        # Calculate accelerate, move, decelerate phases.
        phases = (m.dt0, m.a0),(m.dtm, m.am),(m.dt1, m.a1)
      pe, eb, el, h = self.pe, self.eb, 0.0, self.h
      for t,a in phases:
        while t>0:
          dt = min(0.001, t)  # Use this iteration time interval.
          dv = a*dt
          dl = (v + dv/2)*dt  # Do trapezoidal integration for dl.
          assert dl >= 0.0
          v += dv
          assert v >= -0.001
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
      assert abs(v - m.v1) < 0.001, f'{m} {v=} != {m.v1=}'
    else:
      # With Kf and Kb zero, there is no nozzle or bead pressure, so all
      # filament gets extruded and pe<=0 (can be retracted). The bead size is
      # max the target move bead size, but can be less if insuficient filament
      # is extruded. The filament added to the line is all the filament extruded
      # minus the bead.
      ex = self.eb + self.pe + m.de
      pe, eb, el = min(ex, 0), min(m.eb, max(ex, 0)), max(0, ex - m.eb)
    #print(f'{self.Kf=} {self.Kb=} {self.Re=} {self.pe=} {self.eb=} {m} -> {pe=} {eb=} {el=}')
    assert abs((pe + eb + el)-(self.pe + self.eb + m.de)) < 0.00001, f'{pe=}+{eb=}+{el=} != {self.pe=}+{self.eb=}+{m.de=}'
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
    if w is None: w = self.layer.w
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
    if self.en_dynret and not de: de = -0.00001
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
    if self.en_dynret and not de: de = 0.00001
    self.move(e=e, de=de, v=vb, s=s)

  def up(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      vt=None, vz=None, ve=None, vb=None,
      s=1.0, h=None):
    """ Do a retract, raise, and move. """
    # default hop up is to Zh above layer height.
    if h is None: h = self.layer.h + self.Zh
    self.retract(e=e, de=de, ve=ve, s=s)
    self.move(z=z, dz=dz, v=vz, s=s, h=h)

  def dn(self,
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

  def getCode(self):
    if self.en_optmov:
      # Join together any moves that can be joined.
      ljoin(self.gcode)
    # Fix all the velocities.
    lfixv(self.gcode)
    if self.en_dynext:
      # Partition all the moves into phases
      self.gcode = lsplit(self.gcode)
    # finalize and return the resulting gcode.
    gcode = self.gcode
    self.resetfile()
    for i,c in enumerate(gcode):
      # Adjust Move's to include dynamic retraction and extrusion.
      if isinstance(c, Move):
        if c.isretract and self.en_dynret:
          # Adjust retraction to Re and also relieve advance pressure.
          c = c.set(de=-self.Re - self.pe)
          # If de is zero, skip adding this.
          if not c.de:
            self.flog('skipping empty retract.')
            continue
        elif c.isrestore and self.en_dynret:
          # Find the next move to get the pe and eb it needs.
          m1 = next((m for m in gcode[i+1:] if isinstance(m, Move)), None)
          if m1 and m1.isdraw:
            # For calculating the pressure pe needed, use middle phase ve.
            eb, pe = m1.eb, self._calc_pe(m1.db, m1.vem)
            #self.log(f'adjusting {c} for next {m1}, {m1.ve=:.4f}({m1.ve0:.4f}<{m1.vem:.4f}>{m1.ve1:.4f}) {m1.isaccel=}')
          else:
            eb, pe = 0.0, self._calc_pe(0.0, 0.0)
          # Adjust existing retraction and add the starting bead.
          c = c.set(de=pe - self.pe + eb)
          # If de is zero, skip adding this.
          if not c.de:
            self.flog('skipping empty restore.')
            continue
        elif c.isdraw and self.en_dynext:
          # For calculating the pressure pe needed, use the ending ve1.
          pe = self._calc_pe(c.db, c.ve1)
          # Adjust de to include required change in pe over the move.
          #self.log(f'adjusting {c} {c.ve=:.4f}({c.ve0:.4f}<{c.vem:.4f}>{c.ve1:.4f}) {c.isaccel=}')
          # TODO: This tends to oscillate, change to a PID or PD controller.
          c = c.set(de=c.de + pe - self.pe)
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

  def preExt(self,x0,y0,x1,y1,le=120.0, lw=100.0, m=10.0, ve=20.0, vw=10, r=4):
    """Preextrude around a box with margin m.

    This extrudes for `le` mm at speed `ve` mm/s, then wipes for another `lw`
    mm at speed `vw` mm/s without extruding, and finally another 10mm at 1mm/s
    to stick the last bit of drool to the plate. This should leave the nozzle
    primed with no residual advance pressure, and prevent the wipe from
    becoming a big piece of stringing that messes with the print.
    """
    self.cmt('preExtrude:{layer.h:.2f}')
    self.cmt('structure:pre-extrude')
    x0,x1 = sorted([x0,x1])
    y0,y1 = sorted([y0,y1])
    x0-=m
    x1+=m
    y0-=m
    y1+=m
    x, y = x0, y0
    l, v = le, ve
    self.dn(x=x,y=y)
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
      elif l == 0 and r == 0 and v>1:
        # We have finished the wipe, go another 10mm at 1mm/s to stick the
        # final drool to the plate.
        l, v, r = 10,1,0
    self.up()

  def fbox(self,x0,y0,x1,y1,**kwargs):
    # Fill a box with frame lines.
    w = kwargs.get('w', self.layer.w)
    # make sure x0<x1 and y0<y1
    x0,x1 = sorted((x0,x1))
    y0,y1 = sorted((y0,y1))
    x0+=w/2
    x1-=w/2
    y0+=w/2
    y1-=w/2
    self.dn(x0,y0)
    while x1-x0>2*w and y1-y0>2*w:
      self.draw(x0,y1,**kwargs)
      self.draw(x1,y1,**kwargs)
      self.draw(x1,y0,**kwargs)
      self.draw(x0+w,y0,**kwargs)
      x0+=w
      x1-=w
      y0+=w
      y1-=w
    if y1-y0 > 2*w:
      fw=x1-x0
      self.move(x=(x0+x1)/2,y=y0+(fw-w)/2)
      self.draw(y=y1-(fw-w)/2, w=fw)
    else:
      fw=y1-y0
      self.move(x=x0+(fw-w)/2,y=(y1+y0)/2)
      self.draw(x=x1-(fw-w)/2, w=fw)
    self.up()

  def dbox(self,x0,y0,x1,y1,**kwargs):
    """ fill a box with diagonal lines. """
    w = kwargs.get('w', self.layer.w)
    x0,x1 = sorted((x0,x1))
    y0,y1 = sorted((y0,y1))
    x0+=w/2
    x1-=w/2
    y0+=w/2
    y1-=w/2
    self.dn(x0,y0)
    # initialize left and right diagonal ends at bottom left corner.
    lx,ly=rx,ry=x0,y0
    sqrt2w=w*2**0.5
    while ry<y1:
      if rx<x1:
        # right diagonal end is following bottom edge
        rx+=sqrt2w
        if rx>x1:
          # hit bottom right corner, start going up right side.
          ry+=x-x1
          rx=x1
      else:
        # right diagonal end is following right side.
        ry+=sqrt2w
        if ry>y1: ry=y1
      if ly<y1:
        # left diagonal end is following left side.
        ly+=sqrt2w
        if ly>y1:
          # hit top left corner, start following top.
          lx+=y-y1
          ly=y1
      else:
        # left diagonal end is following the top.
        lx+=sqrt2w
        if lx:x1: lx=x1
      if self.y == y0 or self.x == x1:
        # starting diagonal from right point.
        self.draw(x=rx,y=ry,**kwargs)
        self.draw(x=lx,y=ly,**kwargs)
      else:
        # starting diagonal from left point.
        self.draw(x=lx,y=ly,**kwargs)
        self.draw(x=rx,y=ry,**kwargs)

  def dot(self, x, y, r=1.0, **kwargs):
    # Note this doesn't seem to render anything in FlashPrint.
    self.dn(x,y)
    # Draw a tiny line so FlashPrint renders something.
    self.draw(dy=-0.2,r=r, **kwargs)
    self.up()

  def line(self, l, r=1.0, **kwargs):
    """Draw a line from a sequence of points."""
    x0,y0 = l[0]
    # If there is only one point, draw a dot instead.
    if len(l) == 1 or (len(l) == 2 and l[0] == l[1]):
      return self.dot(x0, y0, r, **kwargs)
    self.dn(x0,y0)
    for x,y in l[1:]:
      self.draw(x,y,**kwargs)
    self.up()

  def text(self, t, x0, y0, x1=None, y1=None, fsize=5, angle=0, **kwargs):
    w = kwargs.get('w', self.layer.w)
    v = vtext.ptext(t,x0,y0,x1,y1,fsize,angle,w)
    self.log(f'text {x0=} {y0=} {fsize=} {angle=} {t!r}')
    for l in v:
      self.line(l,**kwargs)


class ExtrudeTest(GCodeGen):
  """

  Note that nozzle vs fillament diameters and/or track areas mean that de=1mm
  translates to nearly 20mm of nozzle thread, or possibly even more of track.
  This means 5mm of uncompensated pressure advance translates to more 100mm of
  track smear or stringing.

  Linear Advance or Pressure Advance is mm of extrusion advance per mm/s of
  filament extrusion speed, and has typical values in the range 0.0-2.0. Note
  that 100mm/s print speed for a h=0.3 layer with w=0.4 width using 1.75mm
  diameter filament is 12mm^3/s or 5mm/s fillament rate. This suggests the
  advance could get as high as 10mm, or 200mm of track!

  For Linear Advance we assume flow rate out the nozzle is linear with
  pressure, and pressure is linear with 'R' extrusion advance of filament
  (filament compression). This gives the equation for how R updates every dt
  time is ;

  R' = ve*dt + R - dt/Kf*R

  Note that this means R reaches a steady state when R = Kf*ve. It is also an
  exponential decay equation, with R exponentially approaching the steady
  state value for a fixed ve rate with a timeconstant of Kf. This exponential
  decay means extrusion effectively lags by Kf seconds! So for Kf=2.0 it will
  take about 2 seconds of extruding before the extrusion rate "catches up" and
  is extruding at about the right rate. At 100mm/s, 2 seconds is 200mm worth
  of line. Note that this matches the 200mm of track ooze you will get when
  stopping.

  """
  vx0 = 5  # test min speed
  vx1 = 100  # test max speed
  le0 = 0  # test min retraction length
  le1 = 10  # test max retraction length
  ve0 = 15  # test min retraction speed
  ve1 = 45  # test max retraction speed
  ty = 2    # height for each test run.
  tn = 10   # number of increments for each test (tn+1 test lines).
  t0x = 5  # test line warmup distance
  t1x = 20  # test line phase distance.
  rdx = 2*t0x + 4*t1x # test ruler length.
  rdy = 4  # ruler height including space
  sdy = 30 # test text settings box width.
  tdx = rdx + sdy  # test box total width.
  tdy = rdy + tn*ty + 5 # test box total height.

  def _fval(self,k,v):
    """ Format a single test argument value or range."""
    v=f'{v[0]}-{v[1]}' if isinstance(v, tuple) else f'{round(v,4)}'
    return f'{k}={v}'

  def _fset(self, sep=' ', **kwargs):
    """ Format several test argument values or ranges. """
    return sep.join(self._fval(k,v) for k,v in kwargs.items())

  def _getstep(self,tn,vr):
    """ Convert a value-range into start, stop, delta. """
    if not isinstance(vr, tuple):
      return vr, vr, 0
    v0, v1 = vr
    dv = (v1-v0)/tn
    return v0, v1, dv

  def pushconf(self, **kwargs):
    """Push and change config attributes on self."""
    if not hasattr(self, 'confstack'):
      self.confstack = []
    # Filter out configs set to None.
    kwargs = {k:v for k,v in kwargs.items() if v is not None}
    old_conf = {k:getattr(self,k) for k in kwargs}
    self.confstack.append(old_conf)
    if kwargs:
      self.log(f'setting config {self._fset(**kwargs)}')
    self.add(kwargs)

  def popconf(self):
    """Pop and restore config attributes on self."""
    conf = self.confstack.pop()
    if conf:
      self.log(f'restoring config {self._fset(**conf)}')
    self.add(conf)

  def incconf(self, conf):
    """Apply a config change to the current state."""
    for k, v in conf.items():
      setattr(self, k, v)

  def fadd(self, code):
    # Test config args don't emit anything in the final output.
    if isinstance(code, dict):
      self.incconf(code)
    else:
      super().fadd(code)

  def add(self, code):
    # Test config args need to be set on self.
    if isinstance(code, dict):
      self.incconf(code)
    super().add(code)

  def preextrude(self,n=0,x0=-50,y0=-60,x1=70,y1=60):
    self.preExt(x0-2*n,y0,x1,y1+2*n,m=5)

  def brim(self, x0, y0, x1, y1):
    w = self.layer.w
    self.cmt("structure:brim")
    self.dn(x0,y0+5)
    self.draw(y=y1-5)
    self.draw(dx=-w)
    self.draw(y=y0+5)
    self.up()
    self.dn(x1,y0+5)
    self.draw(y=y1-5)
    self.draw(dx=w)
    self.draw(y=y0+5)
    self.up()

  def ruler(self, x0, y0, rl=rdx, rs=t0x):
    """ Draw a ruler rl long and 3mm high with tics offset rs from each end."""
    ml,tl=3,1.5  # tick marker lengths
    self.cmt('structure:shell-outer')
    # draw ruler
    self.line([(x0, y0),(x0+rl,y0)])
    self.move(x=x0+rs-2, y=y0)
    for d in range(0,rl - 2*rs + 1,2):
      l = tl if d % 10 else ml
      self.dn(dx=2,y=y0)
      self.draw(dy=-l)
      self.up()

  def settings(self, x0, y0, **kwargs):
    """ Draw the test arguments at (x0,y0)."""
    self.cmt("structure:brim")
    self.text(self._fset(sep='\n', **kwargs), x0, y0, fsize=5)

  def doTests(self, name, linefn, tests, n=None):
    """Run up to 4 tests.

    Args:
      name: the text name of the tests
      linefn: The function that executes a single line of a test.
      tests: A sequence of up to 4 dicts containing args for each test.
      n: optional index of a single test to run, otherwise run all.
    """
    assert len(tests) <= 4, 'Can only fit max 4 tests.'
    # set alltests for running all tests and initialize n=0 if needed.
    alltests, n = n is None, n or 0
    t4x, t4y = self.tdx, 4*self.tdy - 5
    x0, y0 = -self.rdx/2, 55
    x1, y1 = x0+t4x, y0-t4y
    self.preextrude(n=n)
    self.startLayer(Vp=10)
    # Only draw the brim if running test 0.
    if n==0:
      self.brim(x0,y0,x0+self.rdx, y1)
    for i,t in enumerate(tests):
      # only run test n if not running all tests.
      if alltests or i == n:
        self.doTest(name, x0, y0, linefn, **t)
      y0 -= self.tdy
    self.endLayer()

  def doTest(self, name, x0, y0, linefn, **kwargs):
    testargs = {k:self._getstep(self.tn, v) for k,v in kwargs.items()}
    assert any(dv for _,_,dv in testargs.values()), 'At least one arg in {kwargs} must be a range.'
    self.log(f'Test {name} {self._fset(**kwargs)}')
    self.ruler(x0, y0)
    y = y0 - self.rdy
    self.cmt(f'structure:infill-solid')
    e=0.00001
    while all(v <= v1+e for (v,v1,dv) in testargs.values()):
      lineargs={k:v for k,(v,_,_) in testargs.items()}
      confargs={k:v for k,v in lineargs.items() if k in ('Kf','Kb','Re')}
      funcargs={k:v for k,v in lineargs.items() if k not in confargs}
      self.log(f'Test Line {name} {self._fset(**lineargs)}')
      self.pushconf(**confargs)
      linefn(x0, y, **funcargs)
      self.popconf()
      for k,(v,v1,dv) in testargs.items():
        testargs[k] = (v+dv,v1,dv)
      y-=self.ty
    self.settings(x0 + self.rdx+1, y0, **kwargs)

  def testRetract(self, x0, y0, vx=None, ve=None, h=None, w=None, r=1.0, vr=10, re=4):
    """ Test retract and restore for different settings.

    This test is to try and figure out how much retraction is needed to stop
    extruding after moving or figure out how retract/restore speed matters.
    The pre/post phases are;

      * 5mm: warmup draw at 1mm/s to prime the nozzle.
      * 10mm: move at 1mm/s to drain actual and estimated pressure.
      * 0mm: restore default Re to prepare for draw.
      * 30mm: draw at vx mm/s to build pressure at that speed.
      * 40mm: move at vr mm/s while retracting <re>mm.
      * 0mm: restore re mm.
      * 5mm: cooldown draw at 5mm/s

    The the length `l` of the drool line after the draw indicates the amount
    of retraction required to stop extruding of roughly `l*(0.025+re/40)`,
    assuming the drool is average half the thickness of of a normal w=0.4 x
    h=0.3 line which is 0.5mm of filament per 10mm of line.

    Args:
      vx: draw speed before/after retract/restore.
      ve: optional extrude speed while drawing.
      h: optional line height.
      w: optional line width.
      vr: move speed while retracting/restoring.
      re: retract/restore distance.
    """
    assert ve is not None or vx is not None, 'must specify at least ve or vx'
    if vx is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out vx from ve, w, and h.
      vx = ve*self.Fa/(h*w*r)
    if ve is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out ve from vx, w, and h.
      ve = h*w*r*vx/self.Fa
    # Set h and w
    if h is None: h = self.layer.h if w is None else ve*self.Fa/(w*r*vx)
    if w is None: w = ve*self.Fa/(h*r*vx)
    assert abs(ve*self.Fa - h*w*r*vx) < 0.0001
    # A printer with a=500mm/s^2 can go from 0 to 100mm/s in 0.2s over 10mm.
    # Retracting de=10mm at vr=100mm/s over 40mm means retracting at ve=25mm/s
    # which is probably fine.
    self.dn(x0, y0, h=h)
    # This slow draw then move is to ensure the nozzle is primed and wiped, and
    # any advance pressure, actual and estimated, has decayed away.
    self.draw(dx=5,v=1,w=w)
    self.move(dx=10,v=1)
    # This restore pre-loads advance pressure for drawing the line. If the
    # start of the line is too thin, there was not enough pressure, and Kf
    # probably needs to be increased. If the start is too thick, pressure was
    # too high and Kf should probably be dropped. The line should be
    # consistently the same width as all the other lines.
    self.restore()
    self.draw(dx=30,v=vx,w=w,r=r)      # draw at v=vx.
    self.move(dx=40,v=vr,de=-re)          # retract while moving at v=vr.
    self.restore(de=re)                   # restore to recover from retract.
    self.draw(dx=self.t0x,v=self.vx0)     # cooldown draw.
    self.up()

  def testStartStop(self, x0, y0, vx):
    """ Test starting and stopping extrusion for different settings.

    This test is to check dynamic retracting for different Kf/Kb/Re
    settings.

    Test phases are;

      * 0mm: default drop and restore.
      * 5mm: warmup draw at 5mm/s to prime nozzle with fillament.
      * 10mm: move at 1mm/s to drain actual and estimated pressure.
      * 0mm: default restore to pre-apply pressure before draw.
      * 50mm: draw at <vx>mm/s to check pre-applied pressure and build pressure.
      * 0mm: default retract to relieve pressure.
      * 10mm: move at 1mm/s to check and drain any vestigial pressure.
      * 0mm: default restore.
      * 15mm: draw at 5mm/s to check slow draw after retract/restore.
      * 0mm: default retract and raise.

    Args:
      vx: movement speed to check against.
    """
    # Set different line-phase lengths.
    dx = (5, 10, 50, 10, 15)
    self.dn(x0, y0)
    # This slow draw then move is to ensure the nozzle is primed and wiped, and
    # any advance pressure, actual and estimated, has decayed away.
    self.draw(dx=dx[0],v=5)
    self.move(dx=dx[1],v=1)
    # This restore pre-loads advance pressure for drawing the line. If the
    # start of the line is too thin, there was not enough pressure, and Kf
    # probably needs to be increased. If the start is too thick, pressure was
    # too high and Kf should probably be dropped. The line should be
    # consistently the same width as all the other lines.
    self.restore()
    self.draw(dx=dx[2],v=vx)
    # This retract and slow move is to relieve the advance pressure, and see
    # if any remaining pressure drools off. If there is any trailing drool,
    # the amount of pressure was underestimated and Kf or Kb probably needs
    # to be increased. If there is still stringing when otherwise Kf/Kb seems
    # right, adding some additional retraction with Re could help.
    self.retract()
    self.move(dx=dx[3],v=1)
    # This restore applies advance pressure for a final slow line after
    # the earlier retraction. It should be the same width as all the other
    # lines. If it starts too thin, it suggests the earlier retraction
    # over-estimated the pressure and over-retracted, so Kf should be reduced.
    # If it starts too thick it suggests the earlier retraction under
    # estimated the pressure and extra retraction was actually relieving
    # pressure, so Kf should be increased and Re could possibly be reduced.
    self.restore()
    self.draw(dx=dx[4],v=5)
    self.up()

  def testBacklash(self, x0, y0, vx, re):
    """ Test retract and restore distances.

    This test is to test retraction and restore distances for different vx
    velocities.

    Test phases are;

      * 0mm: default drop and restore to pre-apply pressure before draw.
      * 45mm: draw at <vx>mm/s to stabilize pressure at that speed.
      * 20mm: moving retract of <re>mm at 5mm/s to measure required retraction.
      * 20mm: moving restore of <re>mm at 5mm/s to measure required restore.
      * 5mm: draw at 5mm to finalize line and stabilize pressure.
      * 0mm: default retract and raise to relieve any vestigial pressure.

    Args:
      vx: movement speed to check against.
      re: moving retract/restore distance.
    """
    self.dn(x0, y0)
    self.draw(dx=45,v=vx)
    self.move(dx=20,de=-re, v=5)
    self.move(dx=20,de=+re, v=5)
    self.draw(dx=5,v=5)
    self.up()

  def testKf(self, x0, y0, vx0, vx1):
    """ Test linear advance changing speed different speeds.

    This is to test that dynamic extrusion for different Kf/Kb/Re settings for
    different speed changes is right.

    The test phases are;

      * 5mm: warmup draw at 5mm/s
      * 20mm: draw at <vx0>mm/s base speed before changing speed.
      * 40mm: draw at <vx1>mm/s to check changing to and from this speed.
      * 20mm: draw at <vx>mm/s base speed after changing speed.
      * 5mm: cooldown draw at 5mm/s

    Args:
      vx0: the slow speed to use.
      vx1: the fast speed to use.
    """
    self.log(f'testKfline vx0={vx0} vx1={vx1}')
    self.dn(x0, y0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.draw(dx=self.t1x,v=vx0)
    self.draw(dx=2*self.t1x,v=vx1)
    self.draw(dx=self.t1x,v=vx0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.up()


class TextTest(GCodeGen):

  def texttest(self,h=5,x0=-70,y0=70,x1=70,y1=None):
    chars=''.join(c for c in vtext.glyphs)
    self.text(chars,x0,y0,x1,y1,fsize=h)


inf=float('inf')
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
  cmdline.add_argument('-Fe', type=RangeType(0.0,1.0), default=0.0,
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
  cmdline.add_argument('-v', action='store_true',
      help='Enable verbose logging output.')


if __name__ == '__main__':
  import sys

  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc,
      Kf=args.Kf, Kb=args.Kb, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve, Vb=args.Vb,
      Lh=args.Lh, Lw=args.Lw, Lr=args.Lr,
      en_relext=args.E,
      en_dynret=args.R, en_dynext=args.P, en_optmov=args.O, en_verb=args.v)

  backpressureargs=dict(name="Backpressure", linefn=gen.testRetract, tests=(
    dict(w=(0.3, 0.8), vx=10, ve=0.3*0.3*10/Move.Fa),
    dict(h=(0.1,0.35), ve=0.5*0.35*10/Move.Fa, w=0.5),
    dict(w=(0.3, 0.8), ve=0.2*0.8*10/Move.Fa, h=0.2),
    dict(vx=(10, 60), h=0.2, w=0.5)))

  retractargs=dict(name="Retract", linefn=gen.testRetract, tests=(
    dict(vx=(20,100),vr=(20,100),re=8.0),
    dict(vx=(20,100),vr=10,re=8.0),
    dict(vx=(1,10),vr=1,re=4.0)))

  startstopargs=dict(name="StartStop", linefn=gen.testStartStop, tests=(
    dict(Kf=(0.0,1.0), Kb=4.0, Re=1.0, vx=60),
    dict(Kf=0.4, Kb=(2.0,6.0), Re=1.0, vx=60),
    dict(Kf=0.4, Kb=4.0, Re=(0.0,2.0), vx=60),
    dict(Kf=0.4, Kb=4.0, Re=1.0, vx=(20,100))))

  backlashargs=dict(name="Backlash", linefn=gen.testBacklash, tests=(
    dict(Kf=0.5, Kb=2.0, Re=1.0, vx=(20,100), re=5.0),
    ))

  gen.startFile()
  gen.doTests(n=args.n, **startstopargs)
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
