#!/bin/python3
"""
gcodegen.py gcode generator.

This attempts to define a pressure to flow rate model that is better than the
Linear Advance model at accounting for all the non-linear behaviour people
have encountered. It is not intended to be perfect, just better enough to
justify the additional complexity.

The Linear Advance model is...

https://klipper.discourse.group/t/modification-of-pressure-advance-for-high-speed-bowden-printers/13053?u=dbaarda

The new model

The nozzle extrudes into a bead of cooling filament that is pressed against
the print surface. This puts back pressure against the nozzle. The back
pressure will be a function of be bead size; the more pressure, the bigger the
bead.

At zero-to-low speeds the pressure does not define the flow rate, just the
bead size. What does define the flow rate is the nozzle velocity. As the
nozzle moves it smears the bead, leaving a bead-diameter-wide track of
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
  Pb = Kb*Db + Cb   # (3)
  vn = v*Db*h/Af    # (4)
  Db = r*w          # (5)
  Af = pi*(Df/2)^2  # (6)

Where:

  Pe is the extruder pressure in mm of filament advance.
  Pn is the nozzle flow pressure in mm of filament advance.
  Pb is the bead back-pressure in mm of filament advance.
  Kf is the nozzle pressure factor.
  Kb is the bead backpressure factor.
  Cb is the bead backpressure constant.
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
nozle flow velocity `nv` this can be simplified to;

  Pe = Kf*ve + Cf
  Cf = Kb*Db + Cb

Where

  Cf is the pressure linear-advance constant.

We can calculate outflow vn as a function of extruder pressure and nozzle
velocity like this;

  Db = vn/v*Af/h                  # (7) from re-arranging (5)
  Pe = Kf*vn + Kb*Db + Cb         # from (1) substituting in (2), (3), and (4)
     = Kf*vn + Kb*vn/v*Af/h + Cb   # substituting in (7)
     = (Kf + Kb/v*Af/h)*vn + Cb
  (Kf + Kb/v*Af/h)*vn = Pe - Cb
     vn = (Pe - Cb) /(Kf + Kb/v*Af/h)
     vn = (Pe - Cb) * v/(Kf*v + Kb*Af/h) # (8)

This means we only get nozzle outflow when `Pe > Cb` and `v > 0`. As nozzle
velocity v approaches infinity the outflow approaches `(Pe-Cb)/Kf`. The outflow
is half that max rate when `v = Kb/Kf*Af/h`.

Note the bead diameter as a function of extruder pressure Pe and nozzle velocity v is;

  Db = vn/v*Af/h                                # from (4) substituting (7)
     = ((Pe - Cb) * v/(Kf*v + Kb*Af/h))/v * Af/h  # substituting in (8)
     = (Pe - Cb)/(Kf*v + Kb*Af/h) * Af/h
     = (Pe - Cb)/(Kf/Af*h*v + Kb))

This means at nozzle velocity v=0 the bead diameter is `(Pe-Cb)/Kb`, and as
the nozzle velocity approaches infinity it tends towards zero. The bead
diameter is half the v=0 max when `v = Kb/Kf*Af/h`.

The rate of change in extruder pressure is;

dPe/dt = ve - vn
       = de/dt - vn
   dPe = de - dt * vn
       = de - dt * (Pe - Cb) * v/(Kf*v + Kb*Af/h)
       = de - (Pe - Cb) * dt/(Kf + Kb*Af/(v*h))

This is an exponentially decaying process with time-constant `T = = Kf +
Kb*Af/(v*h)`. This means with no additional extrusion ve=0,
extruder pressure Pe will decay towards Cb over a time roughly of `Kf +
Kb*Af/(h*v)`. Note that for v=0 that time is infinity, which means it will not
decay at all.

However, that assumes there is already a bead with backflow pressure Pb to
match Pe. When you first start extruding after a restore, there is no bead and
no backpressure, so there will be outflow till the bead is formed.

If there is no bead yet, there is no backpressure, and there will be outflow
to form a bead even if nozzle velocity v is zero.

Some Adv3 numbers from earlier testing with h=0.3 w=0.6 (db=0.6,eb=0.035mm);

de=+2mm seems to be the minimum for a normal starting dot.
de=+1mm barely seems to register a dot at all.
de=+0.5mm doesn't seem to show any difference from minimum drool.

Suggests roughly Cb=0.8, Kb=2.0

de=+2.3mm? for vl=10mm/s, ve=0.75mm/s
de=+3mm for vl=30mm/s, ve=2.25mm/s
de=+4mm for vl=60mm/s, ve=4.5mm/s
de=+5mm seems about right advance for vl=100mm/s ve=7.5mm/s

Suggests roughly Kf=0.4 with Cf=Pb=Kb*0.6+Cb=2.0.

"""
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
  h: float | None = None # The line height.
  w: float | None = None # The line width.
  # Note the extrusion ratio r is calculated as a property.

  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Ap = 500  # printer acceleration in mm/s^2.
  Vp = 100  # printer max velocity in mm/s.
  Gp = 0.1  # printer cornering deviation distance in mm.
  Fa = acircle(Fd)  # Fillament area in mm^2.
  Na = acircle(Nd)  # Nozzle area in mm^2.
  a0 = Ap  # acceleration for first phase of movements.
  am = 0.0 # acceleration for mid phase of movements.
  a1 = -Ap # acceleration for last phase of movememnts.
  log = None  # optional log(msg) function to use.

  def __str__(self):
    dx,dy,dz,de,v,v0,v1,h,w = self
    dl, r = self.dl, self.r
    v,v0,v1 = round(v),round(v0),round(v1)
    if self.isdraw:
      return f'draw {dl=:.2f}@{v0}<{v}>{v1} l={h:.2f}x{w:.2f}x{r:.2f}'
    elif self.isup:
      return f'hopup {dz=:.2f}@{v} {h=}'
    elif self.isdn:
      return f'hopdn {dz=:.2f}@{v} {h=}'
    elif self.isretract:
      return f'retract {de=:.4f}@{v}'
    elif self.isrestore:
      return f'restore {de=:.4f}@{v}'
    else:
      return f'move {dl=:.2f}@{v0}<{v}>{v1} {h=}'
    #return f'{self.__class__.__name__}({dx=:.2f}, {dy=:.2f}, {dz=:.2f}, {de=:.4f}, {v=:.2f}, {v0=:.2f}, {v1=:.2f}, {h=:.2f}, {w=:.2f})'

  def set(self, *args, s=None, **kwargs):
    """ Change or scale a move. """
    if s is not None:
      # add scaled fields to kwargs if not in args or kwargs.
      fields = self._fields[len(args):4]
      kwargs.update({n:s*getattr(self,n) for n in fields if n not in kwargs})
    return self._replace(*args,**kwargs)

  @classmethod
  def _dvdl(cls, v0, v1):
    """ The distance travelled accelerating between two velocities. """
    return abs(v0**2 - v1**2)/(2*cls.Ap)

  @classmethod
  def _dvdt(cls, v0, v1):
    """ The time taken accelerating between two velocities. """
    return abs(v1-v0)/cls.Ap

  @classmethod
  def _eb(cls, h, w, r=1.0):
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
    return min(self._dvdl(self.v0, self.v), self.dl) if self.isgo else 0.0

  @property
  def dl1(self):
    return min(self._dvdl(self.v, self.v1), self.dl) if self.isgo else 0.0

  @property
  def dlm(self):
    dlm = self.dl - self.dl0 - self.dl1
    assert 0 <= dlm <= self.dl
    return dlm

  @property
  def dt0(self):
    return self._dvdt(self.v0, self.v) if self.isgo else 0.0

  @property
  def dt1(self):
     return self._dvdt(self.v,self.v1) if self.isgo else 0.0

  @property
  def dtm(self):
    dtm =  abs(self.dlm if self.isgo else self.dz if self.dz else self.de)/self.v
    assert 0 <= dtm
    return dtm

  @property
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
  def isdraw(self):
    return self.isgo and self.de

  @property
  def ismove(self):
    return self.isgo and not self.de

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
    """ Get cos(a) of the angle a between two moves. """
    # cos(th) = v0.v1/(|v0|*|v1|)
    dx0, dy0, dl0 = self.dx, self.dy, self.dl
    dx1, dy1, dl1 = m.dx, m.dy, m.dl
    # if a move has no horizontal movement, treat as 0deg.
    return (dx0*dx1 + dy0*dy1)/(dl0 * dl1) if dl0 and dl1 else 1.0

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

    if m0 is set the v0 start velocity is set to the m0.v1 end velocity.
    If m1 is set the v1 end velocity is set to the corner velocity.
    if v is not set it uses m.v.
    If v0, or v1 are not otherwise set they default to 0.0.
    """
    if v is None: v = self.v
    if m0 is not None: v0 = m0.v1
    if m1 is not None: v1 = self.cornerv(m1)
    dl = self.dl
    if not dl:
      # A stop move, v0 and v1 are zero.
      m =self.set(v=v, v0=0, v1=0)
    else:
      if self.log and self._dvdl(v0, v1) > dl:
        self.log(f'WARNING:cannot accelerate {v0}mm/s to {v1}mm/s over {dl}mm with a={self.Ap}mm/s^2.')
      #assert self._dvdl(v0, v1) <= dl, f'Cannot accelerate {v0=}mm/s to {v1=}mm/s over {dl=}mm.'
      # This is the limit velocity assuming constant acceleration at a/2.
      # This implements smoothed look-ahead to reduce spikey velocity.
      vlim = sqrt((dl*self.Ap + v0**2 + v1**2)/2)
      # If vlim is less than v0 or v1, we don't need acceleration at both ends.
      vlim = max(v0, v1, vlim)
      v = min(v, vlim)  # don't set v higher than requested.
      assert v0 <= v <= self.Vp
      assert v1 <= v <= self.Vp
      m = self.set(v=v,v0=v0,v1=v1)
    return m

  def join(self, m):
    """ Add two moves together. """
    v, v0, v1 = (self.v+m.v)/2, self.v0, m.v1
    return self._replace(dx=self.dx+m.dx, dy=self.dy+m.dy, dz=self.dz+m.dz,
        de=self.de+m.de, v=v, v0=v0, v1=v1, h=m.h, w=m.w)

  def split(self):
    """ Partition a move into a list of moves for the acceleration phases. """
    v,v0,v1,dl = self.v, self.v0, self.v1, self.dl
    if not dl or v0 == v == v1 or self.de <= 0:
      # A non-moving, constant speed, or non-extruding move, don't partition it.
      return [self]
    # Get distances for acceleration phases.
    dl0, dl1 = self.dl0, self.dl1
    # initialize middle phase v0 and v1 to just v.
    vm0, vm1 = v, v
    # if acceleration distances are small don't bother partitioning them.
    if dl0 < 1.0: dl0, vm0 = 0.0, v0
    if dl1 < 1.0: dl1, vm1 = 0.0, v1
    # Get the middle phase distance.
    dlm = dl - dl0 - dl1
    # Make a list of the s,v0,v1 values for each phase.
    phases = [(dl0/dl, v0, v), (dlm/dl, vm0, vm1), (dl1/dl, v, v1)]
    # Return a list moves for the non-zero phases.
    l=[self.set(v=v, v0=v0, v1=v1, s=s) for s,v0,v1 in phases if s]
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
    cls.a0, cls.a1 = cls.Ap, -cls.Ap

  def __new__(cls, *args, **kwargs):
    """ Make v0 and v1 default to v for moves and 0.0 for non-moves. """
    m = super().__new__(cls, *args, **kwargs)
    if m.v0 is None or m.v1 is None:
      v0 = m.v if m.v0 is None else m.v0
      v1 = m.v if m.v1 is None else m.v1
      if m.isgo:
        return m._replace(v0=v0, v1=v1)
      else:
        return m._replace(v0=0.0, v1=0.0)
    return m


def ljoin(l):
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
  h : float
  w : float
  r : float


class GCodeGen(object):
  """ GCodeGen gcode generator.

  Attributes:
    Te: Extruder temp (degC).
    Tp: Platform temp (degC).
    Fe: Extruder fan speed (0.0->1.0).
    Fc: Case fan speed (0.0->1.0).
    Kf: Linear Advance factor (0.0->4.0 mm/mm/s).
    Kb: Bead backpressure factor (??? mm/mm).
    Cb: Bead backpressure constant (??? mm).
    Re: Retraction distance (0.0->10.0 mm).
    Vp: Speed when printing (mm/s).
    Vt: Speed when traveling (mm/s).
    Vz: Speed when raising/lowering (mm/s).
    Ve: Speed when retracting (mm/s).
    Vb: Speed when restoring (mm/s).
    r: Extrusion ratio (0.1 -> 10.0).
    h: Line height (0.1->0.4mm).
    w: Line width (0.2->0.8mm).
    en_dynret: linear advance dynamic retraction enabled.
    en_dynext: linear advance dynamic extrusion enabled.
    en_optmov: path optimizing enabled.
    pe: current advance/retract length for extruder pressure.
    eb: current extruded bead volume for nozzle backpressure.
    t: current total execution time.
    l: current total printed line length.
    l_t: current layer execution time.
    l_l: current layer printed line length.
    l_e: current layer extruded filament length.
    vl: move/draw velocity at the end of the last move.
    ve: extruder velocity at the end of the last move.
    vn: nozzle flow velocity at the end of the last move.
    db: extruded bead diameter at the end of the last move.
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
;layer_height: {h:.2f}
;line_width: {w:.2f}
;extrusion_ratio: {r:.2f}
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
      Te=210, Tp=50, Fe=1.0, Fc=1.0, Kf=0.0, Kb=0.0, Cb=0.0, Re=5,
      Vp=60, Vt=100, Vz=7, Ve=40, Vb=30,
      h=0.2, w=Nd, r=1.0,
      en_dynret=False, en_dynext=False, en_optmov=False, en_verb=False):
    # Temp and Fan settings.
    self.Te, self.Tp, self.Fe, self.Fc = Te, Tp, Fe, Fc
    # Linear advance settings.
    self.Kf, self.Kb, self.Cb, self.Re = Kf, Kb, Cb, Re
    # Default velocities.
    self.Vp, self.Vt, self.Vz, self.Ve, self.Vb = Vp, Vt, Vz, Ve, Vb
    # Default line height, width, and extrusion ratio.
    self.h, self.w, self.r = h, w, r
    # Processing mode options.
    self.en_dynret, self.en_dynext, self.en_optmov, self.en_verb = en_dynret, en_dynext, en_optmov, en_verb
    class GMove(Move, Fd=self.Fd, Nd=self.Nd, log=self.log): pass
    self.GMove = GMove
    self.resetfile()
    self.startFile()

  def resetfile(self):
    """ Reset all state variables. """
    # Zero all positions so the first move's deltas are actually absolute.
    self.x, self.y, self.z, self.e, self.f = 0.0, 0.0, 0.0, 0.0, 0
    self.t = self.l = 0.0
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
    self.t += dt
    self.l += dl if en else 0.0
    self.l_t += dt
    self.l_l += dl if en else 0.0
    self.l_e += en
    self.vl = v1
    self.ve = de/dl*v1 if dl else de/dt if dt else 0.0
    self.vn = self._vn()
    self.db = self._db()

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
    dx,dy,dz,de,v,v0,v1,*_ = m
    # Don't include args in cmd if they are unchanged.
    # Include both x and y if either are changed.
    x = self.x + dx if dx or dy else None
    y = self.y + dy if dy or dx else None
    z = self.z + dz if dz else None
    e = self.e + de if de else None
    f = round(60 * v)
    if f == self.f: f = None
    cmd = self.fcmd('G1', X=x, Y=y, Z=z, E=e, F=f)
    if self.en_verb:
      cmd = f'{cmd:40s}; {m}'
    return cmd

  def fadd(self, code):
    """ Do a final format to string and add. """
    if isinstance(code, Layer):
      self.inclayer(code)
    elif isinstance(code, Move):
      out = self.fmove(code)
      self.incmove(code)
      self.add(out)
      if self.en_verb:
        self.log('{pe=:.4f} {eb=:.4f} {ve=:.4f} {vn=:.4f} {vl=:.3f} {db=:.2f}')
        # We need to fstr() reformat that last log entry...
        self.gcode[-1] = self.fstr(self.gcode[-1])
    else:
      self.add(self.fstr(code))

  def add(self, code):
    """ Add code Layer, Move or format-string instances. """
    if isinstance(code, Layer):
      self.inclayer(code)
    elif isinstance(code, Move):
      self.incmove(code)
    else:
      code = code.strip()
    self.gcode.append(code)

  def cmd(self, cmd, **kwargs):
    self.add(self.fcmd(cmd,**kwargs))

  def cmt(self, cmt):
    self.add(f';{cmt}')

  def log(self, txt):
    if self.en_verb:
      self.cmt('{ftime(t)}: ' + txt)

  def startFile(self):
    self.add(self.fstr(self.startcode))
    self.startLayer(0, 0.00)

  def endFile(self):
    # Retract to -2.5mm for the next print.
    self.retract(de=-2.5)
    self.add(self.fstr(self.endcode))

  def startLayer(self, n=None, z=None,
      Vp=None, Vt=None, Vz=None, Ve=None, Vb=None,
      h=None, w=None, r=1.0):
    # Note we always start layers from a high "hop" position.
    l = Layer(
        self.layer.n + 1 if n is None else n,
        self.layer.z + self.layer.h if z is None else z,
        Vp or self.Vp,
        Vt or self.Vt,
        Vz or self.Vz,
        Ve or self.Ve,
        Vb or self.Vb,
        h or self.h,
        w or self.w,
        r*self.r)
    self.add(l)
    self.cmt('{"preExtrude" if not layer.n else "layer"}:{layer.h:.2f}')
    self.log('layer={layer.n} h={layer.h:.2f} w={layer.w:.2f} r={layer.r:.2f}')

  def endLayer(self, **upargs):
    self.up(**upargs)
    self.layerstats()

  def _calc_pe(self, db, ve, dt=0.0):
    """ Get steady-state pressure advance pe for db and ve.

    If dt is provided, the pe is a reduced starting pressure where the move
    should get close to the target pressure within dt time.
    """
    pb = self.Kb*db + self.Cb
    pn = self.Kf*ve * e**(-dt/self.Kf)
    return pn+pb

  def _db(self, eb=None, h=None):
    """ Get bead diameter from bead volume."""
    if h is None: h = self.layer.h
    if eb is None: eb = self.eb
    assert eb >= 0.0
    return 2*sqrt(eb*self.Fa/(h*pi))

  def _vn(self,pe=None,eb=None,h=None):
    """Get nozzle flow velocity from pe, eb, and h."""
    if pe is None: pe = self.pe
    if eb is None: eb = self.eb
    if h is None: h = self.layer.h
    db = self._db(eb,h)                # diameter of bead
    pb = max(0, self.Kb*db + self.Cb)  # bead backpressure.
    pn = max(0, pe - pb)               # nozzle pressure.
    return pn/self.Kf if self.Kf else self.ve

  def _calc_pe_eb_en(self, m):
    """ Calculate updated pe, eb, and en values after a move."""
    dx,dy,dz,de,v,v0,v1,h,w = m
    if self.Kf:
      # Calculate change in pressure advance 'pe' and bead volume `eb` by
      # iterating over the time interval with a dt that is small relative to
      # the Kf timeconstant.
      dt, dl = m.dt, m.dl
      if not dl:
        # A hop, drop, retract, restore, or set speed.
        assert v0 == v1
        v, ve, de_dl = 0.0, de/dt, None
        # There is only one phase.
        phases =((dt,0),)
      else:
        v, ve, de_dl = v0, None, de/dl
        # Calculate accelerate, move, decelerate phases.
        phases = (m.dt0, m.a0),(m.dtm, m.am),(m.dt1, m.a1)
      pe, eb, el = self.pe, self.eb, 0.0
      for t,a in phases:
        while t>0:
          dt = min(0.001, t)   # This iteration's time interval.
          dv = a*dt           # do trapezoidal integration for dl.
          dl = (v + dv/2)*dt
          v += dv
          # Note filament doesn't get sucked back into the nozzle if the pressure
          # advance is negative, and the bead cannot give negative backpressure.
          pe += de_dl * dl if dl else ve*dt   # Add filament extruded into the nozzle.
          db = self._db(eb)                   # diameter of bead
          pb = max(0, self.Kb*db + self.Cb)   # bead backpressure.
          pn = max(0, pe - pb)                # nozzle pressure.
          nde = pn*dt/(self.Kf+dt)            # Get filament extruded out the nozzle.
          lde = min(db*h*dl/self.Fa, eb+nde)  # Get volume removed from bead to line,
          pe -= nde
          eb = max(0.0, eb + nde - lde)
          el += lde
          t -= dt
      assert abs(v - m.v1) < 0.001, f'{m} {v=} != {m.v1=}'
    else:
      ex = self.pe + de
      pe, eb, el = min(ex, 0), 0, max(0, ex)
    assert abs((pe+eb+el)-(self.pe+self.eb+de)) < 0.00001, f'{pe+eb+el} != {self.pe+self.eb+de}'
    # Note we return pe,eb,en not pe,eb,el, because "what came out the nozzle"
    # is more useful than "what came out the nozzle not including the bead".
    # We also recalculate en from the raw inputs to avoid accumulated integration errors.
    return pe, eb, self.pe + de - pe

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

    Note this is the simple raw movement. Any pressure advance compensation
    will be added in post-processing.
    """
    if z is None and dz is None and h is not None: z = self.layer.z + h
    dx = x - self.x if x is not None else dx if dx is not None else 0.0
    dy = y - self.y if y is not None else dy if dy is not None else 0.0
    dz = z - self.z if z is not None else dz if dz is not None else 0.0
    x,y,z = self.x + dx, self.y + dy, self.z + dz
    dl = sqrt(dx**2 + dy**2)
    # Note: setting h and z or dz means explicit h overrides z height for de calcs.
    if h is None: h = z - self.layer.z
    if w is None: w = self.layer.w
    r = self.layer.r * r
    el = dl * w * h / self.Fa  # line filament volume not including extrusion ratio.
    de = e - self.e if e is not None else de if de is not None else r * el
    if v is None:
      if de:
        v = s * (self.layer.Vp if dl else self.layer.Ve if de < 0 else self.layer.Vb)
      else:
        v = s * (self.layer.Vt if dl else self.layer.Vz)
    m = self.GMove(dx,dy,dz,de,v,h=h,w=w)
    # Only add it if it does anything.
    # Note we don't add speed changes, as moves all keep their own speeds.
    if any(m[:4]):
      self.add(m)
    else:
      self.log(f'skipping empty move {m}')

  def draw(self, x=None, y=None, r=1.0, **kwargs):
    """ Draw a line.

    This is the same as move() but with the default r=1.0 so it prints instead
    of travels by default.
    """
    self.move(x, y, r=r, **kwargs)

  def retract(self, e=None, de=None, ve=None, s=1.0):
    """Do a retract.

    This is a simple retraction with default self.Re. Note compensation for
    advance pressure is added in post processing if enabled.
    """
    if de is None: de = -self.Re
    self.move(e=e, de=de, v=ve, s=s)

  def restore(self, e=None, de=None, vb=None, s=1.0):
    """ Do a restore.

    This is a simple restore with default self.Re. Note adding starting dots
    and compensating for advance pressure is added in post-processing if
    enabled.
    """
    if de is None: de = self.Re
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
    # default hop down is to layer height.
    if h is None: h = self.layer.h
    if (x, y, dx, dy) != (None, None, None, None):
      self.move(x=x, y=y, dx=dx, dy=dy, v=vt, s=s)
    self.move(z=z, dz=dz, v=vz, h=h, s=s)
    self.restore(e=e, de=de, vb=vb, s=s)

  def getCode(self):
    if self.en_optmov:
      # Join together any moves that can be joined.
      ljoins(self.gcode)
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
      if isinstance(c,Move):
        if c.isretract and self.en_dynret:
          # Adjust retraction to also relieve advance pressure.
          c = c.set(de=c.de - self.pe)
        elif c.isrestore and self.en_dynret:
          # Find the next move to get the pe and eb it needs.
          m1 = next((m for m in gcode[i+1:] if isinstance(m, Move)), None)
          # For calculating the pressure pe needed, use the starting ve0 for
          # dynamic extrusion accelerations, otherwise use average ve.
          if m1 and m1.isdraw:
            ve = m1.ve0 if self.en_dynext and m1.isaccel else m1.ve
            eb, pe = m1.eb, self._calc_pe(m1.db, ve)
            #self.log(f'adjusting {c} for next {m1}, {m1.ve=:.4f}({m1.ve0:.4f}<{m1.vem:.4f}>{m1.ve1:.4f}) {m1.isaccel=}')
          else:
            eb, pe = 0.0, self._calc_pe(0.0, 0.0)
          # Adjust existing retraction and add the starting bead.
          c = c.set(de=pe - self.pe + eb)
        elif c.isdraw and self.en_dynext:
          # For calculating the pressure pe needed, use the ending ve1 for
          # accelerations, otherwise use average ve.
          ve = c.ve1 if c.isaccel else c.ve
          pe = self._calc_pe(c.db, ve)
          # Adjust de to include required change in pe over the move.
          #self.log(f'adjusting {c} {c.ve=:.4f}({c.ve0:.4f}<{c.vem:.4f}>{c.ve1:.4f}) {c.isaccel=}')
          c = c.set(de=c.de + pe - self.pe)
      self.fadd(c)
    return '\n'.join(self.gcode) + '\n'

  def layerstats(self):
    h='{layer.h:.2f}'
    w='{layer.w:.2f}'
    r='{(Fa * l_e) / (layer.h * layer.w * l_l):.2f}'
    v='{l_l/l_t if l_t else 0.0:.1f}'
    ve='{l_e/l_t if l_t else 0.0:.3f}'
    self.log('layer={layer.n} finished.')
    self.log(f'  avg:{h}x{w}x{r}@{v} ve={ve}')
    self.log('  t={l_t:.1f} l={l_l:.1f} e={l_e:.1f}')

  def preExt(self,x0,y0,x1,y1,le=120.0, lw=100.0, m=10.0, ve=20.0, vw=10, r=4):
    """Preextrude around a box with margin m.

    This extrudes for `le` mm at speed `ve` mm/s, then wipes for another `lw`
    mm at speed `vw` mm/s without extruding, and finally another 10mm at 1mm/s
    to stick the last bit of drool to the plate. This should leave the nozzle
    primed with no residual advance pressure, and prevent the wipe from
    becoming a big piece of stringing that messes with the print.
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
  tn = 10   # number of test runs for each test.
  dvx = 2.5  # test speed step
  dle = 0.2  # test extraction length step
  dve = 1  # test extraction speed step
  t0x = 5  # test line warmup distance
  t1x = 20  # test line phase distance.
  rdx = 2*t0x + 4*t1x # test ruler length.
  rdy = 4  # ruler height including space
  sdy = 30 # test text settings box width.
  tdx = rdx + sdy  # test box total width.
  tdy = rdy + tn*ty + 5 # test box total height.

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

  def _fval(self,k,v):
    """ format a single gcode command argument."""
    v=f'{v[0]}-{v[1]}' if isinstance(v, tuple) else f'{round(v,4)}'
    return f'{k}={v}'

  def _fset(self, sep='\n', **kwargs):
    return sep.join(self._fval(k,v) for k,v in kwargs.items())

  def settings(self, x0, y0, **kwargs):
    self.cmt("structure:brim")
    self.text(self._fset(**kwargs), x0, y0, fsize=5)

  def doRetractTests(self):
    t4x, t4y = self.tdx, 4*self.tdy - 5
    x0 = -self.rdx/2
    y0 = 55
    x1 = x0+t4x
    y1 = y0-t4y
    self.preextrude()
    self.startLayer(z=0.0, Vp=10)
    self.brim(x0,y0,x0+self.rdx, y1)
    self.testRetract(x0,y0,Kf=0.0,vxr=20,ler=(-1,4),ver=self.Ve,vex=0.0)
    self.testRetract(x0,y0-self.tdy,Kf=1.0,vxr=20,ler=(-1,4),ver=self.Ve,vex=0.0)
    self.testRetract(x0,y0-2*self.tdy,Kf=2.0,vxr=20,ler=(-1,4),ver=self.Ve,vex=0.0)
    self.testRetract(x0,y0-3*self.tdy,Kf=3.0,vxr=20,ler=(-1,4),ver=self.Ve,vex=0.0)
    #self.testRetract(x0,y0,vxr=(self.vx0, self.vx1),ler=self.Re,ver=self.Ve)
    #self.testRetract(x0,y0-self.tdy,vxr=self.Vp,ler=(self.le0,self.le1),ver=self.Ve)
    #self.testRetract(x0,y0-2*self.tdy,vxr=self.Vp,ler=self.Re, ver=(self.ve0,self.ve1))
    self.endLayer()

  def _getstep(self,tn,vr):
    """ Convert a value-range into start, stop, delta. """
    if not isinstance(vr, tuple):
      return vr, vr, 0
    v0, v1 = vr
    dv = (v1-v0)/tn
    return v0, v1, dv

  def testRetract(self, x0=None, y0=None, Kf=None,
      vxr=vx0, ler=le0, ver=ve0, vex=None):
    vx, vx1, dvx = self._getstep(self.tn,vxr)
    le, le1, dle = self._getstep(self.tn,ler)
    ve, ve1, dve = self._getstep(self.tn,ver)
    logset = self._fset(sep=' ', vx=vxr, le=ler, ve=ver)
    self.log(f'testRetract {logset}')
    self.ruler(x0, y0)
    y = y0 - self.rdy
    self.cmt(f'structure:infill-solid')
    while vx <= vx1 and le <= le1 and ve <= ve1:
      self.testRetractLine(x0, y, Kf, vx, le, ve, vex)
      vx+=dvx
      le+=dle
      ve+=dve
      y-=self.ty
    if vex is None: vex='2*ve'
    self.settings(x0 + self.rdx+1, y0, Kf=Kf, vx=vxr, le=ler, ve=ver, vex=vex)

  def testRetractLine(self, x0, y0, Kf, vx, le=0, ve=35, vex=None):
    """ Test starting and stopping extrusion for different settings. """
    # A printer with a=500mm/s^2 can go from 0 to 100mm/s in 0.2s over 10mm.
    # Retracting de=10mm at ve=20mm/s will take 0.5s, which is 50mm at 100mm/s
    # So retracting over a distance of 25mm is max le=5mm@20mm/s or le=10mm@40mm/s.
    # Alternatively we can reduce our speed to 2*ve while extracting.

    # For vex use vx to retract without slowing down, 2*ve to retract at reduced
    # speed, or 0 to stop for retraction.
    oldKf = self.Kf
    if Kf is not None:
      self.Kf = Kf
    if vex is None: vex = 2*ve  # moving speed while retracting.
    self.log(f'testRetractLine {Kf=} v={vx} le={le} ve={ve} vex={vex}')
    dt = le/ve  # time to do retraction
    ex = vex*dt # distance travelled while retracting.
    assert ex <= self.t1x
    self.dn(x0, y0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.draw(dx=self.t1x,v=vx)
    if vex == 0.0:
      self.retract(de=-le, ve=ve)  # retract while stopped.
      self.move(dx=2*self.t1x,v=vx)    # move gap between moving retract/restore.
      self.restore(vb=ve)   # restore while stopped.
    else:
      self.move(dx=ex,v=vex,de=-le)  # retract while moving at v=vex.
      self.move(dx=2*(self.t1x-ex),v=vex) # move gap between moving retract/restore.
      self.move(dx=ex,v=vex,de=le)   # restore while moving at v=vx
    self.draw(dx=self.t1x,v=vx)
    self.draw(dx=self.t0x,v=self.vx0)
    self.up()
    self.Kf = oldKf


  def doStartStopTests(self):
    t4x, t4y = self.tdx, 4*self.tdy - 5
    x0, y0 = -self.rdx/2, 55
    x1, y1 = x0+t4x, y0-t4y
    self.preextrude()
    self.startLayer(z=0.0, Vp=10)
    self.brim(x0,y0,x0+self.rdx, y1)
    self.testStartStop(x0,y0,Kfr=0.8,vxr=(10,100),ler=0.0,tn=9)
    self.testStartStop(x0,y0-self.tdy,Kfr=0.8,vxr=(10,100),ler=1.0,tn=9)
    self.testStartStop(x0,y0-2*self.tdy,Kfr=0.8,vxr=60,ler=(0.0,2.0))
    self.testStartStop(x0,y0-3*self.tdy,Kfr=0.8,vxr=100,ler=(0.0,2.0))
    self.endLayer()

  def testStartStop(self, x0=None, y0=None, Kfr=0.0,
      vxr=vx0, ler=le0, tn=tn):
    Kf, Kf1, dKf = self._getstep(tn,Kfr)
    vx, vx1, dvx = self._getstep(tn,vxr)
    le, le1, dle = self._getstep(tn,ler)
    assert dKf or dvx or dle
    logset = self._fset(sep=' ', Kf=Kfr, vx=vxr, le=ler)
    self.log(f'testStartStop {logset}')
    self.ruler(x0, y0)
    y = y0 - self.rdy
    self.cmt(f'structure:infill-solid')
    e=0.00001
    while Kf <= Kf1+e and vx <= vx1+e and le <= le1+e:
      self.testStartStopLine(x0, y, Kf, vx, le)
      Kf+=dKf
      vx+=dvx
      le+=dle
      y-=self.ty
    self.settings(x0 + self.rdx+1, y0, Kf=Kfr, vx=vxr, le=ler)

  def testStartStopLine(self, x0, y0, Kf, vx, le=0):
    """ Test starting and stopping extrusion for different settings. """
    # if le is set to 0, set it to a tiny value so that retract/restore don't
    # get optimized away.
    if not le: le=0.001
    # A printer with a=500mm/s^2 can go from 0 to 100mm/s in 0.2s over 10mm.
    oldKf,self.Kf = self.Kf, Kf
    self.log(f'testStartStopLine {Kf=:.2f} vx={round(vx)} {le=:.2f}')
    self.dn(x0, y0)
    # This slow draw and move is to ensure the nozzle is primed and wiped, and
    # any advance pressure, actual or estimated, has decayed away.
    self.draw(dx=5,v=1)
    self.move(dx=10,v=1)
    # This restore pre-loads advance pressure for drawing the line. If the
    # start of the line is too thin, there was not enough pressure, and Kf
    # probably needs to be increased. If the start is too thick, pressure was
    # too high and Kf should probably be dropped. The line should be
    # consistently the same width as all the other lines.
    self.restore(de=le)
    self.draw(dx=50,v=vx)
    # This retract and slow move is to relieve the advance pressure, and see
    # if any remaining pressure drools off. If there is any trailing drool,
    # the amount of pressure was underestimated and Kf probably needs to be
    # increased. If there is still stringing when otherwise Kf seems right,
    # adding some additional retraction with de could help.
    self.retract(de=-le)
    self.move(dx=10,v=1)
    # This restore applies low advance pressure for a final slow line after
    # the earlier retraction. It should be the same width as all the other
    # lines. If it starts too thin, it suggests the earlier retraction
    # over-estimated the pressure and over-retracted, so Kf should be reduced.
    # If it starts too thick it suggests the earlier retraction under
    # estimated the pressure and extra retraction was actually relieving
    # pressure, so K should be increased and de could possibly be reduced.
    self.restore(de=le)
    self.draw(dx=15,v=5)
    # This starts drawing at 1mm/sec extruding at 0.1mm/sec for 15mm to see
    # how much over-retraction there was. When the retraction is restored the
    # line should be 34% thicker than the previous fast line. Every mm without
    # without a thick line minus 0.1*Kf is 0.1mm of overretraction.
    #self.draw(dx=15, de=1.5, v=1)
    self.up()
    self.Kf = oldKf


  def doKfTests(self, n=1, Kf=0.0):
    n -= 1 # test number for offset indexes starts at 0.
    x0,y0 = -self.rdx/2, 55 - n*self.tdy
    self.preextrude(n)
    self.startLayer(z=0, Vp=10)
    # Only draw the brim on the first test.
    if n == 0:
      x1,y1 = x0+self.rdx, y0-4*self.tdy+5
      self.brim(x0,y0,x1,y1)
    self.testKf(x0,y0,Kf)
    self.endLayer()

  def testKf(self, x0, y0, Kf=0.0, vxr=(vx0,vx1)):
    vx, vx1, dvx = self._getstep(self.tn, vxr)
    logset = self._fset(sep=' ', Kf=Kf, vx=vxr)
    self.log(f'testKf {logset}')
    self.ruler(x0, y0)
    y = y0 - self.rdy
    self.cmt(f'structure:infill-solid')
    while vx <= vx1:
      self.testKfLine(x0, y, vx0=20, vx1=vx)
      vx+=dvx
      y-=1
    self.settings(x0 + self.rdx+1, y0, Kf=self.Kf, vx0=20, vx=vxr)

  def testKfLine(self, x0, y0, vx0, vx1):
    """ Test linear advance changing speed different speeds. """
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


if __name__ == '__main__':
  import argparse, sys

  inf=float('inf')
  def RangeType(min=-inf, max=inf):
    def init(s):
      v=type(min)(s)
      if v<min or max<v:
        raise argparse.ArgumentTypeError(f'must be between {min} and {max}')
      return v
    return init

  cmdline = argparse.ArgumentParser(description='Generate test gcode.')
  cmdline.add_argument('-Te', type=RangeType(0, 265), default=245,
      help='Extruder temperature.')
  cmdline.add_argument('-Tp', type=RangeType(0, 100), default=100,
      help='Platform temperature.')
  cmdline.add_argument('-Fe', type=RangeType(0.0,1.0), default=0.0,
      help='Extruder fan speed between 0.0 to 1.0.')
  cmdline.add_argument('-Fc', type=RangeType(0.0,1.0), default=0.0,
      help='Case fan speed between 0.0 to 1.0.')
  cmdline.add_argument('-Kf', type=RangeType(0.0,4.0), default=0.4,
      help='Linear Advance factor between 0.0 to 4.0 in mm/mm/s.')
  cmdline.add_argument('-Kb', type=RangeType(0.0,10.0), default=2.0,
      help='Bead backpressure factor between 0.0 to 10.0 in mm/mm.')
  cmdline.add_argument('-Cb', type=RangeType(-5.0,5.0), default=0.8,
      help='Bead backpressure offset between -5.0 to 5.0 in mm.')
  cmdline.add_argument('-Re', type=RangeType(0.0,10.0), default=1.0,
      help='Retraction distance between 0.0 to 10.0 in mm.')
  cmdline.add_argument('-Vp', type=RangeType(5,100), default=50,
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
  cmdline.add_argument('-Lw', type=RangeType(0.2,0.8), default=0.6,
      help='Line width in mm.')
  cmdline.add_argument('-Lr', type=RangeType(0.0,10.0), default=1.0,
      help='Line extrusion ratio between 0.0 to 10.0.')
  cmdline.add_argument('-R', action='store_true',
      help='Enable Linear Advance dynamic retract/restore.')
  cmdline.add_argument('-P', action='store_true',
      help='Enable Linear Advance dynamic extrusion rates.')
  cmdline.add_argument('-O', action='store_true',
      help='Enable optimize, merging draw/move cmds where possible.')
  cmdline.add_argument('-v', action='store_true',
      help='Enable verbose logging output.')
  cmdline.add_argument('-n', type=RangeType(1,5), default=1,
      help='Test number for linear advance tests.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc,
      Kf=args.Kf, Kb=args.Kb, Cb=args.Cb, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve, Vb=args.Vb,
      h=args.Lh, w=args.Lw, r=args.Lr,
      en_dynret=args.R, en_dynext=args.P, en_optmov=args.O, en_verb=args.v)
  #gen.doKfTests(args.n, args.Kf)
  #gen.doRetractTests(args.n, args.Kf)
  gen.doStartStopTests()
  gen.endFile()
  data=gen.getCode()
  sys.stdout.write(data)
