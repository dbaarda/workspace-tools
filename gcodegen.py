#!/bin/python3
"""

"""
import re
from math import pi, inf
import vtext

nl = '\n'  # for use inside fstrings.

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
  return eval(f'f"""{s}"""',globals,locals)


class GCodeGen(object):
  """ GCodeGen gcode generator.

  Attributes:
    Te: Extruder temp (degC).
    Tp: Platform temp (degC).
    Fe: Extruder fan speed (0.0->1.0).
    Fc: Case fan speed (0.0->1.0).
    Kf: Linear Advance factor (0.0->4.0 mm/mm/s).
    Re: Retraction distance (0.0->10.0 mm).
    Vp: Speed when printing (mm/s).
    Vt: Speed when traveling (mm/s).
    Vz: Speed when raising/lowering (mm/s).
    Ve: Speed when retracting/recovering (mm/s).
    r: Extrusion ratio (0.1 -> 10.0).
    h: Line height (0.1->0.4mm).
    w: Line width (0.2->0.8mm).
  """

  Fd = 1.75  # Fillament diameter.
  Nd = 0.4   # Nozzle diameter.
  Fa = pi * (Fd/2)**2  # Fillament area.
  Na = pi * (Nd/2)**2  # Nozzle area.
  Zh = Nd    # z-hop height.

  argfmt = dict(X='.2f',Y='.2f',Z='.3f',E='.4f',F='d')

  startcode = """\
;gcode_flavor: flashforge
;machine_type: Adventurer 3 Series
;filament_diameter0: {self.Fd}
;right_extruder_temperature: {self.Te}
;platform_temperature: {self.Tp}
;extruder_fan: {int(self.Fe*100)}%
;case_fan: {int(self.Fc*100)}%
;base_print_speed: {int(self.Vp)}
;travel_speed: {int(self.Vt)}
;vertical_speed: {int(self.Vz)}
;extruder_speed: {int(self.Ve)}
;layer_height: {self.f_h:.2f}
;line_width: {self.f_w:.2f}
;extrusion_ratio: {self.f_r:.2f}
;start gcode
M118 X{75.0:.2f} Y{75.0:.2f} Z{75.0:.2f} T0
M140 S{int(self.Tp)} T0
M104 S{int(self.Te)} T0
M104 S0 T1
M107
M900 K{self.Kf:.3f} T0
G90
G28
M132 X Y Z A B
G1 Z50.000 F420
G161 X Y F3300
M7 T0
M6 T0
{self._setFc(self.Fc)}
M108 T0
{self._setFe(self.Fe)+nl if self.Fe else ""}\
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
    return f'M652 S{int(fc*255)}' if fc else 'M651'

  def _setFe(self, fe):
    return f'M106' if fe == 1.0 else f'M106 S{int(fe*255)}' if fe else 'M107'

  def __init__(self,
      Te=210, Tp=50, Fe=1.0, Fc=1.0, Kf=0.0, Re=5,
      Vp=60, Vt=100, Vz=7, Ve=35,
      h=0.2, w=Nd, r=1.0):
    self.Te, self.Tp, self.Fe, self.Fc = Te, Tp, Fe, Fc
    self.Kf, self.Re = Kf, Re
    self.Vp, self.Vt, self.Vz, self.Ve = Vp, Vt, Vz, Ve
    self.f_h, self.f_w, self.f_r = h, w, r  # Default line height, width, and extrusion ratio.
    self.gcode=[]
    self.t = self.l = 0.0
    self.x = self.y = self.z = self.e = 0.0
    self.f = 0
    self.v, self.h, self.w, self.r = 0,0,0,0 # Last write's velocity, height, width, and ratio
    self.re = 0  # current retracted de amount.
    self.startFile()

  def _fmt(self,k,v):
    return f'{v:{self.argfmt.get(k,"")}}'

  def add(self, code):
    self.gcode.append(code.strip())

  def cmd(self, cmd, **kwargs):
    args = ' '.join(f'{k}{self._fmt(k,v)}' for k,v in kwargs.items() if v is not None)
    cmd = f'{cmd} {args}' if args else cmd
    self.add(cmd)

  def cmt(self, cmt):
    self.add(f';{cmt}')

  def log(self, txt):
    self.cmt(f'{ftime(self.t)}: {txt}')

  def startFile(self):
    self.add(fstr(self.startcode,locals=locals()))
    self.startLayer('P')

  def endFile(self):
    # Restore any remaining retraction minus 2.5mm for the next print.
    self.restore(de=self.re-2.5,r=0)
    self.add(fstr(self.endcode, locals=locals()))

  def startLayer(self, n=None, Vp=None, h=None, w=None, r=1.0, **dnargs):
    # Note we always start layers from a high "hop" position.
    try:
      self.l_n = n if n else self.l_n + 1
      self.l_z += self.l_h
    except (TypeError, AttributeError):
      # The previous layer was pre-emissions or uninitialized.
      self.l_n = n if n else 1
      self.l_z = 0.0
    self.l_Vp = self.Vp if Vp is None else Vp
    self.l_h = self.f_h if h is None else h
    self.l_w = self.f_w if w is None else w
    self.l_r = self.f_r * r
    self.l_t = self.l_l = self.l_e = 0.0
    self.cmt(f'{"preExtrude" if self.l_n == "P" else "layer"}:{self.l_h:.2f}')
    self.log(f'layer={self.l_n} h={self.l_h:.2f} w={self.l_w:.2f} r={self.l_r:.2f}')
    if dnargs:
      self.dn(h=h, w=w, **dnargs)

  def endLayer(self, **upargs):
    self.up(**upargs)
    self.layerstats()

  def move(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      v=None, s=1.0, r=0.0, w=None, h=None):
    """ Move the extruder.

    Movement can be absolute using `x,y,z,e` arguments, or relative using
    `dx,dy,dz,de` arguments. The z axis can also be specified with `h` height
    above the base of the current layer. The extrusion can also be specified
    using an `r` extrusion ratio muliplier for a `h*w` line, where `h` and `w`
    default to the current layer's if not specified. Speed can be specified
    directly using the `v` velocity argument, or the `s` speed multiplier of
    the default `Vp,Vt,Vz,Ve` speeds for print, travel, vertical, and extruder
    movement.
    """
    if x is None:
      x = None if dx is None else self.x + dx
    if y is None:
      y = None if dy is None else self.y + dy
    if z is None:
      z = None if dz is None else self.z + dz
    if z is None:
      z = None if h is None else self.l_z + h
    dx = 0.0 if x is None else x - self.x
    dy = 0.0 if y is None else y - self.y
    dz = 0.0 if z is None else z - self.z
    dl=(dx**2 + dy**2 + dz**2)**0.5
    if h is None: h = self.l_h
    if w is None: w = self.l_w
    r = self.l_r * r
    a = w * h
    if e is None:
      e = None if de is None else self.e + de
    if e is None:
      e = None if not r else self.e + r * dl * a / self.Fa
    de = 0.0 if e is None else e - self.e
    r = round(de * self.Fa / (dl * a), 4) if dl else None
    if v is None:
      if de:
        v = s * (self.l_Vp if dl else self.Ve)
      else:
        v = s * (self.Vt if abs(dz) < self.Vz/self.Vt*dl else self.Vz)
    f = int(60 * v)
    dt=(dl if dl else abs(de))/v
    # Don't include args in cmd if they are unchanged.
    if x == self.x: x = None
    if y == self.y: y = None
    if z == self.z: z = None
    if e == self.e: e = None
    if f == self.f: f = None
    # Include both X and Y if either are included.
    if x is None and y is not None: x = self.x
    if y is None and x is not None: y = self.y
    if (x,y,z,e,f) != (None, None, None, None, None):
      # if extruding with changed settings, log and update.
      if r and (v,h,w,r) != (self.v, self.h, self.w, self.r):
        self.cmt(f'draw:{h:.2f}x{w:.2f}x{r:.2f}@{int(dl/dt)} ve={de/dt:.2f}')
        self.v, self.h, self.w, self.r = v, h, w, r
      self.cmd('G1', X=x, Y=y, Z=z, E=e, F=f)
      self.x += dx
      self.y += dy
      self.z += dz
      self.e += de
      self.f = self.f if f is None else f
      # Get the amount actually extruded after restoring.
      xe = self.re + de
      self.re, xe = min(0 , xe), max(0, xe)
      # Only count de and dl for extruding moves in stats.
      dl = dl if xe else 0.0
      self.l_e += xe
      self.l_l += dl
      self.l_t += dt
      self.l += dl
      self.t += dt

  def draw(self, x=None, y=None, r=1.0, **kwargs):
    """ Draw a line.

    This is the same as move() but with the default r=1.0 so it prints instead
    of travels by default.
    """
    self.move(x, y, r=r, **kwargs)

  def retract(self, e=None, de=None, v=None, s=1.0):
    if de is None: de = -self.Re
    self.move(e=e, de=de, v=v, s=s)

  def restore(self, e=None, de=None, v=None, s=1.0, r=1.0, w=None, h=None):
    # default restore is default retraction plus enough to half-fill starting dot.
    # Set r=0 to turn this off.
    if de is None: de = self.Re
    if h is None: h = self.l_h
    if w is None: w = self.l_w
    r = self.l_r * r
    sv = r * 0.5 * pi*(w/2)**2 * h # volume to half-fill starting dot.
    se = sv/self.Fa # extra extrusion to half-fill starting dot.
    de += se
    self.move(e=e, de=de, v=v, s=s)

  def up(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      v=None, s=1.0, h=None):
    """ Do a retract, raise, and move. """
    # default hop up is to Zh above layer height.
    if h is None: h = self.l_h + self.Zh
    self.retract(e=e, de=de, s=s)
    self.move(z=z, dz=dz, h=h, s=s)
    if (x, y, dx, dy) != (None, None, None, None):
      self.move(x=x, y=y, dx=dx, dy=dy, v=v, s=s)

  def dn(self,
      x=None, y=None, z=None, e=None,
      dx=None, dy=None, dz=None, de=None,
      v=None, s=1.0, r=1.0, w=None, h=None):
    # default hop down is to layer height.
    if h is None: h = self.l_h
    if (x, y, dx, dy) != (None, None, None, None):
      self.move(x=x, y=y, dx=dx, dy=dy, v=v, s=s)
    self.move(z=z, dz=dz, h=h, s=s)
    self.restore(e=e, de=de, s=s, r=r, w=w, h=h)

  def getCode(self):
    return '\n'.join(self.gcode + [''])

  def layerstats(self):
    h=f'{self.l_h:.2f}'
    w=f'{self.l_w:.2f}'
    r=f'{(self.Fa * self.l_e) / (self.l_h * self.l_w * self.l_l):.2f}'
    v=f'{self.l_l/self.l_t if self.l_t else 0.0:.1f}'
    ve=f'{self.l_e/self.l_t if self.l_t else 0.0:.3f}'
    self.log(f'layer={self.l_n} finished.')
    self.cmt(f' avg:{h}x{w}x{r}@{v} ve={ve}')
    self.cmt(f' t={self.l_t:.1f} l={self.l_l:.1f} e={self.l_e:.1f}')

  def preExt(self,x0,y0,x1,y1,l=120.0, m=10.0, v=20.0, r=4):
    self.cmt('structure:pre-extrude')
    xw,yw = x1-x0 + 2*m, y0-y1 + 2*m
    dy = min(l, yw) / 2
    x, y = x0 - m, y0 + m - dy
    self.dn(x=x,y=y,de=0,r=0) # drop down but don't restore.
    self.draw(dy=dy,v=v,r=r)
    l -= dy
    while l:
      dx = min(xw,l)
      self.draw(dx=dx,v=v,r=r)
      l-=dx
      dy = min(yw,l)
      self.draw(dy=-dy,v=v,r=r)
      l-=dy
      dx = min(xw,l)
      self.draw(dx=-dx,v=v,r=r)
      l-=dx
      dy = min(yw,l)
      self.draw(dy=dy,v=v,r=r)
      l-=dy
      xw+=2
      yw+=2
    self.up()

  def fbox(self,x0,y0,x1,y1,**kwargs):
    # Fill a box with frame lines.
    w = kwargs.get('w', self.l_w)
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
    w = kwargs.get('w', self.l_w)
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
    # double r so the restore extrudes a full dot.
    # Note this doesn't seem to render in FlashPrint.
    self.dn(x,y, r=r*2, **kwargs)
    # Draw a tiny line so FlashPrint renders something.
    self.draw(dy=-0.2,r=r, **kwargs)
    self.up()

  def line(self, l, r=1.0, **kwargs):
    """Draw a line from a sequence of points."""
    x0,y0 = l[0]
    # If there is only one point, draw a dot instead.
    if len(l) == 1 or (len(l) == 2 and l[0] == l[1]):
      return self.dot(x0, y0, r, **kwargs)
    self.dn(x0,y0,r=r,**kwargs)
    for x,y in l[1:]:
      self.draw(x,y,**kwargs)
    self.up(**kwargs)

  def text(self, t, x0, y0, x1=None, y1=None, fsize=5, angle=0, **kwargs):
    w = kwargs.get('w', self.l_w)
    #print(f't={t!r} {x0=} {y0=} {x1=} {y1=} {fsize=} {angle=} {w=}')
    v = vtext.ptext(t,x0,y0,x1,y1,fsize,angle,w)
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
  rdy = 5  # ruler height
  vx0 = 5  # test min speed
  vx1 = 100  # test max speed
  le0 = 0  # test min retraction length
  le1 = 10  # test max retraction length
  ve0 = 15  # test min retraction speed
  ve1 = 45  # test max retraction speed
  tn = 20   # number of test runs for each test.
  dvx = 2.5  # test speed step
  dle = 0.2  # test extraction length step
  dve = 1  # test extraction speed step
  t0x = 5  # test line warmup distance
  tlx = 20  # test line phase distance.

  def rdx(self,tlx):
    return 4*tlx + 2*self.t0x

  def dotests(self, tn=tn, tlx=tlx):
    tlx = self.rdx(tlx)
    tdy = self.rdy + tn + 5
    tly = 3*tdy - 5
    x0 = -tlx/2
    y0 = tly/2
    x1 = x0+tlx
    y1 = y0-tly
    dvx = (self.vx1 - self.vx0)/tn
    dle = (self.le1 - self.le0)/tn
    dve = (self.ve1 - self.ve0)/tn
    self.preExt(x0,y0,x1,y1)
    self.startLayer(Vp=10)
    self.brim(x0,y0,x1,y1)
    self.testrun(x0,y0,dvx=dvx,le0=self.Re,ve0=self.Ve)
    self.testrun(x0,y0-tdy,dle=dle,ve0=self.Ve,vx0=self.Vp)
    self.testrun(x0,y0-2*tdy,dve=dve,le0=self.Re,vx0=self.Vp)
    self.endLayer()

  def testrun(self, x0=None, y0=None, tlx=tlx, dvx=0, dle=0, dve=0,
      vx0=vx0, vx1=vx1, le0=le0, le1=le1, ve0=ve0, ve1=ve1, vex=None):
    assert dvx+dle+dve > 0
    if dvx:
      self.log(f'extrude-test speed={vx0}-{vx1} step={dvx}, ext-len={le0}, ext-vel={ve0}')
    if dle:
      self.log(f'extrude-test ext-len={le0}-{le1} step={dle} speed={vx0} ext-vel={ve0}')
    if dve:
      self.log(f'extrude-test ext-vel={ve0}-{ve1} step={dve} speed={vx0} ext-len={le0}')
    self.ruler(x0,y0,tlx)
    y,vx,le,ve = y0 - self.rdy, vx0, le0, ve0
    self.cmt(f'structure:infill-solid')
    while vx <= vx1 and le <= le1 and ve <= ve1:
      self.testline(x0,y,tlx,vx,le,ve,vex)
      vx+=dvx
      le+=dle
      ve+=dve
      y-=1
    self.settings(x0 + self.rdx(tlx)+1, y0-5, dvx, dle, dve, vx0, vx1, le0, le1, ve0, ve1, vex)

  def brim(self, x0, y0, x1, y1):
    w = self.l_w
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

  def ruler(self, x0, y0, tlx=tlx):
    t0x, rl = self.t0x, self.rdx(tlx)
    ml,tl=3,1.5  # tick marker lengths
    self.cmt('structure:shell-outer')
    # draw ruler
    self.dn(x0, y0)
    self.draw(dx=rl)
    self.up()
    for d in (t0x,tlx,tlx,tlx,tlx):
      self.dn(dx=-d,y=y0)
      self.draw(dy=-ml)
      self.up()
    self.move(x=x0+self.t0x + tlx, y=y0)
    for d in range(2,2*tlx,2):
      l = tl if d % 10 else ml
      self.dn(dx=2,y=y0)
      self.draw(dy=-l)
      self.up()

  def settings(self, x0, y0, dvx=0, dle=0, dve=0,
        vx0=vx0, vx1=vx1, le0=le0, le1=le1, ve0=ve0, ve1=ve1, vex=None):
    self.cmt("structure:brim")
    Kf=f'Kf={self.Kf:.3f}'
    vx=f'vx={vx0}-{vx1}' if dvx else f'vx={vx0}'
    le=f'le={le0}-{le1}' if dle else f'le={le0}'
    ve=f've={ve0}-{ve1}' if dve else f've={ve0}'
    vex=f'vex={vex}' if vex is not None else f'vex=2*ve'
    t= '\n'.join([Kf,vx,le,ve,vex])
    self.text(t, x0, y0, fsize=5)

  def testline(self, x0, y0, tlx, vx, le=0, ve=35, vex=None):
    # A printer with a=500mm/s^2 can go from 0 to 100mm/s in 0.2s over 10mm.
    # Retracting de=10mm at ve=20mm/s will take 0.5s, which is 50mm at 100mm/s
    # So retracting over a distance of 25mm is max le=5mm@20mm/s or le=10mm@40mm/s.
    # Alternatively we can reduce our speed to 2*ve while extracting.

    # For vex use vx to retract without slowing down, 2*ve to retract at reduced
    # speed, or 0 to stop for retraction.
    if vex is None: vex = 2*ve  # moving speed while retracting.
    self.cmt(f'testline v={vx} le={le} ve={ve} vex={vex}')
    dt = le/ve  # time to do retraction
    ex = vex*dt # distance travelled while retracting.
    assert ex <= tlx
    self.dn(x0, y0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.draw(dx=tlx,v=vx)
    if not le:
      self.move(dx=2*tlx,v=vx)    # move without retract/restore.
    elif vex == 0.0:
      self.retract(de=-le, v=ve)  # retract while stopped.
      self.move(dx=2*tlx,v=vx)    # move gap between moving retract/restore.
      self.restore(de=le, v=ve)   # restore while stopped.
    else:
      self.move(dx=ex,v=vex,de=-le)  # retract while moving at v=vex.
      self.move(dx=2*(tlx-ex),v=vex) # move gap between moving retract/restore.
      self.move(dx=ex,v=vex,de=le)   # restore while moving at v=vx
    self.draw(dx=tlx,v=vx)
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
  cmdline.add_argument('-Kf', type=RangeType(0.0,4.0), default=0.0,
      help='Linear Advance factore between 0.0 to 4.0 in mm/mm/s.')
  cmdline.add_argument('-Re', type=RangeType(0.0,10.0), default=5.0,
      help='Retraction distance between 0.0 to 10.0 in mm.')
  cmdline.add_argument('-Vp', type=RangeType(5,100), default=50,
      help='Base printing speed in mm/s.')
  cmdline.add_argument('-Vt', type=RangeType(5,100), default=100,
      help='Base travel speed in mm/s.')
  cmdline.add_argument('-Vz', type=RangeType(1,10), default=7,
      help='Base raise/lower speed in mm/s.')
  cmdline.add_argument('-Ve', type=RangeType(1,50), default=30,
      help='Base retract/restore speed in mm/s.')
  cmdline.add_argument('-Lh', type=RangeType(0.1,0.4), default=0.3,
      help='Layer height in mm.')
  cmdline.add_argument('-Lw', type=RangeType(0.2,0.8), default=0.5,
      help='Line width in mm.')
  cmdline.add_argument('-Lr', type=RangeType(0.0,10.0), default=1.0,
      help='Line extrusion ratio between 0.0 to 10.0.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc,
      Kf=args.Kf, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve,
      h=args.Lh, w=args.Lw, r=args.Lr)
  gen.dotests()
  gen.endFile()
  data=gen.getCode()
  sys.stdout.write(data)
