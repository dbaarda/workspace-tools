#!/bin/pypy3
"""
Printer Bridge Flow Calibration Tests

These tests assume an older printer with Vp=100mm/sec, Ap=500mm/sec^2, and a
0.4mm nozzle. This means it takes 0.2sec to reach full speed over 10mm of
distance.

Bridge lines are between D0.2@10-D0.6@60, with areas 0.031-0.283mm^2,
flows 0.314-16.965mm^3/sec, and ve=0.1-7.1mm/s.

Normal lines are between 0.1x0.4@10-0.3x0.6@100 with areas 0.040-0.180mm^2,
flows 0.4-18.0mm^3/sec, and ve=0.17-7.5mm/sec.

If we assume PA has Kf=2.0, it will take about 6.3sec for flow to stabilize.
For @100mm/sec speeds that's 630mm of line, or 13 spans of a 50mm bridge, 7
traversals of a 100mm line, or 1.7 times around a 12cm diameter spiral.

We need a warmup line to restore and stabilize the nozzle pressure to a
repeatable sane state before each test. If we assume the worse possible
over-retraction remaining from the previous test is 5mm, we need to extrude
5mm and then draw a continuous line for 6.3 secs before doing a retract and
hopup to be ready to hopdon/restore and start the next test. Drawing a 0.2x0.5
line needs to be 24mm long per mm of fillament used, or 12cm for 5mm. Drawing
at 10mm/sec will take 6secs or 6cm to stabilize for Kf=2.0. So worst case we
need to draw at least 18cm of line before each test. We can do 2 12cm lines
next to each other. If we do 11 tests that will require 120x11mm of space.

"""
import argparse
import gcodegen
from pprint import *
from math import *

def roundup(v,m=1):
  return m*ceil(v/m)


class GCodeCfgMixin(gcodegen.GCodeGenBase):
  """GCode mixin to dynamically changing config settings."""

  # These are printer configuration test arguments.
  _CONF_ARGS = ('Kf','Kb','Re', 'fKf', 'Zh', 'dynret', 'dynext')

  def resetfile(self,*args,**kwargs):
    """Initialize confstack whenever we resetfile()."""
    super().resetfile(*args,**kwargs)
    self.confstack = []

  def _fval(self,k,v):
    """ Format a single test argument value or range."""
    v=f'{v[0]}-{v[1]}' if isinstance(v, tuple) else f'{round(v,4)}' if isinstance(v,(float,int)) else f'{v.__name__}'
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

  def _incstep(self, v0, v1, dv):
    return v0+dv, v1, dv

  def _endstep(self, v0, v1, dv):
    e = 0.0001
    return (dv>0 and v0>v1+e) or (dv<0 and v0<v1-e)

  def setconf(self, **kwargs):
    # Filter out configs set to None or already set.
    kwargs = {k:v for k,v in kwargs.items() if k in self._CONF_ARGS and v is not None and getattr(self,k) != v}
    if kwargs:
      self.log(f'Setting config {self._fset(**kwargs)}')
    self.add(kwargs)

  def pushconf(self, **kwargs):
    """Push and change config attributes on self."""
    old_conf = {k:getattr(self,k) for k in self._CONF_ARGS}
    self.confstack.append(old_conf)
    self.log(f'Pushed config {self._fset(**old_conf)}')
    self.setconf(**kwargs)

  def popconf(self):
    """Pop and restore config attributes on self."""
    conf = self.confstack.pop()
    self.log(f'Popped config {self._fset(**conf)}')
    self.setconf(**conf)

  def isconf(self, code):
    return isinstance(code, dict)

  def incconf(self, conf):
    """Apply a config change to the current state."""
    for k, v in conf.items():
      setattr(self, k, v)

  def fconf(self, conf):
    """Format a config change."""
    # Test config args only emit fKf settings.
    fKf = conf.get('fKf')
    if fKf is not None:
      return self.fcmd('M900', K=f'{fKf:.3f}', T=0)

  def fmt(self, code):
    if self.isconf(code):
      return self.fconf(code)
    else:
      return super().fmt(code)

  def inc(self, code):
    if self.isconf(code):
      self.incconf(code)
    else:
      super().inc(code)


class ExtrudeTest(GCodeCfgMixin, gcodegen.GCodeGen):
  """ Generate GCode Extrusion Tests.

  Note that nozzle vs fillament diameters and/or track areas mean that de=1mm
  translates to nearly 20mm of nozzle thread, or possibly even more of track.
  This means 1mm of uncompensated pressure advance translates to more 20mm of
  track smear or stringing.
  """
  #Zh = 1.2 # set Zh high enough to clear bridge-tests.
  vx0 = 5  # test min speed
  vx1 = 100  # test max speed
  le0 = 0  # test min retraction length
  le1 = 10  # test max retraction length
  ve0 = 15  # test min retraction speed
  ve1 = 45  # test max retraction speed
  t0x = 10  # test line warmup x distance.
  t1x = 10  # test line cooldn x distance.
  tex = 60  # test line execute x distance.
  tmx = 60  # test line measurement x distance.
  tdx = t0x+tex+tmx+t1x    #  total test line x distance.
  tn = 10  # number of increments for each test (tn+1 test lines).
  sdy = 5  # test settings y space.
  rdy = 4  # test ruler y space.
  ldy = 5  # test line y distance.
  tdy = (tn+1)*ldy  # test runs y distance.
  cdy = sdy+rdy+tdy+2  # test case y distance.
  rulerbox = (tdx, rdy)   # test ruler box size including space.
  testsbox = (tdx, tdy)   # test run box size including space.
  tcasebox = (tdx, cdy)   # test case box size including space.
  rtics = (0, t0x, (t0x+tex, tmx), tdx) # Ruler ticks.
  bnh = 4 # number of layers for bridge test supports.
  bnw = 2 # number of lines for bridge test supports.

  def preextrude(self,n=0,x0=-40,y0=-74,x1=40,y1=-74,v=20):
    self.cmt('TYPE:Custom')
    self.line(l=((x0,y0-n),(x1,y0-n)), v=v, r=2)

  def title(self,x0,y0, title='', **kwargs):
    self.settings(x0,y0,sep=' ', title=title, **kwargs)

  def settings(self, x0, y0, sep=' ', title='', **kwargs):
    """ Draw the test arguments at (x0,y0)."""
    if title: title += sep
    self.text(title + self._fset(sep=sep, **kwargs), x0, y0, fsize=5)

  def brim(self, x0=-70, y0=70, x1=70, y1=-70):
    w = self.layer.w
    self.cmt("TYPE:Brim")
    for x in (x0-2*w, x1+w):
      self.line(((x,y0),(x,y1),(x+w,y1), (x+w,y0)))

  def ruler(self, x0, y0, box=rulerbox, tics=rtics):
    """ Draw a ruler that fits in box with tics.

    ticks is a sequence of either single pos values (for major tics), or a
    (pos,len) tuple which will draw a major tics every 1cm and minor ticks
    every 2mm for length len.
    """
    ml = box[1] - 1 # marker length
    self.cmt('TYPE:Outer wall')
    # draw ruler
    self.line([(x0, y0),(x0+box[0],y0)])
    for tic in tics:
      p,d = tic if isinstance(tic, tuple) else (tic,0)
      self.move(x=x0+p-2, y=y0)
      for d in range(0,d + 1,2):
        l = ml/2 if d % 10 else ml
        self.hopdn(dx=2,y=y0)
        self.draw(dy=-l)
        self.hopup(dz=0.1)

  def stabilize(self, x0, y0, dx=t0x, dy=ldy, v=20, w=None):
    """ Stabilize nozzle pressure by drawing filled rectangle."""
    if w is None: w = self.layer.w
    x0,dx = x0 + w, dx - 3*w
    n = floor(dy / w)
    self.hopdn(x0, y0)
    self.cmt('TYPE:Bottom surface')
    for i in range(n):
      if i: self.draw(dy=-w, v=v, w=w)
      self.draw(dx=dx, v=v, w=w)
      dx = -dx
    self.hopup()

  def grid(self, x0, y0, dx=tmx, dy=-tdy, gx=2, gy=1, v=None, h=0.2, w=0.5):
    """ Draws a grid at `(x0,y0)` of size `(dx,dy)` with lines every gx,gy.

    The grid starts at `(x0,y0)` and extends in the direction of `(dx,dy)`, so
    inverted and reversed grids can be drawn with negative dx and/or dy.
    """
    # get number of x and y lines and adjust dx and dy to fit.
    gx,gy = copysign(gx,dx), copysign(gy,dy)
    nx, ny = floor(dx/gx), floor(dy/gy)
    dx,dy = nx*gx, ny*gy
    x1,y1 = x0+dx, y0+dy
    lx0,lx1 = x0,x1
    for iy in range(ny+1):
      y = y0 + iy*gy
      self.line([(lx0,y),(lx1,y)], v=v, h=h, w=w)
      lx0, lx1 = lx1, lx0
    ly0, ly1 = y0, y1
    for ix in range(nx+1):
      x = x0 + ix*gx
      ty = 0 if ix % 5 else copysign(1, ly1-ly0)
      self.line([(x,ly0-ty), (x, ly1+ty)], v=v, h=h, w=w)
      ly0, ly1 = ly1, ly0

  def mgrid(self, x0, y0, dx=tmx, dy=-tdy, ly=ldy, tx=2, v=None, h=0.2, w=0.5):
    """ Draws a pressure measurement grid at `(x0,y0)` of size `(dx,dy)`.

    There are vertical lines every tx, and horizontal lines with a
    gap of `w` between them every `ly`.
    """
    # get number of x and y lines and adjust dx and dy to fit.
    assert ly >= 3*w, "grid test lines packed too dense"
    ly,gy = copysign(ly,dy), copysign(2*w, dy)
    ny = floor(dy/ly)
    x1,y1 = x0+dx, y0+dy
    self.cmt('TYPE:Outer wall')
    self.line([(x0-w,y0), dict(dy=dy)], h=h, w=w)
    for iy in range(ny):
      y = y0 + iy*ly
      self.grid(x0,y,dx,gy, v=v, h=h, w=w)

  def testprep(self, x0, y0, dx=tdx, dy=tdy, tl=(0,t0x-0.5, t0x+tex-0.5,t0x+tex+tmx+0.5, tdx), mx0=t0x+tex, mdx=tmx, h=0.2, w=0.5):
    self.Zh=0.1
    # draw two horizontal lines just outside test area.
    fx0,fy0,fx1,fy1 = x0,y0+w,x0+dx, y0-dy-w
    self.line([(fx0,fy0), (fx1,fy0)], h=h, w=w)
    self.line([(fx1,fy1), (fx0,fy1)], h=h, w=w)
    # draw vertical lines at each tl offset.
    for x in tl:
      fx = fx0 + x
      self.line([(fx,fy0),(fx,fy1)], h=h, w=w)
      fy0,fy1=fy1,fy0
    # draw the measurement grid at the mx0 offset.
    self.mgrid(x0+mx0,y0,dx=mdx,dy=-dy)

  def measure(self, x0, y0, de0=0.0, dx=tmx, de=-0.1*tmx, vx=10, ve=None, h=0.2, w=0.5, tx=2):
    """Measure pressure required to start or stop flow.

    This works by adjusting pressure by `de0@ve`, hopping down to `(x0,y0)`,
    then moving `dx@vx` while applying `de`. This should leave a line where
    the length shows how much of de was required to start/stop extruding after
    the de0 adjustment. The default arguments will retract `de=6mm` over
    `dx=6cm`, or 1mm every 1cm.

    At the end it restores enough that the whole measure, excluding any
    extruded fillament, is equivalent to retracting by `Re`.

    Note that a 0.2x0.5x1.0 line uses 1mm^3/cm. Assuming filament flow is
    limited by a standard measurement grid with 0.2x0.5x1.5mm channel space
    every 2mm, the extruded fillament should be limited to de=0.312mm/cm, or
    1.87mm over the whole 6cm measurement line.
    """
    # Adjust y0 to align with measurement grid.
    y0 -= w
    # To minimize drool when moving into position for measure, if de0 is
    # positive, do a normal hopdn() to apply pressure after move, otherwise
    # retract() before hopdn() to relieve pressure before move.
    if de0>=0:
      self.hopdn(x0, y0, h=h, de=de0, vb=ve)
    else:
      self.retract(de=de0, ve=ve)
      self.hopdn(x0, y0, h=h, de=0)
    self.cmt('WIPE_START')
    self.move(dx=dx, de=de, v=vx)
    self.cmt('WIPE_END')
    # Restore any retraction that exceeds a normal -Re retract.
    self.restore(de=-(self.Re+de0+de))

  def bridgeprep(self, x0, y0, sx=t0x, dx=tex, dy=tdy, nh=bnh, nw=bnw, h=None, w=None):
    """Setup Bridge Test Base

    This draws two walls nw lines wide starting at x offset sx that are dy
    long, dx-10 apart, and nh layers high. It also draws a brim line to attach
    the measure lines to, since the bridge line is not attached to the measure
    line.

    Note that this does increment layers, and ends hopped-up Zh above the last
    layer of the second wall, and then sets the next layer back down at the
    starting layer.
    """
    if h is None: h=self.layer.h
    if w is None: w=self.layer.w
    # Save the current layer details for restoring at the end.
    base = self.layer
    # Draw a brim line for the measure lines to attach to.
    x0 += sx
    self.hopdn(x0+dx-w, y0)
    self.cmt('TYPE:Brim')
    self.draw(dy=-dy)
    self.draw(dx=w)
    self.draw(dy=dy)
    self.hopup(dz=0.1)
    for x,y in ((x0+5,y0),(x0+dx+(1-nw)*w-5, y0)):
      self.hopdn(x,y)
      self.cmt('TYPE:Support')
      dyi = -dy
      for l in range(nh):
        # before each layer except the first, increment the layer and hopdn at
        # the next start.
        if l:
          self.nextLayer()
          self.hopdn(dx=(1-nw)*w, de=0)
        for i in range(nw):
          # before each line except the first, move to the next line.
          if i:
            self.draw(dx=w, dz=0, w=w, h=h)
          self.draw(dy=dyi, dz=0, w=w, h=h)
          dyi = -dyi
      self.hopup(dz=0.1)
      # Set the layer back to the starting layer/height.
      self.nextLayer(*base)
    self.Zh = nh * self.Lh

  def bridgeline(self, x0, y0, dx=tex, dy=ldy, Bd=None, Bfw=0.0, Bfh=1.0, vl=None, ve=None, n=None, nh=bnh):
    """ Draw a bridge starting at x0,y0 of length dx with n crossings on top of supports with nh layers. """
    assert (Bd,ve,vl).count(None) <= 1, "Must specify at least 2 of Bd, ve, and vl."
    # Take into account 5mm gaps before/after bridge supports.
    x0, dx = x0+5, dx-10
    if Bd:
      Ba = gcodegen.acirlcle(Bd)
    elif ve and vl:
      Ba = self.Fa * ve/vl
    else:
      Ba = self.Na
    if not Bd: Bd = gcodegen.dcircle(Ba)
    Bl = sqrt(Ba)
    Bs = Bd - Bfw*(Bd - Bl)
    Bh = Bd - Bfh*(Bd - Bl)
    Bz = (Bd - Bh)/2
    vl,ve,h,w,r = self.getVlVehwr(vl=vl, ve=ve, h=Bh, w=Bs, r=Ba/(Bh*Bs))
    base = self.layer
    if not n:
      # default n for 5 secs of extrusion, up to as many bridges as will fit.
      n = min(dy/Bs, ceil(vl*5/dx))
    # for even number of crossings, start on the opposite end.
    if not n%2:
      x0, dx = x0+dx,-dx
    self.nextLayer(n=base.n + nh + 1, z=base.z + nh*self.Lh, h=h)
    self.hopdn(x0, y0, h=h+Bz)
    self.cmt('TYPE:Bridge')
    for i in range(n):
      # Before each line except the first, move to the next line position.
      if i: self.draw(dy=-Bs, dz=0, v=vl, h=h, w=w, r=r)
      self.draw(dx=dx, dz=0, v=vl, h=h, w=w, r=r)
      dx = -dx
    # Set the layer back to the original layer.
    self.nextLayer(*base)

  ########
  # Grid Tests
  #

  # This is designed to test and calibrate the grid pressure measurement
  # system. The test is designed to drain any accumulated pressure and wipe
  # away any drool. It uses a grid similar to the measurement grid but does
  # multiple slow passes, which should ensure there is no pressure or drool
  # left at the end. The measurement phase can then pre-apply known de0
  # pressures that should be directly reflected in the measurement result.

  def gridprep(self, x0, y0, dx=tex, dy=tdy, h=0.2, w=0.5):
    """Setup Grid Test Base"""
    self.grid(x0,y0,dx-5,-dy, gx=5, gy=dy, h=h, w=w)

  def gridline(self, x0, y0, dx=tex, dy=ldy, dt=10, vx=10.0, h=0.2, w=0.5):
    # Adjust start and size to fit grid.
    x0, y0, dx, dy = x0, y0-w, dx-5, dy-w
    # calc number of line passes required for dt secs that will fit.
    n = min(floor(dy/w), ceil(dt*vx/dx))
    # Start from the other end if there are an even number of lines.
    x0, dx = (x0, dx) if (n % 2) else (x0+dx,-dx)
    self.hopdn(x0,y0,h=h, w=w)
    self.cmt('WIPE_START')
    for i in range(n):
      if i: self.move(dy=-w)
      self.move(dx=dx, v=vx)
      dx = -dx
    self.cmt('WIPE_END')

  def doTests(self, name, tests, linefn, prepfn=None, n=None, **kwargs):
    """Run up to as many tests as will fit.

    This calls doTest() with arguments from kwargs overridden by arguments for
    each test case in tests.

    Args:
      name: the text name of the tests
      linefn: The function that executes a single line of a test.
      prepfn: Optional function that prepares for the test.
      tests: A sequence of up to 4 dicts containing args for each test.
      n: optional index of a single test to run, otherwise run all.
      **kwargs: key-value args common to all tests and to put in the title.
    """
    maxn = floor((140-self.sdy)/self.tcasebox[1])
    assert len(tests) <= maxn, f'Can only fit max {maxn} tests.'
    # set alltests for running all tests and initialize n=0 if needed.
    alltests, n = (n is None, n or 0)
    tnx, tny = self.tcasebox[0], self.sdy + self.tcasebox[1]*maxn
    x0, y0 = -70, 70
    x1, y1 = x0+tnx, y0-tny
    self.preextrude(n=n)
    self.startLayer(Vp=10)
    # Only draw the title and brim if running test 0.
    if n==0: self.title(x0,y0, title=name, **kwargs)
    y0 -= self.sdy
    if n==0: self.brim(x0,y0,x1, y1)
    for i,t in enumerate(tests):
      # only run test n if not running all tests.
      if alltests or i == n:
        self.doTest(f'{name}:{i}', x0, y0, linefn, prepfn, runargs=kwargs, **t)
      y0 -= self.tcasebox[1]
    self.endLayer()

  def doTest(self, name, x0, y0, linefn, prepfn=None, runargs=None, **kwargs):
    allargs = (runargs or {}) | kwargs
    allsteps = {k:self._getstep(self.tn, v) for k,v in allargs.items()}
    assert any(dv for _,_,dv in allsteps.values()), 'At least one arg in {allargs} must be a range.'
    self.log(f'Test {name} {self._fset(**allargs)}')
    self.settings(x0, y0, sep=' ', **kwargs)
    y0 -= self.sdy
    #self.ruler(x0, y0)
    #y0 -= self.rdy
    self.testprep(x0,y0)
    if prepfn is not None:
      prepfn(x0+self.t0x, y0)
    self.cmt(f'TYPE:Inner wall')
    y = y0
    e = 0.00001
    while all(not self._endstep(*s) for s in allsteps.values()):
      de0,_,_ = allsteps.get('de0', (0.0,0.0,0.0))
      lineargs={k:v for k,(v,_,_) in allsteps.items() if k != 'de0'}
      confargs={k:v for k,v in lineargs.items() if k in self._CONF_ARGS}
      funcargs={k:v for k,v in lineargs.items() if k not in self._CONF_ARGS}
      self.log(f'Test Line {name} {self._fset(**lineargs)}')
      x = x0
      self.pushconf(**confargs)
      self.stabilize(x, y)
      x += self.t0x
      linefn(x, y, **funcargs)
      x += self.tex
      self.popconf()
      self.measure(x, y, de0=de0)
      x += self.tmx + self.layer.w
      self.stabilize(x, y)
      for k, s in allsteps.items():
        allsteps[k] = self._incstep(*s)
      y-=self.ldy
    # Restore the config arguments after the test.
    #self.popconf()


if __name__ == '__main__':
  import sys, argparse

  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  gcodegen.GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=gcodegen.RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(**gcodegen.GCodeGetArgs(args))

  gridargs=dict(name="GrCal", linefn=gen.gridline, prepfn=gen.gridprep, Kf=gen.Kf, Kb=gen.Kb, Re=gen.Re, tests=[
    dict(de0=(-0.2, 5.8), vx=5),
    dict(de0=(-0.2, 5.8), vx=60)])

  bridgeargs=dict(name="BrPA2", linefn=gen.bridgeline, prepfn=gen.bridgeprep, Kf=gen.Kf, Kb=gen.Kb, Re=gen.Re, tests=(
    dict(ve=(0.1, 2.6), vl=10),
    dict(ve=(2.5, 5.0), vl=60, de0=(0.0,-1.0)),))

  gen.startFile()
  gen.doTests(**gridargs)
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
