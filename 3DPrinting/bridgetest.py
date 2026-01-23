#!/bin/pypy3
"""
Printer Bridge Flow Calibration Tests

These tests assume an older printer with Vp=100mm/sec, Ap=500mm/sec^2, and a
0.4mm nozzle. This means it takes 0.2sec to reach full speed over 10mm of
distance.

Bridge lines are between 0.2@10-0.6@60, with areas 0.031-0.283mm^2,
flows 0.314-16.965mm^3/sec, and ve=0.1-7.1mm/s.

Normal lines are between 0.1x0.4@10-0.3x0.6@100 with areas 0.040-0.180mm^2,
flows 0.4-18.0mm^3/sec, and ve=0.17-7.5mm/sec.

If we assume PA has Kf=2.0, it will take about 6.3sec for flow to stabilize.
For @100mm/sec speeds thats 630mm of line, or 13 spans of a 50mm bridge, 7
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
  _CONF_ARGS = ('Kf','Kb','Re', 'fKf', 'dynret', 'dynext')

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
  rdy = 3  # test ruler y space.
  ldy = 5  # test line y distance.
  tdy = (tn+1)*ldy  # test runs y distance.
  cdy = sdy+rdy+1+tdy+2  # test case y distance.
  rulerbox = (tdx, rdy+1)   # test ruler box size including space.
  testsbox = (tdx, tdy+2)   # test run box size including space.
  tcasebox = (tdx, cdy+2)   # test case box size including space.
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
    for x,y in ((x0-w,y0), (x1+w,y0)):
      self.line(((x,y0),(x,y1),(x-w,y1), (x-w,y0)))

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
        self.hopup()

  def stabilize(self, x0, y0, dx=t0x, dy=ldy, v=20, w=None):
    """ Stabilize nozzle pressure by drawing filled rectangle."""
    if w is None: w = self.layer.w
    n = floor(dy / w)
    self.hopdn(x0, y0)
    self.cmt('TYPE:Bottom surface')
    for i in range(n):
      if i: self.draw(dy=-w, v=v, w=w)
      self.draw(dx=dx, v=v, w=w)
      dx*=-1
    self.hopup()

  def measure(self, x0, y0, de0=0.0, dx=1.0, de=-0.1, n=60,
              vx=10, ve=None, h=None):
    """Measure pressure required to start or stop flow.

    This works by hopping down, adjusting pressure by de0, then doing n steps
    of move dx, retract/restore de. This should leave a line where the length
    shows how many de's were required to start/stop extruding after the de0
    adjustment. The default arguments will do retract 5mm over 5cm of 1mm
    moves.
    """
    self.hopdn(x0, y0, h=h, de=de0, vb=ve)
    for i in range(n):
      self.move(dx=dx,v=vx)
      self.retract(de=de,ve=ve)
    self.hopup(de=0)

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
    ln, lz = self.layer.n, self.layer.z
    # Draw a brim line for the measure lines to attach to.
    x0 += sx
    self.hopdn(x0+dx-w, y0)
    self.cmt('TYPE:Brim')
    self.draw(dy=-dy)
    self.draw(dx=w)
    self.draw(dy=dy)
    self.hopup()
    for x,y in ((x0+5,y0),(x0+dx+(1-nw)*w-5, y0)):
      self.hopdn(x,y)
      self.cmt('TYPE:Support')
      dyi = -dy
      for l in range(nh):
        # before each layer except the first, increment the layer and hopdn at the next start.
        if l:
          #self.move(dz=h)
          self.nextLayer()
          self.hopdn(dx=(1-nw)*w, de=0)
        for i in range(nw):
          # before each line except the first, move to the next line.
          if i:
            self.draw(dx=w, dz=0, w=w, h=h)
          self.draw(dy=dyi, dz=0, w=w, h=h)
          dyi*=-1
      self.hopup()
      # Set the layer back to the starting layer/height.
      self.nextLayer(n=ln, z=lz)

  def bridgeline(self, x0, y0, dx=tex, dy=ldy, Bd=None, Bfw=1.0, Bfh=1.0, vl=None, ve=None, n=None, nh=bnh):
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
    vl,ve,h,w,r = self._getVlVehwr(vl=vl, ve=ve, h=Bh, w=Bs, r=Bh*Bs/Ba)
    ln,lz = self.layer.n, self.layer.z
    if not n:
      # default n for 5 secs of extrusion, up to as many bridges as will fit.
      n = min(dy/Bs, ceil(vl*5/dx))
    # for even number of crossings, start on the opposite end.
    if not n%2:
      x0, dx = x0+dx,-dx
    self.nextLayer(n=ln+nh+1,z=lz + nh*self.Lh, h=h)
    self.hopdn(x0,y0,h=h+Bz)
    self.cmt('TYPE:Bridge')
    for i in range(n):
      # Before each line except the first, move to the next line position.
      if i: self.draw(dy=-Bs, dz=0, v=vl, h=h, w=w, r=r)
      self.draw(dx=dx, dz=0, v=vl, h=h, w=w, r=r)
      dx*=-1
    # Set the layer back to the original layer.
    self.nextLayer(n=ln,z=lz)

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
    maxn = 140/self.tcasebox[1]
    assert len(tests) <= maxn, f'Can only fit max {maxn} tests.'
    # set alltests for running all tests and initialize n=0 if needed.
    alltests, n = n is None, n or 0
    tnx, tny = self.tcasebox[0],self.tcasebox[1]*maxn
    x0, y0 = -70, 70
    x1, y1 = x0+tnx, y0-tny
    self.preextrude(n=n)
    self.startLayer(Vp=10)
    # Only draw the title and brim if running test 0.
    if n==0:
      self.title(x0,y0, title=name, **kwargs)
      self.brim(x0,y0,x1, y1)
    y0 -= 5
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
    self.settings(x0, y0, title=name, sep=' ', **kwargs)
    y0 -= self.sdy
    self.ruler(x0, y0)
    y0 -= self.rdy
    if prepfn is not None:
      prepfn(x0,y0)
    self.cmt(f'TYPE:Inner wall')
    y = y0
    e = 0.00001
    while all(v <= v1+e for (v,v1,dv) in allsteps.values()):
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
      x += self.tmx
      self.stabilize(x, y)
      for k,(v,v1,dv) in allsteps.items():
        allsteps[k] = (v+dv,v1,dv)
      y-=self.ldy
    # Restore the config arguments after the test.
    #self.popconf()

  def _getLineArgs(self, dl=None, de=None, dt=None, vl=None, ve=None, h=None, w=None, r=1.0):
    """ Get dl, de, dt, vl, ve, h, w, and r from each other. """
    # Use default layer vl if there is not enough to derive it.
    if (vl,ve) == (None,None) and None in (dl,dt) and None in (de,dt):
      vl = self.layer.Vp if r else self.layer.Vt
    # Derive vl,dl, and dt from each other.
    if None not in (dl,dt):
      assert vl is None, "cannot specify vl with dl and dt"
      vl = dl/dt if dl else 0.0
    elif None not in (vl,dt):
      assert dl is None, "cannot specify dl with vl and dt"
      dl = vl*dt
    elif None not in (vl,dl):
      assert dt is None, " cannot specify dt with vl and dl"
      dt = dl/vl if dl else 0.0
    # Derive ve, de, and dt from each other.
    if None not in (de,dt):
      assert ve is None, "cannot specify ve with de and dt"
      ve = de/dt if de else 0.0
    elif None not in (ve,dt):
      assert de is None, "cannot specify de with ve and dt"
      de = ve*dt
    elif None not in (ve,de):
      assert dt is None, "cannot specify dt with ve and de"
      dt = de/ve if de else 0.0
      # Also set derived vl and dl from dt.
      if vl is None and dl is not None:
        vl = dl/dt
      elif dl is None and vl is not None:
        dl = vl*dt
    # Derive h,w,r,vl,ve from each other.
    if vl and ve is not None:
      assert None in (h,w,r), 'cannot specify all h,w,r with vl and ve'
      if h is None: h = self.h if None in (w,r) else ve*self.Fa/(vl*w*r)
      if w is None: w = self.w if None in (h,r) else ve*self.Fa/(vl*h*r)
      if r is None: r = ve*self.Fa/(vl*h*w)
    # Set h,w to defaults because we don't have enough to derive them.
    else:
      assert r is not None, 'require r if without ve and vl'
      if h is None: h = self.h
      if w is None: w = self.wl
    # Derive vl,de,dt from ve and dl.
    if vl is None:
      assert ve is not None, 'should have ve if without vl'
      assert r, 'require r!=0 to derive vl from ve'
      vl = ve*self.Fa/(h*w*r)
    if ve is None:
      assert vl is not None, 'should have vl if without ve'
      ve = vl*h*w*r/self.Fa
    assert abs(ve*self.Fa - vl*h*w*r) < 0.0001
    # Set dt from dl/vl or de/ve if we still don't have it.
    if dt is None:
      if vl and dl is not None:
        dt=dl/vl
      elif ve and de is not None:
        dt=de/ve
    # Set dl and de from dt if we still don't have them.
    if dt is not None:
      if dl is None:
        dl = vl*dt
      if de is None:
        de = ve*dt
    assert (dl,de,dt) == (None, None, None) or abs(de*self.Fa - dl*h*w*r) < 0.0001
    return dl,de,dt,vl,ve,h,w,r

  def _getVlVehwr(self, vl=None, ve=None, h=None, w=None, r=1.0):
    """ Get vl, ve, h, w, and r from each other. """
    assert any((ve,vl)), 'must specify at least ve or vl'
    if vl is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out vl from ve, w, h, and r.
      vl = ve*self.Fa/(h*w*r)
    if ve is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out ve from vl, w, h, and r.
      ve = h*w*r*vl/self.Fa
    # Set h and w
    if h is None: h = self.layer.h if w is None else ve*self.Fa/(w*r*vl)
    if w is None: w = ve*self.Fa/(h*r*vl)
    assert abs(ve*self.Fa - h*w*r*vl) < 0.0001
    return vl,ve,h,w,r


if __name__ == '__main__':
  import sys, argparse

  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  gcodegen.GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=gcodegen.RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(**gcodegen.GCodeGetArgs(args))

  bridgeargs=dict(name="BrPA", linefn=gen.bridgeline, prepfn=gen.bridgeprep, Kf=gen.Kf, Kb=gen.Kb, Re=gen.Re, tests=(
    dict(ve=(0.1, 2.6), vl=10),
    dict(ve=(2.5, 7.5), vl=40),))

  gen.startFile()
  gen.doTests(**bridgeargs)
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
