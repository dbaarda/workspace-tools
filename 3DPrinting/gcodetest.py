#!/bin/pypy3
"""
"""
import argparse
import gcodegen
from pprint import *
from math import *

def solvequad(a,b,c):
  """ Solve a quadratic equation, returning the two roots."""
  d = sqrt(b**2 - 4*a*c)
  return (-b - d)/(2*a), (-b + d)/(2*a)

def roundup(v,m=1):
  return m*ceil(v/m)

class ExtrudeTest(gcodegen.GCodeGen):
  """ Generate GCode Extrusion Tests.

  Note that nozzle vs fillament diameters and/or track areas mean that de=1mm
  translates to nearly 20mm of nozzle thread, or possibly even more of track.
  This means 1mm of uncompensated pressure advance translates to more 20mm of
  track smear or stringing.
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
  # Constants for spiral tests.
  R0 = 70
  R1 = 42
  C0 = 0.0
  dRdC = -2
  LN = 5

  # These are printer configuration test arguments.
  _CONF_ARGS = ('Kf','Kb','Re', 'fKf', 'dynret', 'dynext')

  def resetfile(self,*args,**kwargs):
    """Initialize confstack whenever we resetfile()."""
    super().resetfile(*args,**kwargs)
    self.confstack = []

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

  def setconf(self, **kwargs):
    # Filter out configs set to None or already set.
    kwargs = {k:v for k,v in kwargs.items() if v is not None and getattr(self,k) != v}
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

  def preextrude(self,n=0,x0=-50,y0=-60,x1=70,y1=60):
    self.preExt(x0-2*n,y0,x1,y1+2*n,m=5,lw=0)

  def title(self,x0,y0, **kwargs):
    self.settings(x0+2,y0+5,sep=' ', **kwargs)

  def brim(self, x0, y0, x1, y1):
    w = self.w
    self.cmt("structure:brim")
    self.hopdn(x0,y0+5)
    self.draw(y=y1-5)
    self.draw(dx=-w)
    self.draw(y=y0+5)
    self.hopup()
    self.hopdn(x1,y0+5)
    self.draw(y=y1-5)
    self.draw(dx=w)
    self.draw(y=y0+5)
    self.hopup()

  def ruler(self, x0, y0, rl=rdx, rs=t0x):
    """ Draw a ruler rl long and 3mm high with tics offset rs from each end."""
    ml,tl=3,1.5  # tick marker lengths
    self.cmt('structure:shell-outer')
    # draw ruler
    self.line([(x0, y0),(x0+rl,y0)])
    self.move(x=x0+rs-2, y=y0)
    for d in range(0,rl - 2*rs + 1,2):
      l = tl if d % 10 else ml
      self.hopdn(dx=2,y=y0)
      self.draw(dy=-l)
      self.hopup()

  def settings(self, x0, y0, sep='\n', title='', **kwargs):
    """ Draw the test arguments at (x0,y0)."""
    self.text(self._fset(sep=sep, **kwargs), x0, y0, fsize=5)

  def doTests(self, name, linefn, tests, n=None, **kwargs):
    """Run up to 4 tests.

    Args:
      name: the text name of the tests
      linefn: The function that executes a single line of a test.
      tests: A sequence of up to 4 dicts containing args for each test.
      n: optional index of a single test to run, otherwise run all.
      **kwargs: key-value args to put in the title.
    """
    assert len(tests) <= 4, 'Can only fit max 4 tests.'
    # set alltests for running all tests and initialize n=0 if needed.
    alltests, n = n is None, n or 0
    t4x, t4y = self.tdx, 4*self.tdy - 5
    x0, y0 = -self.rdx/2, 55
    x1, y1 = x0+t4x, y0-t4y
    self.preextrude(n=n)
    self.startLayer(Vp=10)
    # Only draw the title and brim if running test 0.
    if n==0:
      self.title(x0,y0,**kwargs)
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
    self.settings(x0 + self.rdx+1, y0, **kwargs)
    y = y0 - self.rdy
    # Push the config arguments without applying changes before the test.
    self.pushconf()
    self.cmt(f'structure:shell-inner')
    e=0.00001
    while all(v <= v1+e for (v,v1,dv) in testargs.values()):
      lineargs={k:v for k,(v,_,_) in testargs.items()}
      confargs={k:v for k,v in lineargs.items() if k in self._CONF_ARGS}
      funcargs={k:v for k,v in lineargs.items() if k not in self._CONF_ARGS}
      self.log(f'Test Line {name} {self._fset(**lineargs)}')
      self.setconf(**confargs)
      linefn(x0, y, **funcargs)
      for k,(v,v1,dv) in testargs.items():
        testargs[k] = (v+dv,v1,dv)
      y-=self.ty
    # Restore the config arguments after the test.
    self.popconf()

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
      if h is None: h = self.layer.h if None in (w,r) else ve*self.Fa/(vl*w*r)
      if w is None: w = self.layer.w if None in (h,r) else ve*self.Fa/(vl*h*r)
      if r is None: r = ve*self.Fa/(vl*h*w)
    # Set h,w to defaults because we don't have enough to derive them.
    else:
      assert r is not None, 'require r if without ve and vl'
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
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

  def rc2xy(self, R=None, C=None, dR=None, dC=None, x0=0.0, y0=0.0):
    """Get x,y from R,C or dR,dC.

    Note R and C default to the current R,C position with optional dR,dC offsets.
    """
    if dR is None: dR = 0.0
    if dC is None: dC = 0.0
    if R is None: R = dist((x0, y0), (self.x, self.y)) + dR
    if C is None: C = atan2(self.y-y0, self.x-x0)/(2*pi) + dC
    a = C*2*pi
    return x0 + R*cos(a), y0 + R*sin(a)

  def xy2rc(self, x=None, y=None, dx=None, dy=None, x0=0.0, y0=0.0):
    """Get R,C from x,y or dx, dy.

    Note x and y default to the current x,y position with optional dx,dy offsets.
    """
    if dx is None: dx = 0.0
    if dy is None: dy = 0.0
    if x is None: x = self.x + dx
    if y is None: y = self.y + dy
    return dist((x0, y0), (x, y)), atan2(y-y0, x-x0)/(2*pi)

  def drawrc(self, R=None, C=None, dR=None, dC=None, x0=0.0, y0=0.0, **kwargs):
    x,y=self.rc2xy(R,C,dR,dC,x0,y0)
    self.draw(x=x,y=y,**kwargs)

  def hopdnrc(self, R=None, C=None, dR=None, dC=None, x0=0.0, y0=0.0, **kwargs):
    x, y = self.rc2xy(R, C, dR, dC, x0, y0)
    self.hopdn(x=x,y=y,**kwargs)

  def spiral(self, R=None, C=None, dR=None, dC=None, dRdC=None, x0=0.0, y0=0.0,
      dl=None, de=None, dt=None, vl=None, ve=None, h=None, w=None, r=1.0, pe=None):
    """ Draw a spiral around x0,y0.

    This draws a spiral from the current R0,C0 position to a position R,C. The
    position R,C can be specified using any sufficient combination of
    R,C,dR,dC,dRdC,dl.
    """
    R0,C0 = self.xy2rc(x0=x0,y0=y0)
    if dl == 0:
      # This is explicitly not a draw or move.
      if de is not None:
        # dl=0, de=<value> means apply a specific retract/restore distance.
        self.move(de=de, v=ve)
      elif pe is not None:
        # dl=0, pe=<value> means do a pressure compensating restore.
        self.restore(pe=pe)
      elif h is not None:
        # dl=0, h=<value> means apply a generic hopup() or hopdn().
        self.hopup() if h > 0 else self.hopdn()
      elif r is not None:
        # dl=0, r=<value> means apply a generic retract() or restore().
        self.retract() if r<0 else self.restore()
      return R0, C0
    R1,C1 = R,C
    if R1 is not None: dR = R1 - R0
    if C1 is not None: dC = C1 - C0
    if dR is None and None not in (dC, dRdC): dR = dC*dRdC
    if dC is None and dR is not None and dRdC: dC = dR/dRdC
    if dR is not None and dC is not None: dRdC = dR/dC
    # Get dl from dR and dC if provided.
    if None not in (dR,dC):
      assert dl is None, 'cannot specify dl with dR,dC'
      dl = dC*pi*(2*R0 + dR)
    dl, de, dt, vl, ve, h, w, r = self._getLineArgs(dl,de,dt,vl,ve,h,w,r)
    # Get dC from dl.
    if dC is None:
      assert dl is not None, 'need to specify enough for dl without dC'
      assert dRdC is not None, 'need to specify dRdC with dl'
      # Solve dRdC*dC^2 + 2*R0*dC - dl/pi = 0 and take the smaller result.
      dC = min(solvequad(a=dRdC, b=2*R0, c=-dl/pi))
    if dR is None: dR = dC*dRdC
    self.log(f'spiral RC0=({R0:.1f},{C0:.3f}) dRC=({dR:.1f},{dC:.3f}) {dl=:.1f}@{vl:.0f} {de=:.3f}@{ve:.1f} l={h:.2f}x{w:.2f}x{r:.2f}')
    R, C, C1, dC = R0, C0, C0+dC, 1/72
    while C < C1:
      dC = min(dC, C1 - C)
      C += dC
      R += dC*dRdC
      self.drawrc(R, C, x0=x0, y0=y0, v=vl, w=w, r=r)

  def dial(self, R=R0, x0=0.0, y0=0.0, tl=3):
    self.hopdnrc(R, 0, x0=x0, y0=y0)
    self.spiral(R, 1, x0=x0, y0=y0)
    self.hopup()
    for i in range(0,72):
      C = i/(72)
      self.hopdnrc(R,C,x0=x0,y0=y0)
      l = tl/2 if i % 6 else tl
      self.drawrc(dR=l,x0=x0,y0=y0)
      self.hopup()

  def dialStart(self,C, R0=R0, R1=R1):
    self.cmt('structure:shell-outer')
    self.hopdnrc(R0,C)
    self.drawrc(R1)
    self.hopup()

  def dialRuler(self, R0=R0, R1=R1):
    self.cmt('structure:shell-outer')
    self.dial(R0,tl=3)
    self.dial(R1,tl=-3)

  def doSpiralTests(self, name, lines, tests, ln=LN, **kwargs):
    """ Do tests in a spiral pattern instead of a grid for longer path lengths.

    tests is a list of test kwargs for lines, and  'lines' is a
    list of spiral kwargs. The spiral kwargs can have string name values
    referring to values in the test kwargs.
    """
    # Note for Kf=1.0,Kb=0.0 the decay to steady state takes about 3sec, which
    # for v=100mm/s is 300mm of line! A circle of D=100mm have about 300mm of
    # circumference, so a whole circuit of the spiral is needed. However for
    # v=25mm only 1/4 of a circle is needed. So we should either adjust the
    # diameter for the test print speed, or adjust the fraction of the circle
    # for each phase based on speed.
    # PreExtrude
    R0,R1,C0,dRdC = self.R0, self.R1, self.C0, self.dRdC
    self.hopdnrc(R=R0+5, C=1/16, h=0.2)
    self.spiral(R=R0+5,C=3/16,vl=10, w=0.4, r=4)
    self.hopup()
    # Prepare tests
    self.startLayer(Vp=10)
    self.settings(-R0,R0,**kwargs)
    self.dialRuler(R0,R1)
    R0+=self.dRdC; R1-=self.dRdC
    # Do tests
    Rl,Cl = R0,C0
    for tn,t in enumerate(tests):
      self.settings([0,-30,-30,0][tn],[25,25,0,0][tn], **t)
      # Draw the dial start if this is the outer radius test.
      if Rl == R0:
        self.dialStart(Cl)
      self.doSpiralTest(tn, name, Rl, Cl, lines, ln, **(kwargs|t))
      # Get the R,C position of the end of the last test and updated it for the next test.
      Rf,Cf = self.xy2rc()
      if Rl - Rf < Rf - R1:
        # Start next test on the next inside track.
        Rl = roundup(Rf, self.dRdC)
      else:
        # Start next test on the next quarter.
        Rl,Cl = R0, roundup(Cf,1/4)
        if Cl >= 1: break
    self.endLayer()

  def doSpiralTest(self, tn, name, R, C, lines=None, ln=5, **kwargs):
    """ Do a single spiral test.

    This does a sequence of spiral test lines with varying arguments. Each
    test line starts at the same C angle and changes R by self.dRdC for each track.

    Args:
      tn: the number of this test.
      name: the name of this test.
      R,C: the R,C coordinates where this test starts.
      dRdC: the change in R per rotation.
      lines: a list of spiral kwargs where some values can be string names refering to (possibly range-valued) kwargs arguments.
      ln: the number of test lines in this test.
      **kwargs: the (possibly range-valued) arguments for this test.
    """
    R1 = self.R1 - self.dRdC
    warmup = [
      dict(dC=1/72, vl=5, h=self.layer.h, w=self.layer.w, r=1.0),
      dict(dC=3/72, vl=1, h=self.layer.h, w=self.layer.w, r=0.0)]
    cooldn = [
      dict(dl=0, pe=2.5),
      dict(dC=3/72, vl=5,  h=self.layer.h, w=self.layer.w, r=1.0)]
    lines = warmup + lines + cooldn
    testargs = {k:self._getstep(ln, v) for k,v in kwargs.items()}
    assert any(dv for _,_,dv in testargs.values()), 'At least one arg in {kwargs} must be a range.'
    self.log(f'Test{tn} {name} {self._fset(**kwargs)}')
    # Adjust the retraction from the last hopup to what it would be for a
    # de=-5.0 retraction from the last pre-test draw. This is the retraction
    # performed after each line cooldown, which we try to match for the start
    # of the first line. For vl=5mm/s ve=0.5mm/s l=0.3x0.5x1.0 and Kf=1.6,
    # Kb=3.2 we get pe=2.5, so a de=-0.5 retract should have pe=-2.5.
    self.retract(pe=-2.5)
    # Push the config arguments without applying changes before the test.
    self.pushconf()
    self.cmt(f'structure:shell-inner')
    e=0.00001
    while all(v <= v1+e for (v,v1,dv) in testargs.values()) and R>R1:
      lineargs={k:v for k,(v,_,_) in testargs.items()}
      confargs={k:v for k,v in lineargs.items() if k in self._CONF_ARGS}
      self.log(f'Test Line {name} {self._fset(**lineargs)}')
      self.setconf(**confargs)
      self.hopdnrc(R,C, h=lineargs.get('h',self.layer.h), de=5.0)
      for l in lines:
        # Replace any name values in l with the corresponding value in lineargs.
        sargs={k:(eval(v,locals=lineargs) if isinstance(v,str) else v) for k,v in l.items()}
        self.spiral(dRdC=self.dRdC,**sargs)
        R, _ = self.xy2rc()
        if R <= R1: break
      self.hopup(de=-5.0)
      # Get the R,C position of the start of the next 1/4 spiral
      R = roundup(R,self.dRdC)
      # Update the testargs for the next lines.
      for k,(v,v1,dv) in testargs.items():
        testargs[k] = (v+dv,v1,dv)
    # Restore the config arguments after the test.
    self.popconf()

  def _getVxVehwr(self, vx=None, ve=None, h=None, w=None, r=1.0):
    """ Get vx, ve, h, w, and r from each other. """
    assert ve is not None or vx is not None, 'must specify at least ve or vx'
    if vx is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out vx from ve, w, h, and r.
      vx = ve*self.Fa/(h*w*r)
    if ve is None:
      if h is None: h = self.layer.h
      if w is None: w = self.layer.w
      # figure out ve from vx, w, h, and r.
      ve = h*w*r*vx/self.Fa
    # Set h and w
    if h is None: h = self.layer.h if w is None else ve*self.Fa/(w*r*vx)
    if w is None: w = ve*self.Fa/(h*r*vx)
    assert abs(ve*self.Fa - h*w*r*vx) < 0.0001
    return vx,ve,h,w,r

  def begLine(self, x0, y0, h):
    """ This is a standard warmup for a test line. """
    # Go to the start of the line and restore for the warmup.
    self.hopdn(x0, y0, h=h)
    # This slow draw then move is to ensure the nozzle is primed and wiped, and
    # any advance pressure, actual and estimated, has decayed away.
    self.draw(dx=self.t0x, v=self.vx0, w=self.layer.w, r=0.5)
    self.move(dx=10,v=1)

  def endLine(self, de=None):
    """ This is a standard cooldown for a test line."""
    self.restore(de=None)
    self.draw(dx=self.t0x, v=self.vx0, w=self.layer.w, r=1.0)
    self.hopup()

  def testRetract(self, x0, y0, vx=None, ve=None, h=None, w=None, r=1.0, vr=10, re=4):
    """ Test retract and restore for different settings.

    This test is to try and figure out how much retraction is needed to stop
    extruding after moving or figure out how retract/restore speed matters.
    The pre/post phases are;

      * 5mm: warmup draw at v=5mm/s,r=0.5 to prime the nozzle.
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
    vx,ve,h,w,r = self._getVxVehwr(vx,ve,h,w,r)
    # A printer with a=500mm/s^2 can go from 0 to 100mm/s in 0.2s over 10mm.
    # Retracting de=10mm at vr=100mm/s over 40mm means retracting at ve=25mm/s
    # which is probably fine.
    self.begLine(x0, y0, h)
    # This restore pre-loads advance pressure for drawing the line. If the
    # start of the line is too thin, there was not enough pressure, and Kf
    # probably needs to be increased. If the start is too thick, pressure was
    # too high and Kf should probably be dropped. The line should be
    # consistently the same width as all the other lines.
    self.restore()
    self.draw(dx=30,v=vx,w=w,r=r)      # draw at v=vx.
    self.move(dx=40,v=vr,de=-re)       # retract while moving at v=vr.
    self.endLine(de=re)

  def testStartStop(self, x0, y0, vx=None, ve=None, h=None, w=None, r=1.0):
    """ Test starting and stopping extrusion for different settings.

    This test is to check dynamic retracting for different Kf/Kb/Re
    settings.

    Test phases are;

      * 0mm: default drop and restore.
      * 5mm: warmup draw at v=5mm/s,r=0.5 to prime nozzle with fillament.
      * 10mm: move at 1mm/s to drain actual and estimated pressure.
      * 0mm: default restore to pre-apply pressure before draw.
      * 50mm: draw at <vx>mm/s to check pre-applied pressure and build pressure.
      * 0mm: default retract to relieve pressure.
      * 20mm: move at 10mm/s to check and drain any vestigial pressure.
      * 0mm: default restore.
      * 5mm: draw at v=5mm/s,r=1.0 to check slow draw after retract/restore.
      * 0mm: default retract and raise.

    Args:
      vx: draw speed to test.
      ve: optional extrude speed while drawing.
      h: optional line height.
      w: optional line width.
      r: optional extrusion ratio.
    """
    vx,ve,h,w,r = self._getVxVehwr(vx,ve,h,w,r)
    # Set different line-phase lengths.
    dx = (5, 10, 50, 20, 5)
    self.begLine(x0,y0,h)
    # This restore pre-loads advance pressure for drawing the line. If the
    # start of the line is too thin, there was not enough pressure, and Kf
    # probably needs to be increased. If the start is too thick, pressure was
    # too high and Kf should probably be dropped. The line should be
    # consistently the same width as all the other lines.
    self.restore()
    self.draw(dx=dx[2],v=vx,w=w,r=r)
    # This retract and slow move is to relieve the advance pressure, and see
    # if any remaining pressure drools off. If there is any trailing drool,
    # the amount of pressure was underestimated and Kf or Kb probably needs
    # to be increased. If there is still stringing when otherwise Kf/Kb seems
    # right, adding some additional retraction with Re could help.
    self.retract()
    self.move(dx=dx[3],v=10)
    # This restore applies advance pressure for a final slow line after
    # the earlier retraction. It should be the same width as all the other
    # lines. If it starts too thin, it suggests the earlier retraction
    # over-estimated the pressure and over-retracted, so Kf should be reduced.
    # If it starts too thick it suggests the earlier retraction under
    # estimated the pressure and extra retraction was actually relieving
    # pressure, so Kf should be increased and Re could possibly be reduced.
    self.endLine()

  def testBacklash(self, x0, y0, vx=None, ve=None, h=None, w=None, r=1.0, r0=0.0, re=5.0):
    """ Test retract and restore distances.

    This test is to test retraction and restore distances for different vx and/or ve
    velocities. Note draw speeds can be set with any combination of vx and/or
    ve optionally together with h, w, and r.

    Test phases are;

      * 0mm: hopdn and restore to h=0.2.
      * 5mm: draw at 5mm/s to prime nozzle.
      * 10mm: move at 1mm/s to drain nozzle.
      * 0mm: hopdn and restore to h=h.
      * 30mm: draw at <vx>mm/s to stabilize pressure at that speed.
      * 0mm: move to h=0.2mm for start of pressure measurement.
      * 0mm: retract initial de=re0.
      * 20mm: moving retract of <re>mm at 10mm/s to measure required retract.
      * 0mm: restore initial de=re0.
      * 20mm: moving restore of <re>mm at 10mm/s to measure required restore.
      * 0mm: default retract and hopup, then hopdn and restore to h=0.2mm to remove backlash.
      * 5mm: draw at 5mm/s to finalize line and stabilize pressure.
      * 0mm: default retract and hopup to relieve any vestigial pressure.

    Args:
      vx: movement speed to check.
      ve: extrusion speed to check.
      h,w,r: line height, width, and ratio to check.
      r0: initial retract/restore distance.
      re: moving retract/restore distance.
    """
    # For Kf=0.5,Kb=4.0 we expect values to be in the ranges;
    # Ve = 0.01 -> 6.0, or 1.0 for v=30 l=0.2x0.4x1.0.
    # Db = 0.1 -> 2.0, or 0.5 for w=0.5.
    # Pn = 0.005 -> 3.0 or 0.5 for ve=1.0.
    # Pb = 0.4 -> 8.0 or 2.0 for w=0.5.
    # Pe = 0.4 -> 11.0 or 2.5 for ve=1,w=0.5, or 4.5 for ve=1,w=1.0.
    # le = 0.008 -> 0.25 mm/mm (l=0.2x0.1x1.0 -> 0.3x2.0x1.0) or 0.05 for l=0.2x0.6x1.0.
    # For v=10mm/s, we need at least ve=0.2*0.1*10/Fa=0.08mm/s for a w=0.1mm line, or Pe=0.44mm for Kf=0.5,Kb=4.0.
    vx,ve,h,w,r = self._getVxVehwr(vx,ve,h,w,r)
    self.begLine(x0, y0, h=0.2)
    self.hopdn(h=h)
    self.draw(dx=30,v=vx,w=w,r=r)
    self.move(h=0.2)
    self.retract(de=r0)
    self.move(dx=20,de=-re, v=10, w=self.layer.w)
    self.restore(de=r0)
    self.move(dx=20,de=+re, v=10, w=self.layer.w)
    self.hopup()
    self.hopdn(h=0.2)
    self.endLine()

  def testKf(self, x0, y0, vx0, vx1):
    """ Test linear advance changing speed different speeds.

    This is to test that dynamic extrusion for different Kf/Kb/Re settings for
    different speed changes is right.

    The test phases are;

      * 5mm: warmup draw at 5mm/s
      * 20mm: draw at <vx0>mm/s base speed before changing speed.
      * 40mm: draw at <vx1>mm/s to check changing to and from this speed.
      * 20mm: draw at <vx0>mm/s base speed after changing speed.
      * 5mm: cooldown draw at 5mm/s

    Args:
      vx0: the slow speed to use.
      vx1: the fast speed to use.
    """
    self.log(f'testKfline {vx0=} {vx1=}')
    self.hopdn(x0, y0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.draw(dx=self.t1x,v=vx0)
    self.draw(dx=2*self.t1x,v=vx1)
    self.draw(dx=self.t1x,v=vx0)
    self.draw(dx=self.t0x,v=self.vx0)
    self.hopup()


class TextTest(gcodegen.GCodeGen):

  def texttest(self,h=5,x0=-70,y0=70,x1=70,y1=None):
    chars=''.join(c for c in vtext.glyphs)
    self.text(chars,x0,y0,x1,y1,fsize=h)


if __name__ == '__main__':
  import sys, argparse

  cmdline = argparse.ArgumentParser(description='Generate test gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  gcodegen.GCodeGenArgs(cmdline)
  cmdline.add_argument('-n', type=gcodegen.RangeType(0,4), default=None,
      help='Test number when running individual tests.')
  args=cmdline.parse_args()

  gen=ExtrudeTest(**gcodegen.GCodeGetArgs(args))

  backpressureargs=dict(name="Backpressure", linefn=gen.testRetract, tests=(
    dict(w=(0.3, 0.8), vx=10, ve=0.3*0.3*10/gen.Fa),
    dict(h=(0.1,0.35), ve=0.5*0.35*10/gen.Fa, w=0.5),
    dict(w=(0.3, 0.8), ve=0.2*0.8*10/gen.Fa, h=0.2),
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

  backlashargs=dict(name="Backlash", linefn=gen.testBacklash,
    Kf=args.Kf, Kb=args.Kb, Re=args.Re,
    # Note ve=1.0mm/s is the same as vx=30mm/s with h=0.2,w=0.4
    # or vx=10mm/s with h=0.3,w=0.8.
    tests=(
      dict(ve=1.0, h=0.3, w=(0.3,0.8), r0=0.0, re=4.0),
      dict(ve=1.0, h=0.3, w=(0.8,1.6), r0=2.0, re=4.0),
      dict(ve=1.0, h=0.2, w=(0.3,0.8), r0=0.0, re=4.0),
      dict(ve=1.0, h=0.2, w=(0.8,1.6), r0=2.0, re=4.0),
      ))

  Kf,Kb,Re = 0.4,1.0,1.0  # "best" default values.
  Kfargs=dict(name="Kf", linefn=gen.testKf, tests=(
    dict(Kf=(0.0,1.0), Kb=Kb, Re=Re, vx0=20, vx1=100),
    dict(Kf=Kf, Kb=(0.0,4.0), Re=Re, vx0=20, vx1=100),
    dict(Kf=Kf, Kb=Kb, Re=(0.0,5.0), vx0=20, vx1=100),
    dict(Kf=Kf, Kb=Kb, Re=Re, vx0=(10,50), vx1=(20,100))))

  # Pre Marlin v1.5 Kf was in extruder advance steps per mm/sec, with numbers
  # for PLA in the range 30-130. In v1.5 it was changed to mm of advance per
  # mm/s, with numbers in the range 0.1=2.0. We don't know what the
  # Adventurer3 supports. It's also possible setting pressure advance works
  # the same as turning the extruder fan on/off and it's executed a
  # parse-time, not pushed into the execution queue. This would mean it can
  # only realistically be set once at the start of a print, and cannot
  # reliably be adjusted mid-print.
  fKfargs=dict(name="fKf", linefn=gen.testKf,
    Kf=args.Kf, Kb=args.Kb, Re=args.Re,
    tests=(
      dict(fKf=(0.0,2.0), vx0=20, vx1=100),
      dict(fKf=0.0, vx0=(10,50), vx1=(20,100)),
      dict(fKf=1.0, vx0=(10,50), vx1=(20,100)),
      dict(fKf=2.0, vx0=(10,50), vx1=(20,100))))

  # Use the StartStopTest to try and characterise backpressure.
  backpressure1args=dict(name="Backpressure1", linefn=gen.testStartStop,
    Kf=args.Kf, Kb=args.Kb, Re=args.Re,
    tests=(
      dict(Kb=(0.0,4.0), Re=0.5, ve=1.0, h=0.2, w=0.3),
      dict(Kb=(0.0,4.0), Re=0.5, ve=1.0, h=0.2, w=0.4),
      dict(Kb=(0.0,4.0), Re=0.5, ve=1.0, h=0.2, w=0.8),
      dict(Kb=(0.0,4.0), Re=0.5, ve=1.0, h=0.2, w=1.6)))

  # We want to figure out the ideal restore and retract distances for
  # different fixed ve values and various w values.
  # ve=1.0 w=0.3, h=0.2 vl=1
  startstopspirl = [
      dict(dl=0, de='Re'),
      dict(dt=4, ve='ve', h='h', w='w', r=1),
      dict(dl=0, de='-(Re+Be)'),
      dict(dl=20, vl=5, h='h', w='w', r=0)]
  spiralStartStopargs = dict(name="startstop", lines=startstopspirl,
    ve=1.0, h=(0.18,0.38), dynret=False,
    tests=(
      dict(Re=2.9, Be=1.3, w=0.3),
      dict(Re=3.2, Be=1.1, w=0.4),
      dict(Re=4.0, Be=0.4, w=0.6),
      dict(Re=4.5, Be=0.0, w=0.8)))

  gen.startFile()
  #gen.doTests(n=args.n, **backpressure1args)
  gen.doSpiralTests(**spiralStartStopargs)
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
