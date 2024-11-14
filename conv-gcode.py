#!/bin/python3
"""
For this model we have;

t = time
e = extruded fillament length
l = extruded line length
z = extruded height/layer height
w = line width
h = line height/layer thickness
T = temperature of extruded material
Te = temperature of extruder
H = heat of extruded material
A = area of the layer
V = volume of extruded material

Temp is roughly heat divided by volume;

T = H/V

heat/temp decays at a rate proportional to temp and area, multiplied by a `hc`
convection constant that depends on air velocity;

dHc = hc * A * T * dt
    = hc * A * H/V * dt
    = hc * A/V * H * dt

Note that we care about temp for meltyness so we actually have;

dTc = dHc/V
    = hc * A * T * dt / V
    = hc * A/V * T * dt

So temperature decays at exactly the same rate as heat. This is exponential
decay with timeconstant hT which is;

dHc = H * dt/(hT+dt)
    ~= 1/hT * H * dt   # for dt small relative to hT

which gives us;

1/hT = hc * A/V
hT   = (V/A) / hc
    = z / hc

Note this suggests the heat decay rate decreases linearly as you go up the
layers. This would be true if the material was a perfect conductor with the
heat evenly distributed through the whole model and only loosing heat through
the top. However, most of the heat will be concentrated near the top layer,
with the lower layers mostly "cooled to equilibrium" with a constant
temperature and no heat transfer. This probably means only the top couple of
layers count, so z is effectively constant, making hT only proportional to 1/hc.

To calculate hT, a rough formula for hc for air velocities v in the range 2m/s
to 20m/s is;

# From https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
hc(v) = 10.45 - v + 10*sqrt(v)
hc(0) = 10.45
hc(5) = 27.81
hc(10) = 32.07
hc(20) = 35.17

Note the Adv3 extruder fan does 7.1m^/h (static pressure 58Pa), which
translates to about 19m/s through a 1cmx1cm nozzle. I'd assume it can only
squeeze about half of that through the blower nozzle, so about 10m/s would be
max blower speed?

So hc with fan on should be about 3x as hc with fan off, so fan on h1T should
be 1/3 of fan off h0T.

We extrude at constant temp Te, so heat extruded is proportional to volume. So
heat output ~= volume ~= length of extruded fillament.

dHe = Te*de

So combining heat added by extruded fillament

dH = dHe - dHc
   = Te*de - H*dt/(hT+dt)
H = H + dH
  = H + Te*de - H*dt/(hT+dt)
  = Te*de + H*hT/(hT+dt)  # heat added plus previous heat decayed.

Note if all the heat is in the last laid fillament and it is all at a
temperature considered "hot" then H effectively gives us the volume of laid
fillament that is still "hot". If we assume "hot" is Te, then H/Te is the
volume of fillament that the extruded heat could heat to the extruder
temperature. We can normalize and simpliy with Te=1 so it's just H. We could
just drive the fan to keep that matching a target constant.

However, the number of layers that are "hot" is H/layer_e. If we assume we
want to keep a constant number of layers hot, we should drive to keep that
constant.

Note A is the surface area of the top layer, which should be layer_l*w, but
that assumes constant line width, so perhaps using layer_e/h is better. If
only the top couple of layers count, then V should also only include those
layers. Perhaps A and V could be calculated using an exponentially decaying
sum, decaying as a function of layers instead of time?

V = layer_e + V*lT/(lT+1)
A = V/h

where lT is the "number of layers that count" decay rate. Note V is measured
in "mm of fillament" and A is in "mm of filament per mm of layer height".

Note Flashforge Adv3 automatically adjusts the extruder fan speed. I've seen
comments online claiming it "adjusts the speed based on the temperature", but I
assume it just adjusts it based on extrusion rate. Lets assume it varies the
speed so hc linearly matches extrusion rate.
"""

import re


def RelativeToAbsoluteE(data):
  lines=[]
  E=0.0
  for line in data.splitlines(keepends=True):
    g1 = re.match(r'(G1 .*E)(-?\d+(\.\d+)?)', line)
    if g1:
      E+=float(g1[2])
      line = re.sub(r'E(-?\d+(\.\d+)?)', 'E%.4f' % E, line, 1)
    lines.append(line)
  return ''.join(lines)

def FixFlashPrintLayerStartHop(data):
  """ At the start of each layer move the z-hop or layer-step immediately after the retract."""
  data=re.sub(r'(\n;layer:.*\n(?:M.*\n)*G1 E.*\n)((?:.*\n)*?)(G1 Z.*\n)'.encode(),r'\1\3\2'.encode(),data)
  return data

def ftime(t):
  h,s = divmod(t,3600)
  m,s = divmod(s,60)
  return f'{int(h):02d}:{int(m):02d}:{s:06.3f}'

def fbool(b):
  """ Format a boolean as 'Y', 'N', or '?' for True, False, None (unknown). """
  return '?' if b is None else 'Y' if b else 'N'

class Printer(object):
  fT = 5.0    # fan speed timeconstant.
  pT = 1.0    # fan pwm cycle time.
  h0T = 15.0  # decay rate of heat extruded with 0% fan.
  h1T = 5.0   # decay rate of heat extruded with 100% fan.
  ve1 = 1.0   # extrusion rate for max fan in mm/s.
  # the dynamic fan control giving fan output fo in the range 0->1 is;
  # fo = Kf*de/dt + Cf
  Cf = 0.2  # min fan while on is 20%.
  Kf = (1-Cf)/ve1
  # The hT decay rate is assumed to vary linearly with fan output.
  # hT = Kht*fo + ChT
  KhT = (h1T - h0T)
  ChT = h0T

  def __init__(self, heatext=0.0, fanspeed=1.0, layerpause=0.0):
    self.hs = heatext  # heat extruded setting.
    self.fanspeed=fanspeed  # fan speed scale setting.
    self.layerpause=layerpause
    self.gcode=[]
    self.fan_on = None
    self.t = 0.0
    self.x = self.y = self.z = self.e = self.f = 0.0
    self.ho = 0.0  # heat extruded output.
    self.fs = 0.0  # dynamicly adjusted fanspeed.
    self.fo = 0.0  # fan speed output.
    self.fd = 0.0  # fan speed drive.
    self.de = 0.0  # heat extruded since last fan cycle.
    self.dt = 0.0  # time since last fan cycle.
    self.de_dt = 0.0 # extrusion rate of last move command.
    self.layer_n = -1
    self._init_layer()

  def _init_layer(self):
    self.layer_n += 1
    self.layer_t = self.layer_l = self.layer_e = 0.0

  def getCode(self):
    return b'\n'.join(self.gcode + [b''])

  def runCode(self, code):
    code = FixFlashPrintLayerStartHop(code)
    for line in code.splitlines():
      try:
        self.runline(line.decode())
      except UnicodeDecodeError:
        # If it's not text, just append it.
        self.gcode.append(line)

  def _addline(self, line):
    self.gcode.append(line.encode())

  def fanstats(self, de, dt, fs, fd):
    fo, ho, hs = self.fo, self.ho, self.hs
    self._addline(
        f';{ftime(self.t)}: fan={fbool(self.fan_on)} '
        f'de/dt={de:.2f}/{dt:.2f} '
        f'fo/fs={fo:.2f}/{fs:.2f} '
        f'ho/hs={ho:.2f}/{hs:.2f} '
        f'fd={fd:.1f}')

  def setfan(self, fs):
    de, dt = self.de, self.dt
    if self.dt > self.pT and fs > 0.0:
      self.dt = self.de = 0.0
    fd = 1.0 if self.dt < fs*self.pT else 0.0
    if fd != self.fd or fs != self.fs:
      self.fanstats(de, dt, fs, fd)
    if fd != self.fd:
        self._addline('M106' if fd else 'M107')
    self.fs, self.fd = fs,fd

  def dofan(self, de, dt):
    self.dt += dt
    self.de += de
    if dt:
      self.de_dt = de/dt
    fo = min(1.0, self.fd*(self.Kf*self.de_dt + self.Cf))
    self.fo = (dt*fo + self.fT*self.fo)/(self.fT + dt)
    hT = self.KhT*fo + self.ChT
    self.ho = (de + self.ho)*hT/(hT + dt)
    # Turn on if there is excess heat extruded.
    on = self.fan_on and self.ho > self.hs
    self.setfan(self.fanspeed if on else 0.0)

  def layerstats(self):
    self._addline(
        f';{ftime(self.t)}: layer={self.layer_n} '
        f't={self.layer_t:.1f} l={self.layer_l:.1f} e={self.layer_e:.1f} '
        f'l/t={self.layer_l/self.layer_t:.3f} '
        f'e/t={self.layer_e/self.layer_t:.3f} '
        f'e/l={self.layer_e/self.layer_l:.3f} '
        f'ho={self.ho:.1f} fo={self.fo:.2f}')
    self._init_layer()

  def runline(self,line):
    dt=de=0.0
    if line.startswith(';start gcode'):
      self._addline(f';heat_ext: {self.hs}')
      self._addline(f';fan_speed: {self.fanspeed}')
      self._addline(f';layer_pause: {self.layerpause}')
      self._addline(line)
    elif line.startswith(';layer:'):
      #if self.layer_n > 0 and self.layerpause > 0.0:
      #  self.runline(f'G4 P{int(self.layerpause*1000)}')
      self.layerstats()
      self._addline(line)
    elif line.startswith(';end gcode'):
      self.layerstats()
      self._addline(line)
    elif line.startswith('M106'):
      dt = self.M106(line)
    elif line.startswith('M107'):
      dt = self.M107(line)
    elif line.startswith('G4 '):
      dt = self.G4(line)
    elif line.startswith('G1 '):
      dt, de = self.G1(line)
    else:
      self._addline(line)
    self.t += dt
    self.layer_t += dt
    self.dofan(de, dt)

  def M106(self, line):
    """ Fan on."""
    self.fan_on=True
    self.setfan(self.fanspeed)
    return 0.0

  def M107(self, line):
    """ Fan off. """
    self.fan_on=False
    self.fd=1.0  # Trick setfan() into always emitting an M107.
    self.setfan(0.0)
    return 0.0

  def G1(self, line):
    m = re.match(r'G1(?: X([-+0-9.]+))?(?: Y([-+0-9.]+))?(?: Z([-+0-9.]+))?(?: E([-+0-9.]+))?(?: F([-+0-9]+))?', line)
    x,y,z,e,f = m.groups()
    x = self.x if x is None else float(x)
    y = self.y if y is None else float(y)
    z = self.z if z is None else float(z)
    e = self.e if e is None else float(e)
    f = self.f if f is None else float(f)/60
    dx,dy,dz,de,df=x-self.x,y-self.y,z-self.z,e-self.e,f-self.f
    dl=(dx*dx+dy*dy+dz*dz)**0.5
    dt=(dl if dl else abs(de))/f
    # only count de and dl for extruding moves.
    de = de if dl else 0.0
    dl = dl if de else 0.0
    self.x,self.y,self.z,self.e,self.f = x,y,z,e,f
    self.layer_l += dl
    self.layer_e += de
    self._addline(line)
    #self._addline(f';{ftime(self.t)}: dl={dl:.3f} de={de:.3f} dt={dt:.3f} f={f:.3}')
    return dt, de

  def G4(self, line):
    """ wait. """
    # Turn fan full on during waits if fan_on.
    if self.fan_on:
      fs, self.fanspeed = self.fanspeed, 1.0
      self.setfan(1.0)
      # move around to blow the fan around for layerpause seconds.
      dx = self.layerpause / 1.01
      self.runline(f'G1 X{self.x+dx:.2f} Y{self.y:.2f} F300')
      self.runline(f'G1 X{self.x-2*dx:.2f} Y{self.y:.2f} F300')
      self.runline(f'G1 X{self.x+2*dx:.2f} Y{self.y:.2f} F300')
      self.fanspeed = fs
      self.runline(f'G1 X{self.x-dx:.2f} Y{self.y:.2f} F6000')
      dt = 0.0
    elif self.layerpause:
      # Change all layer waits to wait for layerpause duration.
      dt = self.layerpause
      self._addline(f'G4 P{int(dt*1000)}')
    else:
      m = re.match(r'G4 P([0-9]+)', line)
      dt = float(m[1])/1000.0
      self._addline(line)
    return dt


if __name__ == '__main__':
  import argparse, sys

  inf=float('inf')
  def FloatRangeType(min=-inf, max=inf):
    def init(s):
      f=float(s)
      if f<min or max<f:
        raise argparse.ArgumentTypeError(f'must be between {min} and {max}')
      return f
    return init

  cmdline = argparse.ArgumentParser(description='Postprocess FlashPrint gcode.')
  cmdline.add_argument('-e', type=FloatRangeType(0.0), default=0.0, help='Heat extruded length target.' )
  cmdline.add_argument('-f', type=FloatRangeType(0.0,1.0), default=1.0, help='Fan speed between 0.0 -> 1.0.' )
  cmdline.add_argument('-p', type=FloatRangeType(0.0), default=0.0, help='Seconds to pause between each layer.')
  cmdline.add_argument('infile', nargs='?', type=argparse.FileType('rb'), default=sys.stdin.buffer,
                    help='Input gcode file to postprocess.')
  args=cmdline.parse_args()

  data = args.infile.read()
  p=Printer(heatext=args.e, fanspeed=args.f, layerpause=args.p)
  p.runCode(data)
  data=p.getCode()
  sys.stdout.buffer.write(data)
