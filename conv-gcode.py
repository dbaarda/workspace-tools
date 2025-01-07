#!/bin/pypy3
"""
Post-process gcode files to fix various problems.

This was primarily written to fix FlashPrint files to improve the printed
results. It was created because it doesn't seem possible to configure
FlashPrint to give good results with high temperature ASA filament. It has the
following features;

* Fixes some minor FlashPrint bugs.
* Implements a better heat model for dynamically adjusting cooling in response
  to extruded heat.
* Can vary fan output for printers with fixed or automatic fan speed by
  inserting PWM fan on/off instructions into the gcode.
* Can adjust or extend wait times between layers beyond the slicer's limits.
* Can dynamically optimize retract/restore distances based on Linear Advance
  for printers that don't implement that feature, using nozzle pressures to
  relieve/restore the right amount of pressure for the extrusion rates
  before/after.

# FlashPrint Bugs

The main bug this fixes is flashprint moves the extruder between layers at the
extrusion height of the previous layer, and doesn't z-hop or raise it to the
new layer height until just before it starts extruding. This means it creates
drag-marks on top of the layers, causes drag-stringing artifacts when crossing
outlines, and can even cause layer-shifts from hitting outline edges. This
fixes the gcode by shifting the z-hop or layer height change command to just
after the previous end-of-layer retraction.

# Heat/Cooling improvements.

This uses a heat model to track how much heat is in the extruded material to
control the amount of cooling to apply. It can reduce the cooling fan output
by inserting PWM like fan on/off commands into the gcode, or adjust the wait
time between layers.

For the heat  model we have;

t = time (seconds)
e = extruded fillament length (mm of fillament)
l = extruded line length (mm of layer-line)
z = extruded height/layer height (Z mm)
w = line width (mm)
h = line height/layer thickness (mm)
T = temperature of extruded material (degC)
Te = temperature of extruder (degC)
H = heat of extruded material
A = area of the layer (mm^2)
V = volume of extruded material (mm^3)

Heat is energy and Temp is roughly heat divided by volume;

T = H/V

heat/temp cooling through convection decays at a rate proportional to temp
difference to ambient and surface area, multiplied by a `hc` convection
constant that depends on air velocity;

dHc = hc * A * T * dt
    = hc * A * H/V * dt
    = hc * A/V * H * dt

Note that we care about temp for meltyness so we actually have;

dTc = dHc/V
    = hc * A * T * dt / V
    = hc * A/V * T * dt

So temperature decays at exactly the same rate as heat. This is exponential
decay with timeconstant Th which is;

dHc = H * dt/Th

which gives us;

1/Th = hc * A/V
Th   = (V/A) / hc
     = z / hc

Note this suggests the heat decay rate decreases linearly as you go up the
layers. This would be true if the material was a perfect heat conductor with
the heat evenly distributed through the whole model and only loosing heat
through the top. However, most of the heat will be concentrated near the top
layer, with the lower layers mostly "cooled to equilibrium" with a constant
temperature and no heat transfer. This probably means only the top couple of
layers count, so z is effectively constant, making Th only proportional to
1/hc.

To calculate Th, a rough formula for hc for air velocities v in the range 2m/s
to 20m/s is;

# From https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
hc(v) = 10.45 - v + 10*sqrt(v)
hc(0) = 10.45
hc(1) = 19.45
hc(5) = 27.81
hc(10) = 32.07
hc(20) = 35.17

Note the Adv3 extruder fan does 7.1m^/h (static pressure 58Pa), which
translates to about 19m/s through a 1cmx1cm nozzle. I'd assume it can only
squeeze about half of that through the blower nozzle, so about 10m/s would be
max blower speed?

So hc with fan on should be about 3x as hc with fan off, or fan on h1T should
be 1/3 of fan off Th0. But note the fan at only 10% with v=1m/s has hc 2x Th0,
so a fan at only 10% will double the rate of cooling.

We extrude at constant temp Te, so heat extruded is proportional to extruded
volume. So heat output is proportional to extruded volume which is
proportional to length of extruded fillament.

dHe = Te*de

So combining heat added by extruded fillament with heat lost through
convection we can get the curren heat H from the previous heat H' after time
dt with;

dH = dHe - dHc
   = Te*de - H*dt/Th

           H = H' + dH
           H = H' + Te*de - H*dt/Th
(Th+dt)/Th*H = H' + Te*de
           H = (H' + Te*de)*Th/(Th+dt)  # heat added plus previous heat decayed.

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

The Flashforge Adventurer3 reputedly automatically adjusts the extruder fan
speed. I've seen comments online claiming it "adjusts the speed based on the
temperature", but I don't think it has extra sensors for that, so maybe it
just adjusts it based on extrusion rate? Observing it in action it doesn't
seem that straight-forward; it does blow hard when first turned on after lots
of hot extrusion, but it doesn't seem to vary much in response to extrusion
speed.

Note for electric fans, rotational speed is linearly proportional to PWM duty
cycle, and airflow is linearly proportional to rotational speed. This means
airflow velocity v should vary linearly with "percent on". With option -A we
assume that fan output is automatically adjusted linearly between 10% to 100%
for extrusion rates between 0mm/sec to 1mm/sec for calculating the amount of
cooling when the fan is on.

The worst thing is the Adv3 executes fan-on and fan-off commands at gcode
parse time, not at execution time. This means fan toggling on/off in gcode
happens before the gcode before it is completely executed. The FlashPrint "fan
on after first layer" feature puts a fan-on command after the first layer in
the gcode, but you will notice the fan comes on before the first layer is
completely printed. This makes it impossible to implement reduced fan output
by toggling the fan on and off in the gcode.

This means that the gcode PWM fan control this implements with option -P
doesn't work for the Adv3, and the only way we can control cooling is by
adding or changing wait times between layers.

Note that FlashPrint also starts printing the next layer at the same place
where the previous layer finished, which means it starts printing on top of
the least-cooled part of the previous layer. The various shell start point
options give only limited ability to prevent this. It seems using
"Shells;Start Points" set to "use random start points" or "inner recess point"
with "Permit optimize start points" set to "No" works best, but even then it
will sometimes start layers on top of the end of the previous layer. This
means the wait time between layers cannot really include the layer print time
in the wait duration, which is what the "Cooling; Maximum Delay Time" does,
and it also cannot be set higher than 6secs. So the -p option forces the wait
time to always be the specified duration which can be arbitrarily large.

# Linear Advance Retract/Restore.

For printers that don't implement Linear Advance, the correct retract/restore
distance depends on the amount of nozzle pressure (measured in extruder
advance distance) that needs to be relieved on retraction, and applied on
restore. The right amount of pressure to relieve, and thus length of
retraction, depends on the extrusion rate before the retraction, and and the
right amount of pressure to restore, and thus length of restore, depends on
the extrusion rate after the restore.

The Flashforge Adventurer3 doesn't seem to implement Linear Advance.
FlashPrint does have a setting for it and does emit the M900 command to set
it, but tests setting it between 0.0 to 3.0, even doing each different Kf
value in a different print in case it is executed at parse time like fan
on/off, showed no difference. (TODO: maybe it has pre 1.5 linear advance which
is measured in steps with values typically in the range 0-150. Test if I just
didn't set it high enough).

I have experienced problems related to this in prints that have layers with a
mixture of large fill areas that are printed fast, and fiddly fine features
that are printed slow. The retraction after completing the fast fills relieves
a large amount of advance pressure, which is then restored in the restore
before printing the slow fiddly bits. This restores too much pressure which
comes out as a blob at the start of the slow print, making a mess out of it.
This can also contribute to stringing artifacts. At first I mistook this for a
cooling problem which prompted the fancy cooling model included here.

So this can modify the gcode changing all the restores and retracts to use
dynamic lengths based on the linear advance pressure from the print speeds
before and after.
"""

import re
import gcodegen

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

def fbool(b):
  """ Format a boolean as 'Y', 'N', or '?' for True, False, None (unknown). """
  return '?' if b is None else 'Y' if b else 'N'

def hc(v):
  """ Get hc as a function of air velocity v in m/s."""
  assert v>= 0
  return 10.45 - v + 10*(v**0.5)


class Printer(gcodegen.GCodeGen):
  Tf = 5.0    # fan speed timeconstant.
  Tp = 1.0    # fan pwm cycle time.
  # The Th0 heat decay rate is the seconds it takes to cool 63%, or
  # about 2x the time to cool from 245degC to 160degC in a 30degC enclosure.
  Th0 = 15.0  # decay rate of heat extruded with 0% fan.
  # The Th decay rate varies with fan output fo as;
  # Th = Kht/hc(fva1*fo)
  fva1 = 10.0  # airflow velocity for max fan in m/s.
  KTh = Th0*hc(0)
  # For printers with automatic fan speed control we assume fan speed is a
  # linear function of extrusion rate. The fan output fo in the range 0->1 is;
  # fo = fd*(Fk*de/dt + Fc)
  fve1 = 1.0   # extrusion rate for max fan in mm/s.
  Fc = 0.1  # min fan speed while fan is on is 10%.
  Fk = (1-Fc)/fve1

  def _Th(self, fo):
    """Get Th temperature decay rate for fo fan output."""
    assert fo >= 0
    return self.KTh / hc(self.fva1 * fo)

  def _fo(self, fd, de_dt):
    # Get fo fan output for autofan if enabled, else fo = fd fan drive.
    assert de_dt >=0.0
    return fd*min(1.0, self.Fk*de_dt + self.Fc) if self.en_autofan else fd

  def __init__(self, Hs=0.0, Lp=0.0, en_pwmfan=False, en_autofan=True, **kwargs):
    super().__init__(**kwargs)
    self.Hs = Hs  # heat extruded setting.
    self.Lp=Lp # duration to pause after each layer.
    self.en_pwmfan=en_pwmfan  # whether to control fanspeed with gcode pwm.
    self.en_autofan=en_autofan  # whether the printer autocontrols fan speed.

  def resetfile(self):
    """ Reset all state variables. """
    super().resetfile()
    self.fan_on = None
    self.ho = 0.0  # accumulated heat extruded out.
    self.fs = 0.0  # current dynamicly adjusted fanspeed.
    self.fo = 0.0  # filtered fan speed output.
    self.fd = 0.0  # current fan speed drive.
    self.fe = 0.0  # heat extruded since last fan cycle.
    self.ft = 0.0  # time since last fan cycle.
    self.de_dt = 0.0 # extrusion rate of last move command.

  def runCode(self, code):
    self.resetfile()
    code = FixFlashPrintLayerStartHop(code)
    for line in code.splitlines():
      try:
        self.runline(line.decode())
      except UnicodeDecodeError:
        # If it's not text it's *.gx file binary data; just append it.
        self.add(line)

  def fanstats(self, de, dt, fs, fd):
    fo, ho, Hs = self.fo, self.ho, self.Hs
    self.log(
        f'fan={fbool(self.fan_on)} '
        f'de/dt={de:.2f}/{dt:.2f} '
        f'fo/fs={fo:.2f}/{fs:.2f} '
        f'ho/Hs={ho:.2f}/{Hs:.2f} '
        f'fd={fd:.1f}')

  def setfan(self, fs):
    de, dt = self.fe, self.ft
    if self.ft > self.Tp and fs > 0.0:
      self.ft = self.fe = 0.0
    if self.en_pwmfan:
      fd = 1.0 if self.ft < fs*self.Tp else 0.0
    else:
      fd = fs
    if fd != self.fd or fs != self.fs:
      self.fanstats(de, dt, fs, fd)
    if fd != self.fd:
      if 0.0 < fd < 1.0:
        self.cmd('M106', S=f'{int(fd*255)}')
      else:
        self.cmd('M106' if fd else 'M107')
    self.fs, self.fd = fs,fd

  def dofan(self, de, dt):
    self.ft += dt
    self.fe += de
    if dt:
      self.de_dt = de/dt
    # Get and filter fan output fo.
    fo = self._fo(self.fd, self.de_dt)
    self.fo = (dt*fo + self.Tf*self.fo)/(self.Tf + dt)
    # Get Th and calculate extruded heat ho.
    Th = self._Th(fo)
    self.ho = (de + self.ho)*Th/(Th + dt)
    # Turn on if there is excess heat extruded.
    on = self.fan_on and self.ho > self.Hs
    self.setfan(self.Fe if on else 0.0)

  def layerstats(self):
    super().layerstats()
    self.log(f'  ho={self.ho:.1f} fo={self.fo:.2f}')

  def runline(self,line):
    """ Parse, execute, and add a line of gcode. """
    #print(line)
    dt=de=0.0
    if line.startswith(';start gcode'):
      self.cmt('heat_ext: {Hs}')
      self.cmt('layer_pause: {Lp}')
      self.add(line)
      self.startLayer(0, 0.0)
    elif line.startswith(';preExtrude:'):
      self.Layer(line)
    elif line.startswith(';layer:'):
      self.Layer(line)
    elif line.startswith(';end gcode'):
      self.layerstats()
      self.add(line)
    elif line.startswith('M106'):
      dt = self.M106(line)
    elif line.startswith('M107'):
      dt = self.M107(line)
    elif line.startswith('G4 '):
      dt = self.G4(line)
    elif line.startswith('G1 '):
      dt, de = self.G1(line)
    else:
      self.add(line)
    self.dofan(max(de, 0.0), dt)

  def Layer(self, line):
    """ Start a new layer. """
    t, h= re.match(';(preExtrude|layer):(\d+.\d+)', line).groups()
    h = float(h)
    if t == 'preExtrude':
      self.startLayer(n=0, z=0.0, h=h)
      self.add(line)
    else:
      self.layerstats()
      self.startLayer(h=h)

  def M106(self, line):
    """ Fan on."""
    self.fan_on=True
    self.setfan(self.Fe)
    return 0.0

  def M107(self, line):
    """ Fan off. """
    self.fan_on=False
    self.fd=1.0  # Trick setfan() into always emitting an M107.
    self.setfan(0.0)
    return 0.0

  def G1(self, line):
    m = re.match(r'G1(?: X([-+0-9.]+))?(?: Y([-+0-9.]+))?(?: Z([-+0-9.]+))?(?: E([-+0-9.]+))?(?: F([-+0-9]+))?', line)
    x,y,z,e,f = (a if a is None else float(a) for a in m.groups())
    # Always use the previous f value when parsing gcode files.
    if f is None: f = self.f
    self.move(x=x, y=y, z=z, e=e, v=f/60)
    # Get the last added move to get it's dt and de. Note it could be logged
    # as skipped if it didn't do anything.
    m=self.gcode[-1]
    return (m.dt, m.de) if isinstance(m, gcodegen.Move) else (0.0, 0.0)

  def G4(self, line):
    """ wait. """
    if self.pe > self.Re:
      # We should never pause without retracting to relieve nozzle pressure.
      self.log(f';WARNING: pause with residual nozzle pressure, pe={self.pe:.4f}>{self.Re:.4f}')
    if self.Lp:
      # Change all layer waits to wait for Lp duration.
      dt = self.Lp
    else:
      m = re.match(r'G4 P([0-9]+)', line)
      dt = float(m[1])/1000.0
    self.wait(dt)
    return dt


if __name__ == '__main__':
  import argparse, sys
  from gcodegen import RangeType

  cmdline = argparse.ArgumentParser(description='Postprocess FlashPrint gcode.',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  gcodegen.GCodeGenArgs(cmdline)
  cmdline.add_argument('-Hs', type=RangeType(0.0), default=0.0,
      help='Heat extruded target in mm of filament that is hot.' )
  cmdline.add_argument('-Lp', type=RangeType(0.0), default=0.0,
      help='Seconds to pause between each layer.')
  cmdline.add_argument('-F', action='store_true',
      help='Control fan speed by pulsing it on/off.')
  cmdline.add_argument('-A', action='store_true',
      help='Printer automatically adjusts fan speed when on.')
  cmdline.add_argument('infile', nargs='?', type=argparse.FileType('rb'), default=sys.stdin.buffer,
      help='Input gcode file to postprocess.')
  args=cmdline.parse_args()

  data = args.infile.read()
  p=Printer(
      Te=args.Te, Tp=args.Tp, Fe=args.Fe, Fc=args.Fc,
      Kf=args.Kf, Kb=args.Kb, Cb=args.Cb, Re=args.Re,
      Vp=args.Vp, Vt=args.Vt, Vz=args.Vz, Ve=args.Ve, Vb=args.Vb,
      h=args.Lh, w=args.Lw, r=args.Lr,
      en_dynret=args.R, en_dynext=args.P, en_optmov=args.O, en_verb=args.v,
      Hs=args.Hs, Lp=args.Lp, en_pwmfan=args.F, en_autofan=args.A)
  p.runCode(data)
  data=p.getCode()
  sys.stdout.buffer.write(data)
