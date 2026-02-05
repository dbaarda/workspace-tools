#!/bin/pypy3
"""
3D Printer Bed Z Calibration.

This generates gcode for drawing spiral pads on a 3x3 grid that can be
measured with calipers to calibrate the bed's z height.
"""
import argparse
import gcodegen
from math import *

class ZCalibrationTest(gcodegen.GCodeGen):
  """ Generate GCode Z Height Calibration Test."""
  # Note we spiral from outside to inside to minimize initial loop retraction.
  PadR0=10.0
  PadR1=2.5
  lfw=0.5  # intentionally slightly under-extrude for a little space.
  vl0=10.0 # first layer spiral speed.
  xvals=(-60, 0, 60)
  yvals=(60, 0, -60)
  lvals=((0.3, 1, 10), (0.2,0.25,None), (0.2,0.5,None), (0.2,1,None), (0.1,1,None))
  #lvals=((0.2,0.25),)*5

  def preextrude(self, n=0, x0=-40, y0=-74, x1=40, y1=-74, v=20):
    self.cmt('TYPE:Custom')
    self.log('{layer=} {rd=} {rl=} {m=} {d=}')
    self.line(l=((x0, y0+n), (x1, y0+n)), v=v, w=1.0, r=1.0)

  def spiral_infill(self, x0, y0, r0=PadR0, r1=PadR1, d=0.25, vl=None):

    """Do a spiral infill."""
    n = self.layer.n
    w = self.layer.w
    # We alternate spiral direction each layer. For 90 degree crossings of
    # alternating layers we need 45 degree spirals, or `dRdC=r*2*pi`, which
    # varies at different radius. We aim for 45 degree spirals at the radius
    # that includes half the area, which is `rm=sqrt((r1^2+r0^2)/2)`.The density
    # is `d = s*w/dRdC` so for a given density and width we need `s=d/w*dRdC`
    # for a given dRdC, or `dRdC=s*w/d` for a given number of spirals.  This
    # means ideally we want `s=d/w*r*2*pi` spirals.  We can only have an integer
    # number of spirals so we round this, and then adjust dRdC for the exact
    # density.
    rm = sqrt((r1*r1 + r0*r0)/2)
    rm = min(r0,r1)
    spirals=round(d/w * rm * 2*pi)
    # Odd layers have r1<r0 and spiral anti-clockwise for the first spiral.
    # This means dR and dRdC are both negative for odd layers, and both
    # positive for even layers.
    dR = r1 - r0
    dRdC = copysign(spirals*w/d, dR)
    dCdS = 1/spirals        # change in starting rotation per spiral
    dCdL = 9/(7*spirals*pi) # Change in starting rotation per layer.
    dC = dR/dRdC
    ci = n*dCdL # inside edge rotation angle of first spiral for layer.
    c0 = ci if r0<r1 else ci-dC
    # draw single line walls.
    self.cmt('TYPE:Outer wall')
    # place the seam between where infill lines touch the wall.
    for r,c in ((r0, c0 + dCdS/2),(r1, c0 + dCdS/2 + dC)):
      self.hopdnrc(x0=x0, y0=y0, R=r, C=c, w=w)
      self.spiral(x0=x0, y0=y0, R=r, dC=copysign(1, dR), dRdC=0, vl=vl, w=w)
      self.hopup()
    # Adjust r0,r1 inwards by w/2 for infill with a bit of overlap.
    r0+= copysign(w/2, dR)
    r1-= copysign(w/2, dR)
    dR = r1 - r0
    dC = dR/dRdC
    # draw infill.
    self.cmt('TYPE:Sparse infill')
    for s in range(spirals):
      self.hopdnrc(x0=x0, y0=y0, R=r0, C=c0, w=w)
      self.spiral(x0=x0, y0=y0, R=r1, dRdC=dRdC, vl=vl, w=w)
      self.hopup()
      # swap r0,r1 and increment c0 for the next spiral.
      r0,r1,dC,c0=r1,r0,-dC,c0 + dC + dCdS

  def pad(self, x0, y0, r0=PadR0, r1=PadR1, vl=None, h=None, w=None, r=None, lfw=lfw):
    """ Draw a spiral pad."""
    if h is None: h = self.layer.h
    if w is None: w = self.layer.w
    if r is None: r = self.lfw2r(lfw, h, w) if lfw else self.layer.r
    dr = w if r0 < r1 else -w
    self.cmt('TYPE:Bottom surface' if self.layer.n == 1 else "TYPE:Top surface")
    self.hopdnrc(x0=x0, y0=y0, R=r0, C=0, h=h, w=w, r=r)
    self.spiral(x0=x0, y0=y0, R=r1, dRdC=dr, vl=vl, h=h, w=w, r=r)
    self.hopup()

  def doTest(self, xvals=xvals, yvals=yvals, lvals=lvals, r0=PadR0, r1=PadR1, **kwargs):
    """Draw 9 pads on platform."""
    self.preextrude()
    for h,d,v in lvals:
      self.nextLayer(h=h, Vp=v)
      for y in yvals:
        for x in xvals:
          if d == 1:
            self.cmt(f'Draw pad at ({x:.1f}, {y:.1f}).')
            self.pad(x, y, r0, r1, **kwargs)
          else:
            self.spiral_infill(x,y,r0,r1,d)
      r0,r1 = r1, r0
    # Add labels.
    self.nextLayer(h=0.2)
    for r , y in enumerate(yvals):
      for c, x in enumerate(xvals):
        self.text(f'{r},{c}', x-5 , y+max(r0,r1)-2)


if __name__ == '__main__':
  import sys, argparse

  cmdline = argparse.ArgumentParser(description='Generate 3D Printer Bed Z Calibration gcode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  gcodegen.GCodeGenArgs(cmdline)
  args=cmdline.parse_args()

  gen=ZCalibrationTest(**gcodegen.GCodeGetArgs(args))

  gen.startFile()
  gen.doTest()
  gen.endFile()
  data=gen.getCode()
  sys.stdout.buffer.write(data)
  #sys.stdout.write(data)
