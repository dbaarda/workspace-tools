#!/bin/pypy3
"""
3D Printer Bed Z Calibration.

This generates gcode for drawing spiral pads on a 3x3 grid that can be
measured with calipers to calibrate the bed's z height.
"""
import argparse
import gcodegen


class ZCalibrationTest(gcodegen.GCodeGen):
  """ Generate GCode Z Height Calibration Test."""
  # Note we spiral from outside to inside to minimize initial loop retraction.
  PadR0=10.0
  PadR1=2.0
  lfw=0.5  # intentionally slightly under-extrude for a little space.
  vl0=10.0 # spiral speed
  xvals=(-60, 0, 60)
  yvals=(60, 0, -60)

  def preextrude(self, n=0, x0=-40, y0=-74, x1=40, y1=-74, v=20):
    self.cmt('TYPE:Custom')
    self.log('{layer=} {rd=} {rl=} {m=} {d=}')
    self.line(l=((x0, y0+n), (x1, y0+n)), v=v, w=1.0, r=1.0)

  def pad(self, x0, y0, r0=PadR0, r1=PadR1, vl=vl0, h=None, w=None, r=None, lfw=lfw):
    """ Draw a spiral pad."""
    if h is None: h = self.layer.h
    if w is None: w = self.layer.w
    if r is None: r = self.lfw2r(lfw, h, w) if lfw else self.layer.r
    dr = w if r0 < r1 else -w
    self.cmt('TYPE:Brim')
    self.hopdnrc(x0=x0, y0=y0, R=r0, C=0, h=h, w=w, r=r)
    self.spiral(x0=x0, y0=y0, R=r1, dRdC=dr, vl=vl, h=h, w=w, r=r)
    self.hopup()

  def doTest(self, xvals=xvals, yvals=yvals, **kwargs):
    """Draw 9 pads on platform."""
    self.preextrude()
    self.startLayer()
    for y in yvals:
      for x in xvals:
        self.cmt(f'Draw pad at ({x:.1f}, {y:.1f}).')
        self.pad(x, y, **kwargs)
    self.endLayer()


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
