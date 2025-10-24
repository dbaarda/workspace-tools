# 3DPrinting TODO

## Code changes

### gcodegen.py

* Change all Move de volumes to be in mm^3, not mm of fillament. Note this
  also changes ve to mm^3/s. This means translation into E values is only done
  in the final stages of generating gcode, and everything else doesn't need to
  know the fillament diameter to calculate volumes from line dimensions. Note
  however that Ve and Ae velocity and acceleration limits will also need to be
  in mm^3/sec and mm^3/sec^2, so will need to be calculated from the extruder
  limits using Af.

* Change Move v0,vm,v1 to also be start/mid/end speeds for adjusts and hops,
  not just horizontal moves. Add/fix vl0,vlm,vl1,vz0,vzm,vz1,ve0,vem,ve1
  attributes for start/mid/end speeds in the various dimensions, calculating
  them from v0,vm,v1 as necessary depending on the move type.

* Improve pressure compensation to adjust ve0/ve/ve1 of future moves, not just
  try to stay within jerk limits.

* Factor out advanced drawing methods into an optional mixin.

* Add support for getting/setting printer acceleration/velocity limits from
  gcode cmds. OrcaSlicer optionally outputs these so we can set them in the
  Slicer UI.

* Add support for getting/setting Kf (and maybe Kb?) from gcode cmds. This
  means we can support gcode pressure-advance tests or adaptive pressure
  advance that change Kf during the print.

* Add support for modeling hardware PA and applying additional pressure
  adjustments on top of it. This could be useful for leveraging hardware PA
  while still adding additional PA improvements.

Note: I considered adding support to Moves for dh, dw, dr for lines that
change their height/width/ratio, but gcode doesn't have support for changing
the line extrusion rate during a draw, so this can't be done.

### conv-gcode.py

* Factor out heat model, pressure model, and input parsing into optional mixins.

* Improve cmdline options to support optional overriding of temps, fans,
  pressure advance, velocitys, layer/line properties, or just leave/use them
  unmodified from the input gcode.

* improve and test heat model.

* improve ADv3 fan support to include per-layer fan toggling.

* Factor ADv3 vs OrcaSlicer support into mixins and make them optional. Maybe
separate mixins for parsing input vs generating output?

## Testing

* Do spiral tests to measure pressure advance.

* Test pressure advance needed for overhangs and bridging. In theory the
  backpressure should be proportional to the fraction not overhanging.

* Figure out some kind of final test-suite for measuring/calibrating printers.

## OrcaSlicer

Improve OrcaSlicer https://github.com/SoftFever/OrcaSlicer with;

* Get https://github.com/SoftFever/OrcaSlicer/pull/10747 MarchingSquares
  optimizations approved.

* Get https://github.com/SoftFever/OrcaSlicer/pull/10876 FillTpmsFk
optimizations approved.

* Add FillAdaptiveTpms fill methods.

* Fix multiline_fill() to use offset() for pollylines. See 
  https://github.com/SoftFever/OrcaSlicer/pull/10876#issuecomment-3393839823
  
* Improve distance, length, and coordinate code everywhere.

## Other


