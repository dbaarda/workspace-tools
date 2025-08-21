# 3DPrinting TODO

## Code changes

### gcodegen.py

* Fix layer.n = -1 initialization and FlashPrint Preextrude being its own
  layer.n=0. Instead we should just make layer.n=0 the initialized layer and
  make pre-extrudes normal draws the same as it is for OrcaSlicer.

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

* Add Move dd attribute for combined dl,dz distance. Should it include de
  distance too? Note de is typically way less than 0.5% of dl (less than 15%
  if de is in mm^3) so it adds nearly nothing to normal move distances, but
  means retract/restore have non-zero dd which might be useful?

* Improve pressure compensation to adjust ve0/ve/ve1 of future moves, not just
  try to stay within jerk limits.

* Change stats to be added from incstate() instead of incmove(), and have it
  accumulate and collect data at 1mm increments instead of per-move.

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

* improve ADv3 gcode support to strip trailing comments, remove unsuported
cmds, and translate different gcode cmds.

* improve OrcaSlicer support to include moving final layer retract/wipe/hopup
cmds before the start of the next layer.

* Factor ADv3 vs OrcaSlicer support into mixins and make them optional. Maybe
separate mixins for parsing input vs generating output?


## Testing

* Do spiral tests to measure pressure advance.

* Test pressure advance needed for overhangs and bridging. In theory the
  backpressure should be proportional to the fraction not overhanging.

* Figure out some kind of final test-suite for measuring/calibrating printers.
