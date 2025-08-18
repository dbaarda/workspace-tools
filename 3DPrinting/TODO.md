# 3DPrinting TODO

## Code changes

### gcodegen.py

* Fix layer.n = -1 initialization and FlashPrint Preextrude being its own
  layer.n=0. Instead we should just make layer.n=0 the initialized layer and
  make pre-extrudes normal draws the same as it is for OrcaSlicer.

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
