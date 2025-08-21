#!/bin/bash
# Post process gcode to make it compatible with the Flashforge Adventurer3.

# Add a comment indicating how this has been post-processed.
ADD_HEAD="/^; *(start gcode|HEADER_BLOCK_END)/ i ; postprocess_cmd: $(basename $0) \"$@\"\n"

# Strip trailing comments and whitespace.
FIX_TCMT='/^[MG].*;/ s: *;.*$::'

# Replace extruder fan `M106 S<nnn>` cmds with `M106` (on) for S>=100 or `M107` (off) for S<100.
FIX_EFAN='s:^M106 S[0-9]{1,2}$:M107:; s:^M106 S[0-9]+$:M106:'

# Replace case fan `M106 P3 S<nnn>` cmds with `M651 Snnn` or `M652` (off).
FIX_CFAN='s:^M106 P3 S0$:M652:; s:^M106 P3 (S[0-9]+)$:M651 \1:'

# Delete any part fan `M106 P2 S<nnn>` cmds (ADV3 has no part fan).
FIX_PFAN='/^M106 P2/ d'

# Strip trailing decimal point of any temp settings.
FIX_TEMP='/^M(140|104) .*\./ s:\.0::'

# Replace any G0 cmds with G1.
FIX_G0='s:^G0 :G1 :'

# Strip any unsupported commands.
#FIX_CMDS='/(G20)( |$)/ d'

#TODO: Split G1 cmds with XY and Z into two commands? Strip redundant G1 Z cmds?

SCMD="
${ADD_HEAD};
${FIX_TCMT};
${FIX_EFAN};
${FIX_CFAN};
${FIX_PFAN};
${FIX_TEMP};
${FIX_G0};"
#${FIX_CMDS};

sed -ri "$SCMD" "$@"
