;start machine end gcode
M107
M104 S0 T0
M140 S0 T0
G162 Z F1800
G28 X Y
M132 X Y A B
{if activate_air_filtration[current_extruder] && support_air_filtration && complete_print_exhaust_fan_speed[current_extruder]>0}; Setting enclosure exhaust fan speed for completed print.
M651 S{complete_print_exhaust_fan_speed[current_extruder]*2.55}
{else}; Turning enclosure exhaust fan off for completed print.
M652
{endif}G91
M18
; end machine end gcode
