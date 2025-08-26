; start machine start gcode
M118 X75 Y75 Z150 T0
M140 S[bed_temperature_initial_layer] T0
M104 S[nozzle_temperature_initial_layer] T0
M104 S0 T1
M107
M652
M900 K[pressure_advance] T0
G90
G28
M132 X Y Z A B
G1 Z50 F420
G161 X Y F3300
M7 T0
M6 T0
M108 T0
{if close_fan_the_first_x_layers[current_extruder]>0}; Leaving extuder fan off for first [close_fan_the_first_x_layers] layers.
{elsif fan_min_speed[current_extruder]>0}; Turning extruder fan on for all layers.
M106
{else}; Leaving extruder fan off for all layers.
{endif}{if activate_air_filtration[current_extruder] && support_air_filtration}; Setting enclosure exhaust fan speed during print.
M651 S{during_print_exhaust_fan_speed_num[current_extruder]}
{else}; Leaving enclosure exhaust fan off during print.
{endif};end machine start gcode
;start preextrude
G1 X-37.5 Y-74 F6000
G1 Z0.2 F420
G1 X37.5 Y-74 E9.5 F1200
;end preextrude
