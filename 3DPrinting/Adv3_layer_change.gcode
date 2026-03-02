;starting layer {layer_num+1}/[total_layer_count]{if close_fan_the_first_x_layers[current_extruder]>0 && layer_num==close_fan_the_first_x_layers[current_extruder] && fan_min_speed[current_extruder]>0}
;Turning extruder fan on after first [close_fan_the_first_x_layers] layers.
M106{endif}
