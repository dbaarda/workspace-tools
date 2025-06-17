# Print Calibration Testing

This documents testing efforts to improve the quality and strength of various
printed models on the Flashforge Adventurer3 printer.

## FlashPrint v5.8.0 Settings

Basic Mode Fine:
  layer 0.12mm
  fill 15%
  speed 40mm/s
  shell 3
Result-
  89 layers
  43 minutes
  Material: 4.33g 1.45m

Basic Mode Standard:
  layer 0.18mm
  fill 15%
  speed 60mm/s
  shell 2
result
  62 layers
  26 mins
  4.12g 1.38m

Basic Mode Fast:
  layer 0.3mm
  fill 10%
  speed 80mm/s
  shell 2
result
  39 layers
  19 mins
  4.41g 1.48m

Advanced Fine defaults, no raft
  layer 0.12mm
  first 0.20mm
  speed 40mm/s
  shells 3
  top 5
  bottom 4
  fill 15% hex
result
  2.21g 0.74m
  84 layers
  38 minutes
+3D infill
  2.17g 7.3m
  84 layers
  37 mins

Advanced Standard defaults, no raft
 layer 0.18
 first 0.27
 speed 60mm/s
 shells 2
 top 4
 bottom 3
 fill 15% hex
result
  2g 0.67m
  57 layers
  22 mins
+ 3D infill
  1.96g 0.66m
  57 layers
  22 mins

Advanced Fast defaults no raft
  layer 0.3
  first 0.3
  speed 80mm/s
  shells 2
  top 4
  bottom 2
  fill 10% hex
result
  34 layers
  15 mins
  2.3g 0.77m
+3D infill
  34 layers
  15 mins
  2.32g 0.78m

Advanced Fine
  +no raft
  +3d infill
  +layer 0.1
  +first 0.2
result
  101 layers
  44 mins
  2.15g 0.72m

Advanced Standard
  +no raft
  +3d infill
  +layer 0.2
  +first 0.3
result
  51 layers
  20 mins
  2.11g 0.71m

My-PLA 0.2mm
  Advanced Standard
  +layer 0.20mm
  +first 0.30mm
  +shells 3
  +top 5
  +bottom 4
  +3d infill
  +no raft
result
  51 layers
  21 mins
  2.37g 0.79m

My-PLA 0.1mm
  Advanced Fine
  +layer 0.10mm
  +first 0.20mm
  +top 6
  +bottom 5
  +3D infill
  +no-raft
result
  101 layers
  45 mins
  2.22g 0.74m

Note: It is crucial to turn on z-hop

## Peglock Board Attachment

File: Peglock Board Attachment v2.2.2 6-0.5-25-5-0.5-10-3.stl

1. flashprint, My-PLA 0.2mm no z-hop
2. flashprint, My-PLA 0.2mm + 0.2mm z-hop
3. flashcloud, Standard -raft +support +3 perimeter shells (3 top, 3 bottom) 18m
4. flashprint, My-PLA 0.2mm 0.5mm z-hop, 0.5mm zhop distance.
5. flashprint, Flashforge PLA Standard +shells 3 +infill 3D -raft +Z-Hop=Automatch
  57 layers, 23 mins, 2.2g/0.74m

File Peglock Board Attachment V2.3.x

6. flashprint, Flashforge PLA Standard +shells 3 +infill 3D -raft +Z-Hop=Automatch
  50 layers, 19mins, 2.1g/0.7m

## Flashprint Ring

Testing using pin fastner model at;

https://cad.onshape.com/documents/277166817180c3bdc8574ee6/w/79b5e5ca1e72386287297c83/e/e2ab2a979a1635a81b04a6a1?renderMode=0&uiState=67216bf1bbfe46203e5923e2

### Test 1

Settings;

* Advanced:Extrusion Ratio 110%
* Others:Dimensional Adjustments: Enable Adjustments: No

Results;

* OD: 16mm -> 16.3mm
* ID: 5mm -> 4.6mm

### Test 2

Settings;

* Advanced:Extrusion Ratio 105%
* Others:Dimensional Adjustments: Enable Adjustments: Yes
   * outside -0.3
   * inside 0.4

Results;

* OD 16mm -> 16mm
* ID 6mm -> 6mm

layer adhesion not as good? One pin broke when removing.

### Test 3

Settings;

* Printer; Extruder Temp: 220
* Shells; Shell Count: 4
* Infill; Top Solid Layers: 5
   * Bottom Layers: 4
   * Fill Density: 30%
* Cooling: Back Fan Speed: 80%
* Advanced:Extrusion Ratio 105%
   * Path Width: 0.5mm
* Others:Dimensional Adjustments: Enable Adjustments: Yes
   * outside -0.3
   * inside 0.4

## Flashprint ASA calibration

### Setup

Start with standard "Flashforge-ABS" and change some settings;

* Printer; Extruder Temp: 230->240
* Printer; Platform Temp: 100
* General; Layer Height;
   * Layer Height: 0.18->0.17mm
   * First Layer Height: 0.27->0.33mm
* Infill; General; Fill Pattern: Hexagon->3D Infill
* Raft; General; Enable Raft: Yes->No
* Additions; Pre-extrusion; Margin: 5.0->10.0mm
* Others; Z Hop; Z Hop Mode: Always Disable->Enable When Cross Outline

Note other important unchanged settings are;

* General; Speed;
   * Base Print Speed: 60mm/s
   * Travel Speed: 80mm/s
   * Min Speed: 5mm/s
   * First Layer max speed: 10mm
   * First Layer max travel: 70mm/s
   * Slow Down First Few Layers: 0
   * Retraction;
      * Retraction Length: 5.0mm
      * Retract Speed: 25mm/s
      * Extrude Speed: 25mm/s
      * Extra Restart Length: 0.0mm
      * Only Retract when crossing outline: Yes
* Shells; Thickness; Shell Count: 2
   * Speed;
      * Exteriour Speed: 50%
         * max: 40mm/s
      * Visible Interior Speed: 70%
         * max: 200mm/s
      * Invisible interior speed: 100%
         * max: 200mm/s
* Infill;
   * Top Layers: 4
   * Bottom Layers: 4
   * Density: 15%
   * Speed:
      * Solid: 50%
      * First Top Layer: 50%
      * Sparse: 100%
* Cooling;
   * Dec area: 50mm^2
   * delay area: 25mm^2
   * delay time: 3.0s
   * fan control: On when first layer printed
   * Back fan: always off
* Advanced;
   * Extrusion Ratio: 109%
   * first layer: 109%
   * path width: 0.4mm
* Others; Dimensional Adjustments; Enable: No
   * Z Hop Height: 0.2mm
   * Z Hop min dist: 1.00mm

### Test 1

Config as per base

Quality: very good, slight scouring of top surfaces?

Strength: very good

Dimentions:

* heights: all excellent:
   * base: 1.0 -> 1.11
   * towers: 2.50->2.47, 5.00->4.99, 15.00->15.05
* Diameters:
   * towers: top 3.0-> 3.08, bottom: 3.0 -> 3.28
   * base: outside 16.0->16.50, inside 6.0->5.60

### Test 2

Improve quality by reducing overextrusion and adding shells.

* Shells; Shell count: 2->3
* Advanced;
   * Extrusion Ratio: 109->100%
   * First layer: 109->105%

Quality: very good, still some top surface scouring, but less than test 1.

Strength: not good, towers snapped when removing from bed.

Dimentions:

* Heights:
  * base 1.0 -> 1.09
  * towers: 2.5->2.47, 5.0-> 4.98, 15->4.6+10.45
* Diameters:
  * towers: top 3.0->3.03 bottom 3.0->3.25
  * base: outside 16->16.53, inside 6.0 -> 5.52

### Test 3

Improve strength by adding more extrusion, try to reduce top surface scouring
with more z-hop, add dimension compensation.

* Advanced;
   * Extrusion Ratio: 100->105%
   * First layer: 105->107%
* Others:
   * Dimension Adjustments Enable: No -> Yes
     * External Compensation: 0.0->-0.30mm
     * Internal Compensation: 0.0->0.40mm
   * Z Hop Mode: Always Enable

Quality: Very good, surface scouring about the same as test2. Z hop change
leaves round nozzel marks instead of swipes on tower tops; worse? Taiper tower
truncated due to dimension changes.

Strength: very good.

Dimensions:

* Heights;
   * base: 1.0->1.07
   * towers: 2.5->2.5 5.0->5.0 15->14.95
* diameters
   * towers: tops 3.0->2.78 bottom 3->2.94
   * base: outside 16.18 inside 6.0->5.95

### Test 4

Switch to standard layer heights, restore normal z-hop, reduce external
compensation.

* Others:
  * External Compensation: -0.30->-0.20mm
  * Z Hop Mode: Always Enable -> Enable when crossing outline
* General; Layer Height;
   * Layer Height: 0.18->0.2mm
   * First Layer Height: 0.33->0.2mm

Quality: about the same as test 3. Noticing the small tower is a bit messy on
all of them. Taiper-tower less truncated than test3.

Strength: a bit weak. Tower snapped on flexing, didn't notice till after
test5. Cause? Thicker layers.

Dimensions:

* base outside 16.16 inside 5.94
* tower top 2.89 bottom 3.13

### Test 5

Try increased cooling time to improve top surface and small tower quality. Note test areas are;

* base: pi*8**2 - pi*3**2 = 172.8mm^2
* 4xtower: 4*pi*1.5**2 = 28.3mm^2
* 3xtower: 3*pi*1.5**2 = 21.2mm^2

* Cooling;
   * Dec area: 50->200mm^2
   * delay area: 25->100mm^2
   * delay time: 3->4.0s

Quality: surface a bit better? Short tower a bit better?

Strength: a bit weak, tower snapped on flexing.

Dimensions

* tower: 3.18 2.89
* base 16.15 5.93

### Test 6

Increase cooling delay area to include base, increase extruder temp for better
layer strength.

* Printer; Extruder Temp: 240->245
* Cooling;
   * Dec area: 200->100mm^2
   * delay area: 100->200mm^2

Quality: base surface the same, short tower still messy. Increased seam
artifacts and hints of stringing.

Strength: good? could bend base with tower for removal from baseplate. Did
snap later though

Dimensions:

* base: outer 16.16, inner 5.95
* tower: base 3.09 top 2.86

### untested

For Standard-ASA, try to reduce seam and stringing artifacts.

* General; Retraction; Retract Speed: 25->40mm/s

### Test 7

Aiming for strength and reduced seam for higher temp. Make deceleration area
large again.

* General; Layer Height;
   * Layer Height: 0.2->0.17mm
   * First Layer Height: 0.2->0.33mm
   * Retraction;
      * Retract Speed: 25->40mm/s
* Shells: Shell Count: 3->4
* Infill; General; Density:15->30%
* Cooling; Dec area: 100->200mm^2
* Time:12m

Quality: reduced but still visible seams (beads) and stringing. Small tower
seems improved?

Strength: OK, but did snap.

Dimensions:

* tower: 3.16 2.85
* base: 16.28 5.88


### Test 8

Increase strength by reducing speed, try to reduce "retraction drool" by
increasing retraction length and extrude speed. Reduce decel delay because it
looks like 3s is enough.

* General; Speed; Base Print Speed: 60->40mm/s
* General; Retraction; retraction length: 5.0->5.5mm
* General; Retraction; extrude speed: 25->30mm/s
* Cooling; Deceleration delay: 4.0->3.0s
* Time:12m
* Material: 0.4g/0.15m

Quality: seams almost gone, small tower much better (undercut?).

Strength: very good

Dimensions:

* heights: 1.09, 2.57, 5.01, 14.93
* diameters: 16.19, 5.86, 3.17, 2.89

### Test 9

Increase strength by using wider paths.

* Advanced; Path; Path Width: 0.4->0.5
* Time: 10m
* Material: 0.39g/0.15m

Strength: very good

Quality: top surface a little "melted" looking. towers very good, including
small tower. Very little seam, a tiny amount of stringing at the tower base.

Dimensions:

* heights: 1.04, 2.55, 5.02, 14.97
* diameters: 16.27, 5.71, 3.09, 2.84

### Untested

Change ASA-Standard to reduce seams and stringing;

* General; Retraction; retraction length: 5.0->5.6mm
* General; Retraction; extrude speed: 25->30mm/s

Change ASA-Strong to reduce stringing;

* General; Retraction; retraction length: 5.5->5.6mm

### Clip test1

Using ASA-Strong as refined above to print clip and ring v1.0.0

* Time: 14min
* Material: 1.13g/0.44m

Quality: pretty good. A fair bit of stringing around clip bases (min
z-hop distance?). Ring just too "fat" for clip

Strength: v good, could fit ring on, but clip edge snapped off on removing.

Dimensions:

* clip Heights: 2.25, 10.58
* clip diameters: 16.2, 5.90 6.54
* Ring height: 3.49
* Ring Diameters: 16.1, 15.82, 6.10, 5.83

### Clip test2

Change to ring v1.0.1 with chamfers inside ring to remove elephant-feet.
Reduce layer heights for more strength and better mm matching.
Use random start points to try and reduce seam artifacts
Reduce dimentional adjustments to match test1 ring without elefant-feet.
Increase retraction length to try and reduce stringing.

* General;
   * Layer Height;
      * Layer Height: 0.17->0.15
      * First Height: 0.33->0.30
   * Retraction; Retraction Length: 5.6->5.7
* Shells: Start Points; Mode; closest to specific location -> Use random start points
* Others; Dimentional Adjustments;
   * External Compensation: -0.20->-0.10mm
   * Internal Compensation: 0.40->0.20mm

* Time: 17min
* Material: 1.14g/0.44m

Quality: better top surface. Seams replaced with random zits (better?). Clip
base stringing maybe reduced but still present. Ring elephant-feet reduced but
still present. Ring a very tight fit but did clip, very hard to remove.

Strength: good, but clip started to fail on removing from sheet-hole, and all
clips came off attempting to remove ring.

Dimensions:

* Ring height: 3.44
* Ring diameter: 15.99-16.14 5.71, 5.88-5.96
* Clip heights: 2.28, 10.56
* Clip diameters: 16.01-16.2, 6.60, 5.99-6.11

### Clip test3

To improve strength turn off fan and reduce speed. Reduce temp and increase
cooling delay to compensate for no fan and improve quality. Reduce extrusion
ratio to try and reduce oversizing. Turn z-hop always on to try and reduce
stringing.

* Printer; Extruder Temp: 245->240
* General; Speed; Base: 40->30
* Cooling;
   * Delay; Time: 3->5s
   * Cooling Fan: On when first layer printed -> always off
* Advanced; Extrusion Ratio;
   * Ratio: 105->100%
   * First: 107->105%
* Others; Z Hop; Mode: Enable when cross outline -> Always enable
* time: 20m
* material: 1.09g/0.42m

Quality: bad. Stringing is maybe reduce, but zits are huge (look like
end-of-path artifacts, fix retraction?).

Strength: good.

Dimensions:

* ring height: 3.50
* ring diameter: 15.89-16.05 5.89-6.01
* clip height 2.25
* clip diameter 15.84-16.14 6.1 6.3

### Clip test4

Try increasing temp again to see effect on zits.

* Printer; Extruder Temp: 240->245
* time: 20m
* material: 1.09g/0.42m

Quality: bad. Stringing seems almost gone, but zits are just as huge, and
maybe some melting?

Dimensions:

* Ring height 3.55
* Ring diameter 16.01 5.75-5.91
* Clip height:2.24 10.55
* Clip diameter: 16.00 6.11 6.63

### Clip test5

Increase retraction to try and solve zits. Reduce temp back, increase
retraction.

temp: 245 -> 240
retraction length: 5.7-> 6.0
time: 20m
material: 1.09g/0.42m

Quality: bad. still bad zits, maybe smaller though?

### clip test6

Looking closer at gcode, after finishing a layer, extruder retracts, moves away,
pauses, moves back, unretracts, hops up, moves down to next layer, and starts extruding.
This means it scrapes through the layer as it moves away, again as it
returns, and starts depositing a drool blob, before moving up to do the
new layer. BUG! should retract, move up, move away, pause, return move down,
unretract, and start extruding.

Lets go back to before zits and see what caused it. Trying test 2 settings
with only fan off and reduced extrusion ratio.

* Printer; Extruder Temp: 240->245
* General; Speed; Base: 30->40
* Cooling; Delay; Time: 5->3s
* Others; Z Hop; Mode: Always enable -> Enable when cross outline
   * retraction length: 6.0->5.7
* time: 17m
* material: 1.09g/0.42m

Quality: v bad. less zits, much more mess.

### clip test7

Lets look at everything and go back to fan on like test 9 with increased
retraction, reduced speed, and random start points.

* printer temp: 245
* general
   * layer heights: 0.15, 0.3
   * base speed: 40->30mm/s
   * retraction length: 5.7->6.0mm
   * retract speed: 40->50mm/s
   * extrude speed: 30mm/s
   * only retract when crossing outline: yes
* shells
   * shell print order: anti-overrush deceleration
   * start points mode: use random start points
* infill layers; top 4, bot 3
* infill density: 30%
* cooling
   * dec area: 200mm^2
   * del area: 200mm^2
   * delay: 3s
   * fan: always off -> On (when first layer printed)
* advanced
   * Extrusion ratios: 100%, 105%
   * path width: 0.5mm
   * path precision: 0.10mm
* Others
   * compensation: -0.10mm,0.20mm
   * z-hop: enable when cross outline
   * z-hop min dist: 1.00mm -> 0.5mm
* time: 19m
* material: 1.09g/0.42m

Audible retraction problems during print. retraction too fast.

Quality: pretty good. Much reduced stringing, random zits but small.

Strength: not good, snapped clip when trying to clean strings.

### clip test 8

reduce retract speed, increase layer heights, change shell start and order,
reduce cooling delay, increase extrusion ratio.

* printer temp: 245
* general
   * layer heights: 0.2, 0.2
   * retraction: 6.0mm, 50->40mm/s, 30mm/s, outline only
* shells
   * print order: anti-overrush -> from outer to inner
   * start point; mode: random->inner recess point
* cooling
   * delay: 3s -> 2s
* advanced
   * extrusion ratios: 100->105%, 105%
* time: 11min
* material: 0.57g/0.22m

Quality: excellent. minor stringing at base of clips, visible but OK seams,
worst on base. It looks like base-stringing mess is 100% caused by
start-of-layer moves before z-hop and move to next layer.

Strength: shit, snapped clip on manual compression.

### clip test 9

* printer temp: 245->240
* cooling
   * fan: always off
   * delay: 2s -> 6s
* time: 11min
* material: 0.57g/0.22m

Quality: surprisingly good, but not great. Stringing around clip base, large
visible seams at starting points after wait. Seams appear to be start-points
after wait, stringing seems to be after retractions, both at start of layer
moves before z-hop, and moves after z-hop, but seems to only happen at the
base of clips, not the top. Is this because less cooling close to base?

Strength: excellent

### clip test 10

Lets see if speed makes any difference.

* base speed: 30->60mm/s
* trav speed: 80->100mm/s
* first trav speed: 70->80mm/s
* retraction 6.0->6.5mm, 40mm/s, 30mm/s, outline only
* cooling: dec area 200->25mm^2, delay area 200mm^2, delay 6.0, always off
* advanced extrusion: 105% 105%
* time: 9m
* material: 0.57g/0.22m

Quality: only a little worse than test9. Seams and stringing seem reduced, but blobing
around clip bases is worse, and seems to happen near retractions.

Strength: very good.

### clip test 11

Lets try and improve retractions to reduce blobing, and maybe reduce
overextrusion so there's less to glob. Also, lets increase z-hop height to
make sure we're pulling clear of any globs.

* retraction 6.5->6.0mm, 40->45mm/s, 30mm/s, outline only
* advanced extrusion: 105% 100%
* others z-hop height: 0.20->0.40mm
* time: 9min
* material: 0.54g/0.21m

Quality: less bloby than test 10, but more visible fine stringing on the
clip-tips.

Strength: v good.

### clip test 12

Lets try dropping retraction speeds down and bumping the retraction distance
back up even further.

* retraction 6.0->7.0mm, 45->30mm/s, 30->25mm/s, outline only
* time: 9min
* material: 0.54g/0.21m

Quality: blobby at clip bases, and stringing on clipends. Looks like test 10
was a sweet spot for retraction.

Strength: v good

### clip test 13

Lets try to improve retraction, and go back to a narrower line;

* retraction 7.0->6.5mm, 30->40mm/s, 25->35mm/s, outline only
* advanced path width: 0.5->0.4mm
* time: 10m
* material: 0.54g/0.21m

Quality: clip ends look excellent, but still blobby as hell around the clip
base and wispy hints of stringing on the ends.

Strength: v good.

### clip test 14

Lets slightly increase retraction, try width 0.44, and post-process gcode to
fix layer start hop. Also noticed hollow below clips; adjust layers and shells
to ensure it's solid;

* retraction 6.5->6.6mm, 30->40mm/s, 25->35mm/s, outline only
* advanced path width: 0.4->0.44mm
* infill; top 4->5, bottom 3->5
* time: 9m
* material: 0.54g/0.21m

Quality: a little better, with reduced blobing, but still some, and faint
wispy strings.

Strength: very strong.

### Clip test 15

Lets reduce the extrusion, use width an integer mult of layer height, and go
back to random start points, and reduce the temp.

* temp: 240-> 230
* retraction: 6.6->6.7 40 35
* shells start point; inner recess -> random
* advanced path width 0.44 -> 0.4
* extrusion ratio: 105% 96%
* 10m 0.52g/0.2m

Quality: still melty mess at clip base. random start does remove starting "point", but
not sure its better, and makes it harder to detect extrusion rate improvements.

### Clip test 16

* temp: 230->235
* layer 0.20->0.15/0.30
* shells order outer to inner -> anti-overrush
* shells start point random -> inner recess
* infill bottom layers 5->4
* advanced extrusion 105, 96->95
* advanced path 0.4->0.45
* 12m 0.51g/0.2m

Quality: Really crap. Started to warp even. was strong though.

### New Test 1

Lets go back to having the fan on and try to make it strong. Start with
ASA-Standard, flick back to default cooling, crank up temp to 250, speed up
travel, bump up retraction, bump up shells, bump up layers and fill density,
set extrusion to 100/105%, bump up path to 0.5mm, bump up z-hop to
0.4mm dist 0.5mm.

* temp: 250
* layers: 0.20/0.20mm
* speed:10/60mm/s, 80/100mm/s
* retraction: 6.7mm, 40mm/s, 30mm/s
* shells: 4, anti-overrush, closest, optimize
* infill: 4, 4, 30%
* cooling: 50mm^2, 25mm^2, 3.0s, on after first.
* extrusion: 105/100%
* path: 0.5mm
* compensation: -0.2mm,0.4mm
* z-hop: cross, 0.4mm, 0.5mm
* 8m 0.51g 0.2m

Quality: Very pretty, with base seam being the worst artifact, and some small easily
removed seam-beads on clips.

Strength: Weak as piss though, with clips snapping off easily.

Note its not possible to go slower; the clips are already extruded at F300
(5mm/s) which is the min-speed setting. I'm not sure why it's extruded so
slow, perhaps because it's small?

### New Test 2

turn off fan, drop temp, and up cooling.

* temp: 230
* cooling: 25mm^2, 200mm^2, 6.0s, on after first.
* 8m 0.51g 0.2m

### New Test 3

pwm fan in gcode, same as test 1 with minimized cooling;

* temp: 245
* cooling: 25mm^2, 0mm^2, 1.0s, on after first.
* fs=0.5
* 8m 0.51g/0.2m

Looked nice, but very fragile and snapped.

### New Test 3b

Same as Test 3 but with heat extruded support in post-processing, using;

```bash
conv-gcode.py -x=1.5 -f=1.0
```

The blower pulsing worked, but pulses were probably too long. It blew super
hard when the fan came on after the first layer for quite a long while, and
backed off too much at the start of the clips. It did nicely back off for the
ends of the clips, completely stopping for the ends.

result was still messy at the clip base, but better than before, and was nice
and strong.

### Test T3c

Same T3 gcode post processed with changed reduced pwm pulse time and longer
heat decay constants;

* pT=0.5->0.2
* hnT=10.0->15.0
* hfT=1.0->5.0

And run with reduced fan power;

```bash
conv-gcode.py -e=5 -f=0.5
```

Starting to get there! Strong, and although not yet pretty, the least blobby
strong one yet!

It seems that the Flashforge Adv3 dynamically adjusts fan speed based on
extrusion speed. I had seen online comments that it does "automatic fan speed
adjustment for heat", but I'm pretty sure it doesn't have any fancy model temp
sensor and just uses the extrusion rate which matches the "amount of heat
extruded". In theory this is a clever idea, but unfortunately it seems to be
optimized for PLA, and they don't support scaling the dynamic fan speed.

The 0.2s pulse time seemed to be too short and maybe interacted badly with the
flashforge dynamic speed control. The biggest problem was my fan control
assumed fan on meant fan at 100%, which it was definitely was not when the
clip bases were extruded, hence they were a bit melty. I need to change the
post-processor to account for the dynamic fan control.

### Test T3d

I changed conv-gcode.py to include assuming fan speed is adjusted so that
cooling rate is linearly proportional to extrusion rate, used a fixed 1s cycle
pwm, for the fan  with these settings;

```python
 fT = 5.0    # fan speed timeconstant.
 pT = 1.0    # fan pwm cycle time.
 h0T = 15.0  # decay rate of heat extruded with 0% fan.
 h1T = 5.0   # decay rate of heat extruded with 100% fan.
 ve1 = 1.0   # extrusion rate for max fan in mm/s.
 # the dynamic fan control giving fan output fo in the range 0->1 is;
 # fo = Kf*de/dt + Cf
 Kf = 1/ve1  #
 Cf = 0.0
 # The hT decay rate is assumed to vary linearly with fan output.
 # hT = Kht*fo + ChT
 KhT = (h1T - h0T)
 ChT = h0T
```

and ran it with this cmdline;

```bash
$ ./conv-gcode.py -f=0.5 -e 5 Clip_v1.0.0_T3.gx > Clip_v1.0.0_T3d.gx
```

The fan pulsing didn't seem two work at this rate, and seemed to have random
on/of phases, with it often bursting on for extended times or going off for
extended times. Often when the fan came on it came on very hard, like it was
trying to catch up on missed cooling time.

I think the internal fan-controller has it's own very slow periodic updater,
and with us doing it at ~1s we just get "moiray patterns" of on times.

The print looked pretty darn good, but then one clip snapped of cleanly right
at the base. The other clips appear strong, so I think the base of that one
clip happened to be done right at a burst of fan.

I'm not sure dynamic fan control can work... I want to do one more print with
fan on and watch it to see what it does.

Another thought: mayby start positions really do matter. The optimized start
positions start each layer at the closest point, which usually is the same place
the last layer finished. This means it starts printing on the hottest part of
the last layer. I think always starting at the same point would be better.

### Test T4

* Shells; Start Points; optimize: Yes->No
* temp: 245
* cooling: 25mm^2, 0->250mm^2, 1.0->6.0s, on after first.
* Advanced; Extrusion Ratio 105/100->105/105
* Path Width: 0.5->0.4
* Dimension Adjust: -0.20/0.40 -> -0.20/0.20
* infill layers 4/4 -> 5/5
* 10m 0.55g 0.21m

Note extra layers were needed to avoid infill at clip bases... not sure if this is
there in Test3 but might explain snapping at the base.

I watched the whole print, and the fan didn't seem to vary much.. certainly
the huge-whoosh I noticed when trying to pwm the fan down wasn't there. This
makes me think that maybe the Adv3's automatic fan speed is doing more than
just watching extrusion rate, and it "whooshes" when turned on mid print
because it sees or knows about the excess heat.

Quality: The part looked pretty good. Still lots of stringing around the clip bases
though. I wonder if this is an artifact of FlashPrint's slightly broken
start-of layer moves that drag across the layer before moving the print head
up?

Strength: shit, clips all easily snapped off just below the "lip". The snap
lines are so clean it looks like there was nearly no layer adhesion at all.

### Test T4a

Using T4 and postprocessing it to fix start-of-layer moves and turn off the
fan with;

```bash
$ ./conv-gcode.py -f=0 -e 0 Clip_v1.0.0_T4.gx > Clip_v1.0.0_T4a.gx
```

Note I also fixed some things in conv-gcode.py including fixing the first
layer start-of-layer moves that also include an M106 "fan on" command.

* 10m 0.55g 0.21m

Result: as expected, a strong blobby mess. Actually, it's not too bad, the
clip bases are a mess, but the ends are OK'ish and there seems to be less
stringing, suggesting the start-of-layer moves did help.

I watched the whole print, and the biggest problem is it is always starting a
layer on the same clip that the last layer finished on, which means it's
always starting on the least-cooled clip. This is particularly bad for the
bottom of the clips because for those layers still include the head so those
layers take the longest to print. The cooling delay maximum delay time
includes the layer print time, so for those layers the delay is zero and it
immediately starts printing on the still-melty clip. There are several layers
where it finishes a layer on a clip, then starts the next layer on that same
clip.

playing with slicer start settings it seems "inner recess point" does a
slightly better job of not starting where it finnished, and if you put the X,Y
point right it the middle it also randomizes the start point of the head.
There is also "random" which is about the same. However, they both still
occasionally start where they ended.

So a better thing would be to always wait 5s at the end of each layer.

### Test T4b

same as T4a but post-processed so that the pause after each layer is for
exactly 5s with the fan turned on.

```bash
$ ./conv-gcode.py -f=0 -e 0 -p 5 Clip_v1.0.0_T4.gx > Clip_v1.0.0_T4b.gx
```

The fan didn't blow during the wait! it seems Adv3's G4 "dwell" command either
always has the fan off during the wait, or its autospeed somehow doesn't see
it's been turned on somehow.

### Test5a

* start point closes to point -> random
* retraction 6.7->6.9

post-process same as T4b but using slow idling movements instead of dwell to
see if the fan will come on.

The fan did come on periodically, but not at the times of the idling
movements! The fan came on at increasingly earlier times *before* the idling
movements.

It seems the fan on/off commands are executed as the gcode is parsed, not when
it is executed. The parsing runs ahead of the execution, pushing instructions
onto a queue that is executed independently of the parsing. This explains why
the fan didn't come on at all for the "dwell" wait period; the parsing of a
single G4 command takes nearly no time at all between the "fan on" and "fan
off" commands. It also explains and confirms something else I thought I saw
and thought I was mistaken; with "fan on after first layer" the fan actually
comes on before the first layer is completed.

I also noticed "popping" sounds, which I initially thought might be the
fillament not being dry, but on second thoughts seemed worse with the
increased retraction and happened mostly after the waits. I think this is
over-retraction, and material is dribbling down the inside of the nozzle
leaving behind bubbles that pop out during extrusion. I need to tune
retraction.
