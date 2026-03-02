# Printer calibration and Control

## Printer characteristics.

Low end printers like my Flashforge Adventurer 3 have X and Y speeds up to
100mm/s and acceleration of about 500mm/s^2.

A printer with 500mm/s^2 acceleration can go from 0 to 100mm/s in 0.2s over
10mm of distance. For 0 to 60mm/s it is 0.12s over 3.6mm of distance.

Note that diameter 0.4mm nozzle vs 1.75mm filament means that extruding
`de=1mm` of filament translates to nearly 20mm of nozzle thread. For layer
`h=0.2mm` with line `w=0.4mm` for 0.2x0.4mm line areas, `de=1mm` is 30mm of
printed line. For layer 0.1x0.4mm lines lines it is 60mm, and for 0.3x0.4mm
it's 20mm. So line print velocity `vl` is between 20x to 60x, or typically
30x, extruder velocity `ve`.

## Line width and spacing.

OrcaSlicer uses a line-model that assumes "bulges" on the edge of the lines,
as described at;

https://manual.slic3r.org/advanced/flow-math 

But this ends up almost exactly the same as simpler "rectangular" lines
because they always use overlap_factor=1. See;

https://github.com/SoftFever/OrcaSlicer/issues/8313#issuecomment-3082244171

This means Orca effectively uses thinner lines than other slicers by
layer_height * (1 - PI/4), or roughly layer_height/5. So add that much to your
line_width setting. Note that the adjustment depends on layer_height, so if
you have a different first_layer_height, have enabled Precise_Z_height (which
can adjust layer heights to give a more accurate model height), or enabled
infill combination (which can merge infill layers into one thick layer), then
the adjustment will be off for those layers. Note that this also means Orca
effectively uses narrower lines for thicker layers with the same line-width
settings, which is IMHO the opposite of what you would want. For
nozzle-width=0.4mm and layer_heights=0.1, 0.2, 0.3 I set Orca line-widths as a
% of the nozzle width and add 5%, 11% and 16% respectively.

Orca effectively shifts the outer-wall inwards by layer_height * (1-PI/4)/2,
or roughly layer_height/10. This can be compensated for by adding it to
"Quality; Precision; X-Y contour compensation" and subtracting it from
"Quality; Precision; X-Y hole compensation". Note in Orca the "X-Y contour
compensation" and "X-Y hole compensation" distance settings are offsets
(radius). I believe FlashPrint's equivalent "External Compensation" and
"Internal Compensation" settings are edge-to-edge (diameters) so you would
need to halve those settings to map them to Orca.

My thoughts on how line-width and extrusion-ratio settings should work is at;

https://github.com/SoftFever/OrcaSlicer/pull/11105#issuecomment-3443404410

## Seam Gap.

OrcaSlicer has a "Quality; Seam; Seam Gap" setting that leaves a space between the seam join
ends.

If we assume the start and end of the lines include the circular "bead" bulge
beyond the start and stop points, we can calculate the gap needed for the
bulge overlap to exactly match the volume of no seam.

gap = line_spacing * pi/4 ~= 0.7854 * line_spacing

The setting can be set as a % of nozzle diameter or length, with the default
being 10%. For a 0.4mm nozzle this is 0.04mm, which is much less than our
calculated 0.3142mm for a 0.4mm width line.

## Linear Advance

Linear Advance or Pressure Advance is mm of extrusion advance per mm/s of
filament extrusion speed, and has typical values in the range 0.0-2.0. Note
that 100mm/s print speed for a 0.3x0.4 lines using 1.75mm diameter filament is
12mm^3/s or 5mm/s filament rate. This suggests the advance could get as high
as 10mm, or 200mm of line!

This means 5mm of uncompensated pressure advance translates to more 100mm of
line smear or stringing.

Linear advance assumes flow rate out the nozzle is linear with pressure, and
pressure is linear with `pe` extrusion advance of filament (how many mm of
filament the extruder has pushed that has not yet come out the nozzle and is
compressed filament). This is based on assuming the filament and boden tube
between the extruder and nozzle behaves like a spring with force linear with
compression distance as per [Hooke's
law](https://en.wikipedia.org/wiki/Hooke%27s_law), and nozzle flow rate is
linear with pressure in the nozzle as per [Poiseuille's
law)(https://en.wikipedia.org/wiki/Poiseuille_law). The filament force divided
by filament cross-section area (nozzle input area) is the pressure in the
nozzle, so nozzle flow is linear with filament compression distance, or the
extruder advance distance. See this for some thoughts on these assumptions;

https://klipper.discourse.group/t/modification-of-pressure-advance-for-high-speed-bowden-printers/13053/18?u=dbaarda

Also more thoughts posted in OrcaSlicer discussions;

https://github.com/SoftFever/OrcaSlicer/discussions/10832#discussioncomment-14603955

Note the linear advance factor `Kf` for compensating for this is the mm of
advance `pe` needed per mm/sec of extruder velocity for the nozzle output flow
rate to match the extruder velocity.

```python
  pe = Kf * de/dt
```

Where

```
  de/dt is the extruder velocity in mm of filament/sec.
  Kf is the linear advance factor in mm per mm/sec of extruder velocity.
  pe is the mm of filament linear advance "pressure" distance.
```

Note that at the steady state, the filament rate into the nozzle `de/dt`
equals the filament rate out of the nozzle `dn/dt`. This means you can replace
`de/dt` with `dn/dt` in the above equation and get the flow rate out of the
nozzle as a function of advance `pe`;

```python
  pe = Kf * dn/dt
  dn/dt = pe/Kf
```

The rate of change in `pe` is `dpe/dt` and is equal to the flow rate in the
extruder minus the flow rate out the nozzle;

```python
  dpe/dt = de/dt - dn/dt
        = de/dt - pe/Kf
```

Which gives us the change in `pe` of `dpe` over time `dt`;

```python
  dpe = de - pe * dt/Kf
```

Note that this is an exponential decay equation, with `pe` advance
exponentially approaching the steady state value for a fixed `de/dt` rate with
a timeconstant of `Kf`. Exponential decay processes cause both lag and a
low-pass-filter "frequency cutoff" effects that depends on the timeconstant.

The lag effect means that changes in the extruder input rate will take roughly
`Kf` time in seconds to appear at the nozzle output. More specifically if you
drive the extruder velocity with a sin-wave, the peaks and troughs in the
sin-wave nozzle output will be exactly `Kf' seconds behind the peaks and
troughs in the sin-wave input.

The low-pass-filter "frequency cutoff" effect means that for a sin-wave
extruder input with amplitude `Ae`, if the frequency `fe` is higher than
`fc=1/(2*pi*Kf)`, then the nozzle output amplitude `Az` will be reduced by
`Az=(fc/fe)*Ae`. This roughly means changes in extruder speed that happen
every 2x faster than every `pi*Kf` seconds will only get half way to the max
and min target nozzle output rates (but note those output peaks and troughs
will also lag by `Kf` seconds).

So for `Kf=1.0` it will take about 1 seconds of extruding before the extrusion
rate "catches up" and is extruding at about the right rate, and switching
speeds faster than every 1.5sec will mean it only gets about half way there.
At 100mm/s, 1 seconds is 100mm worth of line. Note that also means 100mm of
track ooze you would get if you keep moving without retracting after you stop
extruding a line at 100mm/s.

So the pressure can take some time to build up enough for the right extrusion
rate when extruding lines at constant velocity. It also can take some time to
decay back enough when you change to a lower extrusion velocity.

Printers that implement linear advance dynamically adjust the extrusion rates
to take this into account, effectively compensating for the lag and frequency
cutoff effects by adjusting extrusion rates in advance. See the following for
how it's implemented;

https://mmone.github.io/klipper/Kinematics.html

## Retract and Restore

Retraction withdraws the extruder when you stop extruding, and restore
advances it again before you start extruding. This helps avoid stringing when
moving between extrusions. The idea is that retraction sucks the filament back
into the nozzle to prevent "drool", and restore pushes it back to the tip of
the nozzle again. There are lots of calibration tests to help you figure out
the right amount of retraction to minimize stringing.

However, when you look at Linear Advance works and how the nozzle output rate
varies with pressure, it becomes pretty clear that what retraction actually
does is relieve the accumulated advance pressure, and even more importantly,
restore re-applys that pressure! This means when you print at constant speed,
the advance pressure initially starts low giving under-extruded lines, but
builds up over the `Kf` timeconstant duration to give the right extrusion
rate. When the printer stops and retracts, it quickly relieves that pressure
so it can move without extruding, and when it restores, it re-applys that
pressure so the printer is ready to print at the previous speed without
under-extruding at the start.

What this means is the ideal retraction distance is the accumulated pressure
advance distance. If the printer implements linear advance, and it is
implemented and tuned perfectly, this means you shouldn't need any retraction
or restore at all, but perhaps a tiny token amount to compensate for tiny
miss-calibration errors and backlash would be wise. However, if the printer
doesn't implement linear advance, then the ideal retraction distance depends
on the previous extrusion speed, and the ideal restore distance depends on the
**next** extrusion speed.

If the printer always printed at constant extrusion rates, using a constant
sufficiently large retraction would be fine, since the pressure to releave and
reapply would be the same. However, slicers don't use constant extrusion
rates, they slow down for outer edges and corners. Also printers cannot
instantaneously reach the requested line-speed, and have to accelerate and
decelerate, adjusting extrusion rates along the way.

I have been battling under/over extrusion problems and retraction tuning for a
print that has quite a lot of speed variation due to having a large base (fast
print) with small raised clips (slow print). The small clips had lots of
stringing and over extrusion for the layers that also included large base
areas, but the clips were fine for the higher layers without the large base
areas. No amount of tuning retraction and extrusion rates could fix this.
These pressure advance affects totally explain why.

It should be possible to implement linear advance in gcode for printers that
don't implement it in firmware. It would involve dynamically adjusting
retraction/restore distances and explicitly inserting
acceleration/deceleration line-segments with extrusion rates adjusted to
build/decay the linear advance amount.

Even without implementing full linear-advance including acceleration in gcode,
it would probably be worth implementing dynamic retract/restore based on the
previous/next extrusion rates. This would probably fix my clip print.

## Print Layer BackPressure

After playing around with trying to figure out the `Kf` value for my Adv3
printer, it's become apparent that the simple linear advance model is
insufficient to accurately reflect the real nozzle flow characteristics. It
seems the relationship between nozzle flow rate and advance pressure is not
linear. This is also discussed in various places online like;

https://klipper.discourse.group/t/modification-of-pressure-advance-for-high-speed-bowden-printers/13053/26

The reported experiences vary, with some reporting `Kf` seems to drop with
increasing velocity, and others that it increases. In the speed ranges
possible with the Adv3 it seems `Kf` drops with increasing velocity.

After much fiddling around and experimenting, I've formed a theory that this
is because extruding against the print surface generates back pressure. The
nozzle is extruding into a bead of cooling filament that is pressed against
the print surface. Below a certain pressure threshold it cannot push past that
bead and the flow is effectively zero. So what actually drives the flow rate
at low pressures is nozzle movement, moving away from the flow-blocking bead.
At a constant head speed and with sufficient feed rate it will draw a pretty
consistent line, with the bead smearing behind he nozzle. The advance pressure
doesn’t make much difference because the flow rate is dominated by the
back-pressure and nozzle speed moving away from that back-pressure.

At higher nozzle speeds and flow rates, the pressure for the flow rate through
the nozzle starts to dominate. I suspect the back-pressure might be pretty
constant at any speed for a given material, layer height, track width, and
extrusion ratio, because the smearing bead shape will be pretty much the same
for the same track shape regardless of the speed. It might drop a bit as speed
increases because the smearing bead has less time to cool and will flow
better, but i’m not sure how much difference that would make. I don’t think
the bead cools much under the nozzle at even super slow speeds. that is
relieved by moving the print head, smearing the bead away. This would explain
why `Kf` appears to decrease with velocity.

At even higher speeds I suspect the linear [Poiseuille's
equation](https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation)
assumption starts to fail and [Bernoulli’s
equation](https://en.wikipedia.org/wiki/Bernoulli%27s_principle) with pressure
as a function of velocity squared starts to dominate, which could explain why
some people are seeing `Kf` increase with velocity.

We will assume the back-pressure is a function of only the bead diameter,
which is the track width (and extrusion ratio), and is independent of the
nozzle velocity. Note this could be very wrong, with faster speeds giving the
bead less time to cool and harden, lowering the back-pressure, but we assume
this for now. We also assume that the back-pressure for a given material and
layer height is a linear function of the bead diameter. I could, and have,
invent all sorts of rationalizations for this assumption but I'm basicly
making stuff up so I won't include it here. We'll assume this for simplicity
for now.

The flow through the nozzle also requires pressure. Although I suspect that
this should be Bernoulli’s equation, the existing PA model assumes
Poiseuille’s equation, so for now we will too.

This gives us;

```python
  Pe = Pn + Pb      # (1) total extruder advance pressure.
  Pn = Kf*vn        # (2) nozzle pressure from flow rate.
  Pb = Kb*Db        # (3) back pressure from extruded bead.
  vn = v*Db*h/Af    # (4) nozzle flow velocity.
  Db = r*w          # (5) extruded bead diameter.
  Af = pi*(Df/2)^2  # (6) cross section area of filament.
```

Where:

```
  Pe is the extruder pressure in mm of filament advance.
  Pn is the nozzle flow pressure in mm of filament advance.
  Pb is the bead back-pressure in mm of filament advance.
  Kf is the nozzle pressure factor.
  Kb is the bead backpressure factor.
  vn is the nozzle outflow velocity in mm/s of filament.
  ve is the extruder velocity in mm/s of filament.
  Df is the filament diameter in mm.
  Db is the bead diameter in mm.
  Af is the filament cross-section area in mm^2.
  r is the extrusion ratio.
  w is the nominal line width in mm.
  h is the layer and line height in mm.
  v is the nozzle velocity in mm/s
```

Note for a constant track width at steady state where extruder velocity `ve`
equals nozzle flow velocity `nv` we get;

```python
  Pe = Kf*ve + Cf
  Cf = Kb*Db
```

This is a Linear Advance model with a constant offset that depends on the line
width (and material, and nozzle temp, and ...), which so far appears to match
what I've seen.

Note in earlier work on this I included a `Cb` constant offset for the
backpressure, but I've adbandoned it because it makes no difference to, and is
impossible to distinguish from, additional `Re` retraction.

To validate and calibrate this model we need to measure `Pe` to see how it
varies with other variables. Measuring `Pe` can be done by pre-applying
different `Pe` values with a restore before printing a line, and seeing which
Pe gives a good consistent line, or applying different retract values after a
printed line has stabilized and seeing how much retraction is required to
relieve the pressure.

Note that at very low speeds `Pn` should be very small, so we can almost
directly measure `Pb` without it being significantly impacted by `Pn`,
assuming it doesn't vary with speeds. At low speeds we should be able to get
several start/stop cycles per line with different retract/restore distances.

### Backpressure effects for overhangs, bridging, and middle-walls.

For overhangs and bridging I'd expect the backpressure to be reduced by the
percentage overhang, so for bridging it would be zero. Note Orca has settings
for different pressure-advance for bridging.

### Backpressure Testing

We need to measure the `Pb` values to see how they vary with;

* print speeds: 20-100mm? 10mm? 1mm? Note higher speeds will also be
  affected by `Kf` and nozzle flow rates, and it will be hard to separate these
  affects. However, it's also very possible that print speeds affect `Pb` in
  ways that cannot just be included with the linear `Kf` affects, and make `Pb`
  not purely dependent on line width.
* layer heights: Try values 0.1,0.2,0.3,0.4?
* track widths: 0.3,0.4.0.5,0.6,0.7,0.8?

Can we measure `Pb` on a single line using the RetractTest technique of doing
a moving retraction and measuring the smear length? We need to make sure the
test is run without `-P` so it doesn't pressure-compensate in the deceleration
of the draw. The RetractTest1 results below show that the `vr` retracting-move
speed will have an impact on the measurement. I think 10mm/s is probably a
sweet spot between too fast to get a readable track and so slow it just drools
forever. Any pressure too low to leave a track at 10mm/s can be considered
zero. Note drool with an average thinkness of half a h=0.3 x w=0.4 x r=1.0 is
about 0.025mm of filament per mm, and we should include that in the
end-of-draw pressure.

The other problem is should we try to measure the restore pressure needed to
start the line? The StartStopTests measure restore and retract pressure by
trying a range of retract/restore lengths so the ideal pressure can be found
by picking the best line. Unfortunately this requires 10 lines to figure out
the right pressure, instead of the single line the RetractTest uses.

So I think we use a RetractTest to get the initial readings, and then maybe
validate it with a StartStopTest.

If we keep `ve` constant we can avoid nozzle flow rate variations. We can do
this by varying the line dimensions while keeping the line area and print
speed constant. We can also vary the line area and print speed to keep the
extrusion rate constant.

For all 4 tests we should use;

* Use args `-Kf=0.4 -Kb=3.0 -Cb=0.0 -Re=1.0 -R`. These settings seem to be
   pretty close and give good retract/restore pressure-handling to help ensure
   there is no pressure accumulation artifacts between tests and give good
   ruler/brim/etc results.

* Change the RetractTest to include the StartStopTest initial line
   prime/drain/restore steps to ensure the lines start without carryover
   pressure artifacts from earlier tests and get a good line start without
   `-P` enabled.

* Extend the retraction distance to 40mm, removing the 20m move. The move was
  originally added to ensure there was no deceleration at the end of the
  retracting move to mess up the drool length measure. However, at 10mm/s the
  deceleration is so short it's negilgable.

* Only retract 4mm, making each 10mm of ooze be 1mm of retraction, noting that
  each 10mm of line is also about `0.6*0.2*10/Fa=0.5mm` of fillament, so
  assuming the average ooze width is half that, each 10mm of ooze is about
  1.25mm of total ooz+retract Pe.

* Use a fixed ve where possible for the tests to avoid nozzle flow rate
  differences.

For the 4 tests use;

1. Constant `vx`, constant `ve`, to test how `Pe` varies with `w` and `h`. We
  expect `Pe` to vary linearly with `w`;

  * `w=(0.3, 0.8)` for 10 steps every 0.05mm width.
  * `h=0.09/w` giving h=0.3 to 0.1125.
  * `vx=10`
  * `ve=0.3*0.3*10/Fa` = 0.3742mm/s.
  * Expected Pn=0.1, Pb=0.9 to 2.4, Pe=1.0 to 2.5.

1. Constant `w`, constant `ve` to test that `Pe` doesn't vary with `h` and
  `vx`. We expect that `Pe` should be the same for all these lines;

   * `w=0.5`
   * `h=(0.1, 0.35)` for 10 steps every 0.025mm.
   * `vx=3.5/h` giving vx=10 to 35mm/s.
   * `ve=0.5*0.35*10/Fa` = 0.7276mm/s.
   * Expected Pn=0.3, Pb=1.5, Pe=1.8 mm for all tests.

1. Constant `h`, constant `ve` to test how `Pe` varies with `w` and `vx`. We
  expect `Pe` to vary linearly with `w`;

   * `w=(0.3, 0.8)` for 10 steps every 0.05mm width.
   * `h=0.2`
   * `vx=8/w` giving vx=26.7 to 10mm/s.
   * `ve=0.2*0.8*10/Fa` = 0.6652mm/s
   * Expected Pn=0.3, Pb=0.9 to 2.4, Pe=1.2 to 2.7 mm.

1. Constant 'h', constant 'w' to test how `Pe` varies with `vx` and `ve`. We
  expect `Pe` to increase linearly with `ve`.

   * w=0.5
   * h=0.2
   * vx=(10,60) for 10 steps every 5mm/s
   * ve=0.2*0.5/Fa*vx = 0.0416*vx giving ve=0.416 to 2.495 mm/s.
   * Expected Pn=0.2 to 1.0, Pb=1.5, Pe=1.7 to 2.5 mm.

The results of this are in [Retracttest2](#RetractTest2) below.


# Printer Tolerance

## Mesh approximation affects

Note with outer walls CAD programs typically export by creating an approximate
polygon mesh out of points on the model surface. There is usually a tollerance
resolution setting (onshape calls it "Chordal tolerance") which is the max
distance between the model surface and the closest polygon surface. This means
that for concave surfaces (eg holes), the mesh "fills over curves" compared to
the model, and for convex surfaces (eg cylinders) it "cuts off curves", by as
much as the tolerance setting. For fitting parts, you care about the maximum
dimensions for convex surfaces and the minimum dimensions for concave
surfaces. With convex surfaces the maximum dimensions are the points, which
are correct and don't need adjustment, but for concave surfaces the minimum
dimension is the middle of the polygons which is too small by the tolerance
distance and needs to be adjusted outwards by that much. So the diameter of
cylinders should not need adjusting, but the diameter of holes will need to be
increased by 2*tolerance. For OnShape "fine" settings tolerance=0.06mm, so the
diameter of holes need to be increased by 0.12mm.

## Slicer approximation affects

When the slicer converts the mesh into layers, another similar approximation
is applied which can further reduce hole sizes. This is controlled by the
slicer precision setting (Orca has "Quality; Resolution; Precision"). This is
approximating polygon cross-sections with simpler (less sides) polygons that
deviate by less than precision from the cross-sections. This can reduce
hole-sizes similar to the mesh approximation, and shouldn't adversely affect
external contours much. However, in this case the max deviation would be from
the cross-section polygon points when the simplification "cuts a polygon point
off", and the mesh points are accurately on the model surface. Note the
slice-polygon points correspond with a mesh edge, and not necessarily with a
mesh point, but they are closer to the real model surface than the middle of
the polygon slice lines. So the larger out of precision and tolerance has the
most effect, but the smaller can also contribute a bit.

A rough guess of the combined reduction in hole diameter would be
`2*sqrt(precision^2 + tolerance^2)`. Note that setting precision super small
will make the slicing approximation closer and closer to the mesh
approximation, and cannot go finer. So in practice setting precision to
something really small like 0.01mm means your quality and gcode size is bound
by the tolerance used to export the model. For precision=0.01mm and
tolerance=0.06mm the precision effect is negligable and the diameter increase
needed for holes is still about 0.12mm.

### Line shape compensation

https://manual.slic3r.org/advanced/flow-math
https://github.com/SoftFever/OrcaSlicer/issues/8313#issuecomment-3082244171

This attempts to take into account that the line cross-section is not
rectangular, but bulges out the sides giving a wider line-width. It then
"overlaps" the bulges between adjacent lines to fill the gaps. Unfortunately
after all this maths feedback shows the required overlap is to fill the gaps
100%, which means it's equivalent to just using a simple rectangular cross
section.

At least everywhere except the outer wall. For the outer wall it does give you
a perimiter that tries to match the outer-edge of the layer-bulges to the
perimiter, which could be good for non-interference fits. However, for other
applications the average outer-edge is a more useful solution.

In orca this manifests as the line-spacing (which is the effective line width)
being less than the line-width setting by an amount that varys with the layer
height. So effective line-widths are `layer_height * (1 - PI/4)` or roughly
`layer_height/5` narrower than the line-width setting used. It also moves the
outer wall inward by `layer_height * (1-PI/4)/2` or roughly `layer_height/10`
compared to "rectangular" slicers.

### Polyholes

Orca also has "convert holes to polyholes", which is an interesting idea
because it attempts to correct the "holes are smaller" problem by using a
polygon where the middle of the edges is on the model boundary instead of the
points. However, in this case the hole is already a mesh approximation of the
real model boundry, and I think polyholes also has been tuned to try and also
compensate for that. It does some strange things like reducing the number of
polygon edges as the hole gets smaller, which effectively reduces the amount
extruded at the hole edge in a way that is a bit hard to analyse and possibly
varies in hard-to-characterize ways with other dimensions like layer-height.

### Arc Compensation

https://reprap.org/wiki/ArcCompensation

The theory is that curved lines extrude even volumes on the inside and outside
of the line, but because the inside is a shorter path than the outside, it
ends up over-extruded on the inside and under extruded on the outside,
resulting in an overall reduction in the arc radius. This affect depends on
the hole radius, and increases as the hole gets smaller.

Even in theory this affect is usually small, with a hole diameter larger than
13mm not needing compensation, and a 2mm hole only needing a 0.06mm increase
in diameter, but does get significant for smaller holes, with a 1mm diameter
hole needing to be increased by 0.48mm.

In practice I suspect the backpressure effect causes some of the inside
over-extrusion to be pushed back inside, reducing the affect, and probably
quite significantly. With backpressure linearly proportional to bead width, I
would guess the affect is at least halved.

The estimated required radius `r` or diameter 'd' for a hole of radius `R` or
diameter `D` and line-thickness `t` is;

```python
  r = (t + sqrt(t^2 + 4*R^2))/2
  d = (t + sqrt(t^2 + D^2))
```

Another possible related effect is "corner pulling" where the filament is
pulled inwards on a bend by (cooling?) tension.

### Thermal expansion/shrinkage effects

As extruded fillament cools it shrinks, which will tend to reduce outside and
inside dimensions. I'm not sure how much difference this makes. Note that
extruded volumes are normally measured in input fillament volume which is
cool. This means the thermal expansion effects on the extruded fillament are an
over-extrusion while hot that shrinks back to the desired volume when cool,
So in theory it should kinda cancel out.

Bed and enclosure temperatures can also have an effect, as these mean the
filament doesn't cool all the way back to ambient until after the print is
finished and removed. Perhaps on those cases it should be taken into account?

Exceptions are special filaments that have some permanent expansion after
heating, like foaming filaments, or perhaps some with tempering affects?

### Arc Fitting

https://github.com/SoftFever/OrcaSlicer/wiki/quality_settings_precision#arc-fitting

This attempts to map sliced polygons back into circular arcs. This is of
questionable value because;

* Not many printers support the G2 and G3 arc move commands.

* Even printers that support G2 and G3 ususally implement them by
  approximating the curve with straight lines, which might have worse
  resolution than the polygons the arcs were approximated from.

* Approximating polygons back to arcs can get it wrong; maybe the original
  model didn't even have curves? Even if the original model did slice to arcs,
  the mesh approximation can distort and shift them so the generated arc
  doesn't match them.

* Most curved surfaces don't project into circular arcs when sliced. Even
  circular holes become elipses when cut diagonally.

However it could be useful if;

* When you have a very course mesh, this might do a good job of improving the
  resolution of the output.

* It can significantly reduce gcode size while retaining resolution for a very
  fine resolution mesh.

* It lets the printer use it's own optimal curve resolution instead of
  under/over doing it in the gcode.

* If the printer's internal resolution is higher than the mesh-model then arcs
  will probably be closer to the original model dimensions.

### layer slice height tweaks

I can't remember the details, but I read that some slicers have settings that
let you adjust the slice-plane height between the top/middle/bottom of the
layer depending
on the wall slope.

Normally the middle of the layer is used as the slice-plane
for calculating the edge polygons, but for sloping sides there is a
"stair-case" effect that gives outermost and innermost edges that deviate from
the slope, increasing the outer-edges, reducing internal hole sizes and
increasing external
contours.

By using the the layer bottom height as the slice plane you shift the
outer-wall in for outward-sloping surfaces, and out for inward sloping
surfaces, and the opposite when using the layer top height, where the amount
shifted depends on the slope.

Slicers with this setting usually offer the option of using the combination
that shifts the layers inwards, improving fits for non-vertical holes.
