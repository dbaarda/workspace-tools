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

