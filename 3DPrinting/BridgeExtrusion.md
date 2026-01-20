# Bridge Extrusion Maths

This defines a model for FDM printing's extrusion lines for both normal
printed lines and bridge lines. They are based on the model used by OrcaSlicer
that comes from Slic3r, but extends it in a few ways that can be used to
improve OrcaSlicer's bridging implementation, and shows how we can map it as
close as possible to the current OrcaSlicer settings.

## References

This is based on ideas and discussions at;

* https://manual.slic3r.org/advanced/flow-math
* https://github.com/OrcaSlicer/OrcaSlicer/pull/11255

There is also the following supporting documents

* [spreadsheet for translating to/from OrcaSlicer](https://docs.google.com/spreadsheets/d/1q8VEJtvlpofTINNaKx3oTjjiWdCzXQhwqQcuNmZ9_eM/edit?usp=sharing)
* [Onshape CAD model](https://cad.onshape.com/documents/d6c789081e4718027433abc0/w/fab8c0f8c830995797647553/e/4c161cfdaa646b1caebb929e?configuration=Advanced_Settings%3Dtrue%3BAfw%3D1.0%3BAtl%3D4.0E-4%2Bmeter%3BBd%3D3.8E-4%2Bmeter%3BBdr%3D0.95%3BBfh%3D1.0%3BBfp%3D0.0%3BBfw%3D1.0%3BBh%3D3.37E-4%2Bmeter%3BBs%3D3.37E-4%2Bmeter%3BDn%3D4.0E-4%2Bmeter%3BLfp%3D0.0%3BLfw%3D1.0%3BLh%3D2.0E-4%2Bmeter%3BList_10u8p1SlZapMRd%3DDefaults%3BList_vwkf0slImDvp20%3DSet_Line_Width%3BLs%3D4.0E-4%2Bmeter%3BLw%3D5.0E-4%2Bmeter&renderMode=0&uiState=6960a1006d2f50394174a759)

## Intro

There are three different kinds of lines related to bridging;

* Normal perimeter and solid fill lines the bridge sits on and under.

* Bridge lines that span the bridge. These lines are not printed squished onto
  a lower layer, but string through space with minimal suport from adjacent
  lines.

* Anchor lines that sit on top of perimeter lines at the ends of the bridge
  lines. Since they are extensions of bridging lines, they share the same
  `line_spacing`. However, because they are printed onto a supporting layer,
  they are printed squished, so still need different line characteristics.

The transition point between Anchor and Bridge lines requires changing the
line characteristics and print height (as explained below). Sudden line
characteristic changes and hops can have undesirable artefacts from print
speed and pressure changes. Additionally, the transition from printing onto a
supporting surface and printing over free-air alone can have significant
pressure change artefacts. For this reason the transition between the
different line characteristics and print height should probably be done
gradually over a small distance on the anchor-side of the transition edge.

## Settings and parameters

In OrcaSlicer, line cross-sections are based on the [flow
maths](https://manual.slic3r.org/advanced/flow-math) that models the extruded
line cross-section as a rounded-ended rectangle, which is used to define the
`line_width` and `line_area`, but the lines are packed together with
`overlap_factor=1`, which gives a `line_spacing` where `line_spacing =
line_area / layer_height`. For the outer-wall lines the external perimeter is
at the outer-most-edge of the rounded-edge "bulge", so there is no perimeter
overlap, or `perimeter_factor=0`.

So OrcaSlicer doesn't (yet) have exposed settings for `overlap_factor` which
is effectively hard-coded internally to 1, or for `perimeter_factor`, which is
effectively internally hard-coded to 0. For completeness we include these
factors in the model, but note they cannot (yet) be changed in OrcaSlicer from
the defaults.

For bridge lines, they are not squished against a layer below so we model them
as circular sections that overlap with adjacent bridge (or perimeter) lines,
and the layer above. The `line_diameter` defines both the `line_width` and
`line_thickness`, and is typically the same or slightly smaller (due to
stretching) than the nozzle diameter. Trying to use larger diameters typically
results in sagging instead (what about foaming fillaments?).

In addition to rounded sides that overlap horizontally, they have rounded top
and bottom surfaces which overlap with the layer above and bulge towards the
external perimeter below. This makes them like normal perimeter lines, where
the rounded bulge makes the external perimeter distance ambiguous. They also
need to be printed from slightly above the layer top surface height to account
for the top overlap. The `bridge_height` is the vertical equivalent of
`line_spacing`, with the same `overlap_factor` on the top and bottom, so it is
the theoretical `layer_height` of the line if printing a bridge line on top of
something, like perhaps a second bridge layer.

The `anchor_spacing` must be the same as `bridge_spacing` since they are
continuations of the bridge lines, and `anchor_height` must be the same as
`layer_height` since they are wedged between layer lines above and below, so
there are no separate values for those.

The parameters that define these lines use the following input settings. Note
either `line_width` or `line_spacing` can be set, and the one that is set will
be used to derive the other.

* Printer and layer settings
   * Dn - `nozzle_diameter` (range: [0.1mm,1.0mm], default: 0.4mm)
   * Lh - normal `layer_height` (range: [0.1mm,1.0mm], default: 0.2mm)
* Normal line settings
   * Lw - normal `line_width` (range: [Dn/2+Lh,Dn*4+Lh], default: Dn+Lh)
   * Ls - normal `line_spacing` (range: [Dn/2,Dn*4], default: derived from Lw)
   * Lfw - normal side `overlap_factor` (range: [0,1], default: 1)
   * Lfp - normal side `perimeter_factor` (range: [0,1], default: 0)
* Bridge line settings
   * Bd - bridge `line_diameter` (range: [Dn/2,Dn], default: Dn)
   * Bh - bridge `layer_height` (range: [Dn/2,Dn], default: derived from Bd)
   * Bs - bridge `line_spacing` (range: [Dn/2,Dn], default: derived from Bd)
   * Bfw - bridge side `overlap_factor` (range: [0,1], default: 1)
   * Bfh - bridge vertical `overlap_factor` (range: [0,1], default: 1)
   * Bfp - bridge bottom `perimeter_factor` (range: [0,1], default: 0)
* Bridge anchor line settings;
   * Afw - anchor `overlap_factor` (range: [0,1], default: Bfw)
   * Atl - anchor `transition_length` (range:[0mm,10mm], default: Ls)

And from these settings the following attributes can be derived;

* Printer and layer settings
   * An - `nozzle_area` or `Dn^2 * pi/4`.
* Normal line attributes
   * Lp - normal `perimeter_offset` from line center.
   * La - normal `line_area`
* Bridge line attributes
   * Bl - height/width of a square with the same area as Ba.
   * Bz - bridge `print_offset` distance from bridge-layer top up to print height.
   * Bp - bridge `perimeter_offset` from external perimeter up to print height.
   * Ba - bridge `line_area`
* Anchor line attributes
   * Aw - anchor `line_width`
   * Aa - anchor `line_area`

Note that any `line_area` settings is is also the extrusion rate in mm^3/mm.
This can be calculated from OrcaSlicers reported flow rates in mm^3/sec by
dividing it by the reported speed in mm/sec.

## Line calculations

Line characteristics can be specified using `layer_height`, `overlap_factor`,
and either `line_width` or `line_spacing`. The relationship between
`line_width`, `line_spacing`, and `line_area` for normal and anchor lines is;

```
Ls = Lw - Lfw*Lh*(1-pi/4)      ; Ls from Lw, Lh, and Lfw
Lw = Ls + Lfw*Lh*(1-pi/4)      ; Lw from Ls, Lh, and Lfw

La = Lh*(Lw - Lh*(1 - pi/4))
Lp = (Lw - Lfp*Lh*(1-pi/4))/2
```

Note OrcaSlicer currently hard-codes `Lfw=1.0` and `Lfp=0.0` which simplifys
this to;

```
Ls = Lw - Lh*(1-pi/4)      ; Ls from Lw, Lh, and Lfw
Lw = Ls + Lh*(1-pi/4)      ; Lw from Ls, Lh, and Lfw

La = Lh*(Lw - Lh*(1 - pi/4)) = Lh*Ls
Lp = Lw/2
```

But also note that OrcaSlicer has other settings for some types of lines
modifying the "density" which will adjust the spacing without adjusting the
width. When it does this the area is always calculated from the width.

## Bridge Calculations

For bridges, overlap is calculated assuming that a nominal `line_spacing` with
horizontal and vertical `overlap_factor=1` is equivalent to a square
cross-section with the same area as the circular extrusion. This gives the
following relationships between settings and attributes.

```
Bl = Bd*sqrt(pi/4)                ; Bl from Bd.
Bs = Bd*(1 - Bfw*(1-sqrt(pi/4)))  ; Bs from Bd and Bfw
Bh = Bd*(1 - Bfh*(1-sqrt(pi/4)))  ; Bh from Bd and Bfh
Bd = Bs/(1 - Bfw*(1-sqrt(pi/4)))  ; Bd from Bs and Bfw
Bd = Bh/(1 - Bfh*(1-sqrt(pi/4)))  ; Bd from Bh and Bfh

Bfw = (1 - Bs/Bd)/(1-sqrt(pi/4))  ; Bfw from Bs and Bd.
Bfh = (1 - Bh/Bd)/(1-sqrt(pi/4))  ; Bfh from Bh and Bd.

Ba = Bd^2*pi/4
Bz = (Bd - Bh)/2 = Bd*Bfh/2*(1-sqrt(pi/4))
Bp = Bd*(1 - Bfp/2*(1-sqrt(pi/4)))
```

Since `line_diameter` interacts with both `line_spacing` and `layer_height`
the interaction between these complicates which options can be set vs derived.
Bridge quality is highly dependent on the `line_diameter` and
`overlap_factor`, so it's recommended that `line_diameter` be set with
`line_spacing` and `layer_height` being derived. However all the following
setting combinations could also work, but it would require careful limiting to
keep all settings within reasonable limits;

1. Derive `line_spacing` and `layer_height`. Set Bd, Bfw, Bfh. Derive Bs, Bh.
1. Specify `layer_height` and derive `line_spacing`. Set Bd, Bh, Bfw. Derive Bs, Bfh.
1. Specify `line_spacing` and derive `layer_height`. Set Bs, Bfw, Bfh. Derive Bd, Bh.
1. Specify `line_spacing` and `layer_height`. Set Bs, Bh, Bfw. Derive Bd, Bfh.

## Anchor Calculations

The anchor lines are extensions of the bridge lines, but they are printed over
supporting perimeter and infill lines, so they are printed squished like
normal lines and require different line characteristics to the unsupported
bridge lines. Without this there will be nasty artifacts at the start and end
of the bridge lines.

Anchor line spacing is the same as bridge-spacing `Bs`, and `layer_height` is
the same as normal lines `Lh`.

```
Aw = Bs + Afw*Lh*(1-pi/4)         ; Aw from Bs, Lh, and Afw
Aa = Lh*(Aw - Lh*(1 - pi/4))      ; Aa from Lh and Aw.
Aa = Lh*(Bs - (1-Afw)*Lh*(1-pi4)) ; Aa from Lh, Bs, and Afw.
Afw = (Aw - Bs) / (Lh*(1-pi/4))   ; Afw from Aw, Bs, and Lh.
```

Note these need to be printed with `print_height` at the normal print height
of the layer, which is below the `print_height` for the bridge lines they are
extensions of. This means there should be a little hop-up at the transition
from anchor to bridge line.

## Pressure Advance

A significant factor in Pressure Advance is back-pressure from the print
surface. Bridge lines are not printed against a surface, which means the
pressure advance required is much less.

However, the anchor lines **are** printed against the supporting lines, which
means more Pressure Advance is required. In practice this means pressure will
build up while the anchor lines are being printed, and then be suddenly released
when the line moves past the edge of the supporting lines and starts bridging.
This will contribute extra slight over-extrusion at the start of each bridge
line that will decay away until the pressure stabilizes. This is probably a
contributing factor for the staggered line ends affect highlighted at the end
of this video. You can even see the slighly larger diameter of the
lines leaving the support lines vs approaching in the staggered cross-section
vs the more uniform diameters in the smooth central bridge;

https://www.youtube.com/watch?v=Mrs2kAuRCBk

This problem is made even worse by OrcaSlicer using the same extrusion rate
for the anchor lines as the bridging lines, despite the completely different
cross section area and thus flow rate required for those lines. This usually
means that extrusion rates and line spacing that works perfectly for the
bridge lines is hugely over-extruded for the anchor lines. For the short
anchor lines at the beginning and end of the bridges this usually just
contributes more excess pressure that is relieved as more over-extrusion at
the start of the bridge lines. However OrcaSlicer also often has extra bridge
lines at the sides of the bridges that overlap supporting perimeters and are
thus actually anchor lines for their whole length. For those the
over-extrusion can really build up and make a mess, as seen in this video;

https://www.youtube.com/watch?v=eaasEkFULKE

The only way to minimize the problems caused by different bridge and anchor
extrusion rates for slicers like OrcaSlicer that don't differentiate between
them is to pick settings that produce bridge and anchor lines with almost the
same cross section area. Assuming anchor overlap-factor Afw=1 is reqired to
maximise anchor area and the optimal bridge overlap-factor is `Bfw=1`, the
bridge and anchor areas are the same when `Bs=Lh`.

This is normally significantly less than the nozzle diameter, which might be
OK for very stringy materials, but could cause "snapping" for some materials.

## Thoughts.

The bonding strength between adjacent lines and layers probably depends more
linearly on the contact length than the overlap distance. This suggests it
might be better to use contact-factor settings that map to `overlap_factor`
settings. The maths for this could be done by assuming the overlapped profile
is a chamfered-rectangle with the same area as the nominal profile, with the
flat contact faces being the contact length. A more complicated model might
assume bridge lines have eliptical chamfers where the elipse distances depend
on the overlap factors for vertical vs horizontal overlaps.

However, I'm not sure how much difference that would make, and is probably not
worth it. It also will not translate easily to `overlap_factor`s slightly
beyond the [0,1] range that might be useful for deliberate extra spacing or
overlap.

I need to figure out exactly how the best bridging from empirical testing
results translates into overlap. It's possible that the square for
`overlap_factor=1` assumption doesn't match the best settings in practice, in
which case maybe another target rectangle would be better.


## OrcaSlicer v2.3.2-dev

The current implementation in OrcaSlicer doesn't use this model, and is a bit
messy. However, we can translate between settings for our model and equivalent
OrcaSlicer settings as closely as possible.

### Settings

OrcaSlicer distinguishes between internal and external bridges, and duplicates
settings for both. We will focus on external bridges, but this also applys for
internal bridges. It also has two different modes for different bridging models;

* "legacy" bridges models bridges as modified internal solid infill lines and
  aims to keep the same `line_spacing` or `layer_height`. It uses the rounded
  box cross-section and adjusts the height while keeping the same spacing when
  increasing flows, and adjusts the spacing while keeping the same height when
  reducing flows.

* "thick" bridges models bridges with a circular cross-section where the
  height and width are the same as the diameter. The diameter is adjusted to
  adjust the flow. The default spacing is a tiny constant `0.05mm` larger than
  the diameter

The relevant settings OrcaSlicer exposes are;

* Internal solid infill `line_width` - The bridge `line_width` in "legacy" mode.
  (default: `Dn`, range: `[Dn/2,Dn*5)`

* `thick_bridge` - a checkbox toggling "legacy" and "thick" bridge modes.
  (default: `N`).

* `bridge_flow` - multiplier for the flow (default:`1.0`, range: `[0,2]`).

* `bridge_density` - a percentage scaling factor for reduction in line-spacing.
  Higher values means closer spacing (default:`100%`, range: `[10%,120%]`).

### Implementation

The class `Flows` in `Flow.hpp` have lengths `m_width`, `m_height`,
`m_spacing`, `m_nozzle_diameter`, and a `bool m_bridge`, with other attributes
all derived by methods from these attributes. `Flow::mm3_per_mm()` returns the
area of a circle with diameter=width for bridge flows, otherwise calculates
the rounded-rectange area using width and height.

For thick bridges, `Flow::bridging_flow()` in `Flow.hpp` takes a diameter
and the nozzle_diameter, and creates a flow with `width = height = diameter`
and `spacing = Flow::bridge_extrusion_spacing(diameter) = diameter+0.05mm`.
This suggests that it was once believed briding worked best with small gaps
between the lines, but recent empirical testing suggests the opposite.

For legacy bridges `Flow::with_flow_ratio()` is used to modify the solid
infill flow by the setting `bridge_flow`, which calls
`Flow::with_cross_section()`, implemented in `Flow.cpp`. This increases the
height and keeps the spacing when increasing flows, and keeps the height and
decreases the spacing when decreasing flows provided width is greater than the
height. If the reduced width is less than the height it switches to
calculating height and width as the diameter of a circle with the target area.
Note that height changes when increasing flows also gives you tiny line-width
changes to keep the same line-spacing. Reducing flows to widths less than the
layer height is unlikely to produce decent bridges so we will ignore that in
the following conversions and assume it just reduces spacing.

In `LayerRegion.cpp`, `LayerRegion::bridging_flow(FlowRole role, bool
thick_bridge)` returns different flows for legacy vs thick bridge mode. For
thick bridges it returns `Flow::bridging_flow()` with a diameter scaled by
sqrt of the `bridge_flow` setting. For legacy bridges it calls
`Flow::with_flow_ratio()` to return a non-bridging solid infill flow adjusted
by the `bridge_flow` setting to give the required area.

For both bridge models, the `bridge_density` setting is then used to adjust
the line-spacing without changing the flow to `flow.m_spacing/bridge_density`.

### Observations

A table of the observed bridge characteristics for given settings is below.
Note for `thick_bridge=N`, `Bd` is calculated from the observed flow `Ba`. We
use `Lh=0.2` for the layer height, `Ls=0.5` for the solid-infill line-spacing
(from `Lw=0.543` line width).

|thick_bridge |bridge_flow |bridge_density  |Bd   |Bh   |Bw   |Bs    |Ba    |
|-------------|------------|----------------|-----|-----|-----|------|------|
| Y           |1.0         |100%            |0.40 |0.40 |0.40 |0.450 |0.1255|
| Y           |1.0         |110%            |0.40 |0.40 |0.40 |0.409 |0.1255|
| Y           |1.0         |120%            |0.40 |0.40 |0.40 |0.375 |0.1255|
| Y           |1.1         |100%            |0.42 |0.42 |0.42 |0.470 |0.1380|
| Y           |1.2         |100%            |0.44 |0.44 |0.44 |0.488 |0.1510|
| N           |1.0         |100%            |0.36 |0.20 |0.54 |0.500 |0.1000|
| N           |1.0         |110%            |0.36 |0.20 |0.54 |0.454 |0.1000|
| N           |1.0         |120%            |0.36 |0.20 |0.54 |0.417 |0.1000|
| N           |1.1         |100%            |0.37 |0.22 |0.55 |0.500 |0.1100|
| N           |1.2         |100%            |0.39 |0.24 |0.55 |0.500 |0.1200|
| N           |0.9         |100%            |0.34 |0.20 |0.49 |0.450 |0.0900|
| N           |0.8         |100%            |0.32 |0.20 |0.44 |0.400 |0.0800|

**NOTE:** although the bridge flows usually have a different layer-height than
the normal layer height, Orca doesn't do anything about adjusting the
extrusion heights to compensate for that. It always prints bridges with the
extruder at the same height as the rest of the layer, and it always prints the
bridge on the same layer regardless of it's `Bh` bridge height or how "low it
hangs" below the top of the layer. This means in practice that bridges
normally hang much lower than the external perimeter height.

It a lot of bridge testing done by others that I've seen they have been using
"legacy" bridges that require significant `bridge_flow` increases (order 1.4x)
combined with `bridge_density` increases (around 110% to 125% with modified
OrcaSlicer to bypass settings limits) to produce decent results. These tests
often don't mention the `layer_height` used, but as the modelling shows the
flow rates for legacy bridges are highly dependent on it, so the ideal
`bridge_flow` and `bridge_density` settings will vary significantly with
`layer_height`.

If testing is done using "thick" bridges, it uses a better model and the flow
rates are independent of `layer_height`. The flow rates are significantly
higher, and I guess ideal `bridge_flow` settings would be in the range
`[0.5,1.0]`. It's a bit unfortunate that by default it adds a small space
between bridge lines because it also means `bridge_overlap` settings need to
be higher to overcome that, and the 120% limit is likely too low for ideal
results.

Note that for both kinds of bridges `layer_height` will have a significant
impact on the over-extrusion of anchor lines. I suspect ideal settings for
bridges will always benefit from larger `layer_height` settings to reduce the
anchor over-extrusion.

### Settings Translation

From these observations and looking at the code we can derive the formulas for
translating between our model's settings and OrcaSlicer's. For both bridging
models, the following mappings always apply.

Because it doesn't do any print-height compensation to account for the
normally larger bridge layer-height, the some bridge line hight
characteristics are fixed. For any bridge model, once we've mapped the
settings to our model's `Ba` area and `Bs` spacing settings, the `Bd` diameter
and `Bfw` overlap factor can be calculated. This gives us;

```
Bd = sqrt(Ba/(pi/4))
Bfw = (1 - Bs/Bd)/(1-sqrt(pi/4))
Bfh = 0.0  ; no vertical overlap for bridges.
Bz = 0.0   ; no print-height offset for bridges.
```

With no bridge height compensation at all, we cannot adjust the bridge's
external perimeter height at all, which is normally lower than the bottom of
the current print layer, which is generally the target perimeter. However, for
a given `Lh`, `Bp`, and `Bfp` perimeter factor to define exactly how we define
the perimeter, we can calculate the "error" between the bottom of the layer
and the bridge's actual bottom perimeter.

```
Be = Bp - Lh  ; distance from layer bottom to bridge bottom perimeter.
   ~= Bd - Lh  ; exact for Bfp=0, but still within 5% of Bd for Bfp<=1.
```

**NOTE:** for `Bfp=0` (ie, no perimeter overlap so use lowest point of the
round bridge extrusion) this simplifies to just `Be = Bd - Lh`. Given that for
thick bridges typically `Bd` is near the nozzle diameter and `Lh` is near 1/2
the nozzle diameter, this means the perimeters below bridges are nearly a
whole layer lower than they should be. For legacy bridges the default settings
give a smaller diameter, but it still hangs below the layer bottom, and for
decent bridges you need to increase the flow until it has the same problem.

OrcaSlicer calculates `Bh` height and `Bw` width values for it's bridge lines
that are different for thick and legacy bridges. We include them in the
translation calculations for each model, but these have no impact on the
resulting extrusion, which only depends on `Ba` and `Bs`.

OrcaSlicer does not differentiate between bridge lines and anchor lines, so
anchors use the same extrusion rates as bridges. This means we cannot adjust
anchors, but we can calculate their settings.

```
Aa = Ba
Ah = Lh
Aw = Ba/Ah + Lh*(1-pi/4)
As = Bs
Afw = (Aw - As)/(Lh*(1-pi/4))
```

**NOTE:** the default thick bridge lines have an area equal the nozzle area,
but normal lines have an area closer to half the nozzle area. This means thick
bridge anchor lines are pretty close to 2x over-extruded. The default legacy
bridge lines do not have this problem until you increase the flow enough to
give otherwise decent bridges.

#### Thick external bridges

With Thick external bridges enabled, the bridge characteristics from the
settings `bridge_flow` and `bridge_density` are;

```
Ba = bridge_flow * An
Bh = Bd
Bw = Bd
Bs = (Bd + 0.05)/bridge_density
```

We can use this to calculate required settings from our target model's `Ba`,
`Bd` and `Bs` values with;

```
bridge_flow = Ba/An
bridge_density = (Bd + 0.05)/Bs
```

This gives the following table for Bs/Dn, bridge_flow, and bridge_density for
Bd/Dn and Bfw values;

|Bd/Dn |  Bfw  | Bs/Dn   | bridge_flow | bridge_density|
|-----:|------:|--------:|------------:|--------------:|
|1.00  | 1.000 | 0.88622 | 1.0000      | 1.2694        |
|0.95  | 1.000 | 0.84191 | 0.9025      | 1.2769        |
|0.90  | 1.000 | 0.79760 | 0.8100      | 1.2851        |
|1.00  | 0.549 | 0.93750 | 1.0000      | 1.2000        |

**NOTE:** for Bd/Dn=1.0 we need bridge_density=127%, but the valid settings
range is limited to 120%. As Bd/Dn drops, bridge_flow also drops, but
bridge_density increases even further. So in practice the highest
overlap-factor we can get without setting `Bd` larger than the nozzle diameter
is `Bfw=0.549` with `Bd = nozzle_diameter`, `bridge_flow=1.0`, and
`bridge_density=120%`.

#### Legacy external bridges

Legacy external bridges, use an adjusted solid infill line. In the following
`Ls`, `Lw`, and `La` refer to the solid infill line characteristics.

The translations to our model from OrcaSlicer settings are;

```
Bls = min(bridge_flow, 1)*Ls    ; flow-adjusted solid infill line spacing
Bh = max(bridge_flow, 1)*Lh
Bw = Bls + Bh*(1-pi/4)
Bs = Bls / bridge_density
Ba = La * bridge_flow = Bh * Bls
```

The translation to OrcaSlicer settings from our model is;

```
bridge_flow = Ba / La
Bls = min(bridge_flow, 1)*Ls    ; flow-adjusted solid infill line spacing
bridge_density = Bls / Bs
```

**NOTE:** the typical solid infill line has an area about half the nozzle area
but a line spacing about the nozzle diameter. For defaults `bridge_flow=1` and
`bridge_density=100%` this means the default bridge lines are thin and widely
spaced. However, increasing `bridge_flow` increases the area without
increasing the spacing, allowing you to achieve similar `Bd` bridge diameters
as thick bridges but with higher `Bfw` overlap factors. However, it does make
bridge line characteristics vary sigificantly with `Lh` layer height for the
same `bridge_flow` setting, requiring retweaking them for different layer
heights.
