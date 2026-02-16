# Line and Bridge Extrusion Maths

This defines a model for FDM printing's extrusion lines for both normal lines
and bridge lines. They are based on the model for normal lines used by
OrcaSlicer that comes from Slic3r, but extends it in a few ways that can be
used to improve OrcaSlicer's bridging implementation, shows how we can map it
as close as possible to the current OrcaSlicer settings, and suggests changes
to OrcaSlicer to make it easier to set ideal normal and bridge line settings.

Things like print-speed, nozzle-temperature, and fan-speeds also have an
effect on normal and bridge line quality, but this document focuses on
modeling and optimizing the extrusion volume and spacing of the lines
indepenently of speeds, temperatures, and cooling.

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
  they are printed squished, so they have different line characteristics and
  ideally require different settings.

The transition point between Anchor and Bridge lines ideally requires changing
the line settings and print height (see below). Sudden hops, line
characteristic changes, and settings changes can have undesirable artefacts
from print speed and pressure changes. So it is worth trying to minimize the
changes and perhaps apply the setting changes gradually over a small distance.

## A Bridge Extrusion Model

In OrcaSlicer, normal line cross-sections are based on the [flow
maths](https://manual.slic3r.org/advanced/flow-math) that models the extruded
line cross-section as a rounded-ended rectangle, which is used to define the
`line_width` and `line_area`. However, the lines are packed together with a
hard-coded `overlap_factor=1`, which gives a `line_spacing =
line_area/layer_height`, which is equivalent to a simple rectangular model
with `line_width = line_spacing`. For the outer-wall lines the external
perimeter is at the outer-most-edge of the rounded-edge "bulge", so there is
no perimeter overlap, or `perimeter_factor=0`.

OrcaSlicer doesn't (yet) have exposed settings for `overlap_factor` which is
hard-coded internally to 1, or for `perimeter_factor`, which is internally
hard-coded to 0. For completeness we include these factors in the model, but
note they cannot (yet) be changed in OrcaSlicer from the defaults.

Bridge lines are not squished against a layer below so we model them as
circular sections that overlap with adjacent bridge (or perimeter) lines, and
also with the layer above. The `line_diameter` defines both the `line_width`
and `line_thickness`, and is normally set smaller than the nozzle diameter so
the extruded filament is "stretched" to prevent sagging.

In addition to rounded sides that overlap horizontally, bridge lines have
rounded top and bottom surfaces which overlap with the layer above and bulge
towards (or past) the external perimeter below. This makes them like normal
perimeter lines, where we can use a vertical `overlap_factor` to define the
top surface overlap with the layer above, and a vertical `perimeter_factor` to
define the bottom `perimeter_offset`. This gives us a `print_offset` printing
height slightly above the layer top surface height to account for the vertical
overlap with the layer above, and a `perimeter_offset` distance below the
print height where the bottom external perimeter is. From this we also get a
`perimeter_error` distance that the perimeter protrudes below the bottom of
the current layer.

The `bridge_height` is the vertical equivalent of `line_spacing`, with the
same `overlap_factor` on the top and bottom, so it is the theoretical
`layer_height` of the bridge accounting for overlaps if printing a bridge line
on top of something, like perhaps a second bridge layer.

The `anchor_spacing` must be the same as `bridge_spacing` since they are
continuations of the bridge lines, and `anchor_height` must be the same as
`layer_height` since they are wedged between layer lines above and below.

### Settings and parameters

The parameters that define these lines use the following input settings with
the suggested optimal ranges and default values. Note either `line_width` or
`line_spacing` can be set, and the one that is set will be used to derive the
other.

* Printer and layer settings
   * `Df` - printer `filament_diameter` (range: [1.5mm,3.0mm], default: 1.75mm)
   * `Dn` - printer `nozzle_diameter` (range: [0.1mm,1.0mm], default: 0.4mm)
   * `Lh` - current `layer_height` (range: [0.1mm,1.0mm], default: 0.2mm)
* Normal line settings
   * `Lw` - line `line_width` (range: [`Dn/2+Lh`,`Dn*2+Lh`], default: derived from `Ls`)
   * `Ls` - line `line_spacing` (range: [`Dn/2,Dn*2`], default: `Dn`)
   * `Lfw` - line side `overlap_factor` (range: [0,1], default: 1)
   * `Lfp` - line side `perimeter_factor` (range: [0,1], default: 0)
* Bridge line settings
   * `Bd` - bridge `line_diameter` (range: [`Dn/4`,`Dn`], default: `0.9*Dn`)
   * `Bs` - bridge `line_spacing` (range: [`Dn/4`,`Dn`], default: derived from `Bd`)
   * `Bfw` - bridge side `overlap_factor` (range: [0,1], default: 1)
   * `Bfh` - bridge vertical `overlap_factor` (range: [0,1], default: 1)
   * `Bfp` - bridge bottom `perimeter_factor` (range: [0,1], default: 0)
* Bridge anchor line settings;
   * `Afw` - anchor `overlap_factor` (range: [0,1], default: 1)

And from these settings the following attributes can be derived;

* Printer and layer settings
   * `Af` - printer `filament_area` or `Df^2 * pi/4`.
   * `An` - printer `nozzle_area` or `Dn^2 * pi/4`.
* Normal line attributes
   * `Ll` - width of a rectangle with the same `Lh` and area `La`.
   * `Lp` - line `perimeter_offset` from line center.
   * `Lv` - line `line_space` cross-section area space available.
   * `La` - line `line_area` cross-section area extruded.
   * `Lr` - line `line_ratio` extrusion ratio.
* Bridge line attributes
   * `Bh` - bridge `layer_height` accounting for vertical overlap.
   * `Bl` - height/width of a square with the same area as `Ba`.
   * `Bz` - bridge `print_offset` distance from print height to bridge's top
     surface.
   * `Bp` - bridge `perimeter_offset` distance from print height to the bridge's
     bottom perimeter.
   * `Bv` - bridge `line_space` cross-section area space available.
   * `Ba` - bridge `line_area` cross-section area extruded.
   * `Br` - bridge extrusion ratio.
   * `Be` - bridge `perimeter_error` distance the bridge's bottom perimeter
     protudes past the bottom of the layer.
* Anchor line attributes
   * `Al` - width of a rectangle with the same `Lh` and area `Aa`.
   * `Aw` - anchor `line_width`
   * `Av` - anchor `line_space` cross-section area space available.
   * `Aa` - anchor `line_area` cross-section area extruded.
   * `Ar` - anchor extrusion ratio.

**NOTE:** `line_area` and `line_space` settings are also the extrusion rate
and volume available for the line in mm^3/mm volume per line-length. The
`line_area` can be calculated from OrcaSlicers reported flow rates in mm^3/sec
by dividing it by the reported speed in mm/sec. The `line_space` can be
calculated by getting the distance between adjacent lines and multiplying it
by the layer height.

### Line calculations

First we define some useful constants.

```python
Ka = pi/4 = 0.785       # area ratio of a circle to its enclosing square.
Kda = 1-Ka = 0.215      # area ratio of the difference between a square and its enclosed circle.
Kl = sqrt(pi/4) = 0.886 # length ratio of a square's side to the same-area circle's diameter.
Kdl = 1-Kl = 0.114      # length ratio of difference between a circle's diameter and same-area square's side.
```

Line characteristics can be specified using `layer_height`, `overlap_factor`,
and either `line_width` or `line_spacing`. The relationship between all the
attributes for normal and anchor lines is;

```python
Ll = Lw - Lh*Kda               # Ll from Lw, and Lh
Ll = Ls - (1-Lfw)*Lh*Kda       # Ll from Ls, Lh, and Lfw
Ls = Lw - Lfw*Lh*Kda           # Ls from Lw, Lh, and Lfw
Lw = Ls + Lfw*Lh*Kda           # Lw from Ls, Lh, and Lfw

Lv = Lh*Ls                     # Lv from Lh and Ls
La = Lh*Ll                     # La from Lh and Ll
La = Lh*(Lw - Lh*Kda)          # La from Lh and Lw
Lr = La/Lv = Ll/Ls             # Lr from La and Lv or Ll and Ls
Lr = 1 - (1-Lfw)*(Lh*Kda/Ls)   # Lr from Lh, Ls, and Lfw.
Lfw = 1 - (1-Lr)/(Lh*Kda/Ls)   # Lfw from Lh, Ls, and Lr.
Lfw = (Lw - Ls)/(Lh*Kda)       # Lfw from Lh, Ls, and Lw.
Lp = (Lw - Lfp*Lh*Kda)/2       # Lp from Lw, Lh, and Lfp.
```

Note OrcaSlicer currently hard-codes `Lfw=1.0` and `Lfp=0.0` which limits and
simplifys this to;

```python
Ll = Ls
Ls = Lw - Lh*Kda           # Ls from Lw, and Lh
Lw = Ls + Lh*Kda           # Lw from Ls and  Lh

Lv = Lh*Ls = La
La = Lh*Ll = Lv
Lr = 1
Lp = Lw/2
```

**NOTE:** OrcaSlicer has extra `line_density` settings for some types of lines
which will adjust the `Ls` spacing without adjusting the `Lw` width or `La`
extrusion rate area. When it does this the `La` extrusion rate area is always
calculated from the unmodified `Lw` width. This means these `line_density`
settings are effectively a very coarse way of adjusting the `Lfw` overlap
factor. There are also `extrusion_ratio` settings that adjust the `La`
extrusion rate without changing any other setting, which is also equivalent to
ajusting the `Lfw` overlap factor.

## Bridge Calculations

For bridges, overlap is calculated assuming that a nominal `line_spacing` with
horizontal and vertical `overlap_factor=1` is equivalent to a square
cross-section with the same area as the circular extrusion. This gives the
following relationships between settings and attributes.

```python
Bl = Bd*Kl                       # Bl from Bd.
Bs = Bd*(1 - Bfw*Kdl)            # Bs from Bd and Bfw
Bh = Bd*(1 - Bfh*Kdl)            # Bh from Bd and Bfh
Bd = Bs/(1 - Bfw*Kdl)            # Bd from Bs and Bfw
Bd = Bh/(1 - Bfh*Kdl)            # Bd from Bh and Bfh
Bd = sqrt(Ba/Ka)                 # Bd from Ba
Bfw = (1 - Bs/Bd)/Kdl            # Bfw from Bs and Bd.
Bfh = (1 - Bh/Bd)/Kdl            # Bfh from Bh and Bd.

Bv = Bh*Bs                       # Bv from Bh and Bs.
Ba = Bd^2*Ka                     # Ba from Bd.
Ba = Bl^2                        # Ba from Bl.
Br = Ba/Bv                       # Br from Ba and Bv.
Br = Ka/((1 - Bfh*Kdl)*(1 - Bfw*Kdl)) # Br from Bfh and Bfw
Bz = (Bd - Bh)/2                 # Bz from Bd and Bh.
Bz = Bd*Bfh*Kdl/2                # Bz from Bd and Bfh.
Bp = Bd*(1 - Bfp*Kdl/2)          # Bp from Bd and Bfp.
Be = Bp - Bz - Lh
```

If we copy OrcaSlicer's lines and apply similar hard-coded values of `Bfw=Bfh=1`
and `Bfp=0` this simplifies to;

```python
Bl = Bs = Bh = Bd*Kl
Bd = Bs/Kl
Bv = Ba = Bs^2 = Bd^2*Ka
Br = 1
Bz = (Bd - Bh)/2 = Bd*Kdl/2
Bp = Bd
Be = Bd - Bz - Lh = Bd*(1+Kl)/2 - Lh
```

The `Be` bottom `perimeter_errror` above assumes that the print head is
hopped-up by `Bz` to account for the `Bfh` vertical overlap factor of the
rounded top of the bridge lines with the layer above. If the bridge line is
printed without any hopup, this is equivalent to having no vertical overlap,
or `Bfh=0`. We can hard code for `Bfw=1`, `Bfh=Bfp=0` instead which instead
simplifies to;

```python
Bl = Bs = Bd*Kl
Bh = Bd
Bd = Bs/Kl
Ba = Bs^2 = Bd^2*Ka
Bv = Bs*Bh = Bd^2*Kl
Br = Kl
Bz = 0
Bp = Bd
Be = Bd - Lh
```

**NOTE:** For bridge diameters near the typical nozzle diameter of `Bd=0.4mm`
and optimal `Bfh=1` vertical overlap, the optimal hopup height is near `Bz =
0.057*Bd = 0.023mm`. This is well below the minimum `dz=0.1mm` Z resolution of
most 3DPrinters, so in practice adding a hopup is not worth it unless you have
a printer with very fine Z resolution or a very large nozzle diameter and use
very large diameter bridge lines.

## Anchor Calculations

The anchor lines are extensions of the bridge lines, but they are printed over
supporting perimeter and infill lines, so they are printed squished like
normal lines and have different line characteristics to the unsupported bridge
lines. Without accounting for this there will be nasty artifacts at the start
and end of the bridge lines. Anchor line spacing is the same as bridge-spacing
`Bs`, and `layer_height` is the same as normal lines `Lh`.

```
Al = Aw - Lh*Kda              # Al from Aw and Lh.
Al = Bs - (1-Afw)*Lh*Kda      # Al from Bs, Lh, and Afw.
Aw = Bs + Afw*Lh*Kda          # Aw from Bs, Lh, and Afw.
Av = Lh*Bs                    # Av from Lh and Bs.
Aa = Lh*Al                    # Aa from Lh and Al.
Aa = Lh*(Aw - Lh*Kda)         # Aa from Lh and Aw.
Aa = Lh*(Bs - (1-Afw)*Lh*Kda) # Aa from Lh, Bs, and Afw.
Ar = Aa/Av = Al/Bs            # Ar from Aa and Av or Al and Bs
Ar = 1 - (1-Afw)*(Lh*Kda/Bs)  # Ar from Lh, Bs, and Afw.
Afw = 1 - (1-Ar)/(Lh*Kda/Bs)  # Afw from Lh, Bs, and Ar.
Afw = (Aw - Bs)/(Lh*Kda)      # Afw from Lh, Bs, and Aw.
```

If we set Afw=1 this simplifies to;

```python
Al = Bs
Aw = Bs + Lh*Kda
Av = Aa = Lh*Bs
Ar = 1
```

**NOTE:** these need to be printed with print height at the normal print
height of the layer, which for the ideal `Bfh=1` is below the `Bz` offset
print height for the bridge lines they are extensions of. This means ideally
there should be a little hop-up at the transition from anchor to bridge line,
and hop-down at the transition from bridge to anchor. However, in practice the
ideal hop height is well below the minimum Z resolution of most printers, so
it's not worth including.

## OrcaSlicer v2.3.2-dev

The current implementation in OrcaSlicer doesn't use this model, and is a bit
messy. However, we can translate between settings for our model and equivalent
OrcaSlicer settings as closely as possible.

### OrcaSLicer Available Settings

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
  the diameter.

The relevant settings OrcaSlicer exposes are;

* Internal solid infill `line_width` - The bridge `line_width` in "legacy" mode.
  (default: `Dn`, range: ???)

* `thick_bridge` - a checkbox toggling "legacy" and "thick" bridge modes.
  (default: `N` for external, `Y` for internal). Available for both internal
  and external bridges.

* `bridge_flow` - multiplier for the flow (default:`1.0`, range: `[0,2]`).

* `bridge_density` - a percentage scaling factor for reduction in
  line-spacing. Higher values means closer spacing (default:`100%`, range:
  `[10%,100%]` in v2.3.1, `[10%,120%]` in v2.3.2-dev).

**NOTE:** The current v2.3.1 OrcaSlicer limits bridge_density to <=100%, and
for v2.3.2-dev nightly builds to <=120%. This means it is impossible to get
thick bridge lines to even touch in the current release, and even for the
nightly it is not possible to get decent overlap for smaller bridge diameters.

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
between the lines, but recent testing suggests the opposite.

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

A table of the observed bridge characteristics for given settings with
v2.3.2-dev is below. Note `Bd` is calculated from the observed flow `Ba`
calculated from the flow rate in mm^2/sec divided by the print speed in
mm/sec.  We use `Lh=0.2` for the layer height, `Ls=0.5` for the solid-infill
line-spacing (from `Lw=0.543` line width).

|thick_bridge|bridge_flow|bridge_density|Bd  |Bh  |Bw  |Bs  |Bfw  |Ar  |
|------------|-----------|--------------|----|----|----|----|-----|----|
| Y          |1.0        |100%          |0.40|0.40|0.40|0.45|-1.10|1.40| 
| Y          |1.0        |114%          |0.40|0.40|0.40|0.40| 0.12|1.59|
| Y          |1.0        |120%          |0.40|0.40|0.40|0.38| 0.55|1.68|
| Y          |1.1        |120%          |0.42|0.42|0.42|0.39| 0.59|1.77| 
| Y          |0.4        |120%          |0.25|0.25|0.25|0.25| 0.02|1.00| 
| N          |1.0        |100%          |0.36|0.20|0.54|0.50|-3.53|1.00|
| N          |1.0        |120%          |0.36|0.20|0.54|0.42|-1.47|1.20|
| N          |1.4        |120%          |0.42|0.28|0.56|0.42| 0.12|1.68|
| N          |1.7        |120%          |0.47|0.34|0.57|0.42| 0.92|2.04|
| N          |0.6        |120%          |0.34|0.20|0.34|0.25| 0.84|1.20|
| N          |0.4        |100%          |0.23|0.20|0.24|0.20| 1.00|1.00|

**NOTE:** although the bridge flows usually have a different layer-height than
the normal layer height, Orca doesn't do anything about adjusting the
extrusion heights to compensate for that. It always prints bridges with the
extruder at the same height as the rest of the layer, and it always prints the
bridge on the same layer regardless of it's `Bh` bridge height or how "low it
hangs" below the top of the layer. This means in practice that default bridges
normally hang much lower than the external perimeter height.

For thick bridges this shows for defaults there is large gaps between bridge
lines (`Bfw<0`) and significant over-extrusion of anchor lines (`Ar>1`).
Increasing `bridge_density` above 114% is required to make them overlap, and
at 120% the overlap is still not significant. Also increasing `bridge_flow`
increases overlap a bit more but makes the bridge diameter larger than the
nozzle diameter (`Bd>Dn`). Increasing `bridge_flow` or `bridge_density`
increases anchor over-extrusion even more. Reducing `bridge_flow` to 0.4
removes the anchor over-extrusion (`Ar=1`), but makes the lines barely touch
even with `bridge_density` at 120%.

For legacy bridges the defaults have even larger gaps between bridge lines,
but there is no over-extrusion of anchors. Increasing `bridge_density` to 120%
reduces the gaps, but they are still very large. Also increasing `bridge_flow`
to 1.4 finally closes the gaps, and at 1.7 there is a healthy amount of
overlap, but the bridge diameter is larger than the nozzle diameter and there
is significant over-extrusion of anchors. Reducing `bridge_flow` to 0.6 also
gives healthy overlap with minimal anchor over-extrusion. With `bridge_flow`
reduced to 0.4 we can reduce `bridge_density` back to 100% to give
bridge-spacing equal to the layer height for "perfect" overlap and anchor
extrusion.

A lot of bridge testing done by others that I've seen have been using "legacy"
bridges with significant `bridge_flow` increases (order 1.4x) combined with
`bridge_density` increases (around 110% to 125% with modified OrcaSlicer to
bypass settings limits) to produce decent results. These tests often don't
mention the `layer_height` or internal solid infill `line_width` used, but as
the modelling shows the flow rates for legacy bridges are highly dependent on
them, so the ideal `bridge_flow` and `bridge_density` settings will vary
significantly with `layer_height` and internal solid infill settings.

If testing is done using "thick" bridges, it uses a better model and the flow
rates are more independent of `layer_height`. The flow rates are significantly
higher, and I guess ideal `bridge_flow` settings would be in the range
`[0.5,1.0]`. It's unfortunate that by default it adds a small 0.05mm space
between bridge lines because it also means `bridge_overlap` settings need to
be higher to overcome that, and the 120% limit is too low to get ideal
overlap.

For both kinds of bridges `layer_height` will have a significant impact on the
over-extrusion of anchor lines. Using larger `layer_height` settings will help
reduce anchor over-extrusion for larger bridge diameters. Using smaller
diameter bridges than is optimal for the bridge lines to avoid over-extrusion
of anchor lines will probably be the best compromise for overall bridge
quality.

Also interesting is external bridge anchor lines sit on top of 2 lines, one
outer wall and one inner wall, but internal bridge anchor lines sit on top of
a single internal solid infill line. This means the anchor lines are shorter
for internal bridges than external bridges, which would reduce the impact of
anchor line over-extrusion for internal bridges.

### Settings Translation

From these observations and looking at the code we can derive the formulas for
translating between our model's settings and OrcaSlicer's. For both bridging
models, the following mappings always apply.

Because it doesn't do any print-height compensation to account for the
normally larger bridge layer-height, some bridge line height characteristics
are fixed.

```python
Bfh = 0.0  ; no vertical overlap for bridges.
Bz = 0.0   ; no print-height offset for bridges.
```

With no bridge height compensation at all, we cannot adjust the bridge's
external perimeter height at all, which is normally lower than the bottom of
the current print layer, which is generally the target perimeter. However, for
a given `Lh`, `Bp`, and `Bfp` perimeter factor to define exactly how we define
the perimeter, we can calculate the `perimeter_error` between the bottom of
the layer and the bridge's actual bottom perimeter.

```python
Be = Bp - Lh  ; distance from layer bottom to bridge bottom perimeter.
```

**NOTE:** for `Bfp=0` (ie, no perimeter overlap so use lowest point of the
round bridge extrusion) this simplifies to just `Be = Bd - Lh`. Given that for
thick bridges typically `Bd` is near the nozzle diameter and `Lh` is near 1/2
the nozzle diameter, this means the perimeters below default thick bridges are
nearly a whole layer lower than they should be. For legacy bridges the default
settings give a smaller diameter, but it still hangs below the layer bottom,
and for decent bridges you may need to increase the flow until it has the same
problem.

OrcaSlicer calculates `Blh` `layer_height` and `Blw` `line_width` values for
its bridge lines that are different for thick and legacy bridges. It also
calculates a `Bls` nominal line spacing before applying `bridge_density` to
calculate the `Bs` actual spacing. We include them in the translation
calculations for each model, but these have no impact on the resulting
extrusion, which only depends on `Bd` and `Bs`.

OrcaSlicer does not differentiate between bridge lines and anchor lines, so
anchors use the same extrusion rates as bridges. This means we cannot adjust
anchors, but we can calculate their effective settings.

```python
Av = Bs*Lh
Aa = Ba
Al = Ba/Lh
Aw = Al + Lh*Kda
As = Bs
Afw = (Aw - Bs)/(Lh*Kda)
Ar = Aa/Av                     # Ar from Aa and Av
Ar = Bd*Ka/(Lh*(1 - Bfw*Kdl))  # Ar from Bd, Lh, and Bfw
Bd = Ar*Lh*(1 - Bfw*Kdl)/Ka    # Bd from Ar, Lh, and Bfw
```

For `Ar=1` anchor lines that are not over or under
extruded we need the following;

```python
Aa = Av            # required for Ar=1
Ba = Bs * Lh       # from Aa=Ba and Av=Bs*Lh
Bl^2 = Bl * Lh     # from Ba=Bl^2 and Bs=Bl for Bfw=1
Bl = Lh
Bd*Kl = Lh         # From Bl=Bd*Kl
Bd = Lh/Kl
```

This means to avoid overextruding anchor lines we need bridge diameters
slighly larger than the layer height. Note that this also has the nice
property that the perimeter error `Be=0` for the optimal `Bfh=1` vertical
overlap, and only `Be = Lh*Kdl/Kl = Lh*0.128` for `Bfh=0` no vertical overlap.

**NOTE:** the default thick bridge lines have an area equal the nozzle area,
but normal lines have an area closer to half the nozzle area. This means
default thick bridge anchor lines are pretty close to 2x over-extruded. The
default legacy bridge lines have the same area as solid infill lines so they
don't have this problem`.

#### Thick external bridges

With Thick external bridges enabled, the bridge and anchor characteristics from the
settings `bridge_flow` and `bridge_density` are;

```
Blh = Bd
Blw = Bd
Bls = Bd + 0.05
Ba = bridge_flow * An
Bd = sqrt(bridge_flow)*Dn
Bs = Bls/bridge_density
```

We can use this to calculate required settings from our target model's `Ba`,
`Bd`, `Bs`, and `Bls` values with;

```
bridge_flow = Ba/An = (Bd/Dn)^2
bridge_density = Bls/Bs
```

Since OrcaSlicer limits the `bridge_density` to 120%, there is a maximum
possible overlap factor `Bfw` for a given `Bd` or `Bs`.

```python
Bs = (Bd+0.05)/bridge_density  # Bs from Bd for a given bridge_density
Bd = bridge_density*Bs - 0.05  # Bd from Bs for a given bridge_density

Bfw_max = (1 - Bs/(bridge_density_max*Bs - 0.05))/Kdl # Bfw from Bs and bridge_density
Bfw_max = (1 - (Bd + 0.05)/(bridge_density_max*Bd))/Kdl # Bfw from Bd and bridge_density
```

#### Legacy external bridges

Legacy external bridges, use an adjusted solid infill line. This means it his
heavily dependent on the solid infill line width. In the following
`Ls`, `Lw`, and `La` refer to the solid infill line characteristics.

The translations to our model from OrcaSlicer settings are;

```python
Blh = max(bridge_flow, 1)*Lh    # flow-adjusted solid infill layer height.
Bls = min(bridge_flow, 1)*Ls    # flow-adjusted solid infill line spacing.
Blw = Bls + Blh*Kda             # flow-adjusted solid infill line width.
Ba = La * bridge_flow = Blh * Bls
Bd = sqrt(Ba/Ka)
Bs = Bls / bridge_density
```

The translation to OrcaSlicer settings from our model is;

```
bridge_flow = Ba / La
bridge_density = Bls / Bs
```

**Note: for `bridge_flow<1` this gives `Ar=bridge_density`, and when
`bridge_density=100%` this gives `Bs=Ls*bridge_flow`. 

The affect of the OrcaSlicer `bridge_density<=bridge_density_max` limit
depends on whether `bridge_flow` is greater or less than one. If it is greater
than one, then `Bs=Ls/bridge_density` which means the `bridge_density_max`
limit gives a minimum possible `Bs` independent of `Bd`. If it is less than
one, then `Bs=bridge_flow*Ls/bridge_density`, so the `bridge_density_max`
limit still gives a minimum `Bs` but it also decreases linearly with
`bridge_flow`. Since `Bd` decreases at the square-root while `Bs` decreases
linearly with `bridge_flow`, `Bfw_max` increases as you decrease
`bridge_flow`. This means to get decent overlap you need to either increase or
decrease `bridge_flow` enough.

When specifying the bridge thickness with `Bd`, there is a `Bd_min`
corresponding to `bridge_flow=1` where the calculation of `Bfw_max` flips
between the low and high flow formulas.

So when specifying the bridge thickness with `Bs` there is a `Bs_min` over
which `Bfw_max` is effectively `1/Kdl` because `bridge_flow` can be set
arbitrarily high to increase `Bd` without affecting `Bs`. When `Bs<Bs_min`
`bridge_flow` must be less than one to satisfy the `bridge_density_max` limit,
giving a `Bfw_max` limit calculated using the low-flow formula.

```python
Bs_min = Ls/bridge_density_max   # min Bs for bridge_flow >=1
Bd_min = sqrt(La/Ka)             # min Bd for bridge_flow >=1

# for bridge_flow >= 1
Bs = Ls/bridge_density
Bfw_max = (1- Bs_min/Bd)/Kdl    # from Bd for Bd>=Bd_min
Bfw_max = 1/Kdl                 # from Bs for Bs>=Bs_min

# for bridge_flow <1
Bs = bridge_flow*Ls/bridge_density
Bs = Bd^2*Ka/(Lh*bridge_density)                # Bs from Bd
Bd = sqrt(Bs*Lh*bridge_density/Ka)              # Bd from Bs
Bfw_max = (1-Bd*Ka/(Lh*bridge_density_max))/Kdl         # from Bd for Bd<Bd_min
Bfw_max = (1-sqrt(Bs*Ka/(Lh*bridge_density_max)))/Kdl   # from Bs for Bs<Bs_min
```

**NOTE:** the typical solid infill line has an area about half the nozzle area
but a line spacing about the nozzle diameter. For defaults `bridge_flow=1` and
`bridge_density=100%` this means the default bridge lines are thin and widely
spaced. Increasing `bridge_density` will reduce the spacing, but even at the
120% max they will usually not be touching. However, increasing `bridge_flow`
increases the area without increasing the spacing, allowing you to achieve
similar `Bd` bridge diameters as thick bridges but with higher `Bfw` overlap
factors. Also, reducing `bridge_flow` reduces the spacing at a rate faster
than it reduces the bridge diameter, allowing you to set small diameter
bridges with overlap. However, it does make bridge line characteristics vary
sigificantly with `Lh` layer height and `Ls` solid infill spacing for the same
bridge settings, requiring them to be set specifically for different layer
heights and solid infill settings.

### Recommendations

#### Settings for current OrcaSlicer

The ideal bridge line probably has a diameter `Bd` slightly less than the
nozzle diameter `Dn` and an overlap-factor of `Bfw=1`. It's possible that
smaller diameters also have good quality but will require more lines taking
longer to print. It's possible that `Bfw=1` is too much overlap so `Bfw=0.5`
might be better. Too much overlap would probably cause artifacts like
alternating high/low bridge lines for the whole length of the bridge, or some
kind of smear/tear/mess from damaging the adjacent lines. Too little overlap
would cause insufficient support and poor line bonding, resulting in drooping
loose lines.

The ideal anchor line should have an extrusion ratio of `Ar=1`. Lower
extrusion ratios might also be OK and help reduce pressure artifacts at the
anchor-to-bridge transitions, but risk poor line and layer adhesion. Line and
layer adhesion might not matter for the small anchor lines inside the print.

Since OrcaSlicer doesn't adjust the flow for anchor lines, an ideal large
diameter bridge line results in over-extruding the anchor lines. This means we
must choose a compromise between how thick the bridge lines are and how
over-extruded the anchor lines are.

OrcaSlice also doesn't adjust the height of the bridge to account for the
bottom perimeter based on the thickness of the brige lines, which means
dimensional accuracy favours using bridge lines with a thickness matching the
layer height, which is normally also less than the ideal bridge line diameter.

External bridge lines define the external perimeter and are visible, so it's
better to favour dimensional acuracy and quality over speed. They also have
longer anchor lines that increase the risk of artifacts from anchor line
over-extrusion. This suggests favouring smaller diameter bridge lines with
decent overlap for support.

Internal bridge lines are less important so speed matters more. They also have
shorter anchor lines that have less risk of artifacts from anchor
over-extrusion. This suggests using a larger diameter with less overlap for
faster printing.

In theory OrcaSlicer's thick bridges setting uses a better bridge model than
the legacy bridges, but unfortunately it includes a small hard-coded gap
between lines that makes it hard to get decent overlap for thinner bridges. It
is possible to get OK overlap for large diameters. The legacy bridges can get
better overlap for small diameters. It is also possible to get OK overlap for
larger diameters, but not as good as you can get with thick briges.

The 120% limit on `bridge_density` limits how much overlap you can get for
larger diameters using either bridge type, but limits it a little less for
thick bridges. It doesn't limit the overlap for small diameters with legacy
bridges.

This suggests a few options with various different compromises.

##### Print and Line settings

First select your layer height and preferred line-spacing.

Larger layer heights allow for larger bridge diameters with less anchor
over-extrusion. Something near 1/2 to 3/4 your nozzle diameter is good. I use
0.2mm for a 0.4mm nozzle.

Only the "internal solid infill" line spacing matters for bridges, but the
width of any line type can be set this way to get your desired line-spacing.
Use a spacing equal to or slightly larger than your nozzle diameter. I use
0.4mm or 0.5mm for a 0.4mm nozzle.

```python
Dn = 0.4                # your printers nozzle diameter
Lh = 0.2                # the layer height setting to use in OrcaSlicer.
Ls = 0.5                # your print's line spacing
Lw = Ls + Lh*(1 - pi/4) # the line width setting to use in OrcaSlicer.

# For different Lh layer heights this is;
Lw = Ls + 0.0215        # For Lh=0.1mm
Lw = Ls + 0.0429        # For Lh=0.2mm
Lw = Ls + 0.0644        # For Lh=0.3mm
```

##### Bridge best quality settings.

This assumes there are no problems with very thin bridge lines with a diameter
just larger than the layer height. Use this to optimize dimensional accuracy
of the bottom perimeter with primeter error of only `Be = Lh*Kdl/Kl =
Lh*0.128`, have no over-extrusion of anchor lines with Ar=1, and have optimal
bridge overlap with `Afw=1`. This is what you probably want for exernal briges.

```python
# gett layer height setting used in OrcaSlicer
Lh = layer_height

# get line width setting for internal solid infill from OrcaSlicer.
Lw = solid_infill_line_width

# get solid infill line spacing from line width
Ls = Lw - Lh*(1-pi/4)

# set thick bridges off in OrcaSlicer
thick_bridges = off

# set bridge flow setting to use in OrcaSlicer
bridge_flow = Lh/Ls

# set bridge density setting to use in OrcaSlicer
bridge_density = 100%
```

If the the overlap for Afw=1 is too high you can keep the same diameter but
increase the spacing by reducing `bridge_density` to something below 100% but
above 88.6% (for `Bfw=0`) like 94% (for `Bfw=0.5`). This keeps the dimensional
accuracy but reduces the anchor extrusion ratio to `Ar=bridge_density` which
might also help reduce some kinds of pressure artifacts near the
anchor-to-bridge transition.

Alternatively you can increase `bridge_flow` to increase the diameter and
spacing to preserve `Ar=1` anchor extrusion while reducing `Afw` overlap and
slightly reducing the dimensional accuracy by increasing the `Be` bottom
perimeter error. Setting bridge_flow 13% higher will set `Afw=0.5`.

##### Bridge best fast settings.

This is a compromise middle diameter for bridges to give Ok quality and speed.
It gets Ok overlap of `Afw=0.40` and over-extrudes anchors by 20% (Ar=1.2),
for a bottom perimeter error of `Be=0.46*Lh`.

```python
# set thick bridges off in OrcaSlicer
thick_bridges = off

# set bridge flow setting to use in OrcaSlicer
bridge_flow = 1.67*Lh/Ls

# set bridge density setting to use in OrcaSlicer
bridge_density = 120%
```

##### Bridge rough fast settings.

This gives fast thick `Bd=0.9*Dn` bridges with minimal `Bfw=0.25` overlap,
giving `Be=0.9*Dn-Lh` bottom perimeter error and `Ar=0.728*Dn/Lh` anchor
extrusion. Use this for internal bridges, of maybe for external bridges if you
have some mechanism like a post-processing script to dial down the anchor
extrusion.

```python
# set thick bridges on in OrcaSlicer
thick_bridges = on

# set bridge flow setting to use in OrcaSlicer
bridge_flow = 0.9

# set bridge density setting to use in OrcaSlicer
bridge_density = 120%
```

#### Sugested OrcaSlicer changes

Some possible changes to improve OrcaSlicer changes in increasing order of
complexity are;

1. Remove the 0.05mm spacing from thick bridges. This makes bridge overlap
vary with `bridge_flow` requiring very large `bridge_density` settings for
small diameter bridges to get decent overlap. With that offset removed
`bridge_density` will simply set the bridge spacing to `Bs=Bd/bridge_density`,
which means it will directly set overlap to `Bfw=(1-1/bridge_density)/Kdl` so
it will be independent of `bridge_flow`. A setting of `bridge_density=112.8%`
would be equivalent to the optimal `Afw=1.0` overlap factor. This would make
it possible to use thick bridges with very low `bridge_flow` settings to give
small diameters and still get optimal overlap.

2. Increase the `bridge_density` limit to about 160%. When not using thick
bridges with `bridge_flow=1` to get bridge diameters near the nozzle diameter,
`bridge_density` settings up to 158% are required to get `Bfw=1` optimal
overlap factors. Note with fix 1 above thick bridges will probably always be
the best option so this is not required, but would still be nice-to-have and
simple to implement.

3. Adjust the flow of the anchor part of bridge-lines to prevent
over-extrusion. The anchor lines should have an extrusion cross-section area
and corresponding flow rate of `Aa=Lh*Bs`, which is normaly much less than the
`Ba=Bd^2*Ka` area for bridge lines. This would allow bridge line settings to
be set to the optimum values for bridging without having to worry about
over-extrusion of the anchor lines.

4. Adjust which layer bridges are on to take into account the bridge
thickness. Optimal bridge diameters can be around 2x the layer height, which
means they hang a full layer lower than the bottom of the layer they are
printed on. Ideally the slicer should take into account bridge lines have a
different layer-height and include them on the layer where their
lower-perimeter most closely matches the model's perimeter. This would allow
bridge settings for optimum bridge quality without worrying about trying to
set bridge thickness close to the layer hight for dimensional accuracy.

5. Add Ui settings to allow setting `Bd` `bridge_diameter` or `Bs`
`bridge_spacing` and `Bfw` overlap-factor directly. It would also be nice to
be able to optionally set normal line thickness using `Ls` `line_spacing`
instead of `Lw` `line_width`.

## Other Factors

### Nozzle Temperature

Using higher temperatures makes the fillament more fluid and prone to sag. It
also increases the re-heating of adjacent bridge lines, risking more sag and
reducing the support they can provide for the current line.

Lower temperatures make the filament less fluid and probably more likely to
"snap" than "string", particularly for bridge lines with diameters much
smaller than the nozzle.

### Fan Speeds

Higher fan speeds increase cooling rates and should reduce sag. However, very
strong fans could blow the the bridges away before they "set", potentially
increasing sag, introducing sideways displacement, and maybe even "snapping".

### Print Speeds

Lower print speeds have a heap of advantages;

* Gives the filament more time to flow, reducing snapping.
* Supports the line with string-tension to the nozzle for longer, reducing
  sagging.
* Keeps the line under the extruder fan longer, cooling it faster to
  reduce sag.
* Reduces pressure and flow artifacts large velocity changes when accelerating.

However, there are also some possible risks for extremely slow print speeds.

* Gives the filament more time to flow, risking sagging.
* Doesn't apply as much string-tension, risking sagging.
* Keeps the hot nozzle on the extruded line melting it for longer, risking
  sagging.
* Keeps the hot nozzle close to the adjacent lines for longer, risking
  re-melting them and adding sagging.
* Speeds significantly lower than used for printing other lines risk pressure
  and flow artifacts from large speed changes after/before printing slowly.

### Pressure Advance

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
bridge and anchor areas are the same when the bridge spacing is the same as
the layer height, or `Bs=Lh`.

This is normally significantly less than the nozzle diameter, which might be
OK for very stringy materials, but could cause "snapping" for some materials.

### Line Bonding Strength.

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
