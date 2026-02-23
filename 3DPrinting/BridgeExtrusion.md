# Line and Bridge Extrusion Maths

This defines a model for FDM printing's extrusion lines for both normal lines
and bridge lines. They are based on the model for normal lines used by
OrcaSlicer that comes from Slic3r, but extends it in a few ways that can be
used to improve OrcaSlicer's bridging implementation, shows how we can map it
as close as possible to the current OrcaSlicer settings, and suggests changes
to OrcaSlicer to make it easier to set ideal normal and bridge line settings.

Speed, temperature, and cooling do have an effect on normal and bridge line
quality, but this document focuses mainly on the line geometry. So it focuses
on modeling and optimizing the extrusion volume and spacing of the lines
indepenently of things like print-speed, nozzle-temperature, and fan-speeds.

## References

This is based on ideas and discussions at;

* https://manual.slic3r.org/advanced/flow-math
* https://github.com/OrcaSlicer/OrcaSlicer/pull/11105
* https://github.com/OrcaSlicer/OrcaSlicer/pull/11255
* https://github.com/OrcaSlicer/OrcaSlicer/issues/11209
* https://github.com/OrcaSlicer/OrcaSlicer/issues/11412
* https://github.com/OrcaSlicer/OrcaSlicer/issues/11954

Much of which was spawned by the following videos;

* https://www.youtube.com/watch?v=xQBLv3cPUbo
* https://www.youtube.com/watch?v=eaasEkFULKE

There is also the following supporting documents

* [spreadsheet for translating to/from OrcaSlicer](https://docs.google.com/spreadsheets/d/1fX_NlD-0Lr_l3PO8elObgf09ENHlFkTImmb3xV1XEus/edit?usp=sharing)
* [Onshape CAD model](https://cad.onshape.com/documents/d6c789081e4718027433abc0/w/fab8c0f8c830995797647553/e/4c161cfdaa646b1caebb929e?configuration=Advanced_Settings%3Dtrue%3BAfw%3D1.0%3BAtl%3D4.0E-4%2Bmeter%3BBd%3D3.8E-4%2Bmeter%3BBdr%3D0.95%3BBfh%3D1.0%3BBfp%3D0.0%3BBfw%3D1.0%3BBh%3D3.37E-4%2Bmeter%3BBs%3D3.37E-4%2Bmeter%3BDn%3D4.0E-4%2Bmeter%3BLfp%3D0.0%3BLfw%3D1.0%3BLh%3D2.0E-4%2Bmeter%3BList_10u8p1SlZapMRd%3DDefaults%3BList_vwkf0slImDvp20%3DSet_Line_Width%3BLs%3D4.0E-4%2Bmeter%3BLw%3D5.0E-4%2Bmeter&renderMode=0&uiState=6960a1006d2f50394174a759)

## Intro

The geometry of any kind of line can be specified in many different ways using
different parameters based on different kinds of models. It is possible to
translate parameters for any way of specifying the line geometry into
parameters for any other way. However these all ultimately need to be
translated into the fundimental parameters defining the line that are used by
the slicer and printer.

The geometry of any kind of line depends fundimentally on three parameters;
the vertical spacing, the horizontal spacing, and the flow or extruded volume
per line length. The horizontal and vertical spacing together define another
important derived attribute; the cross-section area of the line.

Normal lines are lines printed in layers stacked on top of each other. They
are printed adjacent to other lines or perimeters and "squished" onto the
surface of the layer below or print bed. Their fundimental parameters and
derived atrribute are;

* `Lh` - the `layer_height` or vertical distance in mm between the top and
  bottom of the line, as defined by the top and bottom of the layer.

* `Ls` - the `line_spacing` the horizontal distance in mm between adjacent
  lines with the same `line_spacing`.

* `Lf` - the `line_flow` extruded volume per line length in mm^3/mm, or
  extruded cross-section area in mm^2.

* `La = Ls * Lh` - the `line_area` line available space cross-section area in
  mm^2, or line volume per line length in mm^3/mm.

Bridge lines span a gap in the layer below with other layers stacked above
them. The are printed adjacent to other lines with nothing underneath them, so
they are not printed "squished" against a print surface but "strung" through
space with only their ends supported by the layer below. This means their top
and bottom surface is not clearly defined by the top surface of the layer
below and the height of the print head. They can dangle down with their
defined top surface below the height of the print head. This introduces
another important parameter for printing bridge lines, giving us the following
fundimental parameters and derived attribute;

* `Bh` - the `bridge height` vertical distance in mm between the defined top and
  bottom of the brige line. Note this is not necesarily defined by the top and
  bottom of the current layer.

* `Bs` - the `bridge_spacing` horizontal distance in mm between adjacent lines
  with the same `bridge_spacing`.

* `Bf` - the `bridge_flow` extruded volume per line length in mm^3/mm, or
  extruded cross-section area in mm^2.

* `Ba = Bs * Bh` - the `bridge_area` line available space cross-section area
  in mm^2, or line volume per line length in mm^3/mm.

* `Bz` - the bridge `print_offset` distance in mm between the print-head and
  the defined top surface of the bridge line.

When printing bridges, there are three different kinds of lines involved;

* Normal perimeter and solid fill lines that provide support under the
  bridge ends and are printed on top of the bridge

* Bridge lines that span the bridge gap.

* Anchor lines that sit on top of supporting lines of the layer below at the
  ends of the bridge lines. They are extensions of the bridge lines so they
  share the same `bridge_spacing`. However they are printed on top of supporting
  lines, so they are like normal lines. This means they have different line
  characteristics and ideally require different settings. Importantly their
  available space cross-section area can be very different, requiring
  different flows to fill that area without over/under extruding.

The transition point between anchor and bridge lines ideally requires changing
the line parameters. Sudden hops, line characteristic changes, and parameter
changes can have undesirable artefacts from print speed and pressure changes,
so it helps to minimize these changes and/or apply them slowly.

These fundimental parameters are not what OrcaSlicer uses for specifying line
and bridge parameters. Instead it uses an underlying extrusion model for the
lines which are defined in terms of parameters for that model. These
parameters are then translated internally into the fundimental parameters for
printing. The models and parameters it uses are IMHO not very intuitive, so
the models could be tweaked and better parameters exposed.

Any extrusion model ultimately needs to be able to express the fundimental
line parameters using model parameters that are easy to understand and
correspond to important physical characteristics of the resulting line.
Ideally they should have minimal coupling between how they map to the
important characteristics, so changing one parameter to change one
characteristic should not require changing other parameters to preserve other
characteristics.

The important characteristics of any line are;

* the vertical `line height` used by the line excluding any kind of modeled
  overlap with lines above or below. This is the vertical space used by the
  line, so `N` lines above each other should have a height of `N*line_height`.

* the horizontal `line_spacing` used by the line excluding any kind of modeled
  overlap with adjacent lines. This is the horizontal space used by the line,
  so N adjacent lines beside each other should be `N*line_spacing` wide. Note
  OrcaSlicer doesn't use `line_spacing` and instead uses the model's internal
  `line_width` which is wider than `line_spacing`. IMHO this is not ideal
  because `line_width` doesn't reflect the line geometry of adjacent lines for
  overall dimensional accuracy, but instead includes an extra "bulge" of
  overlap with the adjacent lines that varies with `layer_height`. This bulge
  really only matters at the external perimeters for outer-wall lines.

* the degree of physical `line_contact` for bonding or space between adjacent
  lines. This is how tightly the extruded lines are packed together or spaced
  apart. Ideally this should directly correspond with how strong the bond and
  support with adjacent lines is, what the surface coverage is for sparse
  grills or internal supports, and/or the what the extrusion ratio is for the
  degree of overfill or underfill of the printed volume. These could be
  specified independently for horizontal and vertical contact but are often
  coupled.

* the `perimeter_edge` for external lines that defines where the line's
  external perimeter is for aligning it to the printed model's perimiter. The
  external surfaces of these lines not actually flat, so the perimeter can be
  defined in different ways like "outer-most edge", "inner-most edge", or
  "average edge" of the extrusion model's cross-section.

The `line_height` and `line_spacing` are pretty straightforward and directly
correspond with the fundimental line parameters.

There are many different ways that `line_contact` could be specified. It
largely depends on the model used, and typically on what special `line_width`
it defines for the width of the extruded line cross-section. Some
possibilities are;

* `extrusion_ratio = line_flow / line_area` - the ratio of extruded
  cross-section area to line available space cross-section area or fill-ratio
  of the printed volume. This is 1 when the extruded flow matches the
  available area, <1 when there are voids or gaps in the line, and >1 when the
  flow exceeds the available area. This directly controls over vs under
  extrusion that produces void or bulge artifacts. It is less useful for
  controlling line bonding or surface coverage.

* `overlap_ratio = (line_width - line_spacing)/line_width` - the ratio of
  overlap distance to line_width. This is 0 when lines just touch, >0 if they
  overlap, <0 if there is a gap, and tends to 1 when line-spacing tends to 0.
  This is useful for controlling the overlap distance relative to `line_width`
  for dimensional accuracy, but doesn't give a clear indication of line
  bonding, surface coverage, or extrusion ratio.

* `line_density = line_width/line_spacing` - the ratio of line_width to
  line_spacing. This is 1 when the lines just touch, >1 when the lines
  overlap, <1 when there are gaps between the lines, tends towards infinity as
  line_spacing tends towards 0, and tends towards 0 as line spacing tends
  towards infinity. This is good for controlling surface coverage for grills
  or internal supports. It is OK-ish for controlling line bonding because it
  has a gradient that reflects increasing sensitivity as the overlap
  increases, but the range of values for reasonable overlap is too small and
  usually varies with `layer_height`. For this model it's between 100% to
  112.8% for bridges and 100% to 104.3% for `line_spacing=0.5mm`,
  `layer_height=0.1mm` lines. These problems also apply to using it for
  extrusion ratio.

* `overlap_factor = <fancy model formula>` - A special model specific factor
  for measuring the amount of overlap that reflects the degree of line
  bonding. For example the [flow
  maths](https://manual.slic3r.org/advanced/flow-math) used by OrcaSlicer,
  which is 0 for "just touching", 1 for "full contact", and <0 for gaps. This
  is best for controlling line bonding. For this model it is nicely
  independent of `layer_height`, and also controls extrusion ratio because
  `overlap_factor=1` also corresponds with `extrusion_ratio=1`. However it's
  not useful for surface coverage.

There are also many different ways to define the `perimeter_edge`, but an easy
option is to also use the same `overlap_factor` formula to give a
`perimeter_factor` for the amount the edge-bulge overlaps the perimeter. This
means `perimeter_factor=0` is equivalent to "outer-most edge" which is useful
for clearance fits, and `perimeter_factor=1` is equivalent to "average edge"
which is useful for interference fits or accurate enclosed volume.

## A Bridge Extrusion Model

In OrcaSlicer, normal line cross-sections are based on the [flow
maths](https://manual.slic3r.org/advanced/flow-math) that models the extruded
line cross-section as a rounded-ended rectangle, which is used to define the
`line_width` and `line_flow`, with `overlap_factor` defining the
`line_contact` between adjacent lines. For the outer-wall lines we add a
`perimeter_factor` which is calculated the same as `overlap_factor` but
specifies the overlap for the external perimeter. We use this normal line
extrusion model for any line that is printed "squished" onto a supporting
surface.

**NOTE:** OrcaSlicer doesn't have exposed settings for `overlap_factor` which is
hard-coded internally to 1, or for `perimeter_factor`, which is effectively
internally hard-coded to 0. Hard-coding `overlap_factor=1` means `line_spacing
= line_flow/layer_height`, which is equivalent to a simple rectangular model
with `line_width = line_spacing`. Hard-coding `perimeter_factor=0` means the
perimeter is defined as at the outer-most-edge of the rounded-edge "bulge".
For completeness we include these factors in the model, but note they cannot
be changed in OrcaSlicer from the defaults.

Bridge lines are not squished against a layer below so we model them as
circular cross-sections that overlap with adjacent bridge lines (or possibly
perimeter lines at the bridge edges), and also with the layer above. The
`bridge_diameter` defines both the horizontal width and vertical thickness.

In addition to rounded sides that overlap horizontally with a
`horizontal_overlap_factor`, bridge lines have rounded top and bottom surfaces
which overlap with the layer above and bulge towards (or past) the external
perimeter below. This makes them like normal perimeter lines, where we can use
a `vertical_overlap_factor` to define the top surface overlap with the layer
above, and a vertical `perimeter_factor` to define the bottom
`perimeter_offset`. This gives us a `print_offset` printing height slightly
above the layer top surface height to account for the vertical overlap with
the layer above, and a `perimeter_offset` distance below the print height
where the bottom external perimeter is. From this we also get a
`perimeter_error` distance that the perimeter protrudes below the bottom of
the current layer.

The `bridge_height` is the vertical equivalent of `bridge_spacing`, with the
same `vertical_overlap_factor` on the top and bottom. It is the theoretical
`layer_height` of the bridge accounting for overlaps if printing a bridge line
on top of something, like perhaps a second bridge layer.

The `anchor_spacing` must be the same as `bridge_spacing` since they are
continuations of the bridge lines, and `anchor_height` must be the same as
`layer_height` since they are wedged between layer lines above and below.

Although `overlap_factor` is the primary parameter for specifying
`line_contact`, there are other ways to define this that might be better in
some cases. We also define `line_density` as the width to spacing ratio, and
`extrusion_ratio` as the flow to area ratio. For bridges these are defined
both horizontaly and vertically.

### Settings and parameters

The parameters that define these lines use the following input settings with
the suggested optimal ranges and default values. Note either `line_width` or
`line_spacing` can be set, and the one that is set will be used to derive the
other.

* Printer and layer settings
   * `Fd` - printer `filament_diameter` (range: [1.5mm,3.0mm], default: 1.75mm)
   * `Nd` - printer `nozzle_diameter` (range: [0.1mm,1.0mm], default: 0.4mm)
   * `Lh` - current `layer_height` (range: [0.1mm,1.0mm], default: 0.2mm)
* Normal line settings
   * `Lw` - `line_width` (range: [`Nd/2+Lh`,`Nd*2+Lh`], default: derived from
     `Ls`, `Lcw`, and `Lh`)
   * `Ls` - `line_spacing` (range: [`Nd/2,Nd*2`], default: `Nd`)
   * `Lcw` - line `overlap_factor` (range: [0,1], default: 1)
   * `Lcp` - line `perimeter_factor` (range: [0,1], default: 0)
* Bridge line settings
   * `Bd` - `bridge_diameter` (range: [`Nd/4`,`Nd`], default: `0.9*Nd`)
   * `Bs` - `bridge_spacing` (range: [`Nd/4`,`Nd`], default: derived from `Bd`
     and `Bcw`)
   * `Bcw` - bridge `horizontal_overlap_factor` (range: [0,1], default: 1)
   * `Bch` - bridge `vertical_overlap_factor` (range: [0,1], default: 1)
   * `Bcp` - bridge `perimeter_factor` (range: [0,1], default: 0)
* Bridge anchor line settings;
   * `Acw` - anchor `overlap_factor` (range: [0,1], default: 1)

And from these settings the following attributes can be derived;

* Printer and layer settings
   * `Fa` - printer `filament_area` or `Fd^2 * pi/4`.
   * `Na` - printer `nozzle_area` or `Nd^2 * pi/4`.
* Normal line attributes
   * `Ll` - width of a rectangle with the same `Lh` and area `Lf`.
   * `Lp` - line `perimeter_offset` from line center to perimeter.
   * `La` - `line_area` line cross-section area.
   * `Lf` - `line_flow` extrusion cross-section area.
   * `Lg` - `line_density` ratio of width to spacing.
   * `Lr` - line `extrusion_ratio` of flow to area.
* Bridge line attributes
   * `Bh` - `bridge_height` between defined top and bottom surfaces.
   * `Bl` - width of a square with the same area as `Bf`.
   * `Bz` - bridge `print_offset` distance from print-head to defined top
     surface.
   * `Bp` - bridge `perimeter_offset` distance from print-head to defined
     bottom perimeter.
   * `Ba` - `bridge_area` bridge cross-section area.
   * `Bf` - `bridge_flow` extrusion cross-section area.
   * `Bgh` - bridge `vertical_density` ratio of diameter to height.
   * `Bgw` - bridge `horizontal_density` ratio of diameter to spacing.
   * `Brh` - bridge `vertical_extrusion_ratio` of flow to height-area.
   * `Brw` - bridge `horizontal_extrusion_ratio` of flow to spacing-area.
   * `Br` - bridge `extrusion_ratio` of flow to area.
   * `Be` - bridge `perimeter_error` distance the bridge's bottom perimeter
     protudes past the bottom of the layer.
* Anchor line attributes
   * `Al` - width of a rectangle with the same `Lh` and area `Af`.
   * `Aw` - anchor `line_width` from extrusion edge-to-edge.
   * `Aa` - anchor `line_area` line cross-section area.
   * `Af` - anchor `line_flow` extrusion cross-section area.
   * `Ag` - anchor `line_density` ratio of width to spacing.
   * `Ar` - anchor `extrusion_ratio` of flow to area.

**NOTE:** The `line_flow` can be calculated from OrcaSlicer's reported flow
rates in mm^3/sec by dividing it by the reported speed in mm/sec. The
`line_area` can be calculated by getting the distance between adjacent lines
and multiplying it by the layer height.

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
Ll = Lw - Kda*Lh               # Ll from Lw, and Lh
Ls = Lw - Kda*Lcw*Lh           # Ls from Lw, Lh, and Lcw
Lw = Ls + Kda*Lcw*Lh           # Lw from Ls, Lh, and Lcw

La = Lh*Ls                     # La from Lh and Ls
Lf = Lh*Ll                     # Lf from Lh and Ll
Lg = Lw/Ls                     # Lg from Lw and Ls
Lr = Lf/La = Ll/Ls             # Lr from Lf and La or Ll and Ls
Lcw = (Lw - Ls)/(Kda*Lh)       # Lcw from Lh, Ls, and Lw.
Lp = (Lw - Kda*Lcp*Lh)/2       # Lp from Lw, Lh, and Lcp.
```

OrcaSlicer currently hard-codes `Lcw=1.0` and `Lcp=0.0` which limits and
simplifys this to;

```python
Ll = Ls
Ls = Lw - Kda*Lh           # Ls from Lw, and Lh
Lw = Ls + Kda*Lh           # Lw from Ls and  Lh

La = Lh*Ls = Lf
Lf = Lh*Ll = La
Lg = Lw/Ls
Lr = 1
Lp = Lw/2
```

**NOTE:** OrcaSlicer has extra `line_density` settings for some types of lines
which will adjust the `Ls` spacing without adjusting the `Lw` width or `Lf`
flow area. When it does this the `Lf` flow area is always calculated from the
unmodified `Lw` width. This means these `line_density` settings are
effectively a very coarse way of adjusting the `Lcw` overlap factor. There are
also `extrusion_ratio` settings that adjust the `Lf` flow without changing any
other setting, which is also equivalent to ajusting the `Lcw` overlap factor.

## Bridge Calculations

For bridges, overlap is calculated assuming that a nominal `line_spacing` with
horizontal and vertical `overlap_factor=1` is equivalent to a square
cross-section with the same area as the circular extrusion. This gives the
following relationships between settings and attributes.

```python
Bl = Kl*Bd                       # Bl from Bd.
Bs = Bd*(1 - Kdl*Bcw)            # Bs from Bd and Bcw
Bh = Bd*(1 - Kdl*Bch)            # Bh from Bd and Bch
Bd = Bs/(1 - Kdl*Bcw)            # Bd from Bs and Bcw
Bd = Bh/(1 - Kdl*Bch)            # Bd from Bh and Bch
Bd = sqrt(Bf/Ka)                 # Bd from Bf
Bcw = (1 - Bs/Bd)/Kdl            # Bcw from Bs and Bd
Bch = (1 - Bh/Bd)/Kdl            # Bch from Bh and Bd

Ba = Bh*Bs
Bf = Ka*Bd^2 = Bl^2              # Bf from Bd or Bl.
Bgw = Bd/Bs
Bgh = Bd/Bh
Brw = Bf/Bs^2 = (Bl/Bs)^2
Brh = Bf/Bh^2 = (Bl/Bh)^2
Br = Bf/Ba = sqrt(Brw*Brh)
Bz = (Bd - Bh)/2                 # Bz from Bd and Bh
Bz = Kdl*Bch*Bd/2                # Bz from Bd and Bch
Bp = Bd*(1 - Kdl*Bcp/2)          # Bp from Bd and Bcp
Be = Bp - Bz - Lh
```

If we copy OrcaSlicer's lines and apply similar hard-coded values of `Bcw=Bch=1`
and `Bcp=0` this simplifies to;

```python
Bl = Bs = Bh = Kl*Bd
Bd = Bs/Kl
Ba = Bf = Bs^2 = Ka*Bd^2
Bgw = Bgh = 1/Kl = 1.128
Brw = Brh = Br = 1
Bz = (Bd - Bh)/2 = Kdl*Bd/2
Bp = Bd
Be = Bd - Bz - Lh = Bd*(1+Kl)/2 - Lh
```

The `Be` bottom `perimeter_error` above assumes that the print head is
hopped-up by `Bz` to account for the `Bch` vertical overlap factor of the
rounded top of the bridge lines with the layer above, and that the rounded top
is at the print height. If the bridge line is printed without any hopup, this
is equivalent to having no vertical overlap, or `Bch=0`. We can hard code for
`Bcw=1`, `Bch=Bcp=0` instead which instead simplifies to;

```python
Bl = Bs = Kl*Bd
Bh = Bd
Bd = Bs/Kl
Bf = Bs^2 = Ka*Bd^2
Ba = Bs*Bh = Kl*Bd^2
Bgw = 1/Kl = 1.128
Bgh = 1
Brw = 1
Brh = Ka = 0.785
Br = Kl = 0.886
Bz = 0
Bp = Bd
Be = Bd - Lh
```

**NOTE:** For bridge diameters near the typical nozzle diameter of `Bd=0.4mm`
and optimal `Bch=1` vertical overlap, the optimal hopup height is `Bz =
0.057*Bd = 0.023mm`. This is well below the minimum `dz=0.1mm` Z resolution of
most 3DPrinters, so in practice adding a hopup is not worth it unless you have
a printer with very fine Z resolution or a very large nozzle diameter and use
very large diameter bridge lines. Or perhaps if the the assumption that the
rounded top of the bridge line is at the print height is wrong, and it
actually hangs significantly lower.

## Anchor Calculations

The anchor lines are extensions of the bridge lines, but they are printed over
supporting perimeter and infill lines, so they are printed squished like
normal lines and have different line characteristics to the unsupported bridge
lines. Without accounting for this there will be nasty artifacts at the start
and end of the bridge lines. Anchor line spacing is the same as bridge-spacing
`Bs`, and `layer_height` is the same as normal lines `Lh`.

```python
Al = Aw - Kda*Lh              # Al from Aw and Lh.
Aw = Bs + Kda*Acw*Lh          # Aw from Bs, Lh, and Acw.
Aa = Lh*Bs                    # Aa from Lh and Bs.
Af = Lh*Al                    # Af from Lh and Al.
Ag = Aw/As
Ar = Af/Aa = Al/Bs            # Ar from Af and Aa or Al and Bs
Acw = (Aw - Bs)/(Kda*Lh)      # Acw from Lh, Bs, and Aw.
```

If we set Acw=1 this simplifies to;

```python
Al = Bs
Aw = Bs + Kda*Lh
Aa = Af = Lh*Bs
Ag = Aw/Bs
Ar = 1
```

**NOTE:** these need to be printed with print height at the normal print
height of the layer, which for the ideal `Bch=1` is below the `Bz` offset
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
internal bridges. It also has two different modes that use different bridge
extrusion models;

* "thin" bridges models bridges as modified internal solid infill lines and
  aims to keep the same `line_spacing` or `layer_height` when modifying it. It
  uses the rounded box cross-section and adjusts the height while keeping the
  same spacing when increasing flows, and adjusts the spacing while keeping
  the same height when reducing flows.

* "thick" bridges models bridges with a circular cross-section where the
  height and width are the same as the diameter. The diameter is adjusted to
  adjust the flow. The default spacing is a small constant `0.05mm` larger than
  the diameter.

The relevant settings that affect bridges that OrcaSlicer exposes are;

* `layer_height` - the layer height for the print in mm. (default: `Nd/2`?,
  range ???)

* Internal solid infill `line_width` - The bridge `line_width` in mm for
  "thin" mode. (default: `Nd`, range: ???)

* `thick_bridge` - a checkbox toggling "thin" and "thick" bridge modes.
  (default: `N` for external, `Y` for internal). Available for both internal
  and external bridges.

* `bridge_flow` - multiplier for the flow (default:`1.0`, range: `[0,2]`).

* `bridge_density` - a percentage scaling factor for reduction in
  line-spacing. Higher values means closer spacing (default:`100%`, range:
  `[10%,100%]` in v2.3.1, `[10%,120%]` in v2.3.2-dev).

The current v2.3.1 OrcaSlicer limits bridge_density to <=100%, and
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

For thin bridges `Flow::with_flow_ratio()` is used to modify the solid
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
thick_bridge)` returns different flows for thin vs thick bridge mode. For
thick bridges it returns `Flow::bridging_flow()` with a diameter scaled by
sqrt of the `bridge_flow` setting. For thin bridges it calls
`Flow::with_flow_ratio()` to return a non-bridging solid infill flow adjusted
by the `bridge_flow` setting to give the required area.

For both bridge models, the `bridge_density` setting is then used to adjust
the line-spacing without changing the flow to `flow.m_spacing/bridge_density`.

### Observations

Tables for thick and thin bridges of the observed bridge characteristics from
OrcaSlicer's preview views for the given settings using v2.3.2-dev are below.
Both tables are using `layer_height=0.2` and `line_width=0.543mm` for the
Internal solid infill. Note `Bd` is calculated from the bridge flow `Bf`
calculated from the preview's "Flow" in mm^2/sec divided by the preview's
"Speed" in mm/sec. The observed spacing `Ls` is calculated from the distance
between adjacent lines. The preview's "Layer Height" gives us `Blh`, and the
preview's `Line Width` gives use `Blw`. We calculate perimeter-error `Be` and
anchor extrusion-ratio `Ar` from the other observed values.

**WARNING:** The OrcaSlicer preview modes showing the "Layer Height" and "Line
Width" of lines are a bit of a lie. They don't really indicate how wide or
high the line is. They don't even reflect the available volume for the line
between the adjacent layers and lines. The "Layer Height" doesn't change the
print height above the layer below from the layer's default, and doesn't
change a bridge's layer to align the bottom perimeter with the perimeter of
the printed model. The "Line Width" doesn't indicate how close adjacent lines
are, but for outer-wall lines it is accounted for and reflects the distance to
the edge perimeter. So what these views actually show is just the height and
width of the theoretical model cross-section of the extruded line that was
used to calculate how much filament to extrude, and the spacing to use between
adjacent lines. For bridges, the theoretical model used is often not even
close to the actual extruded line cross-section at all.

The table of results for thick bridges is;

|Thick bridge settings results                                 |
|bridge_flow|bridge_density|Bd  |Blh |Blw |Bs  |Bcw  |Be  |Ar  |
|-----------|--------------|----|----|----|----|-----|----|----|
|1.0        |100%          |0.40|0.40|0.40|0.45|-1.10|0.20|1.40|
|1.0        |113%          |0.40|0.40|0.40|0.40| 0.04|0.20|1.58|
|1.0        |120%          |0.40|0.40|0.40|0.38| 0.55|0.20|1.68|
|1.1        |120%          |0.42|0.42|0.42|0.39| 0.59|0.22|1.77|
|0.4        |120%          |0.25|0.25|0.25|0.25| 0.02|0.05|1.00|

This shows that for thick bridges the defaults have large gaps between bridge
lines (`Bcw<-1`) significant over-extrusion of anchor lines (`Ar>1.25`), and a
whole layer-height worth of perimeter-error (`Be=Lh`). Increasing
`bridge_density` above 114% is required to make them overlap (`Bcw>0`), but
even at 120% the overlap is still not significant (0.5<Bcw<1.0), and it
increases anchor over-extrusion (`Ar>1.5`). Increasing `bridge_flow` increases
overlap a bit more but makes the bridge diameter larger than the nozzle
diameter (`Bd>Nd`), increases anchor over-extrusion even more (`Ar>1.75`), and
slighly increases perimeter-error. Reducing `bridge_flow` to 0.4 removes the
anchor over-extrusion (`Ar=1`) and gives very little perimeter-error
(`Be=Lh/4`), but makes the lines barely touch even with `bridge_density` at
120%.

The table of observed results for thin bridges are;

|Thin bridge settings results                                  |
|bridge_flow|bridge_density|Bd  |Bh  |Bw  |Bs  |Bcw  |Be  |Ar  |
|-----------|--------------|----|----|----|----|-----|----|----|
|1.0        |100%          |0.36|0.20|0.54|0.50|-3.53|0.16|1.00|
|1.0        |120%          |0.36|0.20|0.54|0.42|-1.47|0.16|1.20|
|1.4        |120%          |0.42|0.28|0.56|0.42| 0.12|0.22|1.68|
|1.7        |120%          |0.47|0.34|0.57|0.42| 0.92|0.27|2.04|
|0.7        |120%          |0.30|0.20|0.39|0.29| 0.20|0.10|1.20|
|0.6        |120%          |0.28|0.20|0.34|0.25| 0.84|0.08|1.20|
|0.4        |100%          |0.23|0.20|0.24|0.20| 1.00|0.03|1.00|

This shows that for thin bridges the defaults have even larger gaps between
bridge lines than there is for thick bridges (`Bcw<-3`), the perimeter-error
is also large but less than one layer height (`Lh/2<Be<Lh`), and there is no
anchor over-extrusion (`Ar=1`). Increasing `bridge_density` to 120% reduces
the gaps but they are still large (`Bcw<-1`), and adds some anchor
over-extrusion (`1<Ar<1.25`). Also increasing `bridge_flow` to 1.4 finally
closes the gaps, and at 1.7 there is a healthy amount of overlap (`Bcw>0.75`),
but the bridge diameter is larger than the nozzle diameter (`Bd>Nd`),
perimeter-error is more than the layer height (`Be>1.25*Lh`) and there is
significant anchor over-extrusion (`Ar>2`). Reducing `bridge_flow` to 0.7
instead of increasing it also closes the gaps but reduces the perimeter-error,
and at 0.6 there is healthy overlap (`Bcw>0.75`), small perimeter-error
(`Be<Lh/2`), and moderate anchor over-extrusion (`1<Ar<1.25`). With
`bridge_flow` reduced down to 0.4 the overlap is good enough that we can
reduce `bridge_density` back to 100% to give bridge-spacing equal to the layer
height for "perfect" overlap (`Bcw=1`) and anchor extrusion (`Ar=1`) with
nearly no perimeter-error (`Be<Lh/4`).

Also interesting is external bridge anchor lines sit on top of 2 lines, one
outer wall and one inner wall, but internal bridge anchor lines sit on top of
a single internal solid infill line. This means the anchor lines are shorter
for internal bridges than external bridges, which would reduce the impact of
anchor line over-extrusion for internal bridges.

**NOTE:** although the bridge lines usually have a different layer-height than
the normal layer height, Orca doesn't do anything about adjusting the
extrusion heights or bridge layer to compensate for that. It always prints
bridges with the extruder at the same height as the rest of the layer, and it
always prints the bridge on the same layer regardless of how "low it hangs"
below the top of the layer. This means in practice that the bottom surface of
default bridges is much lower than the external perimeter height of the
printed model.

A lot of bridge testing done by others that I've seen have been using "thin"
bridges with significant `bridge_flow` increases (order 1.4x) combined with
`bridge_density` increases (around 110% to 125% with modified OrcaSlicer to
bypass settings limits) to produce decent results. These tests often don't
mention the `layer_height` or internal solid infill `line_width` used, but as
the modelling shows the flow rates for thin bridges are highly dependent on
them, so the ideal `bridge_flow` and `bridge_density` settings will vary
significantly with `layer_height` and internal solid infill settings.

If testing is done using "thick" bridges, it uses a better model and the flow
rates are more independent of `layer_height`. The flow rates are significantly
higher, and I guess ideal `bridge_flow` settings would be in the range
`[0.5,1.0]`. It's unfortunate that by default it adds a small 0.05mm space
between bridge lines because it also means `bridge_overlap` settings need to
be higher to overcome that, and the 120% limit is too low to get ideal
overlap. It also means that `bridge_flow` and `bridge_density` are coupled and
both affect overlap `Bcw`, making it harder to test for optimal diameter and
overlap independently.

For both kinds of bridges `layer_height` will have a significant impact on the
bridge perimeter-error and over-extrusion of anchor lines. Using larger
`layer_height` settings will help reduce perimeter-error and anchor
over-extrusion for larger bridge diameters. Using smaller diameter bridges
than is optimal for the bridge lines to minimize perimeter-error and avoid
over-extrusion of anchor lines will probably be the best compromise for
overall bridge quality. This means the `bridge_density` limit excludes using
`thick_bridges` because it's impossible to get decent overlap at smaller
diameters.

### Settings Translation

From these observations and looking at the code we can derive the formulas for
translating between our model's settings and OrcaSlicer's. For both bridging
models, the following mappings always apply.

Because it doesn't do any print-height compensation to account for the
normally larger bridge layer-height, some bridge line height characteristics
are fixed.

```python
Bch = 0.0  # no vertical overlap for bridges.
Bz = 0.0   # no print-height offset for bridges.
```

With no bridge height compensation at all, we cannot adjust the bridge's
external perimeter height at all, which is normally lower than the bottom of
the current print layer, which is generally the target perimeter. However, for
a given `Lh`, `Bp`, and `Bcp` perimeter factor to define exactly how we define
the perimeter, we can calculate the `perimeter_error` between the bottom of
the layer and the bridge's actual bottom perimeter.

```python
Be = Bp - Lh  # distance from layer bottom to bridge bottom perimeter.
```

**NOTE:** for `Bcp=0` (ie, no perimeter overlap so use lowest point of the
round bridge extrusion) this simplifies to just `Be = Bd - Lh`. Given that for
thick bridges typically `Bd` is near the nozzle diameter and `Lh` is near 1/2
the nozzle diameter, this means the perimeters below default thick bridges are
nearly a whole layer lower than they should be. For thin bridges the default
settings give a smaller diameter, but it still hangs below the layer bottom,
and for decent bridges you may need to increase the flow until it has the same
problem.

OrcaSlicer calculates `Blh` `layer_height` and `Blw` `line_width` values for
its bridge lines that are different for thick and thin bridges. It also
calculates a `Bls` nominal line spacing before applying `bridge_density` to
calculate the `Bs` actual spacing. We include them in the translation
calculations for each model, but these have no impact on the resulting
extrusion, which only depends on `Bd` and `Bs`.

OrcaSlicer does not differentiate between bridge lines and anchor lines, so
anchors use the same flows as bridges. This means we cannot adjust anchors,
but we can calculate their effective settings.

```python
Aa = Bs*Lh
Af = Bf
Al = Bf/Lh
Aw = Al + Kda*Lh
As = Bs
Acw = (Aw - Bs)/(Kda*Lh)
Ar = Af/Aa = Al/Bs             # Ar from Af and Aa or Al and Bs
Bd = Ar*Lh*(1 - Kdl*Bcw)/Ka    # Bd from Ar, Lh, and Bcw
```

For `Ar=1` anchor lines that are not over or under extruded we need the
following;

```python
Af = Aa            # required for Ar=1
Bf = Bs * Lh       # from Af=Bf and Aa=Bs*Lh
Bl^2 = Bl * Lh     # from Bf=Bl^2 and Bs=Bl for Bcw=1
Bl = Lh
Kl*Bd = Lh         # From Bl=Kl*Bd
Bd = Lh/Kl
```

This means to avoid overextruding anchor lines we need bridge diameters
slighly larger than the layer height. Note that this also has the nice
property that there is no perimeter error `Be=0` for the optimal `Bch=1`
vertical overlap, and only `Be = Kdl*Lh/Kl = Lh*0.128` for no vertical overlap
`Bch=0`.

**NOTE:** the default thick bridge lines have an area equal the nozzle area,
but normal lines have an area closer to half the nozzle area. This means
default thick bridge anchor lines are pretty close to 2x over-extruded. The
default thin bridge lines have the same area as solid infill lines so they
don't have this problem`.

#### Thick bridges

With Thick bridges enabled for internal or external bridges, the bridge and
anchor characteristics from the settings `bridge_flow` and `bridge_density`
are;

```python
Lh = layer_height
Blh = Bd
Blw = Bd
Bls = Bd + 0.05
Bf = bridge_flow*Na
Bd = sqrt(bridge_flow)*Nd
Bs = Bls/bridge_density
Bcw = (1 - Bs/Bd)/Kdl
Ar = Bf/(Bs*Lh)
Be = Bp - Lh
```

We can use this to calculate required settings from our target model's `Bf`,
`Bd`, `Bs`, and `Bls` values with;

```python
layer_height = Lh
bridge_flow = Bf/Na = (Bd/Nd)^2
bridge_density = Bls/Bs
```

Changing `bridge_flow` changes `Bd` and `Bls`, but the +0.05mm in `Bls` means
`Bs` changes are not directly proportional to `Bd` changes, so it changes
diameter `Bd`, spacing `Bs` , overlap `Bcw`, perimeter-error `Be` and anchor
extrusion-ratio `Ar`, which is everything related to bridge quality.

Changing `bridge_density` changes `Bls` without changing `Bd`, so it changes
spacing `Bs`, overlap `Bcw`, and anchor extrusion-ratio `Ar`, but diameter `Bd` and
perimeter-error `Be` don't change.

Changing layer height `Lh` changes perimeter-error `Be` and anchor
extrusion-ratio `Ar` only, so thick bridges are much less coupled to
layer-height changes than thin bridges.

Since OrcaSlicer limits the `bridge_density` to 120%, there is a maximum
possible overlap factor `Bcw` for a given `Bd` or `Bs`.

```python
Bs = (Bd+0.05)/bridge_density  # Bs from Bd for a given bridge_density
Bd = bridge_density*Bs - 0.05  # Bd from Bs for a given bridge_density

Bcw_max = (1 - Bs/(bridge_density_max*Bs - 0.05))/Kdl # Bcw from Bs and bridge_density
Bcw_max = (1 - (Bd + 0.05)/(bridge_density_max*Bd))/Kdl # Bcw from Bd and bridge_density
```

The +0.05mm extra space combines with the 120% `bridge_density` setting limit
to limit the possible overlap `Bcw` more as the diameter `Bd` gets smaller.
When the bridge diameter drops below 0.25mm, which is still larger than the
typical layer height, it is nolonger possible to make them touch.

#### Thin external bridges

Thin external bridges use an adjusted solid infill line. This means it his
heavily dependent on the solid infill line width and layer height. In the
following `Ls`, `Lw`, and `Lf` refer to the solid infill line characteristics.

The translations to our model from OrcaSlicer settings are;

```python
Lh = layer_height               # layer height of the print.
Blh = max(bridge_flow, 1)*Lh    # flow-adjusted solid infill layer height.
Bls = min(bridge_flow, 1)*Ls    # flow-adjusted solid infill line spacing.
Blw = Bls + Kda*Blh             # flow-adjusted solid infill line width.
Bf = Lf*bridge_flow = Blh*Bls
Bd = sqrt(Bf/Ka)
Bs = Bls/bridge_density
Bcw = (1 - Bs/Bd)/Kdl
Be = Bd - Lh
Ar = Bf/(Bs*Lh)
```

The translation to OrcaSlicer settings from our model is;

```python
layer_height = Lh
bridge_flow = Bf/Lf
bridge_density = Bls/Bs
```

Changing `bridge_flow` for `bridge_flow<1` changes `Bf` and `Bls` at the same
rate, but it means `Bd` and `Bs` change at different rates, so it changes
diameter `Bd`, spacing `Bs`, overlap `Bcw` and perimeter-error `Be`, but
anchor extrusion-ratio `Ar` doesn't change.

Changing `bridge_flow` for `bridge_flow>1` changes `Bf` without changing
`Bls`, so it changes diameter `Bd`, overlap `Bcw`, perimeter-error `Be`, and
anchor extrusion-ratio `Ar`, but spacing `Bs` doesn't change.

Changing `bridge_density` for any `bridge_flow` setting changes `Bs` without
changing `Bf`, so spacing `Bs`, overlap `Bcw`, and anchor extrusion-ratio `Ar`
change, but diameter `Bd` and perimeter-error `Be` don't change.

Changing `layer_height` changes `Bf` proportional to `Lh` but not `Bls`, so it
changes diameter `Bd`, overlap `Bcw` and perimeter-error `Be`, but not spacing
`Bs` or anchor extrusion-ratio `Ar`.

The affect of the OrcaSlicer `bridge_density<=bridge_density_max` limit
depends on whether `bridge_flow` is greater or less than one. If it is greater
than one, then `Bs=Ls/bridge_density` which means the `bridge_density_max`
limit gives a minimum possible `Bs` independent of `Bd`. If it is less than
one, then `Bs=bridge_flow*Ls/bridge_density`, so the `bridge_density_max`
limit still gives a minimum `Bs` but it also decreases linearly with
`bridge_flow`. Since `Bd` decreases at the square-root while `Bs` decreases
linearly with `bridge_flow`, `Bcw_max` increases as you decrease
`bridge_flow`. This means to get decent overlap you need to either increase or
decrease `bridge_flow` enough.

When specifying the bridge thickness with `Bd`, there is a `Bd_min`
corresponding to `bridge_flow=1` where the calculation of `Bcw_max` flips
between the low and high flow formulas.

So when specifying the bridge thickness with `Bs` there is a `Bs_min` over
which `Bcw_max` is effectively `1/Kdl` because `bridge_flow` can be set
arbitrarily high to increase `Bd` without affecting `Bs`. When `Bs<Bs_min`
`bridge_flow` must be less than one to satisfy the `bridge_density_max` limit,
giving a `Bcw_max` limit calculated using the low-flow formula.

```python
Bs_min = Ls/bridge_density_max   # min Bs for bridge_flow >=1
Bd_min = sqrt(Lf/Ka)             # min Bd for bridge_flow >=1

# for bridge_flow >= 1
Bs = Ls/bridge_density
Bcw_max = (1- Bs_min/Bd)/Kdl    # from Bd for Bd>=Bd_min
Bcw_max = 1/Kdl                 # from Bs for Bs>=Bs_min

# for bridge_flow <1
Bs = bridge_flow*Ls/bridge_density
Bs = Ka*Bd^2/(Lh*bridge_density)                # Bs from Bd
Bd = sqrt(Bs*Lh*bridge_density/Ka)              # Bd from Bs
Bcw_max = (1-Ka*Bd/(Lh*bridge_density_max))/Kdl         # from Bd for Bd<Bd_min
Bcw_max = (1-sqrt(Ka*Bs/(Lh*bridge_density_max)))/Kdl   # from Bs for Bs<Bs_min
```

**NOTE:** the typical solid infill line has an area about half the nozzle area
but a line spacing about the nozzle diameter. For defaults `bridge_flow=1` and
`bridge_density=100%` this means the default bridge lines are thin and widely
spaced. Increasing `bridge_density` will reduce the spacing, but even at the
120% max they will usually not be touching. However, increasing `bridge_flow`
increases the area without increasing the spacing, allowing you to achieve
similar `Bd` bridge diameters as thick bridges but with higher `Bcw` overlap
factors. Also, reducing `bridge_flow` reduces the spacing at a rate faster
than it reduces the bridge diameter, allowing you to set small diameter
bridges with overlap. However, it does make bridge line characteristics vary
sigificantly with `Lh` layer height and `Ls` solid infill spacing for the same
bridge settings, requiring them to be set specifically for different layer
heights and solid infill settings.

### Recommendations

#### Settings for current OrcaSlicer

The ideal bridge line probably has a diameter `Bd` slightly less than the
nozzle diameter `Nd` and an overlap-factor of `Bcw=1`. It's possible that
smaller diameters also have good quality but will require more lines taking
longer to print. It's possible that `Bcw=1` is too much overlap so `Bcw=0.5`
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
over-extrusion. This suggests using a larger diameter with less overlap, and
maybe even spacing, for faster printing.

In theory OrcaSlicer's thick bridges setting uses a better bridge model than
the thin bridges, but unfortunately it includes a small hard-coded gap between
lines that makes it hard to get decent overlap for thinner bridges. It is
possible to get OK overlap for large diameters. The thin bridges can get
better overlap for small diameters. It is also possible to get OK overlap for
larger diameters, but not as good as you can get with thick briges.

The 120% limit on `bridge_density` limits how much overlap you can get for
larger diameters using either bridge type, but limits it a little less for
thick bridges. It doesn't limit the overlap for small diameters with thin
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
Nd = 0.4                # your printer's nozzle diameter.
Lh = 0.2                # your print's layer height to use.
Ls = 0.5                # your print's line spacing to use.
Lw = Ls + Lh*(1-pi/4)   # get the line width from line spacing.
Ls = Lw - Lh*(1-pi/4)   # or get line spacing from line width.

# For different Lh layer heights this is;
Lw = Ls + 0.0215        # For Lh=0.1mm
Lw = Ls + 0.0429        # For Lh=0.2mm
Lw = Ls + 0.0644        # For Lh=0.3mm

# Set the OrcaSlicer settings to these;
layer_height = Lh       # The print's layer height setting.
line_width = Lw         # The internal solid infill line width setting.
```

##### Bridge best quality settings.

This assumes there are no problems with very thin bridge lines with a diameter
just larger than the layer height. Use this to optimize dimensional accuracy
of the bottom perimeter with primeter error of only `Be = Kdl*Lh/Kl =
Lh*0.128`, have no over-extrusion of anchor lines with `Ar=1`, and have
optimal bridge overlap with `Acw=1`. This is what you probably want for
exernal briges.

```python
thick_bridges = off
bridge_flow = Lh/Ls
bridge_density = 100%
```

If the the overlap for `Acw=1` is too high you can keep the same diameter but
increase the spacing by reducing `bridge_density` to something below 100% but
above 88.6% (`Bcw=0`) like 94% (`Bcw=0.5`). This keeps the dimensional
accuracy but reduces the anchor extrusion ratio to `Ar=bridge_density` which
might also help reduce some kinds of pressure artifacts near the
anchor-to-bridge transition.

Alternatively you can increase `bridge_flow` to increase the diameter `Bd` and
spacing `Bs` while keeping `Ar=1` anchor extrusion, reducing overlap `Acw`,
and slightly increasing the perimeter-error `Be`. Setting bridge_flow 13%
higher will give `Acw=0.5`.

##### Bridge best fast settings.

This is a compromise middle diameter for bridges to give Ok quality and speed.
It gets Ok overlap of `Acw=0.40` and over-extrudes anchors by 20% (Ar=1.2),
for a bottom perimeter error of `Be=0.46*Lh` or less than half a layer.

```python
thick_bridges = off
bridge_flow = 1.67*Lh/Ls
bridge_density = 120%
```

##### Bridge fat fast settings.

This gives fast thick `Bd=0.9*Nd` bridges with OK overlap `Bcw=0.45`,
fairly large perimeter-error `Be=0.9*Nd-Lh` and a fair bit of anchor
over-extrusion `Ar=0.745*Nd/Lh`. Use this for external bridges if having the
diameter near the nozzle diameter turns out to be more important than anchor
over-extrusion for bridge quality, small overlaps give good results, and you
care more about bridge appearance than dimensional accuracy. It's probably
also good for internal bridges. For even better results it could be used
together with some mechanism like a post-processing script to dial down the
anchor over-extrusion.

```python
thick_bridges = on
bridge_flow = 0.81
bridge_density = 120%
```

##### Bridge very fast and sparse settings.

This gives very fast thick `Bd=0.9*Nd` bridges with small gaps `Bcw=-2.21` for
80% coverage (`Bgw=0.80`), fairly large perimeter-error `Be=0.9*Nd-Lh` and
very little anchor over-extrusion `Ar=0.565*Nd/Lh`. It has spacing `Bs=0.45`
so it covers areas with less lines and prints fast. Use this for internal
bridges that just need to provide support for layers above.

```python
thick_bridges = on
bridge_flow = 0.81
bridge_density = 91%
```

#### Sugested OrcaSlicer changes

Some possible changes to improve OrcaSlicer changes in order of increasing
complexity and decreasing priority are;

1. Remove the 0.05mm spacing from thick bridges. This makes bridge overlap
   vary with `bridge_flow` requiring very large `bridge_density` settings for
   small diameter bridges to get decent overlap. With that offset removed
   `bridge_density` will simply set the bridge spacing to
   `Bs=Bd/bridge_density`, which means it will directly set the
   `overlap_factor` to `Bcw=(1-1/bridge_density)/Kdl` making it independent of
   `bridge_flow`. A setting of `bridge_density=112.8%` would be equivalent to
   the optimal `Acw=1.0` overlap factor for any `bridge_flow` setting. This
   would make it possible to use thick bridges with any `bridge_flow` setting
   and get optimal overlap, even for very small diameters.

2. Increase the `bridge_density` limit to about 160%. When not using thick
   bridges with `bridge_flow=1` to get bridge diameters near the nozzle
   diameter, `bridge_density` settings up to 158% are required to get `Bcw=1`
   optimal overlap factors. Note with fix 1 above thick bridges will probably
   always be the best option so this is not required, but would still be
   nice-to-have and simple to implement, at least until fix 4 below is
   implemented.

3. Change the UI settings to allow directly setting the bridge width. Perhaps
   replace the `bridge_flow` setting with a `bridge_width` setting similar to
   the `line_width` settings that can be a percentage of nozzle width or a
   number distance in mm. For thick bridges this would be the
   `bridge_diameter` for the circular model, and for thin bridges it would be
   the `line_width` for the normal line model. This would remove the confusing
   and arbitrary coupling of thin bridges and the internal solid infill line
   settings.

4. Adjust the flow of the anchor part of bridge-lines to prevent
   over-extrusion. The anchor lines should have a flow of `Af=Lh*Bs`, which is
   much less than the typical flow of `Bf=pi/4*Bd^2` for bridge lines. Perhaps
   remove the `thick_bridges` option and always use the circular model for the
   bridge parts and the normal line model for anchor parts, changing the
   reasonable limit for `bridge_density` back to about 120%. Note that
   `bridge_density` should adjust the spacing without affecting the flow for
   the bridge parts so that it controls bridge overlap, but it should also
   adjust the flow for the anchor parts so they are not over-extruded when
   there is overlap or under-extruded when there is gaps. This would allow
   bridge line settings to be set to optimum values for bridging without
   having to worry about over-extrusion of the anchor lines.

5. Adjust which layer bridges are on to take into account the bridge
   thickness. Optimal bridge diameters can be around 2x the layer height,
   which means they hang a full layer lower than the bottom of the layer they
   are printed on. Ideally the slicer should take into account the
   `perimeter_offset` (`Bp`) using a `perimeter_factor` (`Bcp`), but assuming
   a hard-coded `perimeter_factor=0` like there currently is for lines would
   be fine until fix 9 below is implemented. Note it already sets the bridge's
   `layer_height` to the bridge diameter (`Bd`) for thick lines, and using
   this to define the lower perimeter would be equivalent to
   `perimeter_factor=0` (`Bcp=Bch=0`). It should put the bridge lines on the
   layer where their perimeter most closely matches the model's perimeter.
   This would allow bridge settings for optimum bridge quality without
   worrying about trying to set bridge thickness close to the `layer_height`
   for dimensional accuracy.

6. It would be nice to be able to optionally set line and bridge sizes using
   spacing instead of width, since this better reflects the actual geometry
   and makes the `line_spacing` and thus overall wall thickness independent of
   `layer_height`.

7. It would be nice to add another preview view mode that shows `line_spacing`
   similar to the one that shows `line_width`.

8. It would be nice to be able to set the `overlap_factor` for both lines and
   bridges, instead of it being effectively hard-coded to `overlap_factor=1`.
   This would give fine-grained control over the amount of overlap for
   controlling adjacent line bonding and bridge support. Currently the overlap
   can be controlled by adjusting the `flow_ratio` for lines and
   `bridge_density` for bridges, but it is hard to map these into values that
   correlate with line overlap, the range of reasonable values is tiny, and
   they vary with `line_width` and `layer_height`. The `flow_ratios` settings
   are better for compensating for systematic under/over extrusion or filament
   expansion/contraction after the desired `overlap_factor` is applied.

9. It would be nice to be able to set the `perimeter_factor` for both lines
   and bridges instead of it being effectively hard-coded to
   `perimeter_factor=0`. This would give fine-grained control over how the
   perimeter is defined in terms of the edge-bulge overlap, with
   `perimeter_factor=0` being outer-most-edge-of-bulge (important for
   clearance fits), and `perimeter_factor=1` being average-edge-of-bulge
   (important for interference fits or enclosed volume). For lines, this can
   be controlled by adjusting the "X-Y hole compensation" and "X-Y contour
   compensation" settings, but it is hard to map these into values that
   correlate with perimeter overlap, the range of reasonable values is tiny,
   and they vary with `line_width` and `layer_height`. These settings are
   better for compensation for systematic dimension offets before
   `perimeter_factor` is applied. For bridges, without fix 5 above the bridge
   bottom perimeters are not taken into account at all and can even be on the
   wrong layer, but you can control it by carefully setting the bridge
   diameter based on the `layer_height`, but this is not easy and forces the
   use of bridge diameters that might not be optimal for bridging.


## Other Factors

### Nozzle Temperature

Using higher temperatures makes the fillament more fluid and prone to sag. It
also increases the re-heating of adjacent bridge lines, risking more sag and
reducing the support they can provide for the current line.

Lower temperatures make the filament less fluid and probably more likely to
"snap" than "string", particularly for bridge lines with diameters much
smaller than the nozzle.

In general it's probably best to just set the nozzle temperature to whatever
is normally used for that fillament. For heated chambers, it may be worth it
to slightly lower the chamber temperature.

### Fan Speeds

Higher fan speeds increase cooling rates and should reduce sag. However, very
strong fans could blow the the bridges away before they "set", potentially
increasing sag, introducing sideways displacement, and maybe even "snapping".

In general using max fan speeds is probably best, unless the extruder fan is
particularly powerful.

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

Tests by others suggest a very slow 10mm/sec gives the best results, but its
possible these very low speeds were required to mitigate affects from
over-extruded anchor-lines and bad Pressure Advance, so getting those right
might mean higher speeds can be used without adverse results.

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

This problem is made even worse by OrcaSlicer using the same flow
for the anchor lines as the bridging lines, despite the completely different
cross section area and thus flow rate required for those lines. This usually
means that flows and line spacing that works perfectly for the
bridge lines is hugely over-extruded for the anchor lines. For the short
anchor lines at the beginning and end of the bridges this usually just
contributes more excess pressure that is relieved as more over-extrusion at
the start of the bridge lines. However OrcaSlicer also often has extra bridge
lines at the sides of the bridges that overlap supporting perimeters and are
thus actually anchor lines for their whole length. For those the
over-extrusion can really build up and make a mess, as seen in this video;

https://www.youtube.com/watch?v=eaasEkFULKE

The only way to minimize the problems caused by different bridge and anchor
flows for slicers like OrcaSlicer that don't differentiate between them is to
pick settings that produce bridge and anchor lines with almost the same cross
section area. Assuming anchor overlap-factor `Acw=1` is reqired to maximise
anchor area and the optimal bridge overlap-factor is `Bcw=1`, the bridge and
anchor areas are the same when the bridge spacing is the same as the layer
height, or `Bs=Lh`.

This is normally significantly less than the nozzle diameter, which might be
OK for very stringy materials, but could cause "snapping" for some materials.

Better PA handling would take into account the different advance needed when
there is no backpressure from the print surface based on the layer below. Even
better would be to take this into account when there is partial backpressure
for overhangs. Failing that, PA problems can nearly always be mitigated by
using slower speeds.

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

**NOTE:** the actual overlap distances are tiny, and the difference between
`Bcw=0` "just touching" and `Bcw=1` "optimal/full overlap" is just 5.69% of
`Bd` bridge diameter on each side, which is 0.046mm for a full 0.4mm nozzle
diameter wide bridge, or 0.026mm for an optimal
`layer_height=bridge_spacing=0.2mm` bridge. For lines it's even less, and for
a `layer_height=0.1mm,line_spacing=0.5mm` line it's a very tiny 0.011mm.
