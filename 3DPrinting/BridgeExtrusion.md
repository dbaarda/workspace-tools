# Bridge Extrusion Maths

There are three different kinds of lines around bridging;

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
different line characteristics and print height should be done gradually over
a small (1mm?) distance on the anchor-side of the transition edge.

## Settings and parameters

Line cross-sections are based on the flow maths at;

https://manual.slic3r.org/advanced/flow-math

This models the extruded line cross-section as a rounded-ended rectangle,
which is used to define the `line_width` and `line_area`, but the lines are
packed together with `overlap_factor=1`, which gives a `line_spacing` where
`line_spacing = line_area / layer_height`. For the outer-wall lines the
external perimeter is at the outer-most-edge of the rounded-edge "bulge", so
there is no perimeter overlap, or `perimeter_factor=0`.

We model bridging lines similarly but with a circular cross section, with
overlap on the horizontal edges with adjacent bridging lines or perimeter-wall
lines, and with the layer above.

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

* Normal line attributes
   * Lp - normal `perimeter_offset` from line center.
   * La - normal `line_area`
* Bridge line attributes
   * Bz - bridge `print_offset` distance from bridge-layer top up to print height.
   * Bp - bridge `perimeter_offset` from external perimeter up to print height.
   * Ba - bridge `line_area`
* Anchor line attributes
   * Aw - anchor `line_width`
   * Aa - anchor `line_area`

Note that the `line_area` is also the extrusion rate in mm^3/mm.

OrcaSlicer doesn't (yet) have settings for `overlap_factor` which is
effectively hard-coded internally to 1, or for `perimeter_factor`, which is
effectively internally hard-coded to 0. For completeness we include these
factors in the model, but note they cannot (yet) be changed in OrcaSlicer from
the defaults.

Bridge lines are not squished against a layer below so they are extruded as
circular sections that overlap with adjacent bridge (or perimeter) lines. The
`line_diameter` defines both the `line_width` and `line_thickness`, and is
typically the same or slightly smaller (due to stretching) than the nozzle
diameter. Trying to use larger diameters typically results in sagging instead
(what about foaming fillaments?).

In addition to rounded sides that overlap horizontally, they have rounded top
and bottom surfaces which overlap with the layer above and bulge towards the
external perimeter below. This makes them like normal perimeter lines, where
the rounded bulge makes the actual external perimeter ambiguous. They also
need to be printed from slightly above the layer top surface height to account
for the top overlap. The `layer_height` is the vertical equivalent of
`line_spacing` with the same Bfh overlap on the top and bottom, so it's the
theoretical `layer_height` if printing a bridge line on top of something, like
perhaps a second bridge layer.

Anchor-spacing must be the same as bridge-spacing, and anchor-height must be
the same as `line_height`, so there are no separate values for those.

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

## Bridge Calculations

For bridges, overlap is calculated assuming that a nominal `line_spacing` with
horizontal and vertical `overlap_factor=1` is equivalent to a square
cross-section with the same area as the circular extrusion. This gives the
following relationships between settings and attributes.

```
Bs = Bd*(1 - Bfw*(1-sqrt(pi/4)))   ; Bs from Bd and Bfw
Bh = Bd*(1 - Bfh*(1-sqrt(pi/4)))  ; Bh from Bd and Bfh
Bd = Bs/(1 - Bfw*(1-sqrt(pi/4)))   ; Bd from Bs and Bfw
Bd = Bh/(1 - Bfh*(1-sqrt(pi/4)))  ; Bd from Bh and Bfh

Bfw = (1 - Bs/Bd)/(1-sqrt(pi/4))   ; Bfw from Bs and Bd.
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
setting combinations could also work, but it would require careful limits to
keep all settings within reasonable limits;

1. Derive `line_spacing` and `layer_height`. Set Bd, Bfw, Bfh. Derive Bs, Bh.
1. Specify `layer_height` and derive `line_spacing`. Set Bd, Bh, Bfw. Derive Bs, Bfh.
1. Specify `line_spacing` and derive `layer_height`. Set Bs, Bfw, Bfh. Derive Bd, Bh.
1. Specify `line_spacing` and `layer_height`. Set Bs, Bh, Bfw. Derive Bd, Bfh.

## Anchor Calculations

The anchor lines are extensions of the bridge lines, but they are printed over
supporting perimeter and infill lines, so they are printed squished like
normal lines and require different line characteristics to the unsupported
bridge lines.

Anchor line spacing is the same as bridge-spacing `Bs`, and `line_height` is
the same as normal lines `Lh`.

```
Aw = Bs + Afw*Lh*(1-pi/4)         ; Aw from Bs, Lh, and Afw
Aa = Lh*(Aw - Lh*(1 - pi/4))      ; Aa from Lh and Aw.
Aa = Lh*(Bs - (1-Afw)*Lh*(1-pi4)) ; Aa from Lh, Bs, and Afw.
```

Note these need to be printed with `print_height` at the top of the layer,
which is below the `print_height` for the bridge lines they are extensions of.
This means there should be a little hopup at the transition from anchor to
bridge line.


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

Flows in Flow.hpp have lengths m_width, m_height, m_spacing,
m_nozzle_diameter, and bool m_bridge. `Flow::mm3_per_mm()` returns the area of
a circle with diameter=width for bridge flows, otherwise does the
rounded-rectange area using width and height.

For Bridging flows, `Flow::bridging_flow()` takes a diameter and the
nozzle_diameter, and creates a flow with `width = height = diameter` and
`spacing = Flow::bridge_extrusion_spacing(diameter) = diameter+0.05mm`

In LayerRegion.cpp, `LayerRegion::bridging_flow(FlowRole role, bool
thick_bridge)` returns `Flow::bridging_flow()` with a diameter scaled by sqrt
of the bridge_flow ratio, but only if thick_bridge is true, otherwise it
returns a non-bridging flow the same as the solid infill adjusted by the
flow_ratio. When increasing the flow ratio, the adjustment bumps the height,
when decreasing the flow, it adjusts the width and spacing... strange.

A table of the bridge characteristics for given settings is below. Note for
thick_bridge=N, Bd is calculated from the observed flow Ba;

Lh  Ls  thick_bridge bridge_flow bridge_density  Bd     Bh     Bw     Bs      Ba
0.2 0.5  Y           1.0         100%            0.40mm 0.40mm 0.40mm 0.450mm 0.1255mm^2
0.2 0.5  Y           1.0         110%            0.40mm 0.40mm 0.40mm 0.409mm 0.1255mm^2
0.2 0.5  Y           1.0         120%            0.40mm 0.40mm 0.40mm 0.375mm 0.1255mm^2
0.2 0.5  Y           1.1         100%            0.42mm 0.42mm 0.42mm 0.470mm 0.1380mm^2
0.2 0.5  Y           1.2         100%            0.44mm 0.44mm 0.44mm 0.488mm 0.1510mm^2
0.2 0.5  N           1.0         100%            0.36mm 0.20mm 0.54mm 0.500mm 0.1000mm^2
0.2 0.5  N           1.0         110%            0.36mm 0.20mm 0.54mm 0.454mm 0.1000mm^2
0.2 0.5  N           1.0         120%            0.36mm 0.20mm 0.54mm 0.417mm 0.1000mm^2
0.2 0.5  N           1.1         100%            0.37mm 0.22mm 0.55mm 0.500mm 0.1100mm^2
0.2 0.5  N           1.2         100%            0.39mm 0.24mm 0.55mm 0.500mm 0.1200mm^2

### Thick external bridges

With Thick external bridges enabled, the bridge characteristics from the
settings bridge_flow and bridge_density are;

```
Bd = sqrt(bridge_flow)*Dn
Bh = Bd
Bw = Bd
Bs = (Bd + 0.05)/bridge_density
Ba = Bd^2 * pi/4
Bfw = (1 - ((Bd+0.5)/Bd)/bridge_density)/(1-sqrt(pi/4))
```

We can use this to calculate required settings from our target Bd and overlap
factor Bfw;

```
bridge_flow = (Bd/nozzle_diameter)^2
bridge_density = (Bd + 0.05)/Bs
               = (Bd + 0.05)/(Bd*(1 - Bfw*(1-sqrt(pi/4))))
               = ((Bd + 0.05)/Bd)/(1 - Bfw*(1-sqrt(pi/4)))
```

This gives the following table for Bs/Dn, bridge_flow, and bridge_density for
Bd/Dn and Bfw values;

Bd/Dn   Bfw   Bs/Dn    bridge_flow  bridge_density
1.0    1.0    0.88622  1.0000       1.2694
0.95   1.0    0.84191  0.9025       1.2769
0.9    1.0    0.79760  0.8100       1.2851
1.0    0.5493 0.93750  1.0000       1.2000

For Bd/Dn=1.0 we need bridge_density=127%, but the valid settings range is
limited to 120%. As Bd/Dn drops, bridge_flow also drops, but bridge_density
increases even further. So in practice the highest overlap we can get is
`Bfw=0.5493` with `Bd = nozzle_diameter`, `bridge_flow=1.0`, and
`bridge_density=120%`.

Note that although it sets the bridge height to the diameter, it doesn't
actually change the bridging print height even layer to take that extra height
into account, and prints the bridge on the same layer with the same nozzle
height that it would use if it assume the bridge was only the normal layer
height. This means that the external surface of bridges is actually lower than
it should be by `Bd - Lh`.
