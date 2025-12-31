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
   * Dn - nozzle_diameter
   * Lh - `layer_height` of the normal perimeter lines.
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
