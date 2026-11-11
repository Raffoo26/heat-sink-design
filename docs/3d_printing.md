# 3D printing the heat sinks

The two optimised designs are printed in plastic as **physical models**: to see and handle the result of the calculation. A plastic heat sink does not cool anything (PLA conducts about 1000 times less heat than aluminium; see the table in the README).

## Parts

| File | Size [mm] | Suggested colour | Notes |
|---|---|---|---|
| `cad/led_natural_heatsink.stl` | 150 × 150 × 46 | grey / black (anodised look) | needs a bed of at least 150 mm; can be scaled to 70 % |
| `cad/led_source.stl` | 28 × 28 × 3 | yellow / white | the COB LED |
| `cad/cpu_fan_heatsink.stl` | 80 × 100 × 48 | grey / silver | 0.8 mm fins: **never scale it down** |
| `cad/cpu_source.stl` | 37.5 × 37.5 × 3 | red | the processor heat spreader |

All files are in millimetres. To export your own design: `write_stl(hs, "my_heatsink.stl")` from `heatsink.export`.

## Slicer settings

| Setting | Value | Why |
|---|---|---|
| Orientation | base on the bed, fins pointing up | fins are vertical walls: no supports, no overhangs |
| Supports | none | |
| Nozzle / layer height | 0.4 mm / 0.2 mm | 0.8 mm fins = exactly two extrusion lines |
| Walls (perimeters) | 2–3 | fins come out solid anyway |
| Infill | 10–15 % | only the base has infill; do not use 0 %: the top of the base would be printed over nothing |
| Top / bottom layers | 4 / 4 | closed, flat faces |
| Material | PLA or PETG | it is a model: no temperature requirement |
| Brim | 3–5 mm (CPU heat sink) | 40 mm tall, thin fins on a small base: better adhesion |
| Speed | reduce outer-wall speed on the CPU heat sink | thin tall walls wobble at high speed |

Before printing, check in the slicer preview that every fin is a continuous wall (if a fin disappears, the slicer thinks it is too thin: enable "detect thin walls" or the equivalent).

## Assembly

Glue the source block in the centre of the flat face of the base (the side without fins), as in the drawings of the case studies: it marks where the heat enters.

## Check the printed dimensions

Measure with a caliper and compare with the design (FDM printers are typically accurate to about ±0.1–0.2 mm):

| Dimension | Design [mm] | Measured [mm] |
|---|---|---|
| LED heat sink: fin thickness | 1.50 | |
| LED heat sink: fin spacing | 8.40 | |
| LED heat sink: fin height above the base | 40.0 | |
| CPU heat sink: fin thickness | 0.80 | |
| CPU heat sink: fin spacing | 1.75 | |
| CPU heat sink: base thickness | 8.0 | |

For the real aluminium part, the 1.75 mm channels of the CPU design would need bonded or skived fins: standard extrusion cannot make fins that thin and that close (fin height / spacing ≈ 23).
