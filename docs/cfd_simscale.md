# CFD verification in SimScale

Goal: simulate the two optimised designs with a **conjugate heat transfer (CHT)** analysis, where conduction in the metal and flow and convection in the air are solved together. Then compare with the analytical model. SimScale runs in the browser on cloud computers; the free Community plan is enough, and its projects are public, so they can be linked from this README.

Interface contact is modelled as perfect in CFD, so compare with the model **without** the thermal paste resistance.

## Case 2 — CPU heat sink with fan (do this first: simpler domain)

### Geometry (Onshape, free account)

| Part | Dimensions |
|---|---|
| Heat sink (aluminium 6063) | base 80 × 100 × 8 mm; 32 fins 0.8 mm thick, 40 mm high, 100 mm long; fin pitch 2.5548 mm (spacing 1.7548 mm), first and last fin flush with the base edges |
| Heat source (copper block) | 37.5 × 37.5 × 2 mm, centred under the base |
| Air duct (fluid region) | box 80 mm wide × 40 mm high (from the base top to the fin tips), from 50 mm upstream to 100 mm downstream of the fins, **minus** the fins (Boolean subtract, keeping the heat sink) |

The duct walls touch the fin tips and the two outer fins: the flow is fully ducted, as in the model.

### Set-up

| Item | Value |
|---|---|
| Analysis | Conjugate heat transfer, steady state, **laminar** (channel Re ≈ 870), gravity off |
| Materials | air; aluminium k = 209 W/(m·K); copper k = 391 W/(m·K) |
| Inlet | velocity **3.08 m/s** (fan operating point from the model), temperature 35 °C |
| Outlet | pressure 0 Pa (gauge) |
| Duct walls | no-slip, adiabatic |
| Solid outer faces | adiabatic (all heat leaves through the air) |
| Power source | **125 W** in the copper block |
| Mesh | at least 5 cells across the 1.75 mm channels; run 2–3 meshes (coarse/medium/fine) for mesh independence |
| Result controls | area average temperature: bottom face of the copper block, base top face between the fins, outlet; area average pressure: inlet and outlet |

### What the model predicts (to compare)

| Quantity | Model | CFD |
|---|---|---|
| Pressure drop inlet → outlet | 50.6 Pa | |
| Mean air outlet temperature | 46.3 °C | |
| Mean base temperature (finned side) | 52.3 °C | |
| Mean source temperature (no paste) | 58.4 °C | |
| Heat sink resistance (source → inlet air) | 0.187 K/W | |

## Case 1 — LED heat sink, natural convection

### Geometry

| Part | Dimensions |
|---|---|
| Heat sink (aluminium 6063, black anodised) | base 150 × 150 × 6 mm; 16 fins 1.5 mm thick, 40 mm high, 150 mm long (vertical); pitch 9.9 mm (spacing 8.4 mm) |
| LED (aluminium block) | 28 × 28 × 2 mm, centred behind the base |
| Air (external region) | box around the heat sink, about 3 × the heat sink size sideways and in front, 2 × below and 4 × above (the plume rises) |

### Set-up

| Item | Value |
|---|---|
| Analysis | Conjugate heat transfer, steady state, laminar, **gravity on** along the fins (−L direction), **radiation on** |
| Ambient | 35 °C, pressure boundaries open on all sides of the air box |
| Emissivity | 0.85 on heat sink surfaces |
| Back of the base and of the LED | adiabatic (mounted on the lamp body) |
| Power source | **40 W** in the LED block |
| Result controls | area average temperature of the LED back face and of the base; heat flow by convection and radiation |

### What the model predicts

| Quantity | Model | CFD |
|---|---|---|
| Mean base temperature | 64.8 °C | |
| Mean LED temperature (no paste) | 70.9 °C | |
| Heat removed by radiation | 9.7 W (24 %) | |
| Heat sink resistance | 0.90 K/W | |

### Extra study: short fins

Repeat Case 1 with fins 20 mm high, at the model optimum spacing and at 12 mm. This checks the validation point where the model and the experiments disagree on the optimum spacing for short fins (air entering from the open front of the channels).
