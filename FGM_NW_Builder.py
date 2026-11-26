#!/usr/bin/env python3
"""
FGM Nanostructure Builder  and  Multilayer & Multi-element Nanostructure Builder   (version 1.1)

Copyright (C) 2026 Md. Shabbir Rahman Siam

FGM Nanowire Builder
Version 1.0

Copyright (c) 2026 Md. Shabbir Rahman Siam

FGM Nanowire Builder is a Python application for generating
functionally graded multilayer nanowire structures for atomistic simulations.

The software provides tools for nanowire construction, material
composition control, grading profiles, and exporting atomic
structures for molecular dynamics simulations.

Developed using Python, Tkinter, and NumPy.

This software is released under the MIT License.
See the LICENSE file for details.

Author: Md. Shabbir Rahman Siam
HOW TO USE   (the full manual is inside the program: Help > Manual)
    Open the window ........ double-click this file, or run:   py FGM_NW_Builder.py
    Build without a window . py FGM_NW_Builder.py build --direction axial --bins 10 --out wire.nw
    Check against a file ... py FGM_NW_Builder.py validate --ref reference_file.txt --direction axial
    Multilayer wire ........ py FGM_NW_Builder.py layered --example > spec.json   then   layered --spec spec.json --out wire.nw
    Graded along + across .. py FGM_NW_Builder.py combined --example > spec.json   then   combined --spec spec.json --out wire.nw
    Crystal-layer stacks ... py FGM_NW_Builder.py crystals --example > spec.json   then   crystals --spec spec.json --out wire.nw
    Built-in checks ........ py FGM_NW_Builder.py selftest
    More options ........... py FGM_NW_Builder.py build --help

NEEDS: Python 3 with NumPy. If NumPy is missing the program offers to install it
(internet needed once).
"""
import argparse
import math
import os
import sys
import tempfile
import time
try:  # the window needs Tkinter; the command line and the web version do not
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError:
    tk = filedialog = messagebox = ttk = None
_TK_ROOT = tk.Tk if tk else object
_TTK_FRAME = ttk.Frame if ttk else object

VERSION = "v1.1"
MANUAL_TEXT = r"""# FGM Nanostructure Builder - User Manual

Version 1.1, by Md. Shabbir Rahman Siam

## 1. What this program does

This program builds atom-position files for nanowires in the LAMMPS data-file format. You describe the wire (shape, size, crystal, elements, how the composition changes) and the program writes the file, ready to be read by LAMMPS with `read_data`. It builds STRUCTURE files only. It does not write LAMMPS input scripts and it does not run simulations.

The program has two parts, shown as two tabs:

- FGM Nanostructure Builder (first tab): a functionally graded nanowire, nanowire cuboid or nanoparticle (sphere) made of two materials. The composition changes gradually from the core to the surface (radial grading) or along the wire (axial grading).
- Multilayer & Multi-element Nanostructure Builder (second tab): wires made of stacked layers, wires graded along AND across, and wires made of different crystals (for example hexagonal Ti with cubic Al), with 2 to 5 elements.

Everything is in one file, `FGM_NW_Builder.py`. Its source code is readable on purpose.

## 2. Starting the program

You need Python 3 and the NumPy package.

- Windows: install Python from python.org (tick 'Add python.exe to PATH'). Then open Command Prompt once and type `py -m pip install numpy`.
- Start the window: double-click `FGM_NW_Builder.py`, or open Command Prompt in its folder and type `py FGM_NW_Builder.py`. A black console window may stay behind the program. It is harmless.
- If you see 'Python was not found', use `py` instead of `python`.
- If NumPy is missing, the program offers to install it (internet needed once).

Keep `LICENSE.txt` in the same folder as the program.

## 3. The window

- Left side: the input form (two tabs).
- Top right: the buttons Build, Save data file... and Compare with reference file...
- Below them the Result list. It chooses what you see in the large area: the data file, tables, plots and pictures.
- Bottom: a status line with the last message.
- Help menu (top left): this manual (also key F1).

Typical workflow: fill in the form, press Build, look at the tables and pictures, press Save data file..., then read the saved file in LAMMPS.

## 4. First tab: FGM Nanostructure Builder

### 4.1 The idea

The wire is cut into slices called bins. In every bin a certain fraction of the atoms belongs to material 2, and the rest to material 1. The fraction in each bin comes from a grading function that depends on the position b along the way from 0 (start) to 1 (end).

- Radial grading: b = r / R, from the axis of the wire (0) to its surface (1). The bins are rings.
- Axial grading: b = z / L, from one end of the wire (0) to the other (1). The bins are slices.

Atoms are chosen at random inside every bin, but the NUMBER of atoms of each material in every bin is fixed by the function.

### 4.2 Fields

- Shape: Cylindrical (a wire along z), Cuboid (a square or rectangular wire along z) or Sphere (a nanoparticle).
- X-dimension, Y-dimension (nm): the cross-section sizes of a cuboid. Greyed out for other shapes.
- Radius R (nm): radius of the cylinder or sphere.
- Length L (z) (nm): length of the cylinder or cuboid. A sphere has no separate length.
- Grading Direction: Radial or Axial.
- Grading Function: P-FGM, E-FGM, S-FGM or Custom (see 4.3).
- Custom function f(x): your own formula in x, only for Custom.
- p: the exponent of the function (0 or larger). Not used by Custom.
- Atomic Mass 1 and 2 (amu): masses of the two materials. They are written to the data file, and used to turn mass fractions into atom counts.
- Lattice Type: FCC or BCC.
- Lattice Parameter 1 and 2 (nm): lattice constants of the two materials. The wire uses their average for the whole crystal.
- Number of bins: how many slices or rings.

### 4.3 Grading functions

The value x is the fraction of material 2 (by mass in the standard setting). b runs from 0 to 1.

```
P-FGM:   x = b^p
E-FGM:   x = 1 - exp(-b^p)
S-FGM:   x = 0.5 * (2b)^p            for b <= 0.5
         x = 1 - 0.5 * (2(1-b))^p    for b > 0.5
Custom:  x = f(b)  (your formula)
```

Examples: P-FGM with p = 1 is a straight line from 0 to 1. A larger p keeps material 1 longer and then rises quickly. S-FGM is an S-shaped curve: slow at both ends, fast in the middle.

Custom formulas may use: x, numbers, + - * / ^, comparisons (> < >= <= ==), and / or, the functions exp, log, sqrt, sin, cos, tan, tanh, sinh, cosh, atan, asin, acos, abs, min, max, floor, ceil, pow, step(x,a) and box(x,a,b), and the constants pi and e. Results below 0 or above 1 are cut to 0 and 1. step(x,a) is 0 below a and 1 from a on; box(x,a,b) is 1 from a up to b and 0 elsewhere. These two (and comparisons such as x>0.7) make sharp edges, which give solid regions instead of a mix (see 5.9 in the second tab, where the same formulas work).

### 4.4 Advanced options

- Mode: Standard (validated) or Corrected bins.
  - Standard follows the published algorithm with all its details, so that results agree with it. Two details matter: along the wire the last bin holds only the top atomic layer (the effective number of bins is one less than typed), and the ring centre is about 1 Angstrom off the true axis, so the outermost atoms can form an extra ring.
  - Corrected bins: exactly the number of bins you typed, every atom in exactly one bin, rings centred on the true axis, and the length rounded to whole unit cells so that a periodic wire has no stacking seam.
- Grading sampled at (corrected mode only): Bin centre (default) evaluates the function in the middle of each bin. 'First bin 0 / last bin 1' makes the first bin pure material 1 and the last pure material 2.
- Grading value is a (corrected mode only): Mass fraction (default) or Atom fraction.
- Vacuum x,y and Vacuum z (Angstrom): empty space added around the wire. Zero means a tight box. With a tight box use free (non-periodic) boundaries in x and y in LAMMPS.
- Random seed: a whole number makes the random choice repeatable. Empty means a different random choice every time.

### 4.5 Buttons and results

- Build: makes the structure.
- Save data file...: writes the LAMMPS data file with all atoms.
- Compare with reference file...: see 4.6.
- Result list: LAMMPS DATA FILE (first 5000 atoms shown, the saved file has all), Bin table (atoms of each material in each bin), Grading curves (wanted curve and built points), Comparison with reference file.

### 4.6 Compare with reference file

Use this to check the program against a data file made by another program or an earlier version. Type exactly the same values in the form, press Compare with reference file... and pick the file. The program rebuilds the wire with the form values and checks: the header (atom count, box, masses), every atom position, and the number of atoms of material 2 in every bin. Which atoms are chosen is random, so the individual types cannot match. For axial grading, atoms exactly on a bin edge belong to two bins, so the check there is a statistical test. A neighbour test checks that the two materials are mixed randomly. The report ends with ALL CHECKS PASSED or lists what differs. The comparison always uses the Standard mode and no vacuum.

## 5. Second tab: Multilayer & Multi-element Nanostructure Builder

### 5.1 Structure type

Choose at the top:

- Layers (stacked along the wire): several layers one above the other, each with its own composition.
- Combined: graded along + across: composition changes along the wire and from the core to the surface at the same time.
- Crystal layers: FCC(111) / HCP(0001) mix: every layer is its own crystal, for example hexagonal Ti and cubic Al.

All three share the wire settings and the element table.

### 5.2 Wire and settings

- Shape: Cylindrical or Cuboid. Radius (nm) for a cylinder, X- and Y-dimension (nm) for a cuboid.
- Lattice type: FCC or BCC (not used for crystal layers, where every layer has its own crystal).
- Composition is in: Atom % (default, usual for alloys) or Mass %.
- Grading sampled at: Bin centre or first/last bin at the ends.
- Lattice parameter: Mean of elements (average of the constants in the table) or One value (the same constant for all). This is the spacing in the plane across the wire.
- Spacing along the wire: Same as in-plane, or Follows composition (Vegard), where every layer (or slice) gets the spacing of its own composition. See 5.7.
- Vacuum x,y and Vacuum z (Angstrom), Random seed: as in 4.4.

### 5.3 Elements (2 to 5)

Choose the number of elements. For each element give the symbol (a list of common elements fills in the mass and the lattice constant, and you can type any other), the mass (amu) and the lattice parameter a (nm). IMPORTANT: element 1, element 2, element 3 ... are LAMMPS atom types 1, 2, 3 ... in the data file. The order of the table is also the order of the numbers in every composition you type. Use that same order in `pair_coeff` in LAMMPS.

For elements without a standard FCC or BCC constant (for example Cr as FCC, or Ti), type the value from your interatomic potential.

### 5.4 Layers

Each row is one layer, from the bottom to the top.

- Thickness (nm): asked thickness. It is rounded to whole atomic layers; the table shows asked and built thickness.
- Start % and End %: the composition at the start and at the end of the layer, as percent of each element in the table order, for example `100,0` or `40,30,30`. The numbers must add up to 100. Leave End empty for a uniform layer.
- Function and p or f(x): how the mix moves from Start to End: P-, E- or S-FGM (type p in the box), or Custom (type the formula f(x) in the same box, see 5.9).
- Direction: Axial (along the wire, inside the layer) or Radial (from the core to the surface).
- Bins: into how many cells the layer is cut (a uniform layer needs none).
- a_z (nm): optional lattice constant of this layer along the wire (needs Vegard mode).

A graded interface is simply a thin graded layer between two others. Reversing a grading (Ni to Cu, then Cu to Ni) just means swapping Start and End. The buttons Add layer and Remove last change the number of layers (up to 12). Example buttons fill in ready-made stacks.

### 5.5 Combined: graded along + across

Give four corner compositions: core at the bottom, surface at the bottom, core at the top, surface at the top (percent per element, adding to 100). Give a grading function with p and a number of slices for the direction along the wire, and a function with p and a number of rings for the direction from the core to the surface. Each function may also be Custom: type the formula in the p or f(x) box (see 5.9). Every cell gets a blend of the four corners, so any mix is valid and no value can leave 0 to 100 %. If two corners are equal the wire is graded in one direction only. Length (nm) is the length of the wire.

### 5.6 Crystal layers: FCC(111) / HCP(0001)

Every layer is a close-packed crystal with its atomic planes across the wire: FCC with the (111) planes, HCP with the (0001) planes. Close-packed rows run in the same direction in all layers, so different crystals can be joined.

- Crystal: FCC (111) or HCP (0001).
- a and c (nm): lattice constants of that crystal. For HCP, c blank means the ideal ratio.
- Thickness (nm) and Composition % (one recipe per layer).
- In-plane spacing: Mean of the layers, Match first layer, or One value (typed nearest-neighbour distance in nm). All layers are strained to this common in-plane spacing; the table reports the strain of each layer. Each layer keeps its own spacing between atomic planes.

A single-crystal wire gets extra atomic planes when needed so that the stacking is continuous across the periodic boundary (a multiple of 3 planes for FCC, an even number for HCP). The message says so.

### 5.7 Spacing along the wire (Vegard)

With Same as in-plane the whole wire uses one lattice parameter. With Follows composition (Vegard) each layer (or slice) gets the spacing along the wire of its own composition, calculated as the atom-weighted mean of the element constants. The sideways spacing stays common, so layers are strained to fit each other. The distance between two neighbouring atomic layers is the mean of their own spacings, and the box ends one interface spacing above the top. In Vegard mode every element needs its own lattice constant in the table.

### 5.8 Result views

LAMMPS DATA FILE; Layer & cell table (or Cell table) with the atoms of each element in every cell and green PASS lines for the built-in exactness checks; Composition along the wire (atom percent of each element, measured on the finished atoms, with layer boundaries); Composition core to surface (combined mode); Side view (a thin slice through the wire axis with every element in its own colour).

### 5.9 Custom functions: pure wires and solid layers

P-, E- and S-FGM give smooth curves, so most cells are a RANDOM MIX of the elements. To get solid regions without mixing, choose Custom. In the Function column pick Custom and type the formula f(x) in the p or f(x) box. The letter x is the position inside the layer from 0 to 1: along the layer for Axial, from the core (0) to the surface (1) for Radial, and in the combined mode along the wire or from the core to the surface. The result is how far to go from the Start recipe (0) to the End recipe (1); it is cut to the range 0 to 1.

- f(x) = 0 gives the pure Start material, a pure wire with no mixing. f(x) = 1 gives the pure End material.
- f(x) = x>0.7 (or step(x,0.7)) is 0 below 0.7 and 1 from 0.7 on. With Direction Radial this is a SOLID OUTER SHELL around a pure core. With Axial it is a solid piece at the end of the layer.
- f(x) = box(x,0.3,0.7) is a solid band between 0.3 and 0.7 with the Start material elsewhere.
- f(x) = 0.4 is neither 0 nor 1, so every cell is a random alloy with 40 % End material. Only the values 0 and 1 give solid regions.

The formula is evaluated once per cell (bin), in the middle of the cell, so the edge of a solid region lies on a cell edge. With Bins 10 and f(x) = x>0.7 the edge is at 70 % of the radius and the shell is three rings (30 %) thick. Use more bins for a thinner or more exact shell. A cell cannot be thinner than about half an atomic layer, otherwise it is empty and the program says so. Check the Layer & cell table (every cell of a solid region shows atoms of ONE element only) and the Side view.

In the combined mode both functions can be Custom. The Solid regions example uses the corners Ni, Co, Co, Co (core bottom first) with x>=0.5 in both directions and 2 slices and 2 rings: only the lower core is Ni, everything else is Co.

In the Crystal layers form every layer has one composition, so it is already solid. Graded crystal layers (with Custom) are available only through the JSON description of the command line.

## 6. Using the files in LAMMPS

- Units metal, atom_style atomic, and `boundary s s p` (free surfaces across the wire, periodic along it). If you added vacuum in x and y you may also use `p p p`.
- Atom type k is element k of the table. The masses are in the file. In `pair_coeff` list the elements in the same order.
- The nearest-neighbour distance, the lattice constant and the relaxation depend on your potential. The structure is a starting point: relax it (for example `minimize`, and `fix box/relax z 0.0` for zero stress along the wire) before deforming it.
- Thin wires (2 nm diameter) shorten by a few percent when relaxed because of surface stress.
- The text files `in.geometry_check` and `in.relax_check` that came with the program are small LAMMPS inputs that check a data file: the first needs no potential.

## 7. Command line

```
py FGM_NW_Builder.py build --direction axial --bins 10 --seed 1 --out wire.nw
py FGM_NW_Builder.py layered --example > spec.json
py FGM_NW_Builder.py layered --spec spec.json --out wire.nw
py FGM_NW_Builder.py combined --example > spec.json
py FGM_NW_Builder.py crystals --example > spec.json
py FGM_NW_Builder.py validate --ref reference.nw --direction axial
py FGM_NW_Builder.py selftest
```

`build --help` lists all options. The JSON examples show every field of the other modes. `selftest` runs the built-in checks (a few seconds) and must end with ALL CHECKS PASSED.

## 8. Worked examples

### 8.1 A Ni-Cu graded cylinder (first tab)

Shape Cylindrical, Radius 2, Length 20, Grading Direction Radial, Grading Function P-FGM, p 1, masses 58.69 and 63.55, FCC, lattice parameters 0.352 and 0.361, bins 10. Press Build. Look at Bin table: material 2 rises from the core to the surface. Save the file.

### 8.2 A Ni/Co multilayer with graded interfaces (second tab)

Structure type Layers. Press the Ni/Co example. Layer 1 is pure Ni (4 nm), layer 2 grades from Ni to Co with an S-FGM (6 nm), layer 3 is pure Co, layer 4 grades back. Press Build and open Composition along the wire and Side view.

### 8.3 A Ti/Al stack of different crystals

Structure type Crystal layers. Press the Ti/Al example (Ti hexagonal, a 0.2951 nm, c 0.4683 nm; Al cubic, a 0.405 nm). The table shows strains of about -1.5 % and +1.5 %. Build and save. In LAMMPS use a potential that covers both elements and relax the structure.

### 8.4 A solid cobalt shell around a nickel wire

Second tab, Structure type Layers. Press the Solid shell example. Layer 2 has Start 100,0 (Ni), End 0,100 (Co), Function Custom, formula x>0.7, Direction Radial, Bins 10. Press Build: the Layer & cell table shows rings that contain only Ni and the outer three rings that contain only Co, and the Side view shows a red shell around a blue core. To make the shell thinner use more bins, for example 20 bins and x>0.85. For a completely pure wire type 0 as the formula.

## 9. Messages and troubleshooting

- 'Python was not found': start with `py` instead of `python`.
- 'NumPy is missing': accept the offer to install it, or run `py -m pip install numpy`.
- 'the custom function f(x) = ... could not be evaluated': the formula has a mistake, uses a name that is not allowed, or fails at x = 0, 0.25, 0.5, 0.75 or 1 (for example a division by zero). Fix the formula.
- A red message after Build tells you which field is wrong (for example a composition adding to 90). Fix the field and press Build again.
- 'Length set to N whole unit cells': the corrected mode and the second tab round the length so that a periodic wire has no seam.
- 'Layer thicknesses were rounded': thicknesses are whole atomic layers.
- 'Warning: empty cells ...': you asked for more bins, slices or rings than the wire has atomic layers or atoms. Use fewer.
- 'Note: N cells hold fewer than 20 atoms': the composition of small cells is rounded coarsely. Use fewer bins or a larger wire.
- 'A layer is strained by more than 5 %': the crystals of neighbouring layers differ a lot. Treat the structure as a rough starting point.
- A very large structure can take a while and a lot of memory. The side view shows at most a few thousand atoms.

## 10. License and credits

MIT License

Copyright (c) 2026 Md. Shabbir Rahman Siam

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## 11. Limitations and invalid values

### 11.1 Why there is no sphere in the second tab

The second tab stacks layers along the axis of a wire, or grades along the axis and from the axis outwards. A sphere has no length and no single stacking direction. Close-packed crystal planes, for example, can only be stacked in one chosen direction, and a sphere cut from such a stack would have no well-defined layers. The sphere is therefore available only in the first tab, with radial grading (the bins are spherical shells) or axial grading along the diameter. Core-shell nanoparticles with different crystals would need a different method and are not included.

### 11.2 Other limits of the structures

- One lattice type (FCC or BCC) and one in-plane spacing for the whole wire in the Layers and Combined modes. Different crystals per layer exist only as FCC(111) and HCP(0001) layers, with close-packed rows parallel in all layers. BCC layers, other orientations and a twist between layers are not available.
- The structures are perfect crystals with random mixing inside each cell: no vacancies, dislocations, twins, stacking faults or grain boundaries, and no chemical short-range order. Mismatched layers are strained to fit; real interfaces relax by misfit dislocations. Relax the structure in LAMMPS.
- Composition is rounded to whole atoms in every cell. Cells with few atoms give lumpy compositions.
- Layer thicknesses are whole atomic layers. The wire length is made a whole number of unit cells so that a periodic z boundary is seamless.
- Only cylinders and cuboids in the second tab; wires are built along z.
- The edge of a solid region made with a Custom function lies on a cell edge (ring or slice), so its position is as fine as the number of bins allows, and never finer than the atomic layers or atom rings of the wire. Values strictly between 0 and 1 always give a random mix, never a sharp edge.
- Some check scripts and plots assume metals; potentials are not provided with the program.
- The program builds structures. It does not check that your potential can describe your structure.

### 11.3 Allowed values

```
Radius, X, Y, Length, Thickness      number larger than 0 (nm)
Mass, lattice parameter a, c         number larger than 0
p                                    number 0 or larger
Number of bins (first tab)           whole number, 2 or more
Bins, slices, rings (second tab)     whole number, 1 or more
Vacuum                               0 or larger (Angstrom)
Random seed                          whole number, or empty
Composition (second tab)             one number per element, none negative, adding to 100 (+- 0.01)
Custom formula (f(x) box)            must give a number at x = 0, 0.25, 0.5, 0.75 and 1; allowed names in 4.3
Elements                             2 to 5, symbols not empty and all different
Layers                               1 to 12
a_z (layer lattice constant)         only with Spacing = Follows composition (Vegard)
Vegard mode                          every element needs its own lattice constant
In-plane spacing 'One value'         needs the nearest-neighbour distance in nm
```

What happens otherwise: the program shows a message that names the field and does not build. Other points:

- Custom formulas can only use x, numbers, + - * / ^, comparisons, and / or, and the functions and constants listed in 4.3. Any other name is refused before the build. A formula that gives a value below 0 or above 1 is cut to 0 or 1, and a formula must give a number at x = 0, 0.25, 0.5, 0.75 and 1.
- In the first tab the number of bins is limited by the size of the wire in the standard mode (a bin cannot be thinner than half a lattice constant); the program lowers the number and says so.
- Compare with reference file works only on the first tab, in the standard mode, without vacuum.
- Elements in the second tab are limited to 5 and the layer tables to 12 rows.

### 11.4 What has been checked

The first tab was compared with an independent implementation of the published algorithm on many random settings. The other modes are checked by the built-in tests (exact atom counts in every cell, crystal geometry, spacings, periodic seams) and by spot checks in LAMMPS. They have not been tested with every possible potential or size. Always inspect your own structure and relax it before a long simulation.
"""

try:
    import numpy as np
except ImportError:  # handled in ensure_numpy()
    np = None


def ensure_numpy(interactive):
    """Make sure NumPy is available; in the window version offer to install it."""
    global np
    if np is not None:
        return True
    msg = "NumPy is not installed. It can be installed with:\n    %s -m pip install numpy" % (sys.executable,)
    if not interactive:
        print(msg)
        return False
    root = tk.Tk()
    root.withdraw()
    if messagebox.askyesno("NumPy is missing", "This program needs NumPy.\n\nInstall it now? (internet needed once)"):
        import subprocess
        subprocess.call([sys.executable, "-m", "pip", "install", "numpy"])
        try:
            import numpy as _np
            np = _np
        except ImportError:
            messagebox.showerror("Install failed", "Could not install NumPy automatically.\n\n" + msg)
    root.destroy()
    return np is not None


# =====================================================================================
# PART 1 - the builder 
# =====================================================================================
SHAPES = {"cylindrical": 1, "cuboid": 2, "sphere": 3}


# ----------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------
def round_half_up(v):
    """MATLAB/Octave round(): 0.5 rounds up (Python's round() would round to even)."""
    return int(math.floor(v + 0.5))


def _step(x, a):
    """0 below a, 1 from a upwards: a sharp edge, which makes a SOLID region instead of a mix."""
    return 1.0 if x >= a else 0.0


def _box(x, a, b):
    """1 between a (included) and b (not included), 0 elsewhere: a solid band."""
    return 1.0 if (x >= a and x < b) else 0.0


_CUSTOM_NAMES = {"exp": math.exp, "log": math.log, "sqrt": math.sqrt, "sin": math.sin, "cos": math.cos, "tan": math.tan,
                 "tanh": math.tanh, "sinh": math.sinh, "cosh": math.cosh, "atan": math.atan, "asin": math.asin,
                 "acos": math.acos, "pi": math.pi, "e": math.e, "abs": abs, "min": min, "max": max,
                 "floor": math.floor, "ceil": math.ceil, "pow": pow, "step": _step, "box": _box}


def eval_custom(b, custom):
    """Value of the typed formula f(x) at x = b. Octave-style text such as x^2, x.^2 or && is accepted."""
    expr = str(custom).strip().replace(".^", "**").replace("^", "**").replace(".*", "*").replace("./", "/")
    expr = expr.replace("&&", " and ").replace("||", " or ")
    if not expr:
        raise ValueError("the custom function is empty")
    try:
        x = eval(expr, {"__builtins__": {}}, dict(_CUSTOM_NAMES, x=b))  # your own formula, run on your own PC
        return x if isinstance(x, complex) else float(x)
    except Exception as e:
        raise ValueError("the custom function f(x) = %s could not be evaluated at x = %g: %s" % (custom, b, e))


def grading_fraction(b, kind, p, custom=""):
    """Grading function: position b (0..1) -> mass fraction x (0..1) of material 2."""
    if kind == "power":  # P-FGM
        x = b ** p
    elif kind == "exp":  # E-FGM
        x = 1.0 - math.exp(-(b ** p))
    elif kind == "sigmoid":  # S-FGM
        if b <= 0.5:
            x = 0.5 * (2 * b) ** p
        else:
            x = 1.0 - 0.5 * (2 * (1 - b)) ** p
    elif kind == "custom":
        x = eval_custom(b, custom)
        if x != x:
            raise ValueError("the custom function f(x) = %s gave 'not a number' at x = %g" % (custom, b))
    else:
        raise ValueError("unknown grading kind: %s" % kind)
    if isinstance(x, complex):  # e.g. negative number to a fractional power
        x = 1.0
    return min(max(float(x), 0.0), 1.0)


def check_grading(kind, p, custom, what):
    """Refuse a bad function before the build starts (a typed formula is tried at five positions)."""
    if kind not in ("power", "exp", "sigmoid", "custom"):
        raise ValueError("%s: unknown grading function %r" % (what, kind))
    if kind == "custom":
        for b in (0.0, 0.25, 0.5, 0.75, 1.0):
            try:
                grading_fraction(b, kind, p, custom)
            except ValueError as e:
                raise ValueError("%s: %s" % (what, e))


def fn_label(kind, p, custom=""):
    """Short text for tables."""
    if kind == "custom":
        return "Custom f(x)=%s" % custom
    return "%s p=%g" % ({"power": "P-FGM", "exp": "E-FGM", "sigmoid": "S-FGM"}[kind], p)


def count_type2(x, n_atoms, m1, m2, fraction="mass"):
    """Fraction x of material 2 -> how many of n_atoms atoms must be type 2.
    fraction="mass": x is a mass fraction (nanoHUB). fraction="atom": x is an atom fraction."""
    if fraction == "atom":
        ratio = x
    else:
        ratio = (x / m2) / (x / m2 + (1.0 - x) / m1)
    k = round_half_up(ratio * n_atoms)
    return min(max(k, 0), n_atoms)


USE_PREFILTER = True  # only generate the part of the crystal block that is needed


def make_lattice(kind, N, a, ranges=None):
    """All lattice sites of an N x N x N block of cubic cells, in the same order as
    the original MATLAB meshgrid code (so atom IDs line up with nanoHUB).
    ranges = ((xlo,xhi),(ylo,yhi),(zlo,zhi)) skips sites far outside the wire; the sites
    that are kept stay in exactly the same order, so the result is identical but the
    computer needs far less memory for long wires."""
    k = np.arange(N + 1)
    corners = k * a
    half = 0.5 * a + k * a  # face-centre / body-centre positions

    def pick(arr, axis):
        if ranges is None:
            return arr
        lo, hi = ranges[axis]
        return arr[(arr >= lo) & (arr <= hi)]

    def grid(xs, ys, zs):
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="xy")  # same shape rule as MATLAB
        return np.column_stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")])

    Cx, Cy, Cz = (pick(corners, i) for i in range(3))
    Hx, Hy, Hz = (pick(half, i) for i in range(3))
    if kind == "fcc":
        parts = [grid(Cx, Cy, Cz), grid(Hx, Hy, Cz), grid(Cx, Hy, Hz), grid(Hx, Cy, Hz)]
    else:  # bcc
        parts = [grid(Cx, Cy, Cz), grid(Hx, Hy, Hz)]
    return np.vstack(parts)


# ----------------------------------------------------------------------------------
# the builder
# ----------------------------------------------------------------------------------
# ---- FROZEN REFERENCE: validated against the original nanoHUB main.m. Do not edit. ----
# The program uses engine_build(); this function is kept only so that "selftest" can prove
# the engine gives exactly the same structures. Known difference: for an axial S-FGM wire with a
# fractional p the last bin is a complex-number case; the engine reproduces the original there,
# this older function does not (the selftest leaves that one case out of its comparison).
def reference_build_structure(shape="cylindrical", radius_nm=2.0, length_nm=20.0, xdim_nm=2.0,
                    ydim_nm=2.0, direction="radial", grading="power", p=1.0, custom="",
                    mass1=63.55, mass2=58.69, lattice="fcc", a1_nm=0.36, a2_nm=0.36,
                    bins=10, seed=None, mode="nanohub", sample="center", fraction="mass",
                    vacuum_xy=0.0, vacuum_z=0.0):
    """mode="nanohub": exact copy of main.m.   mode="fixed": corrected version.
    sample (fixed only): "center" = grading evaluated at the middle of each bin,
                         "edge"   = first bin gets b=0, last bin gets b=1.
    fraction (fixed only): "mass" (as nanoHUB) or "atom" = what the grading value means."""
    fixed = mode == "fixed"
    rng = np.random.default_rng(seed)
    dim = SHAPES[shape]
    axial = direction == "axial"
    messages = []
    frac_kind = fraction if fixed else "mass"

    # --- sizes (nm -> Angstrom) ---
    leng = length_nm * 10 if dim != 3 else radius_nm * 20
    if dim != 2:
        r_base = radius_nm * 10
        x_len = 2 * r_base
        y_len = 2 * r_base
        z_len = leng
    else:
        x_len = 10 * xdim_nm
        y_len = 10 * ydim_nm
        z_len = leng
        r_base = math.sqrt(x_len ** 2 + y_len ** 2) / 2
    lat_p = 10 * (a1_nm + a2_nm) / 2  # one averaged lattice constant, like nanoHUB

    # --- nanoHUB only: bins cannot be thinner than half a lattice constant ---
    if not fixed:
        limit = (2 * leng if axial else 2 * r_base) / lat_p
        if bins > limit:
            bins = round_half_up(limit)
            messages.append("Warning: number of bins too high, reset to %d." % bins)
    if bins < 2:
        raise ValueError("need at least 2 bins")

    # --- step 1: the crystal block ---
    N = math.ceil(max(x_len, y_len, z_len) / lat_p)
    cmax = 0.5 * lat_p + N * lat_p  # largest coordinate of nanoHUB's cubic block

    z_cut = z_len
    n_cells_z = None
    if fixed and dim != 3:
        # keep a whole number of unit cells along z, so periodic z has no stacking seam
        n_cells_z = max(1, round_half_up(z_len / lat_p))
        z_cut = (2 * n_cells_z - 1) * lat_p / 2 + 1e-6
        if abs(n_cells_z * lat_p - z_len) > 1e-6:
            messages.append("Length set to %d whole unit cells = %.2f A (asked for %.2f A)."
                            % (n_cells_z, n_cells_z * lat_p, z_len))

    if USE_PREFILTER:
        c = cmax / 2
        if dim == 2:
            ranges = ((-lat_p, x_len + lat_p), (-lat_p, y_len + lat_p), (-lat_p, z_cut + lat_p))
        else:
            rr = radius_nm * 10 + lat_p
            ranges = ((c - rr, c + rr), (c - rr, c + rr),
                      (c - rr, c + rr) if dim == 3 else (-lat_p, z_cut + lat_p))
    else:
        ranges = None
    pos = make_lattice(lattice, N, lat_p, ranges)
    X, Y, Z = pos[:, 0], pos[:, 1], pos[:, 2]

    # --- step 2: carve the shape ---
    axis = None  # true centre line / centre point of the carved shape
    if dim == 1:
        cx = cy = cmax / 2
        outside = (X - cx) ** 2 + (Y - cy) ** 2 - (radius_nm * 10) ** 2 > 0
        delete = outside | (Z > z_cut)
        axis = np.array([cx, cy, 0.0])
    elif dim == 3:
        cx = cy = cz = cmax / 2
        delete = (X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2 - (radius_nm * 10) ** 2 > 0
        axis = np.array([cx, cy, cz])
    else:
        delete = (X > x_len) | (Y > y_len) | (Z > z_cut)
    pos = pos[~delete]
    shift = pos.min(axis=0)
    pos = pos - shift  # move the corner of the wire to (0,0,0)
    if dim == 2:
        axis = pos.max(axis=0) / 2 + shift  # centre of the cross-section
    axis = axis - shift
    n_atom = len(pos)
    types = np.ones(n_atom, dtype=int)

    # --- step 3: decide which bin every atom belongs to ---
    if not fixed:
        # nanoHUB way: overlapping closed ranges (an atom on an edge can sit in two bins)
        if axial:
            v = pos[:, 2]
            n_zbox = bins - 1
            step = (v.max() - v.min()) / n_zbox
            lo, hi, denom = v.min(), v.min() + step, n_zbox
        else:
            xc, yc, zc = x_len / 2, y_len / 2, z_len / 2
            if dim == 3:
                v = np.sqrt((pos[:, 0] - xc) ** 2 + (pos[:, 1] - yc) ** 2 + (pos[:, 2] - zc) ** 2)
            else:
                v = np.sqrt((pos[:, 0] - xc) ** 2 + (pos[:, 1] - yc) ** 2)
            step = r_base / (bins - 1)
            lo, hi, denom = 0.0, step, bins
        members = []
        for i in range(1, bins + 1):
            members.append(np.nonzero((v >= lo) & (v <= hi))[0])
            lo = hi
            hi = hi + step
        edges = [(j * step, (j + 1) * step) for j in range(bins)]
        b_list = [i / denom for i in range(1, bins + 1)]
    else:
        # fixed way: each atom in exactly one bin
        if axial:
            half = lat_p / 2
            layer = np.rint(pos[:, 2] / half).astype(int)  # atomic layer number 0..K-1
            K = int(layer.max()) + 1
            bin_id = (layer * bins) // K
            v = pos[:, 2]
            edges = [(j * K / bins * half, (j + 1) * K / bins * half) for j in range(bins)]
            if K < bins:
                messages.append("Warning: only %d atomic layers for %d bins; some bins are empty." % (K, bins))
        else:
            d = pos - axis
            v = np.sqrt((d ** 2).sum(axis=1)) if dim == 3 else np.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2)
            r_max = float(v.max())
            width = r_max / bins
            bin_id = np.minimum(np.floor(v / width).astype(int), bins - 1)
            edges = [(j * width, (j + 1) * width) for j in range(bins)]
        members = [np.nonzero(bin_id == j)[0] for j in range(bins)]
        if sample == "edge":
            b_list = [j / (bins - 1) for j in range(bins)]
        else:
            b_list = [(j + 0.5) / bins for j in range(bins)]

    # --- step 4: repaint atoms of each bin ---
    bin_info = []
    yy = [0.0]  # "modelled" mass-fraction curve (first entry is 0, as in nanoHUB)
    for j in range(bins):
        idx = members[j]
        n_in = len(idx)
        b = b_list[j]
        x = grading_fraction(b, grading, p, custom)
        k = count_type2(x, n_in, mass1, mass2, frac_kind)
        den = (n_in - k) * mass1 + k * mass2
        yy.append(k * mass2 / den if den > 0 else float("nan"))
        if k > 0:
            types[idx[rng.choice(n_in, k, replace=False)]] = 2
        bin_info.append({"bin": j + 1, "lo": edges[j][0], "hi": edges[j][1], "b": b, "x": x,
                         "n_atoms": n_in, "n_type2_requested": k})
    if fixed:
        empty = [bi["bin"] for bi in bin_info if bi["n_atoms"] == 0]
        if empty:
            messages.append("Warning: empty bins %s (too many bins for this size)." % empty)

    atoms = np.column_stack([np.arange(1, n_atom + 1), types, pos])

    def span(c, pad):
        lo_, hi_ = pos[:, c].min(), pos[:, c].max()
        return (lo_ - pad, hi_ + pad) if pad > 0 else (lo_, hi_ + lat_p / 2)

    box = {"x": span(0, vacuum_xy), "y": span(1, vacuum_xy), "z": span(2, vacuum_z)}
    return {"atoms": atoms, "box": box, "lat_p": lat_p, "bins": bins, "bin_info": bin_info,
            "yy": yy, "coord": v, "messages": messages, "mass1": mass1, "mass2": mass2,
            "axial": axial, "mode": mode, "members": members, "axis": axis,
            "n_cells_z": n_cells_z, "lattice": lattice}


def verify(res):
    """Self-checks for the corrected mode (they cannot hold for the nanoHUB copy)."""
    atoms, info = res["atoms"], res["bin_info"]
    types = atoms[:, 1].astype(int)
    out = []
    out.append(("every atom is in exactly one bin",
                sum(len(m) for m in res["members"]) == len(atoms)
                and len(np.unique(np.concatenate(res["members"]))) == len(atoms)))
    n_el = len(info[0]["counts"]) if info and "counts" in info[0] else 2
    if info and "counts" in info[0]:
        ok = all([int((types[m] == e + 1).sum()) for e in range(n_el)] == bi["counts"]
                 for m, bi in zip(res["members"], info))
    else:
        ok = all(int((types[m] == 2).sum()) == bi["n_type2_requested"] for m, bi in zip(res["members"], info))
    out.append(("type-2 count in every bin equals the requested count" if n_el == 2
                else "the count of every element in every bin equals the requested count", ok))
    if res["n_cells_z"] is not None:
        out.append(("wire is a whole number of unit cells along z (no periodic seam)", True))
    return out


def _write_data(f, result, max_atom_lines=None):
    box, atoms = result["box"], result["atoms"]
    masses = result.get("masses") or [result["mass1"], result["mass2"]]
    f.write("%s\n\n" % result.get("comment", FILE_COMMENT))
    f.write("\t\t%4i atoms\n" % len(atoms))
    f.write("\t\t%d atom types\n" % len(masses))
    for ax in "xyz":
        f.write("\t%12.8f %12.8f \t%slo %shi\n" % (box[ax][0], box[ax][1], ax, ax))
    f.write("\nMasses\n\n")
    for i, m in enumerate(masses, 1):
        f.write("\t\t%d %12.8f\n" % (i, m))
    f.write("\nAtoms\n\n")
    np.savetxt(f, atoms if max_atom_lines is None else atoms[:max_atom_lines], fmt="%d %d %f %f %f")


def write_lammps_data(path, result):
    """Write the LAMMPS data file in exactly the layout nanoHUB uses."""
    with open(path, "w", newline="\n") as f:
        _write_data(f, result)


def lammps_data_text(result, max_atom_lines=None):
    """The same file as text (the atom list can be cut short for an on-screen preview)."""
    import io
    buf = io.StringIO()
    _write_data(buf, result, max_atom_lines)
    return buf.getvalue()


# ----------------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------------
def add_arguments(ap):
    ap.add_argument("--shape", choices=list(SHAPES), default="cylindrical")
    ap.add_argument("--radius", type=float, default=2.0, help="nm (cylinder, sphere)")
    ap.add_argument("--length", type=float, default=20.0, help="nm (cylinder, cuboid)")
    ap.add_argument("--xdim", type=float, default=2.0, help="nm (cuboid)")
    ap.add_argument("--ydim", type=float, default=2.0, help="nm (cuboid)")
    ap.add_argument("--direction", choices=["radial", "axial"], default="radial")
    ap.add_argument("--grading", choices=["power", "exp", "sigmoid", "custom"], default="power")
    ap.add_argument("--p", type=float, default=1.0, help="grading parameter")
    ap.add_argument("--custom", default="", help='custom f(x), e.g. "x^2"')
    ap.add_argument("--mass1", type=float, default=63.55)
    ap.add_argument("--mass2", type=float, default=58.69)
    ap.add_argument("--lattice", choices=["fcc", "bcc"], default="fcc")
    ap.add_argument("--a1", type=float, default=0.36, help="lattice parameter 1, nm")
    ap.add_argument("--a2", type=float, default=0.36, help="lattice parameter 2, nm")
    ap.add_argument("--bins", type=int, default=10)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--mode", choices=["standard", "corrected", "nanohub", "fixed"], default="standard",
                    help="standard = validated reference behaviour (default); corrected = exact equal bins")
    ap.add_argument("--sample", choices=["center", "edge"], default="center",
                    help="fixed mode: grading evaluated at bin centre, or first bin b=0 / last bin b=1")
    ap.add_argument("--fraction", choices=["mass", "atom"], default="mass",
                    help="fixed mode: is the grading value a mass fraction or an atom fraction")
    ap.add_argument("--vacuum-xy", type=float, default=0.0, help="empty space on each side in x and y, Angstrom")
    ap.add_argument("--vacuum-z", type=float, default=0.0, help="empty space on each end in z, Angstrom")


def params_from_args(a):
    return dict(shape=a.shape, radius_nm=a.radius, length_nm=a.length, xdim_nm=a.xdim,
                ydim_nm=a.ydim, direction=a.direction, grading=a.grading, p=a.p,
                custom=a.custom, mass1=a.mass1, mass2=a.mass2, lattice=a.lattice,
                a1_nm=a.a1, a2_nm=a.a2, bins=a.bins, seed=a.seed, mode={"standard": "nanohub", "corrected": "fixed"}.get(a.mode, a.mode),
                sample=a.sample, fraction=a.fraction, vacuum_xy=a.vacuum_xy,
                vacuum_z=a.vacuum_z)


def layered_summary_lines(res):
    """Layer table and per-cell table of a layered wire."""
    els = res["elements"]
    names = [str(e.get("symbol", i + 1)) for i, e in enumerate(els)]
    lines = list(res["messages"])
    lines.append("mode: layered   atoms: %d   box: %s" % (len(res["atoms"]),
                 {k: tuple(round(float(x), 4) for x in v) for k, v in res["box"].items()}))
    unit = "atom%" if res.get("basis") == "atom" else "mass%"
    lines.append("")
    lines.append("LAYERS (bottom to top), composition in %s of: %s" % (unit, " / ".join(names)))
    lines.append("%5s %7s %9s %9s %8s %8s %-8s %5s  %s" % ("layer", "planes", "asked nm", "built nm", "a_z nm", "atoms", "grading", "bins", "start -> end"))
    for ly in res["layers"]:
        fmt = lambda r: "/".join("%g" % round(100 * x, 2) for x in r)
        grad = "uniform" if ly["uniform"] else "%s %s" % (fn_label(ly["grading"], ly["p"], ly.get("custom", "")), ly["direction"])
        comp = fmt(ly["start"]) if ly["uniform"] else "%s -> %s" % (fmt(ly["start"]), fmt(ly["end"]))
        lines.append("%5d %7d %9.2f %9.2f %8.4f %8d %-8s %5d  %s" % (ly["layer"], ly["planes"], ly["asked_nm"], ly["built_nm"], ly["a_z_nm"],
                                                              ly["atoms"], "" if ly["uniform"] else "graded", ly["bins"],
                                                              comp if ly["uniform"] else comp + "   (" + grad + ")"))
    lines.append("")
    lines.append("CELLS (atoms of each element)")
    lines.append("%4s %5s %20s %8s " % ("cell", "layer", "range (A)", "atoms") + " ".join("%8s" % nme for nme in names))
    for bi in res["bin_info"]:
        lines.append("%4d %5d %9.2f..%9.2f %8d " % (bi["bin"], bi["layer"], bi["lo"], bi["hi"], bi["n_atoms"])
                     + " ".join("%8d" % c for c in bi["counts"]))
    for name, ok in verify(res):
        lines.append("[%s] %s" % ("PASS" if ok else "FAIL", name))
    return lines


def summary_lines(res):
    """Text summary of a build: warnings, box, per-bin table and (fixed mode) self-checks."""
    if res.get("crystal_layers"):
        return crystal_summary_lines(res)
    if res.get("layers"):
        return layered_summary_lines(res)
    if res.get("combined"):
        return combined_summary_lines(res)
    lines = list(res["messages"])
    lines.append("mode: %s   atoms: %d   box: %s" % ({"nanohub": "standard", "fixed": "corrected"}.get(res["mode"], res["mode"]), len(res["atoms"]),
                 {k: tuple(round(float(x), 4) for x in v) for k, v in res["box"].items()}))
    els = res.get("elements") or []
    if len(els) > 2:
        names = [str(e.get("symbol", i + 1)) for i, e in enumerate(els)]
        lines.append("%4s %20s %8s " % ("bin", "range (A)", "atoms") + " ".join("%8s" % nme for nme in names))
        for bi in res["bin_info"]:
            lines.append("%4d %9.3f..%9.3f %8d " % (bi["bin"], bi["lo"], bi["hi"], bi["n_atoms"])
                         + " ".join("%8d" % c for c in bi["counts"]))
    else:
        lines.append("%4s %20s %8s %8s %10s" % ("bin", "range (A)", "atoms", "type2", "mass%2"))
        for bi, y in zip(res["bin_info"], res["yy"][1:]):
            lines.append("%4d %9.3f..%9.3f %8d %8d %10.2f" % (bi["bin"], bi["lo"], bi["hi"], bi["n_atoms"],
                                                             bi["n_type2_requested"], 100 * y))
    if res["mode"] == "fixed":
        for name, ok in verify(res):
            lines.append("[%s] %s" % ("PASS" if ok else "FAIL", name))
    return lines


def cli_build(argv=None):
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py build", description="FGM Nanostructure Builder (standard and corrected modes)")
    add_arguments(ap)
    ap.add_argument("--out", default="fgm_nw.data")
    a = ap.parse_args(argv)
    res = build_structure(**params_from_args(a))
    write_lammps_data(a.out, res)
    print("\n".join(summary_lines(res)))
    print("written:", a.out)


# =====================================================================================
# PART 1b - the engine
#   One pipeline for every mode:  geometry -> cells -> composition recipe -> atom counts -> atoms
#   * works for any number of elements (2 = binary, 3 = ternary, 4-5 = high-entropy alloys)
#   * the nanoHUB-copy mode runs THROUGH this engine (partition="legacy") and is proven to give
#     exactly the same result as the frozen reference function above (see "selftest").
# =====================================================================================
FILE_COMMENT = "LAMMPS data file generated by FGM Nanostructure Builder - Md. Shabbir Rahman Siam"


def grading_raw(b, kind, p, custom=""):
    """Grading value exactly as the nanoHUB code computes it: no clamping, and it may be a
    complex number (a negative base to a fractional power, which only happens for S-FGM in the
    last, partly empty bin of an axial wire)."""
    if kind == "power":
        return b ** p
    if kind == "exp":
        return 1.0 - math.exp(-(b ** p))
    if kind == "sigmoid":
        if b <= 0.5:
            return 0.5 * (2 * b) ** p
        return 1.0 - 0.5 * (2 * (1 - b)) ** p
    return grading_fraction(b, kind, p, custom)


def oct_round(v):
    """Octave round(): halves round away from zero."""
    return math.copysign(math.floor(abs(v) + 0.5), v)


def legacy_counts_complex(x, n, m1, m2):
    """What the original program does when the grading value x is a complex number.
    Octave rounds the real and imaginary parts separately, compares by size (modulus),
    and the random pick then uses only the real part. Returns (atoms of element 2, mass% curve value)."""
    ratio = (x / m2) / (x / m2 + (1 - x) / m1)
    v = ratio * n
    k = complex(oct_round(v.real), oct_round(v.imag))
    if abs(k) > n:
        k = complex(n, 0.0)
    den = (n - k) * m1 + k * m2
    y = (k * m2 / den).real if abs(den) > 0 else float("nan")
    return min(max(int(k.real), 0), n), y


def sizes_for(shape, radius_nm, length_nm, xdim_nm, ydim_nm):
    """nm -> Angstrom sizes of the wire."""
    dim = SHAPES[shape]
    leng = length_nm * 10 if dim != 3 else radius_nm * 20
    if dim != 2:
        r_base = radius_nm * 10
        x_len = 2 * r_base
        y_len = 2 * r_base
        z_len = leng
    else:
        x_len = 10 * xdim_nm
        y_len = 10 * ydim_nm
        z_len = leng
        r_base = math.sqrt(x_len ** 2 + y_len ** 2) / 2
    return dim, leng, r_base, x_len, y_len, z_len


def make_geometry(shape, radius_nm, length_nm, xdim_nm, ydim_nm, lattice, lat_p, whole_cells_z):
    """Crystal block -> carved shape. Returns the atom positions (wire corner at 0,0,0)."""
    dim, leng, r_base, x_len, y_len, z_len = sizes_for(shape, radius_nm, length_nm, xdim_nm, ydim_nm)
    messages = []
    N = math.ceil(max(x_len, y_len, z_len) / lat_p)
    cmax = 0.5 * lat_p + N * lat_p  # largest coordinate of nanoHUB's cubic block
    z_cut = z_len
    n_cells_z = None
    if whole_cells_z and dim != 3:
        # keep a whole number of unit cells along z, so periodic z has no stacking seam
        n_cells_z = max(1, round_half_up(z_len / lat_p))
        z_cut = (2 * n_cells_z - 1) * lat_p / 2 + 1e-6
        if abs(n_cells_z * lat_p - z_len) > 1e-6:
            messages.append("Length set to %d whole unit cells = %.2f A (asked for %.2f A)."
                            % (n_cells_z, n_cells_z * lat_p, z_len))
    if USE_PREFILTER:
        c = cmax / 2
        if dim == 2:
            ranges = ((-lat_p, x_len + lat_p), (-lat_p, y_len + lat_p), (-lat_p, z_cut + lat_p))
        else:
            rr = radius_nm * 10 + lat_p
            ranges = ((c - rr, c + rr), (c - rr, c + rr),
                      (c - rr, c + rr) if dim == 3 else (-lat_p, z_cut + lat_p))
    else:
        ranges = None
    pos = make_lattice(lattice, N, lat_p, ranges)
    X, Y, Z = pos[:, 0], pos[:, 1], pos[:, 2]
    axis = None  # true centre line / centre point of the carved shape
    if dim == 1:
        cx = cy = cmax / 2
        outside = (X - cx) ** 2 + (Y - cy) ** 2 - (radius_nm * 10) ** 2 > 0
        delete = outside | (Z > z_cut)
        axis = np.array([cx, cy, 0.0])
    elif dim == 3:
        cx = cy = cz = cmax / 2
        delete = (X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2 - (radius_nm * 10) ** 2 > 0
        axis = np.array([cx, cy, cz])
    else:
        delete = (X > x_len) | (Y > y_len) | (Z > z_cut)
    pos = pos[~delete]
    shift = pos.min(axis=0)
    pos = pos - shift  # move the corner of the wire to (0,0,0)
    if dim == 2:
        axis = pos.max(axis=0) / 2 + shift  # centre of the cross-section
    axis = axis - shift
    return {"pos": pos, "axis": axis, "messages": messages, "n_cells_z": n_cells_z}


def partition_legacy(pos, dim, axial, bins, r_base, x_len, y_len, z_len):
    """nanoHUB way: overlapping closed ranges (an atom exactly on an edge can sit in two bins)."""
    if axial:
        v = pos[:, 2]
        n_zbox = bins - 1
        step = (v.max() - v.min()) / n_zbox
        lo, hi, denom = v.min(), v.min() + step, n_zbox
    else:
        xc, yc, zc = x_len / 2, y_len / 2, z_len / 2
        if dim == 3:
            v = np.sqrt((pos[:, 0] - xc) ** 2 + (pos[:, 1] - yc) ** 2 + (pos[:, 2] - zc) ** 2)
        else:
            v = np.sqrt((pos[:, 0] - xc) ** 2 + (pos[:, 1] - yc) ** 2)
        step = r_base / (bins - 1)
        lo, hi, denom = 0.0, step, bins
    members = []
    for i in range(1, bins + 1):
        members.append(np.nonzero((v >= lo) & (v <= hi))[0])
        lo = hi
        hi = hi + step
    edges = [(j * step, (j + 1) * step) for j in range(bins)]
    b_list = [i / denom for i in range(1, bins + 1)]
    return members, edges, b_list, v, []


def partition_equal(pos, dim, axial, bins, lat_p, axis, sample):
    """Corrected way: N equal bins, every atom in exactly one bin."""
    messages = []
    if axial:
        half = lat_p / 2
        layer = np.rint(pos[:, 2] / half).astype(int)  # atomic layer number 0..K-1
        K = int(layer.max()) + 1
        bin_id = (layer * bins) // K
        v = pos[:, 2]
        edges = [(j * K / bins * half, (j + 1) * K / bins * half) for j in range(bins)]
        if K < bins:
            messages.append("Warning: only %d atomic layers for %d bins; some bins are empty." % (K, bins))
    else:
        d = pos - axis
        v = np.sqrt((d ** 2).sum(axis=1)) if dim == 3 else np.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2)
        r_max = float(v.max())
        width = r_max / bins
        bin_id = np.minimum(np.floor(v / width).astype(int), bins - 1)
        edges = [(j * width, (j + 1) * width) for j in range(bins)]
    members = [np.nonzero(bin_id == j)[0] for j in range(bins)]
    if sample == "edge":
        b_list = [j / (bins - 1) for j in range(bins)]
    else:
        b_list = [(j + 0.5) / bins for j in range(bins)]
    return members, edges, b_list, v, messages


def recipe_at(g, start, end):
    """Mix of elements at grading value g (0..1): from the start recipe to the end recipe.
    A recipe is a list of fractions that adds to 1, so every mix is valid."""
    return [(1 - g) * s + g * e for s, e in zip(start, end)]


def counts_from_recipe(c, n, masses, basis="mass"):
    """How many of n atoms each element gets. Always adds up to exactly n.
    basis "mass": c are mass fractions (nanoHUB); basis "atom": c are atom fractions."""
    N = len(c)
    if basis == "mass":
        t = [ci / m for ci, m in zip(c, masses)]
        den = t[0]
        for ti in t[1:]:
            den = den + ti
        r = [ti / den for ti in t]
    else:
        r = list(c)
    active = [i for i in range(N) if r[i] > 0]
    counts = [0] * N
    if len(active) == 1:
        counts[active[0]] = n
    elif len(active) == 2:  # same rounding as nanoHUB: the later element is rounded, the other gets the rest
        i, j = active
        kj = min(max(round_half_up(r[j] * n), 0), n)
        counts[j] = kj
        counts[i] = n - kj
    else:  # three or more: largest-remainder rounding (ties go to the later element)
        q = [ri * n for ri in r]
        base = [int(math.floor(x)) for x in q]
        rem = min(max(n - sum(base), 0), N)
        order = sorted(range(N), key=lambda e: (q[e] - base[e], e), reverse=True)
        for e in order[:rem]:
            base[e] += 1
        counts = base
    return counts


def assign_cell(types, idx, counts, rng, overwrite):
    """Give atoms of one cell their element numbers (1, 2, 3, ...) at random.
    overwrite=True is the nanoHUB behaviour (two elements only): the picks ignore atoms
    that were already painted by a neighbouring bin."""
    n = len(idx)
    if overwrite:
        k = counts[1]
        if k > 0:
            types[idx[rng.choice(n, k, replace=False)]] = 2
        return
    for e in range(1, len(counts)):
        k = counts[e]
        if k <= 0:
            continue
        pool = idx if e == 1 else idx[types[idx] == 1]  # atoms of this cell not yet given an element
        types[pool[rng.choice(len(pool), k, replace=False)]] = e + 1


def engine_build(shape="cylindrical", radius_nm=2.0, length_nm=20.0, xdim_nm=2.0, ydim_nm=2.0,
                 lattice="fcc", elements=None, start=None, end=None, direction="radial",
                 grading="power", p=1.0, custom="", bins=10, partition="equal", sample="center",
                 basis="mass", vacuum_xy=0.0, vacuum_z=0.0, seed=None, comment=None,
                 lattice_param_nm=None, mode_label=None):
    """The one build function behind every mode.
    elements: list of dicts {"symbol": "Ni", "mass": 58.693, "a_nm": 0.3524}  (2 to 8 elements)
    start, end: recipes (fractions per element, each adds to 1); the grading function moves
                the mix from the start recipe to the end recipe.
    partition: "legacy" (nanoHUB copy, two elements only) or "equal" (corrected bins)
    lattice_param_nm: use this one lattice parameter instead of the mean of the elements' values."""
    if not elements or len(elements) < 2:
        raise ValueError("need at least 2 elements")
    n_el = len(elements)
    start = [float(x) for x in (start if start is not None else [1.0] + [0.0] * (n_el - 1))]
    end = [float(x) for x in (end if end is not None else [0.0] * (n_el - 1) + [1.0])]
    for name, rec in (("start", start), ("end", end)):
        if len(rec) != n_el or min(rec) < 0 or abs(sum(rec) - 1.0) > 1e-9:
            raise ValueError("%s recipe must have %d non-negative fractions that add to 1" % (name, n_el))
    legacy = partition == "legacy"
    if legacy and n_el != 2:
        raise ValueError("the nanoHUB-copy partition works with two elements only")
    masses = [float(e["mass"]) for e in elements]
    avals = [float(e["a_nm"]) for e in elements]
    rng = np.random.default_rng(seed)
    axial = direction == "axial"
    messages = []

    dim, leng, r_base, x_len, y_len, z_len = sizes_for(shape, radius_nm, length_nm, xdim_nm, ydim_nm)
    if lattice_param_nm is not None:
        lat_p = 10 * float(lattice_param_nm)
    else:
        lat_p = 10 * sum(avals) / n_el  # one averaged lattice constant, like nanoHUB

    if legacy:  # nanoHUB only: bins cannot be thinner than half a lattice constant
        limit = (2 * leng if axial else 2 * r_base) / lat_p
        if bins > limit:
            bins = round_half_up(limit)
            messages.append("Warning: number of bins too high, reset to %d." % bins)
    if bins < 2:
        raise ValueError("need at least 2 bins")

    geo = make_geometry(shape, radius_nm, length_nm, xdim_nm, ydim_nm, lattice, lat_p,
                        whole_cells_z=not legacy)
    pos, axis = geo["pos"], geo["axis"]
    messages += geo["messages"]
    n_atom = len(pos)
    types = np.ones(n_atom, dtype=int)

    if legacy:
        members, edges, b_list, v, pmsg = partition_legacy(pos, dim, axial, bins, r_base, x_len, y_len, z_len)
    else:
        members, edges, b_list, v, pmsg = partition_equal(pos, dim, axial, bins, lat_p, axis, sample)
    messages += pmsg

    bin_info = []
    yy = [0.0]  # "modelled" mass fraction of the LAST element per bin (first entry 0, as in nanoHUB)
    for j in range(bins):
        idx = members[j]
        n_in = len(idx)
        b = b_list[j]
        g = grading_fraction(b, grading, p, custom)
        raw = grading_raw(b, grading, p, custom) if legacy else None
        if isinstance(raw, complex):  # nanoHUB quirk (see grading_raw); reproduced exactly
            k2, ycurve = legacy_counts_complex(raw, n_in, masses[0], masses[1])
            counts = [n_in - k2, k2]
            yy.append(ycurve)
        else:
            counts = counts_from_recipe(recipe_at(g, start, end), n_in, masses, basis)
            den = 0.0
            for cnt, m in zip(counts, masses):
                den = den + cnt * m
            yy.append(counts[-1] * masses[-1] / den if den > 0 else float("nan"))
        assign_cell(types, idx, counts, rng, overwrite=legacy)
        bin_info.append({"bin": j + 1, "lo": edges[j][0], "hi": edges[j][1], "b": b, "x": g,
                         "n_atoms": n_in, "counts": counts, "n_type2_requested": n_in - counts[0]})
    if not legacy:
        empty = [bi["bin"] for bi in bin_info if bi["n_atoms"] == 0]
        if empty:
            messages.append("Warning: empty bins %s (too many bins for this size)." % empty)

    atoms = np.column_stack([np.arange(1, n_atom + 1), types, pos])

    def span(c, pad):
        lo_, hi_ = pos[:, c].min(), pos[:, c].max()
        return (lo_ - pad, hi_ + pad) if pad > 0 else (lo_, hi_ + lat_p / 2)

    box = {"x": span(0, vacuum_xy), "y": span(1, vacuum_xy), "z": span(2, vacuum_z)}
    return {"atoms": atoms, "box": box, "lat_p": lat_p, "bins": bins, "bin_info": bin_info,
            "yy": yy, "coord": v, "messages": messages, "mass1": masses[0], "mass2": masses[1],
            "axial": axial, "mode": mode_label or ("nanohub" if legacy else "fixed"),
            "members": members, "axis": axis, "n_cells_z": geo["n_cells_z"], "lattice": lattice,
            "elements": elements, "masses": masses,
            "comment": comment or "LAMMPS data file, FGM Nanostructure Builder (elements: %s)"
                                  % " ".join(str(e.get("symbol", i + 1)) for i, e in enumerate(elements))}


def build_structure(shape="cylindrical", radius_nm=2.0, length_nm=20.0, xdim_nm=2.0,
                    ydim_nm=2.0, direction="radial", grading="power", p=1.0, custom="",
                    mass1=63.55, mass2=58.69, lattice="fcc", a1_nm=0.36, a2_nm=0.36,
                    bins=10, seed=None, mode="nanohub", sample="center", fraction="mass",
                    vacuum_xy=0.0, vacuum_z=0.0):
    """The two-element builder used by the window and the command line.
    mode="nanohub": nanoHUB copy.  mode="fixed": corrected bins.
    Both run through engine_build; reference_build_structure (frozen) is only used by the self-test."""
    fixed = mode == "fixed"
    elements = [{"symbol": "1", "mass": mass1, "a_nm": a1_nm}, {"symbol": "2", "mass": mass2, "a_nm": a2_nm}]
    return engine_build(shape=shape, radius_nm=radius_nm, length_nm=length_nm, xdim_nm=xdim_nm,
                        ydim_nm=ydim_nm, lattice=lattice, elements=elements, start=[1.0, 0.0],
                        end=[0.0, 1.0], direction=direction, grading=grading, p=p, custom=custom,
                        bins=bins, partition="equal" if fixed else "legacy", sample=sample,
                        basis=fraction if fixed else "mass", vacuum_xy=vacuum_xy,
                        vacuum_z=vacuum_z, seed=seed, comment=FILE_COMMENT, mode_label=mode)


def vegard_constant(counts, avals_nm):
    """Lattice constant (nm) of a mix with these atom counts: the atom-weighted mean of the
    elements' constants (Vegard's law)."""
    tot = float(sum(counts))
    return sum(float(c) * a for c, a in zip(counts, avals_nm)) / tot if tot > 0 else None


def z_positions(plane_s):
    """Height of every atomic layer when layer k has half-constant plane_s[k]: the distance
    between two neighbouring layers is the mean of their two half-constants."""
    zp = np.zeros(len(plane_s))
    zp[1:] = np.cumsum(0.5 * (plane_s[:-1] + plane_s[1:]))
    return zp


# -------------------------------------------------------------------------------------
# layered mode: stacked layers along the wire, each with its own composition recipe(s),
# grading function and grading direction (along the wire or from core to surface)
# -------------------------------------------------------------------------------------
def layered_build(shape="cylindrical", radius_nm=2.0, xdim_nm=2.0, ydim_nm=2.0, lattice="fcc",
                  elements=None, layers=None, basis="atom", sample="center", vacuum_xy=0.0,
                  vacuum_z=0.0, seed=None, lattice_param_nm=None, comment=None, z_spacing="same"):
    """Build a wire made of stacked layers (2 to 8 elements).
    elements: list of {"symbol", "mass", "a_nm"}
    layers: list (bottom to top) of dicts
        thickness_nm : asked thickness (rounded to whole atomic layers)
        start, end   : recipes (fractions per element, each adds to 1); end=None means uniform
        grading      : "power" | "exp" | "sigmoid"  (moves the mix from start to end)
        p            : grading exponent
        direction    : "axial" (along the wire, inside the layer) or "radial" (core to surface)
        bins         : number of cells the layer is cut into (ignored for a uniform layer)
    z_spacing: "same" = every layer uses the one lattice parameter (as nanoHUB);
               "vegard" = the spacing along the wire follows each layer's own composition (the x,y spacing
               stays common, so the layers are strained to fit each other). A layer may also carry
               z_a_nm: its own lattice constant for the spacing along the wire (overrides the composition).
        custom       : with grading "custom", the formula f(x) (x = position 0..1 inside the layer; 0 gives pure
                       start recipe, 1 pure end recipe, a step such as x>0.7 gives solid regions)
    A graded interface is simply a thin graded layer between two others."""
    if shape not in ("cylindrical", "cuboid"):
        raise ValueError("layers are stacked along the wire: use a cylindrical or cuboid wire")
    if not elements or len(elements) < 2:
        raise ValueError("need at least 2 elements")
    if not layers:
        raise ValueError("need at least one layer")
    if z_spacing not in ("same", "vegard"):
        raise ValueError("z_spacing must be same or vegard")
    n_el = len(elements)
    masses = [float(e["mass"]) for e in elements]
    avals = [float(e["a_nm"]) for e in elements]
    lat_p = 10 * float(lattice_param_nm) if lattice_param_nm is not None else 10 * sum(avals) / n_el
    half = lat_p / 2.0
    rng = np.random.default_rng(seed)
    messages = []

    def recipe_of(vals, what):
        r = [float(x) for x in vals]
        if len(r) != n_el or min(r) < 0 or abs(sum(r) - 1.0) > 1e-9:
            raise ValueError("%s must have %d non-negative fractions that add to 1" % (what, n_el))
        return r

    L = []
    for i, ly in enumerate(layers, 1):
        start = recipe_of(ly["start"], "layer %d start recipe" % i)
        end = recipe_of(ly["end"], "layer %d end recipe" % i) if ly.get("end") is not None else list(start)
        uniform = max(abs(a - b) for a, b in zip(start, end)) < 1e-12
        if ly.get("z_a_nm") is not None and not float(ly["z_a_nm"]) > 0:
            raise ValueError("layer %d: lattice constant must be larger than 0" % i)
        if not float(ly["thickness_nm"]) > 0:
            raise ValueError("layer %d: thickness must be larger than 0" % i)
        bins = 1 if uniform else int(ly.get("bins", 10))
        if bins < 1:
            raise ValueError("layer %d: need at least 1 bin" % i)
        direction = ly.get("direction", "axial")
        if direction not in ("axial", "radial"):
            raise ValueError("layer %d: direction must be axial or radial" % i)
        if not uniform:
            check_grading(ly.get("grading", "power"), float(ly.get("p", 1.0)), ly.get("custom", ""), "layer %d" % i)
        L.append({"start": start, "end": end, "uniform": uniform, "bins": bins, "direction": direction,
                  "grading": ly.get("grading", "power"), "p": float(ly.get("p", 1.0)), "custom": str(ly.get("custom", "") or ""),
                  "asked_nm": float(ly["thickness_nm"]),
                  "z_a": (float(ly["z_a_nm"]) if ly.get("z_a_nm") is not None else None),
                  "planes": max(1, round_half_up(float(ly["thickness_nm"]) * 10 / half))})
    if sum(x["planes"] for x in L) % 2 == 1:
        L[-1]["planes"] += 1
        messages.append("Total length made a whole number of unit cells (one atomic layer added to the last layer) "
                        "so periodic z has no stacking seam.")
    K = sum(x["planes"] for x in L)
    if any(abs(x["planes"] * half / 10 - x["asked_nm"]) > 0.01 for x in L):
        messages.append("Layer thicknesses were rounded to whole atomic layers (see the layer table).")

    geo = make_geometry(shape, radius_nm, K * half / 10.0, xdim_nm, ydim_nm, lattice, lat_p, whole_cells_z=True)
    pos, axis = geo["pos"], geo["axis"]
    messages += geo["messages"]
    plane = np.rint(pos[:, 2] / half).astype(int)
    if int(plane.max()) + 1 != K:
        raise RuntimeError("internal error: the crystal has %d atomic layers, expected %d" % (int(plane.max()) + 1, K))
    n_atom = len(pos)
    types = np.ones(n_atom, dtype=int)
    d = pos - axis
    rad = np.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2)
    r_max = float(rad.max())

    members, bin_info, yy = [], [], [0.0]
    P0 = 0
    for li, ly in enumerate(L, 1):
        P1 = P0 + ly["planes"]
        in_layer = np.nonzero((plane >= P0) & (plane < P1))[0]
        nb = ly["bins"]
        if ly["direction"] == "axial" or ly["uniform"]:
            cid = ((plane[in_layer] - P0) * nb) // ly["planes"]
            edges = [((P0 + j * ly["planes"] / nb) * half, (P0 + (j + 1) * ly["planes"] / nb) * half) for j in range(nb)]
        else:
            width = r_max / nb
            cid = np.minimum(np.floor(rad[in_layer] / width).astype(int), nb - 1)
            edges = [(j * width, (j + 1) * width) for j in range(nb)]
        for j in range(nb):
            idx = in_layer[cid == j]
            if ly["uniform"]:
                b, g, recipe = 0.5, 0.0, list(ly["start"])
            else:
                b = (j + 0.5) / nb if sample != "edge" else (j / (nb - 1) if nb > 1 else 0.5)
                g = grading_fraction(b, ly["grading"], ly["p"], ly["custom"])
                recipe = recipe_at(g, ly["start"], ly["end"])
            counts = counts_from_recipe(recipe, len(idx), masses, basis)
            den = 0.0
            for cnt, m in zip(counts, masses):
                den = den + cnt * m
            yy.append(counts[-1] * masses[-1] / den if den > 0 else float("nan"))
            assign_cell(types, idx, counts, rng, overwrite=False)
            members.append(idx)
            bin_info.append({"bin": len(bin_info) + 1, "layer": li, "direction": ly["direction"],
                             "plane_lo": P0 + j * ly["planes"] / nb, "plane_hi": P0 + (j + 1) * ly["planes"] / nb,
                             "lo": edges[j][0], "hi": edges[j][1], "b": b, "x": g, "n_atoms": len(idx),
                             "counts": counts, "n_type2_requested": len(idx) - counts[0]})
        ly["z0"], ly["z1"] = P0 * half, P1 * half
        ly["atoms"] = int(len(in_layer))
        P0 = P1
    empty = [bi["bin"] for bi in bin_info if bi["n_atoms"] == 0]
    if empty:
        messages.append("Warning: empty cells %s (too many bins for this layer)." % empty)
    small = sum(1 for bi in bin_info if 0 < bi["n_atoms"] < 20)
    if small:
        messages.append("Note: %d cells hold fewer than 20 atoms, so their composition is rounded coarsely." % small)

    gap = half
    pos_unstrained = None
    for ly in L:
        ly["a_z"] = lat_p / 10
    if z_spacing == "vegard":
        pos_unstrained = pos.copy()
        plane_s = np.full(K, half)
        P0 = 0
        for li, ly in enumerate(L, 1):
            P1 = P0 + ly["planes"]
            cells = [bi for bi in bin_info if bi["layer"] == li]
            tot = np.sum([bi["counts"] for bi in cells], axis=0)
            a_layer = vegard_constant(tot, avals)
            if ly["z_a"] is not None:
                plane_s[P0:P1] = 10 * ly["z_a"] / 2
            elif ly["direction"] == "axial" or ly["uniform"]:
                local = (np.arange(ly["planes"]) * ly["bins"]) // ly["planes"]
                for j, bi in enumerate(cells):
                    a_cell = vegard_constant(bi["counts"], avals) if bi["n_atoms"] > 0 else a_layer
                    plane_s[P0:P1][local == j] = 10 * (a_cell if a_cell is not None else lat_p / 10) / 2
            else:
                plane_s[P0:P1] = 10 * (a_layer if a_layer is not None else lat_p / 10) / 2
            ly["a_z"] = float(2 * plane_s[P0:P1].mean() / 10)
            P0 = P1
        zp = z_positions(plane_s)
        pos[:, 2] = zp[plane]
        gap = 0.5 * (plane_s[-1] + plane_s[0])
        zp_ext = np.append(zp, zp[-1] + gap)
        P0 = 0
        for ly in L:
            P1 = P0 + ly["planes"]
            ly["z0"], ly["z1"] = float(zp_ext[P0]), float(zp_ext[P1])
            P0 = P1
        for bi in bin_info:
            if bi["direction"] == "axial" or L[bi["layer"] - 1]["uniform"]:
                bi["lo"] = float(np.interp(bi["plane_lo"], np.arange(K + 1), zp_ext))
                bi["hi"] = float(np.interp(bi["plane_hi"], np.arange(K + 1), zp_ext))

    atoms = np.column_stack([np.arange(1, n_atom + 1), types, pos])

    def span(c, pad):
        lo_, hi_ = pos[:, c].min(), pos[:, c].max()
        return (lo_ - pad, hi_ + pad) if pad > 0 else (lo_, hi_ + (gap if c == 2 else lat_p / 2))

    box = {"x": span(0, vacuum_xy), "y": span(1, vacuum_xy), "z": span(2, vacuum_z)}
    layer_info = [{"layer": i, "planes": x["planes"], "asked_nm": x["asked_nm"], "built_nm": (x["z1"] - x["z0"]) / 10, "a_z_nm": x["a_z"],
                   "atoms": x["atoms"], "start": x["start"], "end": x["end"], "uniform": x["uniform"],
                   "grading": x["grading"], "p": x["p"], "custom": x["custom"], "direction": x["direction"], "bins": x["bins"],
                   "z0": x["z0"], "z1": x["z1"]} for i, x in enumerate(L, 1)]
    return {"atoms": atoms, "box": box, "lat_p": lat_p, "bins": len(bin_info), "bin_info": bin_info,
            "yy": yy, "coord": pos[:, 2], "messages": messages, "mass1": masses[0], "mass2": masses[1],
            "axial": all(x["direction"] == "axial" for x in L), "mode": "layered", "members": members,
            "axis": axis, "n_cells_z": K // 2, "lattice": lattice, "elements": elements, "masses": masses,
            "layers": layer_info, "basis": basis, "plane": plane, "z_spacing": z_spacing, "pos_unstrained": pos_unstrained,
            "comment": comment or "LAMMPS data file, FGM Nanostructure Builder, layered wire (elements: %s)"
                                  % " ".join(str(e.get("symbol", i + 1)) for i, e in enumerate(elements))}


# -------------------------------------------------------------------------------------
# combined mode: composition changes along the wire AND from the core to the surface
#   Four corner recipes (core/surface x bottom/top). Two grading functions decide how fast the
#   mix moves from one corner to the next, one along the wire and one across it. A blend of valid
#   recipes is always a valid recipe, so the composition can never leave 0..100 %.
# -------------------------------------------------------------------------------------
def mix_recipes(g, a, b):
    """Recipe at grading value g between recipes a and b (exactly a when they are equal)."""
    return list(a) if list(a) == list(b) else recipe_at(g, a, b)


def combined_build(shape="cylindrical", radius_nm=2.0, length_nm=10.0, xdim_nm=2.0, ydim_nm=2.0,
                   lattice="fcc", elements=None, corners=None, grading_r="power", p_r=1.0, bins_r=5,
                   grading_z="power", p_z=1.0, bins_z=10, basis="atom", sample="center",
                   vacuum_xy=0.0, vacuum_z=0.0, seed=None, lattice_param_nm=None, comment=None, z_spacing="same",
                   custom_r="", custom_z=""):
    """Build a wire whose composition depends on both the position along the wire (z) and the
    distance from the axis (r).
    corners: four recipes [core-bottom, surface-bottom, core-top, surface-top], each fractions per
             element adding to 1.
    z_spacing: "same" or "vegard" (spacing along the wire follows the composition of each slice)
    custom_r/custom_z : formula f(x) when grading_r/grading_z is "custom" (0 = pure start side, 1 = pure end side;
                        a step such as x>0.7 gives solid regions)
    grading_r/p_r/bins_r : function, exponent and number of rings (core to surface)
    grading_z/p_z/bins_z : function, exponent and number of slices (bottom to top)"""
    if shape not in ("cylindrical", "cuboid"):
        raise ValueError("combined grading needs a cylindrical or cuboid wire")
    if not elements or len(elements) < 2:
        raise ValueError("need at least 2 elements")
    n_el = len(elements)
    if not corners or len(corners) != 4:
        raise ValueError("need four corner recipes: core-bottom, surface-bottom, core-top, surface-top")
    cor = []
    for name, rec in zip(("core-bottom", "surface-bottom", "core-top", "surface-top"), corners):
        r = [float(x) for x in rec]
        if len(r) != n_el or min(r) < 0 or abs(sum(r) - 1.0) > 1e-9:
            raise ValueError("%s recipe must have %d non-negative fractions that add to 1" % (name, n_el))
        cor.append(r)
    c00, c10, c01, c11 = cor
    if z_spacing not in ("same", "vegard"):
        raise ValueError("z_spacing must be same or vegard")
    check_grading(grading_r, p_r, custom_r, "core to surface")
    check_grading(grading_z, p_z, custom_z, "along the wire")
    bins_r, bins_z = int(bins_r), int(bins_z)
    if bins_r < 1 or bins_z < 1:
        raise ValueError("need at least 1 ring and 1 slice")
    masses = [float(e["mass"]) for e in elements]
    avals = [float(e["a_nm"]) for e in elements]
    lat_p = 10 * float(lattice_param_nm) if lattice_param_nm is not None else 10 * sum(avals) / n_el
    half = lat_p / 2.0
    rng = np.random.default_rng(seed)
    messages = []

    geo = make_geometry(shape, radius_nm, length_nm, xdim_nm, ydim_nm, lattice, lat_p, whole_cells_z=True)
    pos, axis = geo["pos"], geo["axis"]
    messages += geo["messages"]
    plane = np.rint(pos[:, 2] / half).astype(int)
    K = int(plane.max()) + 1
    n_atom = len(pos)
    types = np.ones(n_atom, dtype=int)
    d = pos - axis
    rad = np.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2)
    r_max = float(rad.max())
    width = r_max / bins_r
    iz = (plane * bins_z) // K
    ir = np.minimum(np.floor(rad / width).astype(int), bins_r - 1)
    cell = iz * bins_r + ir
    order = np.argsort(cell, kind="stable")  # atoms of one cell stay in file order
    starts = np.searchsorted(cell[order], np.arange(bins_z * bins_r + 1))

    def b_of(j, nb):
        return (j + 0.5) / nb if sample != "edge" else (j / (nb - 1) if nb > 1 else 0.5)

    members, bin_info, yy = [], [], [0.0]
    for jz in range(bins_z):
        bz = b_of(jz, bins_z)
        gz = grading_fraction(bz, grading_z, p_z, custom_z)
        core = mix_recipes(gz, c00, c01)
        surf = mix_recipes(gz, c10, c11)
        for jr in range(bins_r):
            br = b_of(jr, bins_r)
            gr = grading_fraction(br, grading_r, p_r, custom_r)
            recipe = mix_recipes(gr, core, surf)
            c = jz * bins_r + jr
            idx = order[starts[c]:starts[c + 1]]
            counts = counts_from_recipe(recipe, len(idx), masses, basis)
            den = 0.0
            for cnt, m in zip(counts, masses):
                den = den + cnt * m
            yy.append(counts[-1] * masses[-1] / den if den > 0 else float("nan"))
            assign_cell(types, idx, counts, rng, overwrite=False)
            members.append(idx)
            bin_info.append({"bin": len(bin_info) + 1, "iz": jz + 1, "ir": jr + 1,
                             "lo": jz * K / bins_z * half, "hi": (jz + 1) * K / bins_z * half,
                             "r_lo": jr * width, "r_hi": (jr + 1) * width, "b": bz, "br": br, "gz": gz, "gr": gr, "x": gr,
                             "n_atoms": len(idx), "counts": counts, "recipe": recipe,
                             "n_type2_requested": len(idx) - counts[0]})
    empty = [bi["bin"] for bi in bin_info if bi["n_atoms"] == 0]
    if empty:
        messages.append("Warning: empty cells %s (too many rings or slices for this wire)." % empty)
    small = sum(1 for bi in bin_info if 0 < bi["n_atoms"] < 20)
    if small:
        messages.append("Note: %d cells hold fewer than 20 atoms, so their composition is rounded coarsely." % small)

    gap = half
    pos_unstrained = None
    for bi in bin_info:
        bi["a_z_nm"] = lat_p / 10
    if z_spacing == "vegard":
        pos_unstrained = pos.copy()
        plane_s = np.full(K, half)
        slice_of_plane = (np.arange(K) * bins_z) // K
        a_all = vegard_constant(np.sum([bi["counts"] for bi in bin_info], axis=0), avals)
        for jz in range(bins_z):
            cells = [bi for bi in bin_info if bi["iz"] == jz + 1]
            a_sl = vegard_constant(np.sum([bi["counts"] for bi in cells], axis=0), avals)
            a_sl = a_all if a_sl is None else a_sl
            plane_s[slice_of_plane == jz] = 10 * a_sl / 2
            for bi in cells:
                bi["a_z_nm"] = float(a_sl)
        zp = z_positions(plane_s)
        pos[:, 2] = zp[plane]
        gap = 0.5 * (plane_s[-1] + plane_s[0])
        zp_ext = np.append(zp, zp[-1] + gap)
        for bi in bin_info:
            jz = bi["iz"] - 1
            bi["lo"] = float(np.interp(jz * K / bins_z, np.arange(K + 1), zp_ext))
            bi["hi"] = float(np.interp((jz + 1) * K / bins_z, np.arange(K + 1), zp_ext))

    atoms = np.column_stack([np.arange(1, n_atom + 1), types, pos])

    def span(c_, pad):
        lo_, hi_ = pos[:, c_].min(), pos[:, c_].max()
        return (lo_ - pad, hi_ + pad) if pad > 0 else (lo_, hi_ + (gap if c_ == 2 else lat_p / 2))

    box = {"x": span(0, vacuum_xy), "y": span(1, vacuum_xy), "z": span(2, vacuum_z)}
    return {"atoms": atoms, "box": box, "lat_p": lat_p, "bins": len(bin_info), "bin_info": bin_info,
            "yy": yy, "coord": pos[:, 2], "messages": messages, "mass1": masses[0], "mass2": masses[1],
            "axial": False, "mode": "combined", "members": members, "axis": axis, "n_cells_z": K // 2,
            "lattice": lattice, "elements": elements, "masses": masses, "basis": basis,
            "plane": plane, "z_spacing": z_spacing, "pos_unstrained": pos_unstrained,
            "combined": {"corners": cor, "grading_r": grading_r, "p_r": p_r, "bins_r": bins_r,
                         "grading_z": grading_z, "p_z": p_z, "bins_z": bins_z, "planes": K, "r_max": r_max,
                         "custom_r": custom_r, "custom_z": custom_z},
            "comment": comment or "LAMMPS data file, FGM Nanostructure Builder, combined axial-radial grading (elements: %s)"
                                  % " ".join(str(e.get("symbol", i + 1)) for i, e in enumerate(elements))}


def combined_summary_lines(res):
    """Cell table (slice x ring) of a combined-grading wire."""
    els = res["elements"]
    names = [str(e.get("symbol", i + 1)) for i, e in enumerate(els)]
    cb = res["combined"]
    unit = "atom%" if res.get("basis") == "atom" else "mass%"
    fmt = lambda r: "/".join("%g" % round(100 * x, 2) for x in r)
    lines = list(res["messages"])
    lines.append("mode: combined   atoms: %d   box: %s" % (len(res["atoms"]),
                 {k: tuple(round(float(x), 4) for x in v) for k, v in res["box"].items()}))
    lines += ["", "CORNER RECIPES in %s of: %s" % (unit, " / ".join(names)),
              "  core,    bottom: %s" % fmt(cb["corners"][0]), "  surface, bottom: %s" % fmt(cb["corners"][1]),
              "  core,    top   : %s" % fmt(cb["corners"][2]), "  surface, top   : %s" % fmt(cb["corners"][3]),
              "  along the wire : %s, %d slices      core to surface: %s, %d rings" % (
                  fn_label(cb["grading_z"], cb["p_z"], cb.get("custom_z", "")), cb["bins_z"],
                  fn_label(cb["grading_r"], cb["p_r"], cb.get("custom_r", "")), cb["bins_r"]),
              "", "CELLS (atoms of each element)",
              "%4s %5s %4s %19s %17s %7s " % ("cell", "slice", "ring", "z range (A)", "r range (A)", "atoms")
              + " ".join("%7s" % nme for nme in names) + "   a_z nm"]
    for bi in res["bin_info"]:
        lines.append("%4d %5d %4d %8.1f..%8.1f %7.1f..%7.1f %7d " % (
            bi["bin"], bi["iz"], bi["ir"], bi["lo"], bi["hi"], bi["r_lo"], bi["r_hi"], bi["n_atoms"])
            + " ".join("%7d" % c for c in bi["counts"]))
    for name, ok in verify(res):
        lines.append("[%s] %s" % ("PASS" if ok else "FAIL", name))
    return lines


# -------------------------------------------------------------------------------------
# crystal layers: each layer is its own close-packed crystal, stacked along the wire
#   FCC with its (111) planes across the wire, or HCP with its (0001) planes across the wire.
#   Both have the same triangular atom pattern in every plane, so layers of different crystals
#   (for example Ti hexagonal and Al cubic) can be joined with matching close-packed rows.
#   The sideways spacing is shared (layers are strained to fit); each layer keeps its own spacing
#   between atomic planes.
# -------------------------------------------------------------------------------------
CRYSTALS = {"fcc111": "FCC (111)", "hcp0001": "HCP (0001)"}
_PLANE_OFFSETS = [(0.0, 0.0), (0.5, 0.5 / math.sqrt(3.0)), (0.0, 1.0 / math.sqrt(3.0))]  # positions A, B, C in units of d


def crystal_params(crystal, a_nm, c_nm=None):
    """(nearest-neighbour distance in the plane, spacing between planes), both in Angstrom."""
    if crystal not in CRYSTALS:
        raise ValueError("crystal must be fcc111 or hcp0001")
    a = 10.0 * float(a_nm)
    if not a > 0:
        raise ValueError("lattice constant a must be larger than 0")
    if crystal == "fcc111":
        return a / math.sqrt(2.0), a / math.sqrt(3.0)
    c = 10.0 * float(c_nm) if c_nm not in (None, "") else a * math.sqrt(8.0 / 3.0)
    if not c > 0:
        raise ValueError("lattice constant c must be larger than 0")
    return a, c / 2.0


def stacking_letters(layers_planes, kinds):
    """Stacking letter (0=A, 1=B, 2=C) of every atomic plane. A layer starts on a different letter than
    the plane below it; FCC runs A B C A B C..., HCP alternates two letters."""
    def seq(first, kind, n):
        return [(first + k) % 3 if kind == "fcc111" else (first + (k % 2)) % 3 for k in range(n)]
    letters, prev = [], None
    for li, (n, kind) in enumerate(zip(layers_planes, kinds)):
        options = [0] if prev is None else [(prev + 1) % 3, (prev + 2) % 3]
        pick = options[0]
        if li == len(layers_planes) - 1 and letters:  # top layer: also avoid the letter of the very first plane
            for o in options:
                if seq(o, kind, n)[-1] != letters[0]:
                    pick = o
                    break
        s_ = seq(pick, kind, n)
        letters += s_
        prev = s_[-1]
    return letters


def crystal_layers_build(shape="cylindrical", radius_nm=1.5, xdim_nm=2.0, ydim_nm=2.0, elements=None, layers=None,
                         inplane="mean", inplane_nm=None, basis="atom", sample="center", vacuum_xy=0.0,
                         vacuum_z=0.0, seed=None, comment=None):
    """Build a wire from stacked close-packed crystal layers (FCC (111) and HCP (0001) can be mixed).
    elements: list of {"symbol", "mass"} (lattice constants are given per layer)
    layers (bottom to top), each a dict:
        crystal      : "fcc111" or "hcp0001"
        a_nm, c_nm   : lattice constants of that crystal (c_nm only for HCP; blank = ideal ratio)
        thickness_nm : asked thickness (rounded to whole atomic planes)
        start, end   : composition recipes (end=None: uniform); grading/p/direction/bins as in layered_build
    inplane: "mean" (mean of the layers' own in-plane spacings), "first" (match the bottom layer) or
             "value" (use inplane_nm as the nearest-neighbour distance)"""
    if shape not in ("cylindrical", "cuboid"):
        raise ValueError("crystal layers need a cylindrical or cuboid wire")
    if not elements or len(elements) < 2:
        raise ValueError("need at least 2 elements")
    if not layers:
        raise ValueError("need at least one layer")
    n_el = len(elements)
    masses = [float(e["mass"]) for e in elements]
    rng = np.random.default_rng(seed)
    messages = []

    def recipe_of(vals, what):
        r = [float(x) for x in vals]
        if len(r) != n_el or min(r) < 0 or abs(sum(r) - 1.0) > 1e-9:
            raise ValueError("%s must have %d non-negative fractions that add to 1" % (what, n_el))
        return r

    L = []
    for i, ly in enumerate(layers, 1):
        d_own, h = crystal_params(ly.get("crystal"), ly.get("a_nm"), ly.get("c_nm"))
        if not float(ly["thickness_nm"]) > 0:
            raise ValueError("layer %d: thickness must be larger than 0" % i)
        start = recipe_of(ly["start"], "layer %d start recipe" % i)
        end = recipe_of(ly["end"], "layer %d end recipe" % i) if ly.get("end") is not None else list(start)
        uniform = max(abs(a - b) for a, b in zip(start, end)) < 1e-12
        direction = ly.get("direction", "axial")
        if direction not in ("axial", "radial"):
            raise ValueError("layer %d: direction must be axial or radial" % i)
        bins = 1 if uniform else int(ly.get("bins", 10))
        if bins < 1:
            raise ValueError("layer %d: need at least 1 bin" % i)
        if not uniform:
            check_grading(ly.get("grading", "power"), float(ly.get("p", 1.0)), ly.get("custom", ""), "layer %d" % i)
        L.append({"crystal": ly["crystal"], "a_nm": float(ly["a_nm"]), "c_nm": (float(ly["c_nm"]) if ly.get("c_nm") not in (None, "") else None),
                  "d_own": d_own, "h": h, "start": start, "end": end, "uniform": uniform, "bins": bins, "direction": direction,
                  "grading": ly.get("grading", "power"), "p": float(ly.get("p", 1.0)), "custom": str(ly.get("custom", "") or ""), "asked_nm": float(ly["thickness_nm"]),
                  "planes": max(1, round_half_up(float(ly["thickness_nm"]) * 10 / h))})
    if inplane == "first":
        d_xy = L[0]["d_own"]
    elif inplane == "value":
        if inplane_nm is None or not float(inplane_nm) > 0:
            raise ValueError("give the in-plane nearest-neighbour distance (nm)")
        d_xy = 10.0 * float(inplane_nm)
    elif inplane == "mean":
        d_xy = sum(x["d_own"] for x in L) / len(L)
    else:
        raise ValueError("inplane must be mean, first or value")

    # stacking letters; if the top plane would sit directly over the bottom plane (periodic z), add planes
    added = 0
    for attempt in range(4):
        letters = stacking_letters([x["planes"] for x in L], [x["crystal"] for x in L])
        closed = letters[-1] != letters[0]
        if len(L) == 1:  # one crystal only: its stacking must carry on without a fault across the periodic boundary
            closed = (L[0]["planes"] % 3 == 0) if L[0]["crystal"] == "fcc111" else (L[0]["planes"] % 2 == 0)
        if closed:
            break
        L[-1]["planes"] += 1
        added += 1
    if added:
        messages.append("%d atomic plane%s added to the top layer so that, with periodic z, the stacking %s." % (
            added, "s" if added > 1 else "",
            "carries on without a fault across the boundary" if len(L) == 1 else "does not put the top plane directly over the bottom plane"))
    K = sum(x["planes"] for x in L)
    plane_layer = np.concatenate([np.full(x["planes"], i) for i, x in enumerate(L)])
    h_plane = np.array([L[i]["h"] for i in plane_layer])
    zp = z_positions(h_plane)
    gap = 0.5 * (h_plane[-1] + h_plane[0])
    zp_ext = np.append(zp, zp[-1] + gap)

    # atoms of every plane: triangular lattice, shifted by the plane's stacking letter, cut to the wire shape
    if shape == "cylindrical":
        rx = ry = 10.0 * float(radius_nm)
    else:
        rx, ry = 10.0 * float(xdim_nm) / 2, 10.0 * float(ydim_nm) / 2
    rmax_box = math.hypot(rx, ry)
    Mj = int(math.ceil(rmax_box / (d_xy * math.sqrt(3.0) / 2))) + 2
    Mi = int(math.ceil(rmax_box / d_xy)) + Mj // 2 + 3
    jj, ii = np.meshgrid(np.arange(-Mj, Mj + 1), np.arange(-Mi, Mi + 1), indexing="ij")
    x0 = (ii * d_xy + jj * d_xy / 2.0).ravel()
    y0 = (jj * d_xy * math.sqrt(3.0) / 2.0).ravel()
    xs, ys, zs, pl = [], [], [], []
    for k in range(K):
        ox, oy = _PLANE_OFFSETS[letters[k]]
        x, y = x0 + ox * d_xy, y0 + oy * d_xy
        keep = (x * x + y * y <= rx * rx + 1e-9) if shape == "cylindrical" else ((np.abs(x) <= rx + 1e-9) & (np.abs(y) <= ry + 1e-9))
        xs.append(x[keep]); ys.append(y[keep]); zs.append(np.full(int(keep.sum()), zp[k])); pl.append(np.full(int(keep.sum()), k))
    X, Y, Z, plane = np.concatenate(xs), np.concatenate(ys), np.concatenate(zs), np.concatenate(pl)
    minx, miny = X.min(), Y.min()
    pos = np.column_stack([X - minx, Y - miny, Z])
    axis = np.array([-minx, -miny, 0.0])
    n_atom = len(pos)
    types = np.ones(n_atom, dtype=int)
    d = pos[:, :2] - axis[:2]
    rad = np.sqrt((d ** 2).sum(axis=1))
    r_max = float(rad.max())

    members, bin_info, yy = [], [], [0.0]
    P0 = 0
    for li, ly in enumerate(L, 1):
        P1 = P0 + ly["planes"]
        in_layer = np.nonzero((plane >= P0) & (plane < P1))[0]
        nb = ly["bins"]
        if ly["direction"] == "axial" or ly["uniform"]:
            cid = ((plane[in_layer] - P0) * nb) // ly["planes"]
            edges = [(float(np.interp(P0 + j * ly["planes"] / nb, np.arange(K + 1), zp_ext)),
                      float(np.interp(P0 + (j + 1) * ly["planes"] / nb, np.arange(K + 1), zp_ext))) for j in range(nb)]
        else:
            width = r_max / nb
            cid = np.minimum(np.floor(rad[in_layer] / width).astype(int), nb - 1)
            edges = [(j * width, (j + 1) * width) for j in range(nb)]
        for j in range(nb):
            idx = in_layer[cid == j]
            if ly["uniform"]:
                b, g, recipe = 0.5, 0.0, list(ly["start"])
            else:
                b = (j + 0.5) / nb if sample != "edge" else (j / (nb - 1) if nb > 1 else 0.5)
                g = grading_fraction(b, ly["grading"], ly["p"], ly["custom"])
                recipe = recipe_at(g, ly["start"], ly["end"])
            counts = counts_from_recipe(recipe, len(idx), masses, basis)
            den = 0.0
            for cnt, m in zip(counts, masses):
                den = den + cnt * m
            yy.append(counts[-1] * masses[-1] / den if den > 0 else float("nan"))
            assign_cell(types, idx, counts, rng, overwrite=False)
            members.append(idx)
            bin_info.append({"bin": len(bin_info) + 1, "layer": li, "direction": ly["direction"], "lo": edges[j][0], "hi": edges[j][1],
                             "b": b, "x": g, "n_atoms": len(idx), "counts": counts, "n_type2_requested": len(idx) - counts[0]})
        ly["z0"], ly["z1"] = float(zp_ext[P0]), float(zp_ext[P1])
        ly["atoms"] = int(len(in_layer))
        P0 = P1
    empty = [bi["bin"] for bi in bin_info if bi["n_atoms"] == 0]
    if empty:
        messages.append("Warning: empty cells %s (too many bins for this layer)." % empty)
    small = sum(1 for bi in bin_info if 0 < bi["n_atoms"] < 20)
    if small:
        messages.append("Note: %d cells hold fewer than 20 atoms, so their composition is rounded coarsely." % small)
    strains = [100.0 * (d_xy / x["d_own"] - 1.0) for x in L]
    if max(abs(x) for x in strains) > 5.0:
        messages.append("Warning: a layer is strained by more than 5 %% in the plane (%.1f %%); treat the structure as a starting point only."
                        % max(strains, key=abs))

    atoms = np.column_stack([np.arange(1, n_atom + 1), types, pos])

    def span(c_, pad):
        lo_, hi_ = pos[:, c_].min(), pos[:, c_].max()
        return (lo_ - pad, hi_ + pad) if pad > 0 else (lo_, hi_ + (gap if c_ == 2 else d_xy / 2))

    box = {"x": span(0, vacuum_xy), "y": span(1, vacuum_xy), "z": span(2, vacuum_z)}
    layer_info = [{"layer": i, "crystal": x["crystal"], "a_nm": x["a_nm"], "c_nm": x["c_nm"], "planes": x["planes"],
                   "asked_nm": x["asked_nm"], "built_nm": (x["z1"] - x["z0"]) / 10, "a_z_nm": 2 * x["h"] / 10 if x["crystal"] == "hcp0001" else x["h"] * math.sqrt(3.0) / 10,
                   "strain_pct": strains[i - 1], "d_own": x["d_own"], "h": x["h"], "atoms": x["atoms"], "start": x["start"], "end": x["end"],
                   "uniform": x["uniform"], "grading": x["grading"], "p": x["p"], "custom": x["custom"], "direction": x["direction"], "bins": x["bins"],
                   "z0": x["z0"], "z1": x["z1"]} for i, x in enumerate(L, 1)]
    return {"atoms": atoms, "box": box, "lat_p": d_xy * math.sqrt(2.0), "bins": len(bin_info), "bin_info": bin_info, "yy": yy,
            "coord": pos[:, 2], "messages": messages, "mass1": masses[0], "mass2": masses[1],
            "axial": all(x["direction"] == "axial" for x in L), "mode": "crystals", "members": members, "axis": axis,
            "n_cells_z": None, "lattice": "mixed", "elements": elements, "masses": masses, "layers": layer_info, "basis": basis,
            "plane": plane, "letters": letters, "d_xy": d_xy, "z_spacing": "crystal", "pos_unstrained": None, "crystal_layers": True,
            "z_period": float(zp_ext[-1]),
            "comment": comment or "LAMMPS data file, FGM Nanostructure Builder, crystal layers (elements: %s)"
                                  % " ".join(str(e.get("symbol", i + 1)) for i, e in enumerate(elements))}


def crystal_summary_lines(res):
    """Layer table of a crystal-layer wire (crystal, strain, spacing) and its cells."""
    els = res["elements"]
    names = [str(e.get("symbol", i + 1)) for i, e in enumerate(els)]
    unit = "atom%" if res.get("basis") == "atom" else "mass%"
    fmt = lambda r: "/".join("%g" % round(100 * x, 2) for x in r)
    lines = list(res["messages"])
    lines.append("mode: crystal layers   atoms: %d   box: %s" % (len(res["atoms"]), {k: tuple(round(float(x), 4) for x in v) for k, v in res["box"].items()}))
    lines += ["", "In-plane nearest-neighbour distance shared by all layers: %.4f A" % res["d_xy"],
              "LAYERS (bottom to top), composition in %s of: %s" % (unit, " / ".join(names)),
              "%5s %-10s %7s %9s %9s %9s %8s %8s  %s" % ("layer", "crystal", "planes", "asked nm", "built nm", "strain %", "atoms", "grading", "composition")]
    for ly in res["layers"]:
        comp = fmt(ly["start"]) if ly["uniform"] else "%s -> %s (%s %s)" % (fmt(ly["start"]), fmt(ly["end"]), fn_label(ly["grading"], ly["p"], ly.get("custom", "")), ly["direction"])
        lines.append("%5d %-10s %7d %9.2f %9.2f %+9.2f %8d %8s  %s" % (ly["layer"], CRYSTALS[ly["crystal"]], ly["planes"], ly["asked_nm"], ly["built_nm"],
                                                                      ly["strain_pct"], ly["atoms"], "uniform" if ly["uniform"] else "graded", comp))
    lines += ["", "CELLS (atoms of each element)", "%4s %5s %20s %8s " % ("cell", "layer", "range (A)", "atoms") + " ".join("%8s" % nme for nme in names)]
    for bi in res["bin_info"]:
        lines.append("%4d %5d %9.2f..%9.2f %8d " % (bi["bin"], bi["layer"], bi["lo"], bi["hi"], bi["n_atoms"]) + " ".join("%8d" % c for c in bi["counts"]))
    for name, ok in verify(res):
        lines.append("[%s] %s" % ("PASS" if ok else "FAIL", name))
    return lines


# -------------------------------------------------------------------------------------
# self-test: proves the engine reproduces the frozen reference and behaves for 3+ elements
# -------------------------------------------------------------------------------------
def _same_result(a, b):
    """Are two build results identical (every atom, type, bin and message)?"""
    keys = ("lat_p", "bins", "messages", "n_cells_z", "axial", "lattice")
    if any(a[k] != b[k] for k in keys):
        return False
    if not np.array_equal(a["atoms"], b["atoms"]) or not np.array_equal(a["coord"], b["coord"]):
        return False
    if any(a["box"][k] != b["box"][k] for k in "xyz") or not np.array_equal(a["axis"], b["axis"]):
        return False
    if len(a["members"]) != len(b["members"]) or not all(np.array_equal(x, y) for x, y in zip(a["members"], b["members"])):
        return False
    for x, y in zip(a["yy"], b["yy"]):
        if not (x == y or (x != x and y != y)):
            return False
    for x, y in zip(a["bin_info"], b["bin_info"]):
        if any(x[k] != y[k] for k in ("bin", "lo", "hi", "b", "x", "n_atoms", "n_type2_requested")):
            return False
    return True


def _random_classic(rnd):
    return dict(shape=rnd.choice(list(SHAPES)), radius_nm=round(rnd.uniform(0.8, 2.0), 2),
                length_nm=round(rnd.uniform(3, 8), 1), xdim_nm=round(rnd.uniform(1.0, 3.0), 2),
                ydim_nm=round(rnd.uniform(1.0, 3.0), 2), direction=rnd.choice(["radial", "axial"]),
                grading=rnd.choice(["power", "exp", "sigmoid", "custom"]), p=round(rnd.uniform(0.3, 3.5), 2),
                custom=rnd.choice(["x^2", "1-x", "sin(x)", "sqrt(x)", "0.5*x"]),
                mass1=round(rnd.uniform(40, 70), 2), mass2=round(rnd.uniform(40, 70), 2),
                lattice=rnd.choice(["fcc", "bcc"]), a1_nm=round(rnd.uniform(0.30, 0.40), 4),
                a2_nm=round(rnd.uniform(0.30, 0.40), 4), bins=rnd.randint(2, 14),
                mode=rnd.choice(["nanohub", "fixed"]), sample=rnd.choice(["center", "edge"]),
                fraction=rnd.choice(["mass", "atom"]), vacuum_xy=rnd.choice([0.0, 8.0]),
                vacuum_z=rnd.choice([0.0, 5.0]))


def selftest(n=40, seed=2026):
    """Returns a list of (name, ok, detail). Takes about a minute for n=40."""
    import random
    rnd = random.Random(seed)
    out = []

    # 1. the engine must give EXACTLY the same structure as the frozen reference function
    bad = []
    for i in range(n):
        kw = _random_classic(rnd)
        if kw["mode"] == "nanohub" and kw["direction"] == "axial" and kw["grading"] == "sigmoid":
            kw["p"] = float(round(kw["p"]))  # a fractional p gives a complex number in the last bin (tested below)
        s = rnd.randint(1, 10 ** 6)
        try:
            a = build_structure(seed=s, **kw)
            b = reference_build_structure(seed=s, **kw)
            if not _same_result(a, b):
                bad.append(kw)
        except Exception as e:
            bad.append((kw, repr(e)))
    out.append(("nanoHUB-copy and fixed modes through the engine == frozen reference (every atom, type, bin)",
                not bad, "%d random settings, %d differences" % (n, len(bad))))

    # 2. counts: every recipe gives counts that add to the cell size and sit within 1 of the exact share
    bad = 0
    for _ in range(2000):
        N = rnd.randint(2, 6)
        w = [rnd.random() ** 3 for _ in range(N)]
        c = [x / sum(w) for x in w]
        masses = [rnd.uniform(25, 200) for _ in range(N)]
        cell = rnd.randint(0, 3000)
        for basis in ("mass", "atom"):
            k = counts_from_recipe(c, cell, masses, basis)
            if basis == "atom":
                exact = [ci * cell for ci in c]
            else:
                t = [ci / m for ci, m in zip(c, masses)]
                exact = [ti / sum(t) * cell for ti in t]
            if sum(k) != cell or min(k) < 0 or any(abs(ki - ei) >= 1.0 + 1e-6 for ki, ei in zip(k, exact)):
                bad += 1
    out.append(("atom counts add up to the cell size and stay within 1 atom of the exact share (2 to 6 elements)",
                bad == 0, "2000 random recipes, %d problems" % bad))

    # 3. three or more elements: exact counts in every bin, valid elements, random mixing
    bad_counts, bad_types, zs = 0, 0, []
    for i in range(12):
        N = rnd.randint(3, 5)
        els = [{"symbol": "E%d" % (e + 1), "mass": round(rnd.uniform(40, 70), 2), "a_nm": 0.355} for e in range(N)]
        w0 = [rnd.random() for _ in range(N)]
        w1 = [rnd.random() for _ in range(N)]
        st, en = [x / sum(w0) for x in w0], [x / sum(w1) for x in w1]
        res = engine_build(shape=rnd.choice(list(SHAPES)), radius_nm=1.5, length_nm=6, lattice=rnd.choice(["fcc", "bcc"]),
                           elements=els, start=st, end=en, direction=rnd.choice(["radial", "axial"]),
                           grading=rnd.choice(["power", "exp", "sigmoid"]), p=rnd.uniform(0.5, 2.5),
                           bins=rnd.randint(2, 8), partition="equal", basis=rnd.choice(["mass", "atom"]), seed=i)
        ty = res["atoms"][:, 1].astype(int)
        if ty.min() < 1 or ty.max() > N:
            bad_types += 1
        for m, bi in zip(res["members"], res["bin_info"]):
            if [int((ty[m] == e + 1).sum()) for e in range(N)] != bi["counts"]:
                bad_counts += 1
        groups = np.full(len(ty), -1, dtype=int)
        for j, m in enumerate(res["members"]):
            groups[m] = j
        if len(ty) > 300:
            zs.append(mixing_zscore(ty, groups, neighbour_pairs(res["atoms"][:, 2:5], res["lat_p"], res["lattice"]), n_perm=60))
    out.append(("3 to 5 elements: the count of every element in every bin is exactly as requested",
                bad_counts == 0 and bad_types == 0, "12 random wires, %d count errors, %d bad element numbers" % (bad_counts, bad_types)))
    out.append(("3 to 5 elements: elements are mixed randomly inside each bin (neighbour test, |z| < 4)",
                bool(zs) and max(abs(z) for z in zs) < 4, "worst |z| = %.1f" % (max(abs(z) for z in zs) if zs else float("nan"))))

    # 4. a third element at 0% must give exactly the two-element wire (both orders)
    ok = True
    for direction in ("radial", "axial"):
        two = engine_build(elements=[{"mass": 58.9, "a_nm": 0.36}, {"mass": 58.7, "a_nm": 0.36}], start=[1, 0], end=[0, 1],
                           direction=direction, bins=7, grading="exp", p=1.5, seed=5, lattice_param_nm=0.36)
        three_a = engine_build(elements=[{"mass": 58.9, "a_nm": 0.36}, {"mass": 58.7, "a_nm": 0.36}, {"mass": 55.8, "a_nm": 0.36}],
                               start=[1, 0, 0], end=[0, 1, 0], direction=direction, bins=7, grading="exp", p=1.5, seed=5, lattice_param_nm=0.36)
        three_b = engine_build(elements=[{"mass": 58.9, "a_nm": 0.36}, {"mass": 55.8, "a_nm": 0.36}, {"mass": 58.7, "a_nm": 0.36}],
                               start=[1, 0, 0], end=[0, 0, 1], direction=direction, bins=7, grading="exp", p=1.5, seed=5, lattice_param_nm=0.36)
        ta, tb, t2 = three_a["atoms"][:, 1], three_b["atoms"][:, 1].copy(), two["atoms"][:, 1]
        tb[tb == 3] = 2
        ok &= bool(np.array_equal(ta, t2) and np.array_equal(tb, t2) and np.array_equal(two["atoms"][:, 2:], three_a["atoms"][:, 2:]))
    out.append(("a third element at 0% gives exactly the two-element wire (any position of the unused element)", ok, "radial and axial"))
    # 5. the complex-number quirk of the original program, against numbers computed by Octave 8.4
    #    running the original main.m expressions (b = 13/12, p = 0.44, 90 atoms, masses 49.57 / 69.98):
    #    Octave gave k = 82 - 27i and a curve value of 95.5317
    b = 13 / 12
    k2, ycurve = legacy_counts_complex(grading_raw(b, "sigmoid", 0.44), 90, 49.57, 69.98)
    out.append(("original program's complex-number case reproduced (Octave gave 82 atoms, 95.5317)",
                k2 == 82 and abs(100 * ycurve - 95.5317) < 1e-3, "got %d atoms, %.4f" % (k2, 100 * ycurve)))

    # 6. layered mode
    els2 = [{"symbol": "A", "mass": 58.9, "a_nm": 0.36}, {"symbol": "B", "mass": 58.7, "a_nm": 0.36}]
    ok, detail = True, 0
    half0 = 10 * 0.36 / 2.0
    for direction in ("axial", "radial"):  # one layer must equal the engine's single-direction wire
        for basis in ("mass", "atom"):
            for sample in ("center", "edge"):
                probe = engine_build(shape="cylindrical", radius_nm=1.5, length_nm=8, elements=els2, start=[1, 0], end=[0, 1],
                                     partition="equal", seed=3, lattice_param_nm=0.36)
                length_nm = 2 * probe["n_cells_z"] * half0 / 10.0  # the length the engine really builds
                eng = engine_build(shape="cylindrical", radius_nm=1.5, length_nm=length_nm, elements=els2, start=[1, 0], end=[0, 1],
                                   direction=direction, grading="exp", p=1.5, bins=7, partition="equal", sample=sample,
                                   basis=basis, seed=3, lattice_param_nm=0.36)
                lay = layered_build(shape="cylindrical", radius_nm=1.5, elements=els2, basis=basis, sample=sample, seed=3,
                                    lattice_param_nm=0.36,
                                    layers=[dict(thickness_nm=length_nm, start=[1, 0], end=[0, 1],
                                                 grading="exp", p=1.5, direction=direction, bins=7)])
                ok &= bool(np.array_equal(eng["atoms"], lay["atoms"]))
                detail += 1
    out.append(("layered mode with one layer == the single-direction engine wire (every atom and type)", ok, "%d variants" % detail))

    els3 = [{"symbol": "A", "mass": 58.9, "a_nm": 0.355}, {"symbol": "B", "mass": 51.9, "a_nm": 0.355}, {"symbol": "C", "mass": 58.7, "a_nm": 0.355}]
    stack = [dict(thickness_nm=3.1, start=[1, 0, 0]), dict(thickness_nm=2.2, start=[0, 1, 0]),
             dict(thickness_nm=4.7, start=[0, 0, 1]), dict(thickness_nm=3.3, start=[1, 0, 0])]
    res = layered_build(radius_nm=1.2, elements=els3, layers=stack, seed=2)
    ty, z = res["atoms"][:, 1].astype(int), res["atoms"][:, 4]
    want = [1, 2, 3, 1]
    ok_layers = all((ty[(z >= ly["z0"] - 1e-9) & (z < ly["z1"] - 1e-9)] == want[i]).all() for i, ly in enumerate(res["layers"]))
    half = res["lat_p"] / 2
    ok_thick = all(abs(ly["built_nm"] - ly["asked_nm"]) <= half / 20 + 1e-9 + (half / 10 if i == 3 else 0) for i, ly in enumerate(res["layers"]))
    ok_even = sum(ly["planes"] for ly in res["layers"]) % 2 == 0 and all(c for _, c in verify(res))
    out.append(("uniform layers: every atom has its layer's element, thickness within half an atomic layer, even number of layers",
                ok_layers and ok_thick and ok_even, "4 layers, 3 elements"))

    gr = [dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="sigmoid", p=2, bins=9),
          dict(thickness_nm=6, start=[0, 1], end=[1, 0], grading="power", p=1, bins=9),
          dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="exp", p=1.5, bins=6, direction="radial")]
    res = layered_build(radius_nm=1.6, elements=els2, layers=gr, seed=4)
    share1 = [bi["counts"][0] / bi["n_atoms"] for bi in res["bin_info"] if bi["layer"] == 1]
    share2 = [bi["counts"][0] / bi["n_atoms"] for bi in res["bin_info"] if bi["layer"] == 2]
    mono1 = all(b <= a + 0.03 for a, b in zip(share1, share1[1:])) and share1[0] - share1[-1] > 0.5
    mono2 = all(b >= a - 0.03 for a, b in zip(share2, share2[1:])) and share2[-1] - share2[0] > 0.5
    ty = res["atoms"][:, 1].astype(int)
    d3 = res["atoms"][:, 2:4] - res["axis"][:2]
    rr = np.sqrt((d3 ** 2).sum(axis=1))
    in3 = np.zeros(len(ty), bool)
    in3[np.concatenate([m for m, bi in zip(res["members"], res["bin_info"]) if bi["layer"] == 3])] = True
    radial_ok = rr[in3 & (ty == 2)].mean() > rr[in3 & (ty == 1)].mean()
    out.append(("graded layers: A to B falls/rises as asked, the reversed layer runs the other way, a radial layer puts B at the surface",
                mono1 and mono2 and radial_ok and all(c for _, c in verify(res)), "3 layers, exact counts in every cell"))

    res = layered_build(radius_nm=1.5, elements=els3, seed=6, basis="atom",
                        layers=[dict(thickness_nm=5, start=[1 / 3, 1 / 3, 1 / 3]),
                                dict(thickness_nm=6, start=[0.6, 0.3, 0.1], end=[0.1, 0.3, 0.6], grading="sigmoid", p=2, bins=6),
                                dict(thickness_nm=5, start=[0.2, 0.2, 0.6], end=[0.6, 0.2, 0.2], grading="power", p=1, bins=5, direction="radial")])
    grp = np.zeros(len(res["atoms"]), dtype=int)
    for j, m in enumerate(res["members"]):
        grp[m] = j
    zl = mixing_zscore(res["atoms"][:, 1].astype(int), grp, neighbour_pairs(res["atoms"][:, 2:5], res["lat_p"], res["lattice"]), n_perm=80)
    out.append(("layered ternary wire: elements are mixed randomly inside each cell (neighbour test, |z| < 4)", abs(zl) < 4, "|z| = %.1f" % abs(zl)))

    errs = 0
    for bad in (dict(layers=[dict(thickness_nm=3, start=[0.7, 0.2])]), dict(layers=[dict(thickness_nm=0, start=[1, 0])]),
                dict(layers=[dict(thickness_nm=3, start=[1, 0], direction="sideways")]), dict(layers=[])):
        try:
            layered_build(elements=els2, **bad)
        except ValueError:
            errs += 1
    out.append(("layered mode refuses bad input (recipe not adding to 100%, zero thickness, bad direction, no layers)", errs == 4, "%d of 4 refused" % errs))

    # 7. combined axial-radial mode
    A, B = [1.0, 0.0], [0.0, 1.0]
    ok, n_var = True, 0
    for basis in ("mass", "atom"):
        for sample in ("center", "edge"):
            probe = engine_build(shape="cylindrical", radius_nm=1.5, length_nm=8, elements=els2, start=A, end=B,
                                 partition="equal", seed=3, lattice_param_nm=0.36)
            length_nm = 2 * probe["n_cells_z"] * half0 / 10.0
            # no radial change and one ring  == the engine's axial wire
            ax_e = engine_build(shape="cylindrical", radius_nm=1.5, length_nm=length_nm, elements=els2, start=A, end=B,
                                direction="axial", grading="sigmoid", p=2, bins=9, partition="equal", sample=sample,
                                basis=basis, seed=7, lattice_param_nm=0.36)
            ax_c = combined_build(shape="cylindrical", radius_nm=1.5, length_nm=length_nm, elements=els2,
                                  corners=[A, A, B, B], grading_z="sigmoid", p_z=2, bins_z=9, bins_r=1, basis=basis,
                                  sample=sample, seed=7, lattice_param_nm=0.36)
            # no axial change and one slice == the engine's radial wire
            rd_e = engine_build(shape="cylindrical", radius_nm=1.5, length_nm=length_nm, elements=els2, start=A, end=B,
                                direction="radial", grading="exp", p=1.5, bins=6, partition="equal", sample=sample,
                                basis=basis, seed=7, lattice_param_nm=0.36)
            rd_c = combined_build(shape="cylindrical", radius_nm=1.5, length_nm=length_nm, elements=els2,
                                  corners=[A, B, A, B], grading_r="exp", p_r=1.5, bins_r=6, bins_z=1, basis=basis,
                                  sample=sample, seed=7, lattice_param_nm=0.36)
            ok &= bool(np.array_equal(ax_e["atoms"], ax_c["atoms"]) and np.array_equal(rd_e["atoms"], rd_c["atoms"]))
            n_var += 2
    out.append(("combined mode reduces exactly to the axial wire (1 ring) and to the radial wire (1 slice)", ok, "%d variants" % n_var))

    bad_counts, bad_recipe, zs = 0, 0, []
    for i in range(10):
        N = rnd.randint(2, 5)
        els = [{"symbol": "E%d" % (e + 1), "mass": round(rnd.uniform(40, 70), 2), "a_nm": 0.355} for e in range(N)]

        def rec():
            w = [rnd.random() ** 2 for _ in range(N)]
            return [x / sum(w) for x in w]
        res = combined_build(shape=rnd.choice(["cylindrical", "cuboid"]), radius_nm=1.4, length_nm=6, elements=els,
                             corners=[rec(), rec(), rec(), rec()], grading_r=rnd.choice(["power", "exp", "sigmoid"]),
                             grading_z=rnd.choice(["power", "exp", "sigmoid"]), p_r=rnd.uniform(0.5, 2.5), p_z=rnd.uniform(0.5, 2.5),
                             bins_r=rnd.randint(1, 5), bins_z=rnd.randint(1, 6), basis=rnd.choice(["mass", "atom"]),
                             lattice=rnd.choice(["fcc", "bcc"]), seed=i)
        if not all(c for _, c in verify(res)):
            bad_counts += 1
        if any(abs(sum(bi["recipe"]) - 1.0) > 1e-9 or min(bi["recipe"]) < -1e-12 for bi in res["bin_info"]):
            bad_recipe += 1
        if len(res["atoms"]) > 400 and N >= 3:
            grp = np.zeros(len(res["atoms"]), dtype=int)
            for j, m in enumerate(res["members"]):
                grp[m] = j
            zs.append(mixing_zscore(res["atoms"][:, 1].astype(int), grp,
                                    neighbour_pairs(res["atoms"][:, 2:5], res["lat_p"], res["lattice"]), n_perm=60))
    out.append(("combined mode, 2 to 5 elements: exact atom counts in every cell and every mix valid (0 to 100 %)",
                bad_counts == 0 and bad_recipe == 0, "10 random wires, %d count errors, %d invalid mixes" % (bad_counts, bad_recipe)))
    out.append(("combined mode: elements are mixed randomly inside each cell (neighbour test, |z| < 4)",
                bool(zs) and max(abs(z) for z in zs) < 4, "worst |z| = %.1f" % (max(abs(z) for z in zs) if zs else float("nan"))))

    res = combined_build(radius_nm=2.0, length_nm=12, elements=els2, corners=[A, [0.5, 0.5], [0.5, 0.5], B], bins_z=6, bins_r=4,
                         grading_z="power", grading_r="power", sample="edge", basis="atom", seed=1)
    tcomb = res["atoms"][:, 1].astype(int)  # measured on the real atoms, not on the requested counts
    share = {(bi["iz"], bi["ir"]): float((tcomb[m] == 1).mean()) for bi, m in zip(res["bin_info"], res["members"]) if bi["n_atoms"]}
    mono_z = all(share[(z + 1, r)] <= share[(z, r)] + 0.04 for z in range(1, 6) for r in range(1, 5))
    mono_r = all(share[(z, r + 1)] <= share[(z, r)] + 0.04 for z in range(1, 7) for r in range(1, 4))
    corner_ok = share[(1, 1)] == 1.0 and share[(6, 4)] == 0.0 and all(c for _, c in verify(res))
    out.append(("combined mode: element A falls along the wire and from core to surface, corner cells hold the corner recipes",
                mono_z and mono_r and corner_ok, "6 slices x 4 rings"))

    errs = 0
    for bad in (dict(corners=[A, A, B]), dict(corners=[A, A, B, [0.3, 0.3]]), dict(corners=[A, A, B, B], bins_r=0),
                dict(corners=[A, A, B, B], shape="sphere")):
        try:
            combined_build(elements=els2, **bad)
        except ValueError:
            errs += 1
    out.append(("combined mode refuses bad input (wrong number of corners, recipe not adding to 100%, no rings, sphere)",
                errs == 4, "%d of 4 refused" % errs))
    # 8. per-layer lattice constants: spacing along the wire follows the composition
    els_nc = [{"symbol": "Ni", "mass": 58.693, "a_nm": 0.3524}, {"symbol": "Cu", "mass": 63.546, "a_nm": 0.3615}]
    stack_nc = [dict(thickness_nm=4, start=[1, 0]), dict(thickness_nm=4, start=[0, 1]),
                dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="sigmoid", p=2, bins=6)]
    same = layered_build(radius_nm=1.5, elements=els_nc, layers=stack_nc, z_spacing="same", seed=9)
    veg = layered_build(radius_nm=1.5, elements=els_nc, layers=stack_nc, z_spacing="vegard", seed=9)
    plv, zv = veg["plane"], veg["atoms"][:, 4]
    Kp = int(plv.max()) + 1
    zpl = np.array([zv[plv == k].mean() for k in range(Kp)])
    dz = np.diff(zpl)
    p1 = veg["layers"][0]["planes"]
    p2 = p1 + veg["layers"][1]["planes"]
    ok_pure = (np.allclose(dz[:p1 - 1], 10 * 0.3524 / 2, atol=1e-9) and np.allclose(dz[p1:p2 - 1], 10 * 0.3615 / 2, atol=1e-9)
               and abs(dz[p1 - 1] - 0.5 * (10 * 0.3524 / 2 + 10 * 0.3615 / 2)) < 1e-9)
    graded = dz[p2:]
    ok_graded = bool(np.all(np.diff(graded) >= -1e-9)) and graded[0] < graded[-1] and graded[0] >= 10 * 0.3524 / 2 - 1e-9 and graded[-1] <= 10 * 0.3615 / 2 + 1e-9
    out.append(("per-layer spacing: pure Ni and pure Cu layers keep their own spacing, the interface sits in between, a graded layer changes smoothly",
                bool(ok_pure and ok_graded), "layer spacings %.4f / %.4f A" % (dz[0], dz[p1])))
    ok_same = (np.array_equal(same["atoms"][:, 1], veg["atoms"][:, 1]) and np.array_equal(same["atoms"][:, 2:4], veg["atoms"][:, 2:4])
               and np.array_equal(same["plane"], veg["plane"]) and all(c for _, c in verify(veg)))
    out.append(("per-layer spacing changes only the heights: the same atoms get the same elements and the same x, y positions",
                bool(ok_same), "exact counts still hold"))
    box_len = veg["box"]["z"][1] - veg["box"]["z"][0]
    a_top = vegard_constant(veg["bin_info"][-1]["counts"], [0.3524, 0.3615])  # composition of the topmost cell
    gap_expected = 0.5 * (10 * a_top / 2 + 10 * 0.3524 / 2)  # mean of the top and bottom layer spacings
    out.append(("per-layer spacing: the box ends one interface spacing above the top layer (periodic z has no seam)",
                abs(box_len - zpl[-1] - gap_expected) < 1e-9, "gap %.4f A" % (box_len - zpl[-1])))

    ovr = layered_build(radius_nm=1.5, elements=els_nc, z_spacing="vegard", seed=9,
                        layers=[dict(thickness_nm=4, start=[1, 0], z_a_nm=0.40), dict(thickness_nm=4, start=[0, 1], z_a_nm=0.33)])
    zo = ovr["atoms"][:, 4]
    po = ovr["plane"]
    zq = np.array([zo[po == k].mean() for k in range(int(po.max()) + 1)])
    a_first = np.diff(zq)[:ovr["layers"][0]["planes"] - 2].mean()
    a_second = np.diff(zq)[ovr["layers"][0]["planes"] + 1:].mean()
    out.append(("per-layer spacing: a lattice constant typed for a layer overrides the composition value",
                abs(a_first - 10 * 0.40 / 2) < 1e-9 and abs(a_second - 10 * 0.33 / 2) < 1e-9, "%.3f / %.3f A" % (a_first, a_second)))

    eq_els = [{"symbol": "A", "mass": 58.9, "a_nm": 0.36}, {"symbol": "B", "mass": 58.7, "a_nm": 0.36}]
    lay_eq = [dict(thickness_nm=5, start=[1, 0], end=[0, 1], grading="power", p=1, bins=5)]
    s_same = layered_build(radius_nm=1.5, elements=eq_els, layers=lay_eq, z_spacing="same", seed=3)
    s_veg = layered_build(radius_nm=1.5, elements=eq_els, layers=lay_eq, z_spacing="vegard", seed=3)
    out.append(("per-layer spacing: when every element has the same constant the result equals the single-constant wire",
                bool(np.array_equal(s_same["atoms"][:, 1], s_veg["atoms"][:, 1]) and np.abs(s_same["atoms"][:, 2:] - s_veg["atoms"][:, 2:]).max() < 1e-9),
                "max height difference %.1e A" % np.abs(s_same["atoms"][:, 4] - s_veg["atoms"][:, 4]).max()))

    cv = combined_build(radius_nm=1.5, length_nm=8, elements=els_nc, corners=[[1, 0], [1, 0], [0, 1], [0, 1]], bins_z=6, bins_r=2,
                        z_spacing="vegard", seed=2)
    az = [bi["a_z_nm"] for bi in cv["bin_info"] if bi["ir"] == 1]
    out.append(("combined mode with per-slice spacing: the spacing grows smoothly from Ni to Cu along the wire",
                all(b >= a - 1e-12 for a, b in zip(az, az[1:])) and az[0] < az[-1] and 0.3524 <= az[0] and az[-1] <= 0.3615,
                "%.4f -> %.4f nm" % (az[0], az[-1])))
    bad = 0
    for kw in (dict(z_spacing="stretch"), dict(layers=[dict(thickness_nm=3, start=[1, 0], z_a_nm=-0.3)], z_spacing="vegard")):
        try:
            layered_build(elements=els_nc, **({"layers": [dict(thickness_nm=3, start=[1, 0])]} | kw))
        except ValueError:
            bad += 1
    out.append(("per-layer spacing refuses bad input (unknown mode, zero or negative lattice constant)", bad == 2, "%d of 2 refused" % bad))
    # 9. crystal layers: FCC (111) and HCP (0001) layers joined along the wire
    def geo_pairs(pos, d):
        I, J = [], []
        for a0 in range(0, len(pos), 400):
            dd = np.sqrt(((pos[a0:a0 + 400, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
            ii, jj = np.nonzero(dd <= 1.15 * d)
            ii = ii + a0
            keep = ii < jj
            I.append(ii[keep]); J.append(jj[keep])
        return np.concatenate(I), np.concatenate(J)

    def bulk_atoms(res, layer_no, margin_planes=2, margin_rad=1.3):
        ly = res["layers"][layer_no]
        P0 = sum(x["planes"] for x in res["layers"][:layer_no])
        pl = res["plane"]
        d_ = res["atoms"][:, 2:4] - res["axis"][:2]
        rr = np.sqrt((d_ ** 2).sum(axis=1))
        ok = (pl >= P0 + margin_planes) & (pl < P0 + ly["planes"] - margin_planes) & (rr <= rr.max() - margin_rad * res["d_xy"])
        return np.nonzero(ok)[0]

    def local_checks(res, layer_no, kind):
        pos = res["atoms"][:, 2:5]
        d_xy = res["d_xy"]
        idx = bulk_atoms(res, layer_no)
        if len(idx) < 5:
            return False
        good = True
        for i in idx[::max(1, len(idx) // 25)]:
            v = pos - pos[i]
            dist = np.sqrt((v ** 2).sum(axis=1))
            near = (dist <= 1.15 * d_xy) & (dist > 1e-9)
            good &= int(near.sum()) == 12
            same_plane = (res["plane"] == res["plane"][i]) & (dist > 1e-9)
            good &= abs(dist[same_plane].min() - d_xy) < 1e-9
            up = np.round(v[near & (v[:, 2] > 0), :2], 3)
            dn = np.round(v[near & (v[:, 2] < 0), :2], 3)
            if len(up) != 3 or len(dn) != 3:
                good = False
                continue
            sort = lambda a: a[np.lexsort((a[:, 1], a[:, 0]))]
            if kind == "fcc111":   # A B C: the triangle above is the triangle below turned round
                good &= bool(np.allclose(sort(up), sort(-dn), atol=1e-3))
            else:                  # A B A: the triangles above and below coincide
                good &= bool(np.allclose(sort(up), sort(dn), atol=1e-3))
        return good

    ce = [{"symbol": "A", "mass": 58.7}, {"symbol": "B", "mass": 51.9}]
    fcc1 = crystal_layers_build(radius_nm=1.1, elements=ce, layers=[dict(crystal="fcc111", a_nm=0.3524, thickness_nm=6, start=[1, 0])], seed=1)
    hcp1 = crystal_layers_build(radius_nm=1.1, elements=ce, layers=[dict(crystal="hcp0001", a_nm=0.2951, c_nm=0.4683, thickness_nm=6, start=[1, 0])], seed=1)
    ok_f, ok_h = local_checks(fcc1, 0, "fcc111"), local_checks(hcp1, 0, "hcp0001")
    ok_f &= (len(fcc1["letters"]) % 3 == 0 and (fcc1["letters"][0] - fcc1["letters"][-1]) % 3 == 1)  # ...A B C | A B C... across the seam
    ok_h &= (len(hcp1["letters"]) % 2 == 0 and hcp1["letters"][0] != hcp1["letters"][-1])           # ...A B | A B... across the seam
    sp_f = np.diff([fcc1["atoms"][fcc1["plane"] == k, 4].mean() for k in range(8)])
    sp_h = np.diff([hcp1["atoms"][hcp1["plane"] == k, 4].mean() for k in range(8)])
    ok_sp = np.allclose(sp_f, 10 * 0.3524 / math.sqrt(3), atol=1e-9) and np.allclose(sp_h, 10 * 0.4683 / 2, atol=1e-9)
    out.append(("crystal layers: FCC (111) and HCP (0001) layers have 12 neighbours, the right row spacing, plane spacing and ABC / ABA stacking",
                bool(ok_f and ok_h and ok_sp), "FCC a/sqrt(3) = %.4f A, HCP c/2 = %.4f A" % (sp_f[0], sp_h[0])))

    ti_al = [dict(crystal="hcp0001", a_nm=0.2951, c_nm=0.4683, thickness_nm=3, start=[1, 0]), dict(crystal="fcc111", a_nm=0.4050, thickness_nm=3, start=[0, 1]),
             dict(crystal="hcp0001", a_nm=0.2951, c_nm=0.4683, thickness_nm=3, start=[1, 0]), dict(crystal="fcc111", a_nm=0.4050, thickness_nm=3, start=[0, 1])]
    ta = crystal_layers_build(radius_nm=1.0, elements=ce, layers=ti_al, seed=2)
    d_ti, d_al = 10 * 0.2951, 10 * 0.4050 / math.sqrt(2)
    ok_mean = abs(ta["d_xy"] - 0.5 * (d_ti + d_al)) < 1e-12 and ta["layers"][0]["strain_pct"] < 0 < ta["layers"][1]["strain_pct"]
    ok_bulk = all(local_checks(ta, k, ta["layers"][k]["crystal"]) for k in range(4))
    pos = ta["atoms"][:, 2:5]
    mind = 1e9
    for a0 in range(0, len(pos), 400):
        blk = pos[a0:a0 + 400]
        for shift in (0.0, ta["z_period"], -ta["z_period"]):
            dd = np.sqrt(((blk[:, None, :] - (pos + [0, 0, shift])[None, :, :]) ** 2).sum(axis=2))
            if shift == 0.0:
                dd[np.arange(len(blk)), np.arange(a0, a0 + len(blk))] = 1e9
            mind = min(mind, float(dd.min()))
    lt = ta["letters"]
    ok_letters = all(a != b for a, b in zip(lt, lt[1:])) and lt[-1] != lt[0]
    out.append(("crystal layers Ti/Al: shared in-plane spacing is the mean (layers strained by opposite amounts), no atoms overlap anywhere "
                "(also across the periodic seam), neighbouring planes never share a stacking letter",
                bool(ok_mean and ok_bulk and mind >= 0.9 * ta["d_xy"] and ok_letters),
                "strain %+.2f / %+.2f %%, closest atoms %.3f A" % (ta["layers"][0]["strain_pct"], ta["layers"][1]["strain_pct"], mind)))
    out.append(("crystal layers: every atom of a pure layer has that layer's element, exact counts in every cell",
                all(c for _, c in verify(ta)) and all((ta["atoms"][(ta["plane"] >= sum(x["planes"] for x in ta["layers"][:k])) &
                                                                     (ta["plane"] < sum(x["planes"] for x in ta["layers"][:k + 1])), 1] == (1 if k % 2 == 0 else 2)).all()
                                                       for k in range(4)), "4 layers"))

    c3 = [{"symbol": "A", "mass": 58.7}, {"symbol": "B", "mass": 51.9}, {"symbol": "C", "mass": 58.9}]
    mixed = crystal_layers_build(radius_nm=1.4, elements=c3, basis="atom", seed=5, layers=[
        dict(crystal="fcc111", a_nm=0.354, thickness_nm=4, start=[0.5, 0.3, 0.2], end=[0.1, 0.3, 0.6], grading="sigmoid", p=2, bins=5),
        dict(crystal="hcp0001", a_nm=0.2951, c_nm=0.4683, thickness_nm=4, start=[1 / 3, 1 / 3, 1 / 3]),
        dict(crystal="fcc111", a_nm=0.354, thickness_nm=4, start=[0.7, 0.2, 0.1], end=[0.2, 0.2, 0.6], grading="power", p=1, bins=4, direction="radial")])
    grp = np.zeros(len(mixed["atoms"]), dtype=int)
    for j, mm in enumerate(mixed["members"]):
        grp[mm] = j
    zc = mixing_zscore(mixed["atoms"][:, 1].astype(int), grp, geo_pairs(mixed["atoms"][:, 2:5], mixed["d_xy"]), n_perm=60)
    out.append(("crystal layers, 3 elements with graded and radial layers: exact counts in every cell, mixed randomly (|z| < 4)",
                all(c for _, c in verify(mixed)) and abs(zc) < 4, "|z| = %.1f" % abs(zc)))

    first = crystal_layers_build(radius_nm=1.0, elements=ce, layers=ti_al, inplane="first", seed=2)
    valu = crystal_layers_build(radius_nm=1.0, elements=ce, layers=ti_al, inplane="value", inplane_nm=0.29, seed=2)
    out.append(("crystal layers: in-plane spacing options (match the bottom layer, or a typed value) give the expected strains",
                abs(first["layers"][0]["strain_pct"]) < 1e-9 and abs(valu["d_xy"] - 2.9) < 1e-12 and abs(valu["layers"][0]["strain_pct"] - 100 * (2.9 / d_ti - 1)) < 1e-9,
                "first layer strain 0, typed 0.29 nm"))

    bad = 0
    good_layer = dict(crystal="fcc111", a_nm=0.35, thickness_nm=3, start=[1, 0])
    for kw in (dict(layers=[dict(good_layer, crystal="bcc001")]), dict(layers=[dict(good_layer, a_nm=0)]),
               dict(layers=[dict(good_layer, crystal="hcp0001", c_nm=-1)]), dict(layers=[]), dict(layers=[good_layer], shape="sphere"),
               dict(layers=[good_layer], inplane="value"), dict(layers=[dict(good_layer, thickness_nm=0)]),
               dict(layers=[dict(good_layer, start=[0.5, 0.2])])):
        try:
            crystal_layers_build(elements=ce, **kw)
        except ValueError:
            bad += 1
    out.append(("crystal layers refuse bad input (unknown crystal, zero constant, negative c, no layers, sphere, missing spacing, zero thickness, bad recipe)",
                bad == 8, "%d of 8 refused" % bad))
    # 10. Custom function: pure wires and solid layers (no mixing) in every mode
    gf = grading_fraction
    ok_fn = (abs(gf(0.5, "custom", 1, "x^2") - 0.25) < 1e-12 and abs(gf(0.5, "custom", 1, "x.^2") - 0.25) < 1e-12
             and gf(0.69, "custom", 1, "step(x,0.7)") == 0.0 and gf(0.7, "custom", 1, "step(x,0.7)") == 1.0
             and gf(0.8, "custom", 1, "x>0.7") == 1.0 and gf(0.2, "custom", 1, "x>0.7") == 0.0
             and gf(0.45, "custom", 1, "box(x,0.3,0.6)") == 1.0 and gf(0.6, "custom", 1, "box(x,0.3,0.6)") == 0.0
             and abs(gf(0.3, "custom", 1, "min(1,2*x)") - 0.6) < 1e-12 and abs(gf(0.5, "custom", 1, "tanh(x)") - math.tanh(0.5)) < 1e-12
             and gf(0.5, "custom", 1, "(x>0.2) and (x<0.8)") == 1.0 and gf(0.5, "custom", 1, "(x>0.2)&(x<0.8)") == 1.0
             and gf(0.5, "custom", 1, "-3") == 0.0 and gf(0.5, "custom", 1, "7") == 1.0 and gf(0.5, "custom", 1, "0") == 0.0)
    refused = 0
    for bad in ("foo(x)", "1/(x-0.5)", "", "x +", "log(x-1)"):
        try:
            gf(0.5, "custom", 1, bad)
        except ValueError:
            refused += 1
    out.append(("custom function: x^2, x.^2, step, box, comparisons, min, tanh give the right values; values outside 0..1 are cut; bad formulas are refused",
                bool(ok_fn and refused == 5), "%d of 5 bad formulas refused" % refused))

    ce = [{"symbol": "A", "mass": 58.9, "a_nm": 0.355}, {"symbol": "B", "mass": 51.9, "a_nm": 0.355}]

    def lay_custom(formula, **kw):
        return layered_build(radius_nm=1.5, elements=ce, seed=4, layers=[dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="custom",
                                                                             custom=formula, **kw)])

    def cells_pure(res):
        return all(sum(1 for c in bi["counts"] if c) <= 1 for bi in res["bin_info"])

    p0 = layered_build(radius_nm=1.5, elements=ce, seed=4, layers=[dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="custom", custom="0", bins=6)])
    u0 = layered_build(radius_nm=1.5, elements=ce, seed=4, layers=[dict(thickness_nm=6, start=[1, 0])])
    p1 = layered_build(radius_nm=1.5, elements=ce, seed=4, layers=[dict(thickness_nm=6, start=[1, 0], end=[0, 1], grading="custom", custom="1", bins=6)])
    u1 = layered_build(radius_nm=1.5, elements=ce, seed=4, layers=[dict(thickness_nm=6, start=[0, 1])])
    out.append(("layers, custom f(x) = 0 or 1: a pure wire, identical to a uniform layer of that element (every atom)",
                bool(np.array_equal(p0["atoms"], u0["atoms"]) and np.array_equal(p1["atoms"], u1["atoms"]) and (p0["atoms"][:, 1] == 1).all() and (p1["atoms"][:, 1] == 2).all()),
                "x = 0 and x = 1"))

    sh = lay_custom("x>0.7", direction="radial", bins=10)
    tys = sh["atoms"][:, 1].astype(int)
    dd = sh["atoms"][:, 2:4] - sh["axis"][:2]
    rr = np.hypot(dd[:, 0], dd[:, 1])
    ring = np.minimum(np.floor(rr / (rr.max() / 10)).astype(int), 9)
    shell = ring >= 7
    out.append(("layers, custom radial step x>0.7: a SOLID outer shell (rings 8 to 10 all element B, the core all element A, no mixed cell)",
                bool((tys[shell] == 2).all() and (tys[~shell] == 1).all() and 0 < shell.sum() < len(tys) and cells_pure(sh) and all(c for _, c in verify(sh))),
                "%d shell atoms, %d core atoms" % (shell.sum(), (~shell).sum())))

    bd = lay_custom("box(x,0.3,0.7)", direction="axial", bins=10)
    tyb = bd["atoms"][:, 1].astype(int)
    ok_band = all(((tyb[m] == 2).all() if 3 <= j <= 6 else (tyb[m] == 1).all()) for j, m in enumerate(bd["members"]))
    out.append(("layers, custom axial band box(x,0.3,0.7): a SOLID band of element B in cells 4 to 7, element A everywhere else",
                bool(ok_band and cells_pure(bd) and all(c for _, c in verify(bd))), "10 slices"))

    def comb_custom(corners, fz, fr, bz, br):
        return combined_build(radius_nm=1.5, length_nm=8, elements=ce, corners=corners, grading_r="custom", custom_r=fr, grading_z="custom",
                              custom_z=fz, bins_r=br, bins_z=bz, seed=5)

    def region_ok(res, want):  # want(slice, ring) -> element number; every cell must be pure and equal to it
        t = res["atoms"][:, 1].astype(int)
        return all(len(m) > 0 and (t[m] == want(bi["iz"], bi["ir"])).all() for m, bi in zip(res["members"], res["bin_info"])) and cells_pure(res) \
            and all(c for _, c in verify(res))

    A, B = [1, 0], [0, 1]
    c_a = comb_custom([A, B, A, B], "0", "step(x,0.5)", 3, 2)     # across only: solid core (A) and solid shell (B) in every slice
    c_b = comb_custom([A, A, B, B], "step(x,0.5)", "0", 2, 3)     # along only: solid bottom (A) and solid top (B) in every ring
    c_c = comb_custom([A, B, B, B], "x>=0.5", "x>=0.5", 2, 2)    # both: only core-bottom is A
    ok_c = region_ok(c_a, lambda z, r: 1 if r == 1 else 2) and region_ok(c_b, lambda z, r: 1 if z == 1 else 2) \
        and region_ok(c_c, lambda z, r: 1 if (z, r) == (1, 1) else 2)
    out.append(("combined, custom in each direction alone and in both: solid core and shell, solid bottom and top, a solid corner region, none mixed",
                bool(ok_c), "3 wires"))

    half_mix = lay_custom("0.4", direction="axial", bins=3)
    th = half_mix["atoms"][:, 1].astype(int)
    ok_mix = all(abs(int((th[m] == 2).sum()) - 0.4 * len(m)) <= 1.0 for m in half_mix["members"])
    out.append(("custom value between 0 and 1 (f(x) = 0.4) gives a random alloy of that composition in every cell (so only 0 and 1 give solid regions)",
                bool(ok_mix and not cells_pure(half_mix)), "40 % element B in every cell"))

    cr = crystal_layers_build(radius_nm=1.2, elements=[{"symbol": "A", "mass": 58.7}, {"symbol": "B", "mass": 51.9}], seed=2,
                              layers=[dict(crystal="fcc111", a_nm=0.3524, thickness_nm=6, start=[1, 0], end=[0, 1], grading="custom", custom="x>=0.5", bins=4)])
    tcr = cr["atoms"][:, 1].astype(int)
    ok_cr = all(((tcr[m] == 2).all() if j >= 2 else (tcr[m] == 1).all()) for j, m in enumerate(cr["members"]))
    out.append(("crystal layers: a graded crystal layer with custom x>=0.5 gives a solid lower half and a solid upper half",
                bool(ok_cr and all(c for _, c in verify(cr))), "4 slices"))

    refused = 0
    for fn in (lambda: lay_custom("foo(x)", bins=4), lambda: lay_custom("", bins=4), lambda: lay_custom("1/(x-0.5)", bins=4),
               lambda: combined_build(elements=ce, corners=[[1, 0]] * 4, grading_r="custom", custom_r="", bins_r=2, bins_z=2),
               lambda: combined_build(elements=ce, corners=[[1, 0], [0, 1], [1, 0], [0, 1]], grading_z="custom", custom_z="bar", bins_r=2, bins_z=2)):
        try:
            fn()
        except ValueError:
            refused += 1
    out.append(("layered and combined modes refuse a missing, unknown or unsafe custom formula before building", refused == 5, "%d of 5 refused" % refused))
    return out


LAYERED_EXAMPLE = {
    "shape": "cylindrical", "radius_nm": 1.5, "lattice": "fcc", "basis": "atom", "sample": "center",
    "elements": [{"symbol": "Ni", "mass": 58.693, "a_nm": 0.3524}, {"symbol": "Co", "mass": 58.933, "a_nm": 0.3544}],
    "layers": [
        {"thickness_nm": 4, "start": [1, 0]},
        {"thickness_nm": 6, "start": [1, 0], "end": [0, 1], "grading": "sigmoid", "p": 2, "bins": 8, "direction": "axial"},
        {"thickness_nm": 4, "start": [0, 1]},
        {"thickness_nm": 6, "start": [0, 1], "end": [1, 0], "grading": "power", "p": 1, "bins": 6, "direction": "axial"}],
    "vacuum_xy": 0, "vacuum_z": 0, "seed": 1}


def cli_layered(argv=None):
    import json
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py layered",
                                 description="Build a layered / multi-element wire from a JSON description")
    ap.add_argument("--spec", help="JSON file describing the wire")
    ap.add_argument("--example", action="store_true", help="print an example JSON description and stop")
    ap.add_argument("--out", default="layered_wire.nw")
    a = ap.parse_args(argv)
    if a.example or not a.spec:
        print(json.dumps(LAYERED_EXAMPLE, indent=2))
        return
    with open(a.spec) as f:
        spec = json.load(f)
    res = layered_build(**spec)
    write_lammps_data(a.out, res)
    print("\n".join(summary_lines(res)))
    print("written:", a.out)


COMBINED_EXAMPLE = {
    "shape": "cylindrical", "radius_nm": 2.0, "length_nm": 10, "lattice": "fcc", "basis": "atom", "sample": "center",
    "elements": [{"symbol": "Ni", "mass": 58.693, "a_nm": 0.3524}, {"symbol": "Co", "mass": 58.933, "a_nm": 0.3544}],
    "corners": [[1, 0], [0.7, 0.3], [0.3, 0.7], [0, 1]],
    "grading_z": "sigmoid", "p_z": 2, "bins_z": 10, "grading_r": "power", "p_r": 1, "bins_r": 5,
    "vacuum_xy": 0, "vacuum_z": 0, "seed": 1}


def cli_combined(argv=None):
    import json
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py combined",
                                 description="Build a wire graded along AND across from a JSON description")
    ap.add_argument("--spec", help="JSON file describing the wire")
    ap.add_argument("--example", action="store_true", help="print an example JSON description and stop")
    ap.add_argument("--out", default="combined_wire.nw")
    a = ap.parse_args(argv)
    if a.example or not a.spec:
        print(json.dumps(COMBINED_EXAMPLE, indent=2))
        return
    with open(a.spec) as f:
        spec = json.load(f)
    res = combined_build(**spec)
    write_lammps_data(a.out, res)
    print("\n".join(summary_lines(res)))
    print("written:", a.out)


CRYSTALS_EXAMPLE = {
    "shape": "cylindrical", "radius_nm": 1.5, "inplane": "mean", "basis": "atom", "sample": "center",
    "elements": [{"symbol": "Ti", "mass": 47.867}, {"symbol": "Al", "mass": 26.982}],
    "layers": [
        {"crystal": "hcp0001", "a_nm": 0.2951, "c_nm": 0.4683, "thickness_nm": 4, "start": [1, 0]},
        {"crystal": "fcc111", "a_nm": 0.405, "thickness_nm": 4, "start": [0, 1]},
        {"crystal": "hcp0001", "a_nm": 0.2951, "c_nm": 0.4683, "thickness_nm": 4, "start": [1, 0]},
        {"crystal": "fcc111", "a_nm": 0.405, "thickness_nm": 4, "start": [0, 1]}],
    "vacuum_xy": 0, "vacuum_z": 0, "seed": 1}


def cli_crystals(argv=None):
    import json
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py crystals",
                                 description="Build a wire from FCC (111) / HCP (0001) crystal layers from a JSON description")
    ap.add_argument("--spec", help="JSON file describing the wire")
    ap.add_argument("--example", action="store_true", help="print an example JSON description and stop")
    ap.add_argument("--out", default="crystal_wire.nw")
    a = ap.parse_args(argv)
    if a.example or not a.spec:
        print(json.dumps(CRYSTALS_EXAMPLE, indent=2))
        return
    with open(a.spec) as f:
        spec = json.load(f)
    res = crystal_layers_build(**spec)
    write_lammps_data(a.out, res)
    print("\n".join(summary_lines(res)))
    print("written:", a.out)


def cli_selftest(argv=None):
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py selftest")
    ap.add_argument("--n", type=int, default=40, help="how many random settings to compare")
    a = ap.parse_args(argv)
    t0 = time.time()
    results = selftest(a.n)
    for name, ok, detail in results:
        print("[%s] %s  (%s)" % ("PASS" if ok else "FAIL", name, detail))
    allok = all(ok for _, ok, _ in results)
    print("\nSELF-TEST: %s  (%.0f s)" % ("ALL CHECKS PASSED" if allok else "SOMETHING FAILED", time.time() - t0))
    sys.exit(0 if allok else 1)


# =====================================================================================
# PART 2 - comparison with a reference data file
# =====================================================================================
def read_data(path):
    with open(path) as f:
        lines = f.read().split("\n")
    k = next(i for i, l in enumerate(lines) if l.strip() == "Atoms")
    header = [l.rstrip() for l in lines[:k]]
    n = int(next(l for l in header if l.strip().endswith("atoms")).split()[0])
    atoms = np.loadtxt(lines[k + 2:k + 2 + n])
    return header, atoms


def neighbour_pairs(pos, lat_p, lattice):
    """All nearest-neighbour atom pairs (each pair once), found on the lattice grid."""
    half = lat_p / 2.0
    ijk = np.rint(pos / half).astype(np.int64)
    ijk -= ijk.min(axis=0)
    K = int(ijk.max()) + 4
    key = (ijk[:, 0] + 1) * K * K + (ijk[:, 1] + 1) * K + (ijk[:, 2] + 1)
    order = np.argsort(key)
    skey = key[order]
    if lattice == "fcc":  # 12 neighbours (a/2)(+-1,+-1,0) in all orders
        offs = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
    else:  # bcc: 8 neighbours (a/2)(+-1,+-1,+-1)
        offs = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]
    I, J = [], []
    for dx, dy, dz in offs:  # one direction of every pair is enough
        nk = key + dx * K * K + dy * K + dz
        at = np.searchsorted(skey, nk)
        at[at >= len(skey)] = len(skey) - 1
        hit = skey[at] == nk
        I.append(np.nonzero(hit)[0])
        J.append(order[at[hit]])
    return np.concatenate(I), np.concatenate(J)


def mixing_zscore(types, groups, pairs, n_perm=300, seed=0):
    """Are the two kinds of atoms mixed randomly inside every bin?
    Count unlike nearest-neighbour pairs, then shuffle the atoms inside each bin many times
    and see where the real count falls (z-score; near 0 = random, large = clumped/ordered)."""
    I, J = pairs
    obs = int((types[I] != types[J]).sum())
    rng = np.random.default_rng(seed)
    idx = [np.nonzero(groups == g)[0] for g in np.unique(groups)]
    null = np.empty(n_perm)
    for r in range(n_perm):
        t = types.copy()
        for ix in idx:
            t[ix] = types[rng.permutation(ix)]
        null[r] = (t[I] != t[J]).sum()
    sd = null.std()
    return (obs - null.mean()) / sd if sd > 0 else 0.0


def run_checks(res, ref_path):
    """Compare a build result (standard mode) with a reference data file.
    Returns (lines, ok). Each line is (text, tag) with tag 'pass', 'fail' or None."""
    out = []
    state = {"ok": True}

    def add(text, tag=None):
        out.append((text, tag))

    def report(name, ok, detail=""):
        state["ok"] &= bool(ok)
        add("[%s] %s %s" % ("PASS" if ok else "FAIL", name, detail), "pass" if ok else "fail")

    mine = res["atoms"]
    fd, tmp = tempfile.mkstemp(suffix=".nw")
    os.close(fd)
    try:
        write_lammps_data(tmp, res)
        h_mine, _ = read_data(tmp)
    finally:
        os.remove(tmp)
    h_ref, ref = read_data(ref_path)

    # 1. header
    report("header identical (count, box, masses)", h_ref[1:] == h_mine[1:])  # the first line is a free comment
    if h_ref[1:] != h_mine[1:]:
        for r, m in zip(h_ref[1:], h_mine[1:]):
            if r != m:
                add("      reference: %r" % r)
                add("      ours   : %r" % m)

    # 2. positions (same order and same ids)
    same_len = len(ref) == len(mine)
    report("same number of atoms", same_len, "(%d vs %d)" % (len(ref), len(mine)))
    if not same_len:
        add("")
        add("Check that every value in the form is exactly what was used to make the reference file.")
        return out, False
    diff = np.abs(ref[:, 2:5] - mine[:, 2:5]).max()
    report("atom positions identical, same order",
           diff < 1e-5 and (ref[:, 0] == mine[:, 0]).all(), "(max difference %.2e A)" % diff)

    # 3. per-bin type-2 counts, measured on the nanoHUB file
    t2_ref = ref[:, 1] == 2
    t2_me = mine[:, 1] == 2
    add("")
    z_worst = 0.0
    exact_ok = True
    if not res["axial"]:
        # Use the atoms each ring really contains (not a re-measured distance range): atoms that
        # sit EXACTLY on a ring edge (for example the corner edges of a cuboid) would otherwise
        # be counted in two neighbouring rings.
        members = res["members"]
        times = np.zeros(len(mine), dtype=int)
        for m in members:
            times[m] += 1
        on_edge = times > 1
        add("%4s %20s %8s %10s %10s %10s" % ("bin", "range (A)", "atoms", "requested", "reference", "ours"))
        for bi, m in zip(res["bin_info"], members):
            c_ref, c_me = int(t2_ref[m].sum()), int(t2_me[m].sum())
            slack = int(on_edge[m].sum())  # atoms shared with a neighbour ring may add to the count
            req = bi["n_type2_requested"]
            exact_ok &= (req <= c_ref <= req + slack) and (req <= c_me <= req + slack)
            add("%4d %9.3f..%9.3f %8d %10d %10d %10d" % (bi["bin"], bi["lo"], bi["hi"], bi["n_atoms"],
                                                         req, c_ref, c_me))
        if on_edge.any():
            add("(%d atoms sit exactly on a ring edge and belong to two rings.)" % int(on_edge.sum()))
    else:
        for bi in res["bin_info"]:
            v = res["coord"]
            strict = (v > bi["lo"]) & (v < bi["hi"])
            n_s = int(strict.sum())
            p_exp = bi["n_type2_requested"] / bi["n_atoms"] if bi["n_atoms"] else 0.0
            f_ref = t2_ref[strict].mean() if n_s else float("nan")
            f_me = t2_me[strict].mean() if n_s else float("nan")
            sd = math.sqrt(max(p_exp * (1 - p_exp), 1e-12) / max(n_s, 1))
            z = (f_ref - p_exp) / sd if n_s else 0.0
            z_worst = max(z_worst, abs(z))
            add("%4d %9.3f..%9.3f %8d   expected %.3f   reference %.3f   ours %.3f   z=%+.1f"
                % (bi["bin"], bi["lo"], bi["hi"], bi["n_atoms"], p_exp, f_ref, f_me, z))
    if not res["axial"]:
        report("radial: type-2 count per ring identical", exact_ok)
    else:
        report("axial: reference fractions agree with the rule (|z| < 4)", z_worst < 4,
               "(worst |z| = %.1f)" % z_worst)
    # 4. are the type-2 atoms mixed randomly inside each bin? (data-only stand-in for a simulation)
    groups = np.full(len(mine), -1, dtype=int)
    for j in reversed(range(len(res["members"]))):
        groups[res["members"][j]] = j
    groups[groups < 0] = len(res["members"])
    if t2_ref.any() and t2_me.any() and not t2_ref.all():
        pr = neighbour_pairs(mine[:, 2:5], res["lat_p"], res["lattice"])
        z_ref = mixing_zscore(ref[:, 1].astype(int), groups, pr)
        z_me = mixing_zscore(mine[:, 1].astype(int), groups, pr)
        add("")
        report("type-2 atoms are mixed randomly inside each bin (neighbour test)",
               abs(z_ref) < 4 and abs(z_me) < 4, "(z: reference %+.1f, ours %+.1f; |z| < 4 is fine)" % (z_ref, z_me))
    add("")
    if state["ok"]:
        add("OVERALL: ALL CHECKS PASSED", "pass")
    else:
        add("OVERALL: SOMETHING DIFFERS - see FAIL lines above", "fail")
        add("Check that every value in the form is exactly what was used to make the reference file.")
    return out, state["ok"]


def cli_validate(argv=None):
    ap = argparse.ArgumentParser(prog="FGM_NW_Builder.py validate")
    add_arguments(ap)
    ap.add_argument("--ref", required=True, help="reference data file made by another program or version")
    a = ap.parse_args(argv)
    params = params_from_args(a)
    params.update(mode="nanohub", vacuum_xy=0.0, vacuum_z=0.0)  # nanoHUB has neither
    res = build_structure(**params)
    lines, _ = run_checks(res, a.ref)
    for text, _tag in lines:
        print(text)


# =====================================================================================
# PART 3 - the window
# =====================================================================================
UI_SHAPES = {"Cylindrical": "cylindrical", "Cuboid": "cuboid", "Sphere": "sphere"}
UI_DIRS = {"Radial": "radial", "Axial": "axial"}
UI_GRADS = {"P-FGM": "power", "E-FGM": "exp", "S-FGM": "sigmoid", "Custom": "custom"}
UI_LATTS = {"FCC": "fcc", "BCC": "bcc"}
UI_MODES = {"Standard (validated)": "nanohub", "Corrected bins": "fixed"}
UI_SAMPLES = {"Bin centre": "center", "First bin 0 / last bin 1": "edge"}
UI_FRACTIONS = {"Mass fraction": "mass", "Atom fraction": "atom"}
VIEWS = ["LAMMPS DATA FILE", "Bin table", "Grading curves", "Comparison with reference file"]
PREVIEW_ATOMS = 5000
V_DATA, V_TABLE, V_CURVES, V_COMPARE = VIEWS
VIEWS_LAYERED = ["LAMMPS DATA FILE", "Layer & cell table", "Composition along the wire", "Side view (slice)"]
V_TABLE_L, V_PROFILE, V_SIDE = VIEWS_LAYERED[1:]
VIEWS_COMBINED = ["LAMMPS DATA FILE", "Cell table (slice x ring)", "Composition along the wire", "Composition core to surface", "Side view (slice)"]
V_TABLE_C, V_RPROFILE = VIEWS_COMBINED[1], VIEWS_COMBINED[3]
CANVAS_VIEWS = (V_CURVES, V_PROFILE, V_SIDE, V_RPROFILE)
UI_KIND = {"Layers (stacked along the wire)": "layered", "Combined: graded along + across": "combined",
           "Crystal layers: FCC(111) / HCP(0001) mix": "crystal"}
UI_CRYSTAL = {"FCC (111)": "fcc111", "HCP (0001)": "hcp0001"}
UI_INPLANE = {"Mean of the layers": "mean", "Match first layer": "first", "One value": "value"}
UI_SHAPES_L = {"Cylindrical": "cylindrical", "Cuboid": "cuboid"}
UI_BASIS = {"Atom %": "atom", "Mass %": "mass"}
UI_LATMODE = {"Mean of elements": "mean", "One value": "one"}
UI_FUNCS = {"P-FGM": "power", "E-FGM": "exp", "S-FGM": "sigmoid", "Custom": "custom"}
UI_ZSP = {"Same as in-plane": "same", "Follows composition (Vegard)": "vegard"}
PALETTE = ["#1f5fbf", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
# mass (amu), lattice parameter in nm for FCC, for BCC (None = not a stable FCC/BCC metal: type your own value)
ELEMENT_DATA = {
    "Ni": (58.693, 0.3524, None), "Co": (58.933, 0.3544, None), "Cu": (63.546, 0.3615, None),
    "Al": (26.982, 0.4050, None), "Ag": (107.868, 0.4086, None), "Au": (196.967, 0.4078, None),
    "Pd": (106.42, 0.3890, None), "Pt": (195.084, 0.3924, None), "Fe": (55.845, None, 0.2866),
    "Cr": (51.996, None, 0.2884), "Mo": (95.95, None, 0.3147), "W": (183.84, None, 0.3165),
    "Ta": (180.948, None, 0.3303), "Nb": (92.906, None, 0.3300), "V": (50.942, None, 0.3024),
    "Mn": (54.938, None, None), "Ti": (47.867, None, None)}


class ScrollFrame(_TTK_FRAME):
    """A frame with a vertical scroll bar (the layered form is taller than most screens)."""

    def __init__(self, parent, width=560):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, highlightthickness=0, width=width)
        vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(
            scrollregion=self.canvas.bbox("all"), width=max(width, self.inner.winfo_reqwidth())))
        self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.configure(yscrollcommand=vsb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.bind_all("<MouseWheel>", self._wheel, add="+")

    def _wheel(self, e):
        w = self.winfo_containing(e.x_root, e.y_root)
        while w is not None:
            if w is self:
                self.canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
                return
            w = w.master



class FormError(Exception):
    pass


class App(_TK_ROOT):
    def __init__(self):
        super().__init__()
        self.title("FGM Nanostructure Builder %s" % VERSION)
        self.geometry("1360x800")
        self.minsize(1100, 680)
        self.res = None
        self.params = None
        self.curves = None
        self.compare_lines = None
        self.tab = "classic"
        self.states = {"classic": (None, None, None, None), "layered": (None, None, None, None), "combined": (None, None, None, None), "crystal": (None, None, None, None)}
        self.v = {}  # form variables
        self.w = {}  # form widgets
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        menubar = tk.Menu(self)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Manual", command=self.on_manual, accelerator="F1")
        self.bind("<F1>", self.on_manual)
        menubar.add_cascade(label="Help", menu=helpmenu)
        self.config(menu=menubar)
        self._build_tabs()
        self._build_right()
        self.nb.bind("<<NotebookTabChanged>>", self._on_tab)
        self.status = tk.StringVar(value="Fill in the form and press Build.")
        ttk.Label(self, textvariable=self.status, relief="sunken", anchor="w", padding=(8, 3)).grid(
            row=1, column=0, columnspan=2, sticky="ew")
        self._update_state()
        self.show_view(V_DATA)

    # ------------------------------------------------------------------ form
    def _build_left(self, left):
        ttk.Label(left, text="FGM Nanostructure Builder", font=("Segoe UI", 17, "bold")).grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 10))
        self._row = 1
        self._left = left

        def combo(key, label, options, default, parent=left):
            ttk.Label(parent, text=label).grid(row=self._row, column=0, sticky="e", padx=(0, 8), pady=2)
            var = tk.StringVar(value=default)
            cb = ttk.Combobox(parent, textvariable=var, values=list(options), state="readonly", width=24)
            cb.grid(row=self._row, column=1, sticky="w", pady=2)
            cb.bind("<<ComboboxSelected>>", lambda e: self._update_state())
            self.v[key], self.w[key] = var, cb
            self._row += 1

        def entry(key, label, default, parent=left):
            ttk.Label(parent, text=label).grid(row=self._row, column=0, sticky="e", padx=(0, 8), pady=2)
            var = tk.StringVar(value=default)
            en = ttk.Entry(parent, textvariable=var, width=26)
            en.grid(row=self._row, column=1, sticky="w", pady=2)
            self.v[key], self.w[key] = var, en
            self._row += 1

        combo("shape", "Shape", UI_SHAPES, "Cylindrical")
        entry("xdim", "X-dimension (nm)", "2")
        entry("ydim", "Y-dimension (nm)", "2")
        entry("radius", "Radius, R (nm)", "2")
        entry("length", "Length, L (z) (nm)", "20")
        combo("direction", "Grading Direction", UI_DIRS, "Radial")
        combo("grading", "Grading Function", UI_GRADS, "P-FGM")
        entry("custom", "Custom function, f(x) =", "")
        entry("p", "p", "1")
        entry("mass1", "Atomic Mass 1 (amu)", "63.55")
        entry("mass2", "Atomic Mass 2 (amu)", "58.69")
        combo("lattice", "Lattice Type", UI_LATTS, "FCC")
        entry("a1", "Lattice Parameter 1 (nm)", "0.36")
        entry("a2", "Lattice Parameter 2 (nm)", "0.36")
        entry("bins", "Number of bins", "10")

        adv = ttk.LabelFrame(left, text="Advanced (our own options)", padding=(8, 4))
        adv.grid(row=self._row, column=0, columnspan=2, sticky="ew", pady=(12, 0))
        self._row = 0
        combo("mode", "Mode", UI_MODES, "Standard (validated)", adv)
        combo("sample", "Grading sampled at", UI_SAMPLES, "Bin centre", adv)
        combo("fraction", "Grading value is a", UI_FRACTIONS, "Mass fraction", adv)
        entry("vac_xy", "Vacuum x,y (A each side)", "0", adv)
        entry("vac_z", "Vacuum z (A each end)", "0", adv)
        entry("seed", "Random seed (blank = random)", "", adv)

    def _set(self, key, on, readonly=False):
        w = self.w[key]
        if on:
            w.state(["!disabled"])
            if readonly:
                w.state(["readonly"])
        else:
            w.state(["disabled"])

    def _update_state(self):
        shape = UI_SHAPES[self.v["shape"].get()]
        self._set("xdim", shape == "cuboid")
        self._set("ydim", shape == "cuboid")
        self._set("radius", shape in ("cylindrical", "sphere"))
        self._set("length", shape in ("cylindrical", "cuboid"))
        custom = UI_GRADS[self.v["grading"].get()] == "custom"
        self._set("custom", custom)
        self._set("p", not custom)
        fixed = UI_MODES[self.v["mode"].get()] == "fixed"
        self._set("sample", fixed, readonly=True)
        self._set("fraction", fixed, readonly=True)

    def _num(self, key, label, positive=False, minimum=None, integer=False):
        raw = self.v[key].get().strip()
        try:
            x = float(raw)
        except ValueError:
            raise FormError("%s must be a number (you typed %r)." % (label, raw))
        if positive and x <= 0:
            raise FormError("%s must be larger than 0." % label)
        if minimum is not None and x < minimum:
            raise FormError("%s must be at least %g." % (label, minimum))
        if integer:
            if x != int(x):
                raise FormError("%s must be a whole number." % label)
            return int(x)
        return x

    def _collect(self):
        shape = UI_SHAPES[self.v["shape"].get()]
        grading = UI_GRADS[self.v["grading"].get()]
        mode = UI_MODES[self.v["mode"].get()]
        p = dict(shape=shape, direction=UI_DIRS[self.v["direction"].get()], grading=grading,
                 lattice=UI_LATTS[self.v["lattice"].get()], mode=mode,
                 sample=UI_SAMPLES[self.v["sample"].get()], fraction=UI_FRACTIONS[self.v["fraction"].get()])
        p["radius_nm"] = self._num("radius", "Radius", positive=True) if shape != "cuboid" else 2.0
        p["length_nm"] = self._num("length", "Length", positive=True) if shape != "sphere" else 20.0
        p["xdim_nm"] = self._num("xdim", "X-dimension", positive=True) if shape == "cuboid" else 2.0
        p["ydim_nm"] = self._num("ydim", "Y-dimension", positive=True) if shape == "cuboid" else 2.0
        p["p"] = self._num("p", "p", minimum=0) if grading != "custom" else 1.0
        p["custom"] = self.v["custom"].get().strip() if grading == "custom" else ""
        if grading == "custom":
            if not p["custom"]:
                raise FormError("Type a custom function, for example x^2.")
            try:
                grading_fraction(0.5, "custom", 1.0, p["custom"])
            except Exception as e:
                raise FormError("The custom function could not be evaluated: %s" % e)
        p["mass1"] = self._num("mass1", "Atomic Mass 1", positive=True)
        p["mass2"] = self._num("mass2", "Atomic Mass 2", positive=True)
        p["a1_nm"] = self._num("a1", "Lattice Parameter 1", positive=True)
        p["a2_nm"] = self._num("a2", "Lattice Parameter 2", positive=True)
        p["bins"] = self._num("bins", "Number of bins", minimum=2, integer=True)
        p["vacuum_xy"] = self._num("vac_xy", "Vacuum x,y", minimum=0)
        p["vacuum_z"] = self._num("vac_z", "Vacuum z", minimum=0)
        seed = self.v["seed"].get().strip()
        if seed:
            p["seed"] = self._num("seed", "Random seed", integer=True)
        else:
            p["seed"] = None
        return p

    # ---------------------------------------------------------------- result
    def _build_right(self):
        right = ttk.Frame(self, padding=(0, 10, 14, 0))
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(2, weight=1)
        bar = ttk.Frame(right)
        bar.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Button(bar, text="Build", command=self.on_build).pack(side="left", padx=(0, 8))
        self.btn_save = ttk.Button(bar, text="Save data file...", command=self.on_save)
        self.btn_save.pack(side="left", padx=(0, 8))
        self.btn_save.state(["disabled"])
        ttk.Button(bar, text="Compare with reference file...", command=self.on_compare).pack(side="left")

        pick = ttk.Frame(right)
        pick.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        ttk.Label(pick, text="Result:").pack(side="left", padx=(0, 6))
        self.view = tk.StringVar(value=VIEWS[0])
        cb = ttk.Combobox(pick, textvariable=self.view, values=VIEWS, state="readonly", width=36)
        cb.pack(side="left")
        self.view_box = cb
        cb.bind("<<ComboboxSelected>>", lambda e: self.show_view(self.view.get()))

        holder = ttk.Frame(right)
        holder.grid(row=2, column=0, sticky="nsew")
        holder.columnconfigure(0, weight=1)
        holder.rowconfigure(0, weight=1)
        self.text_frame = ttk.Frame(holder)
        self.text_frame.grid(row=0, column=0, sticky="nsew")
        self.text_frame.columnconfigure(0, weight=1)
        self.text_frame.rowconfigure(0, weight=1)
        mono = ("Consolas", 10) if sys.platform.startswith("win") else ("Courier", 10)
        self.text = tk.Text(self.text_frame, wrap="none", font=mono, state="disabled", undo=False)
        self.text.grid(row=0, column=0, sticky="nsew")
        sy = ttk.Scrollbar(self.text_frame, orient="vertical", command=self.text.yview)
        sx = ttk.Scrollbar(self.text_frame, orient="horizontal", command=self.text.xview)
        self.text.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.grid(row=0, column=1, sticky="ns")
        sx.grid(row=1, column=0, sticky="ew")
        self.text.tag_configure("pass", foreground="#117a1f")
        self.text.tag_configure("fail", foreground="#b3261e", font=mono + ("bold",))

        self.canvas = tk.Canvas(holder, bg="white", highlightthickness=1, highlightbackground="#999")
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.canvas.bind("<Configure>", lambda e: self._draw_canvas(self.view.get()) if self.view.get() in CANVAS_VIEWS else None)

    def set_text(self, lines, empty_msg=""):
        """lines: list of text or (text, tag)."""
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        if not lines:
            self.text.insert("end", empty_msg)
        for item in lines:
            txt, tag = item if isinstance(item, tuple) else (item, None)
            self.text.insert("end", txt + "\n", tag)
        self.text.configure(state="disabled")

    def show_view(self, name):
        self.view.set(name)
        if name in CANVAS_VIEWS:
            tk.Misc.tkraise(self.canvas)  # (Canvas has its own tkraise for drawn items)
            self.update_idletasks()
            self._draw_canvas(name)
            return
        self.text_frame.tkraise()
        if name == V_DATA:
            if self.res is None:
                self.set_text([], "Nothing built yet. Fill in the form and press Build.")
            else:
                n = len(self.res["atoms"])
                txt = lammps_data_text(self.res, PREVIEW_ATOMS)
                if n > PREVIEW_ATOMS:
                    txt += "\n... preview shows the first %d of %d atoms. The saved file contains all atoms." % (
                        PREVIEW_ATOMS, n)
                self.set_text(txt.split("\n"))
        elif name in (V_TABLE, V_TABLE_L, V_TABLE_C):
            if self.res is None:
                self.set_text([], "Nothing built yet. Fill in the form and press Build.")
            else:
                self.set_text([(l, "pass" if l.startswith("[PASS]") else "fail" if l.startswith("[FAIL]") else None)
                               for l in summary_lines(self.res)])
        else:
            if self.compare_lines is None:
                self.set_text([], "Press \"Compare with reference file...\" and pick a data file made by another program or version.")
            else:
                self.set_text(self.compare_lines)

    # ------------------------------------------------------------ plotting
    def _make_curves(self):
        p, res = self.params, self.res
        xs = [i / 1000 for i in range(1001)]
        ys = [100 * grading_fraction(x, p["grading"], p["p"], p["custom"]) for x in xs]
        yy = res["yy"]
        if res["mode"] == "nanohub":  # drawn the way nanoHUB draws its second curve
            xm = [i / (len(yy) - 1) for i in range(len(yy))]
            ym = [100 * y for y in yy]
        else:
            xm = [bi["b"] for bi in res["bin_info"]]
            ym = [100 * y for y in yy[1:]]
        pts = [(x, y) for x, y in zip(xm, ym) if y == y]
        return [("Desired grading function", xs, ys, "#1f5fbf", "line"),
                ("Built structure (per bin)", [a for a, _ in pts], [b for _, b in pts], "#d62728", "dots")]

    def _draw_curves(self):
        c = self.canvas
        c.delete("all")
        W, H = c.winfo_width(), c.winfo_height()
        if W < 80 or H < 80:
            return
        if self.curves is None:
            c.create_text(W / 2, H / 2, text="Nothing built yet. Fill in the form and press Build.", fill="#555")
            return
        L, R, T, B = 70, 25, 50, 60
        pw, ph = W - L - R, H - T - B
        px = lambda x: L + x * pw
        py = lambda y: T + ph - y / 100 * ph
        c.create_rectangle(L, T, L + pw, T + ph, outline="#444")
        for i in range(6):
            x, y = i / 5, i * 20
            c.create_line(px(x), T + ph, px(x), T + ph + 5, fill="#444")
            c.create_text(px(x), T + ph + 16, text="%.1f" % x)
            c.create_line(L - 5, py(y), L, py(y), fill="#444")
            c.create_text(L - 28, py(y), text="%d" % y)
            if 0 < i < 5:
                c.create_line(px(x), T, px(x), T + ph, fill="#e4e4e4")
                c.create_line(L, py(y), L + pw, py(y), fill="#e4e4e4")
        c.create_text(L + pw / 2, H - 14, text="Normalised position (0 = start, 1 = end)")
        c.create_text(18, T + ph / 2, text="Material 2 (wt %)", angle=90)
        x0 = L
        for name, xs, ys, color, style in self.curves:
            if style == "line":
                c.create_line([v for x, y in zip(xs, ys) for v in (px(x), py(y))], fill=color, width=2)
            else:
                pts = [(px(x), py(y)) for x, y in zip(xs, ys)]
                if len(pts) > 1:
                    c.create_line([v for pt in pts for v in pt], fill=color, width=1, dash=(3, 3))
                for X, Y in pts:
                    c.create_oval(X - 4, Y - 4, X + 4, Y + 4, fill=color, outline=color)
            c.create_line(x0, 22, x0 + 24, 22, fill=color, width=2)
            c.create_text(x0 + 30, 22, text=name, anchor="w")
            x0 += 230

    # -------------------------------------------------------------- actions
    # ------------------------------------------------------- layered tab
    def _build_tabs(self):
        self.nb = ttk.Notebook(self)
        self.nb.grid(row=0, column=0, sticky="ns", padx=(10, 0), pady=(8, 0))
        tab_c = ttk.Frame(self.nb, padding=(14, 10))
        self.scroll_l = ScrollFrame(self.nb, width=560)
        self.nb.add(tab_c, text="FGM Nanostructure Builder")
        self.nb.add(self.scroll_l, text="Multilayer & Multi-element Nanostructure Builder")
        self._build_left(tab_c)
        self._build_layered(self.scroll_l.inner)

    def _tab_for_ui(self):
        if self.nb.index(self.nb.select()) == 0:
            return "classic"
        return UI_KIND[self.vl["kind"].get()]

    def _on_tab(self, event=None):
        self._switch_state(self._tab_for_ui())

    def _switch_state(self, new):
        if new == self.tab:
            return
        self.states[self.tab] = (self.res, self.params, self.curves, self.compare_lines)
        self.res, self.params, self.curves, self.compare_lines = self.states[new]
        self.tab = new
        views = {"classic": VIEWS, "layered": VIEWS_LAYERED, "combined": VIEWS_COMBINED, "crystal": VIEWS_LAYERED}[new]
        self.view_box.configure(values=views)
        self.btn_save.state(["!disabled"] if self.res is not None else ["disabled"])
        self.status.set({"classic": "FGM Nanostructure Builder: fill in the form and press Build.",
                         "layered": "Multilayer: set the elements and layers, then press Build.",
                         "combined": "Combined grading: set elements, corner recipes and functions, then press Build.",
                         "crystal": "Crystal layers: set elements and the crystal of each layer, then press Build."}[new])
        self.show_view(views[0])

    def _build_layered(self, root):
        self.vl, self.wl = {}, {}
        pad = dict(padx=10, pady=(8, 0), sticky="ew")
        ttk.Label(root, text="Multilayer & Multi-element Nanostructure Builder", font=("Segoe UI", 15, "bold")).grid(
            row=0, column=0, sticky="w", padx=10, pady=(10, 0))

        kindf = ttk.Frame(root)
        kindf.grid(row=1, column=0, sticky="w", padx=10, pady=(6, 0))
        ttk.Label(kindf, text="Structure type").pack(side="left", padx=(0, 8))
        self.vl["kind"] = tk.StringVar(value="Layers (stacked along the wire)")
        kc = ttk.Combobox(kindf, textvariable=self.vl["kind"], values=list(UI_KIND), state="readonly", width=34)
        kc.pack(side="left")
        kc.bind("<<ComboboxSelected>>", self._on_kind)
        wire = ttk.LabelFrame(root, text="Wire and settings", padding=(8, 4))
        wire.grid(row=2, column=0, **pad)

        def field(row, col, key, label, kind, default, options=None, width=17):
            ttk.Label(wire, text=label).grid(row=row, column=col * 2, sticky="e", padx=(0, 6), pady=2)
            var = tk.StringVar(value=default)
            if kind == "combo":
                w = ttk.Combobox(wire, textvariable=var, values=list(options), state="readonly", width=width + 6)
                w.bind("<<ComboboxSelected>>", lambda e: self._update_layered_state())
            else:
                w = ttk.Entry(wire, textvariable=var, width=width)
            w.grid(row=row, column=col * 2 + 1, sticky="w", pady=2, padx=(0, 10))
            self.vl[key], self.wl[key] = var, w

        field(0, 0, "shape", "Shape", "combo", "Cylindrical", UI_SHAPES_L, width=10)
        field(0, 1, "lattice", "Lattice type", "combo", "FCC", UI_LATTS, width=16)
        field(1, 0, "radius", "Radius (nm)", "entry", "1.5", width=16)
        field(1, 1, "basis", "Composition is in", "combo", "Atom %", UI_BASIS, width=16)
        field(2, 0, "xdim", "X-dimension (nm)", "entry", "2", width=16)
        field(2, 1, "sample", "Grading sampled at", "combo", "Bin centre", UI_SAMPLES, width=16)
        field(3, 0, "ydim", "Y-dimension (nm)", "entry", "2", width=16)
        field(3, 1, "latmode", "Lattice parameter", "combo", "Mean of elements", UI_LATMODE, width=16)
        field(4, 0, "vac_xy", "Vacuum x,y (A)", "entry", "0", width=16)
        field(4, 1, "one_a", "One value (nm)", "entry", "0.354", width=22)
        field(5, 0, "vac_z", "Vacuum z (A)", "entry", "0", width=16)
        field(5, 1, "seed", "Random seed (optional)", "entry", "", width=22)
        field(6, 0, "zspace", "Spacing along the wire", "combo", "Same as in-plane", UI_ZSP, width=16)
        ttk.Label(wire, foreground="#444", wraplength=330, justify="left",
                  text="Vegard: each layer (or slice) gets the lattice spacing of its own composition, using the element constants below; the x,y spacing stays common, so layers are strained to fit.").grid(row=6, column=2, columnspan=2, sticky="w")

        el = ttk.LabelFrame(root, text="Elements (2 to 5): this order = composition order = LAMMPS atom types 1, 2, 3...",
                            padding=(8, 4))
        el.grid(row=3, column=0, **pad)
        ttk.Label(el, text="Number of elements").grid(row=0, column=0, sticky="e", padx=(0, 6))
        self.vl["n_el"] = tk.StringVar(value="2")
        sp = ttk.Spinbox(el, from_=2, to=5, width=4, textvariable=self.vl["n_el"], state="readonly",
                         command=self._on_n_elements)
        sp.grid(row=0, column=1, sticky="w")
        for c, t in enumerate(("Symbol", "Mass (amu)", "Lattice parameter a (nm)")):
            ttk.Label(el, text=t, font=("Segoe UI", 9, "bold")).grid(row=1, column=c, sticky="w", padx=(0, 10), pady=(6, 0))
        defaults = [("Ni", "58.693", "0.3524"), ("Co", "58.933", "0.3544")] + [("", "", "")] * 3
        self.el_rows = []
        for i, (sym, mass, a) in enumerate(defaults):
            r = {"sym": tk.StringVar(value=sym), "mass": tk.StringVar(value=mass), "a": tk.StringVar(value=a)}
            r["wsym"] = ttk.Combobox(el, textvariable=r["sym"], values=sorted(ELEMENT_DATA), width=8)
            r["wmass"] = ttk.Entry(el, textvariable=r["mass"], width=12)
            r["wa"] = ttk.Entry(el, textvariable=r["a"], width=14)
            for c, wname in enumerate(("wsym", "wmass", "wa")):
                r[wname].grid(row=2 + i, column=c, sticky="w", padx=(0, 10), pady=1)
            r["wsym"].bind("<<ComboboxSelected>>", lambda e, k=i: self._fill_element(k))
            r["wsym"].bind("<FocusOut>", lambda e, k=i: self._fill_element(k, only_if_empty=True))
            r["sym"].trace_add("write", lambda *a: self._refresh_headers())
            self.el_rows.append(r)

        lay = ttk.LabelFrame(root, text="Layers (bottom to top)", padding=(8, 4))
        lay.grid(row=4, column=0, **pad)
        self.frame_layers = lay
        self.layers_frame = ttk.Frame(lay)
        self.layers_frame.grid(row=0, column=0, sticky="w")
        self.layer_head = []
        for c in range(9):
            h = ttk.Label(self.layers_frame, text="", font=("Segoe UI", 9, "bold"))
            h.grid(row=0, column=c, sticky="w", padx=(0, 6))
            self.layer_head.append(h)
        self.layer_rows = []
        bar = ttk.Frame(lay)
        bar.grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Button(bar, text="Add layer", command=lambda: self._add_layer()).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="Remove last", command=self._remove_layer).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="Ni/Co example", command=self._example_binary).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="Co/Cr/Ni example", command=self._example_ternary).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="Ni/Cu example", command=self._example_nicu).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="Solid shell example", command=self._example_shell).pack(side="left")
        ttk.Label(lay, justify="left", wraplength=520, foreground="#444",
                  text="Composition = percent of each element in the order above, for example 100,0 or 40,30,30 (must add to 100). "
                       "Leave 'End' empty for a uniform layer. A graded layer moves from Start to End with the chosen function. "
                       "Direction: Axial = along the wire, Radial = from the core to the surface. "
                       "A graded interface is simply a thin graded layer between two others. "
                       "Layer thicknesses are rounded to whole atomic layers. a_z: optional lattice constant of that layer along the wire (needs Vegard mode).").grid(row=2, column=0, sticky="w", pady=(6, 0))
        # ---- combined grading box (shown instead of the layers box)
        comb = ttk.LabelFrame(root, text="Combined grading: along the wire AND from the core to the surface", padding=(8, 4))
        comb.grid(row=4, column=0, **pad)
        comb.grid_remove()
        self.frame_comb = comb
        ttk.Label(comb, text="Length (nm)").grid(row=0, column=0, sticky="e", padx=(0, 6), pady=2)
        self.vl["c_len"] = tk.StringVar(value="10")
        ttk.Entry(comb, textvariable=self.vl["c_len"], width=10).grid(row=0, column=1, sticky="w", pady=2)
        self.corner_labels = []
        for r_, (key, base, default) in enumerate((("c_cb", "Core, bottom", "100,0"), ("c_sb", "Surface, bottom", "70,30"),
                                                  ("c_ct", "Core, top", "30,70"), ("c_st", "Surface, top", "0,100")), start=1):
            lbl = ttk.Label(comb, text=base)
            lbl.grid(row=r_, column=0, sticky="e", padx=(0, 6), pady=2)
            self.vl[key] = tk.StringVar(value=default)
            ttk.Entry(comb, textvariable=self.vl[key], width=24).grid(row=r_, column=1, columnspan=3, sticky="w", pady=2)
            self.corner_labels.append((lbl, base))
        for r_, (key, text, fdef, pdef, bdef, btxt) in enumerate((("z", "Along the wire (bottom to top)", "S-FGM", "2", "10", "slices"),
                                                                 ("r", "Core to surface", "P-FGM", "1", "5", "rings")), start=5):
            ttk.Label(comb, text=text).grid(row=r_, column=0, sticky="e", padx=(0, 6), pady=(8 if r_ == 5 else 2, 2))
            sub = ttk.Frame(comb)
            sub.grid(row=r_, column=1, columnspan=3, sticky="w", pady=(8 if r_ == 5 else 2, 2))
            self.vl["c_f" + key] = tk.StringVar(value=fdef)
            self.vl["c_p" + key] = tk.StringVar(value=pdef)
            self.vl["c_b" + key] = tk.StringVar(value=bdef)
            ttk.Combobox(sub, textvariable=self.vl["c_f" + key], values=list(UI_FUNCS), state="readonly", width=8).pack(side="left")
            fcb = sub.winfo_children()[-1]
            fcb.bind("<<ComboboxSelected>>", lambda e, k=key: self._func_changed(self.vl["c_f" + k], self.vl["c_p" + k]))
            ttk.Label(sub, text="  p or f(x)").pack(side="left")
            ttk.Entry(sub, textvariable=self.vl["c_p" + key], width=12).pack(side="left", padx=(2, 8))
            ttk.Entry(sub, textvariable=self.vl["c_b" + key], width=5).pack(side="left", padx=(0, 2))
            ttk.Label(sub, text=btxt).pack(side="left")
        cbar = ttk.Frame(comb)
        cbar.grid(row=7, column=0, columnspan=4, sticky="w", pady=(8, 0))
        ttk.Button(cbar, text="Example: Ni/Co", command=self._example_comb_binary).pack(side="left", padx=(0, 6))
        ttk.Button(cbar, text="Solid regions", command=self._example_comb_solid).pack(side="left", padx=(0, 6))
        ttk.Button(cbar, text="Example: Co/Cr/Ni", command=self._example_comb_ternary).pack(side="left")
        ttk.Label(comb, justify="left", wraplength=520, foreground="#444",
                  text="Each corner is a recipe: percent of every element, in the order above, adding to 100. "
                       "The two functions decide how fast the mix moves from one corner to the next, one along the wire and one "
                       "across it, so any mix between the corners is always valid. Equal corners give a wire graded in one "
                       "direction only. Slices cut the wire along its length, rings cut it from the core to the surface.").grid(
                           row=8, column=0, columnspan=4, sticky="w", pady=(6, 0))

        # ---- crystal layers box (FCC (111) / HCP (0001) mix)
        cry = ttk.LabelFrame(root, text="Crystal layers: FCC (111) and HCP (0001) layers joined along the wire", padding=(8, 4))
        cry.grid(row=4, column=0, **pad)
        cry.grid_remove()
        self.frame_crystal = cry
        top = ttk.Frame(cry)
        top.grid(row=0, column=0, sticky="w")
        ttk.Label(top, text="In-plane spacing").pack(side="left", padx=(0, 6))
        self.vl["cr_inplane"] = tk.StringVar(value="Mean of the layers")
        ttk.Combobox(top, textvariable=self.vl["cr_inplane"], values=list(UI_INPLANE), state="readonly", width=18).pack(side="left")
        ttk.Label(top, text="   one value: nearest-neighbour distance (nm)").pack(side="left")
        self.vl["cr_inplane_nm"] = tk.StringVar(value="")
        ttk.Entry(top, textvariable=self.vl["cr_inplane_nm"], width=8).pack(side="left", padx=(6, 0))
        self.crystal_frame = ttk.Frame(cry)
        self.crystal_frame.grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.crystal_head = []
        for c in range(6):
            h = ttk.Label(self.crystal_frame, text="", font=("Segoe UI", 9, "bold"))
            h.grid(row=0, column=c, sticky="w", padx=(0, 6))
            self.crystal_head.append(h)
        self.crystal_rows = []
        cbar2 = ttk.Frame(cry)
        cbar2.grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Button(cbar2, text="Add layer", command=lambda: self._add_crystal()).pack(side="left", padx=(0, 6))
        ttk.Button(cbar2, text="Remove last", command=self._remove_crystal).pack(side="left", padx=(0, 6))
        ttk.Button(cbar2, text="Ti/Al example", command=self._example_crystal_tial).pack(side="left", padx=(0, 6))
        ttk.Button(cbar2, text="Ni/Cu (111) example", command=self._example_crystal_nicu).pack(side="left")
        ttk.Label(cry, justify="left", wraplength=520, foreground="#444",
                  text="Every layer is its own crystal with its close-packed planes across the wire (FCC along [111], HCP along [0001]) "
                       "and close-packed rows parallel in all layers. a (and c for HCP) are the layer's own lattice constants; "
                       "blank c means the ideal ratio. All layers share one in-plane spacing, so they are strained to fit (the table "
                       "reports the strain); each layer keeps its own spacing between planes. Composition is one recipe per layer "
                       "(for example 100,0 for pure Ti). Lattice type, lattice parameter and spacing settings above are not used "
                       "here. This builds a starting structure: misfit relaxation happens in your LAMMPS run.").grid(
                           row=3, column=0, sticky="w", pady=(6, 0))

        self._example_binary()
        self._example_crystal_tial_rows()
        self._on_n_elements()
        self._update_layered_state()

    def _n_el(self):
        return int(self.vl["n_el"].get())

    def _symbols(self):
        return [self.el_rows[i]["sym"].get().strip() or "E%d" % (i + 1) for i in range(self._n_el())]

    def _refresh_headers(self):
        sy = ",".join(self._symbols())
        for lbl, base in getattr(self, "corner_labels", []):
            lbl.configure(text="%s %% (%s)" % (base, sy))
        for h, t in zip(getattr(self, "crystal_head", []), ("Layer", "Crystal", "a (nm)", "c (nm)", "Thickness (nm)", "Composition %% (%s)" % sy)):
            h.configure(text=t)
        for h, t in zip(self.layer_head, ("Layer", "Thickness (nm)", "Start %% (%s)" % sy, "End %% (%s)" % sy,
                                          "Function", "p or f(x)", "Direction", "Bins", "a_z (nm)")):
            h.configure(text=t)

    def _fill_element(self, i, only_if_empty=False):
        r = self.el_rows[i]
        sym = r["sym"].get().strip()
        if sym not in ELEMENT_DATA:
            return
        mass, a_fcc, a_bcc = ELEMENT_DATA[sym]
        a = a_fcc if UI_LATTS[self.vl["lattice"].get()] == "fcc" else a_bcc
        if not only_if_empty or not r["mass"].get().strip():
            r["mass"].set("%g" % mass)
        if (not only_if_empty or not r["a"].get().strip()) and a is not None:
            r["a"].set("%g" % a)
        elif a is None and not only_if_empty:
            r["a"].set("")
            self.status.set("%s has no standard %s lattice parameter here: type the value from your potential." % (
                sym, self.vl["lattice"].get()))

    def _on_n_elements(self):
        n = self._n_el()
        for i, r in enumerate(self.el_rows):
            for wname in ("wsym", "wmass", "wa"):
                r[wname].state(["!disabled"] if i < n else ["disabled"])
        for row in self.layer_rows:  # keep the composition lists the right length
            for key in ("start", "end"):
                txt = row[key].get().strip()
                if txt:
                    parts = [x for x in txt.replace(";", ",").replace(" ", ",").split(",") if x]
                    parts = (parts + ["0"] * n)[:n]
                    row[key].set(",".join(parts))
        for row in getattr(self, "crystal_rows", []):
            parts = [x for x in row["comp"].get().replace(";", ",").replace(" ", ",").split(",") if x]
            row["comp"].set(",".join((parts + ["0"] * n)[:n]))
        for key in ("c_cb", "c_sb", "c_ct", "c_st"):
            if key in self.vl:
                parts = [x for x in self.vl[key].get().replace(";", ",").replace(" ", ",").split(",") if x]
                self.vl[key].set(",".join((parts + ["0"] * n)[:n]))
        self._refresh_headers()

    def _update_layered_state(self):
        shape = UI_SHAPES_L[self.vl["shape"].get()]
        for key, on in (("radius", shape == "cylindrical"), ("xdim", shape == "cuboid"), ("ydim", shape == "cuboid")):
            self.wl[key].state(["!disabled"] if on else ["disabled"])
        crystal = UI_KIND[self.vl["kind"].get()] == "crystal" if "kind" in self.vl else False
        for key in ("lattice", "latmode", "zspace"):
            w = self.wl[key]
            if crystal:
                w.state(["disabled"])
            else:
                w.state(["!disabled"])
                w.state(["readonly"])
        one = UI_LATMODE[self.vl["latmode"].get()] == "one"
        self.wl["one_a"].state(["!disabled"] if (one and not crystal) else ["disabled"])

    def _add_layer(self, thick="5", start=None, end="", func="P-FGM", p="1", direction="Axial", bins="10", az=""):
        if len(self.layer_rows) >= 12:
            messagebox.showinfo("Layers", "Twelve layers is the limit of this form.")
            return
        n = self._n_el()
        if start is None:
            start = ",".join(["100"] + ["0"] * (n - 1))
        k = len(self.layer_rows) + 1
        row = {"thick": tk.StringVar(value=thick), "start": tk.StringVar(value=start), "end": tk.StringVar(value=end),
               "func": tk.StringVar(value=func), "p": tk.StringVar(value=p), "dir": tk.StringVar(value=direction),
               "bins": tk.StringVar(value=bins), "az": tk.StringVar(value=az)}
        ws = [ttk.Label(self.layers_frame, text=str(k)),
              ttk.Entry(self.layers_frame, textvariable=row["thick"], width=9),
              ttk.Entry(self.layers_frame, textvariable=row["start"], width=14),
              ttk.Entry(self.layers_frame, textvariable=row["end"], width=14),
              ttk.Combobox(self.layers_frame, textvariable=row["func"], values=list(UI_FUNCS), state="readonly", width=7),
              ttk.Entry(self.layers_frame, textvariable=row["p"], width=12),
              ttk.Combobox(self.layers_frame, textvariable=row["dir"], values=list(UI_DIRS)[::-1], state="readonly", width=7),
              ttk.Entry(self.layers_frame, textvariable=row["bins"], width=5),
              ttk.Entry(self.layers_frame, textvariable=row["az"], width=7)]
        for c, w in enumerate(ws):
            w.grid(row=k, column=c, sticky="w", padx=(0, 6), pady=1)
        ws[4].bind("<<ComboboxSelected>>", lambda e, row=row: self._func_changed(row["func"], row["p"]))
        row["widgets"] = ws
        self.layer_rows.append(row)

    def _func_changed(self, func_var, p_var):
        """Switching to Custom puts a starting formula in the p box; switching back puts p = 1."""
        txt = p_var.get().strip()
        try:
            float(txt)
            numeric = True
        except ValueError:
            numeric = False
        if UI_FUNCS[func_var.get()] == "custom" and numeric:
            p_var.set("x>0.7")
        elif UI_FUNCS[func_var.get()] != "custom" and not numeric:
            p_var.set("1")

    def _grading_inputs(self, func_text, box_text, what):
        """(function name, p, formula) from a function selector and its p / f(x) box."""
        kind = UI_FUNCS[func_text]
        if kind == "custom":
            formula = box_text.strip()
            if not formula:
                raise FormError("%s: type the custom formula f(x) in the p / f(x) box (for example x>0.7)." % what)
            try:
                check_grading("custom", 1.0, formula, what)
            except ValueError as e:
                raise FormError(str(e))
            return kind, 1.0, formula
        return kind, self._to_float(box_text, "%s p" % what, minimum=0), ""

    def _remove_layer(self):
        if len(self.layer_rows) > 1:
            for w in self.layer_rows.pop()["widgets"]:
                w.destroy()

    def _clear_layers(self):
        while self.layer_rows:
            for w in self.layer_rows.pop()["widgets"]:
                w.destroy()

    def _set_elements(self, items, lattice="FCC", latmode="Mean of elements", one_a="0.354", zspace="Same as in-plane"):
        self.vl["zspace"].set(zspace)
        self.vl["n_el"].set(str(len(items)))
        self.vl["lattice"].set(lattice)
        self.vl["latmode"].set(latmode)
        self.vl["one_a"].set(one_a)
        for i, r in enumerate(self.el_rows):
            sym = items[i] if i < len(items) else ""
            r["sym"].set(sym)
            if sym:
                self._fill_element(i)
            else:
                r["mass"].set("")
                r["a"].set("")
        self._update_layered_state()

    def _example_binary(self):
        self._set_elements(["Ni", "Co"])
        self._clear_layers()
        self._add_layer("4", "100,0", "")
        self._add_layer("6", "100,0", "0,100", "S-FGM", "2", "Axial", "8")
        self._add_layer("4", "0,100", "")
        self._add_layer("6", "0,100", "100,0", "P-FGM", "1", "Axial", "6")
        self._on_n_elements()

    def _example_nicu(self):
        self._set_elements(["Ni", "Cu"], zspace="Follows composition (Vegard)")
        self._clear_layers()
        self._add_layer("4", "100,0", "")
        self._add_layer("1.5", "100,0", "0,100", "S-FGM", "2", "Axial", "4")
        self._add_layer("4", "0,100", "")
        self._add_layer("1.5", "0,100", "100,0", "S-FGM", "2", "Axial", "4")
        self._on_n_elements()

    def _example_shell(self):
        """A pure Ni wire with a solid Co shell in the middle part: Custom function x>0.7, radial."""
        self._set_elements(["Ni", "Co"])
        self._clear_layers()
        self._add_layer("3", "100,0", "")
        self._add_layer("6", "100,0", "0,100", "Custom", "x>0.7", "Radial", "10")
        self._add_layer("3", "100,0", "")
        self._on_n_elements()

    def _example_ternary(self):
        self._set_elements(["Co", "Cr", "Ni"], latmode="One value", one_a="0.354")
        self._clear_layers()
        self._add_layer("4", "33.3,33.3,33.4", "")
        self._add_layer("6", "50,25,25", "10,20,70", "S-FGM", "2", "Axial", "8")
        self._add_layer("4", "10,20,70", "")
        self._add_layer("5", "10,20,70", "50,25,25", "P-FGM", "1", "Radial", "5")
        self._on_n_elements()

    def _to_float(self, raw, label, positive=False, minimum=None, integer=False):
        raw = str(raw).strip()
        try:
            x = float(raw)
        except ValueError:
            raise FormError("%s must be a number (you typed %r)." % (label, raw))
        if positive and x <= 0:
            raise FormError("%s must be larger than 0." % label)
        if minimum is not None and x < minimum:
            raise FormError("%s must be at least %g." % (label, minimum))
        if integer:
            if x != int(x):
                raise FormError("%s must be a whole number." % label)
            return int(x)
        return x

    def _recipe(self, text, n, what):
        parts = [x for x in text.replace(";", ",").replace(" ", ",").split(",") if x]
        if len(parts) != n:
            raise FormError("%s needs %d numbers (one per element) but has %d." % (what, n, len(parts)))
        vals = [self._to_float(x, what) for x in parts]
        if min(vals) < 0:
            raise FormError("%s cannot have a negative percentage." % what)
        tot = sum(vals)
        if abs(tot - 100.0) > 0.01:
            raise FormError("%s adds up to %g %%; it must add up to 100." % (what, tot))
        return [x / tot for x in vals]

    def _collect_layered(self):
        v = self.vl
        shape = UI_SHAPES_L[v["shape"].get()]
        n = self._n_el()
        one = UI_LATMODE[v["latmode"].get()] == "one"
        one_a = self._to_float(v["one_a"].get(), "Lattice parameter (one value)", positive=True) if one else None
        elements, seen = [], set()
        for i in range(n):
            r = self.el_rows[i]
            sym = r["sym"].get().strip()
            if not sym:
                raise FormError("Element %d needs a symbol." % (i + 1))
            if sym in seen:
                raise FormError("Element symbols must differ (%s appears twice)." % sym)
            seen.add(sym)
            mass = self._to_float(r["mass"].get(), "Mass of %s" % sym, positive=True)
            need_own = (not one) or UI_ZSP[v["zspace"].get()] == "vegard"
            a = self._to_float(r["a"].get(), "Lattice parameter of %s" % sym, positive=True) if need_own else one_a
            elements.append({"symbol": sym, "mass": mass, "a_nm": a})
        layers = []
        for k, row in enumerate(self.layer_rows, 1):
            endtxt = row["end"].get().strip()
            gk, gp, gc = self._grading_inputs(row["func"].get(), row["p"].get(), "Layer %d" % k)
            layers.append(dict(
                thickness_nm=self._to_float(row["thick"].get(), "Layer %d thickness" % k, positive=True),
                start=self._recipe(row["start"].get(), n, "Layer %d start composition" % k),
                end=self._recipe(endtxt, n, "Layer %d end composition" % k) if endtxt else None,
                grading=gk, p=gp, custom=gc,
                direction=UI_DIRS[row["dir"].get()],
                bins=self._to_float(row["bins"].get(), "Layer %d bins" % k, minimum=1, integer=True),
                z_a_nm=self._to_float(row["az"].get(), "Layer %d a_z" % k, positive=True) if row["az"].get().strip() else None))
        if any(ly["z_a_nm"] is not None for ly in layers) and UI_ZSP[v["zspace"].get()] != "vegard":
            raise FormError("A lattice constant a_z for a layer needs \"Spacing along the wire\" set to Follows composition (Vegard).")
        if not layers:
            raise FormError("Add at least one layer.")
        seed = v["seed"].get().strip()
        return dict(shape=shape, radius_nm=self._to_float(v["radius"].get(), "Radius", positive=True) if shape == "cylindrical" else 2.0,
                    xdim_nm=self._to_float(v["xdim"].get(), "X-dimension", positive=True) if shape == "cuboid" else 2.0,
                    ydim_nm=self._to_float(v["ydim"].get(), "Y-dimension", positive=True) if shape == "cuboid" else 2.0,
                    lattice=UI_LATTS[v["lattice"].get()], elements=elements, layers=layers,
                    basis=UI_BASIS[v["basis"].get()], sample=UI_SAMPLES[v["sample"].get()],
                    vacuum_xy=self._to_float(v["vac_xy"].get(), "Vacuum x,y", minimum=0),
                    vacuum_z=self._to_float(v["vac_z"].get(), "Vacuum z", minimum=0),
                    seed=self._to_float(seed, "Random seed", integer=True) if seed else None, lattice_param_nm=one_a,
                    z_spacing=UI_ZSP[v["zspace"].get()])

    def on_build_layered(self):
        try:
            kw = self._collect_layered()
        except FormError as e:
            messagebox.showerror("Please check the form", str(e))
            return
        self._busy(True)
        try:
            t0 = time.time()
            res = layered_build(**kw)
            dt = time.time() - t0
        except Exception as e:
            self._busy(False)
            messagebox.showerror("Build failed", "%s: %s" % (type(e).__name__, e))
            return
        self._busy(False)
        self.res, self.params, self.curves = res, ("layered", kw), None
        self.btn_save.state(["!disabled"])
        extra = ("  " + res["messages"][0]) if res["messages"] else ""
        self.status.set("Built %d atoms in %.1f s: %d layers, %d cells.%s" % (
            len(res["atoms"]), dt, len(res["layers"]), len(res["bin_info"]), extra))
        self.show_view(V_DATA)

    # ------------------------------------------------ combined grading
    def _on_kind(self, event=None):
        kind = UI_KIND[self.vl["kind"].get()]
        for frame in (self.frame_layers, self.frame_comb, self.frame_crystal):
            frame.grid_remove()
        {"layered": self.frame_layers, "combined": self.frame_comb, "crystal": self.frame_crystal}[kind].grid()
        self._update_layered_state()
        self._switch_state(self._tab_for_ui())

    # ------------------------------------------------ crystal layers
    def _add_crystal(self, crystal="HCP (0001)", a="0.2951", c="0.4683", thick="4", comp=None):
        if len(self.crystal_rows) >= 12:
            messagebox.showinfo("Layers", "Twelve layers is the limit of this form.")
            return
        n = self._n_el()
        if comp is None:
            comp = ",".join(["100"] + ["0"] * (n - 1))
        k = len(self.crystal_rows) + 1
        row = {"crystal": tk.StringVar(value=crystal), "a": tk.StringVar(value=a), "c": tk.StringVar(value=c),
               "thick": tk.StringVar(value=thick), "comp": tk.StringVar(value=comp)}
        cb = ttk.Combobox(self.crystal_frame, textvariable=row["crystal"], values=list(UI_CRYSTAL), state="readonly", width=11)
        ec = ttk.Entry(self.crystal_frame, textvariable=row["c"], width=8)
        ws = [ttk.Label(self.crystal_frame, text=str(k)), cb, ttk.Entry(self.crystal_frame, textvariable=row["a"], width=8), ec,
              ttk.Entry(self.crystal_frame, textvariable=row["thick"], width=9), ttk.Entry(self.crystal_frame, textvariable=row["comp"], width=18)]
        for col, w in enumerate(ws):
            w.grid(row=k, column=col, sticky="w", padx=(0, 6), pady=1)

        def sync(event=None, ec=ec, row=row):
            ec.state(["!disabled"] if UI_CRYSTAL[row["crystal"].get()] == "hcp0001" else ["disabled"])
        cb.bind("<<ComboboxSelected>>", sync)
        sync()
        row["widgets"] = ws
        self.crystal_rows.append(row)

    def _remove_crystal(self):
        if len(self.crystal_rows) > 1:
            for w in self.crystal_rows.pop()["widgets"]:
                w.destroy()

    def _clear_crystal(self):
        while self.crystal_rows:
            for w in self.crystal_rows.pop()["widgets"]:
                w.destroy()

    def _example_crystal_tial_rows(self):
        for crystal, a, c, comp in (("HCP (0001)", "0.2951", "0.4683", "100,0"), ("FCC (111)", "0.405", "", "0,100")) * 2:
            self._add_crystal(crystal, a, c, "4", comp)

    def _example_crystal_tial(self):
        self._set_elements(["Ti", "Al"])
        self._clear_crystal()
        for crystal, a, c, comp in (("HCP (0001)", "0.2951", "0.4683", "100,0"), ("FCC (111)", "0.405", "", "0,100")) * 2:
            self._add_crystal(crystal, a, c, "4", comp)
        self._on_n_elements()

    def _example_crystal_nicu(self):
        self._set_elements(["Ni", "Cu"])
        self._clear_crystal()
        for a, comp in (("0.3524", "100,0"), ("0.3615", "0,100")) * 2:
            self._add_crystal("FCC (111)", a, "", "4", comp)
        self._on_n_elements()

    def _collect_crystal(self):
        v = self.vl
        shape = UI_SHAPES_L[v["shape"].get()]
        n = self._n_el()
        elements, seen = [], set()
        for i in range(n):
            r = self.el_rows[i]
            sym = r["sym"].get().strip()
            if not sym:
                raise FormError("Element %d needs a symbol." % (i + 1))
            if sym in seen:
                raise FormError("Element symbols must differ (%s appears twice)." % sym)
            seen.add(sym)
            elements.append({"symbol": sym, "mass": self._to_float(r["mass"].get(), "Mass of %s" % sym, positive=True)})
        layers = []
        for k, row in enumerate(self.crystal_rows, 1):
            crystal = UI_CRYSTAL[row["crystal"].get()]
            ctxt = row["c"].get().strip()
            layers.append(dict(
                crystal=crystal, a_nm=self._to_float(row["a"].get(), "Layer %d a" % k, positive=True),
                c_nm=self._to_float(ctxt, "Layer %d c" % k, positive=True) if (crystal == "hcp0001" and ctxt) else None,
                thickness_nm=self._to_float(row["thick"].get(), "Layer %d thickness" % k, positive=True),
                start=self._recipe(row["comp"].get(), n, "Layer %d composition" % k), end=None))
        if not layers:
            raise FormError("Add at least one layer.")
        mode = UI_INPLANE[v["cr_inplane"].get()]
        seed = v["seed"].get().strip()
        return dict(shape=shape,
                    radius_nm=self._to_float(v["radius"].get(), "Radius", positive=True) if shape == "cylindrical" else 2.0,
                    xdim_nm=self._to_float(v["xdim"].get(), "X-dimension", positive=True) if shape == "cuboid" else 2.0,
                    ydim_nm=self._to_float(v["ydim"].get(), "Y-dimension", positive=True) if shape == "cuboid" else 2.0,
                    elements=elements, layers=layers, inplane=mode,
                    inplane_nm=self._to_float(v["cr_inplane_nm"].get(), "In-plane nearest-neighbour distance", positive=True) if mode == "value" else None,
                    basis=UI_BASIS[v["basis"].get()], sample=UI_SAMPLES[v["sample"].get()],
                    vacuum_xy=self._to_float(v["vac_xy"].get(), "Vacuum x,y", minimum=0),
                    vacuum_z=self._to_float(v["vac_z"].get(), "Vacuum z", minimum=0),
                    seed=self._to_float(seed, "Random seed", integer=True) if seed else None)

    def on_build_crystal(self):
        try:
            kw = self._collect_crystal()
        except FormError as e:
            messagebox.showerror("Please check the form", str(e))
            return
        self._busy(True)
        try:
            t0 = time.time()
            res = crystal_layers_build(**kw)
            dt = time.time() - t0
        except Exception as e:
            self._busy(False)
            messagebox.showerror("Build failed", "%s: %s" % (type(e).__name__, e))
            return
        self._busy(False)
        self.res, self.params, self.curves = res, ("crystals", kw), None
        self.btn_save.state(["!disabled"])
        strains = ", ".join("%+.1f%%" % ly["strain_pct"] for ly in res["layers"])
        self.status.set("Built %d atoms in %.1f s: %d crystal layers (in-plane strain %s)." % (len(res["atoms"]), dt, len(res["layers"]), strains))
        self.show_view(V_DATA)

    def _example_comb_binary(self):
        self._set_elements(["Ni", "Co"])
        for key, val in (("c_len", "10"), ("c_cb", "100,0"), ("c_sb", "70,30"), ("c_ct", "30,70"), ("c_st", "0,100"),
                         ("c_fz", "S-FGM"), ("c_pz", "2"), ("c_bz", "10"), ("c_fr", "P-FGM"), ("c_pr", "1"), ("c_br", "5")):
            self.vl[key].set(val)
        self._on_n_elements()

    def _example_comb_solid(self):
        """Solid regions: pure Ni core and bottom, pure Co shell and top (Custom steps in both directions)."""
        self._set_elements(["Ni", "Co"])
        for key, val in (("c_len", "10"), ("c_cb", "100,0"), ("c_sb", "0,100"), ("c_ct", "0,100"), ("c_st", "0,100"),
                         ("c_fz", "Custom"), ("c_pz", "x>=0.5"), ("c_bz", "2"), ("c_fr", "Custom"), ("c_pr", "x>=0.5"), ("c_br", "2")):
            self.vl[key].set(val)
        self._on_n_elements()

    def _example_comb_ternary(self):
        self._set_elements(["Co", "Cr", "Ni"], latmode="One value", one_a="0.354")
        for key, val in (("c_len", "10"), ("c_cb", "80,10,10"), ("c_sb", "10,80,10"), ("c_ct", "10,10,80"), ("c_st", "33.3,33.3,33.4"),
                         ("c_fz", "P-FGM"), ("c_pz", "1"), ("c_bz", "8"), ("c_fr", "S-FGM"), ("c_pr", "2"), ("c_br", "4")):
            self.vl[key].set(val)
        self._on_n_elements()

    def _collect_common(self):
        v = self.vl
        shape = UI_SHAPES_L[v["shape"].get()]
        n = self._n_el()
        one = UI_LATMODE[v["latmode"].get()] == "one"
        one_a = self._to_float(v["one_a"].get(), "Lattice parameter (one value)", positive=True) if one else None
        elements, seen = [], set()
        for i in range(n):
            r = self.el_rows[i]
            sym = r["sym"].get().strip()
            if not sym:
                raise FormError("Element %d needs a symbol." % (i + 1))
            if sym in seen:
                raise FormError("Element symbols must differ (%s appears twice)." % sym)
            seen.add(sym)
            mass = self._to_float(r["mass"].get(), "Mass of %s" % sym, positive=True)
            need_own = (not one) or UI_ZSP[v["zspace"].get()] == "vegard"
            a = self._to_float(r["a"].get(), "Lattice parameter of %s" % sym, positive=True) if need_own else one_a
            elements.append({"symbol": sym, "mass": mass, "a_nm": a})
        seed = v["seed"].get().strip()
        return n, dict(shape=shape,
                       radius_nm=self._to_float(v["radius"].get(), "Radius", positive=True) if shape == "cylindrical" else 2.0,
                       xdim_nm=self._to_float(v["xdim"].get(), "X-dimension", positive=True) if shape == "cuboid" else 2.0,
                       ydim_nm=self._to_float(v["ydim"].get(), "Y-dimension", positive=True) if shape == "cuboid" else 2.0,
                       lattice=UI_LATTS[v["lattice"].get()], elements=elements,
                       basis=UI_BASIS[v["basis"].get()], sample=UI_SAMPLES[v["sample"].get()],
                       vacuum_xy=self._to_float(v["vac_xy"].get(), "Vacuum x,y", minimum=0),
                       vacuum_z=self._to_float(v["vac_z"].get(), "Vacuum z", minimum=0),
                       seed=self._to_float(seed, "Random seed", integer=True) if seed else None, lattice_param_nm=one_a,
                       z_spacing=UI_ZSP[v["zspace"].get()])

    def _collect_combined(self):
        v = self.vl
        n, kw = self._collect_common()
        corners = [self._recipe(v[k].get(), n, "%s composition" % nm) for k, nm in
                   (("c_cb", "Core, bottom"), ("c_sb", "Surface, bottom"), ("c_ct", "Core, top"), ("c_st", "Surface, top"))]
        gz = self._grading_inputs(v["c_fz"].get(), v["c_pz"].get(), "Along the wire")
        gr_ = self._grading_inputs(v["c_fr"].get(), v["c_pr"].get(), "Core to surface")
        kw.update(length_nm=self._to_float(v["c_len"].get(), "Length", positive=True), corners=corners,
                  grading_z=gz[0], p_z=gz[1], custom_z=gz[2],
                  bins_z=self._to_float(v["c_bz"].get(), "Slices", minimum=1, integer=True),
                  grading_r=gr_[0], p_r=gr_[1], custom_r=gr_[2],
                  bins_r=self._to_float(v["c_br"].get(), "Rings", minimum=1, integer=True))
        return kw

    def on_build_combined(self):
        try:
            kw = self._collect_combined()
        except FormError as e:
            messagebox.showerror("Please check the form", str(e))
            return
        self._busy(True)
        try:
            t0 = time.time()
            res = combined_build(**kw)
            dt = time.time() - t0
        except Exception as e:
            self._busy(False)
            messagebox.showerror("Build failed", "%s: %s" % (type(e).__name__, e))
            return
        self._busy(False)
        self.res, self.params, self.curves = res, ("combined", kw), None
        self.btn_save.state(["!disabled"])
        extra = ("  " + res["messages"][0]) if res["messages"] else ""
        self.status.set("Built %d atoms in %.1f s: %d slices x %d rings = %d cells.%s" % (
            len(res["atoms"]), dt, res["combined"]["bins_z"], res["combined"]["bins_r"], len(res["bin_info"]), extra))
        self.show_view(V_DATA)

    def _draw_rprofile(self):
        c = self.canvas
        c.delete("all")
        W, H = c.winfo_width(), c.winfo_height()
        res = self.res
        if W < 80 or H < 80:
            return
        if res is None or res.get("mode") != "combined":
            c.create_text(W / 2, H / 2, text="Nothing built yet. Set up the combined grading and press Build.", fill="#555")
            return
        atoms = res["atoms"]
        ty = atoms[:, 1].astype(int)
        d = atoms[:, 2:4] - res["axis"][:2]
        rad = np.sqrt((d ** 2).sum(axis=1))
        r_max = float(rad.max())
        nb = 12
        width = r_max / nb
        ib = np.minimum(np.floor(rad / width).astype(int), nb - 1)
        n_el = len(res["elements"])
        counts = np.zeros((nb, n_el))
        np.add.at(counts, (ib, ty - 1), 1)
        pct = 100 * counts / np.maximum(counts.sum(axis=1, keepdims=True), 1)
        rc = (np.arange(nb) + 0.5) * width / 10
        px, py = self._frame(c, W, H, 0, r_max / 10, "Distance from the wire axis (nm)", "Atom % of each element in the ring", xfmt="%.2f")
        for e in range(n_el):
            pts = [v_ for i in range(nb) for v_ in (px(rc[i]), py(pct[i, e]))]
            c.create_line(pts, fill=PALETTE[e % len(PALETTE)], width=2)
        self._legend(c, [str(x.get("symbol", i + 1)) for i, x in enumerate(res["elements"])])

    # ---------------------------------------------- layered result views
    def _frame(self, c, W, H, xmin, xmax, xlabel, ylabel, ymax=100.0, xfmt="%.0f"):
        L, R, T, B = 70, 25, 55, 60
        pw, ph = W - L - R, H - T - B
        px = lambda x: L + (x - xmin) / (xmax - xmin) * pw
        py = lambda y: T + ph - y / ymax * ph
        c.create_rectangle(L, T, L + pw, T + ph, outline="#444")
        for i in range(6):
            xv, yv = xmin + (xmax - xmin) * i / 5, ymax * i / 5
            c.create_line(px(xv), T + ph, px(xv), T + ph + 5, fill="#444")
            c.create_text(px(xv), T + ph + 16, text=xfmt % xv)
            c.create_line(L - 5, py(yv), L, py(yv), fill="#444")
            c.create_text(L - 28, py(yv), text="%d" % yv)
            if 0 < i < 5:
                c.create_line(px(xv), T, px(xv), T + ph, fill="#e8e8e8")
                c.create_line(L, py(yv), L + pw, py(yv), fill="#e8e8e8")
        c.create_text(L + pw / 2, H - 14, text=xlabel)
        c.create_text(18, T + ph / 2, text=ylabel, angle=90)
        return px, py

    def _legend(self, c, names, x0=70, y=22):
        for i, nme in enumerate(names):
            col = PALETTE[i % len(PALETTE)]
            c.create_rectangle(x0, y - 6, x0 + 14, y + 6, fill=col, outline=col)
            c.create_text(x0 + 20, y, text=nme, anchor="w")
            x0 += 24 + 9 * len(nme) + 18

    def _draw_canvas(self, name):
        {V_CURVES: self._draw_curves, V_PROFILE: self._draw_profile, V_SIDE: self._draw_side, V_RPROFILE: self._draw_rprofile}[name]()

    def _draw_profile(self):
        c = self.canvas
        c.delete("all")
        W, H = c.winfo_width(), c.winfo_height()
        res = self.res
        if W < 80 or H < 80:
            return
        if res is None or res.get("mode") not in ("layered", "combined", "crystals"):
            c.create_text(W / 2, H / 2, text="Nothing built yet. Set up the wire and press Build.", fill="#555")
            return
        atoms = res["atoms"]
        half = res["lat_p"] / 2
        ty = atoms[:, 1].astype(int)
        plane = res["plane"] if res.get("plane") is not None else np.rint(atoms[:, 4] / half).astype(int)
        K = int(plane.max()) + 1
        per = max(3, int(round(K / 45)))  # planes per slab (a few planes smooth the random noise)
        slab = plane // per
        ns = int(slab.max()) + 1
        n_el = len(res["elements"])
        counts = np.zeros((ns, n_el))
        np.add.at(counts, (slab, ty - 1), 1)
        pct = 100 * counts / np.maximum(counts.sum(axis=1, keepdims=True), 1)
        zc = np.array([atoms[slab == i, 4].mean() for i in range(ns)]) / 10
        zmax = (res["box"]["z"][1] - res["box"]["z"][0]) / 10
        px, py = self._frame(c, W, H, 0, zmax, "Position along the wire (nm)", "Atom % of each element in the slab", xfmt="%.1f")
        for ly in res.get("layers", []):
            x1 = ly["z1"] / 10
            c.create_line(px(x1), 55, px(x1), H - 60, fill="#999", dash=(4, 4))
            c.create_text(px((ly["z0"] + ly["z1"]) / 20), 66, text="L%d" % ly["layer"], fill="#777")
        for e in range(n_el):
            pts = [v for i in range(ns) for v in (px(zc[i]), py(pct[i, e]))]
            c.create_line(pts, fill=PALETTE[e % len(PALETTE)], width=2)
        self._legend(c, [str(x.get("symbol", i + 1)) for i, x in enumerate(res["elements"])])

    def _draw_side(self):
        c = self.canvas
        c.delete("all")
        W, H = c.winfo_width(), c.winfo_height()
        res = self.res
        if W < 80 or H < 80:
            return
        if res is None or res.get("mode") not in ("layered", "combined", "crystals"):
            c.create_text(W / 2, H / 2, text="Nothing built yet. Set up the wire and press Build.", fill="#555")
            return
        a = res["atoms"]
        yc = float(res["axis"][1])
        lat = res["lat_p"]
        thick = lat * 0.55
        sel = np.abs(a[:, 3] - yc) <= thick
        if sel.sum() > 6000:
            thick = lat * 0.30
            sel = np.abs(a[:, 3] - yc) <= thick
        xs, zs, ty = a[sel, 2], a[sel, 4], a[sel, 1].astype(int)
        wd, ht = float(a[:, 2].max()) + lat, float(a[:, 4].max()) + lat
        top, bottom, side = 45, 30, 40
        scale = min((H - top - bottom) / ht, (W - 2 * side) / wd)
        r = max(1.2, scale * lat * 0.30)
        x0 = (W - wd * scale) / 2
        for xv, zv, tv in zip(xs, zs, ty):
            X, Y = x0 + (xv + lat / 2) * scale, H - bottom - (zv + lat / 2) * scale
            col = PALETTE[(tv - 1) % len(PALETTE)]
            c.create_oval(X - r, Y - r, X + r, Y + r, fill=col, outline=col)
        for ly in res.get("layers", []):
            Y = H - bottom - ly["z1"] * scale
            c.create_line(x0 - 14, Y, x0 + wd * scale + 14, Y, fill="#999", dash=(3, 3))
            c.create_text(x0 - 22, H - bottom - (ly["z0"] + ly["z1"]) / 2 * scale, text="L%d" % ly["layer"], fill="#555", anchor="e")
        self._legend(c, [str(x.get("symbol", i + 1)) for i, x in enumerate(res["elements"])], x0=side)
        c.create_text(W - 10, H - 10, anchor="e", fill="#777",
                      text="thin slice through the wire axis (|y - axis| <= %.1f A), wire shown standing up" % thick)

    def _busy(self, on):
        self.config(cursor="watch" if on else "")
        self.update_idletasks()

    def on_build(self):
        if self.tab == "layered":
            self.on_build_layered()
            return
        if self.tab == "combined":
            self.on_build_combined()
            return
        if self.tab == "crystal":
            self.on_build_crystal()
            return
        try:
            params = self._collect()
        except FormError as e:
            messagebox.showerror("Please check the form", str(e))
            return
        self._busy(True)
        try:
            t0 = time.time()
            res = build_structure(**params)
            self.res, self.params = res, params
            self.curves = self._make_curves()
            dt = time.time() - t0
        except Exception as e:
            self._busy(False)
            messagebox.showerror("Build failed", "%s: %s" % (type(e).__name__, e))
            return
        self._busy(False)
        self.btn_save.state(["!disabled"])
        extra = ("  " + res["messages"][0]) if res["messages"] else ""
        self.status.set("Built %d atoms in %.1f s (mode: %s).%s" % (len(res["atoms"]), dt, res["mode"], extra))
        self.show_view(V_DATA)

    def on_manual(self, event=None):
        """Show the manual in a window inside the program."""
        old = getattr(self, "manual_win", None)
        if old is not None and old.winfo_exists():
            old.lift()
            old.focus_force()
            return
        win = tk.Toplevel(self)
        win.title("Manual - FGM Nanostructure Builder")
        win.geometry("1040x740")
        self.manual_win = win
        win.columnconfigure(1, weight=1)
        win.rowconfigure(0, weight=1)
        win.bind("<Escape>", lambda e: win.destroy())
        mono = ("Consolas", 10) if sys.platform.startswith("win") else ("Courier", 10)
        toc = tk.Listbox(win, width=36, exportselection=False, font=("Segoe UI", 10), activestyle="none")
        toc.grid(row=0, column=0, sticky="ns", padx=(8, 0), pady=8)
        txt = tk.Text(win, wrap="word", padx=16, pady=10, font=("Segoe UI", 10), cursor="arrow", spacing1=2, spacing3=3)
        sb = ttk.Scrollbar(win, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        txt.grid(row=0, column=1, sticky="nsew", pady=8)
        sb.grid(row=0, column=2, sticky="ns", pady=8, padx=(0, 4))
        txt.tag_configure("h1", font=("Segoe UI", 19, "bold"), spacing1=4, spacing3=2)
        txt.tag_configure("sub", font=("Segoe UI", 11, "italic"), foreground="#444", spacing3=10)
        txt.tag_configure("h2", font=("Segoe UI", 14, "bold"), spacing1=16, spacing3=6, foreground="#0b3d91")
        txt.tag_configure("h3", font=("Segoe UI", 11, "bold"), spacing1=10, spacing3=3)
        txt.tag_configure("bullet", lmargin1=22, lmargin2=38, spacing3=2)
        txt.tag_configure("mono", font=mono, background="#f1f1f1", lmargin1=18, lmargin2=18, spacing1=0, spacing3=0)
        marks = []
        lines = MANUAL_TEXT.split("\n")
        i = 0
        while i < len(lines):
            ln = lines[i]
            if ln.startswith("```"):
                j = i + 1
                code = []
                while j < len(lines) and not lines[j].startswith("```"):
                    code.append(lines[j])
                    j += 1
                txt.insert("end", "\n".join(code) + "\n", "mono")
                i = j + 1
                continue
            if ln.startswith("# "):
                txt.insert("end", ln[2:] + "\n", "h1")
                if i + 2 < len(lines) and lines[i + 2].strip():
                    txt.insert("end", lines[i + 2].replace("`", "") + "\n", "sub")
                    i += 2
            elif ln.startswith("## "):
                name = "sec%d" % len(marks)
                txt.mark_set(name, "end-1c")
                txt.mark_gravity(name, "left")
                marks.append(name)
                toc.insert("end", ln[3:])
                txt.insert("end", ln[3:] + "\n", "h2")
            elif ln.startswith("### "):
                txt.insert("end", ln[4:] + "\n", "h3")
            elif ln.startswith("- "):
                txt.insert("end", "\u2022  " + ln[2:].replace("`", "") + "\n", "bullet")
            elif ln.startswith("|"):
                rows = []
                while i < len(lines) and lines[i].startswith("|"):
                    cells = [c.strip().replace("`", "") for c in lines[i].strip("|").split("|")]
                    if not all(set(c) <= set("-: ") for c in cells):
                        rows.append(cells)
                    i += 1
                widths = [max(len(r[k]) for r in rows if k < len(r)) for k in range(max(len(r) for r in rows))]
                for r in rows:
                    txt.insert("end", "  ".join(c.ljust(widths[k]) for k, c in enumerate(r)).rstrip() + "\n", "mono")
                continue
            elif ln.strip():
                para = [ln.replace("`", "")]
                while i + 1 < len(lines) and lines[i + 1].strip() and not lines[i + 1].startswith(("#", "- ", "|", "```")):
                    i += 1
                    para.append(lines[i].replace("`", ""))
                txt.insert("end", " ".join(para) + "\n")
            i += 1
        txt.configure(state="disabled")
        toc.bind("<<ListboxSelect>>", lambda e: txt.yview(marks[toc.curselection()[0]]) if toc.curselection() else None)
        self.manual_text_widget, self.manual_toc = txt, toc

    def on_save(self):
        if self.res is None:
            return
        path = filedialog.asksaveasfilename(
            title="Save LAMMPS data file", defaultextension=".nw", initialfile="fgm_nw.nw",
            filetypes=[("LAMMPS data file", "*.nw *.data *.txt"), ("All files", "*.*")])
        if not path:
            return
        try:
            write_lammps_data(path, self.res)
        except Exception as e:
            messagebox.showerror("Save failed", str(e))
            return
        self.status.set("Saved %d atoms to %s" % (len(self.res["atoms"]), path))

    def on_compare(self):
        if self.tab != "classic":
            messagebox.showinfo("Compare with reference file", "The comparison is available on the first tab (FGM Nanostructure Builder) only: there is no standard reference for multilayer, combined or multi-element wires.")
            return
        try:
            params = self._collect()
        except FormError as e:
            messagebox.showerror("Please check the form", str(e))
            return
        path = filedialog.askopenfilename(
            title="Pick the reference data file",
            filetypes=[("LAMMPS data file", "*.nw *.data *.txt *.dat"), ("All files", "*.*")])
        if not path:
            return
        note = []
        if params["mode"] != "nanohub":
            note.append("Note: the comparison always uses the Standard mode.")
        if params["vacuum_xy"] or params["vacuum_z"]:
            note.append("Note: vacuum padding is ignored for the comparison.")
        params.update(mode="nanohub", vacuum_xy=0.0, vacuum_z=0.0)
        self._busy(True)
        try:
            res = build_structure(**params)
            lines, ok = run_checks(res, path)
        except Exception as e:
            self._busy(False)
            messagebox.showerror("Comparison failed",
                                 "%s: %s\n\nIs this a LAMMPS data file?" % (type(e).__name__, e))
            return
        self._busy(False)
        head = ["Compared with: %s" % path, "Using the values currently in the form."] + note + [""]
        self.compare_lines = head + lines
        self.status.set("Comparison finished: %s" % ("ALL CHECKS PASSED" if ok else "something differs"))
        self.show_view(V_COMPARE)


def hide_own_console():
    """Double-clicked on Windows: hide the black console window that opened only for us."""
    if not sys.platform.startswith("win"):
        return
    try:
        import ctypes
        k = ctypes.windll.kernel32
        procs = (ctypes.c_uint * 8)()
        if k.GetConsoleProcessList(procs, 8) <= 2:  # only py.exe/python.exe are attached
            hwnd = k.GetConsoleWindow()
            if hwnd:
                ctypes.windll.user32.ShowWindow(hwnd, 0)
    except Exception:
        pass


def main():
    if tk is None:
        print("The window needs Tkinter, which is not installed here.\nUse the command line instead (see: python FGM_NW_Builder.py --help) or the web version.")
        return
    if sys.platform.startswith("win"):
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)  # sharper text on high-DPI screens
        except Exception:
            pass
    try:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
    except Exception:
        pass
    try:
        if not ensure_numpy(True):
            return
        hide_own_console()
        App().mainloop()
    except Exception:
        import traceback
        text = traceback.format_exc()
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("FGM Nanostructure Builder - unexpected error", text)
            root.destroy()
        except Exception:
            print(text)


def entry():
    args = sys.argv[1:]
    if args and args[0] in ("build", "validate", "selftest", "layered", "combined", "crystals"):
        if not ensure_numpy(False):
            sys.exit(1)
        {"build": cli_build, "validate": cli_validate, "selftest": cli_selftest, "layered": cli_layered, "combined": cli_combined, "crystals": cli_crystals}[args[0]](args[1:])
    elif args and args[0] in ("-h", "--help", "help"):
        print(__doc__)
    else:
        main()


if __name__ == "__main__":
    entry()
