# FGM Nanostructure Builder

A desktop program that builds atom-position files (LAMMPS data format) for functionally graded (FGM) and multilayer / multi-element nanowires, cuboids and nanoparticles. It writes structure files only: it does not write LAMMPS input scripts or run simulations.

Author: Md. Shabbir Rahman Siam. Released under the MIT License (see `LICENSE.txt`).

## Download and run

**Windows, no Python needed:** go to the [Releases](../../releases) page and download `FGM_NW_Builder.exe`. Windows may show a "protected your PC" warning because the file is not code-signed: choose *More info GåÆ Run anyway*.

**Any system with Python 3 (Windows, macOS, Linux):**

    pip install -r requirements.txt
    python FGM_NW_Builder.py

Tkinter is needed for the window. It comes with Python on Windows and macOS; on Ubuntu/Debian install it with `sudo apt install python3-tk`.

The full manual is inside the program: **Help > Manual** (key F1).

## Command line (no window needed)

    python FGM_NW_Builder.py build --direction axial --bins 10 --out wire.nw
    python FGM_NW_Builder.py layered --example > spec.json
    python FGM_NW_Builder.py layered --spec spec.json --out wire.nw
    python FGM_NW_Builder.py combined --example > spec.json
    python FGM_NW_Builder.py crystals --example > spec.json
    python FGM_NW_Builder.py selftest
    python FGM_NW_Builder.py build --help

## Checking the program

`python FGM_NW_Builder.py selftest` runs the built-in checks (exact atom counts, crystal geometry, spacings, periodic seams). They also run automatically on every push through GitHub Actions.

## For the maintainer: publishing a new version

1. Update `VERSION` near the top of `FGM_NW_Builder.py` and add a line to `CHANGELOG.md`.
2. Commit, then tag and push the tag:

       git tag v1.1
       git push origin main --tags

3. GitHub Actions builds `FGM_NW_Builder.exe` and attaches it to a new Release automatically.

## How to cite

See `CITATION.cff` (GitHub shows a "Cite this repository" button).
