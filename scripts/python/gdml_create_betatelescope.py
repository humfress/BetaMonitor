"""Generate the BetaMonitor telescope GDML geometry with pyg4ometry.

The dimensions and placements reproduce ``original_geometry_telescope_export.gdml``.
All dimensions are in millimeters.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pyg4ometry.geant4 as g4
import pyg4ometry.gdml as gdml


@dataclass(frozen=True)
class GeometryParams:
    """Dimensions and placements for the telescope detector assembly."""

    world_xyz: tuple[float, float, float] = (400.0, 400.0, 450.0)

    # T-pipe and flanges.
    pipe_inner_radius: float = 17.399
    pipe_outer_radius: float = 19.05
    pipe_inner_major_length: float = 196.85
    pipe_inner_minor_length: float = 98.425
    pipe_outer_major_length: float = 171.45
    pipe_outer_minor_length: float = 85.725
    pipe_inner_offset_x: float = -49.2125
    pipe_outer_offset_x: float = -42.8625
    pipe_z: float = 111.125
    flange_inner_radius: float = 17.399
    flange_outer_radius: float = 34.544
    flange_thickness: float = 12.7
    front_flange_z: float = 6.35
    rear_flange_z: float = 19.05
    side_flange_offset_x: float = 92.075

    # Thin entrance window.
    window_radius: float = 10.0
    window_thickness: float = 0.05
    window_z: float = 0.025

    # Telescope scintillators. Their aluminum skins are 0.0001 mm larger.
    scintillator_radius: float = 25.0
    aluminum_skin_delta: float = 0.0001
    scintillator_a_thickness: float = 2.0
    scintillator_a_z: float = -9.0001
    scintillator_b_thickness: float = 50.0
    scintillator_b_z: float = -38.0002


PARAMS = GeometryParams()
DEGREES = "deg"


def build_materials(registry: g4.Registry) -> dict[str, object]:
    """Create the standard and experiment-specific materials."""
    materials = {
        name: g4.MaterialPredefined(name, registry)
        for name in (
            "G4_Al",
            "G4_PLASTIC_SC_VINYLTOLUENE",
            "G4_AIR",
        )
    }
    elements = {
        name: g4.ElementSimple(name, symbol, atomic_number, atomic_mass, registry)
        for name, symbol, atomic_number, atomic_mass in (
            ("Hydrogen", "H", 1, 1.008),
            ("Chromium", "Cr", 24, 51.996),
            ("Manganese", "Mn", 25, 54.938),
            ("Iron", "Fe", 26, 55.845),
            ("Nickel", "Ni", 28, 58.693),
        )
    }

    def mixture(
        name: str,
        density: float,
        components: tuple[tuple[float, str], ...],
        state: str,
    ) -> object:
        material = g4.MaterialCompound(name, density, len(components), registry, state=state)
        for mass_fraction, element in components:
            material.add_element_massfraction(elements[element], mass_fraction)
        return material

    materials["interGalactic"] = mixture("interGalactic", 1e-25, ((1.0, "Hydrogen"),), "gas")
    materials["Stainless_Steel"] = mixture(
        "Stainless_Steel",
        8.03,
        ((0.72, "Iron"), (0.18, "Chromium"), (0.08, "Nickel"), (0.02, "Manganese")),
        "solid",
    )
    return materials


def tube(name: str, rmin: float, rmax: float, length: float, registry: g4.Registry) -> object:
    """Create a full cylindrical tube."""
    return g4.solid.Tubs(name, rmin, rmax, length, 0.0, 360.0, registry, aunit=DEGREES)


def transform(
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
    rotation: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> list[list[float | str]]:
    """Return a pyg4ometry boolean-solid transform in degrees and millimeters."""
    return [list(rotation) + [DEGREES], list(position)]


def build_solids(params: GeometryParams, registry: g4.Registry) -> dict[str, object]:
    """Create the T-pipe, flanges, window, and two telescope scintillators."""
    p = params
    solids: dict[str, object] = {}

    def add(solid: object) -> object:
        solids[solid.name] = solid
        return solid

    # The T-pipe is a vertical tube with a perpendicular branch.
    pipe_1i = add(tube("Pipe1i", 0.0, p.pipe_inner_radius, p.pipe_inner_major_length, registry))
    pipe_2i = add(tube("Pipe2i", 0.0, p.pipe_inner_radius, p.pipe_inner_minor_length, registry))
    add(g4.solid.Union("InnerTPipe", pipe_1i, pipe_2i, transform((p.pipe_inner_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))
    pipe_1o = add(tube("Pipe1o", 0.0, p.pipe_outer_radius, p.pipe_outer_major_length, registry))
    pipe_2o = add(tube("Pipe2o", 0.0, p.pipe_outer_radius, p.pipe_outer_minor_length, registry))
    outer_pipe = add(g4.solid.Union("OuterTPipe", pipe_1o, pipe_2o, transform((p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))
    inner_main = add(tube("Pipe1i_tol", 0.0, p.pipe_inner_radius + 1e-6, 184.15, registry))
    steel_pipe = add(g4.solid.Subtraction("TDecayVolume", outer_pipe, inner_main, transform(), registry))
    inner_branch = add(tube("Pipe2i_tol", 0.0, p.pipe_inner_radius + 1e-6, p.pipe_inner_minor_length, registry))
    add(g4.solid.Subtraction("SteelTPipeSolid", steel_pipe, inner_branch, transform((p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))

    add(tube("FlangeScint", 10.0, p.flange_outer_radius, p.flange_thickness, registry))
    add(tube("Flange", p.flange_inner_radius, p.flange_outer_radius, p.flange_thickness, registry))
    add(tube("VacuumWindowDisk", 0.0, p.window_radius, p.window_thickness, registry))

    for name, thickness in (
        ("AScintCore", p.scintillator_a_thickness),
        ("BScintCore", p.scintillator_b_thickness),
    ):
        core = add(tube(name, 0.0, p.scintillator_radius, thickness, registry))
        aluminum = add(tube(f"{name}Aluminum", 0.0, p.scintillator_radius + p.aluminum_skin_delta, thickness + 2 * p.aluminum_skin_delta, registry))
        add(g4.solid.Subtraction(f"{name}InnerAl", aluminum, core, transform(), registry))

    add(g4.solid.Box("World", *p.world_xyz, registry))
    return solids


def place(
    name: str,
    material: object,
    solid: object,
    world: g4.LogicalVolume,
    registry: g4.Registry,
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
    rotation: tuple[float, float, float] = (0.0, 0.0, 0.0),
    copy_number: int = 0,
) -> None:
    """Create a logical volume and place it directly in the world."""
    volume = g4.LogicalVolume(solid, material, name, registry)
    g4.PhysicalVolume(list(rotation) + [DEGREES], position, volume, name, world, registry, copyNumber=copy_number)


def build_geometry(output_path: str | Path, params: GeometryParams = PARAMS) -> Path:
    """Build and write the complete telescope GDML geometry to ``output_path``."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    registry = g4.Registry()
    materials = build_materials(registry)
    solids = build_solids(params, registry)
    world = g4.LogicalVolume(solids["World"], materials["G4_AIR"], "World", registry)
    p = params

    placements = (
        ("Vacuum", "interGalactic", "InnerTPipe", (0.0, 0.0, p.pipe_z), (0.0, 90.0, 0.0), 1),
        ("SteelTPipe", "Stainless_Steel", "SteelTPipeSolid", (0.0, 0.0, p.pipe_z), (0.0, 90.0, 0.0), 0),
        ("tPipeFlangeScint", "Stainless_Steel", "FlangeScint", (0.0, 0.0, p.front_flange_z), (0.0, 0.0, 0.0), 0),
        ("tPipeFlange1", "Stainless_Steel", "Flange", (0.0, 0.0, p.rear_flange_z), (0.0, 0.0, 0.0), 0),
        ("tPipeFlange2", "Stainless_Steel", "Flange", (p.side_flange_offset_x, 0.0, p.pipe_z), (0.0, 90.0, 0.0), 0),
        ("tPipeFlange3", "Stainless_Steel", "Flange", (-p.side_flange_offset_x, 0.0, p.pipe_z), (0.0, 90.0, 0.0), 0),
        ("VacuumWindow", "G4_Al", "VacuumWindowDisk", (0.0, 0.0, p.window_z), (0.0, 0.0, 0.0), 2),
        ("AScintillator", "G4_PLASTIC_SC_VINYLTOLUENE", "AScintCore", (0.0, 0.0, p.scintillator_a_z), (0.0, 0.0, 0.0), 3),
        ("Ali_sq10", "G4_Al", "AScintCoreInnerAl", (0.0, 0.0, p.scintillator_a_z), (0.0, 0.0, 0.0), 0),
        ("BScintillatorLV", "G4_PLASTIC_SC_VINYLTOLUENE", "BScintCore", (0.0, 0.0, p.scintillator_b_z), (0.0, 0.0, 0.0), 4),
        ("Ali_sq20", "G4_Al", "BScintCoreInnerAl", (0.0, 0.0, p.scintillator_b_z), (0.0, 0.0, 0.0), 0),
    )
    for name, material, solid, position, rotation, copy_number in placements:
        place(name, materials[material], solids[solid], world, registry, position, rotation, copy_number)

    registry.setWorld(world.name)
    writer = gdml.Writer()
    writer.addDetector(registry)
    writer.write(path)
    return path


def main(argv: Sequence[str] | None = None) -> int:
    """Parse CLI arguments and write the requested telescope GDML export."""
    parser = argparse.ArgumentParser(
        description="Generate telescope geometry_export.gdml for BetaMonitor."
    )
    parser.add_argument("output", nargs="?", default="./dat/geometry_telescope_export.gdml", help="Output GDML path.")
    args = parser.parse_args(argv)
    print(build_geometry(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
