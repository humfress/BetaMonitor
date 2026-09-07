"""Generate the BetaMonitor GDML geometry with pyg4ometry.

All dimensions are in millimeters.  Scintillator A and B dimensions are
independent so studies can adjust either detector without changing its mate.
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
    """Dimensions and positions for the physical BetaMonitor assembly."""

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
    world_pipe_z: float = 98.425
    flange_inner_radius: float = 17.399
    flange_outer_radius: float = 34.544
    flange_thickness: float = 12.7
    flange_scint_z: float = -6.35
    flange_1_z: float = 6.35
    flange_side_offset_x: float = 92.075
    flange_side_z: float = 98.425

    # Window and wedge.
    vacuum_window_radius: float = 17.399
    vacuum_window_thickness: float = 0.0762
    vacuum_window_z: float = -0.0381
    wedge_zplanes: tuple[tuple[float, float, float], ...] = (
        (17.399, 9.906, 0.0),
        (17.399, 9.906, 1.0),
        (17.399, 17.398999, 7.0),
    )

    # Both scintillators use the same reference solid dimensions.
    scint_a_pos_z: float = -9.9836
    scint_a_xy: float = 25.0
    scint_a_z: float = 3.0
    scint_a_inner_al_delta: float = 0.0002
    scint_a_mylar_delta: float = 0.1522
    scint_a_outer_al_delta: float = 0.1524

    scint_b_pos_z: float = -13.136
    scint_b_xy: float = 25.0
    scint_b_z: float = 3.0
    scint_b_inner_al_delta: float = 0.0002
    scint_b_mylar_delta: float = 0.1522
    scint_b_outer_al_delta: float = 0.1524

    # Rear mechanical parts and calibration-source holder.
    sipm_lid_rmin: float = 5.08
    sipm_lid_rmax: float = 60.325
    sipm_lid_z: float = 8.89
    sipm_source_rmax: float = 6.325001
    sipm_source_z: float = 8.89
    sipm_lid_pos_z: float = -44.3708
    scint_holder_xy: float = 27.2
    scint_holder_z: float = 30.4434
    scint_holder_cutout_20_xy: float = 21.3
    scint_holder_cutout_20_z: float = 50.4434
    scint_holder_cutout_10_xy: float = 25.5
    scint_holder_cutout_10_z: float = 29.4434
    scint_holder_pos_z: float = -23.1291
    scint_holder_base_rmin: float = 10.0
    scint_holder_base_rmax: float = 22.86
    scint_holder_base_z: float = 1.575
    scint_holder_base_pos_z: float = -39.1383
    scint_sleeve_xy: float = 25.5
    scint_sleeve_z: float = 10.0
    scint_sleeve_cutout_xy: float = 18.5
    scint_sleeve_cutout_z: float = 20.0
    scint_sleeve_pos_z: float = -19.7122
    source_alo_rmin: float = 4.694999
    source_alo_rmax: float = 6.325
    source_alo_z: float = 3.1749
    source_alo_pos_z: float = -42.719651
    source_al_rmax: float = 6.325
    source_al_z: float = 12.7
    source_al_cutout_rmax: float = 3.1623
    source_al_cutout_z: float = 6.3498
    source_al_cutout_shift_z: float = -4.76255
    source_al_pos_z: float = -50.657102
    source_my_rmax: float = 4.694999
    source_my_z: float = 0.0064
    source_my_pos_z: float = -42.71965
    source_cal_rmax: float = 4.694999
    source_cal_z: float = 0.001
    source_cal_pos_z: float = -44.3066


PARAMS = GeometryParams()
DEGREES = "deg"


def build_materials(registry: g4.Registry) -> dict[str, object]:
    """Create the standard and experiment-specific materials."""
    materials = {
        name: g4.MaterialPredefined(name, registry)
        for name in (
            "G4_Al",
            "G4_PLASTIC_SC_VINYLTOLUENE",
            "G4_MYLAR",
            "G4_Cu",
            "G4_AIR",
        )
    }
    elements = {
        name: g4.ElementSimple(name, symbol, atomic_number, atomic_mass, registry)
        for name, symbol, atomic_number, atomic_mass in (
            ("Hydrogen", "H", 1, 1.008),
            ("Oxygen", "O", 8, 15.999),
            ("Silicon", "Si", 14, 28.085),
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
    materials["Vetronite"] = mixture(
        "Vetronite",
        2.0,
        ((0.467434920603219, "Silicon"), (0.532565079396781, "Oxygen")),
        "solid",
    )
    return materials


def box(name: str, xyz: tuple[float, float, float], registry: g4.Registry) -> object:
    """Create a box with the full x, y, and z lengths used by GDML."""
    return g4.solid.Box(name, *xyz, registry)


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
    """Create all primitive and boolean solids used by the detector."""
    p = params
    solids: dict[str, object] = {}

    def add(solid: object) -> object:
        solids[solid.name] = solid
        return solid

    # The T-pipe is the union of a vertical tube and a perpendicular branch.
    pipe_1i = add(tube("Pipe1i", 0.0, p.pipe_inner_radius, p.pipe_inner_major_length, registry))
    pipe_2i = add(tube("Pipe2i", 0.0, p.pipe_inner_radius, p.pipe_inner_minor_length, registry))
    inner_pipe = add(g4.solid.Union("InnerTPipe", pipe_1i, pipe_2i, transform((p.pipe_inner_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))
    pipe_1o = add(tube("Pipe1o", 0.0, p.pipe_outer_radius, p.pipe_outer_major_length, registry))
    pipe_2o = add(tube("Pipe2o", 0.0, p.pipe_outer_radius, p.pipe_outer_minor_length, registry))
    outer_pipe = add(g4.solid.Union("OuterTPipe", pipe_1o, pipe_2o, transform((p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))
    inner_main = add(tube("Pipe1i_tol", 0.0, p.pipe_inner_radius + 1e-6, 184.15, registry))
    steel_pipe = add(g4.solid.Subtraction("TDecayVolume", outer_pipe, inner_main, transform(), registry))
    inner_branch = add(tube("Pipe2i_tol", 0.0, p.pipe_inner_radius + 1e-6, p.pipe_inner_minor_length, registry))
    add(g4.solid.Subtraction("TDecayVolume2", steel_pipe, inner_branch, transform((p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0)), registry))

    add(tube("Flange", p.flange_inner_radius, p.flange_outer_radius, p.flange_thickness, registry))
    add(tube("VacuumWindowDisk", 0.0, p.vacuum_window_radius, p.vacuum_window_thickness, registry))
    rmax, rmin, z = zip(*p.wedge_zplanes)
    wedge = add(g4.solid.Polycone("WedgeRevolvedRing", 0.0, 360.0, z, rmin, rmax, registry, aunit=DEGREES))
    # The wedge is placed at the world-origin window interface. Carve its
    # footprint from the rotated vacuum T-pipe so the two world siblings meet
    # at a boundary rather than overlapping.
    add(
        g4.solid.Subtraction(
            "VacuumWithWedgeCutout",
            inner_pipe,
            wedge,
            transform((-p.world_pipe_z, 0.0, 0.0), (0.0, 90.0, 0.0)),
            registry,
        )
    )

    core = add(box("ScintCore", (p.scint_a_xy, p.scint_a_xy, p.scint_a_z), registry))
    inner_geom = add(
        box(
            "ScintInnerGeom",
            (
                p.scint_a_xy + p.scint_a_inner_al_delta,
                p.scint_a_xy + p.scint_a_inner_al_delta,
                p.scint_a_z + p.scint_a_inner_al_delta,
            ),
            registry,
        )
    )
    add(g4.solid.Subtraction("ScintInnerAl", inner_geom, core, transform(), registry))
    mylar_geom = add(
        box(
            "ScintMylarGeom",
            (
                p.scint_a_xy + p.scint_a_mylar_delta,
                p.scint_a_xy + p.scint_a_mylar_delta,
                p.scint_a_z + p.scint_a_mylar_delta,
            ),
            registry,
        )
    )
    add(g4.solid.Subtraction("ScintMylar", mylar_geom, inner_geom, transform(), registry))
    outer_geom = add(
        box(
            "ScintOuterGeom",
            (
                p.scint_a_xy + p.scint_a_outer_al_delta,
                p.scint_a_xy + p.scint_a_outer_al_delta,
                p.scint_a_z + p.scint_a_outer_al_delta,
            ),
            registry,
        )
    )
    add(g4.solid.Subtraction("ScintOuterAl", outer_geom, mylar_geom, transform(), registry))

    lid_base = add(tube("sipmLidBase", p.sipm_lid_rmin, p.sipm_lid_rmax, p.sipm_lid_z, registry))
    lid_source = add(tube("sipmLidSource", 0.0, p.sipm_source_rmax, p.sipm_source_z, registry))
    add(g4.solid.Subtraction("sipmLid", lid_base, lid_source, transform((0.0, 0.0, -1.2064)), registry))
    holder_base = add(box("scintHolderbase", (p.scint_holder_xy, p.scint_holder_xy, p.scint_holder_z), registry))
    holder_cutout_20 = add(box("scintHoldercutout20", (p.scint_holder_cutout_20_xy, p.scint_holder_cutout_20_xy, p.scint_holder_cutout_20_z), registry))
    holder_10 = add(g4.solid.Subtraction("scintHolder10", holder_base, holder_cutout_20, transform(), registry))
    holder_cutout_10 = add(box("scintHoldercutout10", (p.scint_holder_cutout_10_xy, p.scint_holder_cutout_10_xy, p.scint_holder_cutout_10_z), registry))
    add(g4.solid.Subtraction("scintHolder", holder_10, holder_cutout_10, transform(), registry))
    add(tube("scintHolderbaseRing", p.scint_holder_base_rmin, p.scint_holder_base_rmax, p.scint_holder_base_z, registry))
    sleeve_base = add(box("scintSleevebase", (p.scint_sleeve_xy, p.scint_sleeve_xy, p.scint_sleeve_z), registry))
    sleeve_cutout = add(box("scintHoldercutout20b", (p.scint_sleeve_cutout_xy, p.scint_sleeve_cutout_xy, p.scint_sleeve_cutout_z), registry))
    add(g4.solid.Subtraction("scintSleeve", sleeve_base, sleeve_cutout, transform(), registry))
    add(tube("SourceAlo", p.source_alo_rmin, p.source_alo_rmax, p.source_alo_z, registry))
    source_al_1 = add(tube("SourceAl1", 0.0, p.source_al_rmax, p.source_al_z, registry))
    source_al_2 = add(tube("SourceAl2", 0.0, p.source_al_cutout_rmax, p.source_al_cutout_z, registry))
    add(g4.solid.Subtraction("SourceAl", source_al_1, source_al_2, transform((0.0, 0.0, p.source_al_cutout_shift_z)), registry))
    add(tube("SourceMy", 0.0, p.source_my_rmax, p.source_my_z, registry))
    add(tube("SourceCal", 0.0, p.source_cal_rmax, p.source_cal_z, registry))
    add(box("World", p.world_xyz, registry))
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
) -> g4.LogicalVolume:
    """Create a logical volume and place it directly in the world."""
    volume = g4.LogicalVolume(solid, material, name, registry)
    g4.PhysicalVolume(list(rotation) + [DEGREES], position, volume, name, world, registry, copyNumber=copy_number)
    return volume


def build_geometry(output_path: str | Path, params: GeometryParams = PARAMS) -> Path:
    """Build and write the complete GDML geometry to ``output_path``."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    registry = g4.Registry()
    materials = build_materials(registry)
    solids = build_solids(params, registry)
    world = g4.LogicalVolume(solids["World"], materials["G4_AIR"], "World", registry)
    p = params

    placements = (
        ("Vacuum", "interGalactic", "VacuumWithWedgeCutout", (0.0, 0.0, p.world_pipe_z), (0.0, 90.0, 0.0), 1),
        ("SteelTPipe", "Stainless_Steel", "TDecayVolume2", (0.0, 0.0, p.world_pipe_z), (0.0, 90.0, 0.0), 0),
        # ("tPipeFlangeScint", "Stainless_Steel", "Flange", (0.0, 0.0, p.flange_scint_z), (0.0, 0.0, 0.0), 0),
        ("tPipeFlange1", "Stainless_Steel", "Flange", (0.0, 0.0, p.flange_1_z), (0.0, 0.0, 0.0), 0),
        ("tPipeFlange2", "Stainless_Steel", "Flange", (p.flange_side_offset_x, 0.0, p.flange_side_z), (0.0, 90.0, 0.0), 0),
        ("tPipeFlange3", "Stainless_Steel", "Flange", (-p.flange_side_offset_x, 0.0, p.flange_side_z), (0.0, 90.0, 0.0), 0),
        ("VacuumWindow", "G4_Al", "VacuumWindowDisk", (0.0, 0.0, p.vacuum_window_z), (0.0, 0.0, 0.0), 2),
        ("WedgeRevolvedRingLV", "Stainless_Steel", "WedgeRevolvedRing", (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0),
        ("AScintillator", "G4_PLASTIC_SC_VINYLTOLUENE", "ScintCore", (0.0, 0.0, p.scint_a_pos_z), (0.0, 0.0, 0.0), 3),
        ("Ali_sq10", "G4_Al", "ScintInnerAl", (0.0, 0.0, p.scint_a_pos_z), (0.0, 0.0, 0.0), 0),
        ("Mylar_sq10", "G4_MYLAR", "ScintMylar", (0.0, 0.0, p.scint_a_pos_z), (0.0, 0.0, 0.0), 0),
        ("Alo_sq10", "G4_Al", "ScintOuterAl", (0.0, 0.0, p.scint_a_pos_z), (0.0, 0.0, 0.0), 0),
        ("BScintillatorLV", "G4_PLASTIC_SC_VINYLTOLUENE", "ScintCore", (0.0, 0.0, p.scint_b_pos_z), (0.0, 0.0, 0.0), 4),
        ("Ali_sq20", "G4_Al", "ScintInnerAl", (0.0, 0.0, p.scint_b_pos_z), (0.0, 0.0, 0.0), 0),
        ("Mylar_sq20", "G4_MYLAR", "ScintMylar", (0.0, 0.0, p.scint_b_pos_z), (0.0, 0.0, 0.0), 0),
        ("Alo_sq20", "G4_Al", "ScintOuterAl", (0.0, 0.0, p.scint_b_pos_z), (0.0, 0.0, 0.0), 0),
        ("sipm_Lid", "G4_Cu", "sipmLid", (0.0, 0.0, p.sipm_lid_pos_z), (0.0, 0.0, 0.0), 0),
        ("scint_holder", "Vetronite", "scintHolder", (0.0, 0.0, p.scint_holder_pos_z), (0.0, 0.0, 0.0), 0),
        ("scint_holderBase", "Vetronite", "scintHolderbaseRing", (0.0, 0.0, p.scint_holder_base_pos_z), (0.0, 0.0, 0.0), 0),
        ("scint_Sleeve", "G4_Al", "scintSleeve", (0.0, 0.0, p.scint_sleeve_pos_z), (0.0, 0.0, 0.0), 0),
        ("SourceAlo", "G4_Al", "SourceAlo", (0.0, 0.0, p.source_alo_pos_z), (0.0, 0.0, 0.0), 0),
        ("SourceAl", "G4_Al", "SourceAl", (0.0, 0.0, p.source_al_pos_z), (0.0, 0.0, 0.0), 0),
        ("SourceMy", "G4_MYLAR", "SourceMy", (0.0, 0.0, p.source_my_pos_z), (0.0, 0.0, 0.0), 0),
        ("SourceCal", "G4_Al", "SourceCal", (0.0, 0.0, p.source_cal_pos_z), (0.0, 0.0, 0.0), 1),
    )
    logical_volumes = {
        name: place(name, materials[material], solids[solid], world, registry, position, rotation, copy_number)
        for name, material, solid, position, rotation, copy_number in placements
    }

    registry.setWorld(world.name)
    writer = gdml.Writer()
    writer.addDetector(registry)
    writer.write(path)
    return path


def main(argv: Sequence[str] | None = None) -> int:
    """Parse CLI arguments and write the requested GDML export."""
    parser = argparse.ArgumentParser(description="Generate geometry_export.gdml for BetaMonitor.")
    parser.add_argument("output", nargs="?", default="./dat/geometry_export_new.gdml", help="Output GDML path.")
    args = parser.parse_args(argv)
    print(build_geometry(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
