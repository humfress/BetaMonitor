"""Generate geometry_export.gdml for BetaMonitor.

Scintillator A and B dimensions are independently configurable.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import argparse
import xml.etree.ElementTree as ET


GDML_XSI = "http://www.w3.org/2001/XMLSchema-instance"
GDML_SCHEMA = "http://cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd"


@dataclass(frozen=True)
class GeometryParams:
    world_xyz: tuple[float, float, float] = (400.0, 400.0, 450.0)

    # Pipe and flange geometry.
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

    # Scintillator A placement and independent dimensions.
    scint_a_pos_z: float = -9.9836
    scint_a_xy: float = 20.0
    scint_a_z: float = 3.0
    scint_a_inner_al_delta: float = 0.0002
    scint_a_mylar_delta: float = 0.1522
    scint_a_outer_al_delta: float = 0.1524

    # Scintillator B placement and independent dimensions.
    scint_b_pos_z: float = -13.136
    scint_b_xy: float = 20.0
    scint_b_z: float = 3.0
    scint_b_inner_al_delta: float = 0.0002
    scint_b_mylar_delta: float = 0.1522
    scint_b_outer_al_delta: float = 0.1524

    # Rear mechanics and source pieces.
    sipm_lid_rmin: float = 5.08
    sipm_lid_rmax: float = 60.325
    sipm_lid_z: float = 8.89
    sipm_source_rmax: float = 6.325001
    sipm_source_z: float = 8.89
    sipm_lid_pos_z: float = -44.3708

    scint_holder_xy: float = 22.2
    scint_holder_z: float = 30.4434
    scint_holder_cutout_20_xy: float = 16.3
    scint_holder_cutout_20_z: float = 50.4434
    scint_holder_cutout_10_xy: float = 20.5
    scint_holder_cutout_10_z: float = 29.4434
    scint_holder_pos_z: float = -23.1291

    scint_holder_base_rmin: float = 10.0
    scint_holder_base_rmax: float = 22.86
    scint_holder_base_z: float = 1.575
    scint_holder_base_pos_z: float = -39.1383

    scint_sleeve_xy: float = 20.5
    scint_sleeve_z: float = 10.0
    scint_sleeve_cutout_xy: float = 13.5
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


@dataclass(frozen=True)
class MaterialDef:
    name: str
    state: str
    density: float
    mee: float
    temperature: float = 293.15
    pressure: float | None = None
    atomic_number: int | None = None
    atomic_mass: float | None = None
    fractions: tuple[tuple[float, str], ...] = ()


@dataclass(frozen=True)
class VolumeSpec:
    name: str
    material: str
    solid: str
    position: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation: tuple[float, float, float] | None = None
    copy_number: int | None = None


def _fmt(value: float) -> str:
    text = f"{float(value):.12g}"
    return "0" if text == "-0" else text


def _sub(parent: ET.Element, tag: str, **attrs) -> ET.Element:
    element = ET.SubElement(parent, tag)
    for key, value in attrs.items():
        if isinstance(value, float):
            element.set(key, _fmt(value))
        else:
            element.set(key, str(value))
    return element


def _position(parent: ET.Element, name: str, xyz: tuple[float, float, float]) -> ET.Element:
    return _sub(parent, "position", name=name, unit="mm", x=xyz[0], y=xyz[1], z=xyz[2])


def _rotation(parent: ET.Element, name: str, xyz: tuple[float, float, float]) -> ET.Element:
    return _sub(parent, "rotation", name=name, unit="deg", x=xyz[0], y=xyz[1], z=xyz[2])


def volume_spec(
    name: str,
    material: str,
    solid: str,
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
    rotation: tuple[float, float, float] | None = None,
    *,
    copy_number: int | None = None,
) -> VolumeSpec:
    return VolumeSpec(
        name=name,
        material=material,
        solid=solid,
        position=position,
        rotation=rotation,
        copy_number=copy_number,
    )


class GeometryExportBuilder:
    def __init__(self, params: GeometryParams):
        self.p = params

    def write(self, output_path: str | Path) -> Path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tree = ET.ElementTree(self.build())
        ET.indent(tree, space="  ")
        tree.write(path, encoding="UTF-8", xml_declaration=True)
        return path

    def build(self) -> ET.Element:
        root = ET.Element(
            "gdml",
            {
                "xmlns:xsi": GDML_XSI,
                "xsi:noNamespaceSchemaLocation": GDML_SCHEMA,
            },
        )
        ET.SubElement(root, "define")
        self._emit_materials(root)
        solids = ET.SubElement(root, "solids")
        self._emit_solids(solids)
        self._emit_structure(root)
        setup = _sub(root, "setup", name="Default", version="1.0")
        _sub(setup, "world", ref="World")
        return root

    def _emit_materials(self, root: ET.Element) -> None:
        materials = ET.SubElement(root, "materials")

        def isotope(name: str, n: int, z: int, mass: float) -> None:
            node = _sub(materials, "isotope", name=name, N=n, Z=z)
            _sub(node, "atom", unit="g/mole", value=mass)

        def element(name: str, fractions: tuple[tuple[float, str], ...]) -> None:
            node = _sub(materials, "element", name=name)
            for frac, ref in fractions:
                _sub(node, "fraction", n=frac, ref=ref)

        isotopes = (
            ("H1", 1, 1, 1.00782503081372),
            ("H2", 2, 1, 2.01410199966617),
            ("C12", 12, 6, 12.0),
            ("C13", 13, 6, 13.0034),
            ("N14", 14, 7, 14.0031),
            ("N15", 15, 7, 15.0001),
            ("O16", 16, 8, 15.9949),
            ("O17", 17, 8, 16.9991),
            ("O18", 18, 8, 17.9992),
            ("Ar36", 36, 18, 35.9675),
            ("Ar38", 38, 18, 37.9627),
            ("Ar40", 40, 18, 39.9624),
            ("Fe54", 54, 26, 53.9396),
            ("Fe56", 56, 26, 55.9349),
            ("Fe57", 57, 26, 56.9354),
            ("Fe58", 58, 26, 57.9333),
            ("Cr50", 50, 24, 49.946),
            ("Cr52", 52, 24, 51.9405),
            ("Cr53", 53, 24, 52.9407),
            ("Cr54", 54, 24, 53.9389),
            ("Ni58", 58, 28, 57.9353),
            ("Ni60", 60, 28, 59.9308),
            ("Ni61", 61, 28, 60.9311),
            ("Ni62", 62, 28, 61.9283),
            ("Ni64", 64, 28, 63.928),
            ("Mn55", 55, 25, 54.938),
            ("Cu63", 63, 29, 62.9296),
            ("Cu65", 65, 29, 64.9278),
            ("Si28", 28, 14, 27.9769),
            ("Si29", 29, 14, 28.9765),
            ("Si30", 30, 14, 29.9738),
        )
        for item in isotopes:
            isotope(*item)

        element("H", ((0.999885, "H1"), (0.000115, "H2")))
        element("C", ((0.9893, "C12"), (0.0107, "C13")))
        element("N", ((0.99632, "N14"), (0.00368, "N15")))
        element("O", ((0.99757, "O16"), (0.00038, "O17"), (0.00205, "O18")))
        element("Ar", ((0.003365, "Ar36"), (0.000632, "Ar38"), (0.996003, "Ar40")))
        element("Fe", ((0.05845, "Fe54"), (0.91754, "Fe56"), (0.02119, "Fe57"), (0.00282, "Fe58")))
        element("Cr", ((0.04345, "Cr50"), (0.83789, "Cr52"), (0.09501, "Cr53"), (0.02365, "Cr54")))
        element("Ni", ((0.680769, "Ni58"), (0.262231, "Ni60"), (0.011399, "Ni61"), (0.036345, "Ni62"), (0.009256, "Ni64")))
        element("Mn", ((1.0, "Mn55"),))
        element("Cu", ((0.6917, "Cu63"), (0.3083, "Cu65")))
        element("Silicon", ((0.922296077703922, "Si28"), (0.0468319531680468, "Si29"), (0.0308719691280309, "Si30")))

        materials_def = (
            MaterialDef("interGalactic", "gas", 1e-25, 19.2, temperature=2.73, pressure=1.3332e-08, fractions=((1.0, "H"),)),
            MaterialDef("Stainless_Steel", "solid", 8.03, 282.252009421111, fractions=((0.72, "Fe"), (0.18, "Cr"), (0.08, "Ni"), (0.02, "Mn"))),
            MaterialDef("G4_Al", "solid", 2.699, 166.0, atomic_number=13, atomic_mass=26.9815),
            MaterialDef("G4_PLASTIC_SC_VINYLTOLUENE", "solid", 1.032, 64.7, fractions=((0.914708531800025, "C"), (0.0852914681999746, "H"))),
            MaterialDef("G4_MYLAR", "solid", 1.4, 78.7, fractions=((0.625010832340889, "C"), (0.0419607170679432, "H"), (0.333028450591168, "O"))),
            MaterialDef("G4_Cu", "solid", 8.96, 322.0, fractions=((1.0, "Cu"),)),
            MaterialDef("Vetronite", "solid", 2.0, 125.663004061076, fractions=((0.467434920603219, "Silicon"), (0.532565079396781, "O"))),
            MaterialDef("G4_AIR", "gas", 0.00120479, 85.7, fractions=((0.000124000124000124, "C"), (0.755267755267755, "N"), (0.231781231781232, "O"), (0.0128270128270128, "Ar"))),
        )

        for m in materials_def:
            attrs: dict[str, object] = {"name": m.name, "state": m.state}
            if m.atomic_number is not None:
                attrs["Z"] = m.atomic_number
            node = _sub(materials, "material", **attrs)
            _sub(node, "T", unit="K", value=m.temperature)
            if m.pressure is not None:
                _sub(node, "P", unit="pascal", value=m.pressure)
            _sub(node, "MEE", unit="eV", value=m.mee)
            _sub(node, "D", unit="g/cm3", value=m.density)
            if m.atomic_mass is not None:
                _sub(node, "atom", unit="g/mole", value=m.atomic_mass)
            for frac, ref in m.fractions:
                _sub(node, "fraction", n=frac, ref=ref)

    def _emit_tube(self, solids: ET.Element, name: str, rmin: float, rmax: float, z: float) -> None:
        _sub(solids, "tube", name=name, aunit="deg", deltaphi=360, lunit="mm", rmin=rmin, rmax=rmax, startphi=0, z=z)

    def _emit_box(self, solids: ET.Element, name: str, x: float, y: float, z: float) -> None:
        _sub(solids, "box", name=name, lunit="mm", x=x, y=y, z=z)

    def _emit_union(self, solids: ET.Element, name: str, first: str, second: str, pos: tuple[float, float, float], rot: tuple[float, float, float] | None = None) -> None:
        node = _sub(solids, "union", name=name)
        _sub(node, "first", ref=first)
        _sub(node, "second", ref=second)
        _position(node, f"{name}_pos", pos)
        if rot is not None:
            _rotation(node, f"{name}_rot", rot)

    def _emit_subtraction(self, solids: ET.Element, name: str, first: str, second: str, pos: tuple[float, float, float] | None = None, rot: tuple[float, float, float] | None = None) -> None:
        node = _sub(solids, "subtraction", name=name)
        _sub(node, "first", ref=first)
        _sub(node, "second", ref=second)
        if pos is not None:
            _position(node, f"{name}_pos", pos)
        if rot is not None:
            _rotation(node, f"{name}_rot", rot)

    def _emit_scintillator_stack_box(
        self,
        solids: ET.Element,
        prefix: str,
        core_xy: float,
        core_z: float,
        inner_al_delta: float,
        mylar_delta: float,
        outer_al_delta: float,
    ) -> None:
        core = f"{prefix}Core"
        inner_geom = f"{prefix}InnerGeom"
        inner = f"{prefix}InnerAl"
        mylar_geom = f"{prefix}MylarGeom"
        mylar = f"{prefix}Mylar"
        outer_geom = f"{prefix}OuterGeom"
        outer = f"{prefix}OuterAl"

        self._emit_box(solids, core, core_xy, core_xy, core_z)
        self._emit_box(solids, inner_geom, core_xy + inner_al_delta, core_xy + inner_al_delta, core_z + inner_al_delta)
        self._emit_subtraction(solids, inner, inner_geom, core)
        self._emit_box(solids, mylar_geom, core_xy + mylar_delta, core_xy + mylar_delta, core_z + mylar_delta)
        self._emit_subtraction(solids, mylar, mylar_geom, inner_geom)
        self._emit_box(solids, outer_geom, core_xy + outer_al_delta, core_xy + outer_al_delta, core_z + outer_al_delta)
        self._emit_subtraction(solids, outer, outer_geom, mylar_geom)

    def _emit_solids(self, solids: ET.Element) -> None:
        p = self.p

        self._emit_tube(solids, "Pipe1i", 0.0, p.pipe_inner_radius, p.pipe_inner_major_length)
        self._emit_tube(solids, "Pipe2i", 0.0, p.pipe_inner_radius, p.pipe_inner_minor_length)
        self._emit_union(solids, "InnerTPipe", "Pipe1i", "Pipe2i", (p.pipe_inner_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0))
        self._emit_tube(solids, "Pipe1o", 0.0, p.pipe_outer_radius, p.pipe_outer_major_length)
        self._emit_tube(solids, "Pipe2o", 0.0, p.pipe_outer_radius, p.pipe_outer_minor_length)
        self._emit_union(solids, "OuterTPipe", "Pipe1o", "Pipe2o", (p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0))
        self._emit_tube(solids, "Pipe1i_tol", 0.0, p.pipe_inner_radius + 1e-6, 184.15)
        self._emit_subtraction(solids, "TDecayVolume", "OuterTPipe", "Pipe1i_tol")
        self._emit_tube(solids, "Pipe2i_tol", 0.0, p.pipe_inner_radius + 1e-6, p.pipe_inner_minor_length)
        self._emit_subtraction(solids, "TDecayVolume2", "TDecayVolume", "Pipe2i_tol", (p.pipe_outer_offset_x, 0.0, 0.0), (0.0, -90.0, 0.0))

        self._emit_tube(solids, "FlangeScint", 10.0, p.flange_outer_radius, p.flange_thickness)
        self._emit_tube(solids, "Flange", p.flange_inner_radius, p.flange_outer_radius, p.flange_thickness)
        self._emit_tube(solids, "VacuumWindowDisk", 0.0, p.vacuum_window_radius, p.vacuum_window_thickness)

        poly = _sub(solids, "polycone", name="WedgeRevolvedRing", aunit="deg", deltaphi=360, lunit="mm", startphi=0)
        for rmax, rmin, z in p.wedge_zplanes:
            _sub(poly, "zplane", rmax=rmax, rmin=rmin, z=z)

        self._emit_scintillator_stack_box(
            solids,
            prefix="AScint",
            core_xy=p.scint_a_xy,
            core_z=p.scint_a_z,
            inner_al_delta=p.scint_a_inner_al_delta,
            mylar_delta=p.scint_a_mylar_delta,
            outer_al_delta=p.scint_a_outer_al_delta,
        )
        self._emit_scintillator_stack_box(
            solids,
            prefix="BScint",
            core_xy=p.scint_b_xy,
            core_z=p.scint_b_z,
            inner_al_delta=p.scint_b_inner_al_delta,
            mylar_delta=p.scint_b_mylar_delta,
            outer_al_delta=p.scint_b_outer_al_delta,
        )

        self._emit_tube(solids, "sipmLidBase", p.sipm_lid_rmin, p.sipm_lid_rmax, p.sipm_lid_z)
        self._emit_tube(solids, "sipmLidSource", 0.0, p.sipm_source_rmax, p.sipm_source_z)
        self._emit_subtraction(solids, "sipmLid", "sipmLidBase", "sipmLidSource", (0.0, 0.0, -1.2064))

        self._emit_box(solids, "scintHolderbase", p.scint_holder_xy, p.scint_holder_xy, p.scint_holder_z)
        self._emit_box(solids, "scintHoldercutout20", p.scint_holder_cutout_20_xy, p.scint_holder_cutout_20_xy, p.scint_holder_cutout_20_z)
        self._emit_subtraction(solids, "scintHolder10", "scintHolderbase", "scintHoldercutout20")
        self._emit_box(solids, "scintHoldercutout10", p.scint_holder_cutout_10_xy, p.scint_holder_cutout_10_xy, p.scint_holder_cutout_10_z)
        self._emit_subtraction(solids, "scintHolder", "scintHolder10", "scintHoldercutout10")

        self._emit_tube(solids, "scintHolderbaseRing", p.scint_holder_base_rmin, p.scint_holder_base_rmax, p.scint_holder_base_z)
        self._emit_box(solids, "scintSleevebase", p.scint_sleeve_xy, p.scint_sleeve_xy, p.scint_sleeve_z)
        self._emit_box(solids, "scintHoldercutout20b", p.scint_sleeve_cutout_xy, p.scint_sleeve_cutout_xy, p.scint_sleeve_cutout_z)
        self._emit_subtraction(solids, "scintSleeve", "scintSleevebase", "scintHoldercutout20b")

        self._emit_tube(solids, "SourceAlo", p.source_alo_rmin, p.source_alo_rmax, p.source_alo_z)
        self._emit_tube(solids, "SourceAl1", 0.0, p.source_al_rmax, p.source_al_z)
        self._emit_tube(solids, "SourceAl2", 0.0, p.source_al_cutout_rmax, p.source_al_cutout_z)
        self._emit_subtraction(solids, "SourceAl", "SourceAl1", "SourceAl2", (0.0, 0.0, p.source_al_cutout_shift_z))
        self._emit_tube(solids, "SourceMy", 0.0, p.source_my_rmax, p.source_my_z)
        self._emit_tube(solids, "SourceCal", 0.0, p.source_cal_rmax, p.source_cal_z)

        self._emit_box(solids, "World", *p.world_xyz)

    def _volume_specs(self) -> tuple[VolumeSpec, ...]:
        p = self.p
        return (
            volume_spec("World", "G4_AIR", "World"),
            volume_spec("Vacuum", "interGalactic", "InnerTPipe", (0.0, 0.0, p.world_pipe_z), (0.0, 90.0, 0.0), copy_number=1),
            volume_spec("SteelTPipe", "Stainless_Steel", "TDecayVolume2", (0.0, 0.0, p.world_pipe_z), (0.0, 90.0, 0.0)),
            volume_spec("tPipeFlangeScint", "Stainless_Steel", "FlangeScint", (0.0, 0.0, p.flange_scint_z)),
            volume_spec("tPipeFlange1", "Stainless_Steel", "Flange", (0.0, 0.0, p.flange_1_z)),
            volume_spec("tPipeFlange2", "Stainless_Steel", "Flange", (p.flange_side_offset_x, 0.0, p.flange_side_z), (0.0, 90.0, 0.0)),
            volume_spec("tPipeFlange3", "Stainless_Steel", "Flange", (-p.flange_side_offset_x, 0.0, p.flange_side_z), (0.0, 90.0, 0.0)),
            volume_spec("VacuumWindow", "G4_Al", "VacuumWindowDisk", (0.0, 0.0, p.vacuum_window_z), copy_number=2),
            volume_spec("WedgeRevolvedRingLV", "Stainless_Steel", "WedgeRevolvedRing"),
            volume_spec("AScintillator", "G4_PLASTIC_SC_VINYLTOLUENE", "AScintCore", (0.0, 0.0, p.scint_a_pos_z), copy_number=3),
            volume_spec("Ali_sq10", "G4_Al", "AScintInnerAl", (0.0, 0.0, p.scint_a_pos_z)),
            volume_spec("Mylar_sq10", "G4_MYLAR", "AScintMylar", (0.0, 0.0, p.scint_a_pos_z)),
            volume_spec("Alo_sq10", "G4_Al", "AScintOuterAl", (0.0, 0.0, p.scint_a_pos_z)),
            volume_spec("BScintillatorLV", "G4_PLASTIC_SC_VINYLTOLUENE", "BScintCore", (0.0, 0.0, p.scint_b_pos_z), copy_number=4),
            volume_spec("Ali_sq20", "G4_Al", "BScintInnerAl", (0.0, 0.0, p.scint_b_pos_z)),
            volume_spec("Mylar_sq20", "G4_MYLAR", "BScintMylar", (0.0, 0.0, p.scint_b_pos_z)),
            volume_spec("Alo_sq20", "G4_Al", "BScintOuterAl", (0.0, 0.0, p.scint_b_pos_z)),
            volume_spec("sipm_Lid", "G4_Cu", "sipmLid", (0.0, 0.0, p.sipm_lid_pos_z)),
            volume_spec("scint_holder", "Vetronite", "scintHolder", (0.0, 0.0, p.scint_holder_pos_z)),
            volume_spec("scint_holderBase", "Vetronite", "scintHolderbaseRing", (0.0, 0.0, p.scint_holder_base_pos_z)),
            volume_spec("scint_Sleeve", "G4_Al", "scintSleeve", (0.0, 0.0, p.scint_sleeve_pos_z)),
            volume_spec("SourceAlo", "G4_Al", "SourceAlo", (0.0, 0.0, p.source_alo_pos_z)),
            volume_spec("SourceAl", "G4_Al", "SourceAl", (0.0, 0.0, p.source_al_pos_z)),
            volume_spec("SourceMy", "G4_MYLAR", "SourceMy", (0.0, 0.0, p.source_my_pos_z)),
            volume_spec("SourceCal", "G4_Al", "SourceCal", (0.0, 0.0, p.source_cal_pos_z), copy_number=1),
        )

    def _emit_structure(self, root: ET.Element) -> None:
        structure = ET.SubElement(root, "structure")
        specs = self._volume_specs()

        nodes: dict[str, ET.Element] = {}
        # Geant4 GDML reader expects referenced volumes to already exist,
        # so emit all non-world volumes before the world volume.
        ordered_specs = tuple(spec for spec in specs if spec.name != "World") + tuple(spec for spec in specs if spec.name == "World")
        for spec in ordered_specs:
            node = _sub(structure, "volume", name=spec.name)
            _sub(node, "materialref", ref=spec.material)
            _sub(node, "solidref", ref=spec.solid)
            nodes[spec.name] = node

        world = nodes["World"]
        for spec in specs:
            if spec.name == "World":
                continue
            attrs = {"name": spec.name}
            if spec.copy_number is not None:
                attrs["copynumber"] = str(spec.copy_number)
            physvol = _sub(world, "physvol", **attrs)
            _sub(physvol, "volumeref", ref=spec.name)
            _position(physvol, f"{spec.name}_pos", spec.position)
            if spec.rotation is not None:
                _rotation(physvol, f"{spec.name}_rot", spec.rotation)


def build_geometry(output_path: str | Path) -> Path:
    return GeometryExportBuilder(PARAMS).write(output_path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate geometry_export.gdml for BetaMonitor.")
    parser.add_argument("output", nargs="?", default="./dat/geometry_export_new.gdml", help="Output GDML path.")
    args = parser.parse_args(argv)

    path = build_geometry(args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
