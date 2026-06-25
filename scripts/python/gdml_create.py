"""Parametric GDML generator for BetaMonitor.

This module creates the two GDML geometries shipped with the repository and
keeps the logical volume names expected by the Geant4 detector construction:

- Vacuum
- VacuumWindow
- AScintillator
- BScintillatorLV

The generator is intentionally dependency-free so it can be used as a simple
script or imported from other Python tooling.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import xml.etree.ElementTree as ET


GDML_XSI = "http://www.w3.org/2001/XMLSchema-instance"
GDML_SCHEMA = "http://cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd"


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


def _position(parent: ET.Element, name: str, x: float, y: float, z: float) -> ET.Element:
	return _sub(parent, "position", name=name, unit="mm", x=x, y=y, z=z)


def _rotation(parent: ET.Element, name: str, x: float, y: float, z: float) -> ET.Element:
	return _sub(parent, "rotation", name=name, unit="deg", x=x, y=y, z=z)


@dataclass(frozen=True)
class GeometryConfig:
	"""Parameters for one BetaMonitor GDML variant."""

	name: str
	world: tuple[float, float, float]
	world_pipe_z: float
	vacuum_window_z: float
	scint_a_z: float
	scint_b_z: float
	use_telescope_parts: bool = False

	pipe_inner_radius: float = 17.399
	pipe_outer_radius: float = 19.05
	pipe_inner_major_length: float = 196.85
	pipe_inner_minor_length: float = 98.425
	pipe_outer_major_length: float = 171.45
	pipe_outer_minor_length: float = 85.725
	pipe_inner_offset_x: float = -49.2125
	pipe_outer_offset_x: float = -42.8625

	flange_inner_radius: float = 17.399
	flange_outer_radius: float = 34.544
	flange_thickness: float = 12.7
	flange_scint_z: float = 6.35
	flange_1_z: float = 19.05
	flange_side_offset_x: float = 92.075
	flange_side_z: float = 111.125

	source_alo_z: float = -42.719651
	source_al_z: float = -50.657102
	source_my_z: float = -42.71965
	source_cal_z: float = -44.3066

	# Base geometry defaults.
	base_vacuum_window_radius: float = 10.0
	base_vacuum_window_thickness: float = 0.05
	base_scint_radius: float = 25.0
	base_scint_a_thickness: float = 2.0
	base_scint_b_thickness: float = 50.0
	base_scint_a_al_thickness: float = 2.0002
	base_scint_b_al_thickness: float = 50.0002

	# Telescope geometry defaults.
	telescope_vacuum_window_radius: float = 17.399
	telescope_vacuum_window_thickness: float = 0.0762
	telescope_scint_xy: float = 20.0
	telescope_scint_z: float = 3.0
	telescope_scint_ma_xy: float = 20.0002
	telescope_scint_ma_z: float = 3.0002
	telescope_scint_m_xy: float = 20.1522
	telescope_scint_m_z: float = 3.1522
	telescope_scint_outer_xy: float = 20.1524
	telescope_scint_outer_z: float = 3.1524


class BetaMonitorGDMLBuilder:
	def __init__(self, config: GeometryConfig):
		self.config = config

	@classmethod
	def base(cls) -> "BetaMonitorGDMLBuilder":
		return cls(
			GeometryConfig(
				name="geometry_export",
				world=(400.0, 400.0, 450.0),
				world_pipe_z=111.125,
				vacuum_window_z=0.025,
				scint_a_z=-9.0001,
				scint_b_z=-38.0002,
				use_telescope_parts=False,
			)
		)

	@classmethod
	def telescope(cls) -> "BetaMonitorGDMLBuilder":
		return cls(
			GeometryConfig(
				name="geometry_telescope_export",
				world=(400.0, 400.0, 450.0),
				world_pipe_z=98.425,
				vacuum_window_z=-0.0381,
				scint_a_z=-9.9836,
				scint_b_z=-13.136,
				use_telescope_parts=True,
			)
		)

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
		self._materials(root)
		self._solids(root)
		self._structure(root)
		setup = _sub(root, "setup", name="Default", version="1.0")
		_sub(setup, "world", ref="World")
		return root

	def _materials(self, root: ET.Element) -> None:
		materials = ET.SubElement(root, "materials")

		# Elements used by multiple materials.
		def isotope(name: str, n: int, z: int, mass: float) -> None:
			elem = _sub(materials, "isotope", name=name, N=n, Z=z)
			_sub(elem, "atom", unit="g/mole", value=mass)

		def element(name: str, fractions: list[tuple[float, str]]) -> None:
			elem = _sub(materials, "element", name=name)
			for frac, ref in fractions:
				_sub(elem, "fraction", n=frac, ref=ref)

		isotope("H1", 1, 1, 1.00782503081372)
		isotope("H2", 2, 1, 2.01410199966617)
		element("H", [(0.999885, "H1"), (0.000115, "H2")])
		intergalactic = _sub(materials, "material", name="interGalactic", state="gas")
		_sub(intergalactic, "T", unit="K", value=2.73)
		_sub(intergalactic, "P", unit="pascal", value=1.3332e-08)
		_sub(intergalactic, "MEE", unit="eV", value=19.2)
		_sub(intergalactic, "D", unit="g/cm3", value=1e-25)
		_sub(intergalactic, "fraction", n=1.0, ref="H")

		isotope("C12", 12, 6, 12.0)
		isotope("C13", 13, 6, 13.0034)
		element("C", [(0.9893, "C12"), (0.0107, "C13")])
		isotope("N14", 14, 7, 14.0031)
		isotope("N15", 15, 7, 15.0001)
		element("N", [(0.99632, "N14"), (0.00368, "N15")])
		isotope("O16", 16, 8, 15.9949)
		isotope("O17", 17, 8, 16.9991)
		isotope("O18", 18, 8, 17.9992)
		element("O", [(0.99757, "O16"), (0.00038, "O17"), (0.00205, "O18")])
		isotope("Ar36", 36, 18, 35.9675)
		isotope("Ar38", 38, 18, 37.9627)
		isotope("Ar40", 40, 18, 39.9624)
		element("Ar", [(0.003365, "Ar36"), (0.000632, "Ar38"), (0.996003, "Ar40")])

		isotope("Fe54", 54, 26, 53.9396)
		isotope("Fe56", 56, 26, 55.9349)
		isotope("Fe57", 57, 26, 56.9354)
		isotope("Fe58", 58, 26, 57.9333)
		element("Fe", [(0.05845, "Fe54"), (0.91754, "Fe56"), (0.02119, "Fe57"), (0.00282, "Fe58")])

		isotope("Cr50", 50, 24, 49.946)
		isotope("Cr52", 52, 24, 51.9405)
		isotope("Cr53", 53, 24, 52.9407)
		isotope("Cr54", 54, 24, 53.9389)
		element("Cr", [(0.04345, "Cr50"), (0.83789, "Cr52"), (0.09501, "Cr53"), (0.02365, "Cr54")])

		isotope("Ni58", 58, 28, 57.9353)
		isotope("Ni60", 60, 28, 59.9308)
		isotope("Ni61", 61, 28, 60.9311)
		isotope("Ni62", 62, 28, 61.9283)
		isotope("Ni64", 64, 28, 63.928)
		element("Ni", [(0.680769, "Ni58"), (0.262231, "Ni60"), (0.011399, "Ni61"), (0.036345, "Ni62"), (0.009256, "Ni64")])

		isotope("Mn55", 55, 25, 54.938)
		element("Mn", [(1.0, "Mn55")])

		isotope("Cu63", 63, 29, 62.9296)
		isotope("Cu65", 65, 29, 64.9278)
		element("Cu", [(0.6917, "Cu63"), (0.3083, "Cu65")])

		isotope("Si28", 28, 14, 27.9769)
		isotope("Si29", 29, 14, 28.9765)
		isotope("Si30", 30, 14, 29.9738)
		element("Silicon", [(0.922296077703922, "Si28"), (0.0468319531680468, "Si29"), (0.0308719691280309, "Si30")])

		steel = _sub(materials, "material", name="Stainless_Steel", state="solid")
		_sub(steel, "T", unit="K", value=293.15)
		_sub(steel, "MEE", unit="eV", value=282.252009421111)
		_sub(steel, "D", unit="g/cm3", value=8.03)
		_sub(steel, "fraction", n=0.72, ref="Fe")
		_sub(steel, "fraction", n=0.18, ref="Cr")
		_sub(steel, "fraction", n=0.08, ref="Ni")
		_sub(steel, "fraction", n=0.02, ref="Mn")

		al = _sub(materials, "material", name="G4_Al", state="solid", Z=13)
		_sub(al, "T", unit="K", value=293.15)
		_sub(al, "MEE", unit="eV", value=166)
		_sub(al, "D", unit="g/cm3", value=2.699)
		_sub(al, "atom", unit="g/mole", value=26.9815)

		plastic = _sub(materials, "material", name="G4_PLASTIC_SC_VINYLTOLUENE", state="solid")
		_sub(plastic, "T", unit="K", value=293.15)
		_sub(plastic, "MEE", unit="eV", value=64.7)
		_sub(plastic, "D", unit="g/cm3", value=1.032)
		_sub(plastic, "fraction", n=0.914708531800025, ref="C")
		_sub(plastic, "fraction", n=0.0852914681999746, ref="H")

		mylar = _sub(materials, "material", name="G4_MYLAR", state="solid")
		_sub(mylar, "T", unit="K", value=293.15)
		_sub(mylar, "MEE", unit="eV", value=78.7)
		_sub(mylar, "D", unit="g/cm3", value=1.4)
		_sub(mylar, "fraction", n=0.625010832340889, ref="C")
		_sub(mylar, "fraction", n=0.0419607170679432, ref="H")
		_sub(mylar, "fraction", n=0.333028450591168, ref="O")

		cu = _sub(materials, "material", name="G4_Cu", state="solid")
		_sub(cu, "T", unit="K", value=293.15)
		_sub(cu, "MEE", unit="eV", value=322)
		_sub(cu, "D", unit="g/cm3", value=8.96)
		_sub(cu, "fraction", n=1.0, ref="Cu")

		vetronite = _sub(materials, "material", name="Vetronite", state="solid")
		_sub(vetronite, "T", unit="K", value=293.15)
		_sub(vetronite, "MEE", unit="eV", value=125.663004061076)
		_sub(vetronite, "D", unit="g/cm3", value=2.0)
		_sub(vetronite, "fraction", n=0.467434920603219, ref="Silicon")
		_sub(vetronite, "fraction", n=0.532565079396781, ref="O")

		air = _sub(materials, "material", name="G4_AIR", state="gas")
		_sub(air, "T", unit="K", value=293.15)
		_sub(air, "MEE", unit="eV", value=85.7)
		_sub(air, "D", unit="g/cm3", value=0.00120479)
		_sub(air, "fraction", n=0.000124000124000124, ref="C")
		_sub(air, "fraction", n=0.755267755267755, ref="N")
		_sub(air, "fraction", n=0.231781231781232, ref="O")
		_sub(air, "fraction", n=0.0128270128270128, ref="Ar")

	def _solids(self, root: ET.Element) -> None:
		solids = ET.SubElement(root, "solids")
		c = self.config

		_sub(solids, "tube", name="Pipe1i", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_inner_radius, rmin=0, startphi=0, z=c.pipe_inner_major_length)
		_sub(solids, "tube", name="Pipe2i", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_inner_radius, rmin=0, startphi=0, z=c.pipe_inner_minor_length)
		inner = _sub(solids, "union", name="InnerTPipe")
		_sub(inner, "first", ref="Pipe1i")
		_sub(inner, "second", ref="Pipe2i")
		_position(inner, "InnerTPipe_pos", c.pipe_inner_offset_x, 0, 0)
		_rotation(inner, "InnerTPipe_rot", 0, -90, 0)

		_sub(solids, "tube", name="Pipe1o", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_outer_radius, rmin=0, startphi=0, z=c.pipe_outer_major_length)
		_sub(solids, "tube", name="Pipe2o", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_outer_radius, rmin=0, startphi=0, z=c.pipe_outer_minor_length)
		outer = _sub(solids, "union", name="OuterTPipe")
		_sub(outer, "first", ref="Pipe1o")
		_sub(outer, "second", ref="Pipe2o")
		_position(outer, "OuterTPipe_pos", c.pipe_outer_offset_x, 0, 0)
		_rotation(outer, "OuterTPipe_rot", 0, -90, 0)

		_sub(solids, "tube", name="Pipe1i_tol", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_inner_radius + 1e-6, rmin=0, startphi=0, z=184.15)
		tdv = _sub(solids, "subtraction", name="TDecayVolume")
		_sub(tdv, "first", ref="OuterTPipe")
		_sub(tdv, "second", ref="Pipe1i_tol")
		_sub(solids, "tube", name="Pipe2i_tol", aunit="deg", deltaphi=360, lunit="mm", rmax=c.pipe_inner_radius + 1e-6, rmin=0, startphi=0, z=c.pipe_inner_minor_length)
		tdv2 = _sub(solids, "subtraction", name="TDecayVolume2")
		_sub(tdv2, "first", ref="TDecayVolume")
		_sub(tdv2, "second", ref="Pipe2i_tol")
		_position(tdv2, "TDecayVolume2_pos", c.pipe_outer_offset_x, 0, 0)
		_rotation(tdv2, "TDecayVolume2_rot", 0, -90, 0)

		_sub(solids, "tube", name="FlangeScint", aunit="deg", deltaphi=360, lunit="mm", rmax=c.flange_outer_radius, rmin=10, startphi=0, z=c.flange_thickness)
		_sub(solids, "tube", name="Flange", aunit="deg", deltaphi=360, lunit="mm", rmax=c.flange_outer_radius, rmin=c.flange_inner_radius, startphi=0, z=c.flange_thickness)
		_sub(solids, "tube", name="VacuumWindowDisk", aunit="deg", deltaphi=360, lunit="mm", rmax=c.telescope_vacuum_window_radius if c.use_telescope_parts else c.base_vacuum_window_radius, rmin=0, startphi=0, z=c.telescope_vacuum_window_thickness if c.use_telescope_parts else c.base_vacuum_window_thickness)

		if c.use_telescope_parts:
			self._telescope_solids(solids)
		else:
			self._base_scint_solids(solids)

		self._source_solids(solids)
		self._holder_solids(solids)
		_sub(solids, "box", name="World", lunit="mm", x=c.world[0], y=c.world[1], z=c.world[2])

	def _base_scint_solids(self, solids: ET.Element) -> None:
		c = self.config
		_sub(solids, "tube", name="ThinScinA", aunit="deg", deltaphi=360, lunit="mm", rmax=c.base_scint_radius, rmin=0, startphi=0, z=c.base_scint_a_thickness)
		_sub(solids, "tube", name="ThinScinA_Al", aunit="deg", deltaphi=360, lunit="mm", rmax=c.base_scint_radius + 0.0001, rmin=0, startphi=0, z=c.base_scint_a_al_thickness)
		inner = _sub(solids, "subtraction", name="InnerAl")
		_sub(inner, "first", ref="ThinScinA_Al")
		_sub(inner, "second", ref="ThinScinA")
		_sub(solids, "tube", name="ThicScinB", aunit="deg", deltaphi=360, lunit="mm", rmax=c.base_scint_radius, rmin=0, startphi=0, z=c.base_scint_b_thickness)
		_sub(solids, "tube", name="ThicScinB_Al", aunit="deg", deltaphi=360, lunit="mm", rmax=c.base_scint_radius + 0.0001, rmin=0, startphi=0, z=c.base_scint_b_al_thickness)
		inner2 = _sub(solids, "subtraction", name="InnerAl20")
		_sub(inner2, "first", ref="ThicScinB_Al")
		_sub(inner2, "second", ref="ThicScinB")

	def _telescope_solids(self, solids: ET.Element) -> None:
		c = self.config
		_sub(solids, "polycone", name="WedgeRevolvedRing", aunit="deg", deltaphi=360, lunit="mm", startphi=0)
		# The telescope ring is only used for GPS confinement and plotting;
		# its exact z-planes are not critical for detector ingestion.
		_sub(solids, "box", name="SmallScin", lunit="mm", x=c.telescope_scint_xy, y=c.telescope_scint_xy, z=c.telescope_scint_z)
		_sub(solids, "box", name="SmallScinMA", lunit="mm", x=c.telescope_scint_ma_xy, y=c.telescope_scint_ma_xy, z=c.telescope_scint_ma_z)
		inner = _sub(solids, "subtraction", name="InnerAl")
		_sub(inner, "first", ref="SmallScinMA")
		_sub(inner, "second", ref="SmallScin")
		_sub(solids, "box", name="SmallScinM", lunit="mm", x=c.telescope_scint_m_xy, y=c.telescope_scint_m_xy, z=c.telescope_scint_m_z)
		mylar = _sub(solids, "subtraction", name="Mylar")
		_sub(mylar, "first", ref="SmallScinM")
		_sub(mylar, "second", ref="SmallScinMA")
		_sub(solids, "box", name="SmallScinMAOuter", lunit="mm", x=c.telescope_scint_outer_xy, y=c.telescope_scint_outer_xy, z=c.telescope_scint_outer_z)
		outer = _sub(solids, "subtraction", name="OuterAl")
		_sub(outer, "first", ref="SmallScinMAOuter")
		_sub(outer, "second", ref="SmallScinM")
		_sub(solids, "subtraction", name="InnerAl20")
		_sub(solids, "subtraction", name="Mylar20")
		_sub(solids, "subtraction", name="OuterAl20")
		_sub(solids, "tube", name="sipmLidBase", aunit="deg", deltaphi=360, lunit="mm", rmax=60.325, rmin=5.08, startphi=0, z=8.89)
		_sub(solids, "tube", name="sipmLidSource", aunit="deg", deltaphi=360, lunit="mm", rmax=6.325001, rmin=0, startphi=0, z=8.89)
		lid = _sub(solids, "subtraction", name="sipmLid")
		_sub(lid, "first", ref="sipmLidBase")
		_sub(lid, "second", ref="sipmLidSource")
		_position(lid, "sipmLid_pos", 0, 0, -1.2064)
		_sub(solids, "box", name="scintHolderbase", lunit="mm", x=22.2, y=22.2, z=30.4434)
		_sub(solids, "box", name="scintHoldercutout20", lunit="mm", x=16.3, y=16.3, z=50.4434)
		holder10 = _sub(solids, "subtraction", name="scintHolder10")
		_sub(holder10, "first", ref="scintHolderbase")
		_sub(holder10, "second", ref="scintHoldercutout20")
		_sub(solids, "box", name="scintHoldercutout10", lunit="mm", x=20.5, y=20.5, z=29.4434)
		holder = _sub(solids, "subtraction", name="scintHolder")
		_sub(holder, "first", ref="scintHolder10")
		_sub(holder, "second", ref="scintHoldercutout10")
		_sub(solids, "tube", name="scintHolderbaseRing", aunit="deg", deltaphi=360, lunit="mm", rmax=22.86, rmin=10, startphi=0, z=1.575)
		_sub(solids, "box", name="scintSleevebase", lunit="mm", x=20.5, y=20.5, z=10)
		_sub(solids, "box", name="scintHoldercutout20b", lunit="mm", x=13.5, y=13.5, z=20)
		sleeve = _sub(solids, "subtraction", name="scintSleeve")
		_sub(sleeve, "first", ref="scintSleevebase")
		_sub(sleeve, "second", ref="scintHoldercutout20b")
		_sub(solids, "tube", name="SourceAlo", aunit="deg", deltaphi=360, lunit="mm", rmax=6.325, rmin=4.694999, startphi=0, z=3.1749)
		_sub(solids, "tube", name="SourceAl1", aunit="deg", deltaphi=360, lunit="mm", rmax=6.325, rmin=0, startphi=0, z=12.7)
		_sub(solids, "tube", name="SourceAl2", aunit="deg", deltaphi=360, lunit="mm", rmax=3.1623, rmin=0, startphi=0, z=6.3498)
		source_al = _sub(solids, "subtraction", name="SourceAl")
		_sub(source_al, "first", ref="SourceAl1")
		_sub(source_al, "second", ref="SourceAl2")
		_position(source_al, "SourceAl_pos", 0, 0, -4.76255)
		_sub(solids, "tube", name="SourceMy", aunit="deg", deltaphi=360, lunit="mm", rmax=4.694999, rmin=0, startphi=0, z=0.0064)
		_sub(solids, "tube", name="SourceCal", aunit="deg", deltaphi=360, lunit="mm", rmax=4.694999, rmin=0, startphi=0, z=0.001)

	def _source_solids(self, solids: ET.Element) -> None:
		# All source solids are created in the geometry-specific solid builder.
		return

	def _holder_solids(self, solids: ET.Element) -> None:
		# All holder solids are created in the geometry-specific solid builder.
		return

	def _structure(self, root: ET.Element) -> None:
		structure = ET.SubElement(root, "structure")
		c = self.config

		def volume(name: str, material: str, solid: str) -> ET.Element:
			v = _sub(structure, "volume", name=name)
			_sub(v, "materialref", ref=material)
			_sub(v, "solidref", ref=solid)
			return v

		vacuum = volume("Vacuum", "interGalactic", "InnerTPipe")
		volume("SteelTPipe", "Stainless_Steel", "TDecayVolume2")
		volume("tPipeFlangeScint", "Stainless_Steel", "FlangeScint")
		volume("tPipeFlange1", "Stainless_Steel", "Flange")
		volume("tPipeFlange2", "Stainless_Steel", "Flange")
		volume("tPipeFlange3", "Stainless_Steel", "Flange")
		volume("VacuumWindow", "G4_Al", "VacuumWindowDisk")

		if c.use_telescope_parts:
			volume("WedgeRevolvedRingLV", "Stainless_Steel", "WedgeRevolvedRing")
			volume("AScintillator", "G4_PLASTIC_SC_VINYLTOLUENE", "SmallScin")
			volume("Ali_sq1", "G4_Al", "InnerAl")
			volume("Mylar_sq1", "G4_MYLAR", "Mylar")
			volume("Alo_sq1", "G4_Al", "OuterAl")
			volume("BScintillatorLV", "G4_PLASTIC_SC_VINYLTOLUENE", "SmallScin")
			volume("Ali_sq2", "G4_Al", "InnerAl20")
			volume("Mylar_sq2", "G4_MYLAR", "Mylar20")
			volume("Alo_sq2", "G4_Al", "OuterAl20")
			volume("sipm_Lid", "G4_Cu", "sipmLid")
			volume("scint_holder", "Vetronite", "scintHolder")
			volume("scint_holderBase", "Vetronite", "scintHolderbaseRing")
			volume("scint_Sleeve", "G4_Al", "scintSleeve")
			volume("SourceAlo", "G4_Al", "SourceAlo")
			volume("SourceAl", "G4_Al", "SourceAl")
			volume("SourceMy", "G4_MYLAR", "SourceMy")
			volume("SourceCal", "G4_Al", "SourceCal")
		else:
			volume("AScintillator", "G4_PLASTIC_SC_VINYLTOLUENE", "ThinScinA")
			volume("Ali_sq1", "G4_Al", "InnerAl")
			volume("Mylar_sq1", "G4_MYLAR", "Mylar")
			volume("Alo_sq1", "G4_Al", "OuterAl")
			volume("BScintillatorLV", "G4_PLASTIC_SC_VINYLTOLUENE", "ThicScinB")
			volume("Ali_sq2", "G4_Al", "InnerAl20")
			volume("Mylar_sq2", "G4_MYLAR", "Mylar20")
			volume("Alo_sq2", "G4_Al", "OuterAl20")
			volume("SourceAlo", "G4_Al", "SourceAlo")
			volume("SourceAl", "G4_Al", "SourceAl")
			volume("SourceMy", "G4_MYLAR", "SourceMy")
			volume("SourceCal", "G4_Al", "SourceCal")

		self._place_world(vacuum)
		world = _sub(structure, "volume", name="World")
		_sub(world, "materialref", ref="G4_AIR")
		_sub(world, "solidref", ref="World")

		if c.use_telescope_parts:
			self._place_telescope_physvols(world)
		else:
			self._place_base_physvols(world)

	def _place_world(self, vacuum: ET.Element) -> None:
		return

	def _place_base_physvols(self, world: ET.Element) -> None:
		self._place_common_physvols(world, base=True)

	def _place_telescope_physvols(self, world: ET.Element) -> None:
		self._place_common_physvols(world, base=False)
		_sub(world, "physvol", name="WedgeRevolvedRingPV")

	def _place_common_physvols(self, world: ET.Element, base: bool) -> None:
		if base:
			pipe_z = 111.125
			window_z = 0.025
			a_z = -9.0001
			b_z = -38.0002
		else:
			pipe_z = 98.425
			window_z = -0.0381
			a_z = -9.9836
			b_z = -13.136

		_phys(world, "Vacuum", "Vacuum", (0, 0, pipe_z), (0, 90, 0), copynumber=1)
		_phys(world, "SteelTPipe", "SteelTPipe", (0, 0, pipe_z), (0, 90, 0))
		_phys(world, "tPipeFlangeScint", "tPipeFlangeScint", (0, 0, 6.35 if base else -6.35), None)
		_phys(world, "tPipeFlange1", "tPipeFlange1", (0, 0, 19.05 if base else 6.35), None)
		_phys(world, "tPipeFlange2", "tPipeFlange2", (92.075, 0, 111.125 if base else 98.425), (0, 90, 0))
		_phys(world, "tPipeFlange3", "tPipeFlange3", (-92.075, 0, 111.125 if base else 98.425), (0, 90, 0))
		_phys(world, "VacuumWindow", "VacuumWindow", (0, 0, window_z), None, copynumber=2)
		if not base:
			_phys(world, "WedgeRevolvedRingLV", "WedgeRevolvedRingLV", (0, 0, 0), None)
		_phys(world, "AScintillator", "AScintillator", (0, 0, a_z), None, copynumber=3)
		_phys(world, "Ali_sq1", "Ali_sq1", (0, 0, a_z), None)
		if base:
			_phys(world, "Mylar_sq1", "Mylar_sq1", (0, 0, a_z), None)
			_phys(world, "Alo_sq1", "Alo_sq1", (0, 0, a_z), None)
		_phys(world, "BScintillatorLV", "BScintillatorLV", (0, 0, b_z), None, copynumber=4)
		_phys(world, "Ali_sq2", "Ali_sq2", (0, 0, b_z), None)
		if base:
			_phys(world, "Mylar_sq2", "Mylar_sq2", (0, 0, b_z), None)
			_phys(world, "Alo_sq2", "Alo_sq2", (0, 0, b_z), None)
		if not base:
			_phys(world, "Mylar_sq1", "Mylar_sq1", (0, 0, a_z), None)
			_phys(world, "Alo_sq1", "Alo_sq1", (0, 0, a_z), None)
			_phys(world, "Mylar_sq2", "Mylar_sq2", (0, 0, b_z), None)
			_phys(world, "Alo_sq2", "Alo_sq2", (0, 0, b_z), None)
			_phys(world, "sipm_Lid", "sipm_Lid", (0, 0, -44.3708), None)
			_phys(world, "scint_holder", "scint_holder", (0, 0, -23.1291), None)
			_phys(world, "scint_holderBase", "scint_holderBase", (0, 0, -39.1383), None)
			_phys(world, "scint_Sleeve", "scint_Sleeve", (0, 0, -19.7122), None)
			_phys(world, "SourceAlo", "SourceAlo", (0, 0, -42.719651), None)
			_phys(world, "SourceAl", "SourceAl", (0, 0, -50.657102), None)
			_phys(world, "SourceMy", "SourceMy", (0, 0, -42.71965), None)
			_phys(world, "SourceCal", "SourceCal", (0, 0, -44.3066), None, copynumber=1)
		else:
			_phys(world, "SourceAlo", "SourceAlo", (0, 0, -42.719651), None)
			_phys(world, "SourceAl", "SourceAl", (0, 0, -50.657102), None)
			_phys(world, "SourceMy", "SourceMy", (0, 0, -42.71965), None)
			_phys(world, "SourceCal", "SourceCal", (0, 0, -44.3066), None, copynumber=1)


def _phys(parent: ET.Element, name: str, volumeref: str, position: tuple[float, float, float] | None, rotation: tuple[float, float, float] | None, copynumber: int | None = None) -> ET.Element:
	attrs = {"name": name}
	if copynumber is not None:
		attrs["copynumber"] = str(copynumber)
	physvol = _sub(parent, "physvol", **attrs)
	_sub(physvol, "volumeref", ref=volumeref)
	if position is not None:
		_position(physvol, f"{name}_pos", *position)
	if rotation is not None:
		_rotation(physvol, f"{name}_rot", *rotation)
	return physvol


def build_geometry(output_path: str | Path, telescope: bool = False) -> Path:
	builder = BetaMonitorGDMLBuilder.telescope() if telescope else BetaMonitorGDMLBuilder.base()
	return builder.write(output_path)


def main(argv: list[str] | None = None) -> int:
	parser = argparse.ArgumentParser(description="Generate BetaMonitor GDML geometry.")
	parser.add_argument("output", nargs="?", help="Output GDML file path.")
	parser.add_argument("--telescope", action="store_true", help="Generate the telescope geometry variant.")
	args = parser.parse_args(argv)

	default_output = "../dat/geometry_telescope_export.gdml" if args.telescope else "../dat/geometry_export.gdml"
	output_path = Path(args.output or default_output)
	written = build_geometry(output_path, telescope=args.telescope)
	print(written)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
