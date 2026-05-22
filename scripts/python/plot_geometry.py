import pyg4ometry as pg4
import vtk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import copy
from matplotlib.patches import Polygon as MplPolygon


class GeometryPlotter:

    VIEW_CFG = [
        ("XY  (top)",   "X (mm)", "Y (mm)", 0, 1, 2),
        ("XZ  (front)", "X (mm)", "Z (mm)", 0, 2, 1),
        ("YZ  (side)",  "Y (mm)", "Z (mm)", 1, 2, 0),
    ]

    COLOR_CYCLE = list(plt.cm.tab10.colors)

    VTK_COLORS = [
        [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
        [0.5, 0.5, 0.5], [1.0, 0.5, 0.0], [0.5, 0.0, 1.0],
        [0.0, 0.5, 1.0], [0.5, 1.0, 0.0], [0.0, 1.0, 0.5],
        [1.0, 0.0, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
    ]

    def __init__(self, gdml_path):
        gdml_parser = pg4.gdml.Reader(gdml_path)
        self.reg = gdml_parser.getRegistry()
        self._material_colors = {}

    # ------------------------------------------------------------------ helpers

    def _eval_param(self, p):
        try:
            if hasattr(p, 'eval'):
                return float(p.eval())
            return float(p)
        except Exception:
            return None

    def _get_solid_dims(self, solid, prefix=''):
        """Return a list of human-readable dimension strings for a solid.
        For Boolean solids, recurses into both constituents."""
        lines = []
        t = type(solid).__name__
        try:
            if any(k in t for k in ('Subtraction', 'Union', 'Intersection')):
                op = 'sub' if 'Subtraction' in t else ('uni' if 'Union' in t else 'int')
                obj1 = getattr(solid, 'obj1', None)
                obj2 = getattr(solid, 'obj2', None)
                lines.append(f"[{op}]")
                if obj1 is not None:
                    sub = self._get_solid_dims(obj1, prefix='A:')
                    lines.extend([f"  A:{s}" if not s.startswith('  ') else s for s in sub])
                if obj2 is not None:
                    sub = self._get_solid_dims(obj2, prefix='B:')
                    lines.extend([f"  B:{s}" if not s.startswith('  ') else s for s in sub])
            elif 'Tubs' in t or 'Tube' in t:
                rmin = self._eval_param(getattr(solid, 'pRMin', None))
                rmax = self._eval_param(getattr(solid, 'pRMax', None))
                dz   = self._eval_param(getattr(solid, 'pDz',   None))
                if rmin is not None and rmin > 1e-6:
                    lines.append(f"{prefix}r_in={rmin:.3g} mm")
                if rmax is not None:
                    lines.append(f"{prefix}r_out={rmax:.3g} mm")
                if dz is not None:
                    lines.append(f"{prefix}h={2*dz:.3g} mm")
            elif 'Box' in t:
                px = self._eval_param(getattr(solid, 'pX', None))
                py = self._eval_param(getattr(solid, 'pY', None))
                pz = self._eval_param(getattr(solid, 'pZ', None))
                if px is not None: lines.append(f"{prefix}dx={2*px:.3g} mm")
                if py is not None: lines.append(f"{prefix}dy={2*py:.3g} mm")
                if pz is not None: lines.append(f"{prefix}dz={2*pz:.3g} mm")
            elif 'Sphere' in t:
                rmin = self._eval_param(getattr(solid, 'pRmin', None))
                rmax = self._eval_param(getattr(solid, 'pRmax', None))
                if rmin is not None and rmin > 1e-6:
                    lines.append(f"{prefix}r_in={rmin:.3g} mm")
                if rmax is not None:
                    lines.append(f"{prefix}r_out={rmax:.3g} mm")
            elif 'Cons' in t:
                rm1 = self._eval_param(getattr(solid, 'pRmax1', None))
                rm2 = self._eval_param(getattr(solid, 'pRmax2', None))
                dz  = self._eval_param(getattr(solid, 'pDz',    None))
                if rm1 is not None: lines.append(f"{prefix}r1={rm1:.3g} mm")
                if rm2 is not None: lines.append(f"{prefix}r2={rm2:.3g} mm")
                if dz  is not None: lines.append(f"{prefix}h={2*dz:.3g} mm")
            elif 'Trd' in t:
                px1 = self._eval_param(getattr(solid, 'pX1', None))
                px2 = self._eval_param(getattr(solid, 'pX2', None))
                py1 = self._eval_param(getattr(solid, 'pY1', None))
                py2 = self._eval_param(getattr(solid, 'pY2', None))
                pz  = self._eval_param(getattr(solid, 'pZ',  None))
                if px1 is not None: lines.append(f"{prefix}dx1={2*px1:.3g} mm")
                if px2 is not None: lines.append(f"{prefix}dx2={2*px2:.3g} mm")
                if py1 is not None: lines.append(f"{prefix}dy1={2*py1:.3g} mm")
                if py2 is not None: lines.append(f"{prefix}dy2={2*py2:.3g} mm")
                if pz  is not None: lines.append(f"{prefix}dz={2*pz:.3g} mm")
            else:
                lines.append(f"{prefix}[{t}]")
        except Exception:
            lines.append(f"{prefix}[{t}]")
        return lines

    def _get_material_color(self, mat_name):
        if mat_name not in self._material_colors:
            idx = len(self._material_colors) % len(self.COLOR_CYCLE)
            self._material_colors[mat_name] = self.COLOR_CYCLE[idx]
        return self._material_colors[mat_name]

    def _get_position(self, pv, verbose=True):
        try:
            p = pv.position
            if p is None:
                return np.zeros(3)
            if hasattr(p, 'eval'):
                return np.array([float(v) for v in p.eval()])
            if hasattr(p, 'x'):
                return np.array([float(p.x.eval()), float(p.y.eval()), float(p.z.eval())])
        except Exception as e:
            if verbose:
                print(f"    Position eval error: {e}")
        return np.zeros(3)

    def _get_rotation(self, pv, verbose=True):
        try:
            r = pv.rotation
            if r is None:
                return np.eye(3)
            angles = r.eval() if hasattr(r, 'eval') else [float(r.x.eval()), float(r.y.eval()), float(r.z.eval())]
            rx, ry, rz = [float(a) for a in angles]
            ry = -ry
            cx, sx = np.cos(rx), np.sin(rx)
            cy, sy = np.cos(ry), np.sin(ry)
            cz, sz = np.cos(rz), np.sin(rz)
            Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
            Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
            Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
            return Rz @ Ry @ Rx
        except Exception as e:
            if verbose:
                print(f"    Rotation eval error: {e}")
        return np.eye(3)

    def _collect_faces(self, pv, parent_trans, parent_rot, all_faces, verbose=True, vol_infos=None):
        local_pos = self._get_position(pv, verbose)
        local_rot = self._get_rotation(pv, verbose)
        world_rot   = parent_rot @ local_rot
        world_trans = parent_rot @ local_pos + parent_trans

        lv = pv.logicalVolume
        vol_name = lv.name.split('0x')[0]
        mat_name = lv.material.name.split('0x')[0]
        color = self._get_material_color(mat_name)

        if verbose:
            print(f"  {vol_name}  pos={np.round(world_trans, 2)}  mat={mat_name}")

        if vol_infos is not None:
            dim_lines = self._get_solid_dims(lv.solid)
            vol_infos.append((vol_name, world_trans.copy(), dim_lines, color))

        try:
            m = lv.solid.mesh()
            verts, polys, *_ = m.toVerticesAndPolygons()
            verts = np.array([[float(v[0]), float(v[1]), float(v[2])] for v in verts], dtype=float)
            if len(verts) >= 3 and len(polys) > 0:
                verts = (world_rot @ verts.T).T + world_trans
                for vi, (_, _, _, i, j, d) in enumerate(self.VIEW_CFG):
                    for poly in polys:
                        face = verts[list(poly)]
                        depth = float(np.mean(face[:, d]))
                        all_faces[vi].append((depth, face[:, [i, j]], color))
        except Exception as e:
            if verbose:
                print(f"    Mesh failed: {e}")

        for daughter in lv.daughterVolumes:
            self._collect_faces(daughter, world_trans, world_rot, all_faces, verbose, vol_infos)

    def _collect_cross_section_segments(self, pv, parent_trans, parent_rot, segments, plane_axis, plane_val, vol_infos=None, verbose=True):
        """Walk the volume tree and collect 3-D line segments that are the
        intersection of each mesh face with the plane `axis[plane_axis] == plane_val`."""
        local_pos = self._get_position(pv, verbose)
        local_rot = self._get_rotation(pv, verbose)
        world_rot   = parent_rot @ local_rot
        world_trans = parent_rot @ local_pos + parent_trans

        lv = pv.logicalVolume
        vol_name = lv.name.split('0x')[0]
        mat_name = lv.material.name.split('0x')[0]
        color = self._get_material_color(mat_name)

        if verbose:
            print(f"  {vol_name}  pos={np.round(world_trans, 2)}  mat={mat_name}")

        if vol_infos is not None:
            dim_lines = self._get_solid_dims(lv.solid)
            vol_infos.append((vol_name, world_trans.copy(), dim_lines, color))

        try:
            m = lv.solid.mesh()
            verts, polys, *_ = m.toVerticesAndPolygons()
            verts = np.array([[float(v[0]), float(v[1]), float(v[2])] for v in verts], dtype=float)
            if len(verts) >= 3 and len(polys) > 0:
                verts = (world_rot @ verts.T).T + world_trans
                for poly in polys:
                    face_v = verts[list(poly)]
                    n = len(face_v)
                    vals = face_v[:, plane_axis] - plane_val
                    cross_pts = []
                    for k in range(n):
                        a = float(vals[k])
                        b = float(vals[(k + 1) % n])
                        if a == 0.0 and b != 0.0:
                            cross_pts.append(face_v[k])
                        elif (a > 0.0) != (b > 0.0):   # strictly crossing edge
                            t = a / (a - b)
                            cross_pts.append(face_v[k] + t * (face_v[(k + 1) % n] - face_v[k]))
                    if len(cross_pts) >= 2:
                        segments.append((cross_pts[0], cross_pts[-1], color))
        except Exception as e:
            if verbose:
                print(f"    Mesh failed: {e}")

        for daughter in lv.daughterVolumes:
            self._collect_cross_section_segments(daughter, world_trans, world_rot, segments,
                                                  plane_axis, plane_val, vol_infos, verbose)

    # ------------------------------------------------------------------ VTK helpers

    def add_vtk_text(self, viewer, text, pos=(0.05, 0.9), font_size=24, color=(1, 1, 1)):
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(text)
        text_prop = text_actor.GetTextProperty()
        text_prop.SetFontSize(font_size)
        text_prop.SetColor(*color)
        text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        text_actor.SetPosition(pos[0], pos[1])
        viewer.ren.AddActor2D(text_actor)

    def add_vtk_point(self, viewer, position, color=(1, 0, 0), radius=2.0):
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(*position)
        sphere.SetRadius(radius)
        sphere.SetThetaResolution(16)
        sphere.SetPhiResolution(16)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color)
        viewer.ren.AddActor(actor)

    def add_vtk_axes(self, viewer, origin=(0, 0, 0), length=50, shaft_radius=0.01, tip_length=0.1, tip_radius=0.05):
        axes_cfg = [
            ((1, 0, 0), (1, 0, 0)),
            ((0, 1, 0), (0, 1, 0)),
            ((0, 0, 1), (0, 0, 1)),
        ]
        for direction, color in axes_cfg:
            arrow = vtk.vtkArrowSource()
            arrow.SetShaftRadius(shaft_radius)
            arrow.SetTipLength(tip_length)
            arrow.SetTipRadius(tip_radius)
            transform = vtk.vtkTransform()
            if direction == (0, 1, 0):
                transform.RotateZ(90)
            elif direction == (0, 0, 1):
                transform.RotateY(-90)
            transform.Translate(origin)
            transform.Scale(length, length, length)
            transform_filter = vtk.vtkTransformPolyDataFilter()
            transform_filter.SetTransform(transform)
            transform_filter.SetInputConnection(arrow.GetOutputPort())
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transform_filter.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(*color)
            viewer.ren.AddActor(actor)

    # ------------------------------------------------------------------ 3D VTK plots

    def plot_by_material(self):
        v = pg4.visualisation.VtkViewer()
        n = 0
        legend_entries = []
        for _, mat in self.reg.materialDict.items():
            if mat.type == 'composite':
                mat_name = mat.name.split('0x')[0]
                mat_vis = pg4.visualisation.VisualisationOptions()
                mat_vis.colour = self.VTK_COLORS[n % len(self.VTK_COLORS)]
                mat_vis.alpha = 0.2
                v.addMaterialVisOption(mat_name, mat_vis)
                print('Added material vis for', mat_name)
                legend_entries.append((mat_name, self.VTK_COLORS[n % len(self.VTK_COLORS)]))
                n += 1
            else:
                print('Skipped material', mat.name, 'of type', mat.type)
        v.addLogicalVolume(self.reg.getWorldVolume())
        fig, ax = plt.subplots()
        patches = [mpatches.Patch(color=color, label=mat_name) for mat_name, color in legend_entries]
        ax.legend(handles=patches, loc='center')
        ax.axis('off')
        plt.show(block=False)
        return v

    def plot_by_volume(self, vol_select=None):
        v = pg4.visualisation.VtkViewer()
        volumes = self.reg.logicalVolumeDict
        n = 0
        legend_entries = []
        original_mat = None
        for _, vol in volumes.items():
            vol_name = vol.name.split('0x')[0]
            mat_name = vol.material.name.split('0x')[0]
            mat_vis = pg4.visualisation.VisualisationOptions()
            if vol_select:
                if vol_name == vol_select:
                    print('Selected volume', vol_name)
                    original_mat = vol.material
                    unique_mat = copy.copy(vol.material)
                    unique_mat.name = "UniqueMat"
                    vol.material = unique_mat
                    legend_entries.append((vol_name, [1.0, 0.0, 0.0]))
                    mat_vis.colour = [1.0, 0.0, 0.0]
                    mat_vis.alpha = 0.3
                    v.addMaterialVisOption("UniqueMat", mat_vis)
                else:
                    mat_vis.colour = [0.8, 0.8, 0.8]
                    mat_vis.alpha = 0.1
                    v.addMaterialVisOption(mat_name, mat_vis)
            else:
                mat_vis.colour = self.VTK_COLORS[n % len(self.VTK_COLORS)]
                mat_vis.alpha = 0.0
                v.addMaterialVisOption(mat_name, mat_vis)
                print('Added volume vis for', vol_name)
                legend_entries.append((vol_name, self.VTK_COLORS[n % len(self.VTK_COLORS)]))
                n += 1
        self.add_vtk_text(v, f'{vol_select}', pos=(0.05, 0.9), font_size=48, color=(1, 0, 0))
        v.addLogicalVolume(self.reg.getWorldVolume())
        if vol_select and original_mat is not None:
            for _, vol in volumes.items():
                if vol.name.split('0x')[0] == vol_select:
                    vol.material = original_mat
        return v

    # ------------------------------------------------------------------ screenshot utilities

    def save_three_views(self, viewer, base_filename="view"):
        viewer.ren.GetActiveCamera().SetPosition(0, 0, 1000)
        viewer.ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
        viewer.renWin.Render()
        viewer.exportScreenShot(f"{base_filename}_xy.png")

        viewer.ren.GetActiveCamera().SetPosition(0, 1000, 0)
        viewer.ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
        viewer.renWin.Render()
        viewer.exportScreenShot(f"{base_filename}_xz.png")

        viewer.ren.GetActiveCamera().SetPosition(1000, 0, 0)
        viewer.ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
        viewer.renWin.Render()
        viewer.exportScreenShot(f"{base_filename}_yz.png")

        self.combine_views_to_png(
            f"{base_filename}_xy.png",
            f"{base_filename}_xz.png",
            f"{base_filename}_yz.png",
            out_path=f"{base_filename}_three_views.png",
        )

    def combine_views_to_png(self, xy_path, xz_path, yz_path, out_path="three_views.png"):
        xy_img = mpimg.imread(xy_path)
        xz_img = mpimg.imread(xz_path)
        yz_img = mpimg.imread(yz_path)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        axes[0].imshow(xy_img)
        axes[0].set_title("XY View")
        axes[0].axis('off')
        axes[1].imshow(xz_img)
        axes[1].set_title("XZ View")
        axes[1].axis('off')
        axes[2].imshow(yz_img)
        axes[2].set_title("YZ View")
        axes[2].axis('off')
        plt.tight_layout()
        plt.show(block=False)

    # ------------------------------------------------------------------ utilities

    def print_vol_names(self):
        vol_list = []
        for _, vol in self.reg.logicalVolumeDict.items():
            vol_name = vol.name.split('0x')[0]
            mat_name = vol.material.name.split('0x')[0]
            print(f'{vol_name}; (Material: {mat_name})')
            vol_list.append(vol_name)
        return vol_list

    # ------------------------------------------------------------------ 2D matplotlib plot

    def plot_three_view_2d(self, alpha=0.25, linewidth=0.8, figsize=None, verbose=True, annotate=True):
        self._material_colors = {}

        # XZ view is index 1 in VIEW_CFG: ("XZ  (front)", "X (mm)", "Z (mm)", 0, 2, 1)
        XZ_VI = 1
        xz_title, xz_xl, xz_yl, xz_ih, xz_iv, _ = self.VIEW_CFG[XZ_VI]

        if annotate:
            if figsize is None:
                figsize = (20, 9)
            fig = plt.figure(figsize=figsize)
            gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.6], figure=fig, wspace=0.08)
            ax_xz = fig.add_subplot(gs[0, 0])
            ax_info = fig.add_subplot(gs[0, 1])
            ax_info.axis('off')
            ax_info.set_title('Volume Callouts', fontsize=9, pad=4)
            ax_xz.set_title(xz_title)
            ax_xz.set_xlabel(xz_xl)
            ax_xz.set_ylabel(xz_yl)
            ax_xz.set_aspect('equal')
            ax_xz.set_facecolor('white')
            ax_xz.grid(True, linewidth=0.3, alpha=0.4, zorder=0)
            draw_axes = [ax_xz]
        else:
            if figsize is None:
                figsize = (18, 7)
            fig, raw_axes = plt.subplots(1, 3, figsize=figsize)
            draw_axes = list(raw_axes)
            ax_info = None
            for ax, (title, xl, yl, *_) in zip(draw_axes, self.VIEW_CFG):
                ax.set_title(title)
                ax.set_xlabel(xl)
                ax.set_ylabel(yl)
                ax.set_aspect('equal')
                ax.set_facecolor('white')
                ax.grid(True, linewidth=0.3, alpha=0.4, zorder=0)

        world_lv = self.reg.getWorldVolume()
        daughters = world_lv.daughterVolumes
        print(f"World volume has {len(daughters)} direct daughters")

        if annotate:
            # Single Y=0 cross-section for the XZ view
            vol_infos = []
            segments = []
            for pv in daughters:
                self._collect_cross_section_segments(
                    pv, np.zeros(3), np.eye(3), segments,
                    plane_axis=1, plane_val=0.0,
                    vol_infos=vol_infos, verbose=verbose)
            for pt_a, pt_b, color in segments:
                ax_xz.plot([pt_a[xz_ih], pt_b[xz_ih]], [pt_a[xz_iv], pt_b[xz_iv]],
                           color=color, linewidth=linewidth, alpha=0.9, solid_capstyle='round')
            ax_xz.autoscale_view()
        else:
            vol_infos = None
            # Each view cuts at its own depth axis == 0
            # VIEW_CFG: (title, xl, yl, i_h, i_v, i_d)
            for vi, (_, _, _, i_h, i_v, i_d) in enumerate(self.VIEW_CFG):
                segments = []
                for pv in daughters:
                    self._collect_cross_section_segments(
                        pv, np.zeros(3), np.eye(3), segments,
                        plane_axis=i_d, plane_val=0.0,
                        vol_infos=None, verbose=(verbose and vi == 0))
                for pt_a, pt_b, color in segments:
                    draw_axes[vi].plot([pt_a[i_h], pt_b[i_h]], [pt_a[i_v], pt_b[i_v]],
                                       color=color, linewidth=linewidth, alpha=0.9,
                                       solid_capstyle='round')
                draw_axes[vi].autoscale_view()

        if annotate and vol_infos:
            n_vols = len(vol_infos)
            LABEL_X = -80.0  # mm — fixed X column for all markers

            # Force layout so transforms are initialised
            fig.canvas.draw()

            orig_data = [(float(info[1][xz_ih]), float(info[1][xz_iv])) for info in vol_infos]

            # Spread labels evenly in Y (Z in XZ view) sorted by their natural Z
            sorted_idx = sorted(range(n_vols), key=lambda i: orig_data[i][1])
            t_to_disp = ax_xz.transData
            r_px     = 5.5 * (fig.dpi / 72.0)
            min_sep_px = r_px * 2 * 1.2

            # Start with Y equal to each volume's own Z, then push apart
            label_y_disp = [t_to_disp.transform((LABEL_X, orig_data[i][1]))[1]
                            for i in sorted_idx]
            for _ in range(500):
                moved = False
                for k in range(len(label_y_disp) - 1):
                    gap = label_y_disp[k + 1] - label_y_disp[k]
                    if gap < min_sep_px:
                        push = (min_sep_px - gap) * 0.5
                        label_y_disp[k]     -= push
                        label_y_disp[k + 1] += push
                        moved = True
                if not moved:
                    break

            t_inv = ax_xz.transData.inverted()
            label_pos = {}   # original index -> (LABEL_X, y_data)
            for k, orig_i in enumerate(sorted_idx):
                _, y_data = t_inv.transform((0, label_y_disp[k]))
                label_pos[orig_i] = (LABEL_X, y_data)

            # Place numbered markers on XZ view only
            for idx, (vol_name, world_pos, dim_lines, color) in enumerate(vol_infos):
                num_label  = str(idx + 1)
                xc_orig, yc_orig = orig_data[idx]
                xc_lbl, yc_lbl  = label_pos[idx]
                # Leader line from true centroid to label column
                ax_xz.plot([xc_orig, xc_lbl], [yc_orig, yc_lbl], '-',
                           color=color, linewidth=0.6, alpha=0.5, zorder=9)
                ax_xz.plot(xc_lbl, yc_lbl, 'o', color=color, markersize=11, alpha=0.75, zorder=10)
                ax_xz.text(xc_lbl, yc_lbl, num_label, color='white', fontsize=5.5,
                           ha='center', va='center', fontweight='bold', zorder=11)

            # Multi-column callout panel: row height scales with content
            # 2 base lines (name + pos) + 1 per dim line, at ~0.018 per line
            line_h = 0.018
            row_heights = [max(0.055, (2 + len(info[2])) * line_h + 0.01) for info in vol_infos]
            y_start = 0.97
            # Pack rows greedily into columns of height <= 0.95
            cols_assignment = []  # (col, y_top) per entry
            col = 0
            y_cursor = y_start
            for rh in row_heights:
                if y_cursor - rh < 0.02 and cols_assignment:
                    col += 1
                    y_cursor = y_start
                cols_assignment.append((col, y_cursor))
                y_cursor -= rh
            n_cols_info = col + 1
            col_width = 1.0 / n_cols_info

            for idx, (vol_name, world_pos, dim_lines, color) in enumerate(vol_infos):
                num_label = str(idx + 1)
                col, y_top = cols_assignment[idx]
                x_body = col * col_width + 0.01
                pos_str = f"({world_pos[0]:.2f}, {world_pos[1]:.2f}, {world_pos[2]:.2f}) mm"
                body = f"{num_label}. {vol_name}\n  pos: {pos_str}"
                if dim_lines:
                    body += "\n  " + "\n  ".join(dim_lines)
                ax_info.text(
                    x_body, y_top, body,
                    transform=ax_info.transAxes,
                    fontsize=6.0, va='top', color='black', family='monospace',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              edgecolor=color, alpha=0.85, linewidth=0.8),
                )

        patches_legend = [mpatches.Patch(color=c, label=m) for m, c in self._material_colors.items()]
        fig.legend(handles=patches_legend, loc='lower center', ncol=min(len(patches_legend), 8), fontsize=8)
        plt.tight_layout(rect=[0, 0.06, 1, 1])
        plt.show()
        return fig, draw_axes


if __name__ == "__main__":
    plotter = GeometryPlotter("output/geometry_telescope_export.gdml")
    # plotter = GeometryPlotter("/Users/harperumfress/UW/betamonitor_data/original_singlethread_data/geometry_export.gdml")
    plotter.plot_three_view_2d()
    plotter.print_vol_names()

    v = plotter.plot_by_volume(vol_select='SourceCal')
    v.ren.GetActiveCamera().SetPosition(-100, 0, 0)
    v.ren.GetActiveCamera().SetFocalPoint(0, 0, 0)

    points = [(-17.3944,0.399395,-88.9715), (-97.4837,4.2073,16.8826), (17.3948,0.383145,91.7289)
    ]
    for point in points:
        plotter.add_vtk_point(v, point, color=(0, 0, 1), radius=0.5)
    plotter.add_vtk_point(v, (0, 0, 0), color=(0, 0, 0), radius=5.0)

    plotter.add_vtk_text(v, "Axes showing origin", pos=(0.7, 0.1), font_size=48, color=(0, 0, 0))
    plotter.add_vtk_axes(v, origin=(0, 0, 0), length=50)
    v.view()