"""
compile_vtk_to_gif.py
=====================
Three-viewport animation of the welded ring gear / differential case model
(travelling tooth load, DLOAD subroutine), built from the VTK frames exported
by odb_to_vtk.py. The plotted field is the nodal maximum principal stress
(average of the element-nodal values), which carries the surface stress used
for the fatigue evaluation.

  1. Panel 1 (left)   : Pinion view. The camera sits radially outside the ring
                        gear at the loaded tooth and orbits with the load, so
                        the meshing tooth always stays in the centre.
  2. Panel 2 (centre) : Front view along the Z axis. The gear stays still and
                        the loaded zone and its reactions travel the full 360 deg.
  3. Panel 3 (right)  : 90 deg section. Only a quarter of the ring gear and of
                        the differential case is kept, one cut face goes through
                        the meshing tooth, and section and camera follow the load.

Writes a stable animated GIF with one global indexed palette (no colour flicker).

Usage:
    python compile_vtk_to_gif.py                 # all frames -> gear_torque_simulation.gif
    python compile_vtk_to_gif.py --frames 0,8,16 # preview: PNG per selected frame
"""

import os
import sys
import glob
import time
import zipfile
import argparse
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import matplotlib.colors as mcolors

# -------------------------------------------------------------
# 0. Settings
# -------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

VTK_DIR = os.path.join(script_dir, "vtk_results")
VTK_ZIP = os.path.join(script_dir, "gear_helical_vtk.zip")
OUTPUT_GIF = os.path.join(script_dir, "gear_torque_simulation.gif")
FRAME_PATTERN = "*_[0-9][0-9][0-9][0-9].vtk"   # leaves out the fatigue cycle summary file

POINT_FIELD = "S_Max_Principal"    # nodal average, plotted in the three panels
LOAD_FIELD = "Von_Mises"           # element centroid, used only to locate the loaded tooth

T_RAMP = 0.20                 # s, end of the load ramp (load starts travelling after it)
STRESS_MAX = 100.0            # MPa, top of the colour scale (bottom is 0: compression is not coloured)
PANEL = 640                   # px, side of each square viewport
FOOTER = 72                   # px, strip below the viewports (legend, time, caption)
FPS = 6
CAPTION = "Ring gear welded to differential case (30% weld defect)  |  270 Nm as a travelling tooth load"

parser = argparse.ArgumentParser(description="Compile VTK frames into a 3-viewport GIF.")
parser.add_argument("--frames", default=None,
                    help="comma-separated frame indices to render as PNG previews")
args = parser.parse_args()

# -------------------------------------------------------------
# 1. Locate the VTK frames (extract the ZIP on first use)
# -------------------------------------------------------------
vtk_files = sorted(glob.glob(os.path.join(VTK_DIR, FRAME_PATTERN)))
if not vtk_files and os.path.exists(VTK_ZIP):
    print(f"[INFO] Extracting {os.path.basename(VTK_ZIP)} to {VTK_DIR} ...")
    with zipfile.ZipFile(VTK_ZIP, "r") as zf:
        zf.extractall(VTK_DIR)
    vtk_files = sorted(glob.glob(os.path.join(VTK_DIR, FRAME_PATTERN)))

if not vtk_files:
    print("[ERROR] No VTK files found.")
    sys.exit(1)

total_frames = len(vtk_files)
print(f"[INFO] VTK frames found: {total_frames}")

# -------------------------------------------------------------
# 2. Geometry (built once) and per-frame stress fields
# -------------------------------------------------------------
# The mesh is identical in every frame. The quadratic elements (C3D10, C3D20) are
# drawn as linear cells on their corner nodes, which keeps the original point
# numbering, so the nodal fields of every frame can be attached directly.
CORNER_NODES = {
    pv.CellType.QUADRATIC_TETRA: (pv.CellType.TETRA, 4),
    pv.CellType.QUADRATIC_HEXAHEDRON: (pv.CellType.HEXAHEDRON, 8),
    pv.CellType.QUADRATIC_WEDGE: (pv.CellType.WEDGE, 6),
    pv.CellType.TETRA: (pv.CellType.TETRA, 4),
    pv.CellType.HEXAHEDRON: (pv.CellType.HEXAHEDRON, 8),
    pv.CellType.WEDGE: (pv.CellType.WEDGE, 6),
}

print("[INFO] Reading mesh geometry ...")
first = pv.read(vtk_files[0])
points = np.asarray(first.points)

cell_blocks, type_blocks, cell_order = [], [], []
for cell_type, connectivity in first.cells_dict.items():
    linear_type, n_corners = CORNER_NODES[pv.CellType(cell_type)]
    n = len(connectivity)
    cell_blocks.append(np.hstack([np.full((n, 1), n_corners), connectivity[:, :n_corners]]).ravel())
    type_blocks.append(np.full(n, linear_type, dtype=np.uint8))
    cell_order.append(np.flatnonzero(first.celltypes == cell_type))
    print(f"[INFO]   {pv.CellType(cell_type).name:<22s} {n} elements")

mesh = pv.UnstructuredGrid(np.concatenate(cell_blocks), np.concatenate(type_blocks), points)
cell_order = np.concatenate(cell_order)     # original cell index of every cell of `mesh`
centroids = mesh.cell_centers().points
del first


def read_time(fpath):
    with open(fpath, "rb") as f:
        f.readline()
        parts = f.readline().decode("ascii", errors="ignore").split()
    return float(parts[parts.index("Time") + 1])


print("[INFO] Reading stress fields ...")
times = np.array([read_time(f) for f in vtk_files])
point_stress, load_stress = [], []
for fpath in vtk_files:
    frame = pv.read(fpath)
    point_stress.append(np.asarray(frame.point_data[POINT_FIELD], dtype=np.float32))
    load_stress.append(np.asarray(frame.cell_data[LOAD_FIELD], dtype=np.float32)[cell_order])
    del frame

print(f"[INFO] Peak nodal maximum principal stress over the whole run: "
      f"{max(float(s.max()) for s in point_stress):.1f} MPa")

# -------------------------------------------------------------
# 3. Kinematics: angular position of the travelling load
# -------------------------------------------------------------
# The rotation axis is Z. The loaded tooth is located in each frame from the most
# stressed elements of the toothed rim, then fitted with a constant angular speed.
rim = np.hypot(centroids[:, 0], centroids[:, 1]) > 75.0
rim_angle = np.arctan2(centroids[rim, 1], centroids[rim, 0])

raw_phi = np.full(total_frames, np.nan)
for i, vm in enumerate(load_stress):
    vm_rim = vm[rim]
    if vm_rim.max() > 1.0:
        hot = vm_rim >= np.percentile(vm_rim, 99.8)
        raw_phi[i] = np.arctan2(np.sin(rim_angle[hot]).mean(), np.cos(rim_angle[hot]).mean())
del load_stress

rolling = (times >= T_RAMP - 1e-6) & ~np.isnan(raw_phi)
phi_unwrapped = np.unwrap(raw_phi[rolling])
omega, phi_ramp = np.polyfit(times[rolling] - T_RAMP, phi_unwrapped, 1)
phi = phi_ramp + omega * np.clip(times - T_RAMP, 0.0, None)
travel_deg = np.abs(np.rad2deg(phi - phi_ramp))
fit_err = np.rad2deg(np.abs(phi_unwrapped - (phi_ramp + omega * (times[rolling] - T_RAMP)))).max()

print(f"[KINEMATICS] Load starts at {np.rad2deg(phi_ramp):.1f} deg and travels at "
      f"{np.rad2deg(omega):.1f} deg/s ({travel_deg[-1]:.1f} deg in total)")
print(f"[KINEMATICS] Largest deviation from constant speed: {fit_err:.1f} deg")

# -------------------------------------------------------------
# 4. Colour scale and legend
# -------------------------------------------------------------
bg_color = "#e0e2e7"

# 10 bands of 10 MPa, from 0 MPa (blue) to 100 MPa (red)
band_colors = [
    "#0000ff", "#005cff", "#00b9ff", "#00ffe7", "#00ff8b",
    "#2eff00", "#e7ff00", "#ffb900", "#ff5c00", "#ff0000"
]
n_bands = len(band_colors)
stress_cmap = mcolors.ListedColormap(band_colors, name="stress_bands")
mesh_kwargs = dict(scalars="Stress", cmap=stress_cmap, clim=[0.0, STRESS_MAX],
                   n_colors=n_bands, show_edges=False, show_scalar_bar=False)


def load_font(size, bold=False):
    try:
        return ImageFont.truetype("arialbd.ttf" if bold else "arial.ttf", size)
    except OSError:
        return ImageFont.load_default()


font_title = load_font(17, bold=True)
font_sub = load_font(14)
font_legend = load_font(12)
font_legend_bold = load_font(12, bold=True)


def create_legend_strip():
    bw, bh = 36, 16
    strip = Image.new("RGB", (n_bands * bw + 60, FOOTER), bg_color)
    draw = ImageDraw.Draw(strip)
    draw.text((14, 8), "Maximum principal stress (MPa)", fill="black", font=font_legend_bold)
    bx, by = 14, 28
    step = STRESS_MAX / n_bands
    for i, col in enumerate(band_colors):
        x0 = bx + i * bw
        draw.rectangle([x0, by, x0 + bw, by + bh], fill=col, outline="black", width=1)
        label = "≤ 0" if i == 0 else f"{i * step:.0f}"
        draw.text((x0, by + bh + 4), label, fill="black", font=font_legend, anchor="ma")
    draw.text((bx + n_bands * bw, by + bh + 4), f"≥ {STRESS_MAX:.0f}", fill="black",
              font=font_legend, anchor="ma")
    return strip


legend_strip = create_legend_strip()

# -------------------------------------------------------------
# 5. Static geometry shared by panels 1 and 2
# -------------------------------------------------------------
surface = mesh.extract_surface(algorithm="dataset_surface")
surface_point_ids = surface.point_data["vtkOriginalPointIds"]
surface_edges = surface.extract_feature_edges(boundary_edges=True, feature_edges=True,
                                              manifold_edges=False, feature_angle=30.0)
r_out = np.hypot(points[:, 0], points[:, 1]).max()


def render_frame(idx):
    ph = phi[idx]
    e_r = np.array([np.cos(ph), np.sin(ph), 0.0])     # radial direction of the loaded tooth
    e_t = np.array([-np.sin(ph), np.cos(ph), 0.0])    # tangential direction

    surface.point_data["Stress"] = point_stress[idx][surface_point_ids]

    # 90 deg sector for panel 3: cut faces on the loaded tooth and 90 deg away from it
    mesh.point_data["Stress"] = point_stress[idx]
    sector = mesh.clip(normal=e_t, origin=(0, 0, 0), invert=False)
    sector = sector.clip(normal=e_r, origin=(0, 0, 0), invert=False).extract_surface(algorithm="dataset_surface")
    sector_edges = sector.extract_feature_edges(boundary_edges=True, feature_edges=True,
                                                manifold_edges=False, feature_angle=30.0)

    plotter = pv.Plotter(shape=(1, 3), window_size=[3 * PANEL, PANEL], off_screen=True, border=False)
    plotter.set_background(bg_color)

    # ---------------------------------------------------------
    # PANEL 1: PINION VIEW (camera orbits with the load)
    # ---------------------------------------------------------
    plotter.subplot(0, 0)
    plotter.add_mesh(surface, preference="point", **mesh_kwargs)
    plotter.add_mesh(surface_edges, color="black", line_width=1.0)
    focal = r_out * e_r + np.array([0.0, 0.0, 4.0])
    plotter.camera.focal_point = focal
    plotter.camera.position = focal + 400.0 * e_r
    plotter.camera.up = (0.0, 0.0, 1.0)
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = 44.0

    # ---------------------------------------------------------
    # PANEL 2: FRONT VIEW ALONG Z (fixed camera)
    # ---------------------------------------------------------
    plotter.subplot(0, 1)
    plotter.add_mesh(surface, preference="point", **mesh_kwargs)
    plotter.add_mesh(surface_edges, color="black", line_width=1.0)
    plotter.camera.focal_point = (0.0, 0.0, 0.0)
    plotter.camera.position = (0.0, 0.0, 500.0)
    plotter.camera.up = (0.0, 1.0, 0.0)
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = 1.20 * r_out

    # ---------------------------------------------------------
    # PANEL 3: 90 DEG SECTION (follows the load)
    # ---------------------------------------------------------
    plotter.subplot(0, 2)
    plotter.add_mesh(sector, preference="point", **mesh_kwargs)
    plotter.add_mesh(sector_edges, color="black", line_width=1.0)
    azimuth = ph - np.deg2rad(115.0)
    elevation = np.deg2rad(25.0)
    view_dir = np.array([np.cos(azimuth) * np.cos(elevation),
                         np.sin(azimuth) * np.cos(elevation),
                         np.sin(elevation)])
    focal = 0.50 * r_out * e_r + 0.42 * r_out * e_t + np.array([0.0, 0.0, -8.0])
    plotter.camera.focal_point = focal
    plotter.camera.position = focal + 400.0 * view_dir
    plotter.camera.up = (0.0, 0.0, 1.0)
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = 76.0

    raw_img = plotter.screenshot(return_img=True)
    plotter.close()

    # Titles, legend and caption drawn with PIL (identical pixels in every frame)
    img = Image.new("RGB", (3 * PANEL, PANEL + FOOTER), bg_color)
    img.paste(Image.fromarray(raw_img), (0, 0))
    draw = ImageDraw.Draw(img)
    titles = [
        ("1. Pinion view (camera follows the load)", "Maximum principal stress, nodal"),
        ("2. Front view along Z (fixed gear)", "Maximum principal stress, nodal"),
        ("3. 90° section (follows the load)", "Maximum principal stress, nodal"),
    ]
    for k, (title, sub) in enumerate(titles):
        draw.text((k * PANEL + 18, 14), title, fill="black", font=font_title)
        draw.text((k * PANEL + 18, 38), sub, fill="#444444", font=font_sub)
        if k:
            draw.line([(k * PANEL, 0), (k * PANEL, PANEL)], fill="#9a9da5", width=1)
    draw.line([(0, PANEL), (3 * PANEL, PANEL)], fill="#9a9da5", width=1)

    if times[idx] < T_RAMP - 1e-6:
        stage = "load ramp"
    else:
        stage = f"load position {travel_deg[idx]:.1f}°"
    img.paste(legend_strip, (4, PANEL + 1))
    draw.text((1.5 * PANEL, PANEL + 26), f"t = {times[idx]:.2f} s   |   {stage}",
              fill="black", font=font_title, anchor="ma")
    draw.text((3 * PANEL - 18, PANEL + 28), CAPTION, fill="black", font=font_sub, anchor="ra")
    return img


# -------------------------------------------------------------
# 6. Render the frames
# -------------------------------------------------------------
if args.frames:
    for idx in [int(s) for s in args.frames.split(",")]:
        out_png = os.path.join(script_dir, f"preview_{idx:04d}.png")
        render_frame(idx).save(out_png)
        print(f"[INFO] Preview written: {out_png}")
    sys.exit(0)

pil_frames = []
t_proc_start = time.time()
for idx in range(total_frames):
    pil_frames.append(render_frame(idx))
    if (idx + 1) % 5 == 0 or (idx + 1) == total_frames:
        print(f"  Frame {idx + 1}/{total_frames} rendered ({time.time() - t_proc_start:.0f} s)")

# -------------------------------------------------------------
# 7. Global palette and GIF export
# -------------------------------------------------------------
# One palette for the whole animation, sampled from frames spread along the run,
# with the exact band colours locked in the first entries.
print("[INFO] Building the global palette ...")
sample_ids = np.linspace(0, total_frames - 1, min(6, total_frames)).astype(int)
sample_img = Image.new("RGB", (3 * PANEL, (PANEL + FOOTER) * len(sample_ids)))
for row, f_idx in enumerate(sample_ids):
    sample_img.paste(pil_frames[f_idx], (0, row * (PANEL + FOOTER)))

fixed_palette_hex = band_colors + ["#000000", "#ffffff", bg_color, "#444444", "#9a9da5"]
raw_palette = list(sample_img.quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                                       dither=Image.Dither.NONE).getpalette())
for i, hex_color in enumerate(fixed_palette_hex):
    raw_palette[i * 3: i * 3 + 3] = [int(round(255 * c)) for c in mcolors.to_rgb(hex_color)]

master_palette_image = Image.new("P", (1, 1))
master_palette_image.putpalette(raw_palette)
quantized_frames = [f.quantize(palette=master_palette_image, dither=Image.Dither.NONE)
                    for f in pil_frames]

print(f"[INFO] Writing {os.path.basename(OUTPUT_GIF)} at {FPS} FPS ...")
quantized_frames[0].save(OUTPUT_GIF, save_all=True, append_images=quantized_frames[1:],
                         duration=int(1000.0 / FPS), loop=0, optimize=False)

print(f"\n[DONE] {OUTPUT_GIF} ({os.path.getsize(OUTPUT_GIF) / 1e6:.1f} MB, {total_frames} frames)")
