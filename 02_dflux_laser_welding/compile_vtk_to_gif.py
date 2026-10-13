"""
compile_vtk_to_gif.py
=====================
Multi-Viewport 3DEXPERIENCE / SIMULIA Animated Visual Report Compiler
for Circumferential Laser Welding Simulation.

Generates synchronized multi-viewport animated GIFs from exported VTK frames
(or from laser_therm_vtk.zip / laser_mech_vtk.zip).

Layout:
  1. Viewport 1 (Left)   : 3D Isometric Global View (European 1st Angle Projection)
  2. Viewport 2 (Center) : Top View (X-Y Plane along Z axis) - Laser Path & HAZ
  3. Viewport 3 (Right)  : 3D Cross-Section (R-Z Cutaway) - Penetration Depth & Melt Pool
"""

import os
import sys
import glob
import time
import shutil
import zipfile
import tempfile
import argparse
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import matplotlib.colors as mcolors

# -----------------------------------------------------------------------------
# 1. 3DEXPERIENCE Standard Colormap Palette & Color Bands
# -----------------------------------------------------------------------------
COLORMAP_3DX = [
    "#0000ff",  # Blue (Lowest)
    "#0055ff",
    "#00aaff",
    "#00ffff",  # Cyan
    "#00ffaa",
    "#00ff55",
    "#00ff00",  # Green
    "#aaff00",
    "#ffff00",  # Yellow
    "#ffaa00",
    "#ff5500",
    "#ff0000"   # Red (Highest)
]

def parse_args():
    parser = argparse.ArgumentParser(description="Compile VTK frames into multi-viewport GIF")
    parser.add_argument("input", nargs="?", default=None, help="Input VTK file, directory, or ZIP archive")
    parser.add_argument("--output", "-o", default=None, help="Output GIF file path")
    parser.add_argument("--field", "-f", default=None, choices=["auto", "Temperature", "Von_Mises", "PEEQ"],
                        help="Scalar field to visualize (default: auto-detect)")
    parser.add_argument("--fps", type=int, default=15, help="Animation playback FPS")
    parser.add_argument("--stride", "-s", type=int, default=2, help="Frame step stride (default: 2)")
    parser.add_argument("--scale", type=float, default=0.75, help="Output resolution scale factor (default: 0.75)")
    parser.add_argument("--max-val", type=float, default=None, help="Max scalar value for color normalization")
    parser.add_argument("--min-val", type=float, default=None, help="Min scalar value for color normalization")
    return parser.parse_args()

def extract_frame_metadata(vtk_path):
    time_val = 0.0
    step_name = ""
    try:
        with open(vtk_path, "r") as f:
            f.readline()
            line2 = f.readline()
            parts = line2.split()
            if "Time" in parts:
                t_idx = parts.index("Time")
                time_val = float(parts[t_idx + 1])
            if "Step" in parts:
                s_idx = parts.index("Step")
                step_name = parts[s_idx + 1]
    except Exception:
        pass
    return time_val, step_name

def locate_and_extract_vtk(input_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    if input_path:
        candidates.append(input_path)
    candidates.extend([
        os.path.join(script_dir, "disk_heatsource_th_vtk.zip"),
        os.path.join(script_dir, "laser_therm_vtk.zip"),
        os.path.join(script_dir, "disk_heatsource_me_vtk.zip"),
        os.path.join(script_dir, "laser_mech_vtk.zip"),
        os.path.join(script_dir, "gear_helical_vtk.zip")
    ])
    candidates.extend(sorted(glob.glob(os.path.join(script_dir, "*_vtk.zip"))))

    zip_file = None
    for cand in candidates:
        if cand and os.path.exists(cand) and cand.endswith(".zip"):
            zip_file = cand
            break

    temp_dir = None
    if zip_file:
        print("[INFO] Found VTK archive: %s" % zip_file)
        temp_dir = tempfile.mkdtemp(prefix="vtk_weld_")
        with zipfile.ZipFile(zip_file, 'r') as zf:
            zf.extractall(temp_dir)
        vtk_files = sorted(glob.glob(os.path.join(temp_dir, "*.vtk")))
    else:
        search_dir = input_path if (input_path and os.path.isdir(input_path)) else script_dir
        vtk_files = sorted(glob.glob(os.path.join(search_dir, "*.vtk")))

    return vtk_files, temp_dir, zip_file

def main():
    args = parse_args()
    vtk_files, temp_dir, zip_source = locate_and_extract_vtk(args.input)

    if not vtk_files:
        print("[ERROR] No VTK files found to process. Please provide an ODB or VTK zip.")
        sys.exit(1)

    total_available = len(vtk_files)
    if args.stride > 1:
        vtk_files = vtk_files[::args.stride]
        print("[INFO] Applied stride %d: using %d / %d frames." % (args.stride, len(vtk_files), total_available))
    else:
        print("[INFO] Processing all %d VTK frames." % len(vtk_files))

    # Inspect first frame to detect active scalar field
    sample_mesh = pv.read(vtk_files[0])
    point_arrays = list(sample_mesh.point_data.keys())
    cell_arrays = list(sample_mesh.cell_data.keys())
    print("[INFO] Point fields: %s" % point_arrays)
    print("[INFO] Cell fields: %s" % cell_arrays)

    active_field = args.field
    if active_field is None or active_field == "auto":
        if "Temperature" in point_arrays:
            active_field = "Temperature"
        elif "Von_Mises" in cell_arrays or "Von_Mises" in point_arrays:
            active_field = "Von_Mises"
        elif "PEEQ" in cell_arrays or "PEEQ" in point_arrays:
            active_field = "PEEQ"
        else:
            all_fields = point_arrays + cell_arrays
            active_field = all_fields[0] if all_fields else None

    if not active_field:
        print("[ERROR] No visualizable scalar field detected in VTK.")
        sys.exit(1)

    print("[INFO] Visualizing field: '%s'" % active_field)

    # Determine scalar bounds and colormap
    min_val = args.min_val if args.min_val is not None else 25.0
    if args.max_val is not None:
        max_val = args.max_val
    elif active_field == "Temperature":
        max_val = 2100.0  # Peak keyhole temperature
        # Smeared 5-Zone Metallurgical Colormap
        # Cor 1: Azul (Base metal <720 °C)
        # Cor 2: Verde (ZTA / HAZ 720 - 1485 °C)
        # Cor 3: Amarelo (Solidus -> Liquidus / Mushy 1485 - 1530 °C)
        # Cor 4: Laranja (Liquidus -> Vap / Melt pool 1530 - 1850 °C)
        # Cor 5: Vermelho (Keyhole / Vaporization >= 1850 °C)
        t_min, t_max = min_val, max_val
        norm_zta = (720.0 - t_min) / (t_max - t_min)
        norm_sol = (1485.0 - t_min) / (t_max - t_min)
        norm_liq = (1530.0 - t_min) / (t_max - t_min)
        norm_vap = (1850.0 - t_min) / (t_max - t_min)
        colors_list = [
            (0.0, '#1d4ed8'),            # Cor 1: Blue (Base metal <720)
            (norm_zta * 0.7, '#2563eb'), # Smooth blue transition
            (norm_zta, '#16a34a'),       # Cor 2: Green (ZTA onset Ac1 = 720 C)
            (norm_sol * 0.95, '#65a30d'),# Smooth Lime-Green
            (norm_sol, '#eab308'),       # Cor 3: Yellow (Solidus = 1485 C)
            (norm_liq, '#ea580c'),       # Cor 4: Orange (Liquidus = 1530 C)
            (norm_vap * 0.98, '#f97316'),# Smooth Orange
            (norm_vap, '#dc2626'),       # Cor 5: Red (Keyhole / Vap >= 1850 C)
            (1.0, '#991b1b')             # Deep Red Peak Core
        ]
        cmap = mcolors.LinearSegmentedColormap.from_list('smeared_thermal_5zones', colors_list)
    elif active_field == "Von_Mises":
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 850.0   # 16MnCr5 ultimate tensile ceiling
        cmap = COLORMAP_3DX
    elif active_field == "PEEQ":
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 0.15
        cmap = COLORMAP_3DX
    else:
        min_val = 0.0 if args.min_val is None else args.min_val
        max_val = 1.0
        cmap = COLORMAP_3DX

    # Output filename
    if args.output:
        out_gif = args.output
    elif zip_source:
        base = os.path.splitext(os.path.basename(zip_source))[0].replace("_vtk", "")
        out_gif = f"{base}_simulation.gif"
    else:
        out_gif = "laser_welding_simulation.gif"

    # Extract clean CAD boundary feature edges from frame 0 (merges coincident nodes to eliminate Z=15 interface line)
    print("[INFO] Computing CAD boundary feature edges...")
    mesh_0_clean = sample_mesh.clean(tolerance=1e-4)
    edges_cad = mesh_0_clean.extract_surface().extract_feature_edges(feature_angle=45, boundary_edges=True)
    print("[INFO] CAD boundary edges extracted (%d cells, exterior rims only)." % edges_cad.n_cells)

    # PyVista Off-screen Multi-Viewport Renderer
    pv.set_plot_theme("document")
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1800, 600))
    plotter.background_color = '#f1f5f9'

    # Font setup for HUD banner and legend
    banner_h = 44
    try:
        font_bold = ImageFont.truetype("arialbd.ttf", 16)
        font_reg = ImageFont.truetype("arial.ttf", 14)
        font_sm = ImageFont.truetype("arial.ttf", 11)
        font_sm_bd = ImageFont.truetype("arialbd.ttf", 11)
    except Exception:
        font_bold = font_reg = font_sm = font_sm_bd = ImageFont.load_default()

    # Pre-render 5-zone metallurgical legend card
    leg_w, leg_h = 240, 140
    card = Image.new("RGBA", (leg_w, leg_h), (255, 255, 255, 235))
    cdraw = ImageDraw.Draw(card)
    cdraw.rectangle([(0, 0), (leg_w - 1, leg_h - 1)], outline=(203, 213, 225), width=1)
    cdraw.text((10, 8), "Zonas Térmicas / Thermal Zones", fill=(15, 23, 42), font=font_sm_bd)

    zones_meta = [
        ("#dc2626", "T ≥ 1850 °C : Keyhole (Vapor)"),
        ("#ea580c", "1530–1850 °C : Poça Fusão (Melt)"),
        ("#eab308", "1485–1530 °C : Solidus → Liquidus"),
        ("#16a34a", "720–1485 °C : ZTA / HAZ (Ac1)"),
        ("#1d4ed8", "T < 720 °C : Metal Base (<Ac1)")
    ]
    for idx_z, (col, txt) in enumerate(zones_meta):
        yy = 32 + idx_z * 20
        cdraw.rectangle([(12, yy + 2), (26, yy + 14)], fill=col, outline=(51, 65, 85))
        cdraw.text((34, yy), txt, fill=(30, 41, 59), font=font_sm)

    rendered_images = []
    print(f"[INFO] Rendering {len(vtk_files)} frames to '{out_gif}'...")
    t0 = time.time()

    for idx, vtk_p in enumerate(vtk_files):
        sim_time, step_name = extract_frame_metadata(vtk_p)
        mesh = pv.read(vtk_p)
        # Clean mesh to merge coincident duplicate nodes across partition boundaries
        mesh_clean = mesh.clean(tolerance=1e-4)

        # Detect peak position for local melt pool tracking
        scalars_arr = mesh_clean.point_data.get(active_field, None)
        if scalars_arr is not None and len(scalars_arr) > 0:
            max_idx = np.argmax(scalars_arr)
            peak_val = float(scalars_arr[max_idx])
            peak_x, peak_y, peak_z = mesh_clean.points[max_idx]
            theta = np.arctan2(peak_y, peak_x)
            v_rad = np.array([np.cos(theta), np.sin(theta), 0.0])
        else:
            peak_val = 0.0
            peak_x, peak_y, peak_z = 0.0, 0.0, 20.0
            theta = 0.0
            v_rad = np.array([1.0, 0.0, 0.0])

        # Subplot 0: 3D Isometric Global View
        plotter.subplot(0, 0)
        plotter.clear_actors()
        plotter.add_mesh(mesh_clean, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, interpolate_before_map=True, show_scalar_bar=False)
        plotter.add_mesh(edges_cad, color='#1e293b', line_width=1.5)
        plotter.camera_position = [(80.0, -80.0, 60.0), (0.0, 0.0, 8.0), (0.0, 0.0, 1.0)]
        plotter.add_text("1. 3D Isometric Global", font_size=11, position="upper_left", color="black")

        # Subplot 1: Top View (XY)
        plotter.subplot(0, 1)
        plotter.clear_actors()
        plotter.add_mesh(mesh_clean, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, interpolate_before_map=True, show_scalar_bar=False)
        plotter.add_mesh(edges_cad, color='#1e293b', line_width=1.5)
        plotter.view_xy()
        plotter.camera.zoom(0.95)
        plotter.add_text("2. Top View (X-Y Seam)", font_size=11, position="upper_left", color="black")

        # Subplot 2: Lateral Melt Profile (Outer material removed to reveal depth & keyhole)
        plotter.subplot(0, 2)
        plotter.clear_actors()
        cut_outer = mesh_clean.clip(normal=v_rad, origin=[peak_x, peak_y, peak_z], invert=False)
        plotter.add_mesh(cut_outer, scalars=active_field, cmap=cmap, clim=[min_val, max_val],
                         show_edges=False, smooth_shading=True, interpolate_before_map=True,
                         show_scalar_bar=False)

        cam_target = np.array([peak_x, peak_y, 16.5])
        cam_eye = cam_target + v_rad * 10.5 + np.array([0, 0, 0.5])
        plotter.camera_position = [cam_eye, cam_target, (0.0, 0.0, 1.0)]
        plotter.add_text("3. Lateral Melt Profile (Keyhole & HAZ)", font_size=11, position="upper_left", color="black")

        img = plotter.screenshot(return_img=True)
        frame_im = Image.fromarray(img)

        # Composite HUD Header Banner
        banner = Image.new("RGB", (frame_im.width, frame_im.height + banner_h), (245, 247, 250))
        banner.paste(frame_im, (0, banner_h))
        draw = ImageDraw.Draw(banner)
        draw.text((18, 12), "ABAQUS / DFLUX: Circumferential Laser Welding (Thermal Analysis)", fill=(15, 23, 42), font=font_bold)
        unit_str = "°C" if active_field == "Temperature" else ""
        hud_info = f"Time: {sim_time:.3f} s  |  Peak: {peak_val:.1f} {unit_str}  |  Frame: {idx + 1} / {len(vtk_files)}"
        draw.text((frame_im.width - 480, 12), hud_info, fill=(30, 41, 59), font=font_reg)
        draw.line([(0, banner_h - 1), (frame_im.width, banner_h - 1)], fill=(226, 232, 240), width=1)

        # Overlay Metallurgical Legend Card on bottom-right
        leg_x = frame_im.width - leg_w - 20
        leg_y = banner.height - leg_h - 20
        banner.paste(card, (leg_x, leg_y), card)

        # Optional downscale for web/README optimization
        if args.scale and args.scale < 1.0:
            new_w = int(banner.width * args.scale)
            new_h = int(banner.height * args.scale)
            banner = banner.resize((new_w, new_h), Image.Resampling.LANCZOS)

        rendered_images.append(banner)

        if (idx + 1) % 10 == 0 or (idx + 1) == len(vtk_files):
            elapsed = time.time() - t0
            fps_render = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"[INFO] Rendered {idx + 1}/{len(vtk_files)} frames ({elapsed:.1f} s, {fps_render:.1f} fps)...")

    plotter.close()

    # Save animated GIF using Pillow
    if rendered_images:
        print(f"[INFO] Saving animated GIF ({len(rendered_images)} frames) at {args.fps} FPS...")
        duration_ms = int(1000.0 / args.fps)
        rendered_images[0].save(
            out_gif,
            save_all=True,
            append_images=rendered_images[1:],
            duration=duration_ms,
            loop=0,
            optimize=True
        )
        gif_size_mb = os.path.getsize(out_gif) / (1024.0 * 1024.0)
        print(f"[SUCCESS] Multi-viewport GIF generated: {out_gif} ({gif_size_mb:.2f} MB)")

    # Cleanup temp directory
    if temp_dir and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()
