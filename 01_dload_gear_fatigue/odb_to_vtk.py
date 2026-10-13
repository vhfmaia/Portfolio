# -*- coding: utf-8 -*-
"""
Abaqus ODB to VTK Exporter (legacy VTK, binary)
Compatible with Abaqus Python 2.7 and 3.x.

Exports every frame of every step to gear_helical_NNNN.vtk with the fields
needed for the fatigue evaluation of the drive + coast load cycle:

  Point data (nodal, averaged from the element-nodal stress tensor):
    Displacement, S_Mises, S_Max_Principal, S_Min_Principal
  Cell data (element centroid):
    Von_Mises, Max_Principal, Min_Principal, PEEQ

The signed principal stresses are what separates the tensile from the
compressive half of the tooth-root cycle when the torque reverses; the nodal
fields carry the surface values at the root fillet and at the weld.

Von Mises is unsigned, so it cannot give a stress range. The fatigue quantities
are therefore written once for the whole load cycle to
gear_helical_fatigue_cycle.vtk (nodal):

    Peak_Max_Principal      highest maximum principal stress of the cycle
                            (local tooth-root stress, ISO 6336)
    Lowest_Min_Principal    lowest minimum principal stress of the cycle
    Principal_Stress_Range  largest principal stress of the difference between
                            the stress tensor at the peak and the tensor of any
                            other frame (weld notch stress range, IIW). It stays
                            valid when the principal directions rotate as the
                            load travels.

All solid continuum elements are exported (C3D4, C3D6, C3D8, C3D10, C3D15, C3D20
and their R / H / M / I variants); any other element type is reported.
"""
import os
import sys
import glob

import numpy as np

try:
    from odbAccess import openOdb
    from abaqusConstants import CENTROID, ELEMENT_NODAL
except ImportError:  # allows the writer to be imported outside Abaqus
    openOdb = None
    CENTROID = ELEMENT_NODAL = None

# Abaqus element family to VTK cell type (node ordering is identical). Matched as a
# prefix, so C3D20, C3D20R, C3D20RH, C3D10M, C3D8I, ... all resolve to their family.
VTK_CELL_FAMILIES = (
    ('C3D10', 24),   # VTK_QUADRATIC_TETRA
    ('C3D15', 26),   # VTK_QUADRATIC_WEDGE
    ('C3D20', 25),   # VTK_QUADRATIC_HEXAHEDRON
    ('C3D4', 10),    # VTK_TETRA
    ('C3D6', 13),    # VTK_WEDGE
    ('C3D8', 12)     # VTK_HEXAHEDRON
)


def vtk_cell_type(abaqus_type):
    """VTK cell type of an Abaqus element type, or None if it is not a supported solid."""
    name = abaqus_type.upper()
    for family, cell_type in VTK_CELL_FAMILIES:
        if name.startswith(family):
            return cell_type
    return None


def stress_invariants(tensor):
    """Von Mises, maximum and minimum principal stress from (n, 6) rows S11, S22, S33, S12, S13, S23."""
    t = np.asarray(tensor, dtype=np.float64)
    s11, s22, s33, s12, s13, s23 = [t[:, k] for k in range(6)]
    mises = np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                    + 3.0 * (s12 ** 2 + s13 ** 2 + s23 ** 2))
    matrix = np.empty((t.shape[0], 3, 3))
    matrix[:, 0, 0] = s11
    matrix[:, 1, 1] = s22
    matrix[:, 2, 2] = s33
    matrix[:, 0, 1] = matrix[:, 1, 0] = s12
    matrix[:, 0, 2] = matrix[:, 2, 0] = s13
    matrix[:, 1, 2] = matrix[:, 2, 1] = s23
    principal = np.linalg.eigvalsh(matrix)  # ascending
    return mises, principal[:, 2], principal[:, 0]


def cycle_fatigue_fields(tensor_file, num_nodes):
    """Fatigue quantities of the whole load cycle from the nodal stress tensors of every frame.

    `tensor_file` holds one (num_nodes, 6) float32 block per frame, in frame order.
    """
    tensors = np.memmap(tensor_file, dtype=np.float32, mode='r').reshape(-1, num_nodes, 6)

    # Envelope of the principal stresses and the tensor at the tensile peak of each node
    peak_max = -np.inf * np.ones(num_nodes)
    lowest_min = np.inf * np.ones(num_nodes)
    peak_tensor = np.zeros((num_nodes, 6))
    for frame_tensor in tensors:
        _, s_max, s_min = stress_invariants(frame_tensor)
        higher = s_max > peak_max
        peak_max[higher] = s_max[higher]
        peak_tensor[higher] = frame_tensor[higher]
        lowest_min = np.minimum(lowest_min, s_min)

    # Stress range: largest principal stress of (peak tensor - frame tensor) over the cycle
    stress_range = np.zeros(num_nodes)
    for frame_tensor in tensors:
        _, d_max, d_min = stress_invariants(peak_tensor - frame_tensor)
        stress_range = np.maximum(stress_range, np.maximum(np.abs(d_max), np.abs(d_min)))

    del tensors
    return [("Peak_Max_Principal", peak_max), ("Lowest_Min_Principal", lowest_min),
            ("Principal_Stress_Range", stress_range)]


def write_vtk(filename, title, points, cells, cell_types, point_data, cell_data):
    """Writes a legacy binary VTK unstructured grid. `cells` is the flat [n, i0, i1, ...] connectivity."""
    def text(line):
        return (line + "\n").encode('ascii')

    with open(filename, 'wb') as f:
        f.write(text("# vtk DataFile Version 3.0"))
        f.write(text(title[:255]))
        f.write(text("BINARY"))
        f.write(text("DATASET UNSTRUCTURED_GRID"))

        f.write(text("POINTS %d float" % len(points)))
        np.asarray(points, dtype='>f4').tofile(f)
        f.write(text(""))

        f.write(text("CELLS %d %d" % (len(cell_types), len(cells))))
        np.asarray(cells, dtype='>i4').tofile(f)
        f.write(text(""))
        f.write(text("CELL_TYPES %d" % len(cell_types)))
        np.asarray(cell_types, dtype='>i4').tofile(f)
        f.write(text(""))

        for header, count, fields in (("POINT_DATA", len(points), point_data),
                                      ("CELL_DATA", len(cell_types), cell_data)):
            if not fields:
                continue
            f.write(text("%s %d" % (header, count)))
            for name, values in fields:
                values = np.asarray(values, dtype='>f4')
                if values.ndim == 2:
                    f.write(text("VECTORS %s float" % name))
                else:
                    f.write(text("SCALARS %s float 1" % name))
                    f.write(text("LOOKUP_TABLE default"))
                values.tofile(f)
                f.write(text(""))


def block_data(block):
    """Bulk data of a field block for single or double precision output databases."""
    try:
        return np.asarray(block.data)
    except Exception:
        return np.asarray(block.dataDouble)


def in_instance(block, instance):
    owner = getattr(block, 'instance', None)
    return owner is None or owner.name == instance.name


def extract_fields(frame, instance, node_index, elem_index, num_nodes, num_cells):
    """Bulk (numpy) extraction of displacements and stresses of one frame.

    Also returns the nodal stress tensor (None if the frame has no stress output).
    """
    point_data, cell_data = [], []
    nodal_tensor = None
    outputs = frame.fieldOutputs

    if 'U' in outputs:
        disp = np.zeros((num_nodes, 3))
        for block in outputs['U'].bulkDataBlocks:
            if not in_instance(block, instance):
                continue
            idx = node_index[np.asarray(block.nodeLabels)]
            keep = idx >= 0
            disp[idx[keep]] = block_data(block)[keep][:, :3]
        point_data.append(("Displacement", disp))

    if 'S' in outputs:
        stress = outputs['S']

        # Nodal stress: average of the element-nodal tensors, then invariants
        total = np.zeros((num_nodes, 6))
        hits = np.zeros(num_nodes)
        for block in stress.getSubset(position=ELEMENT_NODAL).bulkDataBlocks:
            data = block_data(block)
            if not in_instance(block, instance) or data.shape[1] != 6:
                continue
            idx = node_index[np.asarray(block.nodeLabels)]
            keep = idx >= 0
            idx, data = idx[keep], data[keep]
            hits += np.bincount(idx, minlength=num_nodes)
            for k in range(6):
                total[:, k] += np.bincount(idx, weights=data[:, k], minlength=num_nodes)
        total /= np.maximum(hits, 1.0)[:, None]
        nodal_tensor = total
        mises, s_max, s_min = stress_invariants(total)
        point_data += [("S_Mises", mises), ("S_Max_Principal", s_max), ("S_Min_Principal", s_min)]

        # Element centroid stress
        centroid = np.zeros((num_cells, 6))
        for block in stress.getSubset(position=CENTROID).bulkDataBlocks:
            data = block_data(block)
            if not in_instance(block, instance) or data.shape[1] != 6:
                continue
            idx = elem_index[np.asarray(block.elementLabels)]
            keep = idx >= 0
            centroid[idx[keep]] = data[keep]
        mises, s_max, s_min = stress_invariants(centroid)
        cell_data += [("Von_Mises", mises), ("Max_Principal", s_max), ("Min_Principal", s_min)]

    if 'PEEQ' in outputs:
        peeq = np.zeros(num_cells)
        for block in outputs['PEEQ'].getSubset(position=CENTROID).bulkDataBlocks:
            if not in_instance(block, instance):
                continue
            idx = elem_index[np.asarray(block.elementLabels)]
            keep = idx >= 0
            peeq[idx[keep]] = block_data(block)[keep].reshape(-1)
        cell_data.append(("PEEQ", peeq))

    return point_data, cell_data, nodal_tensor


def extract_fields_legacy(frame, node_index, elem_index, num_nodes, num_cells):
    """Value-by-value fallback: displacements and integration-point averaged von Mises."""
    point_data, cell_data = [], []
    outputs = frame.fieldOutputs

    if 'U' in outputs:
        disp = np.zeros((num_nodes, 3))
        for val in outputs['U'].values:
            idx = node_index[val.nodeLabel]
            if idx >= 0:
                disp[idx] = val.data[:3]
        point_data.append(("Displacement", disp))

    if 'S' in outputs:
        total = np.zeros(num_cells)
        hits = np.zeros(num_cells)
        for val in outputs['S'].values:
            idx = elem_index[val.elementLabel]
            if idx >= 0:
                total[idx] += val.mises
                hits[idx] += 1.0
        cell_data.append(("Von_Mises", total / np.maximum(hits, 1.0)))

    return point_data, cell_data


def export_odb_to_vtk(odb_name, step_names=None, output_prefix="gear_vtk"):
    if not os.path.exists(odb_name):
        print("[ERROR] ODB file not found: %s" % odb_name)
        return

    print("[INFO] Opening ODB: %s" % odb_name)
    odb = openOdb(odb_name, readOnly=True)

    # 1. Select analysis steps (all by default)
    step_keys = list(odb.steps.keys())
    if not step_keys:
        print("[ERROR] No steps found in ODB.")
        odb.close()
        return
    if step_names:
        step_keys = [key for key in step_keys if key in step_names]
    print("[INFO] Processing Steps: %s" % ", ".join(step_keys))

    # 2. Select instance / assembly
    assembly = odb.rootAssembly
    if len(assembly.instances) > 0:
        instance = list(assembly.instances.values())[0]
    else:
        instance = assembly

    nodes = instance.nodes
    elements = instance.elements

    # Node labels to 0-based VTK point indices
    node_labels = np.array([node.label for node in nodes], dtype=np.int64)
    points = np.array([node.coordinates for node in nodes], dtype=np.float64)
    num_nodes = len(node_labels)
    node_index = -np.ones(int(node_labels.max()) + 1, dtype=np.int64)
    node_index[node_labels] = np.arange(num_nodes)

    # Mesh cells and connectivity
    cells = []
    cell_types = []
    cell_labels = []
    type_map = {}
    exported, skipped = {}, {}
    for elem in elements:
        e_type = elem.type
        if e_type not in type_map:
            type_map[e_type] = vtk_cell_type(e_type)
        c_type = type_map[e_type]
        if c_type is None:
            skipped[e_type] = skipped.get(e_type, 0) + 1
            continue
        exported[e_type] = exported.get(e_type, 0) + 1
        conn = elem.connectivity
        cells.append(len(conn))
        cells.extend(int(node_index[label]) for label in conn)
        cell_types.append(c_type)
        cell_labels.append(elem.label)

    num_cells = len(cell_types)
    elem_index = -np.ones(max(cell_labels) + 1, dtype=np.int64)
    elem_index[np.array(cell_labels, dtype=np.int64)] = np.arange(num_cells)
    cells = np.array(cells, dtype=np.int64)
    print("[INFO] Mesh extracted: %d nodes, %d elements." % (num_nodes, num_cells))
    for e_type in sorted(exported):
        print("[INFO]   %-8s %d elements" % (e_type, exported[e_type]))
    for e_type in sorted(skipped):
        print("[WARNING] %-8s %d elements NOT exported (unsupported element type)" % (e_type, skipped[e_type]))

    # 3. Process field outputs frame by frame, across all steps
    f_idx = 0
    use_bulk = True
    tensor_filename = "%s_nodal_tensors.tmp" % output_prefix
    tensor_file = open(tensor_filename, 'wb')
    tensor_frames = 0
    for s_idx, step_name in enumerate(step_keys):
        step = odb.steps[step_name]
        for local_idx, frame in enumerate(step.frames):
            if s_idx > 0 and local_idx == 0:
                continue  # same state as the last frame of the previous step
            total_time = step.totalTime + frame.frameValue
            vtk_filename = "%s_%04d.vtk" % (output_prefix, f_idx)
            print("[INFO] Exporting %s frame %d (Time = %.4f s) -> %s"
                  % (step_name, local_idx, total_time, vtk_filename))

            if use_bulk:
                try:
                    point_data, cell_data, nodal_tensor = extract_fields(frame, instance, node_index,
                                                                         elem_index, num_nodes, num_cells)
                    if nodal_tensor is not None:
                        np.asarray(nodal_tensor, dtype=np.float32).tofile(tensor_file)
                        tensor_frames += 1
                except Exception as error:
                    print("[WARNING] Bulk extraction failed (%s). Falling back to value-by-value export." % error)
                    use_bulk = False
            if not use_bulk:
                point_data, cell_data = extract_fields_legacy(frame, node_index, elem_index, num_nodes, num_cells)

            title = "Abaqus ODB Export Frame %d Time %f Step %s StepTime %f" % (
                f_idx, total_time, step_name, frame.frameValue)
            write_vtk(vtk_filename, title, points, cells, cell_types, point_data, cell_data)
            f_idx += 1

    odb.close()
    tensor_file.close()

    # 4. Fatigue quantities of the whole load cycle
    if use_bulk and tensor_frames > 1:
        cycle_filename = "%s_fatigue_cycle.vtk" % output_prefix
        print("[INFO] Evaluating the load cycle over %d frames -> %s" % (tensor_frames, cycle_filename))
        cycle_data = cycle_fatigue_fields(tensor_filename, num_nodes)
        for name, values in cycle_data:
            print("[INFO]   %-24s %.1f to %.1f MPa" % (name, values.min(), values.max()))
        write_vtk(cycle_filename, "Abaqus ODB Export Fatigue Cycle %d frames" % tensor_frames,
                  points, cells, cell_types, cycle_data, [])
    else:
        print("[WARNING] Fatigue cycle file not written (nodal stress tensor not available).")
    os.remove(tensor_filename)

    print("[INFO] Export completed successfully: %d frames." % f_idx)


if __name__ == "__main__":
    target_odb = "Gear_Shaft_assembly_drive_coast.odb"
    if len(sys.argv) > 1 and sys.argv[1].endswith(".odb"):
        target_odb = sys.argv[1]
    elif not os.path.exists(target_odb):
        found_odbs = [f for f in glob.glob("*.odb") if "assembly" in f.lower()] or glob.glob("*.odb")
        if found_odbs:
            target_odb = found_odbs[0]
            print("[INFO] Target ODB not specified, auto-detected: %s" % target_odb)

    export_odb_to_vtk(target_odb, output_prefix="gear_helical")
