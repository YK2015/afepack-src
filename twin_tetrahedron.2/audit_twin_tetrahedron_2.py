#!/usr/bin/env python3

"""Audit twin_tetrahedron.2 basis functions with patch tests.

The oracle here is not a second implementation of the twin macro basis.
Instead, it uses invariants, global quadratic reproduction, and ordinary
tetrahedron P2 interpolation on each sub-tetrahedron.
"""

import ctypes
import math
import random
from pathlib import Path


DoubleP = ctypes.POINTER(ctypes.c_double)
TOL_VALUE = 2.0e-11
TOL_GRAD = 2.0e-10
TOL_FD = 2.0e-6

VERTICES = [
    (0.0, 0.0, 0.0),  # v0
    (1.0, 0.0, 0.0),  # v1
    (0.5, 0.5, 0.0),  # v2
    (0.0, 1.0, 0.0),  # v3
    (0.0, 0.0, 1.0),  # v4
]

# Current AFEPack DOF/basis order after moving the two internal micro-edge
# DOFs to face 0 and face 3.
BASIS_NAMES = (
    [f"phi_{i}" for i in range(1, 6)]
    + [f"psi_{i}" for i in range(1, 8)]
    + ["chi_2", "chi_1"]
)
GRADIENT_NAMES = (
    [f"gradient_phi_{i}" for i in range(1, 6)]
    + [f"gradient_psi_{i}" for i in range(1, 8)]
    + ["gradient_chi_2", "gradient_chi_1"]
)

DOF_POINTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.5, 0.0, 0.0),
    (0.0, 0.5, 0.0),
    (0.0, 0.0, 0.5),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
    (0.75, 0.25, 0.0),
    (0.25, 0.75, 0.0),
    (0.25, 0.25, 0.5),  # chi_2, micro-edge (v2,v4), face 0
    (0.25, 0.25, 0.0),  # chi_1, micro-edge (v0,v2), face 3
]

SUB_TETS = {
    "T+": (0, 1, 2, 4),
    "T-": (0, 2, 3, 4),
}

# For each sub-tetrahedron, list the global DOF indices in ordinary tetra P2
# order: vertices, then edges (01, 02, 03, 12, 13, 23).
LOCAL_DOF_MAP = {
    "T+": [0, 1, 2, 4, 5, 13, 7, 10, 9, 12],
    "T-": [0, 2, 3, 4, 13, 6, 7, 11, 12, 8],
}

SUPPORT_ZERO = {
    "T+": [3, 6, 8, 11],   # phi_4, psi_2, psi_4, psi_7
    "T-": [1, 5, 9, 10],   # phi_2, psi_1, psi_5, psi_6
}

FACE_ZERO = {
    # face 0: (v1,v2,v3,v4), no v0 influence in the trace
    "F0": [0, 5, 6, 7, 13],
    # face 1: (v0,v3,v4), only vertices 0,3,4 and their edges may remain
    "F1": [1, 2, 5, 9, 10, 11, 12, 13],
    # face 2: (v0,v1,v4), only vertices 0,1,4 and their edges may remain
    "F2": [2, 3, 6, 8, 10, 11, 12, 13],
    # face 3: (v0,v1,v2,v3), no v4 influence in the trace
    "F3": [4, 7, 8, 9, 12],
}

FACE_SAMPLERS = {
    "F0": ((1, 2, 4), (2, 3, 4)),
    "F1": ((0, 3, 4),),
    "F2": ((0, 1, 4),),
    "F3": ((0, 1, 2), (0, 2, 3)),
}


def load_function(lib, name):
    fn = getattr(lib, name)
    fn.argtypes = [DoubleP, ctypes.POINTER(DoubleP), ctypes.c_void_p]
    fn.restype = None
    return fn


def make_vertex_ptrs():
    arrays = [(ctypes.c_double * 3)(*v) for v in VERTICES]
    ptrs = (DoubleP * len(arrays))(*[ctypes.cast(item, DoubleP) for item in arrays])
    return arrays, ptrs


def call_scalar(fn, point, vertex_ptrs):
    point_array = (ctypes.c_double * 3)(*point)
    output = (ctypes.c_double * 1)(0.0)
    fn(ctypes.cast(point_array, DoubleP), vertex_ptrs, ctypes.byref(output))
    return output[0]


def call_vector(fn, point, vertex_ptrs):
    point_array = (ctypes.c_double * 3)(*point)
    output = (ctypes.c_double * 3)(0.0, 0.0, 0.0)
    fn(ctypes.cast(point_array, DoubleP), vertex_ptrs, ctypes.byref(output))
    return (output[0], output[1], output[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(c, a):
    return (c * a[0], c * a[1], c * a[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a):
    return math.sqrt(dot(a, a))


def max_abs_vec(a):
    return max(abs(a[0]), abs(a[1]), abs(a[2]))


def mat_vec(mat, vec):
    return tuple(sum(mat[i][j] * vec[j] for j in range(3)) for i in range(3))


def inverse3(mat):
    aug = [[float(mat[i][j]) for j in range(3)] + [1.0 if i == j else 0.0 for j in range(3)]
           for i in range(3)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1.0e-14:
            raise ValueError("singular 3x3 matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        factor = aug[col][col]
        aug[col] = [x / factor for x in aug[col]]
        for row in range(3):
            if row == col:
                continue
            factor = aug[row][col]
            aug[row] = [aug[row][k] - factor * aug[col][k] for k in range(6)]
    return [row[3:] for row in aug]


def tetra_data(tet_name):
    ids = SUB_TETS[tet_name]
    verts = [VERTICES[i] for i in ids]
    v0, v1, v2, v3 = verts
    mat = [
        [v1[0] - v0[0], v2[0] - v0[0], v3[0] - v0[0]],
        [v1[1] - v0[1], v2[1] - v0[1], v3[1] - v0[1]],
        [v1[2] - v0[2], v2[2] - v0[2], v3[2] - v0[2]],
    ]
    inv = inverse3(mat)
    grad_l1 = (inv[0][0], inv[0][1], inv[0][2])
    grad_l2 = (inv[1][0], inv[1][1], inv[1][2])
    grad_l3 = (inv[2][0], inv[2][1], inv[2][2])
    grad_l0 = scale(-1.0, add(add(grad_l1, grad_l2), grad_l3))
    return verts, (grad_l0, grad_l1, grad_l2, grad_l3), inv


def barycentric(tet_name, point):
    verts, _, inv = tetra_data(tet_name)
    rhs = sub(point, verts[0])
    l1, l2, l3 = mat_vec(inv, rhs)
    l0 = 1.0 - l1 - l2 - l3
    return (l0, l1, l2, l3)


def point_from_bary(tet_name, bary):
    verts = [VERTICES[i] for i in SUB_TETS[tet_name]]
    out = (0.0, 0.0, 0.0)
    for coeff, vertex in zip(bary, verts):
        out = add(out, scale(coeff, vertex))
    return out


def random_bary(rng, n, min_weight=0.08):
    vals = [rng.random() for _ in range(n)]
    total = sum(vals)
    vals = [v / total for v in vals]
    return tuple(min_weight + (1.0 - n * min_weight) * v for v in vals)


def random_tet_point(rng, tet_name, min_weight=0.08):
    return point_from_bary(tet_name, random_bary(rng, 4, min_weight))


def random_triangle_point(rng, tri):
    a = random_bary(rng, 3, 0.08)
    out = (0.0, 0.0, 0.0)
    for coeff, vertex_id in zip(a, tri):
        out = add(out, scale(coeff, VERTICES[vertex_id]))
    return out


def p2_shapes_and_grads(tet_name, point):
    lam = barycentric(tet_name, point)
    _, grad_lam, _ = tetra_data(tet_name)
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

    values = [lam[i] * (2.0 * lam[i] - 1.0) for i in range(4)]
    grads = [scale(4.0 * lam[i] - 1.0, grad_lam[i]) for i in range(4)]
    for i, j in pairs:
        values.append(4.0 * lam[i] * lam[j])
        grads.append(add(scale(4.0 * lam[i], grad_lam[j]), scale(4.0 * lam[j], grad_lam[i])))
    return values, grads


def local_p2_eval(tet_name, global_dof_values, point):
    values, grads = p2_shapes_and_grads(tet_name, point)
    local_ids = LOCAL_DOF_MAP[tet_name]
    value = 0.0
    grad = (0.0, 0.0, 0.0)
    for local_value, local_grad, global_id in zip(values, grads, local_ids):
        coeff = global_dof_values[global_id]
        value += coeff * local_value
        grad = add(grad, scale(coeff, local_grad))
    return value, grad


def c_basis_values(value_fns, point, vertex_ptrs):
    return [call_scalar(fn, point, vertex_ptrs) for fn in value_fns]


def c_basis_grads(grad_fns, point, vertex_ptrs):
    return [call_vector(fn, point, vertex_ptrs) for fn in grad_fns]


def c_expansion(value_fns, grad_fns, coeffs, point, vertex_ptrs):
    vals = c_basis_values(value_fns, point, vertex_ptrs)
    grads = c_basis_grads(grad_fns, point, vertex_ptrs)
    value = sum(c * v for c, v in zip(coeffs, vals))
    grad = (0.0, 0.0, 0.0)
    for coeff, basis_grad in zip(coeffs, grads):
        grad = add(grad, scale(coeff, basis_grad))
    return value, grad


def global_quadratic(coeffs, point):
    x, y, z = point
    a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz = coeffs
    value = (
        a0 + ax * x + ay * y + az * z
        + axx * x * x + ayy * y * y + azz * z * z
        + axy * x * y + axz * x * z + ayz * y * z
    )
    grad = (
        ax + 2.0 * axx * x + axy * y + axz * z,
        ay + 2.0 * ayy * y + axy * x + ayz * z,
        az + 2.0 * azz * z + axz * x + ayz * y,
    )
    return value, grad


def finite_difference_gradient(value_fns, coeffs, point, vertex_ptrs, h=1.0e-6):
    grad = []
    for axis in range(3):
        step = [0.0, 0.0, 0.0]
        step[axis] = h
        plus = add(point, tuple(step))
        minus = sub(point, tuple(step))
        v_plus, _ = c_expansion(value_fns, [], coeffs, plus, vertex_ptrs)
        v_minus, _ = c_expansion(value_fns, [], coeffs, minus, vertex_ptrs)
        grad.append((v_plus - v_minus) / (2.0 * h))
    return tuple(grad)


def max_error_record(current, err, detail):
    best_err, _ = current
    if err > best_err:
        return err, detail
    return current


def print_check(name, err, tol, detail):
    status = "PASS" if err <= tol else "FAIL"
    print(f"[{status}] {name}: max_err={err:.3e} tol={tol:.1e}")
    if status == "FAIL" and detail is not None:
        print(f"       worst={detail}")
    return status == "PASS"


def main():
    rng = random.Random(20260423)
    root = Path(__file__).resolve().parent
    lib = ctypes.CDLL(str(root / "twin_tetrahedron.2.bas_fun.so"))
    _, vertex_ptrs = make_vertex_ptrs()
    value_fns = [load_function(lib, name) for name in BASIS_NAMES]
    grad_fns = [load_function(lib, name) for name in GRADIENT_NAMES]

    all_passed = True

    # 1. Nodal interpolation matrix.
    worst_diag = (0.0, None)
    worst_off = (0.0, None)
    for i, point in enumerate(DOF_POINTS):
        vals = c_basis_values(value_fns, point, vertex_ptrs)
        for j, value in enumerate(vals):
            target = 1.0 if i == j else 0.0
            err = abs(value - target)
            detail = (i, j, BASIS_NAMES[j], point, value, target)
            if i == j:
                worst_diag = max_error_record(worst_diag, err, detail)
            else:
                worst_off = max_error_record(worst_off, err, detail)
    all_passed &= print_check("nodal diagonal", worst_diag[0], TOL_VALUE, worst_diag[1])
    all_passed &= print_check("nodal off-diagonal", worst_off[0], TOL_VALUE, worst_off[1])

    # 2. Partition of unity and gradient sum.
    worst_sum = (0.0, None)
    worst_grad_sum = (0.0, None)
    for tet_name in SUB_TETS:
        for _ in range(40):
            point = random_tet_point(rng, tet_name)
            vals = c_basis_values(value_fns, point, vertex_ptrs)
            grads = c_basis_grads(grad_fns, point, vertex_ptrs)
            err_sum = abs(sum(vals) - 1.0)
            grad_sum = (sum(g[0] for g in grads), sum(g[1] for g in grads), sum(g[2] for g in grads))
            worst_sum = max_error_record(worst_sum, err_sum, (tet_name, point, sum(vals)))
            worst_grad_sum = max_error_record(worst_grad_sum, max_abs_vec(grad_sum), (tet_name, point, grad_sum))
    all_passed &= print_check("partition of unity", worst_sum[0], TOL_VALUE, worst_sum[1])
    all_passed &= print_check("gradient sum is zero", worst_grad_sum[0], TOL_GRAD, worst_grad_sum[1])

    # 3. Support mask.
    worst_support = (0.0, None)
    for tet_name, zero_ids in SUPPORT_ZERO.items():
        for _ in range(40):
            point = random_tet_point(rng, tet_name)
            vals = c_basis_values(value_fns, point, vertex_ptrs)
            grads = c_basis_grads(grad_fns, point, vertex_ptrs)
            for idx in zero_ids:
                value_err = abs(vals[idx])
                grad_err = max_abs_vec(grads[idx])
                worst_support = max_error_record(
                    worst_support,
                    max(value_err, grad_err),
                    (tet_name, BASIS_NAMES[idx], point, vals[idx], grads[idx]),
                )
    all_passed &= print_check("support zero mask", worst_support[0], TOL_GRAD, worst_support[1])

    # 4. Global quadratic reproduction.
    worst_quad_value = (0.0, None)
    worst_quad_grad = (0.0, None)
    for case in range(25):
        coeffs = [rng.uniform(-2.0, 2.0) for _ in range(10)]
        dof_values = [global_quadratic(coeffs, p)[0] for p in DOF_POINTS]
        for tet_name in SUB_TETS:
            for _ in range(12):
                point = random_tet_point(rng, tet_name)
                value_c, grad_c = c_expansion(value_fns, grad_fns, dof_values, point, vertex_ptrs)
                value_ref, grad_ref = global_quadratic(coeffs, point)
                worst_quad_value = max_error_record(
                    worst_quad_value,
                    abs(value_c - value_ref),
                    (case, tet_name, point, value_c, value_ref),
                )
                worst_quad_grad = max_error_record(
                    worst_quad_grad,
                    max_abs_vec(sub(grad_c, grad_ref)),
                    (case, tet_name, point, grad_c, grad_ref),
                )
    all_passed &= print_check("global P2 value reproduction", worst_quad_value[0], TOL_VALUE, worst_quad_value[1])
    all_passed &= print_check("global P2 gradient reproduction", worst_quad_grad[0], TOL_GRAD, worst_quad_grad[1])

    # 5. Full continuous piecewise P2 patch test using ordinary tetra P2 on T+/T-.
    worst_patch_value = (0.0, None)
    worst_patch_grad = (0.0, None)
    for case in range(40):
        dof_values = [rng.uniform(-2.0, 2.0) for _ in DOF_POINTS]
        for tet_name in SUB_TETS:
            for _ in range(20):
                point = random_tet_point(rng, tet_name)
                value_c, grad_c = c_expansion(value_fns, grad_fns, dof_values, point, vertex_ptrs)
                value_ref, grad_ref = local_p2_eval(tet_name, dof_values, point)
                worst_patch_value = max_error_record(
                    worst_patch_value,
                    abs(value_c - value_ref),
                    (case, tet_name, point, value_c, value_ref),
                )
                worst_patch_grad = max_error_record(
                    worst_patch_grad,
                    max_abs_vec(sub(grad_c, grad_ref)),
                    (case, tet_name, point, grad_c, grad_ref),
                )
    all_passed &= print_check("piecewise P2 value patch", worst_patch_value[0], TOL_VALUE, worst_patch_value[1])
    all_passed &= print_check("piecewise P2 gradient patch", worst_patch_grad[0], TOL_GRAD, worst_patch_grad[1])

    # 6. Finite-difference self-consistency for C value/gradient functions.
    worst_fd = (0.0, None)
    for case in range(10):
        dof_values = [rng.uniform(-2.0, 2.0) for _ in DOF_POINTS]
        for tet_name in SUB_TETS:
            for _ in range(8):
                point = random_tet_point(rng, tet_name, min_weight=0.18)
                _, grad_c = c_expansion(value_fns, grad_fns, dof_values, point, vertex_ptrs)
                grad_fd = finite_difference_gradient(value_fns, dof_values, point, vertex_ptrs)
                worst_fd = max_error_record(
                    worst_fd,
                    max_abs_vec(sub(grad_c, grad_fd)),
                    (case, tet_name, point, grad_c, grad_fd),
                )
    all_passed &= print_check("C gradient vs finite-difference value", worst_fd[0], TOL_FD, worst_fd[1])

    # 7. One-sided checks near the internal face (v0,v2,v4).
    worst_interface_value = (0.0, None)
    worst_interface_grad = (0.0, None)
    for case in range(20):
        dof_values = [rng.uniform(-2.0, 2.0) for _ in DOF_POINTS]
        for tet_name, local_bary in (
            ("T+", (0.24, 1.0e-7, 0.38, 0.38)),
            ("T-", (0.24, 0.38, 1.0e-7, 0.38)),
        ):
            point = point_from_bary(tet_name, local_bary)
            value_c, grad_c = c_expansion(value_fns, grad_fns, dof_values, point, vertex_ptrs)
            value_ref, grad_ref = local_p2_eval(tet_name, dof_values, point)
            worst_interface_value = max_error_record(
                worst_interface_value,
                abs(value_c - value_ref),
                (case, tet_name, point, value_c, value_ref),
            )
            worst_interface_grad = max_error_record(
                worst_interface_grad,
                max_abs_vec(sub(grad_c, grad_ref)),
                (case, tet_name, point, grad_c, grad_ref),
            )
    all_passed &= print_check("near-interface one-sided value", worst_interface_value[0], TOL_VALUE, worst_interface_value[1])
    all_passed &= print_check("near-interface one-sided gradient", worst_interface_grad[0], TOL_GRAD, worst_interface_grad[1])

    # 8. Boundary trace support checks.
    worst_face = (0.0, None)
    for face_name, zero_ids in FACE_ZERO.items():
        for tri in FACE_SAMPLERS[face_name]:
            for _ in range(20):
                point = random_triangle_point(rng, tri)
                vals = c_basis_values(value_fns, point, vertex_ptrs)
                for idx in zero_ids:
                    worst_face = max_error_record(
                        worst_face,
                        abs(vals[idx]),
                        (face_name, tri, BASIS_NAMES[idx], point, vals[idx]),
                    )
    all_passed &= print_check("boundary face trace zero mask", worst_face[0], TOL_VALUE, worst_face[1])

    print()
    if all_passed:
        print("AUDIT PASS")
        return 0
    print("AUDIT FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
