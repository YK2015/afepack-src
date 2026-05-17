#!/usr/bin/env python3

"""Verify the nodal interpolation matrix of four_tetrahedron.2."""

import ctypes
from pathlib import Path


DoubleP = ctypes.POINTER(ctypes.c_double)


def load_function(lib, name):
    fn = getattr(lib, name)
    fn.argtypes = [DoubleP, ctypes.POINTER(DoubleP), ctypes.c_void_p]
    fn.restype = None
    return fn


def main():
    root = Path(__file__).resolve().parent
    lib = ctypes.CDLL(str(root / "four_tetrahedron.2.bas_fun.so"))

    vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.5, 0.5),
        (0.5, 0.0, 0.5),
        (0.5, 0.5, 0.0),
    ]
    vertex_arrays = [(ctypes.c_double * 3)(*v) for v in vertices]
    vertex_ptrs = (DoubleP * len(vertex_arrays))(
        *[ctypes.cast(item, DoubleP) for item in vertex_arrays]
    )

    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.5, 0.5),
        (0.5, 0.0, 0.5),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (0.0, 0.0, 0.5),
        (0.0, 0.75, 0.25),
        (0.0, 0.25, 0.75),
        (0.25, 0.0, 0.75),
        (0.75, 0.0, 0.25),
        (0.75, 0.25, 0.0),
        (0.25, 0.75, 0.0),
        (0.25, 0.5, 0.25),
        (0.25, 0.25, 0.5),
        (0.5, 0.25, 0.25),
        (0.0, 0.25, 0.25),
        (0.25, 0.0, 0.25),
        (0.25, 0.25, 0.0),
    ]

    names = (
        [f"phi_{i}" for i in range(1, 8)]
        + [f"psi_{i}" for i in range(1, 13)]
        + [f"chi_{i}" for i in range(1, 4)]
    )
    functions = [load_function(lib, name) for name in names]

    max_diag_error = 0.0
    max_off_error = 0.0
    worst = None

    for i, point in enumerate(points):
        point_array = (ctypes.c_double * 3)(*point)
        point_ptr = ctypes.cast(point_array, DoubleP)
        for j, (name, fn) in enumerate(zip(names, functions)):
            output = (ctypes.c_double * 1)(0.0)
            fn(point_ptr, vertex_ptrs, ctypes.byref(output))
            value = output[0]
            target = 1.0 if i == j else 0.0
            err = abs(value - target)
            if i == j:
                if err > max_diag_error:
                    max_diag_error = err
                    worst = ("diag", i, j, name, point, value)
            else:
                if err > max_off_error:
                    max_off_error = err
                    worst = ("off", i, j, name, point, value)

    print(f"max_diag_error={max_diag_error:.3e}")
    print(f"max_off_error={max_off_error:.3e}")
    print(f"worst={worst}")


if __name__ == "__main__":
    main()
