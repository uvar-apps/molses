#!/usr/bin/env python3
"""Self-tests for meshview. Run: python3 tools/meshview/selftest.py

These exist because an early version of this renderer back-face-culled with a
geometric test that disagreed with the stored normals of the meshes under
inspection. It drew a checkerboard of background-coloured speckles on every
mesh — including the reference that was being compared against — and for a
while that was mistaken for a defect in the meshes. A closed surface must not
have unreachable background inside its silhouette, so test 2 below asserts
exactly that, and test 4 asserts it still holds when normals are inconsistent.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from render import Camera, Mesh, load_obj, render  # noqa: E402


def uv_sphere(radius: float = 5.0, nu: int = 96, nv: int = 48) -> Mesh:
    """Closed sphere with outward analytic normals and consistent winding."""
    pts, nrm = [], []
    for j in range(nv + 1):
        theta = math.pi * j / nv
        for i in range(nu):
            phi = 2 * math.pi * i / nu
            n = (
                math.sin(theta) * math.cos(phi),
                math.sin(theta) * math.sin(phi),
                math.cos(theta),
            )
            nrm.append(n)
            pts.append((radius * n[0], radius * n[1], radius * n[2]))
    idx = []
    for j in range(nv):
        for i in range(nu):
            a = j * nu + i
            b = j * nu + (i + 1) % nu
            c = (j + 1) * nu + i
            d = (j + 1) * nu + (i + 1) % nu
            idx.append((a, c, b))
            idx.append((b, c, d))
    return Mesh(np.asarray(pts, float), np.asarray(nrm, float), np.asarray(idx, np.int64))


def enclosed_uncovered(covered: np.ndarray) -> int:
    """Count uncovered pixels the outside cannot reach — the speckle signal.

    Takes the renderer's coverage mask, never the pixel values: a surface
    aligned with the light reaches 255, which is also the background value, so
    an image-based test would report bright pixels as holes. Flood over
    4-neighbours from every border pixel; written out instead of using PIL's
    floodfill, which does not behave on grayscale, so the test cannot silently
    pass.
    """
    img = np.where(covered, 0, 255)
    h, w = img.shape
    bg = img >= 250
    seen = np.zeros((h, w), bool)
    stack = [(0, y) for y in range(w)] + [(h - 1, y) for y in range(w)]
    stack += [(x, 0) for x in range(h)] + [(x, h - 1) for x in range(h)]
    dq = __import__("collections").deque()
    for x, y in stack:
        if bg[x, y] and not seen[x, y]:
            seen[x, y] = True
            dq.append((x, y))
    while dq:
        x, y = dq.popleft()
        for nx_, ny_ in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if 0 <= nx_ < h and 0 <= ny_ < w and bg[nx_, ny_] and not seen[nx_, ny_]:
                seen[nx_, ny_] = True
                dq.append((nx_, ny_))
    return int((bg & ~seen).sum())


def silhouette_area(covered: np.ndarray) -> int:
    return int(covered.sum())


fails: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not cond:
        fails.append(name)


def main() -> int:
    print("meshview selftest")
    sphere = uv_sphere()
    size, zoom = 400, 30.0
    cam = Camera(focus=np.zeros(3), zoom=zoom)
    img, mask = render(sphere, cam, size=size, return_mask=True)

    # 1. Silhouette matches the analytic disc: proves projection + winding.
    want = math.pi * (5.0 * zoom) ** 2
    got = silhouette_area(mask)
    check(
        "silhouette area matches analytic disc",
        abs(got - want) / want < 0.05,
        f"{got} vs {want:.0f}",
    )

    # 2. THE REGRESSION TEST: no isolated uncovered pixel inside a closed
    #    surface. This is the checkerboard failure mode.
    holes = enclosed_uncovered(mask)
    check("no enclosed uncovered pixels", holes == 0, f"{holes} speckles")

    # 3. Shading must vary smoothly across the sphere, not in steps.
    core = img[size // 2 - 60 : size // 2 + 60, size // 2 - 60 : size // 2 + 60]
    check(
        "interior shading varies smoothly",
        2.0 < float(core.std()) < 60.0,
        f"std {core.std():.1f}",
    )
    grad = np.abs(np.diff(core.astype(float), axis=0)).max()
    check("no single-pixel shading discontinuity", grad < 40.0, f"max row step {grad:.0f}")

    # 4. Same robustness when normals point the wrong way: coverage must not
    #    depend on normal orientation (real risk when comparing meshes whose
    #    normals came from different passes).
    bad_img, bad_mask = render(Mesh(sphere.positions, -sphere.normals, sphere.indices), cam, size=size, return_mask=True)
    check("flipped normals cover identically", np.array_equal(bad_mask, mask), "")
    check("no enclosed uncovered pixels (flipped)", enclosed_uncovered(bad_mask) == 0, "")

    # 5. Determinism: identical inputs, identical bytes.
    again = render(sphere, cam, size=size)
    check("render is deterministic", np.array_equal(again, img), "")

    # 6. Loader: explicit fixture with v//vn indices and a polygon face.
    fixture = "\n".join(
        [
            "# fixture",
            "v 0 0 0",
            "v 1 0 0",
            "v 1 1 0",
            "v 0 1 0",
            "v 0 1 1",
            "vn 0 0 1",
            "vn 0 -1 0",
            "vn 1 0 0",
            "f 1//1 2//1 3//1 4//1",
            "f 4//2 5//2 1//2",
            "",
        ]
    )
    with open("/tmp/meshview_probe.obj", "w") as fh:
        fh.write(fixture)
    m = load_obj("/tmp/meshview_probe.obj")
    check("loader reads 5 vertices", len(m.positions) == 5, str(len(m.positions)))
    check("loader fans a quad into triangles", len(m.indices) == 3, str(len(m.indices)))
    check("loader keeps one normal per vertex", len(m.normals) == 5, str(len(m.normals)))
    check(
        "loader normals are unit length",
        float(np.abs(np.linalg.norm(m.normals, axis=1) - 1).max()) < 1e-9,
    )

    print("FAILURES:", len(fails))
    return 1 if fails else 0


def sil_mean(img: np.ndarray) -> float:
    m = img < 250
    return float(img[m].mean() / 255.0) if m.any() else 0.0


if __name__ == "__main__":
    sys.exit(main())
