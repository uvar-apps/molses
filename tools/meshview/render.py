#!/usr/bin/env python3
"""meshview — headless rasterizer for comparing triangle meshes side by side.

No GPU, no display server, no scipy: numpy and Pillow only. Designed for the
one task that matters when tuning a surface generator — putting the *same*
world-space region of several meshes at the *same* scale and orientation in a
single image, so a difference in the mesh is a difference in the picture and
nothing else.

    python3 tools/meshview/render.py --help
    python3 tools/meshview/render.py ref.obj ours.obj -o cmp.png --zoom 110 \
        --focus 4.29,5.05,4.55 --aa 2

Vertices may carry per-vertex normals (`vn`); those are interpolated across
each triangle and used for shading, which is what makes triangulation noise in
the normal field visible. Files without `vn` get area-weighted geometric
normals computed on load.

Known limits are documented in README.md: orthographic only, no specular or
shadows, and a Python per-face loop that is fine up to a few hundred thousand
visible triangles.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass

import numpy as np
from PIL import Image

__all__ = ["Mesh", "load_obj", "Camera", "render", "compare"]

FLOAT = np.float64


# --------------------------------------------------------------------------- #
# geometry input
# --------------------------------------------------------------------------- #


@dataclass
class Mesh:
    """Triangle soup: positions, outward unit normals, 0-based indices."""

    positions: np.ndarray  # (n, 3) float64
    normals: np.ndarray  # (n, 3) float64, unit length
    indices: np.ndarray  # (m, 3) int64

    @property
    def centroid(self) -> np.ndarray:
        return self.positions.mean(axis=0)


def load_obj(path: str) -> Mesh:
    """Read positions, `vn` normals and triangular faces from a Wavefront OBJ.

    Vertex normals are matched through the `v//vn` indices when present and
    fall back to the vertex order otherwise. Normals are re-normalized because
    exported files are not always unit length.
    """
    verts: list[tuple[float, float, float]] = []
    vns: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("v "):
                p = line.split()[1:4]
                verts.append((float(p[0]), float(p[1]), float(p[2])))
            elif line.startswith("vn "):
                p = line.split()[1:4]
                vns.append((float(p[0]), float(p[1]), float(p[2])))
            elif line.startswith("f "):
                idx = []
                for tok in line.split()[1:]:
                    parts = tok.split("/")
                    idx.append(int(parts[0]) - 1)
                # Fan-triangulate polygons so quad sources work too.
                for k in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[k], idx[k + 1]))

    if not verts:
        raise ValueError(f"{path}: no vertices")
    positions = np.asarray(verts, dtype=FLOAT)
    indices = np.asarray(faces, dtype=np.int64).reshape(-1, 3)
    if vns and len(vns) == len(verts):
        normals = np.asarray(vns, dtype=FLOAT)
    else:
        normals = _geometric_normals(positions, indices)
    norm = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / np.maximum(norm, 1e-12)
    return Mesh(positions, normals, indices)


def _geometric_normals(positions: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """Area-weighted vertex normals, used when a file ships without `vn`."""
    out = np.zeros_like(positions)
    a, b, c = positions[indices[:, 0]], positions[indices[:, 1]], positions[indices[:, 2]]
    face = np.cross(b - a, c - a)
    for k in range(3):
        np.add.at(out, indices[:, k], face)
    return out


# --------------------------------------------------------------------------- #
# camera
# --------------------------------------------------------------------------- #


def _rot_x(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=FLOAT)


def _rot_y(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=FLOAT)


def _rot_z(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=FLOAT)


@dataclass
class Camera:
    """Orthographic camera looking down -z, centred on `focus`.

    `focus` is a world-space point mapped to the image centre, which is what
    keeps two different meshes comparable: pass the same focus and zoom to both
    and any pixel difference is geometry, not framing.
    """

    focus: np.ndarray
    zoom: float = 100.0  # pixels per angstrom
    yaw: float = 0.0
    pitch: float = 0.0
    roll: float = 0.0

    @property
    def matrix(self) -> np.ndarray:
        return _rot_z(self.roll) @ _rot_y(self.yaw) @ _rot_x(self.pitch)

    def project(self, points: np.ndarray, size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return screen x, screen y and depth (larger = nearer the camera)."""
        v = points @ self.matrix.T
        f = self.matrix @ self.focus
        x = (v[:, 0] - f[0]) * self.zoom + size / 2.0
        y = -(v[:, 1] - f[1]) * self.zoom + size / 2.0
        return x, y, v[:, 2]


# --------------------------------------------------------------------------- #
# rasterizer
# --------------------------------------------------------------------------- #


def render(
    mesh: Mesh,
    cam: Camera,
    size: int = 900,
    light: tuple[float, float, float] = (0.35, 0.45, 0.82),
    ambient: float = 0.30,
    background: float = 255.0,
    return_mask: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Shade `mesh` into a grayscale `size`x`size` image.

    Visibility is resolved by a z-buffer, and no face is culled on orientation.
    That combination is deliberate: culling on a *geometric* test removes
    edge-on triangles and leaves gaps with a checkerboard pattern behind them,
    while the normals being displayed point the other way back there. Faces are
    instead kept and their interpolated normal is flipped toward the viewer, so
    a crease where stored and geometric normals disagree shows as itself rather
    than as a hole that gets mistaken for a meshing artifact.

    Returns `(size, size)` uint8, or `(image, coverage)` when `return_mask` is
    set. Coverage is tracked explicitly rather than inferred from the pixel
    value: a surface facing the light reaches exactly 255, which is the
    background value, so colour alone cannot tell "brightest pixel" from "no
    surface here". Every measurement in selftest.py goes through the mask.
    """
    sx, sy, depth = cam.project(mesh.positions, size)
    nview = mesh.normals @ cam.matrix.T
    ldir = np.asarray(light, dtype=FLOAT)
    ldir = ldir / max(np.linalg.norm(ldir), 1e-12)

    img = np.full((size, size), float(background))
    zbuf = np.full((size, size), -np.inf)
    covered = np.zeros((size, size), bool)

    faces = mesh.indices
    fx0, fx1 = sx[faces[:, 0]], sx[faces[:, 1]]
    fx2 = sx[faces[:, 2]]
    fy0, fy1, fy2 = sy[faces[:, 0]], sy[faces[:, 1]], sy[faces[:, 2]]
    # Cheap window reject: only faces touching the viewport are rasterized.
    inside = (
        (np.maximum(np.maximum(fx0, fx1), fx2) >= 0)
        & (np.minimum(np.minimum(fx0, fx1), fx2) < size)
        & (np.maximum(np.maximum(fy0, fy1), fy2) >= 0)
        & (np.minimum(np.minimum(fy0, fy1), fy2) < size)
    )
    order = np.nonzero(inside)[0]
    order = order[np.argsort((depth[faces[order]]).mean(axis=1))]

    for t in order:
        i0, i1, i2 = (int(v) for v in faces[t])
        ax, ay = sx[i1] - sx[i0], sy[i1] - sy[i0]
        bx, by = sx[i2] - sx[i0], sy[i2] - sy[i0]
        den = ax * by - ay * bx
        if abs(den) < 1e-9:
            continue
        x0 = max(int(math.floor(min(sx[i0], sx[i1], sx[i2]))), 0)
        x1 = min(int(math.ceil(max(sx[i0], sx[i1], sx[i2]))), size - 1)
        y0 = max(int(math.floor(min(sy[i0], sy[i1], sy[i2]))), 0)
        y1 = min(int(math.ceil(max(sy[i0], sy[i1], sy[i2]))), size - 1)
        if x1 < x0 or y1 < y0:
            continue
        px = (np.arange(x0, x1 + 1) + 0.5)[None, :] - sx[i0]
        py = (np.arange(y0, y1 + 1) + 0.5)[:, None] - sy[i0]
        w1 = (by * px - bx * py) / den
        w2 = (-ay * px + ax * py) / den
        w0 = 1.0 - w1 - w2
        hit = (w0 >= -1e-7) & (w1 >= -1e-7) & (w2 >= -1e-7)
        if not hit.any():
            continue
        z = w1 * depth[i1] + w2 * depth[i2] + w0 * depth[i0]
        win = z > zbuf[y0 : y1 + 1, x0 : x1 + 1]
        win &= hit
        if not win.any():
            continue
        if int(hit.sum()) == 1:
            # Sub-pixel face: flat-shade from the averaged normal and skip the
            # interpolation math. Most of the triangles in a fine mesh land
            # here, so this is where the frame rate comes from.
            n = (nview[i0] + nview[i1] + nview[i2]) / 3.0
            shade = np.full(w0.shape, ambient + (1 - ambient) * _lambert(n, ldir))
        else:
            nx = w1 * nview[i1, 0] + w2 * nview[i2, 0] + w0 * nview[i0, 0]
            ny = w1 * nview[i1, 1] + w2 * nview[i2, 1] + w0 * nview[i0, 1]
            nz = w1 * nview[i1, 2] + w2 * nview[i2, 2] + w0 * nview[i0, 2]
            nl = np.sqrt(nx * nx + ny * ny + nz * nz)
            inv = 1.0 / np.maximum(nl, 1e-12)
            nx, ny, nz = nx * inv, ny * inv, nz * inv
            # Two-sided lighting: turn the normal toward the camera before
            # shading, so an inverted normal reads as dark instead of leaving
            # the pixel at background colour.
            sgn = np.where(nz < 0, -1.0, 1.0)
            d = np.maximum((nx * ldir[0] + ny * ldir[1] + nz * ldir[2]) * sgn, 0.0)
            shade = ambient + (1 - ambient) * d
        block = img[y0 : y1 + 1, x0 : x1 + 1]
        block[win] = 255.0 * shade[win]
        img[y0 : y1 + 1, x0 : x1 + 1] = block
        cmask = covered[y0 : y1 + 1, x0 : x1 + 1]
        cmask[win] = True
        covered[y0 : y1 + 1, x0 : x1 + 1] = cmask
        zb = zbuf[y0 : y1 + 1, x0 : x1 + 1]
        zb[win] = z[win]
        zbuf[y0 : y1 + 1, x0 : x1 + 1] = zb
    out = img.astype(np.uint8)
    return (out, covered) if return_mask else out


def _lambert(n: np.ndarray, ldir: np.ndarray) -> float:
    nl = float(np.dot(n, ldir))
    if n[2] < 0:
        nl = float(np.dot(-n, ldir))
    return max(nl, 0.0)


# --------------------------------------------------------------------------- #
# output
# --------------------------------------------------------------------------- #


def _supersample(draw_one, size: int, aa: int) -> np.ndarray:
    """Render at aa*scale and box-filter down for antialiased edges."""
    if aa <= 1:
        return draw_one(size)
    big = draw_one(size * aa)
    return np.asarray(Image.fromarray(big).resize((size, size), Image.LANCZOS))


def compare(
    paths: list[str],
    out: str,
    focus: tuple[float, float, float] | None = None,
    zoom: float = 100.0,
    size: int = 620,
    aa: int = 2,
    per_row: int = 0,
    labels: list[str] | None = None,
    yaw: float = 0.0,
    pitch: float = 0.0,
    ambient: float = 0.30,
) -> str:
    """Render each mesh with an identical camera into one labelled contact sheet."""
    meshes = [load_obj(p) for p in paths]
    if focus is None:
        focus_t = meshes[0].centroid
    else:
        focus_t = np.asarray(focus, dtype=FLOAT)
    cam = Camera(focus=focus_t, zoom=zoom, yaw=yaw, pitch=pitch)
    images = [
        _supersample(lambda s, m=m: render(m, cam, size=s, ambient=ambient), size, aa)
        for m in meshes
    ]
    names = labels if labels else [p.rsplit("/", 1)[-1] for p in paths]
    cols = per_row or len(paths)
    rows = int(math.ceil(len(images) / cols))
    pad, head = 12, 26
    sheet = Image.new("RGB", (cols * size + (cols + 1) * pad, rows * (size + head) + pad), (248, 248, 248))
    from PIL import ImageDraw

    dr = ImageDraw.Draw(sheet)
    for k, im in enumerate(images):
        r, c = divmod(k, cols)
        x = pad + c * (size + pad)
        y = head + pad + r * (size + head)
        sheet.paste(Image.fromarray(im), (x, y))
        dr.text((x + 4, y - 18), names[k], fill=(200, 0, 0))
    sheet.save(out)
    return out


def _parse_focus(text: str) -> tuple[float, float, float]:
    parts = [float(v) for v in text.replace(" ", "").split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("--focus expects x,y,z")
    return (parts[0], parts[1], parts[2])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Render meshes for visual comparison.")
    ap.add_argument("meshes", nargs="+", help="OBJ files (positions, optional vn, faces)")
    ap.add_argument("-o", "--out", default="render.png")
    ap.add_argument("--zoom", type=float, default=100.0, help="pixels per angstrom")
    ap.add_argument("--focus", type=_parse_focus, default=None, help="world point at image centre")
    ap.add_argument("--size", type=int, default=620, help="tile edge length in pixels")
    ap.add_argument("--aa", type=int, default=2, help="supersampling factor (0/1 = off)")
    ap.add_argument("--per-row", type=int, default=0, help="grid columns (0 = one row)")
    ap.add_argument("--label", action="append", default=None, help="repeatable, one per mesh")
    ap.add_argument("--yaw", type=float, default=0.0, help="radians about z")
    ap.add_argument("--pitch", type=float, default=0.0, help="radians about x")
    ap.add_argument("--ambient", type=float, default=0.30)
    args = ap.parse_args(argv)
    out = compare(
        args.meshes,
        args.out,
        focus=args.focus,
        zoom=args.zoom,
        size=args.size,
        aa=max(args.aa, 1),
        per_row=args.per_row,
        labels=args.label,
        yaw=args.yaw,
        pitch=args.pitch,
        ambient=args.ambient,
    )
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
