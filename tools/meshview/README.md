# meshview

Headless rasterizer for **comparing triangle meshes side by side**. numpy and
Pillow only — no GPU, no display server, no window system, no scipy. Built for
one specific job: look at the same patch of two or three meshes at the same
scale and orientation and see which one is smoother.

It lives here in `tools/` because it was written to debug molses' surface
quality. It has no dependency on molses and reads any OBJ, so it is intended to
become its own repository and crate later.

```sh
python3 tools/meshview/render.py ref.obj ours.obj -o cmp.png \
    --zoom 110 --focus 4.29,5.05,4.55 --aa 2 \
    --label reference --label ours

python3 tools/meshview/selftest.py     # 11 checks, no files needed
```

## Why the camera is the whole point

Every mesh in one invocation shares a single `Camera`, so `--focus` and
`--zoom` place the *same world-space region* at the *same scale* in each tile.
That is what makes a pixel difference mean a geometry difference. Framing each
mesh on its own centroid instead would silently change magnification between
panels, which is exactly how a comparison starts lying to you.

`--focus x,y,z` takes a world-space point (for molses meshes, coordinates in
Å) and puts it at the centre of every tile. With no `--focus`, the first
mesh's centroid is used for all of them.

## Normals, and why the output is honest about them

Shading interpolates the file's own per-vertex `vn` normals, so triangulation
noise in the **normal field** is directly visible — which is the artifact worth
looking for on a subdivided or projected surface. Files without `vn` get
area-weighted geometric normals at load.

Two deliberate choices, both of which exist because getting them wrong produced
a wrong conclusion once:

- **Nothing is back-face culled.** Visibility is resolved by a z-buffer. An
  earlier version culled on a geometric orientation test, which disagrees with
  the stored normals of a creased surface; it drew a checkerboard of
  background-coloured speckles onto every mesh, **including the reference that
  was being compared against**, and those speckles were briefly blamed on the
  meshes. Orientation belongs to the shading, not to visibility: the
  interpolated normal is rotated toward the viewer instead.
- **Coverage is tracked explicitly**, via `render(..., return_mask=True)`, and
  never inferred from pixel values. A surface facing the light reaches exactly
  255, which is also the background value — colour alone cannot distinguish
  "brightest possible pixel" from "no surface here", and a test that conflated
  them failed on a flawless render of a sphere.

## Self-tests

`selftest.py` builds a UV sphere of known radius and analytic normals and
asserts properties that must hold, rather than diffing stored images:

| check | catches |
|---|---|
| silhouette area vs `πr²` | projection, winding and scale errors (currently within 0.07%) |
| no uncovered pixel enclosed by covered ones | the checkerboard regression above |
| shading varies smoothly, no single-pixel step | normals discontinuous across faces |
| flipped normals cover identically | visibility depending on normal sign |
| byte-identical re-render | nondeterminism, e.g. float-order or set iteration |
| fixture loaded from a hand-written OBJ | the loader itself, incl. `v//vn` and polygon fan-triangulation |

Run it as the pre-commit gate before trusting any picture:

```sh
python3 tools/meshview/selftest.py
```

## Limits

Read these before drawing a conclusion from an image.

- **Orthographic only.** No perspective, and at the zooms used for surface
  inspection that is the right choice, but it is not how a molecule looks in a
  real viewer.
- **Lambert with a constant ambient.** No specular, shadows or ambient
  occlusion, so crevices never look as deep as they are.
- **Grayscale.** Colour carries no information here; brightness is the signal.
- **Performance.** The per-face loop is Python. Frames are first filtered to
  those touching the viewport, and sub-pixel faces take a flat-shading fast
  path, which is what most triangles of a fine mesh hit — but it is still about
  5 s per mesh at zoom 110 on a 1M-vertex mesh. This is the main reason to
  port it to Rust.
- **Not a measurement tool.** It shows smoothness; it does not quantify it.
  Numbers come from `molses`'s own `roughness` and `normal_dev` log fields, or
  from the analytic-residual checks described in the main README. A render can
  look right while the metrics disagree, and the metrics are the ones that
  decide.

## Comparing meshes at different resolutions

Do **not** compare shading statistics across zoom levels or grid spacings.
Per-edge normal angle scales roughly linearly with mesh spacing; angle per unit
edge length is the scale-free quantity, and even then geometric normal
roughness rises at fine spacings because more triangles are near-degenerate
slivers. Same resolution, same camera, same light — or state the difference
explicitly.
