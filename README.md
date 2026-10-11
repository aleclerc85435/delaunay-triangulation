# delaunay

Compute the Delaunay triangulation of a set of 2D points using the
Bowyer–Watson incremental insertion algorithm.

## Usage

```js
import { delaunay } from './src/core.js';

const tris = delaunay([
  { x: 0, y: 0 },
  { x: 4, y: 0 },
  { x: 0, y: 4 },
  { x: 1, y: 1 },
]);
// tris === [[0, 2, 3], [0, 1, 3], [1, 2, 3]]  (order may vary)
```

`delaunay(points)` takes an array of `{ x, y }` and returns an array of
`[i, j, k]` index triples into the input. Points must have numeric `x`/`y`.
With fewer than three points the result is `[]`.

## Why

This exists for callers that need a correct, dependency-free Delaunay
triangulation in pure ESM JavaScript — meshes for terrain, proximity graphs,
interpolation scaffolding — without pulling in a full computational-geometry
package. The trade-off is speed: Bowyer–Watson is O(n²) worst case. For inputs
in the hundreds-to-low-thousands range the constant factor and simplicity win;
for hundreds of thousands of points you want a sweep-line implementation
instead.

## Edges you will hit

- **Collinear input** has no Delaunay triangulation (a circumcircle through
  three collinear points has infinite radius). The library returns `[]`
  rather than fabricating zero-area triangles.
- **Duplicate points** are tolerated: no triangle will join two coincident
  vertices, so the mesh stays valid but the duplicate contributes nothing.
- **Co-circular points** (four or more on a common circle) have multiple
  valid Delaunay triangulations; this library returns one of them. Do not
  assert on which diagonal is chosen.
- **Ties** in the in-circumcircle test are broken with a small EPS so that
  near-coincident points terminate instead of looping.

## Exports

- `delaunay(points)` — the function described above.
- `triangulate` — default export, alias of `delaunay`.

## Design notes

The window stores values eagerly rather than keeping running aggregates. Running
sums drift with floating point over long streams, and recomputing from a small
buffer is cheap enough that the drift is not worth the speed.

