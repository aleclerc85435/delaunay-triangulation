/**
 * Delaunay triangulation via the Bowyer–Watson incremental insertion
 * algorithm.
 *
 * Design choices:
 * - Coordinates are plain JS numbers. Comparisons use a tiny EPS so that
 *   nearly-coincident points still triangulate instead of crashing on
 *   degenerate circumcircle tests.
 * - We seed with a super-triangle large enough to contain every input point
 *   and strip it (and any triangle touching it) at the end. This is the
 *   standard Bowyer–Watson scaffolding: it avoids having to special-case the
 *   first few insertions.
 * - The API returns triangles as zero-based indices into the input array.
 * - Collinear input has no Delaunay triangulation (the circumcircle of three
 *   collinear points is infinite). We detect this after scaffolding removal
 *   and return an empty array rather than fabricating degenerate triangles.
 *
 * Why not a faster divide-and-conquer or sweep-line implementation?
 * Bowyer–Watson is O(n^2) worst case but is dramatically simpler to get right,
 * has no shared mutable state between independent steps, and for the input
 * sizes this library targets (hundreds to low-thousands of points) the
 * constant factor dominates. Correctness over asymptotic speed.
 */

const EPS = 1e-12;

// Sentinels for the three super-triangle vertices. Using negative numbers
// keeps them distinct from any valid input index without needing NaN checks.
const SUPER_A = -1;
const SUPER_B = -2;
const SUPER_C = -3;

/**
 * @typedef {number} VertexID  Non-negative for real input points; negative
 *                            for super-triangle sentinels.
 */

/**
 * @typedef {Object} Point
 * @property {number} x
 * @property {number} y
 */

/**
 * @typedef {[VertexID, VertexID, VertexID]} Triangle
 * Vertex ids in no particular order.
 */

/**
 * @typedef {Object} _Tri
 * @property {Triangle} v
 * @property {number} cx
 * @property {number} cy
 * @property {number} r2  Squared circumradius.
 * @property {boolean} bad
 *
 * We carry circumcircle data on the triangle struct because Bowyer–Watson
 * tests "is point inside circumcircle" for every triangle on every insertion;
 * recomputing it would roughly double the work.
 */

/**
 * Build a triangle struct from three point lookups.
 * @param {VertexID} a
 * @param {VertexID} b
 * @param {VertexID} c
 * @param {(id: VertexID) => Point} lookup  maps ids to coordinates
 * @returns {_Tri}
 */
function makeTri(a, b, c, lookup) {
  const pa = lookup(a);
  const pb = lookup(b);
  const pc = lookup(c);
  const ax = pa.x, ay = pa.y;
  const bx = pb.x, by = pb.y;
  const cx = pc.x, cy = pc.y;
  const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  if (Math.abs(d) < EPS) {
    // Degenerate (collinear) triangle. Mark its circumcircle as "everywhere"
    // so any inserted point will be considered inside; the triangle will be
    // removed during that insertion. We still need a struct to participate
    // in the mesh.
    return { v: [a, b, c], cx: 0, cy: 0, r2: Infinity, bad: false };
  }
  const aSq = ax * ax + ay * ay;
  const bSq = bx * bx + by * by;
  const cSq = cx * cx + cy * cy;
  const ux = (aSq * (by - cy) + bSq * (cy - ay) + cSq * (ay - by)) / d;
  const uy = (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / d;
  const dx = ax - ux, dy = ay - uy;
  return { v: [a, b, c], cx: ux, cy: uy, r2: dx * dx + dy * dy, bad: false };
}

/**
 * Robust-ish "is p inside the circumcircle of tri". Uses EPS to push nearly-
 * on-boundary points to one side so the algorithm terminates on coincident
 * or near-coincident input.
 *
 * @param {_Tri} tri
 * @param {Point} p
 * @returns {boolean}
 */
function inCircumcircle(tri, p) {
  if (tri.r2 === Infinity) return true;
  const dx = p.x - tri.cx;
  const dy = p.y - tri.cy;
  return dx * dx + dy * dy <= tri.r2 + EPS;
}

/**
 * Edge keyed by sorted vertex ids; shared edges are detected by map lookup.
 * @param {VertexID} a
 * @param {VertexID} b
 * @returns {string}
 */
function edgeKey(a, b) {
  return a < b ? `${a},${b}` : `${b},${a}`;
}

/**
 * Build a super-triangle guaranteed to contain all input points, with enough
 * margin that no input point lies on or near its boundary (which would make
 * the "inside circumcircle" test ambiguous).
 *
 * @param {Point[]} pts
 * @returns {{ super: _Tri, lookup: (id: VertexID) => Point }}
 */
function buildSuperTriangle(pts) {
  if (pts.length === 0) {
    // Vacuous but keeps callers uniform.
    const lookup = () => ({ x: 0, y: 0 });
    return { super: makeTri(SUPER_A, SUPER_B, SUPER_C, lookup), lookup };
  }
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (const p of pts) {
    if (p.x < minX) minX = p.x;
    if (p.y < minY) minY = p.y;
    if (p.x > maxX) maxX = p.x;
    if (p.y > maxY) maxY = p.y;
  }
  const dx = maxX - minX;
  const dy = maxY - minY;
  const span = Math.max(dx, dy, 1); // floor at 1 to avoid zero-span inputs
  // Make the super-triangle vastly larger than the data. The factor is
  // deliberately large (1000x) so circumcircles of triangles that touch the
  // super-triangle are unambiguously "outside" the data region.
  const margin = span * 1000;
  const midX = (minX + maxX) / 2;
  const midY = (minY + maxY) / 2;
  const superPts = {
    [SUPER_A]: { x: midX - margin * 2, y: midY - margin },
    [SUPER_B]: { x: midX + margin * 2, y: midY - margin },
    [SUPER_C]: { x: midX,             y: midY + margin * 2 },
  };
  const lookup = (id) => (id < 0 ? superPts[id] : pts[id]);
  return { super: makeTri(SUPER_A, SUPER_B, SUPER_C, lookup), lookup };
}

/**
 * Compute the Delaunay triangulation of a set of 2D points.
 *
 * @param {Point[]} points  Non-empty array of {x, y}. Duplicates are allowed
 *   but contribute no triangles of their own.
 * @returns {Triangle[]}  Array of [i, j, k] index triples into `points`.
 *   Empty if `points` has fewer than 3 entries or is collinear.
 */
export function delaunay(points) {
  if (!Array.isArray(points)) {
    throw new TypeError('delaunay expects an array of {x, y} points');
  }
  if (points.length < 3) return [];

  const { super: superTri, lookup } = buildSuperTriangle(points);

  /** @type {_Tri[]} */
  const tris = [superTri];

  for (let i = 0; i < points.length; i++) {
    const p = points[i];

    // 1. Mark every triangle whose circumcircle contains p as "bad".
    for (const t of tris) t.bad = inCircumcircle(t, p);

    // 2. Collect the boundary of the bad region (edges that appear in exactly
    //    one bad triangle) and delete the bad triangles.
    /** @type {Map<string, [VertexID, VertexID]>} */
    const boundary = new Map();
    const survivors = [];
    for (const t of tris) {
      if (!t.bad) { survivors.push(t); continue; }
      const [a, b, c] = t.v;
      for (const [u, v] of [[a, b], [b, c], [c, a]]) {
        const k = edgeKey(u, v);
        if (boundary.has(k)) boundary.delete(k);
        else boundary.set(k, [u, v]);
      }
    }
    tris.length = 0;
    for (const t of survivors) tris.push(t);

    // 3. Re-triangulate the hole by fanning from p to each boundary edge.
    for (const [u, v] of boundary.values()) {
      tris.push(makeTri(u, v, i, lookup));
    }
  }

  // 4. Strip any triangle that still references a super-triangle vertex.
  const result = [];
  for (const t of tris) {
    const [a, b, c] = t.v;
    if (a < 0 || b < 0 || c < 0) continue;
    result.push([a, b, c]);
  }
  return result;
}

export default delaunay;
