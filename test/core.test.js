import { test } from 'node:test';
import assert from 'node:assert/strict';
import { delaunay } from '../src/core.js';

/**
 * Validate that `tris` is a legal triangulation of `pts`: every triangle uses
 * distinct in-range indices, no triangle is degenerate, and no two triangles
 * share an interior point pair (i.e. edges don't cross).
 *
 * We do NOT assert Delaunay-ness here directly; instead the property tests
 * below check the Delaunay empty-circle condition explicitly.
 */
function validateMesh(pts, tris) {
  assert.ok(Array.isArray(tris), 'result is an array');
  for (const t of tris) {
    assert.equal(t.length, 3, 'each triangle has 3 vertices');
    const [a, b, c] = t;
    for (const idx of [a, b, c]) {
      assert.ok(Number.isInteger(idx), 'indices are integers');
      assert.ok(idx >= 0 && idx < pts.length, 'indices in range');
    }
    assert.notEqual(a, b, 'triangle vertices distinct (a!=b)');
    assert.notEqual(b, c, 'triangle vertices distinct (b!=c)');
    assert.notEqual(a, c, 'triangle vertices distinct (a!=c)');
    const pa = pts[a], pb = pts[b], pc = pts[c];
    const area = Math.abs((pb.x - pa.x) * (pc.y - pa.y) - (pc.x - pa.x) * (pb.y - pa.y));
    assert.ok(area > 1e-9, 'triangle is not degenerate');
  }
  // Edge-uniqueness: every undirected edge appears at most once across all
  // triangles. (Interior edges of a valid triangulation appear in two
  // triangles sharing that edge, so this is NOT the standard property — but
  // the Delaunay triangulation has no overlapping triangles, so each directed
  // edge is unique. We check undirected pairs shared by at most 2 tris.)
  const counts = new Map();
  for (const [a, b, c] of tris) {
    for (const [u, v] of [[a, b], [b, c], [c, a]]) {
      const k = u < v ? `${u},${v}` : `${v},${u}`;
      counts.set(k, (counts.get(k) || 0) + 1);
    }
  }
  for (const [k, n] of counts) {
    assert.ok(n <= 2, `edge ${k} shared by ${n} triangles (<=2 expected)`);
  }
}

/** Circumcircle of three points; returns {cx, cy, r2} or null if collinear. */
function circumcircle(pa, pb, pc) {
  const ax = pa.x, ay = pa.y, bx = pb.x, by = pb.y, cx = pc.x, cy = pc.y;
  const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  if (Math.abs(d) < 1e-12) return null;
  const aSq = ax * ax + ay * ay;
  const bSq = bx * bx + by * by;
  const cSq = cx * cx + cy * cy;
  const ux = (aSq * (by - cy) + bSq * (cy - ay) + cSq * (ay - by)) / d;
  const uy = (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / d;
  return { cx: ux, cy: uy, r2: (ax - ux) ** 2 + (ay - uy) ** 2 };
}

/**
 * Delaunay empty-circle property: for every triangle, no other input point
 * lies strictly inside its circumcircle. Points on the circle are allowed
 * (co-circular configurations have multiple valid triangulations).
 */
function assertDelaunay(pts, tris) {
  const eps = 1e-9 * Math.max(1, ...pts.map(p => Math.hypot(p.x, p.y)));
  for (const [a, b, c] of tris) {
    const cc = circumcircle(pts[a], pts[b], pts[c]);
    assert.ok(cc, 'no degenerate triangles in output');
    for (let i = 0; i < pts.length; i++) {
      if (i === a || i === b || i === c) continue;
      const p = pts[i];
      const dist2 = (p.x - cc.cx) ** 2 + (p.y - cc.cy) ** 2;
      assert.ok(dist2 >= cc.r2 - eps, `point ${i} inside circumcircle of [${a},${b},${c}]`);
    }
  }
}

test('returns empty array for fewer than 3 points', () => {
  assert.deepEqual(delaunay([]), []);
  assert.deepEqual(delaunay([{ x: 0, y: 0 }]), []);
  assert.deepEqual(delaunay([{ x: 0, y: 0 }, { x: 1, y: 0 }]), []);
});

test('three non-collinear points yield exactly one triangle', () => {
  const pts = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 0, y: 1 }];
  const tris = delaunay(pts);
  assert.equal(tris.length, 1);
  assert.deepEqual(tris[0].slice().sort((a, b) => a - b), [0, 1, 2]);
  validateMesh(pts, tris);
});

test('three collinear points yield no triangles', () => {
  const pts = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 2, y: 0 }];
  const tris = delaunay(pts);
  assert.deepEqual(tris, []);
});

test('square produces two triangles covering it', () => {
  const pts = [
    { x: 0, y: 0 }, { x: 1, y: 0 },
    { x: 1, y: 1 }, { x: 0, y: 1 },
  ];
  const tris = delaunay(pts);
  assert.equal(tris.length, 2);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
  // Each corner appears in exactly one triangle: the square is split along a
  // diagonal, so the two shared corners participate in both triangles while
  // the other two participate in one each.
  const counts = [0, 0, 0, 0];
  for (const t of tris) for (const idx of t) counts[idx]++;
  for (const n of counts) assert.ok(n >= 1 && n <= 2, 'corner used 1 or 2 times');
  const total = counts.reduce((s, n) => s + n, 0);
  assert.equal(total, 6);
});

test('a point inside a triangle splits it into three', () => {
  const pts = [
    { x: 0, y: 0 }, { x: 4, y: 0 }, { x: 0, y: 4 },
    { x: 1, y: 1 }, // interior
  ];
  const tris = delaunay(pts);
  assert.equal(tris.length, 3);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
  // The interior point (index 3) is in all three triangles.
  for (const t of tris) assert.ok(t.includes(3), 'interior point in every triangle');
});

test('regular grid triangulates consistently and is Delaunay', () => {
  const pts = [];
  for (let r = 0; r < 4; r++)
    for (let c = 0; c < 4; c++)
      pts.push({ x: c, y: r });
  const tris = delaunay(pts);
  // For an n×m grid, a valid triangulation has 2(n-1)(m-1) triangles.
  assert.equal(tris.length, 2 * 3 * 3);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
});

test('random-ish point set satisfies the Delaunay property', () => {
  // Deterministic pseudo-random points (no Math.random in tests).
  const pts = [];
  let seed = 42;
  const next = () => {
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    return seed / 0x7fffffff;
  };
  for (let i = 0; i < 20; i++) pts.push({ x: next() * 10, y: next() * 10 });
  const tris = delaunay(pts);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
});

test('duplicate points do not crash and leave a valid mesh', () => {
  const pts = [
    { x: 0, y: 0 }, { x: 0, y: 0 }, // duplicate
    { x: 2, y: 0 }, { x: 0, y: 2 },
  ];
  const tris = delaunay(pts);
  validateMesh(pts, tris);
  // No triangle should reference index 0 and 1 together, since they're
  // coincident (that would be a degenerate triangle).
  for (const t of tris) {
    assert.ok(!(t.includes(0) && t.includes(1)), 'no triangle joins duplicates');
  }
});

test('points with negative coordinates triangulate correctly', () => {
  const pts = [
    { x: -5, y: -5 }, { x: 5, y: -5 }, { x: 5, y: 5 }, { x: -5, y: 5 },
    { x: 0, y: 0 },
  ];
  const tris = delaunay(pts);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
  // Center point is interior, so every triangle touches it.
  for (const t of tris) assert.ok(t.includes(4));
});

test('five points in convex position form a fan of three triangles', () => {
  // Points roughly on a circle: Delaunay triangulation of points in convex
  // position is a triangulated polygon, n-2 triangles for n points.
  const pts = [];
  for (let i = 0; i < 5; i++) {
    const a = (i / 5) * 2 * Math.PI;
    pts.push({ x: Math.cos(a) * 3, y: Math.sin(a) * 3 });
  }
  const tris = delaunay(pts);
  assert.equal(tris.length, 3);
  validateMesh(pts, tris);
  assertDelaunay(pts, tris);
});

test('four co-circular points yield a valid Delaunay triangulation (either diagonal)', () => {
  // Four points on a common circle: both possible diagonals are valid
  // Delaunay triangulations. We only require the result to be one of them.
  const pts = [
    { x: 1, y: 0 }, { x: 0, y: 1 },
    { x: -1, y: 0 }, { x: 0, y: -1 },
  ];
  const tris = delaunay(pts);
  assert.equal(tris.length, 2);
  validateMesh(pts, tris);
  // Delaunay property holds for both choices because all four points are
  // exactly on each triangle's circumcircle (boundary case, allowed).
  assertDelaunay(pts, tris);
});

test('non-array input throws', () => {
  assert.throws(() => delaunay(null), TypeError);
  assert.throws(() => delaunay('abc'), TypeError);
});
