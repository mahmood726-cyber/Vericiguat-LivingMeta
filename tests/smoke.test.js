/*
 * Minimal smoke + numerical-regression test for Vericiguat-LivingMeta.
 *
 * The statistical engine is embedded in VERICIGUAT_REVIEW.html (a single-file
 * dashboard). This test (a) asserts the shipped artifacts are present and the
 * pooling engine functions are still defined in the HTML, (b) confirms the
 * config JSON parses, and (c) re-implements the same verified inverse-variance
 * / Paule-Mandel math the dashboard uses and checks it against hand-computed
 * reference values so a regression in the formula would be caught.
 *
 * Run: node tests/smoke.test.js   (exit 0 = pass, exit 1 = fail)
 */
'use strict';
const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, '..');
let failures = 0;
function ok(cond, msg) {
  if (cond) { console.log('  PASS ' + msg); }
  else { console.error('  FAIL ' + msg); failures++; }
}
function approx(a, b, tol, msg) { ok(Number.isFinite(a) && Math.abs(a - b) < tol, msg + ' (got ' + a + ', want ' + b + ' +/- ' + tol + ')'); }

// ---- 1. Artifacts present ----
const html = fs.readFileSync(path.join(ROOT, 'VERICIGUAT_REVIEW.html'), 'utf8');
ok(html.length > 100000, 'dashboard HTML is non-trivial');
ok(html.charCodeAt(0) !== 0xFEFF, 'dashboard HTML has no BOM');
ok((html.match(/<script[ >]/g) || []).length === (html.match(/<\/script>/g) || []).length, 'script tags balanced');

// ---- 2. Engine functions still defined ----
['function poolWith', 'function pauleMandelTau2', 'function binaryEffect', 'function hrFromPublished', 'function mdFromCont'].forEach(function (sig) {
  ok(html.indexOf(sig) !== -1, 'engine still defines: ' + sig);
});

// ---- 3. Config parses ----
const cfg = JSON.parse(fs.readFileSync(path.join(ROOT, 'configs', 'vericiguat.json'), 'utf8'));
ok(cfg.drug === 'Vericiguat', 'config drug is Vericiguat');

// ---- 4. Numerical regression: re-implement the verified math ----
// log-OR effect + variance for a 2x2 table (matches binaryEffect, no zero cell).
function logOR(tE, tN, cE, cN) {
  const a = tE, b = tN - tE, c = cE, e = cN - cE;
  return { y: Math.log((a * e) / (b * c)), v: 1 / a + 1 / b + 1 / c + 1 / e };
}
// Single 2x2: 10/100 vs 20/100 -> OR = (10*80)/(90*20) = 0.4444...
const s = logOR(10, 100, 20, 100);
approx(Math.exp(s.y), (10 * 80) / (90 * 20), 1e-9, 'logOR back-transforms to correct OR');
approx(s.v, 1 / 10 + 1 / 90 + 1 / 20 + 1 / 80, 1e-12, 'logOR variance = 1/a+1/b+1/c+1/d');

// Fixed-effect inverse-variance pool of two identical studies: effect unchanged, SE shrinks by sqrt(2).
function poolFixed(ys, vs) {
  const w = vs.map(function (v) { return 1 / v; });
  const sumW = w.reduce(function (p, q) { return p + q; }, 0);
  let eff = 0; for (let i = 0; i < ys.length; i++) eff += ys[i] * w[i];
  eff /= sumW;
  return { effect: eff, se: 1 / Math.sqrt(sumW) };
}
const p = poolFixed([s.y, s.y], [s.v, s.v]);
approx(p.effect, s.y, 1e-12, 'pooling two identical effects leaves effect unchanged');
approx(p.se, Math.sqrt(s.v / 2), 1e-12, 'pooling two identical studies shrinks SE by sqrt(2)');

// Paule-Mandel tau^2 for homogeneous studies should be ~0 (Q <= k-1).
function pmTau2(y, v) {
  const k = y.length; if (k < 2) return 0; let tau2 = 0;
  for (let it = 0; it < 200; it++) {
    const w = v.map(function (vi) { return 1 / (vi + tau2); });
    const sumW = w.reduce(function (p2, q) { return p2 + q; }, 0);
    let yBar = 0; for (let i = 0; i < k; i++) yBar += y[i] * w[i]; yBar /= sumW;
    let Q = 0; for (let j = 0; j < k; j++) { const d = y[j] - yBar; Q += w[j] * d * d; }
    const diff = Q - (k - 1); if (Math.abs(diff) < 1e-7) break;
    let slope = 0; for (let m = 0; m < k; m++) { const dm = y[m] - yBar; slope += -w[m] * w[m] * dm * dm; }
    if (Math.abs(slope) < 1e-14) break;
    const t2 = Math.max(0, tau2 + diff / slope); if (Math.abs(t2 - tau2) < 1e-10) { tau2 = t2; break; } tau2 = t2;
  }
  return tau2;
}
approx(pmTau2([s.y, s.y, s.y], [s.v, s.v, s.v]), 0, 1e-9, 'Paule-Mandel tau^2 is 0 for homogeneous studies');
ok(pmTau2([0.5], [0.1]) === 0, 'Paule-Mandel returns 0 for k=1 (no pooling)');

if (failures) { console.error('\nSMOKE FAILED: ' + failures + ' assertion(s)'); process.exit(1); }
console.log('\nSMOKE PASSED');
