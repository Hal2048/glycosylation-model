/**
 * draw.js — GlycoGlyph headless drawer.
 *
 * Input: process.argv[2] = IUPAC-condensed glycan string.
 * Output: a JSON object written to stdout — {"svg": "...", "width": N, "height": N, "error": "..."}.
 *        Or, if an output file path is supplied (and --json is NOT set), the SVG is written to that file.
 *
 * Design notes:
 *  - No real browser required; we use jsdom to emulate a DOM.
 *  - d3 needs `document` at module-load time, so we inject the jsdom globals into the Node
 *    global scope BEFORE `require('d3')`. Otherwise d3 fails with "document is not defined".
 *  - d3 only issues HTTP-like calls when `d3.xml(url)` / `d3.text(url)` are invoked. Our drawing
 *    pipeline never calls those (we do not invoke the original savesvg), so we simply stub
 *    `window.fetch` to fail fast if anything accidentally tries a network request.
 *  - The SVG emitted by GlycoGlyph uses `<use xlink:href="about:///.../monos.svg#xxx">` which
 *    is meaningless for standalone files. We mirror the original savesvg.js replaceuse step
 *    and inline each referenced `<symbol>` as a real `<svg viewBox=...>`, producing a fully
 *    self-contained SVG.
 */

const fs = require('fs');
const path = require('path');
const {JSDOM} = require('jsdom');

// ---------- CLI parsing ----------
// Convention: if --json / -j is present, write JSON to stdout. Otherwise write the SVG to the
// positional argument after the IUPAC string (if provided).
const args = process.argv.slice(2);
const jsonMode = args.includes('--json') || args.includes('-j');
const filtered = args.filter(a => a !== '--json' && a !== '-j');
if (filtered.length < 1) {
  console.error('Usage: node draw.js <IUPAC> [outputFile] [width] [height] [orientation] [--json]');
  process.exit(2);
}
const iupac = filtered[0];
const outFile = jsonMode ? null : (filtered[1] || null);
const width = parseInt(filtered[2] || '320', 10);
const height = parseInt(filtered[3] || '500', 10);
const orientation = filtered[4] || 'bottom-to-top'; // 'bottom-to-top' or 'right-to-left'

// ---------- Asset paths ----------
const ROOT = path.resolve(__dirname, '..');
const GLYCOGLYPH_JS = path.join(ROOT, 'js', 'glycoglyph.min.js');
const MONOS_SVG_FILE = path.join(ROOT, 'assets', 'monos.svg');
const CSS_FILE = path.join(ROOT, 'assets', 'glycoglyph.css');

const glycoglyphSrc = fs.readFileSync(GLYCOGLYPH_JS, 'utf8');
const monosSvgText = fs.readFileSync(MONOS_SVG_FILE, 'utf8');
const cssText = fs.readFileSync(CSS_FILE, 'utf8');

// ---------- jsdom bootstrap ----------
// runScripts: 'outside-only' lets us inject scripts via window.eval but stops jsdom from
// trying to auto-execute <script> tags in the HTML.
const html = `<!DOCTYPE html><html><head><meta charset="utf-8"></head>
<body><div id="d3glycanstruc"></div></body></html>`;
const dom = new JSDOM(html, {pretendToBeVisual: true, runScripts: 'outside-only'});
const {window} = dom;

// Expose jsdom's key globals on the Node global. d3 reads them once at module-init time,
// so the injection must happen before `require('d3')`.
global.document = window.document;
global.window = window;
global.Element = window.Element;
global.HTMLElement = window.HTMLElement;
global.SVGElement = window.SVGElement;
global.Node = window.Node;
global.navigator = window.navigator;
global.URL = window.URL;
global.DOMParser = window.DOMParser;

const d3 = require('d3');

// Expose d3 and the assets to glycoglyph. The library reads `window.location.protocol/host` to
// build filePaths; we point it at a sensible `file://.../assets` URL. The exact path does not
// matter because we only consume the #id fragment when inlining symbols later.
window.d3 = d3;
window.URL = window.URL;
try {
  Object.defineProperty(window, 'location', {
    value: {protocol: 'file:', host: '', hostname: '', pathname: '/', href: 'file://' + ROOT.replace(/\\/g, '/')},
    configurable: true
  });
} catch (_) {}

// Stub any accidental network access — drawing must be fully offline.
window.fetch = () => Promise.reject(new Error('fetch disabled'));
window.glycoglyphBase = path.dirname(MONOS_SVG_FILE);

// ---------- Load GlycoGlyph ----------
window.eval(glycoglyphSrc);

// ---------- Invoke drawing ----------
const opts = {
  editingMode: false,            // disable click/mouseover editing handlers
  userOverrideOptions: false,    // do not read settings from a DOM form
  drawingareachoice: 'fixglycansize',
  orientation,
  width,
  height,
  symbsize: 25,
  linkageVisible: true,
  linkRotate: true,
  linkAbbr: false,
  linkFontSize: 12,
  margin: {top: 40, right: 40, bottom: 40, left: 40}
};

let svgEl;
let error = null;
try {
  window.glycoglyph.d3glycanstructure(iupac, opts);
  svgEl = window.document.getElementById('d3glycanstrucSVG');
  if (!svgEl) throw new Error('SVG element not produced');
} catch (e) {
  error = e && (e.stack || e.message) || String(e);
}

// ---------- Inline <symbol> references ----------
// Replace every <use href="...monos.svg#ID"> with a real inline <svg viewBox=...> containing the
// original symbol shapes, so the resulting SVG has no external file dependencies.
function inlineSymbols(svgRoot, monosDoc) {
  const symbols = Array.from(monosDoc.querySelectorAll('symbol'));
  const byId = new Map();
  for (const s of symbols) {
    const id = s.getAttribute('data-abbr');
    if (id) byId.set(id.toLowerCase(), {viewBox: s.getAttribute('viewBox'), inner: s.innerHTML.replace(/\n|\t/g, '').trim()});
  }
  const uses = Array.from(svgRoot.querySelectorAll('use'));
  for (const u of uses) {
    const href = u.getAttribute('href') || u.getAttribute('xlink:href') || '';
    const hashIdx = href.indexOf('#');
    if (hashIdx < 0) continue;
    const id = href.slice(hashIdx + 1).toLowerCase();
    const sym = byId.get(id);
    if (!sym) continue;
    const x = u.getAttribute('x') || '0';
    const y = u.getAttribute('y') || '0';
    const w = u.getAttribute('width') || '25';
    const h = u.getAttribute('height') || '25';
    // Build an inline <svg> in the jsdom document
    const inner = window.document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    inner.setAttribute('x', x);
    inner.setAttribute('y', y);
    inner.setAttribute('width', w);
    inner.setAttribute('height', h);
    if (sym.viewBox) inner.setAttribute('viewBox', sym.viewBox);
    // Inject the symbol's original inner shapes via innerHTML
    inner.innerHTML = sym.inner;
    u.parentNode.insertBefore(inner, u);
    u.parentNode.removeChild(u);
  }
}

// Also inline the small set of CSS rules that actually affect rendered nodes (font, cursor,
// etc.), so the SVG looks identical in viewers that ignore external CSS files. UI-only rules
// are skipped on purpose — they target elements that are not part of the drawing.
function inlineCss(svgRoot, cssRaw) {
  const cleaned = cssRaw.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\r?\n|\r/g, '');
  cleaned.split('}').forEach(rule => {
    const [selRaw, bodyRaw] = rule.split('{');
    if (!selRaw || !bodyRaw) return;
    const selector = selRaw.trim();
    const decls = bodyRaw.split(';').map(s => s.trim()).filter(Boolean);
    if (!decls.length) return;
    // Only consider selectors that target glyph nodes; ignore everything UI-related.
    const allowed = /(^|\s|,)(svg|text|path|use|\.symbol)(\s|$|,|\.|:)/;
    if (!allowed.test(selector)) return;
    svgRoot.querySelectorAll('*').forEach(el => {
      if (el === svgRoot) return;
      // Simple matching by tag or class only.
      const tag = el.tagName.toLowerCase();
      const cls = el.getAttribute('class') || '';
      let match = false;
      selector.split(',').forEach(part => {
        const p = part.trim();
        if (p === tag) match = true;
        if (p.startsWith('.') && cls.split(/\s+/).includes(p.slice(1))) match = true;
      });
      if (!match) return;
      decls.forEach(d => {
        const [prop, ...rest] = d.split(':');
        const value = rest.join(':').trim();
        if (!prop || !value) return;
        el.setAttribute(prop, value);
      });
    });
  });
}

let svgString = null;
if (svgEl) {
  // Add xmlns so the file is a valid standalone SVG when saved.
  svgEl.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  svgEl.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  // Parse monos.svg to look up symbol definitions by id.
  const monosDoc = new window.DOMParser().parseFromString(monosSvgText, 'image/svg+xml');
  inlineSymbols(svgEl, monosDoc);
  // Inline the necessary style rules.
  inlineCss(svgEl, cssText);
  svgString = svgEl.outerHTML;
}

const result = {
  svg: svgString,
  width: svgEl ? parseFloat(svgEl.getAttribute('width')) : width,
  height: svgEl ? parseFloat(svgEl.getAttribute('height')) : height,
  error
};

if (outFile) {
  if (!svgString) {
    process.stderr.write('ERROR: ' + error + '\n');
    process.exit(1);
  }
  fs.writeFileSync(outFile, svgString, 'utf8');
  process.stderr.write(`Wrote ${outFile} (${svgString.length} bytes)\n`);
} else {
  process.stdout.write(JSON.stringify(result));
}