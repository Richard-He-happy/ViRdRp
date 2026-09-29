const $ = (selector, root = document) => root.querySelector(selector);
const app = $("#app");
const state = {
  viruses: [],
  byAbbr: new Map(),
  legend: {},
  mclEdges: [],
  structureMatrix: null,
  sequenceMatrix: null,
  heatmapOrders: null,
  njtree: null,
  sample: null
};

const plddtColors = [
  { label: "Very high (pLDDT > 90)", min: 90, color: "#0d57d3" },
  { label: "High (90 > pLDDT > 70)", min: 70, color: "#6acbf1" },
  { label: "Low (70 > pLDDT > 50)", min: 50, color: "#fed936" },
  { label: "Very low (pLDDT < 50)", min: -Infinity, color: "#fd7d4d" }
];

const HMM_REFERENCE_FASTA = `>XXXX\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`;

const NJTREE_STYLE_CONFIG = {
  svgWidth: 680,
  svgHeight: 680,
  fitPadding: 42,
  showNodeLabels: false,
  leafNodeRadius: 3.6,
  nodeStrokeWidth: 0.8,
  linkStroke: "#9fb5af",
  linkStrokeWidth: 1.15,
  linkStrokeOpacity: 0.82,
  hitAreaRadius: 8,
  minZoom: 0.45,
  maxZoom: 5
};

async function json(path) {
  const res = await fetch(new URL(path, document.baseURI));
  if (!res.ok) throw new Error("Failed to load " + path);
  return res.json();
}

async function loadData() {
  const [viruses, legend, mclEdges, structureMatrix, sequenceMatrix, heatmapOrders, njtree, sample] = await Promise.all([
    json("./data/viruses.json"),
    json("./data/legend.json"),
    json("./data/mcl_edges.json"),
    json("./data/structure_matrix.json"),
    json("./data/sequence_matrix.json"),
    json("./data/heatmap_orders.json"),
    json("./data/njtree.json"),
    json("./data/sample-001.json")
  ]);
  Object.assign(state, { viruses, legend, mclEdges, structureMatrix, sequenceMatrix, heatmapOrders, njtree, sample });
  state.byAbbr = new Map(viruses.map(item => [item.abbreviation, item]));
}

function navigate(path) { location.hash = path; render(); }

document.addEventListener("click", event => {
  const disabled = event.target.closest("[data-disabled]");
  if (disabled) { event.preventDefault(); alert("该功能已被禁止"); return; }
  const link = event.target.closest("[data-route]");
  if (!link) return;
  event.preventDefault();
  navigate(link.getAttribute("href").replace(/^#/, "") || "/mcl");
});
window.addEventListener("hashchange", render);
$("#global-search").addEventListener("submit", event => { event.preventDefault(); alert("该功能已被禁止"); });

function page(title, subtitle, content, compact = false) {
  app.innerHTML = `
    <section class="page ${compact ? "page-compact" : ""}">
      <div class="section-head">
        <div>
          <h1>${title}</h1>
          ${subtitle ? `<p>${subtitle}</p>` : ""}
        </div>
      </div>
      ${content}
    </section>
  `;
}

function escapeHtml(value = "") {
  return String(value).replace(/[<>&"]/g, s => ({ "<": "&lt;", ">": "&gt;", "&": "&amp;", '"': "&quot;" }[s]));
}

function displayFor(abbr) { return state.byAbbr.get(abbr)?.displayLabel || "XXXX"; }

function legendHtml() {
  const renderGroup = group => `
    <section class="legend-group">
      <h3>${group}</h3>
      ${(state.legend[group] || []).map(item => `
        <div class="legend-row">
          <span class="swatch" style="background:${item.color}"></span>
          <span>${escapeHtml(item.family)}</span>
        </div>
      `).join("")}
    </section>
  `;
  return `
    <aside class="legend">
      ${renderGroup("(+)ssRNA")}
      ${renderGroup("(-)ssRNA")}
    </aside>
  `;
}

function graphShell(id, title) {
  return `
    <div class="panel canvas-panel">
      <div class="graph-toolbar">
        <strong>${title}</strong>
        <button class="control-btn" data-reset="${id}">Reset</button>
      </div>
      <div class="graph-host" id="${id}"></div>
    </div>
  `;
}

function colorFor(abbr) {
  return state.byAbbr.get(abbr)?.color || "#8ea4ad";
}

function renderMcl() {
  page("MCL Clustering", "Protein structures were clustered using the Markov Cluster Algorithm (MCL) based on pairwise TM-score similarities.", `
    <div class="workspace">
      ${legendHtml()}
      ${graphShell("mcl-graph", "MCL-based clustering of protein structures")}
    </div>
  `);
  drawMcl($("#mcl-graph"));
}

function mclNodeAbbreviation(data = {}) {
  return String(data.abbreviation || data.id || data.label || data.name || "").trim();
}

function openStructurePageFromMclNode(data) {
  navigate("/structure/Sample-001");
}
function drawMcl(host) {
  const width = host.clientWidth;
  const height = host.clientHeight;
  const labels = state.structureMatrix.labels.map(label => label.trim());
  const nodes = labels.map((id, i) => ({
    id,
    abbreviation: id,
    label: id,
    name: id,
    x: width / 2 + Math.cos(i / labels.length * Math.PI * 2) * Math.min(width, height) * 0.35,
    y: height / 2 + Math.sin(i / labels.length * Math.PI * 2) * Math.min(width, height) * 0.35,
    degree: 0
  }));
  const byId = new Map(nodes.map(node => [node.id, node]));
  const edges = state.mclEdges.filter(edge => byId.has(edge.source) && byId.has(edge.target));
  edges.forEach(edge => {
    byId.get(edge.source).degree++;
    byId.get(edge.target).degree++;
  });
  for (let tick = 0; tick < 180; tick++) {
    for (const edge of edges) {
      const a = byId.get(edge.source), b = byId.get(edge.target);
      const dx = b.x - a.x, dy = b.y - a.y;
      const distance = Math.hypot(dx, dy) || 1;
      const force = (distance - 92) * 0.006;
      a.x += dx / distance * force; a.y += dy / distance * force;
      b.x -= dx / distance * force; b.y -= dy / distance * force;
    }
    for (let i = 0; i < nodes.length; i++) {
      for (let j = i + 1; j < nodes.length; j++) {
        const a = nodes[i], b = nodes[j];
        const dx = b.x - a.x, dy = b.y - a.y;
        const d2 = dx * dx + dy * dy || 1;
        const force = Math.min(28, 4200 / d2);
        const d = Math.sqrt(d2);
        a.x -= dx / d * force; a.y -= dy / d * force;
        b.x += dx / d * force; b.y += dy / d * force;
      }
    }
    for (const node of nodes) {
      node.x += (width / 2 - node.x) * 0.005;
      node.y += (height / 2 - node.y) * 0.005;
    }
  }
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  const g = document.createElementNS(svg.namespaceURI, "g");
  svg.append(g);
  host.replaceChildren(svg);
  let scale = 1, tx = 0, ty = 0, dragging = false, last = null;
  const apply = () => g.setAttribute("transform", `translate(${tx} ${ty}) scale(${scale})`);
  for (const edge of edges) {
    const a = byId.get(edge.source), b = byId.get(edge.target);
    const line = document.createElementNS(svg.namespaceURI, "line");
    line.setAttribute("x1", a.x); line.setAttribute("y1", a.y);
    line.setAttribute("x2", b.x); line.setAttribute("y2", b.y);
    line.setAttribute("stroke", "#c9d7d4");
    line.setAttribute("stroke-width", String(Math.max(0.4, edge.value * 1.2)));
    g.append(line);
  }
  for (const node of nodes) {
    const group = document.createElementNS(svg.namespaceURI, "g");
    group.classList.add("node");
    group.setAttribute("transform", `translate(${node.x} ${node.y})`);
    group.style.pointerEvents = "all";
    group.addEventListener("pointerdown", event => event.stopPropagation());
    const circle = document.createElementNS(svg.namespaceURI, "circle");
    circle.setAttribute("r", String(5 + Math.sqrt(node.degree)));
    circle.setAttribute("fill", colorFor(node.id));
    circle.setAttribute("stroke", "#ffffff");
    circle.setAttribute("stroke-width", "1.5");
    const hitArea = document.createElementNS(svg.namespaceURI, "circle");
    hitArea.setAttribute("r", String(11 + Math.sqrt(node.degree)));
    hitArea.setAttribute("fill", "transparent");
    hitArea.setAttribute("pointer-events", "all");
    hitArea.addEventListener("click", event => {
      event.stopPropagation();
      openStructurePageFromMclNode(node);
    });
    circle.addEventListener("click", event => {
      event.stopPropagation();
      openStructurePageFromMclNode(node);
    });
    const title = document.createElementNS(svg.namespaceURI, "title");
    title.textContent = displayFor(node.abbreviation);
    group.append(title, hitArea, circle);
    group.addEventListener("click", event => {
      event.stopPropagation();
      openStructurePageFromMclNode(node);
    });
    const link = document.createElementNS(svg.namespaceURI, "a");
    link.setAttribute("href", "#/structure/Sample-001");
    link.setAttribute("aria-label", "Open Sample-001");
    link.append(group);
    g.append(link);
  }
  svg.addEventListener("wheel", event => {
    event.preventDefault();
    scale = Math.max(0.35, Math.min(5, scale * (event.deltaY < 0 ? 1.1 : 0.9)));
    apply();
  }, { passive: false });
  svg.addEventListener("pointerdown", event => { dragging = true; last = [event.clientX, event.clientY]; svg.setPointerCapture(event.pointerId); });
  svg.addEventListener("pointermove", event => {
    if (!dragging) return;
    tx += event.clientX - last[0]; ty += event.clientY - last[1];
    last = [event.clientX, event.clientY]; apply();
  });
  svg.addEventListener("pointerup", () => { dragging = false; });
  $("[data-reset='mcl-graph']").onclick = () => { scale = 1; tx = 0; ty = 0; apply(); };
}

function renderNjtree() {
  page("NJtree", "Neighbor-joining (NJ) tree of proteins based on pairwise structural and sequence similarities.", `
    <div class="workspace">
      ${legendHtml()}
      <div class="two-up">
        ${graphShell("tree-structure", "Structure similarity NJ tree")}
        ${graphShell("tree-sequence", "Sequence similarity NJ tree")}
      </div>
    </div>
  `);
  drawUnrootedTree($("#tree-structure"), state.njtree.structure, {
    leafNodeRadius: NJTREE_STYLE_CONFIG.leafNodeRadius / 2
  });
  drawUnrootedTree($("#tree-sequence"), state.njtree.sequence);
}

function treeLeaves(node, out = []) {
  if (!node.children || !node.children.length) out.push(node.name);
  else node.children.forEach(child => treeLeaves(child, out));
  return out;
}

function leafCount(node) {
  return (!node.children || !node.children.length)
    ? 1
    : node.children.reduce((sum, child) => sum + leafCount(child), 0);
}

function layoutTree(node, startAngle = 0, endAngle = Math.PI * 2, x = 0, y = 0) {
  if (!node.children || !node.children.length) return { ...node, x, y, leaf: true };
  const total = node.children.reduce((sum, child) => sum + leafCount(child), 0);
  let cursor = startAngle;
  const children = node.children.map(child => {
    const span = (endAngle - startAngle) * leafCount(child) / Math.max(1, total);
    const childStart = cursor;
    const childEnd = cursor + span;
    cursor = childEnd;
    const theta = (childStart + childEnd) / 2;
    const length = Math.pow(Math.max(0.015, child.length || 0.02), 0.85) * 210 * 1.05;
    return layoutTree(child, childStart, childEnd, x + Math.cos(theta) * length, y + Math.sin(theta) * length);
  });
  return { ...node, x, y, children };
}

function collectTreeNodes(node, out = []) {
  out.push(node);
  for (const child of node.children || []) collectTreeNodes(child, out);
  return out;
}

function computeFitTransform(nodes, width, height, padding) {
  const xs = nodes.map(node => node.x);
  const ys = nodes.map(node => node.y);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const dataWidth = Math.max(maxX - minX, 1);
  const dataHeight = Math.max(maxY - minY, 1);
  const availableWidth = Math.max(width - padding * 2, 1);
  const availableHeight = Math.max(height - padding * 2, 1);
  const scale = Math.min(availableWidth / dataWidth, availableHeight / dataHeight);
  return {
    scale,
    tx: width / 2 - scale * (minX + maxX) / 2,
    ty: height / 2 - scale * (minY + maxY) / 2
  };
}

function drawUnrootedTree(host, tree, styleOverrides = {}) {
  const config = { ...NJTREE_STYLE_CONFIG, ...styleOverrides };
  const width = host.clientWidth || config.svgWidth;
  const height = host.clientHeight || config.svgHeight;
  const layout = layoutTree(tree, -Math.PI / 2, Math.PI * 1.5, 0, 0);
  const allNodes = collectTreeNodes(layout);
  const initial = computeFitTransform(allNodes, width, height, config.fitPadding);
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  const g = document.createElementNS(svg.namespaceURI, "g");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("preserveAspectRatio", "xMidYMid meet");
  svg.append(g);
  host.replaceChildren(svg);

  function addEdges(node) {
    for (const child of node.children || []) {
      const line = document.createElementNS(svg.namespaceURI, "line");
      line.setAttribute("x1", node.x); line.setAttribute("y1", node.y);
      line.setAttribute("x2", child.x); line.setAttribute("y2", child.y);
      line.setAttribute("stroke", config.linkStroke);
      line.setAttribute("stroke-opacity", String(config.linkStrokeOpacity));
      line.setAttribute("stroke-width", String(config.linkStrokeWidth));
      line.setAttribute("vector-effect", "non-scaling-stroke");
      g.append(line);
      addEdges(child);
    }
  }
  function addNodes(node) {
    if (node.leaf) {
      const group = document.createElementNS(svg.namespaceURI, "g");
      const record = state.byAbbr.get(node.name);
      group.classList.add("tree-node");
      group.setAttribute("transform", `translate(${node.x} ${node.y})`);
      group.addEventListener("pointerdown", event => event.stopPropagation());
      group.addEventListener("click", event => {
        event.stopPropagation();
        navigate("/structure/Sample-001");
      });

      const title = document.createElementNS(svg.namespaceURI, "title");
      title.textContent = [displayFor(node.name), "Click to open Sample-001"].join("\n");

      const hitArea = document.createElementNS(svg.namespaceURI, "circle");
      hitArea.setAttribute("r", String(config.hitAreaRadius));
      hitArea.setAttribute("fill", "transparent");
      hitArea.setAttribute("pointer-events", "all");
      hitArea.addEventListener("click", event => {
        event.stopPropagation();
        navigate("/structure/Sample-001");
      });

      const dot = document.createElementNS(svg.namespaceURI, "circle");
      dot.setAttribute("r", String(config.leafNodeRadius));
      dot.setAttribute("fill", colorFor(node.name));
      dot.setAttribute("stroke", "#ffffff");
      dot.setAttribute("stroke-width", String(config.nodeStrokeWidth));
      dot.setAttribute("vector-effect", "non-scaling-stroke");
      dot.addEventListener("click", event => {
        event.stopPropagation();
        navigate("/structure/Sample-001");
      });

      group.append(title, hitArea, dot);
      if (config.showNodeLabels) {
        const text = document.createElementNS(svg.namespaceURI, "text");
        const angle = Math.atan2(node.y, node.x);
        text.textContent = displayFor(node.name);
        text.setAttribute("x", Math.cos(angle) * 10);
        text.setAttribute("y", Math.sin(angle) * 10);
        text.setAttribute("font-size", "9.2");
        text.setAttribute("fill", colorFor(node.name));
        text.setAttribute("text-anchor", Math.cos(angle) >= 0 ? "start" : "end");
        group.append(text);
      }
      const link = document.createElementNS(svg.namespaceURI, "a");
      link.setAttribute("href", "#/structure/Sample-001");
      link.setAttribute("aria-label", "Open Sample-001");
      link.append(group);
      g.append(link);
    } else {
      for (const child of node.children || []) addNodes(child);
    }
  }
  addEdges(layout);
  addNodes(layout);

  let scale = initial.scale, tx = initial.tx, ty = initial.ty, dragging = false, last = null;
  const apply = () => g.setAttribute("transform", `translate(${tx} ${ty}) scale(${scale})`);
  svg.addEventListener("wheel", event => { event.preventDefault(); scale = Math.max(config.minZoom, Math.min(config.maxZoom, scale * (event.deltaY < 0 ? 1.1 : 0.9))); apply(); }, { passive: false });
  svg.addEventListener("pointerdown", event => { dragging = true; last = [event.clientX, event.clientY]; });
  svg.addEventListener("pointermove", event => { if (!dragging) return; tx += event.clientX - last[0]; ty += event.clientY - last[1]; last = [event.clientX, event.clientY]; apply(); });
  svg.addEventListener("pointerup", () => { dragging = false; });
  svg.addEventListener("pointerleave", () => { dragging = false; });
  const reset = document.querySelector(`[data-reset='${host.id}']`);
  if (reset) reset.onclick = () => { scale = initial.scale; tx = initial.tx; ty = initial.ty; apply(); };
  apply();
}

function renderHeatmap() {
  page("Heatmap", "Heatmap visualization of proteins based on pairwise structural and sequence similarities.", `
    <div class="workspace">
      ${legendHtml()}
      <div class="two-up">
        ${graphShell("heat-structure", "Structure similarity heatmap")}
        ${graphShell("heat-sequence", "Sequence similarity heatmap")}
      </div>
    </div>
  `);
  drawHeatmap($("#heat-structure"), state.structureMatrix, state.heatmapOrders.structure);
  drawHeatmap($("#heat-sequence"), state.sequenceMatrix, state.heatmapOrders.sequence);
}

function heatColor(value) {
  const stops = ["#FFFFFF", "#F0EBF4", "#E0D6EB", "#D2C1E1", "#C3ACD8", "#B499CF", "#A584C5", "#9770BB", "#875AB2", "#7946A8", "#6A329F"];
  const v = Math.max(0, Math.min(1, value)) * (stops.length - 1);
  const i = Math.floor(v);
  const t = v - i;
  const parse = hex => [1, 3, 5].map(pos => parseInt(hex.slice(pos, pos + 2), 16));
  const a = parse(stops[i]);
  const b = parse(stops[Math.min(stops.length - 1, i + 1)]);
  const r = Math.round(a[0] + (b[0] - a[0]) * t);
  const g = Math.round(a[1] + (b[1] - a[1]) * t);
  const b2 = Math.round(a[2] + (b[2] - a[2]) * t);
  return `rgb(${r},${g},${b2})`;
}

function reorderMatrix(matrix, order) {
  const rowByName = new Map(matrix.rows.map(row => [row.name, row]));
  const indexByLabel = new Map(matrix.labels.map((label, i) => [label, i]));
  return order.map(name => ({
    name,
    values: order.map(label => rowByName.get(name).values[indexByLabel.get(label)])
  }));
}

function drawHeatmap(host, matrix, order) {
  const canvas = document.createElement("canvas");
  const ctx = canvas.getContext("2d");
  host.replaceChildren(canvas);
  const labels = order;
  const rows = reorderMatrix(matrix, order);
  let scale = 1, ox = 0, oy = 0, dragging = false, last = null;
  function resize() {
    canvas.width = host.clientWidth * devicePixelRatio;
    canvas.height = host.clientHeight * devicePixelRatio;
    canvas.style.width = `${host.clientWidth}px`;
    canvas.style.height = `${host.clientHeight}px`;
    draw();
  }
  function metrics() {
    const w = canvas.width / devicePixelRatio, h = canvas.height / devicePixelRatio;
    const pad = 88;
    const cell = Math.max(3, Math.min((w - pad - 18) / labels.length, (h - pad - 18) / labels.length) * scale);
    return { w, h, pad, cell };
  }
  function draw() {
    const { w, h, pad, cell } = metrics();
    ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, w, h);
    ctx.save();
    ctx.translate(pad + ox, pad + oy);
    rows.forEach((row, y) => row.values.forEach((value, x) => {
      ctx.fillStyle = heatColor(value);
      ctx.fillRect(x * cell, y * cell, cell + 0.35, cell + 0.35);
    }));
    ctx.restore();
    ctx.font = "10px Segoe UI, Arial";
    labels.forEach((label, i) => {
      ctx.fillStyle = colorFor(label);
      const pos = pad + ox + i * cell + cell / 2;
      if (pos > pad - 60 && pos < w + 60) {
        ctx.save();
        ctx.translate(pos, pad - 8);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(displayFor(label), 0, 0);
        ctx.restore();
      }
      const y = pad + oy + i * cell + cell * 0.72;
      if (y > pad - 20 && y < h + 20) ctx.fillText(displayFor(label), 10, y);
    });
  }
  canvas.addEventListener("wheel", event => { event.preventDefault(); scale = Math.max(0.45, Math.min(7, scale * (event.deltaY < 0 ? 1.12 : 0.9))); draw(); }, { passive: false });
  canvas.addEventListener("pointerdown", event => { dragging = true; last = [event.clientX, event.clientY]; });
  canvas.addEventListener("pointermove", event => { if (!dragging) return; ox += event.clientX - last[0]; oy += event.clientY - last[1]; last = [event.clientX, event.clientY]; draw(); });
  canvas.addEventListener("pointerup", () => { dragging = false; });
  canvas.addEventListener("click", () => navigate("/structure/Sample-001"));
  const reset = document.querySelector(`[data-reset='${host.id}']`);
  if (reset) reset.onclick = () => { scale = 1; ox = 0; oy = 0; draw(); };
  resize();
}

function renderStructure() {
  const item = state.sample;
  page("Sample-001 Structure", "", `
    <div class="structure-layout">
      <div class="panel">
        <div class="graph-toolbar">
          <strong>PDBe-Molstar Viewer</strong>
          <div><button class="control-btn" id="reset-view">Reset</button><button class="control-btn primary" id="structure-download">Download</button></div>
        </div>
        <div class="viewer" id="viewer"></div>
        <div class="plddt"><span><i class="swatch" style="background:#0d57d3"></i>Very high (pLDDT &gt; 90)</span><span><i class="swatch" style="background:#6acbf1"></i>High (90 &gt; pLDDT &gt; 70)</span><span><i class="swatch" style="background:#fed936"></i>Low (70 &gt; pLDDT &gt; 50)</span><span><i class="swatch" style="background:#fd7d4d"></i>Very low (pLDDT &lt; 50)</span></div>
      </div>
      <aside class="panel structure-info">
        <h2>Sample-001</h2>
        <p class="lead">ColabFold pLDDT Score: XXXX</p>
        <p>UniProt: XXXX</p>
        <div class="taxonomy">
          <div><b>Type:</b><span>XXX</span></div><div><b>Realm:</b><span>XXXX</span></div><div><b>Phylum:</b><span>XXXX</span></div><div><b>Class:</b><span>XXXX</span></div><div><b>Order:</b><span>XXXX</span></div><div><b>Family:</b><span>XXXX</span></div><div><b>Genus:</b><span>XXXX</span></div><div><b>Virus:</b><span>XXXX</span></div>
        </div>
        <h2>Sequence</h2><div class="sequence">XXXXXXXXXXXXXXXX</div>
      </aside>
    </div>
  `);
  $("#structure-download").onclick = () => alert("该功能已被禁止");
  drawMolstarViewer($("#viewer"), item);
}

function disposeMolstarViewer() {
  const plugin = window.currentMolstarPlugin;
  if (!plugin) return;
  try {
    plugin.plugin?.dispose?.();
  } catch {
    // PDBe-Molstar cleanup APIs vary by version.
  }
  window.currentMolstarPlugin = null;
}

function drawMolstarViewer(host, item) {
  disposeMolstarViewer();
  const shell = document.createElement("div");
  shell.className = "molstar-shell molstar-shell--minimal";
  const container = document.createElement("div");
  container.className = "molstar-host";
  shell.append(container);
  host.replaceChildren(shell);

  if (!window.PDBeMolstarPlugin) {
    container.innerHTML = `<div class="notice">PDBe-Molstar bundle did not load.</div>`;
    return;
  }

  const options = {
    customData: {
      url: new URL(item.pdb, document.baseURI).href,
      format: "pdb",
      binary: false
    },
    alphafoldView: true,
    sequencePanel: true,
    hideStructure: ["water"],
    granularity: "residueInstances",
    bgColor: { r: 255, g: 255, b: 255 },
    hideControls: true,
    hideCanvasControls: [],
    landscape: true
  };

  const plugin = new window.PDBeMolstarPlugin();
  window.currentMolstarPlugin = plugin;
  const renderResult = plugin.render(container, options);
  if (renderResult?.catch) {
    renderResult.catch(error => {
      console.error("PDBe-Molstar render failed", error);
      if (window.currentMolstarPlugin === plugin) {
        container.innerHTML = `<div class="notice">PDBe-Molstar could not load this PDB file.</div>`;
      }
    });
  }

  const reset = $("#reset-view");
  if (reset) {
    reset.onclick = () => {
      const currentPlugin = window.currentMolstarPlugin;
      if (currentPlugin?.visual?.reset) {
        currentPlugin.visual.reset({ camera: true }).catch(() => drawMolstarViewer(host, item));
      } else {
        drawMolstarViewer(host, item);
      }
    };
  }
}

function renderHmm() {
  page("HMM_search", "Interface preview only; execution is disabled in the public demo.", `
    <div class="panel upload"><div class="hmm-reference"><h2>Reference FASTA format</h2><p>The example is intentionally anonymized.</p><pre class="hmm-example-fasta">${escapeHtml(HMM_REFERENCE_FASTA)}</pre></div><form id="hmm-form" class="upload"><input class="field" type="file" accept=".fasta,.fa,.txt" disabled aria-disabled="true" /><button class="control-btn primary" type="submit">Run HMM_search</button></form><div class="notice">该功能已被禁止</div></div>
  `);
  $("#hmm-form").addEventListener("submit", event => { event.preventDefault(); alert("该功能已被禁止"); });
}

function renderAbout() {
  page("About", "", `
    <div class="about-card panel">
      <p>This public version preserves the original visualization framework while replacing research identifiers, labels, sequences, matrices, and metadata with anonymized demonstration data.</p>
      <p>Sample-001 and the X placeholders are presentation-only values and do not represent the underlying research records.</p>
    </div>
  `, true);
}

function currentPath() { return location.hash ? location.hash.slice(1) || "/mcl" : "/mcl"; }
function render() {
  const path = currentPath();
  if (path === "/" || path === "/mcl") return renderMcl();
  if (path === "/njtree") return renderNjtree();
  if (path === "/heatmap") return renderHeatmap();
  if (path === "/about") return renderAbout();
  if (path === "/hmm-search") return renderHmm();
  if (path === "/structure/Sample-001") return renderStructure();
  page("Page Not Found", "", "<div class=\"notice\">This page is not available in the public demo.</div>", true);
}

loadData().then(render).catch(error => {
  app.innerHTML = `<section class="page"><div class="notice">${escapeHtml(error.message)}</div></section>`;
});
