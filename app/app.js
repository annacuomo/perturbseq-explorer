// app.js: everything the page DOES.
//
// The flow:
//   1. init()       runs once when the page opens: downloads data/index.json
//                   (the list of perturbations) and wires up the controls.
//   2. select(id)   runs when you pick a gene: downloads that gene's JSON file.
//   3. draw()       turns the current record + ranking choice into a chart.
//
// Python analogies are in the comments where they help.

// ---------------------------------------------------------------- settings
// "const" = a variable that is never reassigned (like a module-level CONSTANT).
const COLORS = {
  down: "#1B6B73", downFaint: "#8CCBCF",
  up: "#7A2E4D", upFaint: "#C98BA6",
  grid: "#E8E1D5", ink: "#2A2A2A",
};
const SHOW_PER_SIDE = 20; // bars shown in each direction

const METRICS = {
  z:   { label: "change in expression (z-score vs controls)" },
  lfc: { label: "log2 fold change (approximate)" },
  abs: { label: "change in expression per cell (approximate)" },
};

// ---------------------------------------------------------------- state
// "let" = a variable that can change. These hold what the page currently shows.
let index = null;        // contents of data/index.json
let current = null;      // the perturbation record currently displayed
let metric = "z";        // "z" or "lfc"
const byId = new Map();  // like a Python dict: "GATA1" -> index entry

// $("chart") is shorthand for "find the element with id='chart' in index.html"
const $ = (id) => document.getElementById(id);

function showMessage(html) {
  $("chart").innerHTML = `<p class="placeholder">${html}</p>`;
}

// ---------------------------------------------------------------- 1. start-up
// "async" functions can "await" slow things (downloads) without freezing the page.
async function init() {
  try {
    const res = await fetch("data/index.json");   // like requests.get(...)
    if (!res.ok) throw new Error(res.status);
    index = await res.json();                     // like res.json() in requests
  } catch (err) {
    showMessage("Could not load data/index.json. Run the pipeline first, " +
      "and view through a local server (see README).");
    return;
  }

  // fill the search box's dropdown suggestions
  const list = $("perturbation-list");
  for (const p of index.perturbations) {          // like: for p in index["perturbations"]:
    byId.set(p.id.toUpperCase(), p);
    const opt = document.createElement("option");
    opt.value = p.id;
    list.appendChild(opt);
  }

  // `...${x}...` is a template string, like a Python f-string
  $("dataset-note").textContent =
    `${index.perturbations.length.toLocaleString()} perturbations · ` +
    `${index.n_genes_measured.toLocaleString()} genes measured · ${index.dataset}`;
  const fdrText = `${Math.round((index.fdr ?? 0.05) * 100)}%`;   // ?? = "if missing, use"
  $("fdr").textContent = fdrText;
  document.querySelectorAll(".fdr-copy").forEach((el) => (el.textContent = fdrText));

  // event listeners: "when X happens, call this function" (callbacks)
  $("search").addEventListener("change", (e) => select(e.target.value));
  $("search").addEventListener("keydown", (e) => { if (e.key === "Enter") select(e.target.value); });

  for (const btn of document.querySelectorAll(".toggle button")) {
    if (btn.dataset.metric !== "z" && !index.has_fold_change) btn.disabled = true;
    btn.addEventListener("click", () => setMetric(btn.dataset.metric));
  }

  // the part of the URL after # remembers the view, e.g. ...#GATA1&lfc
  const [hashId, hashMetric] = decodeURIComponent(location.hash.slice(1)).split("&");
  if (hashMetric in METRICS) setMetric(hashMetric, false);
  if (hashId) { $("search").value = hashId; select(hashId); }
}

// ---------------------------------------------------------------- 2. pick a gene
async function select(query) {
  const q = query.trim().toUpperCase();
  // exact match, or the first entry for that gene (e.g. "MYC" -> "MYC_P1P2")
  const entry = byId.get(q) ?? index.perturbations.find((p) => p.id.toUpperCase().startsWith(q + "_"));
  if (!entry) {
    showMessage(`No perturbation called “${query}” in this screen.`);
    $("summary").hidden = true;
    return;
  }
  current = await (await fetch("data/" + entry.file)).json();
  showSummary(current);
  draw();
}

function setMetric(m, redraw = true) {
  metric = m;
  for (const btn of document.querySelectorAll(".toggle button")) {
    btn.setAttribute("aria-pressed", String(btn.dataset.metric === m));
  }
  if (redraw && current) draw();
}

function showSummary(rec) {
  $("summary").hidden = false;
  $("stat-knockdown").textContent =
    rec.knockdown_pct != null ? `${rec.knockdown_pct}%` :
    rec.knockdown_z != null ? `z ${rec.knockdown_z.toFixed(1)}` : "n/a";
  $("stat-sig").textContent = rec.n_sig == null ? "n/a" : rec.n_sig.toLocaleString();
  $("stat-cells").textContent = rec.n_cells == null ? "n/a" : rec.n_cells.toLocaleString();
}

// ---------------------------------------------------------------- 3. the chart
function draw() {
  if (typeof Plotly === "undefined") {
    showMessage("The charting library (Plotly) didn’t load. Check your internet " +
      "connection, or open the browser console for details.");
    return;
  }
  const rec = current;
  history.replaceState(null, "", "#" + encodeURIComponent(rec.id) + (metric === "z" ? "" : "&" + metric));

  // pick the top genes by the chosen metric. In pandas this would be roughly:
  //   df = df[df[metric].notna()]
  //   down = df[df[metric] < 0].nsmallest(N, metric)
  //   up   = df[df[metric] > 0].nlargest(N, metric)
  const valid = rec.genes.filter((d) => d[metric] != null);
  const down = valid.filter((d) => d[metric] < 0).sort((a, b) => a[metric] - b[metric]).slice(0, SHOW_PER_SIDE);
  const up = valid.filter((d) => d[metric] > 0).sort((a, b) => b[metric] - a[metric]).slice(0, SHOW_PER_SIDE);

  // Plotly draws the first row at the BOTTOM, so order bottom-to-top:
  // strongest down ... weakest down, weakest up ... strongest up
  const rows = [...down, ...up.reverse()];   // ... "spreads" lists, like down + up[::-1]

  const color = (d) => d[metric] < 0 ? (d.sig ? COLORS.down : COLORS.downFaint)
                                     : (d.sig ? COLORS.up : COLORS.upFaint);
  const fmtLevel = (v) => (v == null ? "?" : v < 10 ? v.toFixed(2) : v.toFixed(1));

  // the text shown when you hover a bar (Plotly understands <b> and <br>)
  const hover = rows.map((d) =>
    `<b>${d.gene}</b>${d.sig ? "" : " (not significant)"}<br>` +
    `z-score: ${d.z.toFixed(2)}<br>` +
    (d.lfc != null ? `log2 fold change: ${d.lfc.toFixed(2)}<br>` +
      `approx. per cell: ${fmtLevel(d.before)} → ${fmtLevel(d.after)}` +
      ` (${d.abs > 0 ? "+" : ""}${fmtLevel(d.abs)})` : ""));

  const trace = {
    type: "bar",
    orientation: "h",
    x: rows.map((d) => d[metric]),   // like [d[metric] for d in rows]
    y: rows.map((d) => d.gene),
    marker: { color: rows.map(color) },
    text: hover,
    textposition: "none",            // use the text for hovering only, not as bar labels
    hovertemplate: "%{text}<extra></extra>",
  };

  const layout = {
    title: { text: `Knocking down <b>${rec.target}</b>`, font: { family: "Fraunces, serif", size: 20 } },
    font: { family: "Manrope, sans-serif", color: COLORS.ink },
    height: Math.max(360, rows.length * 20 + 110),
    margin: { l: 90, r: 24, t: 56, b: 56 },
    xaxis: { title: METRICS[metric].label, zeroline: true, zerolinecolor: COLORS.ink, gridcolor: COLORS.grid },
    // list every gene explicitly, otherwise Plotly hides every other label
    yaxis: { automargin: true, tickfont: { size: 11 }, type: "category",
             tickmode: "array", tickvals: rows.map((d) => d.gene), ticktext: rows.map((d) => d.gene) },
    hoverlabel: { bgcolor: "#fff", bordercolor: COLORS.grid, font: { family: "Manrope, sans-serif", color: COLORS.ink } },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    showlegend: false,
  };

  $("chart").innerHTML = "";
  Plotly.newPlot("chart", [trace], layout, { responsive: true, displaylogo: false });
}

init();
