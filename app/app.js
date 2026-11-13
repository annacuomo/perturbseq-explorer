// app.js: everything the page DOES.
//
// The flow:
//   1. init()        runs once when the page opens: downloads data/index.json
//                    (the list of perturbations) and wires up the controls.
//   2. select(id)    runs when you pick a gene: downloads that gene's JSON file.
//   3. draw()        turns the current record into the bar chart.
//   4. showEnrichment() fills the "What do these genes have in common?" panel.
//   5. openEnrichr() sends a gene list to the Enrichr website.
//
// Python analogies are in the comments where they help.

// ---------------------------------------------------------------- settings
// "const" = a variable that is never reassigned (like a module-level CONSTANT).
const COLORS = {
  down: "#1B6B73", up: "#7A2E4D",
  grid: "#E8E1D5", ink: "#2A2A2A",
};
const SHOW_PER_SIDE = 20; // bars shown in each direction
const ENRICHR = "https://maayanlab.cloud/Enrichr";
const MIN_CHECKED_CELLS = 25; // the smallest control has 25 cells: below this the test is unchecked

// ---------------------------------------------------------------- state
// "let" = a variable that can change. These hold what the page currently shows.
let index = null;        // contents of data/index.json
let current = null;      // the perturbation record currently displayed
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
  for (const btn of document.querySelectorAll("button.enrichr")) {
    btn.addEventListener("click", () => openEnrichr(btn.dataset.direction));
  }

  // the part of the URL after # remembers the gene, e.g. ...#GATA1
  // (older links like #GATA1&lfc still work: the part after & is ignored)
  const [hashId] = decodeURIComponent(location.hash.slice(1)).split("&");
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
    $("few-cells").hidden = true;
    $("enrichment").hidden = true;
    return;
  }
  current = await (await fetch("data/" + entry.file)).json();
  history.replaceState(null, "", "#" + encodeURIComponent(current.id));
  showSummary(current);
  draw();
  showEnrichment(current);
}

function showSummary(rec) {
  $("summary").hidden = false;
  $("stat-knockdown").textContent =
    rec.knockdown_pct != null ? `${rec.knockdown_pct}%` :
    rec.knockdown_z != null ? `z ${rec.knockdown_z.toFixed(1)}` : "n/a";
  $("stat-sig").textContent = rec.n_sig == null ? "n/a" : rec.n_sig.toLocaleString();
  $("stat-cells").textContent = rec.n_cells == null ? "n/a" : rec.n_cells.toLocaleString();
  $("few-cells").hidden = !(rec.n_cells != null && rec.n_cells < MIN_CHECKED_CELLS);
}

// The genes behind the chart: significant ones only, ranked by fold change.
// (Files made before the pipeline change also hold non-significant genes, so filter here too.)
function significantGenes(rec) {
  const sig = rec.genes.filter((d) => d.sig && d.lfc != null);
  return {
    down: sig.filter((d) => d.lfc < 0).sort((a, b) => a.lfc - b.lfc),   // most negative first
    up: sig.filter((d) => d.lfc > 0).sort((a, b) => b.lfc - a.lfc),     // most positive first
  };
}

// ---------------------------------------------------------------- 3. the chart
function draw() {
  if (typeof Plotly === "undefined") {
    showMessage("The charting library (Plotly) didn’t load. Check your internet " +
      "connection, or open the browser console for details.");
    return;
  }
  const rec = current;
  const { down, up } = significantGenes(rec);
  if (down.length + up.length === 0) {
    showMessage(`No gene changes significantly when <b>${rec.target}</b> is knocked down.`);
    return;
  }

  // Plotly draws the first row at the BOTTOM, so order bottom-to-top:
  // strongest down ... weakest down, weakest up ... strongest up
  const rows = [...down.slice(0, SHOW_PER_SIDE), ...up.slice(0, SHOW_PER_SIDE).reverse()];

  const fmtLevel = (v) => (v == null ? "?" : v < 10 ? v.toFixed(2) : v.toFixed(1));
  // the text shown when you hover a bar (Plotly understands <b> and <br>)
  const hover = rows.map((d) =>
    `<b>${d.gene}</b><br>` +
    `log2 fold change: ${d.lfc > 0 ? "+" : ""}${d.lfc.toFixed(2)}<br>` +
    `approx. per cell: ${fmtLevel(d.before)} → ${fmtLevel(d.after)}<br>` +
    `z-score: ${d.z.toFixed(2)}`);

  const trace = {
    type: "bar",
    orientation: "h",
    x: rows.map((d) => d.lfc),       // like [d["lfc"] for d in rows]
    y: rows.map((d) => d.gene),
    marker: { color: rows.map((d) => (d.lfc < 0 ? COLORS.down : COLORS.up)) },
    text: hover,
    textposition: "none",            // use the text for hovering only, not as bar labels
    hovertemplate: "%{text}<extra></extra>",
  };

  const layout = {
    title: { text: `Knocking down <b>${rec.target}</b>`, font: { family: "Fraunces, serif", size: 20 } },
    font: { family: "Manrope, sans-serif", color: COLORS.ink },
    height: Math.max(300, rows.length * 20 + 110),
    margin: { l: 90, r: 24, t: 56, b: 56 },
    xaxis: { title: "log2 fold change vs controls (approximate)", zeroline: true,
             zerolinecolor: COLORS.ink, gridcolor: COLORS.grid },
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

// ---------------------------------------------------------------- 4. enrichment panel
// The gene list sent to Enrichr: the pipeline's top list if present,
// otherwise the significant genes in this file.
function geneList(rec, direction) {
  return rec.top_lists?.[direction] ?? significantGenes(rec)[direction].map((d) => d.gene);
}

// "mitochondrial translation (GO:0032543)" -> ["mitochondrial translation", "GO:0032543"]
function splitTerm(term) {
  const m = term.match(/^(.*?)\s*\((GO:\d+)\)$/);
  return m ? [m[1], m[2]] : [term, null];
}

// Escape text before putting it into HTML (gene and term names come from data files).
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

function showEnrichment(rec) {
  $("enrichment").hidden = false;
  $("enrichr-status").textContent = "";
  const hasPrecomputed = rec.enrichment != null;
  const nTop = index.enrich_top ?? 100;
  $("enrichment-note").textContent = hasPrecomputed
    ? `Gene sets over-represented among the top ${nTop} significant genes in each direction ` +
      `(Hallmark and GO biological process; background: all ${index.n_genes_measured.toLocaleString()} measured genes).`
    : "Send the significant genes to Enrichr to see which pathways and processes they share.";

  for (const direction of ["up", "down"]) {
    const col = $(`enrich-${direction}`);
    const genes = geneList(rec, direction);
    const ol = col.querySelector("ol");
    const btn = col.querySelector("button");
    btn.disabled = genes.length === 0;
    btn.textContent = genes.length ? `Explore these ${genes.length} genes in Enrichr ↗` : "No genes to explore";

    if (!hasPrecomputed) { ol.innerHTML = ""; continue; }
    const terms = rec.enrichment[direction];
    if (genes.length < 5) {
      ol.innerHTML = `<li class="empty">Too few significant genes to test.</li>`;
    } else if (terms.length === 0) {
      ol.innerHTML = `<li class="empty">No gene set passes q &lt; ${Math.round((index.fdr ?? 0.05) * 100) / 100}.</li>`;
    } else {
      // one <li> per enriched gene set: name, source, overlap, q-value, and the genes involved
      ol.innerHTML = terms.map((t) => {
        const [name, goId] = splitTerm(t.term);
        return `<li>
          <span class="term" title="${esc(goId ?? t.term)}">${esc(name)}</span>
          <span class="lib">${esc(t.lib)}</span>
          <span class="meta">${t.k} of ${t.n_list} genes · q = ${t.q.toPrecision(2)}</span>
          <span class="genes">${t.genes.map(esc).join(", ")}${t.k > t.genes.length ? ", …" : ""}</span>
        </li>`;
      }).join("");
    }
  }
}

// ---------------------------------------------------------------- 5. Enrichr
// Enrichr's API: POST the list to /addList, get back a short id, open /enrich?dataset=<id>.
// The new tab is opened straight away (browsers block pop-ups opened later, after a download),
// then pointed at the results once Enrichr answers.
async function openEnrichr(direction) {
  const rec = current;
  const genes = geneList(rec, direction);
  const tab = window.open("about:blank", "_blank");
  $("enrichr-status").textContent = "Sending the list to Enrichr…";
  try {
    const form = new FormData();                  // like the data= of requests.post
    form.append("list", genes.join("\n"));
    form.append("description", `${rec.id} knockdown, genes going ${direction} (Replogle 2022 K562)`);
    const res = await fetch(`${ENRICHR}/addList`, { method: "POST", body: form });
    if (!res.ok) throw new Error(res.status);
    const { shortId } = await res.json();
    if (tab) tab.location = `${ENRICHR}/enrich?dataset=${shortId}`;
    else window.open(`${ENRICHR}/enrich?dataset=${shortId}`, "_blank");
    $("enrichr-status").textContent = "";
  } catch (err) {
    // plan B: put the genes on the clipboard and open Enrichr's front page to paste them
    let copied = false;
    try { await navigator.clipboard.writeText(genes.join("\n")); copied = true; } catch (e) { /* clipboard blocked */ }
    if (tab) tab.location = `${ENRICHR}/`;
    $("enrichr-status").textContent = copied
      ? "Couldn’t send the list automatically. The genes are on your clipboard: paste them into Enrichr’s box."
      : `Couldn’t send the list automatically. Paste these genes into Enrichr: ${genes.join(", ")}`;
  }
}

init();
