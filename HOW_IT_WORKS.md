# How this project works: a guide for a Python person

This explains the non-Python parts: what the browser does, what each web file
is for, why we need a local server, and how GitHub turns the repo into a live
website. Read it alongside the code; `app/app.js` is commented with Python
analogies.

---

## 1. The big picture

```
  YOUR LAPTOP (once, offline)                 ANYONE'S BROWSER (every visit)
  ───────────────────────────                 ──────────────────────────────
  .h5ad  ──►  pipeline/03_summarise.py  ──►  app/data/*.json  ──►  app/ (HTML+CSS+JS)  ──►  chart
  (1 GB)      (numpy, scipy, BH, ...)         (small files)        reads the JSON
```

- **All the science happens in Python, ahead of time.** The output is a folder of
  small JSON files.
- **The website does no statistics.** It downloads one small JSON file and draws it.
- **The JSON format is the contract** between the two halves (documented in the
  README). Change one side → check the other.

This is called a **static site**: no server running Python, no database. Just
files. That's why it can be hosted free on GitHub Pages.

---

## 2. What happens when someone opens the page

1. The browser asks the server for `index.html`.
2. Reading it, it finds links to `style.css`, `app.js`, the Google fonts and
   Plotly, and requests each of those.
3. `app.js` runs. Its `init()` function requests `data/index.json` (the list of
   2,176 perturbations) and fills the search box suggestions.
4. You pick GATA1 → `select()` requests `data/perturbations/GATA1.json` (~8 KB).
5. `draw()` sorts those genes by the chosen metric and hands them to Plotly,
   which draws the SVG chart.

So the visitor only ever downloads the index plus one small file per gene they
look at, never the full dataset.

---

## 3. The three files in `app/`

| File | Role | Python-world analogy |
|---|---|---|
| `index.html` | **Structure**: what's on the page (title, search box, toggle, empty chart box) | the layout of a figure: which panels exist |
| `style.css` | **Appearance**: colours, fonts, spacing, the toggle's look | a matplotlib style sheet / rcParams |
| `app.js` | **Behaviour**: download data, react to clicks, draw the chart | your script |

They connect through **ids**. In `index.html`, `<section id="chart">` is the
empty box; in `app.js`, `$("chart")` finds it and Plotly draws into it. In
`style.css`, `.chart { ... }` styles anything with `class="chart"`.

Things worth knowing in `style.css`: colours are defined once at the top as
variables (`--teal: #1B6B73`) and reused as `var(--teal)`, like constants.

---

## 4. JavaScript for a Python person

| JavaScript | Python |
|---|---|
| `const x = 5;` / `let x = 5;` | `x = 5` (`const` can't be reassigned, `let` can) |
| `[1, 2, 3]`, `{gene: "GATA1", z: -1.2}` | list, dict |
| `rows.map((d) => d.z)` | `[d["z"] for d in rows]` |
| `rows.filter((d) => d.z < 0)` | `[d for d in rows if d["z"] < 0]` |
| `rows.sort((a, b) => a.z - b.z)` | `rows.sort(key=lambda d: d["z"])` |
| `` `${n} genes` `` | `f"{n} genes"` |
| `x ?? 0.05` | `x if x is not None else 0.05` |
| `[...a, ...b]` | `a + b` |
| `await fetch(url)` then `.json()` | `requests.get(url).json()` |
| `null` | `None` |

Two ideas with no exact Python equivalent:

- **Events.** The page doesn't run top to bottom and finish. It sets up
  *listeners* ("when the toggle is clicked, call `setMetric`") and then waits.
- **`async` / `await`.** Downloads are slow, so `fetch` returns a promise of a
  result; `await` pauses that function (not the whole page) until it arrives.

---

## 5. Why `python -m http.server`

Double-clicking `index.html` opens it as `file:///...`. For security, browsers
won't let a page loaded that way `fetch()` other files, so the JSON never loads.

`python -m http.server 8000 --directory app` starts a tiny web server that
serves the `app` folder exactly as GitHub Pages will.

- `localhost` = "this computer". `8000` = the **port**, like a door number.
- `--directory app` is **relative to where your terminal is**. Run it from the
  `perturbseq-explorer` folder. (From inside `pipeline/` it would be `../app`.)
- **404 for `/`**: it's serving the wrong folder (see above).
- **"Address already in use"**: an old server still holds port 8000. Use 8001,
  or stop it: `lsof -ti :8000 | xargs kill`.
- **Ctrl+C** stops the server.
- After editing a file, just **refresh** the browser. If it looks stale,
  hard-refresh: **Cmd+Shift+R**.

**Browser caching.** Browsers keep copies of CSS and JS files to load pages
faster, so after a change you may still see the old version. That's why the
toggle once looked like plain grey buttons: the page was using an old
`style.css`. Two defences:
- hard-refresh (Cmd+Shift+R) while developing;
- in `index.html` the links read `style.css?v=4` and `app.js?v=4`. The browser
  treats a new `?v=` as a new file, so **bump the number when you change
  those files**, and visitors to the live site get the new version too.

---

## 6. Your debugger: the browser's developer tools

This is the single most useful skill for web work.

- **Chrome:** Cmd+Option+J opens the **Console**.
- **Safari:** first Settings → Advanced → tick "Show features for web developers",
  then Cmd+Option+C.

What to look at:
- **Console**: JavaScript errors in red, with file and line number. You can also
  type JS there, e.g. `current.genes.length` to inspect the loaded record.
- **Network** tab: every file the page requested, with its status (200 = OK,
  404 = not found). The blank-chart bug earlier would have shown as a failed
  request for `plotly-basic.min.js`.

When something looks wrong, paste me the red Console lines.

---

## 7. Plotly and the CDN

Plotly is the charting library (the same project as Python's `plotly`). The page
loads it from a **CDN** (content delivery network), a public server that hosts
popular libraries, via the `<script src="https://cdnjs...">` line in
`index.html`. The version number in that URL must exist, which is what broke the
first chart (2.35.2 isn't hosted; 2.35.3 is).

---

## 8. Git and GitHub

**Concepts**
- **Repository (repo)**: the project folder plus its full history.
- **Commit**: a saved snapshot with a message. Local only, until you push.
- **Push**: upload your commits to GitHub.
- **`.gitignore`**: files git should never track. Ours excludes the `.h5ad`,
  the `.cache.pkl` and `__pycache__`.
- **GitHub Actions**: scripts GitHub runs for you on each push. Ours
  (`.github/workflows/pages.yml`) copies `app/` to **GitHub Pages**, which hosts
  it at `https://annacuomo.github.io/perturbseq-explorer/`.

**One-time setup**

1. On github.com: **New repository** → name `perturbseq-explorer` → Public →
   tick **nothing** (no README, .gitignore or licence: we have them, and an
   extra commit on GitHub makes the first push fail) → Create.
2. In the new repo: **Settings → Pages → Source: GitHub Actions**. Do this
   *before* pushing, otherwise the first deploy fails with "Get Pages site
   failed" (if that happens: enable it now, then Actions → the failed run →
   Re-run jobs).
3. Use GitHub's private email for commits: copy the
   `…@users.noreply.github.com` address from github.com/settings/emails, then
   `git config --global user.email "that-address"`. (If you already committed
   with your personal email, GitHub will refuse the push; fix with
   `git commit --amend --reset-author --no-edit`.)
4. In the terminal, from the `perturbseq-explorer` folder:
   ```bash
   git init
   git add .
   git status            # check: no .h5ad or .cache.pkl listed
   git commit -m "Perturb-seq explorer: pipeline and first view"
   git branch -M main
   git remote add origin https://github.com/annacuomo/perturbseq-explorer.git
   git push -u origin main
   ```
   If the push asks for a password, GitHub needs a token or the GitHub CLI:
   `brew install gh && gh auth login` is the easiest route.
5. **Actions** tab: watch the "Deploy app to GitHub Pages" run go green (~1 min).
   The URL appears there and under Settings → Pages.

**Every change after that**
```bash
git add .
git commit -m "what you changed"
git push
```
The site updates about a minute later.

Regenerated the JSON? Those files are in `app/data/`, so they're committed and
deployed like everything else. That's intended: the live site needs them.

---

## 9. Glossary

- **Static site**: just files, no code running on a server.
- **DOM**: the page as a tree of elements that JavaScript can find and change.
- **CDN**: public server for common libraries.
- **localhost / port**: your own machine / which "door" a server listens on.
- **404 / 200**: "not found" / "OK".
- **Deploy**: publish a new version of the site.
