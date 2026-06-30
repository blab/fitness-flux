# Fitness flux — manuscript

`press` (not yet publicly released) renders the Markdown-ish manuscript into a themed, self-contained HTML paper with inline interactive figures, math, and citations.

- Entry: `fitness-flux.md`
- Config: `press.config.yaml` (entry, `figures_dir: ../viz`, `bib: fitness_flux.bib`, `theme: ./theme`)
- Bibliography: `fitness_flux.bib`
- Theme: `theme/tokens.css` (blotter look — Museo via Adobe Typekit + warm grayscale chrome)
- Interactive figure components live in `../viz/<component>/`

## Author locally (live preview)

From this directory:

```sh
press preview .
```

Open http://localhost:8000. The page hot-reloads on edits to the Markdown, the `.bib`, the theme, and the live `../viz` components and their data — the fast edit/refresh loop in place of the LaTeX PDF cycle.

## Build a static bundle

```sh
press build .
```

Writes a self-contained bundle to `dist/` (gitignored). Open `dist/index.html` to check it. The bundle freezes each interactive figure, copies static figures, vendors the KaTeX assets, and includes a `.nojekyll` file so GitHub Pages serves it as-is.

> Note: open the built page over **http** (e.g. `python3 -m http.server -d dist 8099`), not by double-clicking the file — browsers block the figures' module imports/fetches over `file://`.

## PDF

For quick PDF, `press pdf` outputs HTML to PDF. Print/PDF can't use the live interactive figures, so press renders each component to a high-DPI static PNG (committed under `figures/`) and uses those for print.

```sh
press figures .  # render each labelled component → figures/<slug>.png (re-run when a component or its data changes)
press pdf .      # build + headless-print → dist/fitness-flux.pdf (US letter, 0.65in margins)
```

`press figures` writes `figures/<slug>.png` (the `#fig:<slug>` label) as each figure's automatic static fallback, and skips any figure with an explicit `static=` override. `figures`/`pdf` use headless Chromium — install it once with `npx playwright install chromium`.

## LaTeX

For a journal-grade PDF, `press tex` emits a LaTeX source that `\includegraphics` the same `figures/<slug>.png` images and cites `fitness_flux.bib` through `plos.bst`, so references come out in PLOS style.

```sh
press figures .   # only if components or data changed since the last run
press tex .       # → fitness-flux.tex
pdflatex fitness-flux && bibtex fitness-flux && pdflatex fitness-flux && pdflatex fitness-flux
```

The generated `fitness-flux.tex` is a clean PLOS-style article (authblk author/affiliation block, abstract, sectioned body, math, `\cite`/`\ref`, figures via `\includegraphics`) — the natural source-submission artifact. The pdflatex byproducts (`.aux`, `.bbl`, `.log`, `.out`, `.pdf`) are build artifacts to gitignore. Requires a LaTeX install (`pdflatex`, `bibtex`).

## Deploy to GitHub Pages

The paper is published from a `gh-pages` branch that holds **only the generated site** (force-pushed as a single commit; `main` is never touched). The live URL is https://blab.github.io/fitness-flux/.

`publish.sh` does the build + publish in one step:

```sh
./publish.sh
```
