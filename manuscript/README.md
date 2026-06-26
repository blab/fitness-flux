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

## Deploy to GitHub Pages

The paper is published from a `gh-pages` branch that holds **only the generated site** (force-pushed as a single commit; `main` is never touched). The live URL is https://blab.github.io/fitness-flux/.

`publish.sh` does the build + publish in one step:

```sh
./publish.sh
```
