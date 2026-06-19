# Portable Visualization Components — Implementation Spec

Trevor Bedford - 2026-06-19

## Purpose

Define a convention for authoring data visualizations as **portable components** that render unchanged across three contexts:

1. **Dashboard** — a list-of-figures page (the current `fitness-flux.html` shape).
2. **Living paper** — a long, responsive, scroll-y page with figures inline among prose, equations, and citations (the dashboard→blog→manuscript continuum).
3. **Slides** — Reveal.js, generally one interactive figure per slide.

The component — its code and its data — is written once and consumed by all three. This document is the contract for how a component is structured and how each context consumes it.

**In scope:** the figure component contract, the component directory layout (the `viz/<component>/` bundle), the three consumption patterns, data versioning, the freeze-for-talks boundary, and the recipe for extracting a component out of an existing monolithic dashboard. **Host *locations* are out of scope** — where the dashboard, the living paper (a top-level `paper/`), or slide decks live in the repo is the host's concern; they sit outside `/viz/` and consume the atomic `viz/<component>/` units.

**Explicitly deferred (not this spec):** the living-paper authoring *framework*. The figure components defined here are framework-agnostic plain ES modules, so a future living-paper host (MyST is the current front-runner; Quarto is the fallback) can consume them via a raw-HTML or OJS block without rework. Do not adopt a manuscript framework as part of implementing this spec. See "Future host" at the end.

## Core concept: component, dataset, figure

Three nouns, kept distinct:

- A **component** is a rendering capability — a pure function over a data *shape* (e.g. `fitness-flux`: draw frequency-weighted fitness violins for any dataset shaped `{ fitnessFlux, colors }`).
- A **dataset** is a concrete instance conforming to that shape (e.g. `sarscov2_clades`, `sarscov2_lineages`, `h3n2_clades`). A component is **one-to-many** with datasets.
- A **figure** is the binding of `(component, dataset, opts)` — one component rendering one dataset with one set of options. A living paper showing the same component for two datasets contains two figures.

This is the correction to an earlier draft that assumed one dataset per component. The component holds its conformant datasets in a `data/` subdirectory; the binding to a single dataset happens at the host (living/dashboard) or is frozen in (talks).

## Two lifecycles

A figure component is a small, self-contained unit of code plus its versioned data. It has two lifecycles, and keeping them distinct is the key idea:

- **Living.** The dashboard and the living paper track the current version of a component (repo `main`). When the analysis workflow reruns and the data is regenerated and committed, these update. Rollback is git (check out a tag or SHA).
- **Frozen.** A talk takes a *copy* of the component and one dataset at talk-time into a dated, self-contained directory. A talk is an artifact, not a living system — the 2026 version of a figure must keep rendering exactly as it did in 2026, offline, years later.

The same component serves both. Sharing within the living library is fine (it is one versioned thing); freezing happens only at the talk boundary.

## The render contract

A figure is an ES module exporting a single pure function:

```js
/**
 * Render a figure into a container.
 *
 * PURE: render() never fetches. One dataset's already-loaded, parsed object
 * arrives via `data`. This is what decouples the figure from data loading,
 * from the dataset selector, and from its host context.
 *
 * @param {HTMLElement} container  Element to render into (cleared first).
 * @param {Object}      data       One dataset, pre-loaded and parsed. Shape is
 *                                  component-specific and documented in this
 *                                  module's header (see below). The dataset is
 *                                  self-contained: it embeds its own color table.
 * @param {Object}      [opts]
 * @param {"inline"|"slide"|"dashboard"} [opts.mode="inline"]
 *                                  Styling/sizing profile. inline = responsive
 *                                  to container width, modest fonts; slide =
 *                                  larger fonts, sized to viewport; dashboard =
 *                                  as inline but expects sibling figures.
 * @param {number}      [opts.width]   Explicit width; if omitted, measured from
 *                                     container.clientWidth.
 * @param {number}      [opts.height]  Explicit height; if omitted, mode default.
 * @param {Object}      [opts...]      Component-specific options (e.g. which
 *                                     variants to highlight) pass through here.
 * @returns {{
 *   element: HTMLElement,          // the rendered node
 *   resize: (width?, height?) => void,  // re-render at a new size
 *   destroy: () => void            // remove listeners/observers, clear container
 * }}
 */
export function render(container, data, opts = {}) { /* ... */ }
```

Rules the contract enforces:

- **No fetching inside `render`.** The host loads one dataset and passes it in. This is non-negotiable — it is what makes the same function work in a dashboard (selector-driven dataset), a living paper (one figure per dataset), and a frozen slide (snapshot dataset).
- **No `ResizeObserver` inside `render` by default.** `render` does a one-time render at the given or measured size and exposes `resize()`. The *host* decides when to call it (see contexts). This avoids the classic Reveal.js failure where a figure measures zero width while its slide is hidden.
- **Return a controller, not nothing.** `destroy()` must remove every listener/observer it added, so a living-paper page with many figures doesn't leak.
- **The dataset is self-contained.** Each dataset object embeds everything the component needs to render it, including its own `colors` table — see "Colors" below. There is no separate shared color fetch.
- **Each module documents its `data` shape** in a header comment, e.g.:

```js
// fitness-flux/fitness-flux.js
//
// data shape (one dataset):
//   {
//     fitnessFlux: Array<{ variant, date (ISO string), fit (number), freq (number) }>,
//     colors:      Array<{ variant, color, display_name, is_major, order }>
//   }
//
// Adapted from <repo>@<sha>:fitness-flux.html (renderFitnessFlux) on 2026-06-19
```

The provenance comment is required on any component extracted from a dashboard (see Conventions).

## Colors

Colors are **per-dataset, embedded in each dataset file** — not a shared cross-component table. Different datasets have different variant universes (h3n2 clades and sarscov2 clades share no variants; `sarscov2_clades` even sources its colors from the Nextstrain tree while the others use a gradient), so a single shared color file breaks down the moment datasets span organisms. Each dataset JSON therefore carries its own `colors` array alongside its data.

`lib/colors.js` remains — but only as the generic *builder* (`colorScale(data.colors)`, `buildLegend(...)`), shared across components because the logic is identical. The color *data* travels inside each dataset. There is no `shared/colors.json`.

## Directory layout

The component library lives under `/viz/`. **Each component is a self-contained bundle directory directly under `/viz/`** whose kebab-case name carries identity; the module repeats that name, and each component holds its conformant datasets in a `data/` subdirectory:

```
/viz/
  fitness-flux/                  # one component = one directory (the portable unit)
    fitness-flux.js              # pure render() — unit of truth (same name as dir)
    index.html                   # dev/preview harness (standalone view + iframe target)
    meta.json                    # optional: dataset ids, labels, default (see Manifest)
    data/                        # this component's conformant datasets
      sarscov2_clades.json       #   each self-contained: { fitnessFlux, colors }
      sarscov2_lineages.json
      h3n2_clades.json
  frequency-panels/
    frequency-panels.js
    index.html
    meta.json
    data/
      sarscov2_clades.json
      sarscov2_lineages.json
      h3n2_clades.json
  lib/
    plot.js                      # vendored Observable Plot (ESM)
    d3.js                        # vendored d3 (ESM)
    colors.js                    # generic color-scale / legend builder (consumes data.colors)
    data.js                      # helper: fetchJSON / parseTSV etc. (host-side use)
  shared/                        # reserved for anything genuinely cross-component (may be empty)
```

Hosts — the dashboard, the living paper, slide decks — live *outside* this tree and consume these atomic units; their locations are out of scope for this spec (see Purpose).

The component directory is the portable unit: copy, move, archive, or drop it into another repo as one thing — code and all its datasets travel together. `index.html` as the harness name means the directory URL resolves to it with no filename, so embeds are clean (`src="/viz/fitness-flux/"`).

Datasets live *inside* the component (`viz/{component}/data/{dataset}.json`) rather than in a top-level `viz/data/{component}/` tree, because a dataset is shaped for one component's render function — it is component-specific, not cross-component — so it belongs in the bundle. Splitting it out would scatter one component across two trees and break the copy-and-it-works property. The only genuinely shared things are the `lib/` helpers; `shared/` is reserved for cross-component data should any ever arise.

`lib/` is shared *within* this living library. A component reaches up — `../lib/plot.js`. At the talk-freeze boundary, the referenced `lib/` files are copied *into* the frozen bundle so the frozen artifact is self-contained — see "Freeze boundary." This is the "living shares, frozen self-contains" split: the dev-time bundle is not fully self-contained, but the frozen one is, which is the end that matters.

Use `data/{dataset}.json` uniformly, even for a component that currently has a single dataset, so the convention is predictable and adding a second dataset later needs no restructure.

## Manifest (`meta.json`)

A small optional per-component manifest enumerating its datasets, so hosts don't hardcode or rely on directory globbing (static hosts can't list a directory):

```json
{
  "datasets": [
    { "id": "sarscov2_clades",   "label": "SARS-CoV-2 clades" },
    { "id": "sarscov2_lineages", "label": "SARS-CoV-2 lineages" },
    { "id": "h3n2_clades",       "label": "H3N2 clades" }
  ],
  "default": "sarscov2_clades"
}
```

Ideally emitted by the analysis workflow, which already knows the dataset list, so it stays in sync. Consumed by the dashboard selector, the dev harness, and (later) living-paper tooling to list datasets with human labels and a defined order. Optional at three datasets; recommended as the library grows.

## Anatomy of a component (one bundle directory)

**1. The module** (`fitness-flux/fitness-flux.js`) — exports the pure `render` per the contract. Imports Plot/d3/helpers from `../lib/`. Contains only rendering logic; no fetch, no selector, no DOM outside `container`.

**2. The dev harness** (`fitness-flux/index.html`) — ~20 lines. Loads one dataset, imports `render`, calls it. Doubles as (a) the independently-viewable preview page and (b) the iframe target. Reads `?dataset=<id>` to preview any dataset, falling back to the `meta.json` default. Owns its own `ResizeObserver` so it reflows when embedded.

```html
<!doctype html>
<meta charset="utf-8">
<style>html,body{margin:0} #fig{width:100vw;height:100vh}</style>
<div id="fig"></div>
<script type="module">
  import { render } from "./fitness-flux.js";
  const meta = await fetch("./meta.json").then(r => r.json()).catch(() => null);
  const id = new URLSearchParams(location.search).get("dataset")
           ?? meta?.default ?? "sarscov2_clades";
  const data = await fetch(`./data/${id}.json`).then(r => r.json());
  const fig = render(document.getElementById("fig"), data, { mode: "slide" });
  new ResizeObserver(() => fig.resize()).observe(document.getElementById("fig"));
</script>
```

**3. The datasets** (`fitness-flux/data/*.json`) — one self-contained JSON per dataset, each carrying both the figure's data and its embedded `colors` table, emitted by the analysis workflow and committed. JSON (not TSV) so hosts parse with `fetch().json()` and no custom TSV parser is needed at render time; the workflow does the TSV→JSON conversion at generation time.

## The three consumption contexts

**Dashboard.** Imports each module and calls `render` into its slot. A dataset selector is a *host* affordance: it reads the dataset list from `meta.json`, and on change loads that dataset's JSON and re-renders (`controller.destroy()` then `render()` again, or re-`render` into the cleared container). The selector is not part of any component.

**Living paper (inline).** Imports a module once and instantiates it **once per dataset** — each its own `<figure>` — directly, with **no iframes**, so figures flow with the prose column, tooltips overflow correctly, and width is responsive. The host owns a `ResizeObserver` per figure.

```html
<p>…prose about SARS-CoV-2…</p>
<figure id="ff-sarscov2" class="figure"></figure>
<p>…prose about H3N2…</p>
<figure id="ff-h3n2" class="figure"></figure>

<script type="module">
  import { render } from "../viz/fitness-flux/fitness-flux.js";
  for (const [slot, id] of [["ff-sarscov2","sarscov2_clades"], ["ff-h3n2","h3n2_clades"]]) {
    const data = await fetch(`../viz/fitness-flux/data/${id}.json`).then(r => r.json());
    const el = document.getElementById(slot);
    const c = render(el, data, { mode: "inline" });
    new ResizeObserver(() => c.resize()).observe(el);
  }
</script>
```

Same component (`fitness-flux.js`), two datasets, two figures — exactly the one-to-many model.

**Slide (Reveal.js).** Embeds the *frozen* bundle via iframe, lazy-loaded on slide entry using Reveal's `data-src` (which defers the load until the slide is reached, sidestepping the hidden-slide zero-width problem). The frozen directory is self-contained, so the deck works offline on stage.

```html
<section>
  <iframe data-src="/viz/fitness-flux-sarscov2-clades-2026-06-18/"
          width="100%" height="100%" frameborder="0" scrolling="no"></iframe>
</section>
```

The slide points at a **frozen, dated copy** whose directory name encodes the `(component, dataset, date)` binding — not at the living `viz/fitness-flux/` path. The directory URL resolves to its `index.html` with no filename. The living paper points at the living path. Same component, two lifecycles, as designed.

## Freeze boundary (talk snapshots)

Freezing binds one `(component, dataset, opts)` into a self-contained, dated directory. The chosen dataset collapses to a single `data.json`; the directory name carries the identity:

1. Create `/viz/<component>-<dataset>-<YYYY-MM-DD>/` (in the talks site's static pool; kebab-case the dataset id for the URL).
2. Copy in the module (`<component>.js`), copy the chosen dataset to `data.json`, and copy the `lib/` files the module imports (`plot.js`, `d3.js`, `colors.js`, …), flattening everything into the dated directory.
3. Write a minimal frozen `index.html` that fetches `./data.json` unconditionally (no `?dataset=`, no `meta.json`) and imports `./<component>.js`. Rewrite the module's `../lib/...` imports to `./...`. The frozen dir must have **no external references** — no CDN, no `../` escapes.
4. Verify it renders opened directly from disk with the network disabled (offline test).

Resulting frozen layout (flat and self-contained, one dataset only):

```
/viz/fitness-flux-sarscov2-clades-2026-06-18/
  index.html          # fetches ./data.json, imports ./fitness-flux.js
  fitness-flux.js
  data.json           # the chosen dataset, copied + renamed (colors embedded)
  plot.js
  d3.js
  colors.js
```

This is the artifact a slide iframes. It never updates. A later talk — or the same talk with a different dataset — freezes a new dated directory.

## Data and versioning workflow

- Each dataset (`viz/{component}/data/{dataset}.json`) and the optional `meta.json` are **generated, not hand-edited.** A re-runnable workflow in the analysis repo (Snakemake/Make) produces them and they are committed to the repo. Git is the version history, the provenance, and the rollback mechanism. (The fitness-flux workflow already emits to `viz/{component}/data/{analysis}.json` per dataset.)
- The living paper and dashboard, served from `main`, always show current data. To present or cite an earlier state, check out a tag/SHA or deploy a pinned ref.
- Served same-origin from Cloudflare/GitHub Pages with `lib/` vendored, there is **no external fetch at render time** — robust to network conditions and durable in the historical record.
- Talks freeze one dataset into their dated directory (Freeze boundary above), so a deck is hermetic regardless of later workflow reruns.

## Static tiers (don't beat PNGs everywhere)

Not every figure should be interactive. Decide per figure:

- **Tier 1 — SVG (or PNG).** Static figures where hover/zoom add nothing. The `render` module can serialize its Plot/d3 output to a standalone `.svg` (vector, frozen, dependency-free — crisp on any projector, strictly better than a PNG for line/scatter). Use this for most figures.
- **Tier 2 — frozen interactive page (iframe).** Figures where hovering genuinely informs (e.g. the fitness-flux violin: read variant/date/fitness/frequency on hover; the frequency panels: per-point tooltips). Vendored + embedded data makes these nearly as bulletproof as a PNG, but live.
- **Tier 3 — choreographed (deferred).** A figure that must react to slide state (fragment-driven reveals, transitions on advance). Not built yet; would use direct module import into the slide rather than an iframe.

The decision rule: does interactivity or vector-crispness justify the live-viz tax for *this* figure? If not, Tier 1.

## Extraction recipe (monolith → component)

To carve a figure out of an existing dashboard like `fitness-flux.html`:

1. Identify the render function and its transitive dependencies (its data inputs, the color layer, any constants like `FLUX_HALF_HEIGHT`, and helpers it calls).
2. Create the bundle directory `viz/<component>/` and move the render logic into `<component>.js` (same name as the directory) as the pure `render(container, data, opts)`. Lift all `fetch`/parse/selector/`status` machinery *out* — it becomes host responsibility.
3. **Freeze the transformed output where possible.** For fitness-flux, the render joins frequencies and fitness into a single array; emit *that already-joined array* (plus the embedded `colors`) as each dataset's JSON rather than raw inputs, so the component carries no join logic and resize works from the passed-in data.
4. Generate one self-contained JSON **per dataset** to `viz/<component>/data/{dataset}.json` via the analysis workflow, each embedding its own `colors`; optionally emit `meta.json`.
5. Write the dev harness (`index.html`) with `?dataset=` + default.
6. Add the provenance header comment.
7. Verify standalone (each dataset via `?dataset=`), inline (responsive), and iframe (resizes) all render.

The current `fitness-flux.html` is already cleanly factored by render function (`renderFitnessFlux`, `renderFrequencyPanels`, `renderMutations`, `renderFlux`), so extraction is mechanical, not archaeological.

## Conventions

- **No build step.** Plain ES modules, vendored libraries, JSON data. No bundler, no transpile. Everything must work served as static files and (for frozen dirs) opened from disk.
- **No browser storage.** No `localStorage`/`sessionStorage`; keep state in memory.
- **Vendor Plot and d3** into `lib/` from the start (don't depend on a CDN at render time). Dev may import from a CDN for convenience, but the committed components and all frozen dirs use vendored copies.
- **Provenance comment** at the top of every extracted module: `// Adapted from <repo>@<sha>:<source> (<function>) on <YYYY-MM-DD>`. Pure metadata; creates no dependency.
- **Data is JSON, generated, per-dataset, colors embedded.** No hand-edited data files; no TSV parsing at render time; no shared color file.
- **Naming:** a component is a directory directly under `viz/` named in kebab-case (`fitness-flux/`); the directory name carries identity and is what appears in URLs. The module repeats the directory name (`fitness-flux.js`) — kebab-case throughout, matching the prevailing JS/web file convention, so directory and file echo each other with no transformation rule. The harness is always `index.html`; datasets are always `data/{dataset}.json`. Dataset ids follow the analysis pipeline's naming (snake_case, e.g. `sarscov2_clades`). Frozen dirs encode the binding as `<component>-<dataset>-<date>/` (kebab-case the dataset id for the URL) and collapse the dataset to `data.json`.

## First milestone & acceptance criteria

Build the `fitness-flux` component end to end across its three datasets. Done when:

- [ ] `fitness-flux/fitness-flux.js` exports a pure `render(container, data, opts)` per the contract; no fetch, selector, or status logic inside.
- [ ] The workflow generates `fitness-flux/data/{sarscov2_clades,sarscov2_lineages,h3n2_clades}.json`, each self-contained with embedded `colors`; `meta.json` optionally emitted.
- [ ] `fitness-flux/index.html` renders any dataset via `?dataset=<id>` (default from `meta.json`), standalone and inside an iframe (`mode: "slide"`, reflows).
- [ ] Rendering the same component for two datasets as two inline figures between prose, each reflowing on container resize (`mode: "inline"`), works — a verification step, not a committed artifact.
- [ ] A freeze of one binding into `/viz/fitness-flux-<dataset>-<date>/` collapses the dataset to `data.json`, is self-contained (vendored libs, dir-local paths, no `../` or CDN), and renders from disk with the network disabled.
- [ ] Provenance comment present; no build step; no browser storage.

## Open decisions (resolve during implementation)

- **Repo placement of `/viz/`.** Recommended: live in the analysis/dashboard repo initially, with the talks site consuming frozen copies. Revisit if a central viz library proves better.

Resolved: components are atomic directories directly under `/viz/` (`viz/<component>/`, no `figures/` level); datasets live inside the component as `data/{dataset}.json` (not a top-level `viz/data/` tree); colors are per-dataset and embedded (no `shared/colors.json`); the figure receives one dataset object via `data`, with its colors inside it. `meta.json` is workflow-emitted, with `default` = `sarscov2_clades`. Host locations are out of scope — the dashboard is retained for now as `viz/fitness-flux.html`, and the living paper will be a top-level `paper/`.

## Future host (deferred — not this spec)

When the manuscript end matures (numbered figures/equations, citations, cross-references, PDF + HTML from one source), evaluate a host framework — MyST first, Quarto as the proven fallback — to author the living paper. Because components here are framework-agnostic ES modules with a pure `render`, that host consumes them via a raw-HTML or OJS block with no change to the components. The figure components are the durable investment; the host is swappable. Do not let the eventual host choice block or reshape this component work.
