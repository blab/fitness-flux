# Portable Visualization Components — Implementation Spec

Trevor Bedford - 2026-06-19

## Purpose

Define a convention for authoring data visualizations as **portable components** that render unchanged across three contexts:

1. **Dashboard** — a list-of-figures page (the current `fitness-flux.html` shape).
2. **Living paper** — a long, responsive, scroll-y page with figures inline among prose, equations, and citations (the dashboard→blog→manuscript continuum).
3. **Slides** — Reveal.js, generally one interactive figure per slide.

The component — its code and its data — is written once and consumed by all three. This document is the contract for how a component is structured and how each context consumes it.

**In scope:** the figure component contract, directory layout, the three consumption patterns, data versioning, the freeze-for-talks boundary, and the recipe for extracting a component out of an existing monolithic dashboard.

**Explicitly deferred (not this spec):** the living-paper authoring *framework*. The figure components defined here are framework-agnostic plain ES modules, so a future living-paper host can consume them via a raw-HTML or OJS block without rework. Do not adopt a manuscript framework as part of implementing this spec. See "Future host" at the end.

## Core concept: one component, two lifecycles

A figure component is a small, self-contained unit of code plus its versioned data. It has two lifecycles, and keeping them distinct is the key idea:

- **Living.** The dashboard and the living paper track the current version of a component (repo `main`). When the analysis workflow reruns and the data file is regenerated and committed, these update. Rollback is git (check out a tag or SHA).
- **Frozen.** A talk takes a *copy* of the component and its data at talk-time into a dated, self-contained directory. A talk is an artifact, not a living system — the 2026 version of a figure must keep rendering exactly as it did in 2026, offline, years later.

The same component serves both. Sharing within the living artifacts is fine (they are one versioned thing); freezing happens only at the talk boundary.

## The render contract

A figure is an ES module exporting a single pure function:

```js
/**
 * Render a figure into a container.
 *
 * PURE: render() never fetches. All data arrives pre-loaded and parsed via
 * the `data` argument. This is what decouples the figure from data loading,
 * from the dataset selector, and from its host context.
 *
 * @param {HTMLElement} container  Element to render into (cleared first).
 * @param {Object}      data       Pre-loaded, parsed inputs this figure needs.
 *                                  Shape is figure-specific and documented in
 *                                  this module's header (see below).
 * @param {Object}      [opts]
 * @param {"inline"|"slide"|"dashboard"} [opts.mode="inline"]
 *                                  Styling/sizing profile. inline = responsive
 *                                  to container width, modest fonts; slide =
 *                                  larger fonts, sized to viewport; dashboard =
 *                                  as inline but expects sibling figures.
 * @param {number}      [opts.width]   Explicit width; if omitted, measured from
 *                                     container.clientWidth.
 * @param {number}      [opts.height]  Explicit height; if omitted, mode default.
 * @param {Object}      [opts...]      Figure-specific options (e.g. which
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

- **No fetching inside `render`.** The host loads data and passes it in. This is non-negotiable — it is what makes the same function work in a dashboard (selector-driven data), a living paper (canonical data), and a frozen slide (snapshot data).
- **No `ResizeObserver` inside `render` by default.** `render` does a one-time render at the given or measured size and exposes `resize()`. The *host* decides when to call it (see contexts). This avoids the classic Reveal.js failure where a figure measures zero width while its slide is hidden.
- **Return a controller, not nothing.** `destroy()` must remove every listener/observer it added, so a living-paper page with many figures doesn't leak.
- **Each module documents its `data` shape** in a header comment, e.g.:

```js
// figures/fitness-flux/fitness-flux.js
//
// data shape:
//   {
//     fitnessFlux: Array<{ variant, date (Date), fit (number), freq (number) }>,
//     colors:      Array<{ variant, color, display_name, order, is_major }>  // from shared/colors.json
//   }
//
// Adapted from <repo>@<sha>:fitness-flux.html (renderFitnessFlux) on 2026-06-19
```

The provenance comment is required on any component extracted from a dashboard (see Conventions).

## Directory layout

The component library lives under `/viz/`. **Each figure is a self-contained bundle directory** whose kebab-case name carries the figure's identity; the module inside repeats that name (so directory and file echo each other), and the wrapper and data use stable generic names:

```
/viz/
  figures/
    fitness-flux/            # one figure = one directory (the portable unit)
      fitness-flux.js        # pure render() — unit of truth (same name as dir)
      index.html             # thin wrapper: standalone view + iframe target
      data.json              # this figure's data, workflow-regenerated
    frequency-panels/
      frequency-panels.js
      index.html
      data.json
  shared/
    colors.json              # cross-figure color + display-name table
  lib/
    plot.js                  # vendored Observable Plot (ESM)
    d3.js                    # vendored d3 (ESM)
    colors.js                # helper: build color scale / display-name map from colors.json
    data.js                  # helper: fetchJSON / parseTSV etc. (host-side use)
  index.html                 # (later) the living paper; out of scope for first milestone
```

The bundle directory is the portable unit: copy, move, archive, or drop it into another repo as one thing. `index.html` as the wrapper name means the directory URL resolves to it with no filename, so embeds are clean (`src="/viz/figures/fitness-flux/"`).

`shared/` and `lib/` are shared *within* this living library (the dashboard and living paper are one versioned artifact, so sharing is correct here). A figure that needs them reaches up — `../../shared/colors.json`, `../../lib/plot.js`. At the talk-freeze boundary, the referenced `shared/` and `lib/` files are copied *into* the frozen bundle so the frozen artifact is self-contained — see "Freeze boundary." This is the "living shares, frozen self-contains" split: the dev-time bundle is not fully self-contained, but the frozen one is, which is the end that matters.

## Anatomy of a figure (one bundle directory)

**1. The module** (`figures/fitness-flux/fitness-flux.js`) — exports the pure `render` per the contract. Imports Plot/d3/helpers from `../../lib/`. Contains only rendering logic; no fetch, no selector, no DOM outside `container`.

**2. The wrapper** (`figures/fitness-flux/index.html`) — ~15 lines. Loads this figure's data, imports `render`, calls it. Doubles as (a) the independently-viewable page and (b) the iframe target for slides. Owns its own `ResizeObserver` so it reflows when embedded.

```html
<!doctype html>
<meta charset="utf-8">
<style>html,body{margin:0} #fig{width:100vw;height:100vh}</style>
<div id="fig"></div>
<script type="module">
  import { render } from "./fitness-flux.js";
  const [fitnessFlux, colors] = await Promise.all([
    fetch("./data.json").then(r => r.json()),
    fetch("../../shared/colors.json").then(r => r.json()),
  ]);
  const c = render(document.getElementById("fig"), { fitnessFlux, colors }, { mode: "slide" });
  new ResizeObserver(() => c.resize()).observe(document.getElementById("fig"));
</script>
```

**3. The data** (`figures/fitness-flux/data.json`, plus `shared/colors.json`) — the figure's own data lives in the bundle; the cross-figure color table lives in `shared/`. Both are emitted by the analysis workflow and committed. JSON (not TSV) so the wrapper and living paper parse it with `JSON.parse`/`fetch().json()` and no custom TSV parser is needed at render time. The analysis workflow does the TSV→JSON conversion at generation time.

## The three consumption contexts

**Dashboard.** Imports each module and calls `render` into its slot. A dataset selector, if retained, is a *host* affordance: on change, the host loads the other dataset's JSON and re-renders (`controller.destroy()` then `render()` again, or just re-`render` into the cleared container). The selector is not part of any component.

**Living paper (inline).** Imports modules directly — **no iframes** — so figures flow with the prose column, tooltips overflow correctly, and width is responsive. The host owns a `ResizeObserver` per figure.

```html
<p>…prose introducing the figure…</p>
<figure id="fitness-flux" class="figure"></figure>
<p>…prose continuing the argument…</p>

<script type="module">
  import { render } from "./figures/fitness-flux/fitness-flux.js";
  const [fitnessFlux, colors] = await Promise.all([
    fetch("figures/fitness-flux/data.json").then(r => r.json()),
    fetch("shared/colors.json").then(r => r.json()),
  ]);
  const el = document.getElementById("fitness-flux");
  const c = render(el, { fitnessFlux, colors }, { mode: "inline" });
  new ResizeObserver(() => c.resize()).observe(el);
</script>
```

**Slide (Reveal.js).** Embeds the *frozen* wrapper page via iframe, lazy-loaded on slide entry using Reveal's `data-src` (which defers the load until the slide is reached, sidestepping the hidden-slide zero-width problem). The frozen directory is self-contained, so the deck works offline on stage.

```html
<section>
  <iframe data-src="/viz/fitness-flux-2026-06-18/"
          width="100%" height="100%" frameborder="0" scrolling="no"></iframe>
</section>
```

Note the slide points at a **frozen, dated copy** (`/viz/fitness-flux-2026-06-18/`), not at the living `figures/fitness-flux/` path. The directory URL resolves to its `index.html` with no filename. The living paper points at the living path. Same component, two lifecycles, as designed.

## Freeze boundary (talk snapshots)

Freezing a figure for a talk produces a self-contained, dated directory. Because the source is already a bundle, the freeze is close to a directory copy:

1. Copy the figure bundle to `/viz/<figure>-<YYYY-MM-DD>/` (in the talks site's static pool).
2. Copy in the `shared/` and `lib/` files the bundle references (`colors.json`, `plot.js`, `d3.js`, `colors.js`, …), flattening them into the dated directory.
3. Rewrite the now-local references so all imports and fetches are directory-local: the wrapper's `../../shared/colors.json` becomes `./colors.json`, the module's `../../lib/plot.js` becomes `./plot.js`. The frozen dir must have **no external references** — no CDN, no `../` escapes.
4. Verify it renders opened directly from disk (offline test).

Resulting frozen layout (flat and self-contained):

```
/viz/fitness-flux-2026-06-18/
  index.html          # wrapper, fetches ./data.json + ./colors.json, imports ./fitness-flux.js
  fitness-flux.js
  data.json
  colors.json
  plot.js
  d3.js
  colors.js
```

This is the artifact a slide iframes. It never updates. A later talk wanting a newer version freezes a new dated directory.

## Data and versioning workflow

- Each figure's `data.json` (and `shared/colors.json`) is **generated, not hand-edited.** A re-runnable workflow in the analysis repo (Snakemake/Make) produces the JSON and it is committed to the repo. Git is the version history, the provenance, and the rollback mechanism.
- The living paper and dashboard, served from `main`, always show current data. To present or cite an earlier state, check out a tag/SHA or deploy a pinned ref.
- Served same-origin from Cloudflare/GitHub Pages with `lib/` vendored, there is **no external fetch at render time** — robust to network conditions and durable in the historical record.
- Talks freeze the data into their dated directory (step 2 above), so a deck is hermetic regardless of later workflow reruns.

## Static tiers (don't beat PNGs everywhere)

Not every figure should be interactive. Decide per figure:

- **Tier 1 — SVG (or PNG).** Static figures where hover/zoom add nothing. The `render` module can serialize its Plot/d3 output to a standalone `.svg` (vector, frozen, dependency-free — crisp on any projector, strictly better than a PNG for line/scatter). Use this for most figures.
- **Tier 2 — frozen interactive page (iframe).** Figures where hovering genuinely informs (e.g. the fitness-flux violin: read variant/date/fitness/frequency on hover; the frequency panels: per-point tooltips). Vendored + inlined makes these nearly as bulletproof as a PNG, but live.
- **Tier 3 — choreographed (deferred).** A figure that must react to slide state (fragment-driven reveals, transitions on advance). Not built yet; would use direct module import into the slide rather than an iframe.

The decision rule: does interactivity or vector-crispness justify the live-viz tax for *this* figure? If not, Tier 1.

## Extraction recipe (monolith → component)

To carve a figure out of an existing dashboard like `fitness-flux.html`:

1. Identify the render function and its transitive dependencies (its data inputs, the shared color layer, any constants like `FLUX_HALF_HEIGHT`, and helpers it calls).
2. Create the bundle directory `figures/<figure-name>/` and move the render logic into `<figure-name>.js` (same name as the directory) as the pure `render(container, data, opts)`. Lift all `fetch`/parse/selector/`status` machinery *out* — it becomes host responsibility.
3. **Freeze the transformed output where possible.** For fitness-flux, the render currently joins `frequencies` and `mutation_fitness` into a `fitnessFluxData` array; emit *that already-joined array* as `data.json` rather than the two raw inputs, so the component carries no join logic and resize still works from the cached array.
4. Generate `data.json` from the chosen dataset via the analysis workflow; write the color table to `shared/colors.json`.
5. Write the thin wrapper (`index.html`).
6. Add the provenance header comment.
7. Verify standalone, inline (responsive), and iframe (resizes) all render.

The current `fitness-flux.html` is already cleanly factored by render function (`renderFitnessFlux`, `renderFrequencyPanels`, `renderMutations`, `renderFlux`), so extraction is mechanical, not archaeological.

## Conventions

- **No build step.** Plain ES modules, vendored libraries, JSON data. No bundler, no transpile. Everything must work served as static files and (for frozen dirs) opened from disk.
- **No browser storage.** No `localStorage`/`sessionStorage`; keep state in memory.
- **Vendor Plot and d3** into `lib/` from the start (don't depend on a CDN at render time). Dev may import from a CDN for convenience, but the committed components and all frozen dirs use vendored copies.
- **Provenance comment** at the top of every extracted module: `// Adapted from <repo>@<sha>:<source> (<function>) on <YYYY-MM-DD>`. Pure metadata; creates no dependency.
- **Data is JSON, generated.** No hand-edited data files; no TSV parsing at render time.
- **Naming:** each figure is a directory under `figures/` named in kebab-case (`fitness-flux/`); the directory name carries identity and is what appears in URLs. The module inside **repeats the directory name** (`fitness-flux.js`) — kebab-case throughout, matching the prevailing JS/web file convention, so directory and file echo each other with no transformation rule to remember. The wrapper is always `index.html` (so the directory URL resolves with no filename) and the figure's data is always `data.json`. Frozen dirs use the kebab name + date (`fitness-flux-2026-06-18/`) and keep the matching module (`fitness-flux.js`). Kebab everywhere is deliberate: it's the ecosystem-native convention for module files and avoids any case-only filename collisions on case-insensitive filesystems.

## First milestone & acceptance criteria

Build exactly one figure — `fitness-flux` — end to end. Done when:

- [ ] `figures/fitness-flux/fitness-flux.js` exports a pure `render(container, data, opts)` per the contract; no fetch, selector, or status logic inside.
- [ ] `figures/fitness-flux/data.json` and `shared/colors.json` are generated extracts of the chosen dataset; `data.json` holds the already-joined flux array.
- [ ] `figures/fitness-flux/index.html` renders the figure when opened standalone and reflows when embedded in an iframe (`mode: "slide"`).
- [ ] A minimal inline test page renders the figure between two prose paragraphs and reflows on container resize (`mode: "inline"`).
- [ ] A freeze into `/viz/fitness-flux-<date>/` is self-contained (vendored libs, dir-local paths, no `../` or CDN) and renders correctly opened from disk with the network disabled.
- [ ] Provenance comment present; no build step; no browser storage.

## Open decisions (resolve during implementation)

- **Repo placement of `/viz/`.** Recommended: live in the analysis/dashboard repo next to the source dashboard initially, with the talks site consuming frozen copies. Revisit if a central viz library proves better.
- **Colors as `data.colors` vs `opts.colors`.** The color table now has a home (`shared/colors.json`), but the open question remains whether the figure receives it as a data input (`data.colors`, current spec, uniform contract) or as an option. Confirm `data.colors` reads well in practice.
- **Which dataset is canonical** for `figures/fitness-flux/data.json` (the dashboard `<select>` implies several).
- **Whether to retain a dashboard context** for fitness-flux at all, or let the living-paper page subsume it later.

## Future host (deferred — not this spec)

When the manuscript end matures (numbered figures/equations, citations, cross-references, PDF + HTML from one source), evaluate a host framework — MyST first, Quarto as the proven fallback — to author the living paper. Because components here are framework-agnostic ES modules with a pure `render`, that host consumes them via a raw-HTML or OJS block with no change to the components. The figure components are the durable investment; the host is swappable. Do not let the eventual host choice block or reshape this component work.
