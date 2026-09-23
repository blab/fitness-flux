// Clade-distance figure: why clade labels are a lossy stand-in for HA1 similarity.
//
// Two panels:
//
//   * a heatmap of pairwise HA1 amino-acid distance between clade MRCA sequences,
//     with clades ordered root-to-tip so the tree's block structure is visible;
//   * the same pairwise distances split by what the clade LABELS imply about
//     relatedness (ancestor/descendant, sibling, unrelated). The three
//     distributions overlap heavily, which is the point: a label-based forecast
//     score charges one flat penalty across a 0-34 amino-acid range.
//
// data = {
//   label, gene,                          // "H3N2", "HA1"
//   names: string[],                      // clades, root-to-tip order
//   clades: Array<{ clade, n_subs, order }>,
//   matrix: number[][],                   // matrix[i][j] = aa distance, names order
//   pairs: Array<{ a, b, d, rel }>,       // upper triangle, rel-tagged
//   relationships: Array<{ key, label, color }>,
//   stats: { [rel]: { n, mean, min, max } },
//   sequential: string[],                 // single-hue ramp, light -> dark
//   summary: { n_clades, n_variable_positions, max_distance, mean_distance,
//              mean_parent_child_distance, identical_pairs, aliases, ... }
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";

const FRAME = "#333";
const GRID = "#e1e0d9";
const MUTED = "#898781";
const INK = "#0b0b0b";
const GAP = 22;
const MIN_PANEL = 300;

function fontStack(size) {
    return `${size}px system-ui, -apple-system, "Segoe UI", sans-serif`;
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "13px" : "11px";
    const headFont = mode === "slide" ? 14 : 12;

    const names = data.names ?? [];
    const matrix = data.matrix ?? [];
    const gene = data.gene ?? "HA1";
    const ramp = data.sequential ?? ["#eef4fc", "#2a78d6", "#1a4f93"];
    const relationships = data.relationships ?? [];
    const stats = data.stats ?? {};
    const summary = data.summary ?? {};

    // Heatmap cells, long form. Diagonal included so the block structure is anchored.
    const cells = [];
    for (let i = 0; i < names.length; i++) {
        for (let j = 0; j < names.length; j++) {
            const d = matrix[i]?.[j];
            if (d != null) cells.push({ a: names[i], b: names[j], d });
        }
    }

    // Strip plot: one row per label relationship, jittered so density is readable.
    // Jitter is deterministic (hashed from the pair) so the figure is stable
    // across renders rather than reshuffling on every resize.
    const relIndex = Object.fromEntries(relationships.map((r, i) => [r.key, i]));
    const strip = (data.pairs ?? []).map((p) => {
        let h = 0;
        const key = `${p.a}|${p.b}`;
        for (let i = 0; i < key.length; i++) h = (h * 31 + key.charCodeAt(i)) | 0;
        const jitter = ((h % 1000) / 1000 - 0.5) * 0.62;
        return { ...p, row: relIndex[p.rel] ?? 0, y: (relIndex[p.rel] ?? 0) + jitter };
    });
    const means = relationships
        .filter((r) => stats[r.key])
        .map((r) => ({
            rel: r.key, row: relIndex[r.key], mean: stats[r.key].mean,
            n: stats[r.key].n, color: r.color,
        }));

    const colorOf = Object.fromEntries(relationships.map((r) => [r.key, r.color]));
    const labelOf = Object.fromEntries(relationships.map((r) => [r.key, r.label]));
    const maxD = summary.max_distance ?? Math.max(...cells.map((c) => c.d), 1);

    // Label every clade when there is room, otherwise thin them out; the tooltip
    // carries the rest, and the root-to-tip ordering keeps the axis legible.
    // Tick labels need ~12px of axis each; thin them out rather than let them
    // collide, and lean on the tooltip plus the root-to-tip ordering for the rest.
    function tickSubset(plotArea) {
        const perLabel = plotArea / Math.max(names.length, 1);
        const step = Math.max(1, Math.ceil(12 / Math.max(perLabel, 1)));
        return names.filter((_, i) => i % step === 0);
    }

    const root = document.createElement("div");
    container.appendChild(root);

    const grid = document.createElement("div");
    Object.assign(grid.style, { display: "grid", gap: `${GAP}px` });
    root.appendChild(grid);

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PANEL, w || 860);
    }

    function heading(title, note) {
        const head = document.createElement("div");
        Object.assign(head.style, { font: fontStack(headFont), margin: "0 0 3px 56px" });
        const name = document.createElement("span");
        name.textContent = title;
        Object.assign(name.style, { fontWeight: "600", color: INK });
        const extra = document.createElement("span");
        extra.textContent = note ? `  ·  ${note}` : "";
        Object.assign(extra.style, { color: MUTED });
        head.append(name, extra);
        return head;
    }

    function heatmapPlot(panelW) {
        const size = Math.max(220, Math.min(panelW, 560));
        const ticks = tickSubset(size - 64);
        return Plot.plot({
            style: { fontSize: axisFont, background: "transparent" },
            width: size,
            height: size,
            marginLeft: 56,
            marginRight: 14,
            marginTop: 6,
            marginBottom: 58,
            x: { domain: names, ticks, tickRotate: -90, label: null },
            y: { domain: names, ticks, label: null },
            color: {
                type: "linear",
                domain: [0, maxD],
                range: [ramp[0], ramp[ramp.length - 1]],
                interpolate: "rgb",
                legend: true,
                label: `${gene} amino-acid differences`,
                width: Math.min(size, 320),
            },
            marks: [
                Plot.cell(cells, { x: "b", y: "a", fill: "d", inset: 0 }),
                Plot.frame({ stroke: FRAME }),
                Plot.tip(
                    cells,
                    Plot.pointer({
                        x: "b",
                        y: "a",
                        title: (d) => `${d.a} vs ${d.b}\n${d.d} ${gene} amino acid${d.d === 1 ? "" : "s"}`,
                    }),
                ),
            ],
        });
    }

    function stripPlot(panelW, panelH) {
        return Plot.plot({
            style: { fontSize: axisFont, background: "transparent" },
            width: panelW,
            height: panelH,
            marginLeft: 148,
            marginRight: 16,
            marginTop: 20,
            marginBottom: 40,
            x: {
                domain: [0, maxD * 1.02],
                label: `${gene} amino-acid differences between clades`,
                labelAnchor: "center",
                labelArrow: "none",
            },
            y: {
                type: "linear",
                // Reversed so the rows read top-to-bottom in the order the legend
                // lists them: closest label relationship first.
                reverse: true,
                domain: [-0.88, relationships.length - 0.12],
                ticks: relationships.map((_, i) => i),
                tickFormat: (i) => relationships[Math.round(i)]?.label ?? "",
                label: null,
                grid: false,
            },
            marks: [
                Plot.gridX({ stroke: GRID, strokeOpacity: 1 }),
                Plot.dot(strip, {
                    x: "d",
                    y: "y",
                    r: 1.9,
                    fill: (d) => colorOf[d.rel],
                    fillOpacity: 0.45,
                    stroke: "none",
                }),
                // Mean marker: a short vertical rule in the series colour, ringed in
                // surface so it stays visible where the dots are dense. ruleX keeps
                // the y scale continuous, which the jittered dots depend on.
                Plot.ruleX(means, {
                    x: "mean", y1: (d) => d.row - 0.42, y2: (d) => d.row + 0.42,
                    stroke: "#fcfcfb", strokeWidth: 5,
                }),
                Plot.ruleX(means, {
                    x: "mean", y1: (d) => d.row - 0.42, y2: (d) => d.row + 0.42,
                    stroke: (d) => d.color, strokeWidth: 2.5,
                }),
                Plot.text(means, {
                    x: maxD * 1.0,
                    y: "row",
                    text: (d) => `mean ${d.mean} aa  ·  n = ${d.n}`,
                    textAnchor: "end",
                    dy: -13,
                    fontSize: 10,
                    fontWeight: 600,
                    fill: INK,
                    stroke: "#fcfcfb",
                    strokeWidth: 4,
                }),
                Plot.frame({ anchor: "left", stroke: FRAME }),
                Plot.frame({ anchor: "bottom", stroke: FRAME }),
                Plot.tip(
                    strip,
                    Plot.pointer({
                        x: "d",
                        y: "y",
                        title: (d) => `${d.a} vs ${d.b}\n${d.d} ${gene} amino acid${d.d === 1 ? "" : "s"}\n${labelOf[d.rel]}`,
                    }),
                ),
            ],
        });
    }

    function draw(totalWidth) {
        const width = Math.max(MIN_PANEL, Math.floor(totalWidth));
        const cols = width >= 2 * MIN_PANEL + GAP ? 2 : 1;
        const panelW = Math.max(MIN_PANEL, Math.floor((width - (cols - 1) * GAP) / cols));

        const left = document.createElement("div");
        left.append(
            heading(
                `Pairwise ${gene} distance`,
                `${summary.n_clades ?? names.length} clades, ordered root to tip`,
            ),
            heatmapPlot(panelW),
        );

        const right = document.createElement("div");
        const identical = (summary.identical_pairs ?? []).length;
        right.append(
            heading(
                "What the clade labels imply",
                `${(data.pairs ?? []).length} clade pairs`,
            ),
            stripPlot(panelW, Math.round(Math.min(360, Math.max(230, panelW * 0.62)))),
        );

        const note = document.createElement("div");
        Object.assign(note.style, {
            font: fontStack(mode === "slide" ? 12 : 10.5),
            color: MUTED,
            margin: "6px 0 0 132px",
            maxWidth: `${panelW - 140}px`,
            lineHeight: "1.45",
        });
        note.textContent =
            `Related and unrelated labels overlap across the full range` +
            (identical ? `, and ${identical} clade pair${identical === 1 ? " has" : "s have"} identical ${gene} sequences` : "") +
            `. A label-based score charges the same penalty for every one of these pairs.`;
        right.append(note);

        grid.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
        grid.replaceChildren(left, right);
    }

    draw(measureWidth());

    return {
        element: root,
        resize(width) {
            draw(width ? Math.max(MIN_PANEL, Math.floor(width)) : measureWidth());
        },
        destroy() {
            root.remove();
        },
    };
}
