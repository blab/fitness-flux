// Variant frequencies through time: per-season panels comparing empirical weekly
// frequencies (dots) with MLR-modeled frequencies (lines) for each variant.
//
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (renderFrequencyPanels)
// on 2026-06-19.
//
// data = {
//   seasonal: Array<{ timepoint, date (ISO string), variant, empirical (number|null), modeled (number|null) }>,
//   colors:   Array<{ variant, color, display_name, is_major, order }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../../lib/plot.js";
import * as d3 from "../../lib/d3.js";
import { colorScale, buildLegend } from "../../lib/colors.js";

const COL_GAP = 10; // px between panel columns
const YEAR_MS = 365.25 * 24 * 3600 * 1000;

// Seasons span ~1 year (SARS-CoV-2) or ~2 years (H3N2); tick every 6 months or
// yearly accordingly. Derived from the data so the component is dataset-agnostic.
function tickIntervalFor(dated) {
    const spans = [];
    for (const [, rows] of d3.group(dated, (d) => d.timepoint)) {
        const [lo, hi] = d3.extent(rows, (d) => d.date);
        if (lo && hi) spans.push(hi - lo);
    }
    spans.sort((a, b) => a - b);
    const median = spans.length ? spans[Math.floor(spans.length / 2)] : 0;
    return median > 1.4 * YEAR_MS ? d3.utcYear : d3.utcMonth.every(6);
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const panelFont = mode === "slide" ? "14px" : "12px";
    const legendFont = mode === "slide" ? "15px" : "14px";
    const panelW = mode === "slide" ? 300 : 232;
    const panelH = mode === "slide" ? 190 : 150;

    const scale = colorScale(data.colors);
    const dated = data.seasonal.map((d) => ({
        ...d,
        date: d.date instanceof Date ? d.date : new Date(d.date),
    }));
    const tickInterval = tickIntervalFor(dated);
    const bySeason = d3.group(dated, (d) => d.timepoint);
    const seasons = [...bySeason.keys()].sort();

    const root = document.createElement("div");
    container.appendChild(root);

    const legend = buildLegend(scale, { orientation: "horizontal", fontSize: legendFont });
    legend.style.margin = "8px 0 4px";
    root.appendChild(legend);

    const grid = document.createElement("div");
    Object.assign(grid.style, { display: "grid", gap: `6px ${COL_GAP}px`, marginTop: "8px" });
    root.appendChild(grid);

    function panelFor(season) {
        const rows = bySeason.get(season);
        const tipPoints = [];
        for (const d of rows) {
            if (typeof d.modeled === "number")
                tipPoints.push({ variant: d.variant, date: d.date, value: d.modeled, kind: "MLR" });
            if (typeof d.empirical === "number")
                tipPoints.push({ variant: d.variant, date: d.date, value: d.empirical, kind: "empirical" });
        }
        return Plot.plot({
            style: { fontSize: panelFont },
            width: panelW,
            height: panelH,
            marginLeft: 42,
            marginBottom: 22,
            marginTop: 6,
            marginRight: 28,
            color: { type: "identity" },
            x: { type: "utc", ticks: tickInterval, tickFormat: "%b %Y", nice: tickInterval, label: null },
            y: { domain: [0, 1], ticks: [0, 0.5, 1], tickFormat: ".0%", label: null },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.line(rows.filter((d) => typeof d.modeled === "number"), {
                    x: "date",
                    y: "modeled",
                    z: "variant",
                    stroke: (d) => scale.color(d.variant),
                    strokeWidth: 1.3,
                }),
                Plot.dot(rows.filter((d) => typeof d.empirical === "number"), {
                    x: "date",
                    y: "empirical",
                    fill: (d) => scale.color(d.variant),
                    r: 1,
                    fillOpacity: 0.5,
                }),
                Plot.tip(
                    tipPoints,
                    Plot.pointer({
                        x: "date",
                        y: "value",
                        title: (d) => `${scale.name(d.variant)}\n${d.kind} ${(d.value * 100).toFixed(1)}%`,
                    }),
                ),
            ],
        });
    }

    function draw(totalWidth) {
        const cols = Math.max(1, Math.floor((totalWidth + COL_GAP) / (panelW + COL_GAP)));
        grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        grid.replaceChildren(...seasons.map(panelFor));
    }

    function measure() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(panelW, w || panelW * 4);
    }

    draw(measure());

    return {
        element: root,
        resize(width) {
            draw(width ? Math.max(panelW, Math.floor(width)) : measure());
        },
        destroy() {
            root.remove();
        },
    };
}
