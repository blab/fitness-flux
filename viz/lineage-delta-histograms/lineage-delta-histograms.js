// Lineage-delta histograms: across all parent→child Pango lineage branch
// comparisons, the distribution of the change in each chosen mutation predictor
// (one panel per predictor) and, in a final panel, the change in log fitness.
// Proportion histograms sharing a y-axis, so the panels are directly comparable.
//
// data = the lineage-deltas viz JSON:
//   { predictors: { key: { label, group } },
//     points: [{ delta_log_fitness, values: { <key>: number } }, ...] }
// opts = { predictors?: "spike,nonspike" | string[], mode?, width?, height? }
//
// Pure: no fetching, no ResizeObserver. Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const GAP = 16;
const MIN_PANEL = 240; // px below which the grid drops a column
const DEFAULT_PREDICTORS = ["spike", "nonspike"];
const BAR_FILL = "#9ecae1"; // light blue
const FITNESS_BIN = 0.05;

// Histogram of `values` into bins of width `interval` centered on multiples of
// `interval` (so unit bins are centered on integers, matching the notebook's
// {min-0.5, max, 1} / {min-0.025, max, 0.05} bin specs). Returns bars with the
// proportion of all values in each bin.
function histogram(values, interval) {
    const n = values.length;
    const counts = new Map();
    for (const v of values) {
        const idx = Math.round(v / interval);
        counts.set(idx, (counts.get(idx) || 0) + 1);
    }
    const indices = [...counts.keys()];
    const lo = Math.min(...indices), hi = Math.max(...indices);
    const bars = [];
    for (let idx = lo; idx <= hi; idx++) {
        const center = idx * interval;
        bars.push({
            x0: center - interval / 2,
            x1: center + interval / 2,
            center,
            p: (counts.get(idx) || 0) / n,
        });
    }
    return bars;
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    let height = opts.height ?? (mode === "slide" ? 260 : 230);

    const requested =
        typeof opts.predictors === "string"
            ? opts.predictors.split(",").map((s) => s.trim()).filter(Boolean)
            : Array.isArray(opts.predictors)
              ? opts.predictors
              : DEFAULT_PREDICTORS;
    let mutationKeys = requested.filter((k) => data.predictors[k]);
    if (mutationKeys.length === 0) {
        mutationKeys = DEFAULT_PREDICTORS.filter((k) => data.predictors[k]);
    }
    const labelOf = (key) => data.predictors[key]?.label ?? key;

    // One panel per mutation predictor, then the Δ log-fitness panel. Bars are
    // proportions over all branches; the x-axis is focused on the central ~99% of
    // the mass so rare founder/recombinant outliers don't stretch it flat (the
    // off-range bars clip away).
    const finite = (vals) => vals.filter((v) => typeof v === "number");
    const makePanel = (label, values, interval, decimals) => {
        const sorted = values.slice().sort((a, b) => a - b);
        const lo = d3.quantile(sorted, 0.01), hi = d3.quantile(sorted, 0.99);
        return {
            label,
            decimals,
            bars: histogram(values, interval),
            xDom: [
                Math.floor(lo / interval) * interval - interval / 2,
                Math.ceil(hi / interval) * interval + interval / 2,
            ],
        };
    };
    const panels = mutationKeys.map((key) =>
        makePanel(
            `Change in ${labelOf(key)} mutations`,
            finite(data.points.map((p) => p.values?.[key])),
            1,
            0,
        ),
    );
    panels.push(
        makePanel(
            "Change in log fitness",
            finite(data.points.map((p) => p.delta_log_fitness)),
            FITNESS_BIN,
            2,
        ),
    );

    // Shared y-domain (max bar proportion across every panel).
    const maxP = Math.max(...panels.flatMap((panel) => panel.bars.map((b) => b.p)));
    const yDom = [0, maxP * 1.04];

    const root = document.createElement("div");
    container.appendChild(root);
    let lastWidth = 0;

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PANEL, w || 820);
    }

    function panel(spec, index, panelW, panelH) {
        const fmt = (v) => v.toFixed(spec.decimals);
        return Plot.plot({
            style: { fontSize: axisFont },
            width: panelW,
            height: panelH,
            marginLeft: 48,
            marginRight: 10,
            marginTop: 8,
            marginBottom: 40,
            clip: true,
            x: { domain: spec.xDom, label: spec.label, labelAnchor: "center", labelArrow: "none" },
            y: {
                domain: yDom,
                label: index === 0 ? "Frequency" : null,
                labelAnchor: "center",
                labelArrow: "none",
                tickFormat: "%",
            },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.rectY(spec.bars, {
                    x1: "x0",
                    x2: "x1",
                    y: "p",
                    fill: BAR_FILL,
                    stroke: "#fff",
                    strokeWidth: 0.5,
                }),
                Plot.ruleX([0], { stroke: "#888", strokeDasharray: "3,3" }),
                Plot.tip(
                    spec.bars,
                    Plot.pointerX({
                        x: "center",
                        y: "p",
                        title: (d) => `${spec.label}: ${fmt(d.center)}\n${(d.p * 100).toFixed(1)}%`,
                    }),
                ),
            ],
        });
    }

    function draw(totalWidth) {
        lastWidth = totalWidth;
        const cols = Math.max(
            1,
            Math.min(panels.length, Math.floor((totalWidth + GAP) / (MIN_PANEL + GAP))),
        );
        const panelW = Math.max(
            MIN_PANEL,
            Math.floor((totalWidth - (cols - 1) * GAP) / cols),
        );
        const panelH = Math.round(Math.min(height, Math.max(170, panelW * 0.85)));
        const grid = document.createElement("div");
        Object.assign(grid.style, {
            display: "grid",
            gridTemplateColumns: `repeat(${cols}, minmax(0, 1fr))`,
            gap: `${GAP}px`,
        });
        panels.forEach((spec, i) => grid.appendChild(panel(spec, i, panelW, panelH)));
        root.replaceChildren(grid);
    }

    draw(measureWidth());

    return {
        element: root,
        resize(width, newHeight) {
            if (newHeight) height = newHeight;
            draw(width ? Math.max(MIN_PANEL, Math.floor(width)) : measureWidth());
        },
        destroy() {
            root.remove();
        },
    };
}
