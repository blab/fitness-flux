// Lineage-delta trends through time: how the relationship between a per-branch
// predictor (e.g. regional substitution count) and change in log fitness evolves
// season by season. For each predictor and season, the parent→child branches are
// summarized into one statistic — regression slope (default), Pearson r, or
// Spearman ρ (toggled) — and plotted as a line per predictor across seasons.
//
// Reads the same per-branch JSON as the lineage-deltas scatter, so predictors= and
// the xmin/xmax/ymin/ymax crop select the same visible cloud; the statistic is
// computed over those shown points (n >= 3 per season).
//
// data = {
//   seasons:    Array<string>,                 // chronological x-axis
//   predictors: { key: { label, group } },
//   points:     Array<{ timepoint, parent, child, delta_log_fitness, values: { <key>: number } }>
// }
// opts = { predictors?: "s1,rbd,…" | string[], stat?: "slope"|"pearson"|"spearman",
//          mode?, width?, height?, xmin?, xmax?, ymin?, ymax? }
//
// Pure: no fetching, no ResizeObserver. Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";
import { colorScale, buildLegend, linkLegendHighlight } from "../lib/colors.js";

const GAP = 16;
const MIN_PLOT = 360;
const MIN_N = 3; // seasons with fewer shared branches are skipped (cf. linear_models.py)
const DEFAULT_PREDICTORS = ["spike", "rbd", "accessory"];
// Categorical palette for the predictor lines (seaborn "deep"); cycles if needed.
const PALETTE = [
    "#4c72b0", "#dd8452", "#55a868", "#c44e52",
    "#8172b3", "#937860", "#da8bc3", "#8c8c8c",
];
// Fixed colors for the common genome regions (palette used for anything else).
const REGION_COLORS = {
    s1: "#dd8452", // orange
    rbd: "#c44e52", // red
    orf1ab: "#4c72b0", // blue
    accessory: "#55a868", // green
};
const STATS = [
    ["slope", "Slope", "Regression slope", 3],
    ["pearson", "Pearson <i>r</i>", "Pearson r", 2],
    ["spearman", "Spearman <i>ρ</i>", "Spearman ρ", 2],
];

function padExtent([lo, hi], frac = 0.05) {
    if (lo === hi) return [lo - 1, hi + 1];
    const pad = (hi - lo) * frac;
    return [lo - pad, hi + pad];
}

function rankValues(values) {
    const order = values.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
    const ranks = new Array(values.length);
    for (let i = 0; i < order.length; ) {
        let j = i;
        while (j + 1 < order.length && order[j + 1][0] === order[i][0]) j++;
        const avg = (i + j) / 2 + 1;
        for (let k = i; k <= j; k++) ranks[order[k][1]] = avg;
        i = j + 1;
    }
    return ranks;
}

function pearson(xs, ys) {
    const n = xs.length;
    let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    for (let i = 0; i < n; i++) {
        sx += xs[i]; sy += ys[i];
        sxx += xs[i] * xs[i]; syy += ys[i] * ys[i]; sxy += xs[i] * ys[i];
    }
    const dx = n * sxx - sx * sx, dy = n * syy - sy * sy;
    return dx > 0 && dy > 0 ? (n * sxy - sx * sy) / Math.sqrt(dx * dy) : null;
}

function seasonStat(points, stat) {
    if (points.length < MIN_N) return null;
    const xs = points.map((p) => p.x), ys = points.map((p) => p.y);
    if (stat === "pearson") return pearson(xs, ys);
    if (stat === "spearman") return pearson(rankValues(xs), rankValues(ys));
    const n = xs.length;
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) {
        sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i];
    }
    const denom = n * sxx - sx * sx;
    return denom === 0 ? null : (n * sxy - sx * sy) / denom;
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    let height = opts.height ?? (mode === "slide" ? 340 : 300);
    let stat = STATS.some(([k]) => k === opts.stat) ? opts.stat : "slope";

    const requested =
        typeof opts.predictors === "string"
            ? opts.predictors.split(",").map((s) => s.trim()).filter(Boolean)
            : Array.isArray(opts.predictors)
              ? opts.predictors
              : DEFAULT_PREDICTORS;
    let panelKeys = requested.filter((k) => data.predictors[k]);
    if (panelKeys.length === 0) {
        panelKeys = Object.keys(data.predictors)
            .filter((k) => DEFAULT_PREDICTORS.includes(k))
            .slice(0, 3);
    }
    const labelOf = (key) => data.predictors[key]?.label ?? key;
    // Synthetic color scale over the predictors, so we can reuse the shared
    // square-swatch legend and Auspice-style legend-hover highlighting.
    const colors = panelKeys.map((key, i) => ({
        variant: key,
        color: REGION_COLORS[key] ?? PALETTE[i % PALETTE.length],
        display_name: labelOf(key),
        is_major: true,
        order: i,
    }));
    const scale = colorScale(colors);

    const num = (v) => {
        if (v == null) return null;
        const n = typeof v === "number" ? v : Number(v);
        return Number.isFinite(n) ? n : null;
    };
    const xmin = num(opts.xmin), xmax = num(opts.xmax);
    const ymin = num(opts.ymin), ymax = num(opts.ymax);
    const inView = (x, y) =>
        (xmin == null || x >= xmin) && (xmax == null || x <= xmax) &&
        (ymin == null || y >= ymin) && (ymax == null || y <= ymax);

    const seasons = data.seasons;
    // predictor -> season -> [{x, y}] over the cropped points
    const grouped = new Map();
    for (const key of panelKeys) {
        const bySeason = new Map(seasons.map((s) => [s, []]));
        for (const p of data.points) {
            const x = p.values?.[key], y = p.delta_log_fitness;
            if (typeof x !== "number" || typeof y !== "number") continue;
            if (!inView(x, y)) continue;
            bySeason.get(p.timepoint)?.push({ x, y });
        }
        grouped.set(key, bySeason);
    }

    // Series for the current statistic: one point per (predictor, season) with a
    // defined statistic, in chronological order so each predictor's line connects.
    function buildSeries(which) {
        const rows = [];
        for (const key of panelKeys) {
            const bySeason = grouped.get(key);
            for (const season of seasons) {
                const value = seasonStat(bySeason.get(season), which);
                if (value != null && Number.isFinite(value)) {
                    rows.push({ predictor: key, season, value, color: scale.color(key) });
                }
            }
        }
        return rows;
    }

    const root = document.createElement("div");
    container.appendChild(root);
    let lastWidth = 0;

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PLOT, w || 820);
    }

    function buildToggle() {
        const wrap = document.createElement("div");
        Object.assign(wrap.style, {
            display: "inline-flex",
            border: "1px solid #ccc",
            borderRadius: "5px",
            overflow: "hidden",
            fontSize: "13px",
        });
        STATS.forEach(([key, label], i) => {
            const active = key === stat;
            const button = document.createElement("button");
            button.type = "button";
            button.innerHTML = label;
            Object.assign(button.style, {
                border: "none",
                borderLeft: i ? "1px solid #ccc" : "none",
                padding: "4px 12px",
                cursor: active ? "default" : "pointer",
                font: "inherit",
                background: active ? "#333" : "#fff",
                color: active ? "#fff" : "#555",
            });
            button.addEventListener("click", () => {
                if (stat === key) return;
                stat = key;
                animateToggle();
            });
            wrap.appendChild(button);
        });
        return wrap;
    }

    function draw(totalWidth) {
        lastWidth = totalWidth;
        root.replaceChildren();

        const header = document.createElement("div");
        Object.assign(header.style, {
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
            gap: "16px",
            flexWrap: "wrap",
            margin: "0 0 4px",
        });
        header.append(
            buildLegend(scale, { orientation: "horizontal", fontSize: "13px" }),
            buildToggle(),
        );
        root.appendChild(header);

        const series = buildSeries(stat);
        const meta = STATS.find(([k]) => k === stat);
        const yLabel = meta[2];
        const yDom = padExtent(d3.extent([0, ...series.map((d) => d.value)]));

        const holder = document.createElement("div");
        Object.assign(holder.style, { overflowX: "auto", marginTop: "10px" });
        holder.appendChild(
            Plot.plot({
                style: { fontSize: axisFont },
                width: totalWidth,
                height,
                marginLeft: 58,
                marginRight: 14,
                marginTop: 10,
                marginBottom: 52,
                color: { type: "identity" },
                x: { type: "point", domain: seasons, label: null, tickRotate: -40 },
                y: { domain: yDom, label: yLabel, labelAnchor: "center", labelArrow: "none" },
                marks: [
                    Plot.frame({ anchor: "left", stroke: "#333" }),
                    Plot.frame({ anchor: "bottom", stroke: "#333" }),
                    Plot.ruleY([0], { stroke: "#bbb", strokeDasharray: "3,3" }),
                    Plot.line(series, {
                        x: "season",
                        y: "value",
                        z: "predictor",
                        stroke: "color",
                        strokeWidth: 2,
                    }),
                    Plot.dot(series, {
                        x: "season",
                        y: "value",
                        fill: "color",
                        r: 3.5,
                        stroke: "white",
                        strokeWidth: 0.6,
                    }),
                    Plot.tip(
                        series,
                        Plot.pointer({
                            x: "season",
                            y: "value",
                            title: (d) =>
                                `${labelOf(d.predictor)}\nseason ${d.season}\n${meta[1].replace(/<\/?i>/g, "")} = ${d.value.toFixed(meta[3])}`,
                        }),
                    ),
                ],
            }),
        );
        root.appendChild(holder);
    }

    // Toggle the statistic with each predictor's line and markers gliding to their
    // new y-positions: rebuild, then transition path `d` and circle cx/cy from old.
    function animateToggle(duration = 450) {
        const oldLines = [...root.querySelectorAll('g[aria-label="line"] path')].map((p) =>
            p.getAttribute("d"),
        );
        const oldDots = [...root.querySelectorAll('g[aria-label="dot"] circle')].map((c) => [
            c.getAttribute("cx"),
            c.getAttribute("cy"),
        ]);
        draw(lastWidth || measureWidth());
        const lines = [...root.querySelectorAll('g[aria-label="line"] path')];
        const finalLines = lines.map((p) => p.getAttribute("d"));
        lines.forEach((p, j) => oldLines[j] != null && p.setAttribute("d", oldLines[j]));
        if (lines.length)
            d3.selectAll(lines).data(finalLines).transition().duration(duration).attr("d", (d) => d);
        const dots = [...root.querySelectorAll('g[aria-label="dot"] circle')];
        const finalDots = dots.map((c) => [c.getAttribute("cx"), c.getAttribute("cy")]);
        dots.forEach((c, j) => {
            if (oldDots[j]) {
                c.setAttribute("cx", oldDots[j][0]);
                c.setAttribute("cy", oldDots[j][1]);
            }
        });
        if (dots.length) {
            const sel = d3.selectAll(dots).data(finalDots).transition().duration(duration);
            sel.attr("cx", (d) => d[0]).attr("cy", (d) => d[1]);
        }
    }

    draw(measureWidth());
    const unlink = linkLegendHighlight(root, scale);

    return {
        element: root,
        resize(width, newHeight) {
            if (newHeight) height = newHeight;
            draw(width ? Math.max(MIN_PLOT, Math.floor(width)) : measureWidth());
        },
        destroy() {
            unlink();
            root.remove();
        },
    };
}
