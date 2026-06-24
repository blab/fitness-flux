// Lineage deltas: for each parent→child Pango branch, the change in log fitness
// (y, always) against the change in a chosen "predictor" (x) — a mutation-region
// substitution count or a mutation-effect score. One panel per predictor named in
// `opts.predictors`, sharing a y-scale so slopes are comparable across regions.
// Points are colored by season on a blue(early)→red(late) ramp.
//
// A figure-wide All / Early / Late toggle subsets the points by season window
// (Early = Jan 2020–Jun 2022, Late = Jul 2022 on). Axis domains are fixed from the
// full data, so toggling fades the out-of-period points and tilts each panel's
// regression line without the axes jumping — making the "effect of spike change
// erodes over time" comparison honest. This one figure stands in for both the
// across-genome (regions as panels) and across-time (early vs late) views.
//
// data = {
//   seasons:    Array<string>,                 // chronological, for the time color domain
//   predictors: { key: { label, group } },
//   points:     Array<{ timepoint, period: "early"|"late", parent, child,
//                        delta_log_fitness, values: { <key>: number } }>
// }
// opts = { predictors?: "spike,rbd,accessory" | string[], mode?, width?, height?, initial?,
//          xmin?, xmax?, ymin?, ymax? }   // axis overrides applied to every panel
//                                          // (any subset; default shows all points)
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const GAP = 16; // px between panels
const MIN_PANEL = 300; // px below which the grid drops a column
const DEFAULT_PREDICTORS = ["spike", "rbd", "accessory"];
const STATES = [
    ["all", "All"],
    ["early", "Early"],
    ["late", "Late"],
];
const DOT_OPACITY = 0.55;

function padExtent([lo, hi], frac = 0.05) {
    if (lo === hi) return [lo - 1, hi + 1];
    const pad = (hi - lo) * frac;
    return [lo - pad, hi + pad];
}

// Average ranks (ties shared), for Spearman's rho.
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

// Least-squares fit plus Pearson r, Spearman rho, and the "doubling" (ln 2 / slope,
// i.e. the predictor change that doubles fitness) — matching the notebook's panel
// annotations.
function panelStats(points) {
    const n = points.length;
    if (n < 2) return null;
    const xs = points.map((p) => p.x), ys = points.map((p) => p.y);
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) {
        sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i];
    }
    const denom = n * sxx - sx * sx;
    if (denom === 0) return null;
    const slope = (n * sxy - sx * sy) / denom;
    return {
        slope,
        intercept: (sy - slope * sx) / n,
        r: pearson(xs, ys),
        rho: pearson(rankValues(xs), rankValues(ys)),
        doubling: slope !== 0 ? Math.log(2) / slope : null,
    };
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    let height = opts.height ?? (mode === "slide" ? 300 : 260);
    let state = STATES.some(([k]) => k === opts.initial) ? opts.initial : "all";

    // Which predictors get panels (keep only keys present in the data).
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
    const xLabelOf = (key) =>
        data.predictors[key]?.group === "Mutation region"
            ? `Change in ${labelOf(key)} mutations`
            : `Change in ${labelOf(key)}`;

    // Season → color along Mathematica's ColorData["Rainbow"] ramp (earliest
    // purple/blue → latest red), the notebook's color scheme. The 0.3 stop is its
    // exact ColorData["Rainbow"][0.3] = RGBColor[0.298, 0.566, 0.752] (#4c90c0).
    const RAINBOW = [
        [0.0, "#781c87"],
        [0.3, "#4c90c0"],
        [0.5, "#5cba4f"],
        [0.7, "#dbc740"],
        [0.85, "#ef8a2b"],
        [1.0, "#d62f26"],
    ];
    const rainbow = d3
        .scaleLinear()
        .domain(RAINBOW.map((d) => d[0]))
        .range(RAINBOW.map((d) => d[1]))
        .clamp(true);
    const seasonIndex = new Map(data.seasons.map((s, i) => [s, i]));
    const seasonSpan = Math.max(1, data.seasons.length - 1);
    const colorOf = (p) => rainbow((seasonIndex.get(p.timepoint) ?? 0) / seasonSpan);

    // Per-panel point arrays (full set; fixed positions) + fixed domains.
    const panelPoints = new Map();
    for (const key of panelKeys) {
        panelPoints.set(
            key,
            data.points
                .filter(
                    (p) =>
                        typeof p.values?.[key] === "number" &&
                        typeof p.delta_log_fitness === "number",
                )
                .map((p) => ({
                    x: p.values[key],
                    y: p.delta_log_fitness,
                    period: p.period,
                    timepoint: p.timepoint,
                    parent: p.parent,
                    child: p.child,
                    color: colorOf(p),
                })),
        );
    }
    // Optional fixed-axis overrides (any subset) applied to every panel; missing
    // bounds keep the data-derived domain (which shows all points). Clip marks only
    // when an override is active, so an override can crop without points spilling.
    const num = (v) => {
        if (v == null) return null;
        const n = typeof v === "number" ? v : Number(v);
        return Number.isFinite(n) ? n : null;
    };
    const xmin = num(opts.xmin), xmax = num(opts.xmax);
    const ymin = num(opts.ymin), ymax = num(opts.ymax);
    const clip = [xmin, xmax, ymin, ymax].some((v) => v != null);

    const xDomain = new Map(
        panelKeys.map((key) => {
            const xs = panelPoints.get(key).map((p) => p.x);
            const [lo, hi] = padExtent(d3.extent([0, ...xs]));
            return [key, [xmin ?? lo, xmax ?? hi]];
        }),
    );
    const allY = panelKeys.flatMap((key) => panelPoints.get(key).map((p) => p.y));
    const [ylo, yhi] = padExtent(d3.extent([0, ...allY]));
    const yDomain = [ymin ?? ylo, ymax ?? yhi];

    const inState = (p) => state === "all" || p.period === state;

    const root = document.createElement("div");
    container.appendChild(root);
    let lastWidth = 0;

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PANEL, w || 820);
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
        STATES.forEach(([key, label], i) => {
            const active = key === state;
            const button = document.createElement("button");
            button.type = "button";
            button.textContent = label;
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
                if (state === key) return;
                state = key;
                animateToggle();
            });
            wrap.appendChild(button);
        });
        return wrap;
    }

    function buildTimeLegend() {
        const wrap = document.createElement("div");
        Object.assign(wrap.style, {
            display: "flex",
            alignItems: "center",
            gap: "8px",
            fontSize: "12px",
            color: "#555",
        });
        const first = document.createElement("span");
        first.textContent = data.seasons[0];
        const bar = document.createElement("span");
        Object.assign(bar.style, {
            display: "inline-block",
            width: "110px",
            height: "10px",
            borderRadius: "2px",
            background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1]
                .map((t) => rainbow(t))
                .join(", ")})`,
        });
        const last = document.createElement("span");
        last.textContent = data.seasons[data.seasons.length - 1];
        wrap.append(first, bar, last);
        return wrap;
    }

    function panel(key, index, panelW, panelH) {
        const wrap = document.createElement("div");
        wrap.className = "ld-panel";
        const title = document.createElement("div");
        title.textContent = labelOf(key);
        Object.assign(title.style, {
            fontSize: "13px",
            fontWeight: "600",
            color: "#333",
            marginBottom: "2px",
        });

        const pts = panelPoints.get(key);
        const xDom = xDomain.get(key);
        // Fit and stats use only the points actually shown — in the active period
        // and within the (possibly overridden) x/y domains, so they describe the
        // visible cloud rather than cropped-away outliers.
        const inView = (p) =>
            p.x >= xDom[0] && p.x <= xDom[1] && p.y >= yDomain[0] && p.y <= yDomain[1];
        const shown = pts.filter((p) => inState(p) && inView(p));
        const st = panelStats(shown);
        const fit = st
            ? [
                  { x: xDom[0], y: st.intercept + st.slope * xDom[0] },
                  { x: xDom[1], y: st.intercept + st.slope * xDom[1] },
              ]
            : [];
        const fmt = (v, d) => (v == null || !Number.isFinite(v) ? "n/a" : v.toFixed(d));
        const mutsUnit = data.predictors[key]?.group === "Mutation region" ? " muts" : "";
        const statRows = st
            ? [
                  `slope = ${fmt(st.slope, 2)}`,
                  `doubling = ${fmt(st.doubling, 1)}${mutsUnit}`,
                  `r = ${fmt(st.r, 2)}`,
                  `ρ = ${fmt(st.rho, 2)}`,
              ]
            : [];
        const statFont = mode === "slide" ? 15 : 13;
        const statLineH = statFont * 1.5;

        const plot = Plot.plot({
            style: { fontSize: axisFont },
            width: panelW,
            height: panelH,
            marginLeft: 52,
            marginRight: 12,
            marginTop: 8,
            marginBottom: 40,
            clip,
            color: { type: "identity" },
            x: {
                domain: xDom,
                label: xLabelOf(key),
                labelAnchor: "center",
                labelArrow: "none",
            },
            y: {
                domain: yDomain,
                label: index === 0 ? "Change in log fitness" : null,
                labelAnchor: "center",
                labelArrow: "none",
            },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.dot(pts, {
                    x: "x",
                    y: "y",
                    r: 3,
                    fill: "color",
                    fillOpacity: (d) => (inState(d) ? DOT_OPACITY : 0),
                    stroke: "none",
                }),
                ...(fit.length
                    ? [Plot.line(fit, { x: "x", y: "y", stroke: "#222", strokeWidth: 1.5 })]
                    : []),
                ...statRows.map((text, j) =>
                    Plot.text([text], {
                        frameAnchor: "top-left",
                        dx: 6,
                        dy: 6 + j * statLineH,
                        textAnchor: "start",
                        lineAnchor: "top",
                        fontSize: statFont,
                        fill: "#222",
                    }),
                ),
                Plot.tip(
                    shown,
                    Plot.pointer({
                        x: "x",
                        y: "y",
                        title: (d) =>
                            `${d.parent} → ${d.child}\nseason ${d.timepoint}\nΔ ${labelOf(key)} ${d.x}\nΔ log fitness ${d.y.toFixed(3)}`,
                    }),
                ),
            ],
        });
        // Italicize only the leading r / ρ symbol (not its "= value"), per the
        // summary-statistic convention, by splitting that text into two tspans.
        for (const t of plot.querySelectorAll('g[aria-label="text"] text')) {
            const s = t.textContent;
            if (!/^(r|ρ) = /.test(s)) continue;
            const cut = s.indexOf(" = ");
            const sym = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
            sym.setAttribute("font-style", "italic");
            sym.textContent = s.slice(0, cut);
            const rest = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
            rest.textContent = s.slice(cut);
            t.replaceChildren(sym, rest);
        }

        wrap.append(title, plot);
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
        header.append(buildTimeLegend(), buildToggle());
        root.appendChild(header);

        const grid = document.createElement("div");
        const cols = Math.max(
            1,
            Math.min(
                panelKeys.length,
                Math.floor((totalWidth + GAP) / (MIN_PANEL + GAP)),
            ),
        );
        const panelW = Math.max(
            MIN_PANEL,
            Math.floor((totalWidth - (cols - 1) * GAP) / cols),
        );
        const panelH = Math.round(Math.min(height, Math.max(220, panelW * 0.78)));
        Object.assign(grid.style, {
            display: "grid",
            gridTemplateColumns: `repeat(${cols}, minmax(0, 1fr))`,
            gap: `${GAP}px`,
            marginTop: "10px",
        });
        panelKeys.forEach((key, i) => grid.appendChild(panel(key, i, panelW, panelH)));
        root.appendChild(grid);
    }

    // Toggle: rebuild in the new state, then transition each panel's out-of-period
    // dots fading (fill-opacity) and its regression line tilting (path `d`) from the
    // old values to the new ones. Domains are fixed, so points never move.
    function animateToggle(duration = 450) {
        const before = [...root.querySelectorAll(".ld-panel")].map((p) => ({
            dots: [...p.querySelectorAll('g[aria-label="dot"] circle')].map((c) =>
                c.getAttribute("fill-opacity"),
            ),
            line: p.querySelector('g[aria-label="line"] path')?.getAttribute("d") ?? null,
        }));
        draw(lastWidth || measureWidth());
        const panelsNow = [...root.querySelectorAll(".ld-panel")];
        panelsNow.forEach((p, i) => {
            const old = before[i];
            if (!old) return;
            const circles = [...p.querySelectorAll('g[aria-label="dot"] circle')];
            const finalOpacity = circles.map((c) => c.getAttribute("fill-opacity"));
            circles.forEach((c, j) => {
                if (old.dots[j] != null) c.setAttribute("fill-opacity", old.dots[j]);
            });
            if (circles.length)
                d3.selectAll(circles)
                    .data(finalOpacity)
                    .transition()
                    .duration(duration)
                    .attr("fill-opacity", (d) => d);
            const linePath = p.querySelector('g[aria-label="line"] path');
            if (linePath && old.line != null) {
                const finalD = linePath.getAttribute("d");
                linePath.setAttribute("d", old.line);
                d3.select(linePath).transition().duration(duration).attr("d", finalD);
            }
        });
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
