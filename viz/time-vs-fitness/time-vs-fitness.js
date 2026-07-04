// Fitness through time: each variant is a horizontal violin whose half-height is
// proportional to its empirical frequency. Two y-axis modes, toggled in the UI:
//   • "flux"     — y is the variant's fixed scaffolded log fitness (cumulative
//                  fitness flux); each blob sits at a constant height.
//   • "relative" — y is log fitness relative to the daily population mean:
//                  log_fitness_i − log(f̄(t)), with f̄(t) = Σ_i x_i e^{log_fitness_i}.
//                  Blobs start above 0, drift down as the mean rises, and fade past 0.
//
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (renderFitnessFlux,
// drawFitnessFlux) on 2026-06-19.
//
// data = {
//   fitnessFlux: Array<{ variant, date (ISO string), log_fitness (number), emp_freq (number) }>,
//   colors:      Array<{ variant, color, display_name, is_major, order }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", measure?: "flux"|"relative", width?, height? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";
import { colorScale, buildLegend, linkLegendHighlight } from "../lib/colors.js";

// Blob half-height as a fraction of the y-axis span, so the visual thickness is
// uniform across datasets and across both modes. Calibrated so the reference
// (SARS-CoV-2 clades cumulative, span ≈ 2.66) keeps its ~0.13 thickness.
const HALF_HEIGHT_FRACTION = 0.049;
const GAP = 16; // px between plot and legend
const MIN_PLOT = 360; // px below which the legend folds under

function decimalYear(date) {
    const d = new Date(date);
    const y = d.getUTCFullYear();
    const start = Date.UTC(y, 0, 1);
    return y + (d.getTime() - start) / (Date.UTC(y + 1, 0, 1) - start);
}

function linearFitSlope(points) {
    const n = points.length;
    if (n < 2) return null;
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (const [x, y] of points) { sx += x; sy += y; sxx += x * x; sxy += x * y; }
    const denom = n * sxx - sx * sx;
    return denom === 0 ? null : (n * sxy - sx * sy) / denom;
}

// Wave velocity from the dominant variant (freq > 0.5) at each date — flux mode only.
function slopeText(points) {
    const dominant = [];
    for (const [, rows] of d3.group(points, (d) => +d.date)) {
        const top = d3.greatest(rows, (d) => d.empFreq);
        if (top && top.empFreq > 0.5) dominant.push([decimalYear(top.date), top.logFitness]);
    }
    const slope = linearFitSlope(dominant);
    return slope && slope > 0
        ? `slope = ${slope.toFixed(2)} per year   ·   doubling = ${(Math.log(2) / slope).toFixed(1)} years`
        : "";
}

function linearFit(points) {
    const n = points.length;
    if (n < 2) return null;
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (const [x, y] of points) { sx += x; sy += y; sxx += x * x; sxy += x * y; }
    const denom = n * sxx - sx * sx;
    if (denom === 0) return null;
    const slope = (n * sxy - sx * sy) / denom;
    return { slope, intercept: (sy - slope * sx) / n };
}

// "2020-2021.5,2021.5-2023" -> [{lo, hi, below}] as half-open decimal-year ranges
// [lo, hi); bounds are literal decimal years, so 2022 = Jan 1 2022 and 2021.5 =
// ~Jul 1 2021. An optional ":side" suffix flips a segment's label to the other side
// of its line -- "below"/"down"/"right" put it under-and-downstream (default above);
// e.g. "2021-2022.5,2022.5-2026:below".
function parseRegressions(spec) {
    if (!spec) return [];
    return String(spec).split(",").map((s) => s.trim()).filter(Boolean).map((token) => {
        const [range, place] = token.split(":");
        const [a, b] = range.split("-").map((y) => parseFloat(y));
        if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
        const below = /^(below|down|right)$/i.test((place ?? "").trim());
        return { lo: a, hi: b, below };
    }).filter(Boolean);
}

function dateFromDecimalYear(dy) {
    const y = Math.floor(dy);
    const start = Date.UTC(y, 0, 1), end = Date.UTC(y + 1, 0, 1);
    return new Date(start + (dy - y) * (end - start));
}

// data value -> pixel, from a Plot scale descriptor (works for utc + linear).
function scaleApply(scale) {
    if (scale && typeof scale.apply === "function") return (v) => scale.apply(v);
    const [d0, d1] = scale.domain, [r0, r1] = scale.range;
    const n0 = +d0, n1 = +d1;
    return (v) => r0 + ((+v - n0) / (n1 - n0)) * (r1 - r0);
}

// Place each segment's doubling label off the pixel-space center of its regression
// line, along the perpendicular (above by default, or below when the segment opts
// in) with a short leader so it clears the blobs. Drawn as raw SVG after Plot
// computes the scales, which the perpendicular needs.
function placeRegressionLabels(plotEl, segments) {
    if (typeof plotEl.scale !== "function") return;
    const svg = plotEl.tagName && plotEl.tagName.toLowerCase() === "svg"
        ? plotEl
        : plotEl.querySelector("svg");
    if (!svg) return;
    const xa = scaleApply(plotEl.scale("x"));
    const ya = scaleApply(plotEl.scale("y"));
    const NS = "http://www.w3.org/2000/svg";
    const LEADER = 15, GAP = 3;
    for (const s of segments) {
        if (!s.t) continue;
        const p0 = [xa(s.link.x1), ya(s.link.y1)];
        const p1 = [xa(s.link.x2), ya(s.link.y2)];
        const mid = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2];
        const dx = p1[0] - p0[0], dy = p1[1] - p0[1];
        const len = Math.hypot(dx, dy) || 1;
        let nx = -dy / len, ny = dx / len;   // rotate line direction 90°
        if (ny > 0) { nx = -nx; ny = -ny; }   // normalize to the upward side
        if (s.below) { nx = -nx; ny = -ny; }  // ...unless this segment opts below
        const ex = mid[0] + nx * LEADER, ey = mid[1] + ny * LEADER;
        const leader = document.createElementNS(NS, "line");
        leader.setAttribute("x1", mid[0]);
        leader.setAttribute("y1", mid[1]);
        leader.setAttribute("x2", ex);
        leader.setAttribute("y2", ey);
        leader.setAttribute("stroke", "#777");
        leader.setAttribute("stroke-width", "1");
        svg.appendChild(leader);
        // Anchor the label's near edge at the leader tip, so the leader meets the
        // side of the text (right edge when the label sits left of the line).
        const toLeft = nx < 0;
        const text = document.createElementNS(NS, "text");
        text.setAttribute("x", ex + (toLeft ? -GAP : GAP));
        text.setAttribute("y", ey);
        text.setAttribute("text-anchor", toLeft ? "end" : "start");
        text.setAttribute("dominant-baseline", "central");
        text.setAttribute("fill", "#555");
        text.setAttribute("font-size", "11");
        text.textContent = s.t;
        svg.appendChild(text);
    }
}

const MEASURES = [
    ["flux", "Cumulative fitness flux"],
    ["relative", "Relative to population average"],
];

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "16px" : "14px";
    const legendFont = mode === "slide" ? "15px" : "14px";
    let height = opts.height ?? (mode === "slide" ? 480 : 430);
    let measure = opts.measure === "relative" ? "relative" : "flux";

    const scale = colorScale(data.colors);
    const points = data.fitnessFlux.map((d) => ({
        variant: d.variant,
        date: d.date instanceof Date ? d.date : new Date(d.date),
        empFreq: d.emp_freq,
        logFitness: d.log_fitness,
    }));

    // Darker stroke per variant for the blob centerline — a stable hex string per
    // variant so Plot groups the line into one path per variant (not per point).
    const darkerByVariant = new Map(
        [...new Set(points.map((d) => d.variant))].map((v) => [
            v,
            d3.color(scale.color(v)).darker(0.6).formatHex(),
        ]),
    );

    // Per-day population mean log fitness: log( Σ x_i e^{logFitness_i} / Σ x_i ).
    // Normalized by Σ x_i since the >1% data cutoff leaves daily sums a touch under 1.
    const logBarFByDate = new Map();
    for (const [date, rows] of d3.group(points, (d) => +d.date)) {
        let num = 0, den = 0;
        for (const r of rows) {
            num += r.empFreq * Math.exp(r.logFitness);
            den += r.empFreq;
        }
        if (den > 0) logBarFByDate.set(date, Math.log(num / den));
    }
    const yOf = (d) =>
        measure === "relative" ? d.logFitness - logBarFByDate.get(+d.date) : d.logFitness;

    // Half-height tracks the current view's y-axis span (per mode), so blobs are
    // the same visual thickness in both modes and across datasets.
    const [fluxLo, fluxHi] = d3.extent(points, (d) => d.logFitness);
    const [relLo, relHi] = d3.extent(points, (d) => d.logFitness - logBarFByDate.get(+d.date));
    const spanFlux = fluxHi - fluxLo;
    const spanRel = relHi - relLo;
    const halfHeightFor = (m) =>
        HALF_HEIGHT_FRACTION * (m === "relative" ? spanRel : spanFlux);

    // Population mean log fitness f̄(t) = Σ x_v f_v (the fitness-flux trajectory),
    // fit piecewise over the ranges in opts.regressions and drawn as gray lines below.
    const meanLogFitPoints = [];
    for (const [date, rows] of d3.group(points, (d) => +d.date)) {
        let num = 0, den = 0;
        for (const r of rows) { num += r.empFreq * r.logFitness; den += r.empFreq; }
        if (den > 0) meanLogFitPoints.push([decimalYear(date), num / den]);
    }
    const segments = parseRegressions(opts.regressions)
        .map(({ lo, hi, below }) => {
            const seg = meanLogFitPoints.filter(([x]) => x >= lo && x < hi);
            const fit = linearFit(seg);
            if (!fit) return null;
            const xs = seg.map(([x]) => x);
            const x0 = Math.min(...xs), x1 = Math.max(...xs);
            const doubling = fit.slope > 0 ? Math.log(2) / fit.slope : null;
            return {
                link: {
                    x1: dateFromDecimalYear(x0), y1: fit.intercept + fit.slope * x0,
                    x2: dateFromDecimalYear(x1), y2: fit.intercept + fit.slope * x1,
                },
                t: doubling
                    ? `doubling ${doubling < 1 ? `${Math.round(doubling * 12)} mo` : `${doubling.toFixed(1)} yr`}`
                    : "",
                below,
            };
        })
        .filter(Boolean);
    // Explicit regression ranges label each line on the plot; otherwise fall back to
    // the single dominant-variant slope annotation in the header.
    const annotation = segments.length ? "" : slopeText(points);

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
        MEASURES.forEach(([key, label], i) => {
            const active = key === measure;
            const button = document.createElement("button");
            button.type = "button";
            button.textContent = label;
            Object.assign(button.style, {
                border: "none",
                borderLeft: i ? "1px solid #ccc" : "none",
                padding: "4px 11px",
                cursor: active ? "default" : "pointer",
                font: "inherit",
                background: active ? "#333" : "#fff",
                color: active ? "#fff" : "#555",
            });
            button.addEventListener("click", () => {
                if (measure === key) return;
                measure = key;
                animateToggle();
            });
            wrap.appendChild(button);
        });
        return wrap;
    }

    function draw(totalWidth) {
        lastWidth = totalWidth;
        const relative = measure === "relative";
        const halfHeight = halfHeightFor(measure);
        root.replaceChildren();

        // Header row: slope annotation on the left (flux mode only), mode toggle on
        // the right. Always rendered with the toggle so switching modes doesn't
        // shift the plot vertically.
        const header = document.createElement("div");
        Object.assign(header.style, {
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
            gap: "16px",
            margin: "0 0 4px",
        });
        const note = document.createElement("div");
        Object.assign(note.style, { fontSize: "13px", color: "#555" });
        note.textContent = !relative && annotation ? annotation : "";
        header.append(note, buildToggle());
        root.appendChild(header);

        const figRow = document.createElement("div");
        Object.assign(figRow.style, {
            display: "flex",
            alignItems: "flex-start",
            gap: `${GAP}px`,
            flexWrap: "wrap",
        });
        root.appendChild(figRow);

        const plotHolder = document.createElement("div");
        Object.assign(plotHolder.style, {
            flex: "1 1 auto",
            minWidth: "0",
            overflowX: "auto",
            marginTop: "12px",
        });

        // Measure the vertical legend to decide whether it fits beside the plot.
        let legend = buildLegend(scale, { orientation: "vertical", fontSize: legendFont });
        figRow.append(plotHolder, legend);
        const legendWidth = legend.offsetWidth;
        const sideBySide = totalWidth - legendWidth - GAP >= MIN_PLOT;

        let plotWidth;
        if (sideBySide) {
            figRow.style.flexDirection = "row";
            legend.style.marginTop = "20px"; // align legend top with the plot frame
            plotWidth = Math.max(MIN_PLOT, Math.floor(totalWidth - legendWidth - GAP));
        } else {
            figRow.style.flexDirection = "column";
            legend.remove();
            legend = buildLegend(scale, { orientation: "horizontal", fontSize: legendFont });
            legend.style.marginTop = "8px";
            figRow.appendChild(legend);
            plotWidth = Math.max(MIN_PLOT, Math.floor(totalWidth));
        }

        const plotEl = Plot.plot({
                style: { fontSize: axisFont },
                width: plotWidth,
                height,
                marginTop: 20,
                marginLeft: 58,
                marginBottom: 36,
                marginRight: 12,
                color: { type: "identity" },
                x: { type: "utc", label: null },
                y: {
                    label: relative
                        ? "Relative fitness to population average"
                        : "Cumulative fitness flux",
                    labelAnchor: "center",
                    labelArrow: "none",
                },
                marks: [
                    Plot.frame({ anchor: "left", stroke: "#333" }),
                    Plot.frame({ anchor: "bottom", stroke: "#333" }),
                    ...(relative
                        ? [Plot.ruleY([0], { stroke: "#bbb", strokeDasharray: "3,3" })]
                        : []),
                    Plot.areaY(points, {
                        x: "date",
                        y1: (d) => yOf(d) - halfHeight * d.empFreq,
                        y2: (d) => yOf(d) + halfHeight * d.empFreq,
                        z: "variant",
                        fill: (d) => scale.color(d.variant),
                        fillOpacity: 0.85,
                        curve: "basis",
                    }),
                    // Darker centerline tracing each blob's mid-point: flat in
                    // flux mode, drifting in relative mode.
                    Plot.line(points, {
                        x: "date",
                        y: (d) => yOf(d),
                        z: "variant",
                        stroke: (d) => darkerByVariant.get(d.variant),
                        strokeWidth: 1,
                        curve: "basis",
                    }),
                    // Gray piecewise regression lines over the mean-fitness trajectory
                    // (flux mode only). Their doubling labels are placed after render
                    // by placeRegressionLabels, which needs the computed scales.
                    ...(!relative
                        ? segments.map((s) =>
                              Plot.link([s.link], {
                                  x1: "x1",
                                  y1: "y1",
                                  x2: "x2",
                                  y2: "y2",
                                  stroke: "#777",
                                  strokeWidth: 1,
                              }),
                          )
                        : []),
                    Plot.tip(
                        points,
                        Plot.pointerX({
                            x: "date",
                            y: (d) => yOf(d),
                            title: (d) =>
                                `${scale.name(d.variant)}\ndate ${d3.utcFormat("%Y-%m-%d")(d.date)}\n${relative ? "relative fitness" : "fitness"} ${yOf(d).toFixed(2)}\nfrequency ${(d.empFreq * 100).toFixed(1)}%`,
                        }),
                    ),
                ],
            });
        if (!relative) placeRegressionLabels(plotEl, segments);
        plotHolder.replaceChildren(plotEl);
    }

    // Toggle modes with the blobs and centerlines gliding to their new y-positions
    // instead of snapping: rebuild in the new mode (axis/labels/legend snap), then
    // transition each area/line path's `d` from its old shape to the new one.
    function animateToggle(duration = 500) {
        const oldArea = [...root.querySelectorAll('g[aria-label="area"] path')].map((p) => p.getAttribute("d"));
        const oldLine = [...root.querySelectorAll('g[aria-label="line"] path')].map((p) => p.getAttribute("d"));
        draw(lastWidth || measureWidth());
        const areaPaths = [...root.querySelectorAll('g[aria-label="area"] path')];
        const linePaths = [...root.querySelectorAll('g[aria-label="line"] path')];
        const finalArea = areaPaths.map((p) => p.getAttribute("d"));
        const finalLine = linePaths.map((p) => p.getAttribute("d"));
        areaPaths.forEach((p, j) => oldArea[j] != null && p.setAttribute("d", oldArea[j]));
        linePaths.forEach((p, j) => oldLine[j] != null && p.setAttribute("d", oldLine[j]));
        if (areaPaths.length)
            d3.selectAll(areaPaths).data(finalArea).transition().duration(duration).attr("d", (d) => d);
        if (linePaths.length)
            d3.selectAll(linePaths).data(finalLine).transition().duration(duration).attr("d", (d) => d);
    }

    draw(measureWidth());
    const unlink = linkLegendHighlight(root, scale, { darker: darkerByVariant });

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
