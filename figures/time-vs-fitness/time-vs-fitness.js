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

    const annotation = slopeText(points);

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

        plotHolder.replaceChildren(
            Plot.plot({
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
            }),
        );
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
