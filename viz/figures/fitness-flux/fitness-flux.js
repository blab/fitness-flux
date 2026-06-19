// Fitness flux through time: each variant is a horizontal violin centred on its
// scaffolded log fitness, with half-height proportional to frequency. A traveling
// wave of fitness as the population explores successive adaptive variants.
//
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (renderFitnessFlux,
// drawFitnessFlux) on 2026-06-19.
//
// data = {
//   fitnessFlux: Array<{ variant, date (ISO string), fit (number), freq (number) }>,
//   colors:      Array<{ variant, color, display_name, is_major, order }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, height? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../../lib/plot.js";
import * as d3 from "../../lib/d3.js";
import { colorScale, buildLegend } from "../../lib/colors.js";

const FLUX_HALF_HEIGHT = 0.13; // fitness units at frequency 1
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

function slopeText(points) {
    const dominant = [];
    for (const [, rows] of d3.group(points, (d) => +d.date)) {
        const top = d3.greatest(rows, (d) => d.freq);
        if (top && top.freq > 0.5) dominant.push([decimalYear(top.date), top.fit]);
    }
    const slope = linearFitSlope(dominant);
    return slope && slope > 0
        ? `slope = ${slope.toFixed(2)} per year   ·   doubling = ${(Math.log(2) / slope).toFixed(1)} years`
        : "";
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "16px" : "14px";
    const legendFont = mode === "slide" ? "15px" : "14px";
    const height = opts.height ?? (mode === "slide" ? 480 : 430);

    const scale = colorScale(data.colors);
    const points = data.fitnessFlux.map((d) => ({
        variant: d.variant,
        date: d.date instanceof Date ? d.date : new Date(d.date),
        freq: d.freq,
        fit: d.fit,
    }));
    const annotation = slopeText(points);

    const root = document.createElement("div");
    container.appendChild(root);

    function measure() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PLOT, w || 820);
    }

    function draw(totalWidth) {
        root.replaceChildren();

        if (annotation) {
            const note = document.createElement("div");
            Object.assign(note.style, { fontSize: "13px", color: "#555", margin: "4px 0" });
            note.textContent = annotation;
            root.appendChild(note);
        }

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
                y: { label: "Cumulative fitness flux", labelAnchor: "center", labelArrow: "none" },
                marks: [
                    Plot.frame({ anchor: "left", stroke: "#333" }),
                    Plot.frame({ anchor: "bottom", stroke: "#333" }),
                    Plot.areaY(points, {
                        x: "date",
                        y1: (d) => d.fit - FLUX_HALF_HEIGHT * d.freq,
                        y2: (d) => d.fit + FLUX_HALF_HEIGHT * d.freq,
                        z: "variant",
                        fill: (d) => scale.color(d.variant),
                        fillOpacity: 0.85,
                        curve: "basis",
                    }),
                    Plot.tip(
                        points,
                        Plot.pointerX({
                            x: "date",
                            y: (d) => d.fit,
                            title: (d) =>
                                `${scale.name(d.variant)}\ndate ${d3.utcFormat("%Y-%m-%d")(d.date)}\nfitness ${d.fit.toFixed(2)}\nfrequency ${(d.freq * 100).toFixed(1)}%`,
                        }),
                    ),
                ],
            }),
        );
    }

    draw(measure());

    return {
        element: root,
        resize(width) {
            draw(width ? Math.max(MIN_PLOT, Math.floor(width)) : measure());
        },
        destroy() {
            root.remove();
        },
    };
}
