// Frequency–fitness phase portrait: each variant is one trajectory through
// (empirical frequency, relative fitness) space, drawn in time order. A variant
// emerges at low frequency / high relative fitness (top-left), arcs right and
// down to ~0 relative fitness near its peak frequency, then falls back down and
// left as it is outcompeted. Lines only — no blobs, no points.
//
//   • x — empirical frequency on a logit scale (spreads the low/high tails).
//   • y — log fitness relative to the daily population mean:
//         log_fitness_i − log(f̄(t)), with f̄(t) = Σ_i x_i e^{log_fitness_i} / Σ_i x_i.
//
// Same data as the time-vs-fitness component, plotted differently.
//
// data = {
//   fitnessFlux: Array<{ variant, date (ISO string), log_fitness (number), emp_freq (number) }>,
//   colors:      Array<{ variant, color, display_name, is_major, order }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, height? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";
import { colorScale, buildLegend, linkLegendHighlight } from "../lib/colors.js";

const GAP = 16; // px between plot and legend
const MIN_PLOT = 360; // px below which the legend folds under

// Logit x-axis: clamp frequency into (1%, 99%) so the transform stays finite,
// with gridline ticks at human-readable frequencies.
const FREQ_MIN = 0.01;
const FREQ_MAX = 0.99;
const LOGIT_TICK_FREQS = [0.01, 0.1, 0.5, 0.9, 0.99];
const toLogit = (p) => Math.log(p / (1 - p));
const fromLogit = (t) => 1 / (1 + Math.exp(-t));
const clampFreq = (p) => Math.min(FREQ_MAX, Math.max(FREQ_MIN, p));
const pct = d3.format(".0%");

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "16px" : "14px";
    const legendFont = mode === "slide" ? "15px" : "14px";
    let height = opts.height ?? (mode === "slide" ? 480 : 430);

    const scale = colorScale(data.colors);
    // Sort ascending by date so each variant's line is drawn in time order: the
    // trajectory loops right-then-left, so input order (not an x-sort) matters.
    const points = data.fitnessFlux
        .map((d) => ({
            variant: d.variant,
            date: d.date instanceof Date ? d.date : new Date(d.date),
            empFreq: d.emp_freq,
            logFitness: d.log_fitness,
        }))
        .sort((a, b) => a.date - b.date);

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
    const relOf = (d) => d.logFitness - logBarFByDate.get(+d.date);
    const xOf = (d) => toLogit(clampFreq(d.empFreq));

    const root = document.createElement("div");
    container.appendChild(root);

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PLOT, w || 820);
    }

    function draw(totalWidth) {
        root.replaceChildren();

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
                marginBottom: 40,
                marginRight: 12,
                color: { type: "identity" },
                x: {
                    domain: [toLogit(FREQ_MIN), toLogit(FREQ_MAX)],
                    ticks: LOGIT_TICK_FREQS.map(toLogit),
                    tickFormat: (t) => pct(fromLogit(t)),
                    label: "Clade frequency",
                    labelAnchor: "center",
                    labelArrow: "none",
                },
                y: {
                    label: "Relative fitness to population average",
                    labelAnchor: "center",
                    labelArrow: "none",
                },
                marks: [
                    Plot.frame({ anchor: "left", stroke: "#333" }),
                    Plot.frame({ anchor: "bottom", stroke: "#333" }),
                    // Population mean: variants above 0 are fitter than average.
                    Plot.ruleY([0], { stroke: "#bbb", strokeDasharray: "3,3" }),
                    Plot.line(points, {
                        x: xOf,
                        y: relOf,
                        z: "variant",
                        stroke: (d) => scale.color(d.variant),
                        strokeWidth: 1.5,
                        strokeOpacity: 0.85,
                        curve: "basis",
                    }),
                    Plot.tip(
                        points,
                        Plot.pointer({
                            x: xOf,
                            y: relOf,
                            title: (d) =>
                                `${scale.name(d.variant)}\ndate ${d3.utcFormat("%Y-%m-%d")(d.date)}\nfrequency ${(d.empFreq * 100).toFixed(1)}%\nrelative fitness ${relOf(d).toFixed(2)}`,
                        }),
                    ),
                ],
            }),
        );
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
