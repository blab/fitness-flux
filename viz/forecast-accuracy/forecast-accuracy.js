// Forecast-accuracy figure (Abousamra et al. Fig 2B adaptation): mean absolute
// error of clade-frequency predictions vs forecast lead time, MLR vs a naive
// (persistence) model. A single panel with:
//
//   * a shared grey line over the hindcast (lead <= 0), where both models are the
//     in-window MLR fit and therefore coincide;
//   * two diverging lines over the forecast (lead > 0) — MLR (blue) and naive
//     (red, flat-forward from the final date);
//   * a divider + shaded band at lead 0 separating hindcast from forecast, and a
//     dashed 5% reference line;
//   * faint per-window-pair curves behind the aggregate (opts.showPairs).
//
// data = {
//   curve:     Array<{ lead (days), mlr (%|null), naive (%|null) }>,   // aggregate MAE
//   pairs:     Array<{ label, points: Array<{ lead, mlr, naive }> }>,  // per-window-pair
//   reference: number,                                                 // e.g. 5 (%)
//   models:    Array<{ key: "mlr"|"naive", label, color }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, height?, showPairs? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const SHARED = "#52514e"; // neutral ink: the hindcast fit shared by both models
const FRAME = "#333";
const GRID = "#e1e0d9";
const MUTED = "#898781";
const MIN_W = 300;
const MIN_H = 220;

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    const showPairs = opts.showPairs !== false;

    const models = data.models ?? [
        { key: "mlr", label: "MLR", color: "#2a78d6" },
        { key: "naive", label: "Naïve", color: "#e34948" },
    ];
    const colorOf = Object.fromEntries(models.map((m) => [m.key, m.color]));
    const labelOf = Object.fromEntries(models.map((m) => [m.key, m.label]));
    const reference = data.reference ?? 5;

    // Hindcast: lead <= 0, where both models are the in-window fit → one grey line.
    const shared = data.curve
        .filter((d) => d.lead <= 0 && d.mlr != null)
        .map((d) => ({ lead: d.lead, error: d.mlr }));
    // Forecast: lead > 0, where the two models diverge.
    const forecastLong = [];
    for (const d of data.curve) {
        if (d.lead <= 0) continue;
        for (const m of models) {
            if (d[m.key] != null) forecastLong.push({ lead: d.lead, model: m.key, error: d[m.key] });
        }
    }
    // Faint per-pair background: hindcast shared, forecast split, same as aggregate.
    const pairShared = [];
    const pairForecast = [];
    for (const p of showPairs ? data.pairs ?? [] : []) {
        for (const pt of p.points) {
            if (pt.lead <= 0) {
                if (pt.mlr != null) pairShared.push({ z: p.label, lead: pt.lead, error: pt.mlr });
            } else {
                for (const m of models) {
                    if (pt[m.key] != null)
                        pairForecast.push({ z: `${p.label}|${m.key}`, model: m.key, lead: pt.lead, error: pt[m.key] });
                }
            }
        }
    }

    const leads = data.curve.map((d) => d.lead);
    const xLo = Math.min(-90, d3.min(leads) ?? -90);
    const xHi = Math.max(180, d3.max(leads) ?? 180) + 4;
    const allErr = data.curve.flatMap((d) => [d.mlr, d.naive]).filter((v) => v != null);
    const yHi = 1.06 * Math.max(reference, d3.max(allErr) ?? reference);

    // --- DOM: header (legend) + plot ---
    const root = document.createElement("div");
    container.appendChild(root);

    const legend = document.createElement("div");
    Object.assign(legend.style, {
        display: "flex",
        flexWrap: "wrap",
        gap: "14px",
        alignItems: "center",
        font: `${mode === "slide" ? 14 : 12}px system-ui, -apple-system, "Segoe UI", sans-serif`,
        color: "#0b0b0b",
        margin: "2px 0 6px 8px",
    });
    const chip = (color, label) => {
        const el = document.createElement("span");
        Object.assign(el.style, { display: "inline-flex", alignItems: "center", gap: "6px" });
        const sw = document.createElement("span");
        Object.assign(sw.style, { width: "16px", height: "3px", borderRadius: "2px", background: color, display: "inline-block" });
        const tx = document.createElement("span");
        tx.textContent = label;
        el.append(sw, tx);
        return el;
    };
    legend.append(
        chip(colorOf.mlr, `${labelOf.mlr} forecast`),
        chip(colorOf.naive, `${labelOf.naive} forecast`),
        chip(SHARED, "In-window fit (both)"),
    );
    root.appendChild(legend);

    const plotBox = document.createElement("div");
    root.appendChild(plotBox);

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_W, w || 760);
    }

    function draw(totalWidth) {
        const width = Math.max(MIN_W, Math.floor(totalWidth));
        const height = opts.height ?? Math.round(Math.min(420, Math.max(MIN_H, width * 0.52)));
        const fig = Plot.plot({
            style: { fontSize: axisFont, background: "transparent" },
            width,
            height,
            marginLeft: 54,
            marginRight: 16,
            marginTop: 14,
            marginBottom: 42,
            x: {
                domain: [xLo, xHi],
                label: "Forecast lead (days)",
                labelAnchor: "center",
                labelArrow: "none",
            },
            y: {
                domain: [0, yHi],
                label: "Mean absolute error (%)",
                labelAnchor: "center",
                labelArrow: "none",
            },
            marks: [
                // Faint shading of the forecast half so the eye reads hindcast vs forecast.
                Plot.rect([{ x1: 0, x2: xHi, y1: 0, y2: yHi }], {
                    x1: "x1", x2: "x2", y1: "y1", y2: "y2", fill: "#000", fillOpacity: 0.03,
                }),
                Plot.gridY({ stroke: GRID, strokeOpacity: 1 }),
                Plot.ruleX([0], { stroke: MUTED, strokeWidth: 1 }),
                Plot.ruleY([reference], { stroke: MUTED, strokeDasharray: "4,3" }),
                Plot.text([`${reference}%`], {
                    frameAnchor: "right", dx: -6, dy: -6, textAnchor: "end", lineAnchor: "bottom",
                    fontSize: 11, fill: MUTED,
                }),
                Plot.text(["Hindcast"], {
                    frameAnchor: "top-left", dx: 8, dy: 6, textAnchor: "start", lineAnchor: "top",
                    fontSize: 11, fontStyle: "italic", fill: MUTED,
                }),
                Plot.text(["Forecast"], {
                    frameAnchor: "top-right", dx: -8, dy: 6, textAnchor: "end", lineAnchor: "top",
                    fontSize: 11, fontStyle: "italic", fill: MUTED,
                }),

                // Faint per-window-pair curves behind the aggregate.
                Plot.line(pairShared, { x: "lead", y: "error", z: "z", stroke: SHARED, strokeOpacity: 0.12, strokeWidth: 1 }),
                Plot.line(pairForecast, { x: "lead", y: "error", z: "z", stroke: (d) => colorOf[d.model], strokeOpacity: 0.12, strokeWidth: 1 }),

                // Aggregate: shared hindcast, then diverging forecast lines.
                Plot.line(shared, { x: "lead", y: "error", stroke: SHARED, strokeWidth: 2 }),
                Plot.line(forecastLong, { x: "lead", y: "error", z: "model", stroke: (d) => colorOf[d.model], strokeWidth: 2 }),
                Plot.dot(shared, { x: "lead", y: "error", fill: SHARED, r: 2 }),
                Plot.dot(forecastLong, { x: "lead", y: "error", fill: (d) => colorOf[d.model], r: 2 }),

                Plot.frame({ anchor: "left", stroke: FRAME }),
                Plot.frame({ anchor: "bottom", stroke: FRAME }),

                Plot.tip(
                    forecastLong,
                    Plot.pointer({
                        x: "lead", y: "error",
                        title: (d) => `${labelOf[d.model]}\nlead ${Math.round(d.lead)} d\nMAE ${d.error.toFixed(1)}%`,
                    }),
                ),
            ],
        });
        plotBox.replaceChildren(fig);
    }

    draw(measureWidth());

    return {
        element: root,
        resize(width) {
            draw(width ? Math.max(MIN_W, Math.floor(width)) : measureWidth());
        },
        destroy() {
            root.remove();
        },
    };
}
