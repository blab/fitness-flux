// Forecast-accuracy figure (Abousamra et al. Fig 2B adaptation): a small-multiples
// grid of MLR-vs-naive clade-frequency forecast error, one panel per virus lineage.
// Each panel plots mean absolute error against forecast lead time, with:
//
//   * a shared grey line over the hindcast (lead <= 0), where both models are the
//     in-window MLR fit and therefore coincide;
//   * two diverging lines over the forecast (lead > 0) — MLR (blue) and naive
//     (red, flat-forward from the final date);
//   * a divider + shaded band at lead 0 separating hindcast from forecast, and a
//     dashed 5% reference line;
//   * faint per-window-pair curves behind the aggregate (opts.showPairs).
//
// All panels share x and y axes for comparability. One legend serves the grid.
//
// data = {
//   panels: Array<{
//     key, label,                       // lineage id + display name (panel title)
//     n_pairs,                          // window pairs formed
//     curve: Array<{ lead, mlr, naive }>,           // aggregate MAE (%)
//     pairs: Array<{ label, points: Array<{ lead, mlr, naive }> }>  // per-pair
//   }>,
//   reference: number,                  // e.g. 5 (%)
//   models: Array<{ key: "mlr"|"naive", label, color }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, showPairs? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const SHARED = "#52514e"; // neutral ink: the hindcast fit shared by both models
const FRAME = "#333";
const GRID = "#e1e0d9";
const MUTED = "#898781";
const GAP = 18;
const MIN_PANEL = 280; // px below which the grid folds to one column

function reshapePanel(panel, models, showPairs) {
    const shared = panel.curve
        .filter((d) => d.lead <= 0 && d.mlr != null)
        .map((d) => ({ lead: d.lead, error: d.mlr }));
    const forecastLong = [];
    for (const d of panel.curve) {
        if (d.lead <= 0) continue;
        for (const m of models) {
            if (d[m.key] != null) forecastLong.push({ lead: d.lead, model: m.key, error: d[m.key] });
        }
    }
    const pairShared = [];
    const pairForecast = [];
    for (const p of showPairs ? panel.pairs ?? [] : []) {
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
    return { shared, forecastLong, pairShared, pairForecast };
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "13px" : "11px";
    const showPairs = opts.showPairs !== false;

    const panels = data.panels ?? [];
    const models = data.models ?? [
        { key: "mlr", label: "MLR", color: "#2a78d6" },
        { key: "naive", label: "Naïve", color: "#e34948" },
    ];
    const colorOf = Object.fromEntries(models.map((m) => [m.key, m.color]));
    const labelOf = Object.fromEntries(models.map((m) => [m.key, m.label]));
    const reference = data.reference ?? 5;

    // Shared axes across panels. y is taken from the AGGREGATE curves only, so a
    // sparse per-pair spike (clipped) cannot blow up the scale.
    const allLeads = panels.flatMap((p) => p.curve.map((d) => d.lead));
    const xLo = Math.min(-90, d3.min(allLeads) ?? -90);
    const xHi = Math.max(180, d3.max(allLeads) ?? 180) + 4;
    const aggMax = d3.max(panels.flatMap((p) => p.curve.flatMap((d) => [d.mlr, d.naive])).filter((v) => v != null));
    const yHi = 1.06 * Math.max(reference, aggMax ?? reference);

    const shaped = panels.map((p) => reshapePanel(p, models, showPairs));

    // --- DOM: header (legend) + grid ---
    const root = document.createElement("div");
    container.appendChild(root);

    const legend = document.createElement("div");
    Object.assign(legend.style, {
        display: "flex", flexWrap: "wrap", gap: "14px", alignItems: "center",
        font: `${mode === "slide" ? 14 : 12}px system-ui, -apple-system, "Segoe UI", sans-serif`,
        color: "#0b0b0b", margin: "2px 0 8px 4px",
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

    const grid = document.createElement("div");
    Object.assign(grid.style, { display: "grid", gap: `${GAP}px` });
    root.appendChild(grid);

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PANEL, w || 820);
    }

    function panelPlot(panel, shp, panelW, panelH, showRegionText) {
        return Plot.plot({
            style: { fontSize: axisFont, background: "transparent" },
            width: panelW,
            height: panelH,
            marginLeft: 46,
            marginRight: 12,
            marginTop: 6,
            marginBottom: 34,
            clip: true, // clip faint per-pair lines to the shared y-domain
            x: { domain: [xLo, xHi], label: "Forecast lead (days)", labelAnchor: "center", labelArrow: "none" },
            y: { domain: [0, yHi], label: "Mean absolute error (%)", labelAnchor: "center", labelArrow: "none" },
            marks: [
                Plot.rect([{ x1: 0, x2: xHi, y1: 0, y2: yHi }], {
                    x1: "x1", x2: "x2", y1: "y1", y2: "y2", fill: "#000", fillOpacity: 0.03,
                }),
                Plot.gridY({ stroke: GRID, strokeOpacity: 1 }),
                Plot.ruleX([0], { stroke: MUTED, strokeWidth: 1 }),
                Plot.ruleY([reference], { stroke: MUTED, strokeDasharray: "4,3" }),
                ...(showRegionText
                    ? [
                          Plot.text(["Hindcast"], { frameAnchor: "top-left", dx: 6, dy: 5, textAnchor: "start", lineAnchor: "top", fontSize: 10, fontStyle: "italic", fill: MUTED }),
                          Plot.text(["Forecast"], { frameAnchor: "top-right", dx: -6, dy: 5, textAnchor: "end", lineAnchor: "top", fontSize: 10, fontStyle: "italic", fill: MUTED }),
                      ]
                    : []),
                Plot.line(shp.pairShared, { x: "lead", y: "error", z: "z", stroke: SHARED, strokeOpacity: 0.12, strokeWidth: 1 }),
                Plot.line(shp.pairForecast, { x: "lead", y: "error", z: "z", stroke: (d) => colorOf[d.model], strokeOpacity: 0.12, strokeWidth: 1 }),
                Plot.line(shp.shared, { x: "lead", y: "error", stroke: SHARED, strokeWidth: 2 }),
                Plot.line(shp.forecastLong, { x: "lead", y: "error", z: "model", stroke: (d) => colorOf[d.model], strokeWidth: 2 }),
                Plot.dot(shp.shared, { x: "lead", y: "error", fill: SHARED, r: 1.6 }),
                Plot.dot(shp.forecastLong, { x: "lead", y: "error", fill: (d) => colorOf[d.model], r: 1.6 }),
                Plot.frame({ anchor: "left", stroke: FRAME }),
                Plot.frame({ anchor: "bottom", stroke: FRAME }),
                Plot.tip(
                    shp.forecastLong,
                    Plot.pointer({
                        x: "lead", y: "error",
                        title: (d) => `${labelOf[d.model]}\nlead ${Math.round(d.lead)} d\nMAE ${d.error.toFixed(1)}%`,
                    }),
                ),
            ],
        });
    }

    function panelCell(panel, shp, panelW, panelH, showRegionText) {
        const cell = document.createElement("div");
        const head = document.createElement("div");
        Object.assign(head.style, {
            font: `${mode === "slide" ? 14 : 12}px system-ui, -apple-system, "Segoe UI", sans-serif`,
            margin: "0 0 2px 46px",
        });
        const name = document.createElement("span");
        name.textContent = panel.label;
        Object.assign(name.style, { fontWeight: "600", color: "#0b0b0b" });
        const n = document.createElement("span");
        n.textContent = `  ·  n = ${(panel.pairs ?? []).length} pairs`;
        Object.assign(n.style, { color: MUTED });
        head.append(name, n);
        cell.append(head, panelPlot(panel, shp, panelW, panelH, showRegionText));
        return cell;
    }

    function draw(totalWidth) {
        const width = Math.max(MIN_PANEL, Math.floor(totalWidth));
        const cols = width >= 2 * MIN_PANEL + GAP ? 2 : 1;
        const panelW = Math.max(MIN_PANEL, Math.floor((width - (cols - 1) * GAP) / cols));
        const panelH = Math.round(Math.min(300, Math.max(200, panelW * 0.72)));
        grid.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
        grid.replaceChildren(
            ...panels.map((panel, i) => panelCell(panel, shaped[i], panelW, panelH, i === 0)),
        );
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
