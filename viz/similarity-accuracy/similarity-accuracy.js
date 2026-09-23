// Similarity-aware forecast accuracy: the same H3N2 forecasts, scored two ways.
//
// Two panels:
//
//   * forecast error in HA1 amino acids (MLR vs naive) against lead time;
//   * the identical forecasts scored by clade label, as mean absolute frequency
//     error (%), where the two models nearly converge.
//
// The panels carry different units, so they are small multiples with their own
// axes — never a dual axis.
//
// data = {
//   label, gene, horizon, n_pairs, max_distance,
//   metrics: Array<{
//     key: "aa"|"mae", label, ylabel, unit,
//     curve: Array<{ lead, mlr, naive }>,
//     pairs: Array<{ label, points: Array<{ lead, mlr, naive }> }>,
//     mlr_mean, naive_mean, gap_pct
//   }>,
//   models: Array<{ key: "mlr"|"naive", label, color }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, showPairs? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const SHARED = "#52514e"; // neutral ink: the hindcast fit both models share
const FRAME = "#333";
const GRID = "#e1e0d9";
const MUTED = "#898781";
const INK = "#0b0b0b";
const GAP = 20;
const MIN_PANEL = 290;

function fontStack(size) {
    return `${size}px system-ui, -apple-system, "Segoe UI", sans-serif`;
}

function reshape(metric, models, showPairs) {
    const shared = metric.curve
        .filter((d) => d.lead <= 0 && d.mlr != null)
        .map((d) => ({ lead: d.lead, error: d.mlr }));
    const forecast = [];
    for (const d of metric.curve) {
        if (d.lead <= 0) continue;
        for (const m of models) {
            if (d[m.key] != null) forecast.push({ lead: d.lead, model: m.key, error: d[m.key] });
        }
    }
    const pairShared = [];
    const pairForecast = [];
    for (const p of showPairs ? metric.pairs ?? [] : []) {
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
    return { shared, forecast, pairShared, pairForecast };
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "13px" : "11px";
    const headFont = mode === "slide" ? 14 : 12;
    const showPairs = opts.showPairs !== false;

    const metrics = data.metrics ?? [];
    const models = data.models ?? [
        { key: "mlr", label: "MLR", color: "#2a78d6" },
        { key: "naive", label: "Naïve", color: "#e34948" },
    ];
    const colorOf = Object.fromEntries(models.map((m) => [m.key, m.color]));
    const labelOf = Object.fromEntries(models.map((m) => [m.key, m.label]));
    const horizon = data.horizon ?? 365;

    const shaped = metrics.map((m) => reshape(m, models, showPairs));

    const allLeads = metrics.flatMap((m) => m.curve.map((d) => d.lead));
    const xLo = Math.min(-90, d3.min(allLeads) ?? -90);
    const xHi = Math.max(horizon, d3.max(allLeads) ?? horizon) + 5;

    const root = document.createElement("div");
    container.appendChild(root);

    const legend = document.createElement("div");
    Object.assign(legend.style, {
        display: "flex", flexWrap: "wrap", gap: "14px", alignItems: "center",
        font: fontStack(mode === "slide" ? 14 : 12), color: INK, margin: "2px 0 8px 4px",
    });
    const chip = (color, label, dashed) => {
        const el = document.createElement("span");
        Object.assign(el.style, { display: "inline-flex", alignItems: "center", gap: "6px" });
        const sw = document.createElement("span");
        Object.assign(sw.style, {
            width: "16px", height: "3px", borderRadius: "2px", display: "inline-block",
            background: dashed
                ? `repeating-linear-gradient(90deg, ${color} 0 4px, transparent 4px 7px)`
                : color,
        });
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
        return Math.max(MIN_PANEL, w || 880);
    }

    function heading(title, note) {
        const head = document.createElement("div");
        Object.assign(head.style, { font: fontStack(headFont), margin: "0 0 2px 54px" });
        const name = document.createElement("span");
        name.textContent = title;
        Object.assign(name.style, { fontWeight: "600", color: INK });
        const extra = document.createElement("span");
        extra.textContent = note ? `  ·  ${note}` : "";
        Object.assign(extra.style, { color: MUTED });
        head.append(name, extra);
        return head;
    }

    function metricPlot(metric, shp, panelW, panelH, annotate) {
        const values = metric.curve.flatMap((d) => [d.mlr, d.naive]).filter((v) => v != null);
        const yHi = 1.1 * (d3.max(values) ?? 1);
        const digits = 2;
        return Plot.plot({
            style: { fontSize: axisFont, background: "transparent" },
            width: panelW,
            height: panelH,
            marginLeft: 54,
            marginRight: 12,
            marginTop: 6,
            marginBottom: 36,
            clip: true,
            x: { domain: [xLo, xHi], label: "Forecast lead (days)", labelAnchor: "center", labelArrow: "none" },
            y: { domain: [0, yHi], label: metric.ylabel, labelAnchor: "center", labelArrow: "none" },
            marks: [
                Plot.rect([{ x1: 0, x2: xHi, y1: 0, y2: yHi }], {
                    x1: "x1", x2: "x2", y1: "y1", y2: "y2", fill: "#000", fillOpacity: 0.03,
                }),
                Plot.gridY({ stroke: GRID, strokeOpacity: 1 }),
                Plot.ruleX([0], { stroke: MUTED, strokeWidth: 1 }),
                ...(annotate
                    ? [
                          Plot.text(["Hindcast"], { frameAnchor: "top-left", dx: 6, dy: 5, textAnchor: "start", lineAnchor: "top", fontSize: 10, fontStyle: "italic", fill: MUTED }),
                          Plot.text(["Forecast"], { frameAnchor: "top-right", dx: -6, dy: 5, textAnchor: "end", lineAnchor: "top", fontSize: 10, fontStyle: "italic", fill: MUTED }),
                      ]
                    : []),
                Plot.line(shp.pairShared, { x: "lead", y: "error", z: "z", stroke: SHARED, strokeOpacity: 0.1, strokeWidth: 1 }),
                Plot.line(shp.pairForecast, { x: "lead", y: "error", z: "z", stroke: (d) => colorOf[d.model], strokeOpacity: 0.1, strokeWidth: 1 }),
                Plot.line(shp.shared, { x: "lead", y: "error", stroke: SHARED, strokeWidth: 2 }),
                Plot.line(shp.forecast, { x: "lead", y: "error", z: "model", stroke: (d) => colorOf[d.model], strokeWidth: 2 }),
                Plot.frame({ anchor: "left", stroke: FRAME }),
                Plot.frame({ anchor: "bottom", stroke: FRAME }),
                Plot.tip(
                    shp.forecast,
                    Plot.pointer({
                        x: "lead", y: "error",
                        title: (d) =>
                            `${labelOf[d.model]}\nlead ${Math.round(d.lead)} d\n${d.error.toFixed(digits)}${metric.unit === "%" ? "%" : " " + metric.unit}`,
                    }),
                ),
            ],
        });
    }

    function draw(totalWidth) {
        const width = Math.max(MIN_PANEL, Math.floor(totalWidth));
        const cols = width >= 2 * MIN_PANEL + GAP ? 2 : 1;
        const panelW = Math.max(MIN_PANEL, Math.floor((width - (cols - 1) * GAP) / cols));
        const panelH = Math.round(Math.min(300, Math.max(210, panelW * 0.7)));

        const cells = metrics.map((metric, i) => {
            const cell = document.createElement("div");
            const gap = metric.gap_pct;
            cell.append(
                heading(
                    metric.label,
                    gap == null ? "" : `MLR ${gap >= 0 ? "better by" : "worse by"} ${Math.abs(gap).toFixed(1)}%`,
                ),
                metricPlot(metric, shaped[i], panelW, panelH, i === 0),
            );
            return cell;
        });

        grid.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
        grid.replaceChildren(...cells);
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
