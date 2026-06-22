// Variant frequencies through time: per-season panels comparing empirical weekly
// frequencies (dots) with MLR-modeled frequencies (lines) for each variant. A
// "Logit" switch toggles the y-axis between a linear frequency scale and a logit
// transform (log(p/(1-p))), which spreads out the low- and high-frequency tails.
//
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (renderFrequencyPanels)
// on 2026-06-19.
//
// data = {
//   seasonal: Array<{ timepoint, date (ISO string), variant, empirical (number|null), modeled (number|null) }>,
//   colors:   Array<{ variant, color, display_name, is_major, order }>
// }
// opts = { mode?: "inline"|"slide"|"dashboard", logit?: boolean, width? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";
import { colorScale, buildLegend } from "../lib/colors.js";

const COL_GAP = 10; // px between panel columns
const YEAR_MS = 365.25 * 24 * 3600 * 1000;

// Logit y-axis: clamp frequency into (1%, 99%) so the transform stays finite,
// with gridline ticks at human-readable frequencies.
const FREQ_MIN = 0.01;
const FREQ_MAX = 0.99;
const LOGIT_TICK_FREQS = [0.01, 0.1, 0.5, 0.9, 0.99];
const toLogit = (p) => Math.log(p / (1 - p));
const fromLogit = (t) => 1 / (1 + Math.exp(-t));
const clampFreq = (p) => Math.min(FREQ_MAX, Math.max(FREQ_MIN, p));
const pct = d3.format(".0%");

// Seasons span ~1 year (SARS-CoV-2) or ~2 years (H3N2); tick every 6 months or
// yearly accordingly. Derived from the data so the component is dataset-agnostic.
function tickIntervalFor(dated) {
    const spans = [];
    for (const [, rows] of d3.group(dated, (d) => d.timepoint)) {
        const [lo, hi] = d3.extent(rows, (d) => d.date);
        if (lo && hi) spans.push(hi - lo);
    }
    spans.sort((a, b) => a - b);
    const median = spans.length ? spans[Math.floor(spans.length / 2)] : 0;
    return median > 1.4 * YEAR_MS ? d3.utcYear : d3.utcMonth.every(6);
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const panelFont = mode === "slide" ? "14px" : "12px";
    const legendFont = mode === "slide" ? "15px" : "14px";
    const panelW = mode === "slide" ? 300 : 232;
    const panelH = mode === "slide" ? 190 : 150;

    let useLogit = opts.logit === true;
    let lastWidth = 0;

    const scale = colorScale(data.colors);
    const dated = data.seasonal.map((d) => ({
        ...d,
        date: d.date instanceof Date ? d.date : new Date(d.date),
    }));
    const tickInterval = tickIntervalFor(dated);
    const bySeason = d3.group(dated, (d) => d.timepoint);
    const seasons = [...bySeason.keys()].sort();

    const root = document.createElement("div");
    container.appendChild(root);

    // Controls row: variant legend (left) + Logit switch (right). Built once so
    // toggling the y-axis only redraws the panels, not this row.
    const header = document.createElement("div");
    Object.assign(header.style, {
        display: "flex",
        alignItems: "flex-start",
        justifyContent: "space-between",
        gap: "16px",
        margin: "8px 0 4px",
    });
    const legend = buildLegend(scale, { orientation: "horizontal", fontSize: legendFont });
    Object.assign(legend.style, { flex: "1 1 auto", minWidth: "0", margin: "0" });
    header.append(legend, buildLogitSwitch());
    root.appendChild(header);

    const grid = document.createElement("div");
    Object.assign(grid.style, { display: "grid", gap: `6px ${COL_GAP}px`, marginTop: "8px" });
    root.appendChild(grid);

    function buildLogitSwitch() {
        const wrap = document.createElement("span");
        Object.assign(wrap.style, {
            display: "inline-flex",
            alignItems: "center",
            gap: "6px",
            fontSize: "13px",
            color: "#555",
            cursor: "pointer",
            userSelect: "none",
            flex: "0 0 auto",
            marginTop: "2px",
        });
        const label = document.createElement("span");
        label.textContent = "Logit";
        const track = document.createElement("span");
        Object.assign(track.style, {
            position: "relative",
            width: "30px",
            height: "16px",
            borderRadius: "8px",
            background: "#ccc",
            transition: "background 0.15s",
            flex: "0 0 auto",
        });
        const knob = document.createElement("span");
        Object.assign(knob.style, {
            position: "absolute",
            top: "2px",
            left: "2px",
            width: "12px",
            height: "12px",
            borderRadius: "50%",
            background: "#fff",
            boxShadow: "0 1px 2px rgba(0, 0, 0, 0.3)",
            transition: "left 0.15s",
        });
        track.appendChild(knob);
        wrap.append(label, track);
        const sync = () => {
            track.style.background = useLogit ? "#333" : "#ccc";
            knob.style.left = useLogit ? "16px" : "2px";
        };
        wrap.setAttribute("role", "switch");
        wrap.addEventListener("click", () => {
            useLogit = !useLogit;
            wrap.setAttribute("aria-checked", String(useLogit));
            sync();
            animateToggle();
        });
        wrap.setAttribute("aria-checked", String(useLogit));
        sync();
        return wrap;
    }

    function panelFor(season) {
        const rows = bySeason.get(season);
        const tipPoints = [];
        for (const d of rows) {
            if (typeof d.modeled === "number")
                tipPoints.push({ variant: d.variant, date: d.date, value: d.modeled, kind: "MLR" });
            if (typeof d.empirical === "number")
                tipPoints.push({ variant: d.variant, date: d.date, value: d.empirical, kind: "empirical" });
        }
        const yOf = (v) => (useLogit ? toLogit(clampFreq(v)) : v);
        const y = useLogit
            ? {
                  domain: [toLogit(FREQ_MIN), toLogit(FREQ_MAX)],
                  ticks: LOGIT_TICK_FREQS.map(toLogit),
                  tickFormat: (t) => pct(fromLogit(t)),
                  label: null,
              }
            : { domain: [0, 1], ticks: [0, 0.5, 1], tickFormat: ".0%", label: null };
        return Plot.plot({
            style: { fontSize: panelFont },
            width: panelW,
            height: panelH,
            marginLeft: 42,
            marginBottom: 22,
            marginTop: 6,
            marginRight: 28,
            color: { type: "identity" },
            x: { type: "utc", ticks: tickInterval, tickFormat: "%b %Y", nice: tickInterval, label: null },
            y,
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.line(rows.filter((d) => typeof d.modeled === "number"), {
                    x: "date",
                    y: (d) => yOf(d.modeled),
                    z: "variant",
                    stroke: (d) => scale.color(d.variant),
                    strokeWidth: 1.3,
                }),
                Plot.dot(rows.filter((d) => typeof d.empirical === "number"), {
                    x: "date",
                    y: (d) => yOf(d.empirical),
                    fill: (d) => scale.color(d.variant),
                    r: 1,
                    fillOpacity: 0.5,
                }),
                Plot.tip(
                    tipPoints,
                    Plot.pointer({
                        x: "date",
                        y: (d) => yOf(d.value),
                        title: (d) => `${scale.name(d.variant)}\n${d.kind} ${(d.value * 100).toFixed(1)}%`,
                    }),
                ),
            ],
        });
    }

    function draw(totalWidth) {
        lastWidth = totalWidth;
        const cols = Math.max(1, Math.floor((totalWidth + COL_GAP) / (panelW + COL_GAP)));
        grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        grid.replaceChildren(...seasons.map(panelFor));
    }

    function measure() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(panelW, w || panelW * 4);
    }

    // Swap a panel to the new y-mode with the marks gliding: render the new panel,
    // start its dots/lines at the old panel's positions, replace (the axis snaps to
    // the new scale), then transition cy / path-d to their final positions.
    function animatePanelSwap(oldPanel, newPanel, duration = 500) {
        const oldDots = oldPanel.querySelectorAll('g[aria-label="dot"] circle');
        const oldLines = oldPanel.querySelectorAll('g[aria-label="line"] path');
        const newDots = [...newPanel.querySelectorAll('g[aria-label="dot"] circle')];
        const newLines = [...newPanel.querySelectorAll('g[aria-label="line"] path')];
        const finalCy = newDots.map((c) => c.getAttribute("cy"));
        const finalD = newLines.map((p) => p.getAttribute("d"));
        newDots.forEach((c, j) => oldDots[j] && c.setAttribute("cy", oldDots[j].getAttribute("cy")));
        newLines.forEach((p, j) => oldLines[j] && p.setAttribute("d", oldLines[j].getAttribute("d")));
        oldPanel.replaceWith(newPanel);
        if (newDots.length)
            d3.selectAll(newDots).data(finalCy).transition().duration(duration).attr("cy", (d) => d);
        if (newLines.length)
            d3.selectAll(newLines).data(finalD).transition().duration(duration).attr("d", (d) => d);
    }

    // Toggle the y-axis with points/lines gliding to their new positions.
    function animateToggle() {
        const oldPanels = [...grid.children];
        const n = Math.min(oldPanels.length, seasons.length);
        for (let i = 0; i < n; i++) animatePanelSwap(oldPanels[i], panelFor(seasons[i]));
    }

    draw(measure());

    return {
        element: root,
        resize(width) {
            draw(width ? Math.max(panelW, Math.floor(width)) : measure());
        },
        destroy() {
            root.remove();
        },
    };
}
