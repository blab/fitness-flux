// Lineage-deltas multiple-regression view: a predicted-vs-observed scatter beside
// a coefficient table. The OLS fits the per-branch change in (log) fitness on the
// change in substitution count across four non-overlapping regions — spike RBD,
// spike S1 outside the RBD, ORF1ab, accessory — so each term's partial estimate
// strips out the confounding in the marginal scatters (a more-evolved lineage
// gains substitutions everywhere at once). The scatter plots each branch's model
// prediction (x) against its observed change in fitness (y); for in-sample OLS
// fitted values the calibration line is exactly y = x.
//
// An All / Early / Late toggle (matching the sibling lineage-deltas component)
// switches the active fit: the table shows that group's coefficients and the
// scatter shows that group's branches with predictions recomputed from those
// coefficients, so table and scatter always agree.
//
// data = {
//   seasons:      Array<string>,                         // color domain (early→late)
//   linear_model: { <group>: { terms: [{ term, key, label, estimate, se, t, p }],
//                              r_squared, n } },          // group in {all, early, late}
//   points:       Array<{ timepoint, period, parent, child, delta_log_fitness,
//                          values: { rbd, s1_nonrbd, orf1ab, accessory, ... } }>
// }
// opts = { mode?, width?, height?, initial? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?, height?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const GAP = 16; // px between scatter and table
const TABLE_W = 230; // px table column width when laid out side by side
const TWO_COL_MIN = 546; // = TABLE_W + GAP + 300 min scatter; below this, panels stack
const DOT_OPACITY = 0.55;
const GROUP_LABELS = [
    ["all", "All"],
    ["early", "Early"],
    ["late", "Late"],
];

// Same blue(early)→red(late) ramp as the sibling lineage-deltas component.
const RAINBOW = [
    [0.0, "#781c87"],
    [0.3, "#4c90c0"],
    [0.5, "#5cba4f"],
    [0.7, "#dbc740"],
    [0.85, "#ef8a2b"],
    [1.0, "#d62f26"],
];

function padExtent([lo, hi], frac = 0.05) {
    if (lo === hi) return [lo - 1, hi + 1];
    const pad = (hi - lo) * frac;
    return [lo - pad, hi + pad];
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

// Fixed 4 decimal places so every estimate shares the same precision.
function fmtEstimate(v) {
    if (v == null || !Number.isFinite(v)) return "—";
    return v.toFixed(4);
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    let height = opts.height ?? (mode === "slide" ? 340 : 300);

    // Normalize the model map: accept either { group: {...} } or a bare { terms, ... }.
    const rawModels = data.linear_model ?? {};
    const models = rawModels.terms ? { all: rawModels } : rawModels;
    const states = GROUP_LABELS.filter(([k]) => models[k]);
    if (states.length === 0) states.push(["all", "All"]);
    let state = states.some(([k]) => k === opts.initial) ? opts.initial : states[0][0];

    const predictorTerms = (m) => (m?.terms ?? []).filter((t) => t.key);
    const interceptOf = (m) =>
        (m?.terms ?? []).find((t) => t.key == null)?.estimate ?? 0;
    const predictFor = (point, group) => {
        const m = models[group];
        let v = interceptOf(m);
        for (const t of predictorTerms(m)) v += t.estimate * (point.values[t.key] ?? 0);
        return v;
    };

    const rainbow = d3
        .scaleLinear()
        .domain(RAINBOW.map((d) => d[0]))
        .range(RAINBOW.map((d) => d[1]))
        .clamp(true);
    const seasonIndex = new Map((data.seasons ?? []).map((s, i) => [s, i]));
    const seasonSpan = Math.max(1, (data.seasons ?? []).length - 1);
    const colorOf = (p) => rainbow((seasonIndex.get(p.timepoint) ?? 0) / seasonSpan);

    // Need every predictor key the "all" model uses to recompute fitted values.
    const requiredKeys = predictorTerms(models.all ?? models[state]).map((t) => t.key);
    const pts = (data.points ?? [])
        .filter(
            (p) =>
                typeof p.delta_log_fitness === "number" &&
                requiredKeys.every((k) => typeof p.values?.[k] === "number"),
        )
        .map((p) => ({
            y: p.delta_log_fitness,
            values: p.values ?? {},
            period: p.period,
            timepoint: p.timepoint,
            parent: p.parent,
            child: p.child,
            color: colorOf(p),
        }));

    const activePoints = (group) =>
        group === "all" ? pts : pts.filter((p) => p.period === group);

    // Fixed domains across every toggle state so axes never jump. Pad x (predicted)
    // and y (observed) independently — as the lineage-deltas panels do — so the
    // x-axis fits the prediction range instead of the wider observed spread.
    const xVals = [0];
    const yVals = [0];
    for (const p of pts) {
        yVals.push(p.y);
        for (const [g] of states) {
            if (g === "all" || p.period === g) xVals.push(predictFor(p, g));
        }
    }
    const xDomain = padExtent(d3.extent(xVals));
    const yDomain = padExtent(d3.extent(yVals));

    const root = document.createElement("div");
    container.appendChild(root);
    let lastWidth = 0;

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(300, w || 820);
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
        first.textContent = data.seasons?.[0] ?? "";
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
        last.textContent = data.seasons?.[data.seasons.length - 1] ?? "";
        wrap.append(first, bar, last);
        return wrap;
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
        states.forEach(([key, label], i) => {
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
                draw(lastWidth || measureWidth());
            });
            wrap.appendChild(button);
        });
        return wrap;
    }

    function buildTable(group) {
        const m = models[group] ?? { terms: [], r_squared: null, n: null };
        const wrap = document.createElement("div");
        const table = document.createElement("table");
        table.className = "ldm-table";
        // Override press.css `table { margin: 1.5rem 0 }` so the table top-aligns with
        // the scatter, and shrink the font slightly.
        Object.assign(table.style, {
            marginTop: "0",
            marginBottom: "0",
            fontSize: "0.82rem",
        });

        const thead = document.createElement("thead");
        const headRow = document.createElement("tr");
        const headers = [
            ["Predictor", ""],
            ["Estimate", "num"],
        ];
        for (const [text, cls] of headers) {
            const th = document.createElement("th");
            th.textContent = text;
            if (cls) th.className = cls;
            headRow.appendChild(th);
        }
        thead.appendChild(headRow);
        table.appendChild(thead);

        const tbody = document.createElement("tbody");
        for (const term of m.terms ?? []) {
            const tr = document.createElement("tr");
            const cells = [
                [term.label ?? term.term, ""],
                [fmtEstimate(term.estimate), "num"],
            ];
            for (const [text, cls] of cells) {
                const td = document.createElement("td");
                td.textContent = text;
                if (cls) td.className = cls;
                tr.appendChild(td);
            }
            tbody.appendChild(tr);
        }
        table.appendChild(tbody);

        wrap.append(table);
        return wrap;
    }

    function buildScatter(group, panelW, panelH) {
        const m = models[group];
        const shown = activePoints(group).map((p) => ({
            x: predictFor(p, group),
            y: p.y,
            color: p.color,
            timepoint: p.timepoint,
            parent: p.parent,
            child: p.child,
        }));
        const lineLo = Math.max(xDomain[0], yDomain[0]);
        const lineHi = Math.min(xDomain[1], yDomain[1]);
        const line = [
            { x: lineLo, y: lineLo },
            { x: lineHi, y: lineHi },
        ];
        const rVal = pearson(shown.map((p) => p.x), shown.map((p) => p.y));
        const r2 = m?.r_squared;
        const statRows = [
            `r = ${rVal == null ? "n/a" : rVal.toFixed(2)}`,
            `R² = ${r2 == null ? "n/a" : r2.toFixed(2)}`,
            `n = ${m?.n ?? shown.length}`,
        ];
        const statFont = mode === "slide" ? 15 : 13;
        const statLineH = statFont * 1.5;

        const plot = Plot.plot({
            style: { fontSize: axisFont },
            width: panelW,
            height: panelH,
            marginLeft: 56,
            marginRight: 12,
            marginTop: 8,
            marginBottom: 42,
            color: { type: "identity" },
            x: {
                domain: xDomain,
                label: "Predicted change in fitness",
                labelAnchor: "center",
                labelArrow: "none",
            },
            y: {
                domain: yDomain,
                label: "Observed change in fitness",
                labelAnchor: "center",
                labelArrow: "none",
            },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.line(line, {
                    x: "x",
                    y: "y",
                    stroke: "#999",
                    strokeDasharray: "4,3",
                }),
                Plot.dot(shown, {
                    x: "x",
                    y: "y",
                    r: 3,
                    fill: "color",
                    fillOpacity: DOT_OPACITY,
                    stroke: "none",
                }),
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
                            `${d.parent} → ${d.child}\nseason ${d.timepoint}\npredicted ${d.x.toFixed(3)}\nobserved ${d.y.toFixed(3)}`,
                    }),
                ),
            ],
        });
        // Italicize the leading r / R² / n symbol of each stat line.
        for (const t of plot.querySelectorAll('g[aria-label="text"] text')) {
            const s = t.textContent;
            if (!/^(R²|r|n) = /.test(s)) continue;
            const cut = s.indexOf(" = ");
            const sym = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
            sym.setAttribute("font-style", "italic");
            sym.textContent = s.slice(0, cut);
            const rest = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
            rest.textContent = s.slice(cut);
            t.replaceChildren(sym, rest);
        }
        return plot;
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
        header.append(buildTimeLegend());
        if (states.length > 1) header.append(buildToggle());
        root.appendChild(header);

        // CSS grid (not wrapping flex): scatter in a shrinkable left column,
        // table fixed-width on the right. Grid never folds the table below even
        // when the two columns sum to exactly the container width.
        const twoCol = totalWidth >= TWO_COL_MIN;
        const grid = document.createElement("div");
        Object.assign(grid.style, {
            display: "grid",
            gridTemplateColumns: twoCol ? `minmax(0, 1fr) ${TABLE_W}px` : "minmax(0, 1fr)",
            gap: `${GAP}px`,
            alignItems: "start",
            marginTop: "8px",
        });

        const scatterW = twoCol
            ? Math.max(300, totalWidth - TABLE_W - GAP)
            : totalWidth;
        const panelH = Math.round(Math.min(height, Math.max(240, scatterW * 0.8)));
        const scatterWrap = document.createElement("div");
        scatterWrap.style.minWidth = "0";
        scatterWrap.appendChild(buildScatter(state, scatterW, panelH));

        grid.append(scatterWrap, buildTable(state));
        root.appendChild(grid);
    }

    draw(measureWidth());

    return {
        element: root,
        resize(width, newHeight) {
            if (newHeight) height = newHeight;
            draw(width ? Math.max(300, Math.floor(width)) : measureWidth());
        },
        destroy() {
            root.remove();
        },
    };
}
