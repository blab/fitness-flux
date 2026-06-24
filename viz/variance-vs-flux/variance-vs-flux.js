// The "fitness wave" figure: a 2×2 grid of four panels built from the daily
// fitness-wave timeseries (population fitness variance and flux/velocity per day):
//
//   [ variance through time ]     [ variance vs flux (Fisher) + best fit ]
//   [ flux through time     ]     [ per-year flux points + yearly means   ]
//
// Variance and flux are shown ×1000 ("× 10⁻³"). The scatter illustrates Fisher's
// fundamental theorem (the rate of fitness gain tracks fitness variance); the
// best-fit slope/r come from the precomputed regression.
//
// Port of the four-panel figure in fitness-flux-analysis/fitness-flux.nb
// (figVariance / figVarianceVsVelocity / figVelocity / figYearlyMeans).
//
// data = {
//   points: Array<{ date (ISO string), variance (number), velocity (number|null) }>,
//   fit:    { slope, intercept, r_squared, n }
// }
// opts = { mode?: "inline"|"slide"|"dashboard", width?, height? }
//
// Pure: no fetching, no ResizeObserver. The host owns data loading and resize.
// Returns { element, resize(width?), destroy() }.

import * as Plot from "../lib/plot.js";
import * as d3 from "../lib/d3.js";

const COL = "#4c90c0"; // notebook's ColorData["Rainbow"][0.3]
const GAP = 16;
const MIN_PANEL = 300; // px below which the grid folds to one column

// Deterministic horizontal jitter in [-0.3, 0.3] from the date string, so the
// per-year point cloud is stable across redraws/resizes (no Math.random).
function jitter(dateStr) {
    let h = 0;
    for (let i = 0; i < dateStr.length; i++) h = (h * 31 + dateStr.charCodeAt(i)) >>> 0;
    return ((h % 1000) / 1000 - 0.5) * 0.6;
}

export function render(container, data, opts = {}) {
    const mode = opts.mode ?? "inline";
    const axisFont = mode === "slide" ? "14px" : "12px";
    // Optional cap on the shared log-fitness value scale (variance + flux, ×10⁻³),
    // applied to every value axis; passed via the figure directive as `scalemax=N`.
    const yMax = typeof opts.scalemax === "number" ? opts.scalemax : null;

    const fit = data.fit ?? {};
    const pts = data.points.map((p) => ({
        dateStr: p.date,
        date: p.date instanceof Date ? p.date : new Date(p.date),
        year: +String(p.date).slice(0, 4),
        variance: p.variance * 1000,
        velocity: typeof p.velocity === "number" ? p.velocity * 1000 : null,
    }));
    const velPts = pts.filter((p) => p.velocity != null);
    // Shared x-domain for the two timeseries panels so their axes align (velocity
    // starts later than variance, after the velocity window).
    const dateExtent = d3.extent(pts, (p) => p.date);
    // Capped value domains (top bounded by scalemax; bottom kept from the data).
    const varDom = yMax == null ? null : [Math.min(0, d3.min(pts, (p) => p.variance) ?? 0), yMax];
    const fluxDom = yMax == null ? null : [Math.min(0, d3.min(velPts, (p) => p.velocity) ?? 0), yMax];

    // Square, shared domain for the scatter so the variance↔flux relationship
    // reads against the diagonal.
    const hi =
        Math.max(d3.max(velPts, (p) => p.variance), d3.max(velPts, (p) => p.velocity)) * 1.03;
    const lo = Math.min(0, d3.min(velPts, (p) => p.velocity));
    const sq = [lo, yMax == null ? hi : yMax];
    // Best-fit line across the domain: y = slope·x + 1000·intercept (in ×10⁻³ units).
    const fitLine = [
        { x: sq[0], y: fit.slope * sq[0] + (fit.intercept ?? 0) * 1000 },
        { x: sq[1], y: fit.slope * sq[1] + (fit.intercept ?? 0) * 1000 },
    ];
    const annotation = `slope = ${(+fit.slope).toFixed(2)}\nr = ${Math.sqrt(fit.r_squared ?? 0).toFixed(2)}`;

    // Per-year flux means: a horizontal segment at each year's mean, plus a label.
    const byYear = d3.group(velPts, (p) => p.year);
    const years = [...byYear.keys()].sort((a, b) => a - b);
    const meanSegs = [];
    const meanLabels = [];
    for (const y of years) {
        const m = d3.mean(byYear.get(y), (p) => p.velocity);
        meanSegs.push({ x: y - 0.3, y: m, year: y }, { x: y + 0.3, y: m, year: y });
        meanLabels.push({ year: y, y: m, label: m.toFixed(1) });
    }
    const yearlyPts = velPts.map((p) => ({
        x: p.year + jitter(p.dateStr),
        y: p.velocity,
        dateStr: p.dateStr,
    }));

    const VAR_LABEL = "Log fitness variance (× 10⁻³)";
    const FLUX_LABEL = "Log fitness flux\nper gen (× 10⁻³)";
    const fmtDate = d3.utcFormat("%Y-%m-%d");

    const root = document.createElement("div");
    container.appendChild(root);
    const grid = document.createElement("div");
    Object.assign(grid.style, { display: "grid", gap: `${GAP}px`, marginTop: "8px" });
    root.appendChild(grid);

    function measureWidth() {
        const w = opts.width ?? Math.floor(container.clientWidth);
        return Math.max(MIN_PANEL, w || 820);
    }

    const base = (panelW, panelH, y) => ({
        style: { fontSize: axisFont },
        width: panelW,
        height: panelH,
        marginLeft: 56,
        marginRight: 14,
        marginTop: 16,
        // Uniform across all four panels so their data regions are the same size
        // and axes line up; sized for the scatter's x-axis label below its ticks.
        marginBottom: 46,
        clip: yMax != null,
        y: { labelAnchor: "center", labelArrow: "none", ...y },
    });

    function timeseriesPanel(panelW, panelH, points, accessor, label, what, yDom) {
        return Plot.plot({
            ...base(panelW, panelH, yDom ? { label, domain: yDom } : { label }),
            x: { type: "utc", domain: dateExtent, label: null },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.dot(points, { x: "date", y: accessor, fill: COL, r: 2, fillOpacity: 0.5 }),
                Plot.tip(
                    points,
                    Plot.pointerX({
                        x: "date",
                        y: accessor,
                        title: (d) => `${fmtDate(d.date)}\n${what} ${accessor(d).toFixed(2)}`,
                    }),
                ),
            ],
        });
    }

    function scatterPanel(panelW, panelH) {
        return Plot.plot({
            ...base(panelW, panelH, { label: FLUX_LABEL, domain: sq }),
            x: { label: VAR_LABEL, domain: sq, labelAnchor: "center", labelArrow: "none" },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.dot(velPts, {
                    x: "variance",
                    y: "velocity",
                    fill: COL,
                    r: 2,
                    fillOpacity: 0.5,
                }),
                Plot.line(fitLine, { x: "x", y: "y", stroke: "#000", strokeDasharray: "4,3" }),
                Plot.text([annotation], {
                    frameAnchor: "top-left",
                    dx: 6,
                    dy: 6,
                    textAnchor: "start",
                    lineAnchor: "top",
                    fontSize: 11,
                }),
                Plot.tip(
                    velPts,
                    Plot.pointer({
                        x: "variance",
                        y: "velocity",
                        title: (d) =>
                            `${fmtDate(d.date)}\nvariance ${d.variance.toFixed(2)}\nflux ${d.velocity.toFixed(2)}`,
                    }),
                ),
            ],
        });
    }

    function yearlyPanel(panelW, panelH) {
        return Plot.plot({
            ...base(panelW, panelH, fluxDom ? { label: FLUX_LABEL, domain: fluxDom } : { label: FLUX_LABEL }),
            x: {
                domain: [years[0] - 0.5, years[years.length - 1] + 0.5],
                ticks: years,
                tickFormat: (y) => String(y).slice(-2),
                label: null,
            },
            marks: [
                Plot.frame({ anchor: "left", stroke: "#333" }),
                Plot.frame({ anchor: "bottom", stroke: "#333" }),
                Plot.dot(yearlyPts, { x: "x", y: "y", fill: COL, r: 1.8, fillOpacity: 0.3 }),
                Plot.line(meanSegs, { x: "x", y: "y", z: "year", stroke: "#000", strokeWidth: 2 }),
                Plot.text(meanLabels, {
                    x: "year",
                    y: "y",
                    text: "label",
                    dy: -8,
                    fontSize: 12,
                    fontWeight: "bold",
                    fill: "#000",
                }),
            ],
        });
    }

    function draw(totalWidth) {
        const cols = totalWidth >= 2 * MIN_PANEL + GAP ? 2 : 1;
        const panelW = Math.max(MIN_PANEL, Math.floor((totalWidth - (cols - 1) * GAP) / cols));
        const panelH = Math.round(Math.min(240, Math.max(160, panelW * 0.6)));
        grid.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
        grid.replaceChildren(
            timeseriesPanel(panelW, panelH, pts, (d) => d.variance, VAR_LABEL, "variance", varDom),
            scatterPanel(panelW, panelH),
            timeseriesPanel(panelW, panelH, velPts, (d) => d.velocity, FLUX_LABEL, "flux", fluxDom),
            yearlyPanel(panelW, panelH),
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
