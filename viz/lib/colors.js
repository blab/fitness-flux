// Variant color layer shared by the figure components.
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (buildColors, variantColor,
// majorLegend) on 2026-06-19.
//
// A "colors" array is a list of { variant, color, display_name, is_major, order }
// rows (the JSON form of {dataset}_colors.tsv).

// Auspice-style legend hover: dim everything but the hovered variant.
const DIM = 0.3;
const DURATION = 350; // ms, matches Auspice's fastTransitionDuration

export function colorScale(colors) {
    const colorByVariant = new Map(colors.map((d) => [d.variant, d.color]));
    const displayName = new Map(colors.map((d) => [d.variant, d.display_name]));
    const majorVariants = colors
        .filter((d) => d.is_major === true || `${d.is_major}` === "true")
        .sort((a, b) => a.order - b.order)
        .map((d) => d.variant);
    return {
        majorVariants,
        color: (v) => colorByVariant.get(v) ?? "#808080",
        name: (v) => displayName.get(v) ?? v,
    };
}

// Build a self-styled legend of the major variants (no external CSS, so it
// survives a freeze). Orientation "vertical" stacks; "horizontal" wraps inline.
export function buildLegend(scale, { orientation = "vertical", fontSize = "14px", swatch = 13 } = {}) {
    const vertical = orientation === "vertical";
    const div = document.createElement("div");
    Object.assign(div.style, {
        display: "flex",
        flexDirection: vertical ? "column" : "row",
        flexWrap: vertical ? "nowrap" : "wrap",
        gap: vertical ? "4px 0" : "4px 14px",
        alignItems: vertical ? "flex-start" : "center",
        fontSize,
        lineHeight: "1.3",
    });
    for (const v of scale.majorVariants) {
        const item = document.createElement("span");
        item.dataset.variant = v;
        Object.assign(item.style, {
            display: "inline-flex",
            alignItems: "center",
            gap: "5px",
            whiteSpace: "nowrap",
            cursor: "pointer",
        });
        const sw = document.createElement("span");
        Object.assign(sw.style, {
            width: `${swatch}px`,
            height: `${swatch}px`,
            borderRadius: "2px",
            flex: "0 0 auto",
            background: scale.color(v),
        });
        item.append(sw, document.createTextNode(scale.name(v)));
        div.appendChild(item);
    }
    return div;
}

let fxCounter = 0;

// Wire Auspice-style hover highlighting: hovering a legend swatch fades every
// plot mark and legend item except the hovered variant's to DIM opacity.
//
// Driven entirely by ONE per-instance CSS rule (no per-element work on hover), so
// it stays O(1) even for the dot-dense multi-panel time-vs-frequency grid — an
// earlier per-mark version stalled the browser and starved Plot's tooltip pointer
// there. Lines/areas are matched per series: a mark is "the hovered variant's"
// when its stroke/fill equals that variant's color (or, for the time-vs-fitness
// centerline, its darker color in `darker`); major-variant colors are unique, so
// this is order-independent. The empirical dots, however, number in the tens of
// thousands — per-circle opacity (let alone a transition over them) is hopeless —
// so the whole dot LAYER is dimmed at the group level (a handful of <g>s) and
// without a transition, receding the cloud in one composite so the hovered
// variant's modeled line stands out.
//
// `root` is the component's stable root element (it contains both the legend and
// the plot SVG[s], even across redraws that rebuild its children). Call once per
// render; returns a teardown that removes the injected <style>.
export function linkLegendHighlight(root, scale, { darker } = {}) {
    const id = ++fxCounter;
    root.dataset.fxScope = id;
    const scope = `[data-fx-scope="${id}"]`;

    // Persistent: the per-series marks (few) + legend items transition smoothly.
    // The dot LAYER is intentionally excluded — it dims instantly (see dimRule).
    const transition =
        `${scope} g[aria-label="area"] path,` +
        `${scope} g[aria-label="line"] path,` +
        `${scope} [data-variant]{transition:opacity ${DURATION}ms ease-in-out}`;

    const style = document.createElement("style");
    style.textContent = transition;
    document.head.appendChild(style);

    // Dim everything whose color/variant is NOT the hovered one; matched marks
    // keep their default opacity. The hovered variant's line keeps stroke=color
    // (modeled/trajectory line) and stroke=darker (time-vs-fitness centerline).
    // The dot layer (`g[aria-label="dot"]`) is dimmed wholesale at the group level
    // — cheap regardless of how many thousands of circles it holds.
    function dimRule(v) {
        const c = scale.color(v);
        const dc = darker?.get(v);
        const lineNot = dc ? `:not([stroke="${c}"]):not([stroke="${dc}"])` : `:not([stroke="${c}"])`;
        return (
            `${scope} g[aria-label="area"] path:not([fill="${c}"]),` +
            `${scope} g[aria-label="line"] path${lineNot},` +
            `${scope} g[aria-label="dot"],` +
            `${scope} [data-variant]:not([data-variant="${v}"]){opacity:${DIM}}`
        );
    }

    let current = null;
    const highlight = (v) => {
        if (v === current) return;
        current = v;
        style.textContent = transition + dimRule(v);
    };
    const reset = () => {
        if (current === null) return;
        current = null;
        style.textContent = transition;
    };

    // Delegated on the stable root: legend items carry data-variant, marks do not,
    // so hovering marks is ignored (and never interferes with Plot's tooltip).
    const onOver = (e) => {
        const item = e.target.closest?.("[data-variant]");
        if (item && root.contains(item)) highlight(item.dataset.variant);
    };
    const onOut = (e) => {
        const item = e.target.closest?.("[data-variant]");
        if (item && !item.contains(e.relatedTarget)) reset();
    };
    root.addEventListener("mouseover", onOver);
    root.addEventListener("mouseout", onOut);

    return () => {
        root.removeEventListener("mouseover", onOver);
        root.removeEventListener("mouseout", onOut);
        style.remove();
    };
}
