// Variant color layer shared by the figure components.
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (buildColors, variantColor,
// majorLegend) on 2026-06-19.
//
// A "colors" array is a list of { variant, color, display_name, is_major, order }
// rows (the JSON form of {dataset}_colors.tsv).

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
        Object.assign(item.style, {
            display: "inline-flex",
            alignItems: "center",
            gap: "5px",
            whiteSpace: "nowrap",
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
