// Host-side data helpers. Used by the dashboard for the figures that still load
// TSVs inline (mutation–fitness, flux wave); the extracted components receive
// pre-parsed JSON and never fetch.
// Adapted from fitness-flux@da79fa9:viz/fitness-flux.html (parseTSV, fetch*) on 2026-06-19.

export function parseTSV(text) {
    const lines = text.trim().split(/\r?\n/);
    const headers = lines[0].split("\t");
    return lines.slice(1).map((line) => {
        const values = line.split("\t");
        const obj = {};
        headers.forEach((header, i) => {
            const value = values[i];
            obj[header] =
                value !== undefined && value !== "" && !isNaN(value) ? +value : value;
        });
        return obj;
    });
}

export async function fetchText(url) {
    const response = await fetch(`${url}?t=${Date.now()}`);
    if (!response.ok) throw new Error(`${response.status} ${url}`);
    return response.text();
}

export async function fetchTSV(url) {
    return parseTSV(await fetchText(url));
}

export async function fetchJSON(url) {
    return JSON.parse(await fetchText(url));
}
