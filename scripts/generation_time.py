"""Per-variant generation time for SARS-CoV-2 (pre- vs post-Omicron).

Generation time tau converts a per-day MLR growth rate into per-generation fitness
(fitness = tau * per-day-rate). SARS-CoV-2 shortened its generation interval with
Omicron, so we assign tau by variant: pre-Omicron variants get `tau_pre` (5.0 days)
and Omicron variants `tau_post` (3.2 days). Because scaffolding stitches windows by
shared variants, tau must depend on variant identity (not the window), so this
classifier is the single source of truth used by both the MLR export (run-mlr-model.py)
and the fitness-flux velocity (fitness_wave.py).

`is_omicron_*` return True (Omicron), False (pre-Omicron), or None (unclassifiable,
e.g. the "other" catch-all or the MLR pivot), which callers map to the default tau.
"""

import re
from typing import Optional

# Nextstrain clade name, e.g. "20A", "21K", "22B".
_CLADE_RE = re.compile(r"^(\d+)([A-Z])")

# The only pre-Omicron Pango recombinants; every later X* recombinant is Omicron-era.
_PRE_OMICRON_RECOMBINANTS = {"XA", "XB", "XC"}


def is_omicron_clade(clade: str) -> Optional[bool]:
    if clade in ("other", "", "WT"):
        return None if clade != "WT" else False
    match = _CLADE_RE.match(clade)
    if match is None:
        return None
    year, letter = int(match.group(1)), match.group(2)
    if year >= 22:
        return True
    if year == 21:
        return letter >= "K"  # 21A-21J pre-Omicron; 21K/21L/21M Omicron
    return False  # 19*, 20*


def is_omicron_lineage(lineage: str, aliasor) -> Optional[bool]:
    if lineage in ("other", ""):
        return None
    try:
        full = aliasor.uncompress(lineage)
    except Exception:
        full = lineage
    if full.startswith("B.1.1.529"):
        return True
    if lineage.startswith("X"):
        return lineage.split(".")[0] not in _PRE_OMICRON_RECOMBINANTS
    return False


def load_aliasor(alias_path: Optional[str] = None):
    """Construct a Pango Aliasor. `alias_path=None` downloads alias_key.json from
    GitHub, matching scripts/collapse-lineage-counts.py."""
    from pango_aliasor.aliasor import Aliasor

    return Aliasor(alias_path)


def generation_time_for(
    variant: str,
    kind: str,
    tau_pre: float,
    tau_post: float,
    aliasor=None,
) -> float:
    """tau for a single variant. `kind` is "clades" or "lineages"."""
    if kind == "clades":
        omicron = is_omicron_clade(variant)
    elif kind == "lineages":
        omicron = is_omicron_lineage(variant, aliasor)
    else:
        omicron = None
    if omicron is False:
        return tau_pre
    return tau_post  # Omicron, "other"/pivot, and unclassifiable -> default
