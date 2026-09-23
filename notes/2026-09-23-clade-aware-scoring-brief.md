# Similarity-aware scoring for clade frequency forecasts

## The problem

The normalized energy score and the cross-entropy loss in the Methods both treat clade
labels as exchangeable. Putting forecast mass on clade B when B.1 sweeps is scored
exactly like putting it on A. For vaccine composition that is the wrong accounting: B and
B.1 share nearly all their HA1 amino acids, so a B-matched vaccine is a near hit, while an
A-matched vaccine is a miss. The metric should track that.

This bites in a specific regime. The <30-sequence collapse already folds small descendants
into their parents, and future observed sequences are remapped onto the clade set defined
at the forecast timepoint — so many near-misses are silently absorbed before scoring. What
survives is exactly the case where a parent and its descendant both clear the threshold and
co-circulate, which is also the hardest call for strain selection.

*(The hierarchy-smoothed cross-entropy from the first draft is dropped. It wasn't a proper
scoring rule, and it added a second tuning knob for little gain.)*

## The distance between two clades

One definition, used by both options below.

**Represent each clade by the reconstructed HA1 amino acid sequence at its MRCA node, and
define the distance between two clades as the number of HA1 positions at which those two
sequences differ.** No curated site list, no weighting, no antigenic model.

`augur ancestral` already writes these node sequences, so this costs nothing new. It is
defensible precisely because there is nothing to argue about: a reviewer can recompute it
from the alignment, and the only upstream assumption is the alignment itself, which is
already in the paper.

Use the MRCA node sequence rather than the consensus of the sequences assigned to the
clade. Collapsed clades have swallowed their small descendants, so their consensus is a
mixture that shifts whenever the collapse threshold changes; the MRCA sequence doesn't.

## Option A — score the mutation profile instead of the clade label

### The idea, in plain terms

Right now the scorer sees a forecast as a list of numbers, one per clade, and compares it
to the observed list. It has no way of knowing that B.1 sits right next to B on the tree.

The fix is to translate both lists into something the scorer can compare meaningfully:
**for each HA1 mutation, what fraction of the viral population carries it?** A forecast
that confuses B with B.1 produces almost the same answer to that question, because those
two clades share almost all their mutations. A forecast that confuses B.1 with a distant
clade produces a very different answer. So near-misses are cheap and real misses are
expensive, automatically, without inventing a new scoring rule.

It also happens to be the question a vaccine strain choice actually turns on.

### Step 1 — build the clade x mutation table

Take the clades present at this timepoint. List every HA1 position where any two of them
differ; each becomes a column. (If a position has more than two amino acids across clades,
give it one column per amino acid.) Each clade is a row, with 1 where it carries that
amino acid and 0 where it doesn't.

A toy example with four clades and five variable positions:

| clade | 145K | 159Y | 186S | 193F | 276E |
|-------|------|------|------|------|------|
| A     | 0    | 0    | 0    | 0    | 0    |
| B     | 1    | 1    | 0    | 0    | 0    |
| B.1   | 1    | 1    | 1    | 0    | 0    |
| C     | 0    | 0    | 0    | 1    | 1    |

Real timepoints will have roughly 10-40 columns rather than five.

The distance between two rows is the square root of the number of columns where they
differ. So if B.1 is the clade that actually wins:

| forecast put its mass on | columns differing from B.1 | distance |
|---|---|---|
| B (the parent)  | 1 | 1.00 |
| A (the root)    | 3 | 1.73 |
| C (a distant clade) | 5 | 2.24 |

Today all three of those cost the same. That table is the whole point of the exercise, and
it is small enough to check by hand for any given timepoint.

The square root isn't a choice — it falls out of using the standard energy score in this
space. It has the mild side effect of compressing large differences, so the fifteenth
amino acid difference counts for less than the second, which is probably what you want
anyway.

### Step 2 — add one identity column per clade

Append four more columns, one per clade, with a small value `eps` on the diagonal and 0
elsewhere. This guarantees that two clades never collapse onto the same row even if their
HA1 sequences happen to be identical.

`eps` is the dial between the old metric and the new one:

- **`eps` large** — the identity columns swamp the mutation columns, every clade sits the
  same distance from every other, and you recover today's MNES exactly.
- **`eps = 0`** — pure mutation profile, maximum credit for near-misses.
- **In between** — adding `eps` raises every squared distance by `2 x eps^2`, so setting
  `eps ≈ 0.7` puts a floor of about one amino acid difference between any two clades.

I'd report `eps = 0` and one moderate value, and show the curve in supplement.

### Step 3 — normalize per timepoint

Divide the whole table (mutation columns and identity columns together) by the largest
clade-to-clade distance at that timepoint, so the biggest possible distance is 1
everywhere. Without this, B lineages score better than A lineages for a reason that has
nothing to do with forecast skill: shallower trees, fewer variable positions, smaller
distances, smaller error. That would manufacture a fake result given the adaptive-rate
comparison in Section 3.

### Step 4 — translate every forecast into a mutation profile

For a single posterior draw giving clade frequencies of, say, A=0.1, B=0.3, B.1=0.5,
C=0.1, multiply each clade's row by its frequency and add them up. The result has one
number per mutation: for 145K you'd get 0.3 + 0.5 = 0.8, i.e. this draw predicts that 80%
of circulating viruses will carry 145K.

Do this for all 100 posterior draws and for the observed frequencies.

### Step 5 — run the existing energy score, unchanged

Feed those profile vectors to the same scorer instead of the clade-frequency vectors.
Same posterior draws, same observed counts, same normalization by N, same code path. The
only new line is one matrix multiply in front of it. Nothing is refit.

### Why this is still a legitimate score

The energy score is *proper*, meaning you can't improve your score by reporting a
distribution you don't actually believe. That property survives here because multiplying
by a fixed table is just a change of coordinates — like scoring in Fahrenheit instead of
Celsius — and a change of coordinates can't create a way to game the score. The one thing
that would break it is if two different clade mixtures mapped onto the identical profile,
and the identity columns from Step 2 are exactly what rules that out.

### Sanity checks before trusting any number

1. Set `eps` to something huge and confirm the new score reproduces the current MNES up to
   a constant scale factor. That's the regression test.
2. Print the distance table for one timepoint and eyeball it against the tree. Parent and
   child should be the closest pair; anything else means the node sequences are misaligned.
3. Take a real timepoint, artificially shift forecast mass from the winning clade onto its
   parent, then onto a distant clade, and confirm the scores come out in that order.
4. Count the mutation columns per timepoint. If that number swings wildly across timepoints,
   normalization is doing a lot of work and deserves a supplementary figure.

## Option C — expected HA1 distance

The number for the vaccine composition audience:

```
L(t, h) = sum over clades k of  p_k(t + h) * d(k, k*)
```

where `p` is the forecast frequency at horizon `h`, `k*` is the clade observed at highest
frequency at `t + h`, and `d` is the raw count of HA1 amino acid differences.

Read it as: *if you picked a strain by sampling from the forecast, how many HA1 amino acid
differences would separate it from the clade that actually dominated?* Report at h = 270
and 365 days.

Note that this uses raw counts, not the square root from Option A. That's deliberate — a
plain expectation reads better in whole amino acids, and nothing here depends on the
geometry. Worth a footnote so nobody reads it as an inconsistency.

If the single-winner framing feels too crude, the fuller version averages `d` over both the
forecast and the observed population rather than collapsing the truth to one clade. Same
distance matrix either way.

## Things to watch

1. **The naive baseline should improve.** Persistence puts mass on the currently dominant
   clade, which is usually the parent of whatever sweeps. A phylogeny-aware metric is kind
   to persistence, so the "naive global is consistently best" result may shift. If the gap
   at 12 months narrows under Option A, that is a finding in itself, and it sharpens the
   Discussion point about overconfident GAs.
2. **Tables 1 and 2 may reorder.** Check whether the geo-resolution and window-size
   recommendations are stable under the new score before committing to them.
3. **`eps` is a free parameter and a reviewer will say so.** Show the sensitivity curve,
   pre-register one value, and anchor it to something interpretable — the one-amino-acid
   floor above is the natural choice.
4. **Sparse clades give noisy MRCA reconstructions**, which puts noise in the distance
   matrix rather than in the forecast. The threshold of 50 helps; worth confirming the
   distances are stable at the country resolution before reporting country-level numbers.

## What to run first

H3N2 only, regional resolution, model m11, HA1 distances, Option A at three values of
`eps`. Recompute Figure 3A. If the MLR and naive curves separate differently than they do
now, that justifies rolling it out to the other lineages and resolutions.
