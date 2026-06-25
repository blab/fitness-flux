# Mutation counts

Per-Pango-lineage amino-acid mutation counts for SARS-CoV-2, consumed by the
lineage-deltas analysis (`sarscov2_lineages_mut_counts.tsv`: columns
`nuc/spike/s1/s2/ntd/rbd/orf1ab/accessory_muts`).

These are produced by the Snakemake `mutation_counts` rules (`rules/mutation_counts.smk`):

1. **Fetch** the Nextclade SARS-CoV-2 reference tree:
   ```
   nextstrain remote download https://nextstrain.org/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs data/nextclade_sarscov2_pango
   ```
2. **Count** with `count_lineage_mutations.py`, which walks each tip's root→leaf path
   and counts the **net** amino-acid substitutions per gene (positions whose final
   state differs from Wuhan-Hu-1), splitting spike by sub-region (S1 14–685, S2
   686–1273, NTD 14–305, RBD 319–541) and taking the median over the tips of each
   Pango lineage. Nucleotide count is the tree's `div`. Output:
   `results/sarscov2_lineages_mut_counts.tsv`.

Counting **net state** rather than summing mutation *events* along the path avoids
the ~5× inflation that Nextclade's recombinant graft nodes (`rec_parent` /
`internal_X`) otherwise produce for the entire XBB/recombinant radiation (e.g.
XBB.1 spike 226 → 42).

Refresh as new Pango lineages are defined:
```
snakemake --forcerun fetch_nextclade_pango mutation-counts/results/sarscov2_lineages_mut_counts.tsv
```

The Mathematica notebook `mutation-counts-sarscov2-lineages.nb` (which summed events
along the path) is retained for reference but is superseded by the Python script.
