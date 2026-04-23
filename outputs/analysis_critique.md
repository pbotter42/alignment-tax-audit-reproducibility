# Alignment-Tax Audit: Key Diagnostics

- Sample after filtering: **591 models x 12061 items**.
- PC1 variance explained: **36.26%**.
- Raw cosine coherence vs ability correlation: **0.994** (near-collinearity diagnostic).
- Alternative coherence (1-factor R2) vs ability correlation: **-0.333**.
- Jackknife impact vs ability correlation: **0.542**.
- Jackknife impact vs 1-factor R2 correlation: **-0.922**.
- Quadratic fit significantly improves over linear for impact~ability: p = **2.55e-292**.

## Interpretation Guardrails

- The original raw-cosine coherence metric is highly coupled to PC1 ability and should be treated as a diagnostic, not a standalone construct-validity endpoint.
- The stronger supported pattern is a **middle-ability fracture** (inverted-U impact profile), not a monotonic 'smarter = more fractured' trend.
- High-ability models in this snapshot are often near-neutral or negative impact; many strongest polluters cluster in the middle ability range.
- Model-type heuristics (merge/instruction keywords) are noisy labels and should be interpreted as exploratory.
