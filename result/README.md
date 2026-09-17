# Local result archive

`result/` is the local home for generated numerical outputs. It follows the
same four-way layout as the maintained workflows:

```text
result/
├─ gas/0d/
├─ gas/0d-cstr/
├─ surface/0d/
└─ surface/0d-cstr/
```

Create one dated campaign directory under the matching calculation type, for
example `result/surface/0d-cstr/campaign_YYYYMMDD/`. Keep the run summary,
case outputs, figures, logs, and the calculation handoff document together in
that directory.

Numerical outputs are intentionally local and ignored by Git. The short
`README.md` files that describe this layout remain tracked.

## Current archived surface campaigns

- `surface/0d/campaign_hong_closed0d_100s_production_20260912_r02/`:
  completed production scan (972/972 cases).
- `surface/0d-cstr/campaign_hong_surface_cstr_tau10ms_20260915/`:
  completed campaign execution with 647 accepted cases; the remaining cases
  are explicitly retained with their failed, invalid-output, or
  needs-extension status.
