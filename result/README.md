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
