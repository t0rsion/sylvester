# 2026-08-19 v0.2 record

This directory freezes the pre-release comparison at commit `f413b0b`.
`REPORT.md` and `results.csv` are unchanged.

The raw `results.json` is 626,092,556 bytes. The repository stores its exact
zstd encoding to stay below host file limits:

```text
zstd -d results.json.zst -o results.json
sha256sum -c SHA256SUMS
```

`SHA256SUMS` includes both the stored archive and its decoded JSON.
