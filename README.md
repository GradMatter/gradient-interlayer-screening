# Gradient Interlayer Screening

Code for screening compatible intermediate layers between two end materials
using Miedema thermodynamic descriptors and a coefficient-of-thermal-expansion
filter.

Requires Python 3.8+; no third-party packages.

## Run

```bash
python screening.py
```

The script reads `parameters.json` and writes `results.csv`.

## Files

- `screening.py`: self-contained screening code and element parameters.
- `parameters.json`: end materials, thresholds, candidate filters, and CTE values.
- `results.csv`: integrated screening results.

## Results

`results.csv` contains four `section` values:

- `end_a` and `end_b`: retained candidate elements.
- `cross_paths`: compatible intermediate element pairs.
- `cte_ranked`: monotonically increasing CTE paths ranked by `elin`.

`worst_dh` is in kJ mol^-1, `worst_delta` is in %, and `elin` uses the same
CTE units as `parameters.json`.

## Code Availability

The screening code used in this study is available at GitHub
(https://github.com/GradMatter/gradient-interlayer-screening).

## License

MIT.
