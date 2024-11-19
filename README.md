# Revun

Reconstruction script for VTK/VTU data format based on numerical mode decomposition

## Version

v1.0.1

# Code description

`Revun` drives dynamic mode decomposition based on VTK/VTU files.


# How to start mode decomposition

## Execution

```console
python3 src/mode_decomposition.py
```

Tutorial case: `testcase/work*`

## Configuration file

Mode decomposition by `Revun` is controled by the configuration file: `config.yml`.

## Requirements

`Revun` requires the following packages:

- numpy (>=1.22.3)
- yaml (>= 5.3.1)
- VTK  (>= 9.2.6)


# Contact:

Yusuke Takahashi, Hokkaido University

ytakahashi@eng.hokudai.ac.jp


# References

- Yusuke Takahashi, "Spatiotemporal mode extraction for fluid–structure interaction using mode decomposition", Journal of Sound and Vibration, Volume 597, Part A, 2025, pp.118804. https://doi.org/10.1016/j.jsv.2024.118804.
