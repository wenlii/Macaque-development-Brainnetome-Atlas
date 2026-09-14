# Code and input requirements

`validate_analysis_release.py` can validate an extracted analysis ZIP: run it with the extracted package directory as its argument. Its `--hash-only` mode uses only the Python standard library; full checks need numpy and nibabel.

`analysis_pipeline/` retains the exact v4 repair/projection/audit source files. They require the full native reconstruction/registration layout in the original research archive, plus Workbench, ANTsPy and Python dependencies. They are provided for method inspection and are not a Run-All analysis using the minimal ZIP alone.

`display_only/` retains the exact cosmetic processing/rendering source. Use the identical scripts inside the extracted display-only ZIP to render its included surfaces; its `code/` location is part of that package's relative-path contract. Rebuilding cosmetic labels additionally requires the v4 source packages and full source-caution inputs documented in that ZIP.

No reconstruction binaries or credentials are included. These source files are preserved rather than silently rewritten for a different directory layout.
