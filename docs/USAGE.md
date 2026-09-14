# Using the release assets

Keep the primary and cosmetic packages in separate directories. Both include relative Workbench specifications and their own label dictionaries. For the analysis package, `0` can mean background, excluded support or unassigned cortex; use its cortex and uncertainty masks rather than treating every zero as noncortical tissue.

Verify the downloaded ZIPs against `SHA256SUMS.txt`. On Windows, use `Get-FileHash -Algorithm SHA256 <archive.zip>`. After extraction, run `python validate_release.py` from the primary package, or `python code/finalize_display.py verify` from the display-only package. Source-preservation checks for the display companion require the original source packages when available; geometric/label checks use its included files.

Midthickness uses the corresponding age-template physical space. Inflated surfaces are display geometry. The cosmetic dense grid is interpolated and is not a new anatomical standard or a new MRI resolution.

The full manuals are included inside the ZIPs. The repository provides a compact usage guide and method/limitation summary; package-relative data links are meant to be followed from the extracted packages.
