# Conda packaging

`meta.yaml` builds the current checkout locally. Run `./build_and_upload.sh`
with conda-build already installed; the script builds and tests without uploading.
Python packaging metadata remains the source of the runtime dependencies.

For Bioconda, generate a recipe from the final source distribution:

```bash
python scripts/make_bioconda_recipe.py dist/circleseeker-1.5.2.tar.gz /tmp/circleseeker-meta.yaml
```

The generated recipe pins the release sdist URL and its actual SHA256. Submit it
to `bioconda/bioconda-recipes` only when the matching GitHub release asset is
publicly available. Update the existing circleseeker recipe; retain old packages.
The obsolete 0.9.8 alternate recipe has been removed from this checkout and remains
available in Git history.
