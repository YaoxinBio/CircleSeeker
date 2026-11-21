# Conda packaging

Bioconda distributes CircleSeeker v1.6.1 as build `1.6.1-0`. The published recipe uses the GitHub tag archive and its SHA256 checksum. [Installation instructions](../README.md#install-from-bioconda) · [Bioconda package](https://bioconda.github.io/recipes/circleseeker/README.html)

To build and test the current checkout locally, install `conda-build` and run the following from the repository root. The script uses `conda-recipe/meta.yaml` with uploads disabled.

```bash
bash conda-recipe/build_and_upload.sh
```

For a release that distributes a Python source distribution (sdist), the helper below generates a recipe using its filename and checksum. With Python's `build` package installed, run these commands from the repository root.

```bash
sdist_dir=$(mktemp -d)
python -m build --sdist --outdir "$sdist_dir"
python scripts/make_bioconda_recipe.py "$sdist_dir"/circleseeker-*.tar.gz "$sdist_dir/bioconda-meta.yaml"
```

The generated URL points to an uploaded release asset named `circleseeker-<version>.tar.gz`. Check that the URL serves the exact file used to generate the checksum before submitting an sdist-based update to the [Bioconda recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/circleseeker/meta.yaml).
