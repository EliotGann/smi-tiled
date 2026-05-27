(installation)=
# Installation

`smi-tiled` works on Linux, macOS, and Windows under Python 3.10–3.13.
It has **no dependency on PyHyperScattering** at runtime; the two
packages can be installed alongside each other, but neither requires
the other.

## With pixi (recommended for development)

[Pixi](https://pixi.sh) is the most reliable way to set up a working
environment because it manages both the conda and PyPI sides of the
dependency tree from a single lockfile.

```bash
git clone https://github.com/EliotGann/smi-tiled.git
cd smi-tiled
pixi install         # solves and installs the default environment
pixi shell           # drops you into it
```

Useful tasks (see `pyproject.toml` `[tool.pixi.tasks]`):

```bash
pixi run test         # run the test suite
pixi run test-quick   # short-form summary
pixi run lint         # flake8 over src/ and tests/
pixi run docs-build   # build these docs locally
pixi run docs-serve   # serve at http://localhost:8765
```

## With pip

```bash
pip install smi-tiled[tiled]
```

The `[tiled]` extra pulls in `tiled[client]` + `bluesky-tiled-plugins`,
which is what the loader actually talks to.  Other optional extras:

| Extra      | Pulls in | What it enables |
|------------|---|---|
| `tiled`    | `tiled[client]`, `bluesky-tiled-plugins` | Loading from a Tiled server (always needed) |
| `upload`   | `httpx`, `zarr` | Writing reduced data back to a sandbox Tiled server |
| `dev`      | `pytest`, `coverage`, `flake8` | Running the test suite |
| `docs`     | `sphinx`, `pydata-sphinx-theme`, `myst-parser`, … | Building this documentation |
| `all`      | everything | Convenience meta-extra |

Example: full development install with documentation:

```bash
pip install -e ".[tiled,upload,dev,docs]"
```

## Verifying the install

```python
import smi_tiled
print(smi_tiled.__version__)
print(smi_tiled.TiledSMISWAXSLoader)
print(smi_tiled.defaults.default_saxs_mask_path())
```

If the bundled mask path resolves and the loader imports cleanly,
you're ready to {ref}`load a scan <quickstart>`.

## Authenticating with the Tiled server

`smi-tiled` defaults to the public NSLS-II Tiled catalog at
`https://tiled.nsls2.bnl.gov` under `smi/migration`.  Public read access
is usually available without credentials, but if you need to read a
private catalog or write to a sandbox, set the `TILED_API_KEY` environment
variable or pass `api_key=...` to {class}`smi_tiled.TiledSMISWAXSLoader`.

```bash
export TILED_API_KEY="…"
```

```python
loader = smi_tiled.TiledSMISWAXSLoader(api_key="…")
```

Interactive logins are also supported:

```python
loader = smi_tiled.TiledSMISWAXSLoader()
loader.login()              # opens the OAuth flow in your terminal
```

## Optional: PyHyperScattering interop

If you also use [PyHyperScattering](https://github.com/EliotGann/PyHyperScattering)
for downstream analysis (`da.rsoxs.*` slicing, `da.fit.*` curve fitting,
etc.), install it separately:

```bash
pip install pyhyperscattering
```

`smi-tiled` emits DataArrays carrying the PyHyper geometry attrs
(`dist, poni1, poni2, pixel1, pixel2, energy, wavelength`), so the SAXS
output is consumable by `PFGeneralIntegrator(geomethod='template_xr')`.
WAXS is **not** — its 3-panel arc geometry requires our own integrator
(see {doc}`user-guide/reduction`).
