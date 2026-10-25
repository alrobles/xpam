# cran-comments.md — xpam submission notes

## First submission (1.0.0)

This is the first CRAN submission of `xpam`.

The package was developed under the name `pam`, but `pam` is already on
CRAN (verified via `available.packages()`), so the package was renamed
`xpam` throughout: DESCRIPTION, `useDynLib`, Rcpp exports, tests,
vignette, man pages and CITATION.  The `pam` 0.1.0--0.2.0 code line was
never submitted to CRAN; the GitHub-only history is documented in
`NEWS.md`.

The package supersedes the GitHub-only `alrobles/PAM` repository, which
is retired in favour of the public `alrobles/xpam` repository.

## Check status

Local check: R 4.6.1, Ubuntu 24.04, GDAL 3.11.4, on the built
`xpam_1.0.0.tar.gz` with `--as-cran`: **0 ERROR, 0 WARNING, 2 NOTEs**.

1. `CRAN incoming feasibility ... NOTE` — "New submission".  An
   earlier check of this same tree reported the package URLs as 404
   because the public repository did not exist yet; the repository and
   its pkgdown site have since been created from this tree and the URLs
   now resolve (verified live before submission).
2. `HTML version of manual ... NOTE` — "Skipping checking HTML
   validation: no command 'tidy' found" and "Skipping checking math
   rendering: package 'V8' unavailable".  A property of the local check
   environment, not of the package.

The package was additionally checked without GDAL
(`--configure-vars='--with-gdal-config=/nonexistent'`) on the same
`xpam_1.0.0.tar.gz`: 0 ERROR, 0 WARNING, 2 NOTEs (the same two),
confirming that machines without GDAL install, load and check cleanly.

## Tests

The test suite was run twice.  Under `R CMD check --as-cran` (network
restricted) the network-dependent tests are skipped by design:
`FAIL 0 | WARN 0 | SKIP 38 | PASS 31`.  With `NOT_CRAN=true` and the
network available, the whole suite runs:
`FAIL 0 | WARN 0 | SKIP 0 | PASS 247`.  The 38 skipped tests are gated
with `skip_on_cran()` + `skip_if_offline()`, not disabled.

## Documentation and vignettes

Following the convention used by our other CRAN packages (`xsdm`,
`xnicher`), the vignette is **not** shipped in the CRAN tarball: it is
built by pkgdown and published on the package website
(<https://alrobles.github.io/xpam/>).  `vignettes/` and `inst/doc/` are
excluded via `.Rbuildignore`, and the `VignetteBuilder` field is
therefore absent so the source package carries no prebuilt vignette
index.  `README.md` and `NEWS.md` are shipped in the tarball.

## Examples and network access

`R CMD check --run-donttest` passes with the network unreachable.
`get_rodentia_mexico()` and `get_world_pca()` document a
`"Can't access data."` character return when their data server is
unreachable, and their examples test for it with `is.character()`; the
`get_geo_pam()` example runs a fully local `SpatVector` example that
needs no download and only then exercises the guarded remote path.
Network-dependent tests are gated with `skip_on_cran()` +
`skip_if_offline()`.

## GDAL is optional

`configure` probes for GDAL; the native tiled raster path links only
when GDAL is found.  Without GDAL the package builds, loads and the
terra-based fallback works; engine tests skip via `skip_without_gdal()`.
`configure` probes are hermetic: `LIBS` is cleared at the top and the
GDAL link test uses only `gdal-config --libs`, so an ambient `LIBS` in
the build environment cannot silently disable the engine.  Verified on
`xpam_1.0.0.tar.gz` under a poisoned `LIBS` (`-lbogus_poison`, GDAL
still detected and `HAVE_GDAL` defined) and a clean one.

## Resolved during pre-submission

* `.hpp` -> `.h` in `src/` (unlikely file names WARNING).
* `test_check("pam")` -> `test_check("xpam")` (rename residual, test
  ERROR).
* GitHub URLs pointed at the private devel repository -> now
  `alrobles/xpam` (public production repository).
* `.aspell/` dictionary committed (`defaults.R` + `xpam.rds`) and
  `inst/WORDLIST` added (89 words); `spelling::spell_check_package()`
  returns 0 flagged words across DESCRIPTION, Rd files, README, NEWS
  and the vignette.
* OpenMP regions honour `OMP_THREAD_LIMIT` via safe thread counts.
* Vignettes excluded from the tarball; `VignetteBuilder` dropped.
* `README.md` and `NEWS.md` now shipped in the tarball.
