AC_INIT([xpam], [1.0.0])
AC_CONFIG_AUX_DIR([tools])

# Hermetic probes.  R CMD INSTALL exports the user's ~/.R/Makevars into
# the environment, so ambient LIBS (e.g. -lRcppParallel appended there)
# would poison both the compiler sanity check and the GDAL link test on
# machines whose Makevars adds libraries configure knows nothing about.
# No probe in this file needs ambient libraries — the package links via
# PKG_LIBS in src/Makevars, never through configure's LIBS.
LIBS=""

# The C++ compiler must be resolved before the GDAL probe: AC_LANG_PUSH
# below would otherwise trigger AC_PROG_CXX lazily *inside* the link test.
AC_PROG_CXX

# --- GDAL detection (optional, 5-layer anti-friction pattern from xbioclim) --
# The engine compiles unconditionally; HAVE_GDAL only gates whether the native
# tiled raster path actually links. Without GDAL the package still builds,
# loads, and the R fallback (pam_env + terra) works.
AC_ARG_WITH([gdal-config],
  AS_HELP_STRING([--with-gdal-config=PATH],
    [full path to the gdal-config script]),
  [GDAL_CONFIG="${withval}"],
  [GDAL_CONFIG=""]
)

GDAL_CPPFLAGS=""
GDAL_LIBS=""
HAVE_GDAL=no

if test -z "${GDAL_CONFIG}"; then
  AC_PATH_PROG([GDAL_CONFIG], [gdal-config], [not-found])
else
  AC_MSG_CHECKING([for gdal-config at ${GDAL_CONFIG}])
  if test -x "${GDAL_CONFIG}"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    GDAL_CONFIG="not-found"
  fi
fi

if test "${GDAL_CONFIG}" != "not-found"; then

  GDAL_CPPFLAGS=`${GDAL_CONFIG} --cflags 2>/dev/null`
  GDAL_LIBS=`${GDAL_CONFIG} --libs 2>/dev/null`

  GDAL_VERSION=`${GDAL_CONFIG} --version 2>/dev/null`
  AC_MSG_CHECKING([GDAL version >= 2.0.1])
  GDAL_VERSION_MAJOR=`echo "${GDAL_VERSION}" | cut -d. -f1`
  GDAL_VERSION_MINOR=`echo "${GDAL_VERSION}" | cut -d. -f2`
  GDAL_VERSION_PATCH=`echo "${GDAL_VERSION}" | cut -d. -f3`
  if test -z "${GDAL_VERSION_PATCH}"; then
    GDAL_VERSION_PATCH=0
  fi

  if test "${GDAL_VERSION_MAJOR}" -ge 3; then
    GDAL_OK=yes
  elif test "${GDAL_VERSION_MAJOR}" -eq 2 && test "${GDAL_VERSION_MINOR}" -ge 1; then
    GDAL_OK=yes
  elif test "${GDAL_VERSION_MAJOR}" -eq 2 && test "${GDAL_VERSION_MINOR}" -eq 0 && test "${GDAL_VERSION_PATCH}" -ge 1; then
    GDAL_OK=yes
  else
    GDAL_OK=no
  fi

  if test "${GDAL_OK}" = "yes"; then
    AC_MSG_RESULT([yes (${GDAL_VERSION})])

    # Real link test: "gdal-config exists" is not enough — the flags must
    # actually compile and link a GDAL call.  The probe is hermetic in
    # LIBS: ambient libraries (e.g. a user ~/.R/Makevars that exports
    # -lRcppParallel via PKG_LIBS into the environment R CMD INSTALL
    # gives configure) must not be able to fail this test — a missing
    # ambient library is not a missing GDAL.  gdal-config --libs already
    # carries GDAL's own dependencies.
    SAVED_CPPFLAGS="${CPPFLAGS}"
    CPPFLAGS="${GDAL_CPPFLAGS} ${CPPFLAGS}"
    LIBS="${GDAL_LIBS}"    # hermetic probe: ambient LIBS stays out

    AC_LANG_PUSH([C++])
    AC_MSG_CHECKING([whether GDAL links correctly])
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM(
        [#include <gdal.h>],
        [GDALAllRegister();]
      )],
      [AC_MSG_RESULT([yes])
       HAVE_GDAL=yes],
      [AC_MSG_RESULT([no - GDAL link test failed; building without GDAL])
       GDAL_CPPFLAGS=""
       GDAL_LIBS=""]
    )
    AC_LANG_POP([C++])

    CPPFLAGS="${SAVED_CPPFLAGS}"
    LIBS=""                                  # back to the cleared baseline

    if test "${HAVE_GDAL}" = "yes"; then
      GDAL_CPPFLAGS="${GDAL_CPPFLAGS} -DHAVE_GDAL"
    fi
  else
    AC_MSG_RESULT([no (${GDAL_VERSION}) - need >= 2.0.1; building without GDAL])
    GDAL_CPPFLAGS=""
    GDAL_LIBS=""
  fi

else
  AC_MSG_WARN([gdal-config not found; building without GDAL (native engine disabled)])
fi

AC_SUBST([GDAL_CPPFLAGS])
AC_SUBST([GDAL_LIBS])

# Object list enumerated explicitly — no GNU Make extensions (CRAN portable).
# xpam_gdal_io.o / pam_rasterise.o are always compiled: without -DHAVE_GDAL
# their bodies are informative-error stubs.
ALL_OBJECTS="PamEngine.o RcppExports.o pam_kernel.o pam_rasterise.o xpam_gdal_io.o"

AC_SUBST([ALL_OBJECTS])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
