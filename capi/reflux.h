/*
 * reflux.h: reflux's C API (M10a), ABI version 1.
 *
 * The library is built from the V module `capi` by `v run tools/python.vsh
 * build`, with gcc, and exports exactly the functions declared here.
 *
 * Conventions (docs/m10a-plan.md, D2):
 * - Every function returns an int32_t: REFLUX_OK, or a negative code for a
 *   failure of the whole call. Where a function takes a buffer (err, err_len),
 *   it writes the reason there, NUL-terminated and cut to fit at a character
 *   boundary. A NULL buffer, or one of no length, gets nothing.
 * - Integers are int32_t or int64_t, never int. An array is a pointer and a
 *   length, and the caller allocates every array a function fills.
 * - Text going in is NUL-terminated UTF-8, and a list of texts an array of
 *   pointers and a count. The library never writes through them: they lack
 *   `const` only because V, which defines these functions, can't say it.
 * - Text coming out goes into a buffer the caller passes (buf, len),
 *   NUL-terminated, and the room it takes, its length and the NUL, goes to
 *   *needed unless needed is NULL. Text that doesn't fit is REFLUX_TOO_SMALL,
 *   with nothing written to buf: call again with *needed bytes.
 * - A handle is an int64_t id that is never used twice, so a stale or freed
 *   one is REFLUX_BAD_HANDLE rather than a crash. A database and a fluid are
 *   independent: freeing one leaves the other as it was.
 * - Any thread may call, and several at once. But the library starts itself
 *   on its first call, and that start must not race: make one call, such as
 *   reflux_abi_version(), from one thread before calling from several.
 */
#ifndef REFLUX_H
#define REFLUX_H

#include <stdint.h>

#define REFLUX_ABI_VERSION 1

/* A call's result. */
#define REFLUX_OK 0
#define REFLUX_BAD_ARGUMENT (-1) /* an argument it can't use */
#define REFLUX_BAD_HANDLE (-2)   /* an id that isn't a live handle of its kind */
#define REFLUX_NOT_FOUND (-3)    /* no compound by that name, or more than one */
#define REFLUX_DATA_ERROR (-4)   /* data it can't read, or can't build from */
#define REFLUX_TOO_SMALL (-5)    /* text that doesn't fit: see *needed */

/* Equations of state. */
#define REFLUX_VDW 0
#define REFLUX_RK 1
#define REFLUX_SRK 2
#define REFLUX_PR 3

/* α functions. */
#define REFLUX_ALPHA_EQUATION 0        /* the equation's own: Soave's, for SRK and PR */
#define REFLUX_ALPHA_MATHIAS_COPEMAN 1 /* ChemSep's Mathias–Copeman constants, for PR */

/* reflux_abi_version returns REFLUX_ABI_VERSION, for a caller to check before
   it calls anything else. */
int32_t reflux_abi_version(void);

/* reflux_version writes reflux's version, such as "0.1.0", into buf. A buffer
   too small for it is REFLUX_BAD_ARGUMENT. */
int32_t reflux_version(char *buf, int32_t len);

/* reflux_free drops a handle. */
int32_t reflux_free(int64_t id, char *err, int32_t err_len);

/* reflux_collect runs a full garbage collection, and writes the collector's
   heap size in bytes to heap_bytes unless it's NULL: for tests of memory. */
int32_t reflux_collect(int64_t *heap_bytes);

/* ---- Component data (docs/data-format.md) ---- */

/* reflux_data_root writes the data root a database opens when given none:
   $REFLUX_DATA if that's set, or else the data directory of the source tree
   the library was built from. A data root holds components/, a file per
   compound, and vendor/chemsep/, ChemSep's kij banks and Mathias–Copeman
   constants. */
int32_t reflux_data_root(char *buf, int32_t len, int32_t *needed);

/* reflux_database_open opens a data root, or the default one if root is NULL
   or "", and writes the database's handle to *db. The n_extra directories of
   component files in extra are read after the root's compounds, and replace
   those with the same CAS number. */
int32_t reflux_database_open(char *root, char **extra, int32_t n_extra, int64_t *db,
                             char *err, int32_t err_len);

/* reflux_database_index writes every compound the database holds, in the
   format of index.toml: each one's CAS number, name, InChIKey and aliases. */
int32_t reflux_database_index(int64_t db, char *buf, int32_t len, int32_t *needed,
                              char *err, int32_t err_len);

/* reflux_compound finds a compound by CAS number, name, alias or InChIKey,
   ignoring case, and writes its component file, TOML. A query that finds
   none, or more than one, is REFLUX_NOT_FOUND, whose reason names the
   closest. */
int32_t reflux_compound(int64_t db, char *query, char *buf, int32_t len, int32_t *needed,
                        char *err, int32_t err_len);

/* ---- Fluids ---- */

/* reflux_fluid_new builds a fluid of the nc compounds named, for an equation
   of state (REFLUX_VDW, REFLUX_RK, REFLUX_SRK or REFLUX_PR), and writes its
   handle to *fluid.
   - alpha is REFLUX_ALPHA_EQUATION or REFLUX_ALPHA_MATHIAS_COPEMAN (PR only);
     a compound mc-PR.lib lacks keeps the equation's own.
   - kij is 1 for ChemSep's kij for the equation (PR and SRK have a bank),
     taking a pair's first listing, or 0 for none.
   - n_overrides pairs set by hand follow: the compounds pair_a[k] and
     pair_b[k], both in the fluid and named any way the database knows, get
     the kij pair_k[k]. */
int32_t reflux_fluid_new(int64_t db, int32_t equation, char **names, int32_t nc, int32_t alpha,
                         int32_t kij, char **pair_a, char **pair_b, double *pair_k,
                         int32_t n_overrides, int64_t *fluid, char *err, int32_t err_len);

/* reflux_fluid_setup writes how a fluid was built, as TOML: the equation, the
   kij bank's file, then each compound's CAS number, name, Tc (K), Pc (Pa),
   acentric factor, molar mass (kg/mol) and alpha function, then each pair
   i < j's kij and where it came from: "file" (with its line and comment, and
   the lines of any listings not used), "override" or "missing". Floats are
   written so that they read back exactly. */
int32_t reflux_fluid_setup(int64_t fluid, char *buf, int32_t len, int32_t *needed, char *err,
                           int32_t err_len);

#endif
