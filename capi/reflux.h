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
 *   reflux_abi_version(), from one thread before calling from several. That
 *   thread may end afterwards.
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

/* A point's status, in a status array: 0 for an answer, or why there's none.
   A point without an answer has NaN for every value. */
#define REFLUX_POINT_OK 0
#define REFLUX_POINT_BAD_INPUT 1        /* a value it can't use: see the call */
#define REFLUX_POINT_NO_CONVERGENCE 2   /* an iteration ran out of steps */
#define REFLUX_POINT_NO_VOLUME_ROOT 3   /* no volume root at a state it needed */
#define REFLUX_POINT_MORE_PHASES 4      /* a flash's phase is itself unstable */
#define REFLUX_POINT_PHASES_SWAPPED 5   /* a vapour-fraction flash's phases swapped */
#define REFLUX_POINT_OUTSIDE_ENVELOPE 6 /* no such point on the phase envelope */
#define REFLUX_POINT_SUPERCRITICAL 7    /* at or above Tc: no vapour pressure */
#define REFLUX_POINT_NOT_FINITE 8       /* an answer that isn't finite: beyond the model */

/* A phase's label, by its phase identification parameter. */
#define REFLUX_NO_PHASE (-1)
#define REFLUX_LIQUID 0
#define REFLUX_VAPOR 1

/* Volume roots, for reflux_props. Where there's one root, it's that one. */
#define REFLUX_ROOT_LIQUID 0 /* the smallest volume */
#define REFLUX_ROOT_VAPOR 1  /* the largest */
#define REFLUX_ROOT_STABLE 2 /* the one with the lower Gibbs energy */

/* reflux_props' properties, per mole of the phase: out[k * n + i] is property
   k at point i. */
#define REFLUX_PROP_Z 0      /* compressibility factor, PV/(nRT) */
#define REFLUX_PROP_V 1      /* molar volume, m3/mol */
#define REFLUX_PROP_H 2      /* enthalpy, J/mol */
#define REFLUX_PROP_S 3      /* entropy, J/(mol K) */
#define REFLUX_PROP_G 4      /* Gibbs energy, J/mol */
#define REFLUX_PROP_U 5      /* internal energy, J/mol */
#define REFLUX_PROP_A 6      /* Helmholtz energy, J/mol */
#define REFLUX_PROP_CP 7     /* isobaric heat capacity, J/(mol K) */
#define REFLUX_PROP_CV 8     /* isochoric heat capacity, J/(mol K) */
#define REFLUX_PROP_W 9      /* speed of sound, m/s */
#define REFLUX_PROP_MU_JT 10 /* Joule-Thomson coefficient, K/Pa */
#define REFLUX_PROPS 11      /* how many there are */

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
   acentric factor, molar mass (kg/mol), dhf and s0 (its ideal gas's ΔHf in
   J/mol and S° in J/(mol K), at 298.15 K and 1 bar, where it has them) and
   alpha function, then each pair i < j's kij and where it came from: "file"
   (with its line and comment, and the lines of any listings not used),
   "override" or "missing". Floats are written so that they read back
   exactly. */
int32_t reflux_fluid_setup(int64_t fluid, char *buf, int32_t len, int32_t *needed, char *err,
                           int32_t err_len);

/* ---- Many points at once ----

   The functions below evaluate n points in one call, n >= 0. An input of
   per-point values has 1 entry, used for every point, or n. Each point gets
   a status (REFLUX_POINT_*), and a point that fails leaves the others as
   they would be alone. err gets the first failed point's reason, or "" if
   none failed, so the reason for one point is the same call made for it
   alone. A call fails as a whole only for a bad handle or an argument it
   can't use: a NULL array, a count other than 1 or n, an unknown root. */

/* reflux_props evaluates one phase of the fluid at n points (T in K, P in
   Pa, and the mole amounts z, nc per point, normalised by their sum), taking
   the volume root asked for (REFLUX_ROOT_*). It writes the REFLUX_PROPS
   properties to out (REFLUX_PROPS * n values: see REFLUX_PROP_*), ln phi to
   ln_phi (n * nc), the phase's label to phase (n), and each point's status
   to status (n). A point's T or P must be finite and above zero, and its z
   finite, non-negative and not all zero: anything else is
   REFLUX_POINT_BAD_INPUT. H, S, G, U and A are on the ideal gas's formation
   basis (ΔHf and S° at 298.15 K and 1 bar), which needs every compound's
   ΔHf and S°: see reflux_fluid_setup. */
int32_t reflux_props(int64_t fluid, int32_t root, int64_t n, double *t, int64_t t_len, double *p,
                     int64_t p_len, double *z, int64_t z_rows, double *out, double *ln_phi,
                     int32_t *phase, int32_t *status, char *err, int32_t err_len);

/* reflux_psat writes one compound's vapour pressure, in Pa, at n
   temperatures, in K: out[i] at t[i], with its status in status[i]. The
   compound is the fluid's component'th, counting from 0, and its vapour
   pressure is the fluid's equation and alpha for it alone, where its liquid
   and vapour fugacities are equal. A temperature at or above its Tc is
   REFLUX_POINT_SUPERCRITICAL, and one that isn't finite and above zero
   REFLUX_POINT_BAD_INPUT. */
int32_t reflux_psat(int64_t fluid, int32_t component, int64_t n, double *t, double *out,
                    int32_t *status, char *err, int32_t err_len);

/* ---- Stability and flashes ---- */

/* reflux_stability tests whether each of n feeds (T in K, P in Pa, and z,
   nc per point) is stable as one phase, by Michelsen's tangent-plane
   analysis: Wilson's two trial phases and a nearly pure one of each
   compound. It writes stable (n: 1, 0, or -1 for a point without an
   answer), the lowest tangent-plane distance tm any trial reached away from
   the feed (n), negative where the feed is unstable, and that trial's
   composition w (n * nc), the phase that would appear. If every trial came
   back to the feed itself, tm is 0 and w the feed. */
int32_t reflux_stability(int64_t fluid, int64_t n, double *t, int64_t t_len, double *p,
                         int64_t p_len, double *z, int64_t z_rows, int32_t *stable, double *tm,
                         double *w, int32_t *status, char *err, int32_t err_len);

/* Flashes, for reflux_flash: what a point's two values a and b are. H, S
   and U are per mole of feed, on the basis reflux_props describes, and V is
   the molar volume. */
#define REFLUX_FLASH_PT 0  /* a: T (K), b: P (Pa) */
#define REFLUX_FLASH_PH 1  /* a: P (Pa), b: H (J/mol) */
#define REFLUX_FLASH_PS 2  /* a: P (Pa), b: S (J/(mol K)) */
#define REFLUX_FLASH_TV 3  /* a: T (K), b: V (m3/mol) */
#define REFLUX_FLASH_UV 4  /* a: U (J/mol), b: V (m3/mol) */
#define REFLUX_FLASH_TVF 5 /* a: T (K), b: vapour fraction, 0 for a bubble point, 1 for a dew point */
#define REFLUX_FLASH_PVF 6 /* a: P (Pa), b: vapour fraction */

/* A flash's values per point, out[k * n + i]: per mole of feed. */
#define REFLUX_FLASH_BETA 0   /* vapour fraction */
#define REFLUX_FLASH_T 1      /* K */
#define REFLUX_FLASH_P 2      /* Pa */
#define REFLUX_FLASH_H 3      /* J/mol */
#define REFLUX_FLASH_S 4      /* J/(mol K) */
#define REFLUX_FLASH_V 5      /* m3/mol */
#define REFLUX_FLASH_VALUES 6 /* how many there are */

/* A phase's values per point, phase[(k * n + i) * 2 + j], per mole of the
   phase except the fraction, which is of the feed's moles. j is the phase's
   label, REFLUX_LIQUID (0) or REFLUX_VAPOR (1): a split's phase of lower mass
   density is the vapour, and a single phase takes the place of its label,
   the other's values being NaN. */
#define REFLUX_PHASE_FRACTION 0
#define REFLUX_PHASE_V 1      /* m3/mol */
#define REFLUX_PHASE_H 2      /* J/mol */
#define REFLUX_PHASE_S 3      /* J/(mol K) */
#define REFLUX_PHASE_VALUES 4 /* how many there are */

/* reflux_flash flashes n feeds of mole amounts z (nc per point), each at
   its two values a and b, whose meaning spec gives (REFLUX_FLASH_*). For
   each point it writes the REFLUX_FLASH_VALUES values to out, the
   REFLUX_PHASE_VALUES values of each phase to phase (REFLUX_PHASE_VALUES *
   n * 2), each phase's composition to x, as mole fractions (x[(i * 2 + j)
   * nc + c]), and the number of phases, 1 or 2, to n_phases. A flash that
   fails at a point gives that point the status of its kind: no
   convergence, no volume root, more phases than two, and for the
   vapour-fraction flashes phases that came out swapped or a point outside
   the phase envelope. A value it can't reach, such as an H no temperature
   from 10 K to 10000 K gives, is REFLUX_POINT_BAD_INPUT. A vapour-fraction
   flash gives both phases, even one with no moles at a vapour fraction of 0
   or 1: its composition is the incipient phase's. */
int32_t reflux_flash(int64_t fluid, int32_t spec, int64_t n, double *a, int64_t a_len, double *b,
                     int64_t b_len, double *z, int64_t z_rows, double *out, double *phase,
                     double *x, int32_t *n_phases, int32_t *status, char *err, int32_t err_len);

#endif
