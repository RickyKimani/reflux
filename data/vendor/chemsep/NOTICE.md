# ChemSep data: source, licence, and what reflux does with it

## Source

The files in this folder are unmodified copies from **ChemSep LITE 8.50**,
installed from chemsep.org on 2026-10-26. Pure-component data came from
the installation's `pcd` folder, and interaction parameters from its `ipd`
folder. They are byte-for-byte identical to the installed files.

| File | From | Version (as stated in the file) | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| `chemsep1.xml` | `pcd/` | ChemSep v8.49 pure component data 1 (2024), 433 compounds | 3818631 | `453b0f61ec0ab185eff2b1f8096aa1e96b05feb1c67f9b03e42a15dcd800ab14` |
| `chemsep2.xml` | `pcd/` | ChemSep v8.50 pure component data 2 (2025, modified 2026), 211 records of 209 compounds | 1604493 | `37c7086eb4c4b89a61406fe565c66720b7b31d5ae5cdeef784e143fa5031607e` |
| `pr.ipd` | `ipd/` | DECHEMA Peng-Robinson parameters | 21929 | `48440f2838e6fa8475e32b59756891c23ea0ea0a58ad56d4a3d79e32089e9d27` |
| `srk.ipd` | `ipd/` | DECHEMA Soave-Redlich-Kwong parameters | 22449 | `d2a94c1566b99bc08847c007f6f4e45dcce9b334098f433d2ef1d37908665d44` |
| `mc-PR.lib` | `pcd/` | Matthias-Copeman Peng-Robinson (2023) | 33138 | `47b498a583037c752da4a74ac1924b73d6046a95ea4926b294b85be9221c461f` |
| `nrtl.ipd` | `ipd/` | DECHEMA NRTL parameters at 1 atm, cal/mol (2008) | 30351 | `6b8c3404ca0ebfddbc99ccbe97bd88df2dfbdc979dc31af160994b405659f8ae` |
| `uniquac.ipd` | `ipd/` | DECHEMA UNIQUAC parameters at 1 atm, cal/mol (2008) | 22763 | `8eb18d212d2fec49cb9cea65314df0a5ca71cf658aa4022b290a3eccf78f7a5b` |
| `wilson.ipd` | `ipd/` | Wilson parameters at 1 atm, cal/mol, Andrew Mason's compilation (2008) | 31229 | `b29365910fab307eaa7ae102e4647055aa2cc69b04db46638e9c16014227f6a9` |
| `henrys_law.ipd` | `ipd/` | Henry's law constants in water, from the API Technical Data Book (2011) | 5308 | `ab632f345ceeff04bbbc7b44641e75ccff9dfb83658beae8db3af8f3fb659974` |
| `unifacvl.gct` | `ipd/` | UNIFAC (VLE) main groups and interaction parameters | 18995 | `7cc0fb2790837df1887e14e1f777eae82895e72b5d9e89e7493c7d3541867132` |
| `unifacvl_ext.gct` | `ipd/` | its extension: new main groups and parameters | 4636 | `906c3afdf30375dc4868b860f01eebd2cb56c254bdb9791bd39668a387b8f699` |
| `unifacrq.gct` | `ipd/` | UNIFAC (VLE) subgroups: main group, R, Q, M | 4405 | `06b4cbfaf8376da2d69b5b390bc70ab97cbee16e2f953c383fa65c33bb10ae7a` |
| `unifacrq_ext.gct` | `ipd/` | its extension: new subgroups | 1044 | `dbba6796d829c1633e1eca8da84061800632c51d8a4b67a3eb3bb3d676e4b702` |
| `unifacll.gct` | `ipd/` | UNIFAC (LLE) main groups and interaction parameters | 7607 | `163de3ebfde5363ce560bd9a5cccd23864ae1f844ca4b738ee9f64143b630ce7` |
| `unifacll_ext.gct` | `ipd/` | its extension | 111 | `edb7b589608e9f3dd3db8888f3d7987598fe82fdfcf45147c864c4aaf2619a80` |
| `llefacrq.gct` | `ipd/` | UNIFAC (LLE) subgroups | 2287 | `1f8a6cf99b82b7e09e60fc3b5e9c408b86ea890c73d155fc8075e433911f49bd` |
| `llefacrq_ext.gct` | `ipd/` | its extension | 100 | `7e472c3e92c30060bdacb1b3a7f26a8cbc93f6415a88dc5e989e0e5a8aa210c1` |
| `modfac_a.gct` | `ipd/` | modified UNIFAC (Dortmund) main groups and a parameters | 23971 | `241dea971728a79cd8e54404c7951c4a25e6ed47dc87cc4ab2138358d3f42b0d` |
| `modfac_a_ext.gct` | `ipd/` | its extension | 400 | `473169418c45b91feb9f0d749d723d4d5489fc2693fed458a90e7a0a2469b583` |
| `modfac_b.gct` | `ipd/` | modified UNIFAC (Dortmund) b parameters | 19302 | `81c7ec45b113f2d550123af3f5e40ef1a57240d02cff45052abfafcf6785e253` |
| `modfac_b_ext.gct` | `ipd/` | its extension | 387 | `65dfad5b0a20959ea5c622a56e2ae42260fe46db4b52014b77f554ab3205b65a` |
| `modfac_c.gct` | `ipd/` | modified UNIFAC (Dortmund) c parameters | 7050 | `f9bb69da0967aba7ecb5364ab783100f7e397ac701c885868a532c65590d86d2` |
| `modfac_c_ext.gct` | `ipd/` | its extension | 345 | `18a1532bcdabb94723a463501c3ca5a7d41f0d819a96864471854a90ac73d6a8` |
| `modfacrq.gct` | `ipd/` | modified UNIFAC (Dortmund) subgroups | 5066 | `d6d53e221970579bf75fb88078dc2b929dd328538522eb361ecab21ba19dd21a` |
| `modfacrq_ext.gct` | `ipd/` | its extension | 175 | `bc6a0a923c2bb116a4ca8087a5fab92206ae71b2c12224b2c0462028ac5e571a` |
| `Artistic_license_2_0.txt` | `pcd/` | the licence text ChemSep ships with its data | 8720 | `da4a665b492e71021b516e2d9e261f77f0c62ac60bfac089b4a4a33ddecd42a7` |

`.gitattributes` marks this folder `-text`, so Git never rewrites the line
endings and these checksums hold on every checkout.

## Licence

**Pure-component data** (`chemsep1.xml`, `chemsep2.xml`). Each file states
its copyright and licence in its root element:

> ChemSep v8.49 pure component data 1 - Copyright (c) Harry Kooijman and
> Ross Taylor (2024) - http://www.perlfoundation.org/artistic_license_2_0

> ChemSep v8.50 pure component data 2 - Copyright (c) Harry Kooijman and
> Ross Taylor (2025) - http://www.perlfoundation.org/artistic_license_2_0

The Artistic License 2.0 is in `Artistic_license_2_0.txt`, as distributed
with the data. These copies are the Standard Version, redistributed
verbatim under §2 of that licence. reflux's converted database (in
`data/components/`) is a Modified Version under §4, distributed under the
same licence, with its changes documented next to it. reflux's code is MIT
licensed (`LICENSE` at the repository root). It's a separate work that
reads this data (§9), aggregated with it under §7.

**Model parameters** (`pr.ipd`, `srk.ipd`, `mc-PR.lib`, and since M6 step 1
the NRTL, UNIQUAC, Wilson and Henry's-law banks and the UNIFAC tables).
These files say "Copyright 2009 Harry Kooijman and Ross Taylor" (2008 or
2011 for some, 2023 for `mc-PR.lib`; the `.gct` tables say nothing) but
**name no licence**.
reflux assumes they're part of the same ChemSep data package, and so
under the same Artistic License 2.0, as DWSIM also does when it
redistributes them. **This is an assumption, not something the files
state.** It should be confirmed with the authors, and the files removed if
it turns out to be wrong.

## Known defects, left as shipped

The XML files are **not well-formed**: 43 SMILES attribute values end in a
stray `<` (1 in `chemsep1.xml`, 42 in `chemsep2.xml`; for example
`CC(C)(C)C1CCCCC1<`). The copies here keep them, because they're verbatim.
reflux's importer drops the trailing `<`, which is invalid SMILES in any
case, and lists every repair in its import report.

Two compounds are in `chemsep2.xml` twice: trifluoromethane (75-46-7) and
3-chloropropene (107-05-1), with the same library index each time. In both,
the first copy is the more complete one, and it's the one reflux keeps,
with one exception. R-23's first-copy vapour pressure is 94–100% below
CoolProp's across its range, where the second copy's is within 3%. So
reflux takes its vapour pressure and Antoine fit from the second copy. The
import report lists that, and everything else the importer noticed: summed
ASOG groups, impossible correlation ranges and similar.

The parameter files have defects of their own:
- **Repeated pairs.** In `pr.ipd` and `srk.ipd`, about a quarter of the
  pairs are listed more than once, each line a separate DECHEMA fit with
  its own value. ChemSep asks you to pick one; reflux takes the first
  listed and reports the others.
- **`mc-PR.lib` lists 60 compounds more than once,** some four times.
  reflux takes the first.
- **A line with no CAS number.** One `mc-PR.lib` line (propane) is missing
  it. reflux skips that line and reports it.
- **`*` in place of numbers.** reflux reads it as a missing value.
- **The activity banks (M6 step 1)** state their units on a `Units=` line,
  cal/mol, and their formulas nowhere. Their A's are over DECHEMA's own R,
  1.98721 cal/(mol·K): DWSIM's NRTL on the same bank is τ = A/(1.98721·T)
  to 4e-10 in ln γ (docs/spike-dwsim-activity.md).
- **UNIFAC's numbering** follows the 1987 *Properties of Gases and Liquids*:
  ChemSep's subgroup 9 is C=C, so the subgroups others number 9 to 69 are
  its 10 to 70, as `unifacrq.gct`'s first line says.
- **`modfac_a_ext.gct`** names its one new main group, `cy-CHNH2`, with no
  number. reflux skips that line and reports it.

## Compared with the copies in DWSIM's repository

DWSIM's GitHub repository ships older exports (v8.3 and v8.31, 2022), with
431 and 127 compounds. Compared with those:

- **Compounds:** `chemsep1` adds 2 compounds; `chemsep2` adds 84 and drops
  2. The two files no longer overlap; there are 642 distinct compounds in
  all.
- **New identifiers:** `InChIKey` and `PubChemID` elements. The SMILES
  element is now spelled `SMILES`, not `Smiles`.
- **Values:** most differences are reformatting, with a median relative
  change of about 1e-6 (e.g. 58.1222 written as 58.122). A few dozen
  compounds have real revisions to Tc, Pc or ω. For example, acrylic acid's
  Tc goes from 654.2 K to 609.5 K, maleic anhydride's Pc from 72.8 to 59.9
  bar, and several glycols' ω roughly double.

`mc-PR.lib` (vendored in M2 step 6) holds Mathias–Copeman parameters
fitted for Peng–Robinson. For the 573 compounds that also have
Mathias–Copeman constants in the XML, none of the sets agree. The XML's are
a different, unlabelled fit, which is why reflux doesn't use them by
default.

## Also in the ChemSep installation, not vendored yet

- `pcd/chemsep.syn` holds name synonyms, which would help lookup by alias.
- `ipd/apisrk.ipd` has API-SRK kᵢⱼ.
- `ipd/` also has Margules, van Laar and Hayden–O'Connell parameters, the
  q′ form of UNIQUAC, ASOG's and PSRK's tables, and UMR's.

Like the `.ipd` files, these name no licence.
