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

**Model parameters** (`pr.ipd`, `srk.ipd`, `mc-PR.lib`). These files say
"Copyright 2009 Harry Kooijman and Ross Taylor" (2023 for `mc-PR.lib`) but
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
the first copy is the more complete and more recently edited one, and it's
the one reflux keeps. The import report lists everything else the importer
noticed: summed ASOG groups, impossible correlation ranges and similar.

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
- `ipd/` also has NRTL, UNIQUAC, Wilson and UNIFAC parameter tables, for M6.

Like the `.ipd` files, these name no licence.
