# Component data

One TOML file per compound, named by its CAS number, in reflux's component
format ([docs/data-format.md](../../docs/data-format.md)). `index.toml` lists
them all, so lookup by name doesn't have to parse every file.

**Everything here is generated. Don't edit it.**
[tools/import_chemsep.vsh](../../tools/import_chemsep.vsh) writes this whole
directory from the ChemSep XML in [data/vendor/chemsep/](../vendor/chemsep/),
and a test fails if it's out of date. To change a compound, fix the importer,
or add your own component file in a directory of your own.

## Licence

These files are a Modified Version of the ChemSep pure-component database,
Copyright (c) Harry Kooijman and Ross Taylor. They are distributed under the
same licence as the original, the Artistic License 2.0
([text](../vendor/chemsep/Artistic_license_2_0.txt)). The Standard Version is
vendored, unmodified, in [data/vendor/chemsep/](../vendor/chemsep/); see its
[NOTICE](../vendor/chemsep/NOTICE.md).

How these files differ from it is documented in [CHANGES.md](CHANGES.md).
