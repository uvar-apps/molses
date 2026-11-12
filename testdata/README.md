# Test data

All structure files in this directory are **CC0 1.0** (public domain) — they
come from the [wwPDB archive](https://www.wwpdb.org/about/usage-policies),
which dedicates PDB data files to the public domain. Please still attribute the
original structure authors when publishing results (wwPDB requests this as
scientific courtesy).

| File | PDB ID | Description | Why |
|---|---|---|---|
| `1crn.pdb` / `1crn.cif` | 1CRN | Crambin (46 residues, 327 atoms) | Classic MSMS benchmark [Sanner et al. 1996] |
| `4hhb.pdb` | 4HHB | Hemoglobin | Multi-chain, altLocs |
| `synthetic.xyz` | — | 3-atom water | Synthetic (generated), no license needed |
| `spheres.xyzr` | — | 2 spheres | Synthetic MSMS interchange fixture |

Download with:

```console
$ curl -L -o testdata/1crn.pdb https://files.rcsb.org/download/1CRN.pdb
$ curl -L -o testdata/1crn.cif https://files.rcsb.org/download/1CRN.cif
$ curl -L -o testdata/4hhb.pdb https://files.rcsb.org/download/4HHB.pdb
```

If a file is missing from this directory, tests that need it skip with a
message rather than failing.
