# MAVISp Mechanistic Indicators 
This script identifies the mechanistic indicators identified in MAVISp for a given variant.

The script reads a list of proteins and their respective variants, and accesses the final MAVISp CSV of each protein across all available modes of analysis (simple/ensemble) and reports the predicted identified mechanisms associated with each variant in broad categories. 
The aim of the script is to streamline the analysis of a larger pool of variants, to be able to have an overview of their effects without consulting the full MAVISp CSV. 

## Features
- Searches all available MAVISp modes (simple, ensemble) in the given MAVISp database version.
- Maps MAVISp annotations to standardized mechanistic indicator categories.
- Reports DeMaSk-based **Loss-of-Function (LoF)** or **Gain-of-Function (GoF)** predictions when no other mechanism is identified.
- Produces a single CSV containing mechanistic annotations from every available MAVISp mode.

## Requirements
- Python

## Input
The input must be a CSV containing at least the following columns:

| Column | Description |
|---------|-------------|
| `protein` | Protein identifier (e.g. `ARID3A`) |
| `variant` | Protein variant in HGVSp format (e.g. `p.Arg175His`) |

Example:
```csv
protein,variant
TP53,p.Arg175His
BRCA1,p.Cys61Gly
PTEN,p.Arg130Gln
```

---

## MAVISp Database Structure
The script expects a MAVISp database organized as:

```
mavisp_db/
│
├── simple_mode/
│   └── dataset_tables/
│       ├── TP53-simple_mode.csv
│       ├── BRCA1-simple_mode.csv
│       └── ...
│
├── ensemble_mode/
│   └── dataset_tables/
│       ├── TP53-ensemble_mode.csv
│       ├── BRCA1-ensemble_mode.csv
│       └── ...
```

Each dataset table must contain an `HGVSp` column.

---

## Usage
```bash
python mechanistic_indicators.py \
    -i variants.csv \
    -m /path/to/mavisp_db \
    -o results.csv
```

### Arguments
| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Input CSV containing proteins and variants |
| `-m`, `--mavisp_db` | Root directory of the MAVISp database |
| `-o`, `--output` | Output CSV (default: `mechanistic_indicators.csv`) |
| `-d`, `--demask_threshold` | Threshold for assigning LoF/GoF from DeMaSk scores (default: `0.25`) |

Example:

```bash
python mechanistic_indicators.py \
    -i variants.csv \
    -m ./MAVISp_DB \
    -o mechanisms.csv \
    -d 0.3
```

---

## Output

The output CSV contains one row per input variant.

Example:

| protein | variant | simple_mode | ensemble_mode |
|----------|----------|-------------|---------------|
| TP53 | p.Arg175His | STABILITY | STABILITY_LONG RANGE |
| BRCA1 | p.Cys61Gly | SS-BOND LOSS | SS-BOND LOSS |
| PTEN | p.Arg130Gln | LoF | LoF |

If a variant is not found in a MAVISp mode, the value will be:

```
NA
```

If a variant is found but no mechanism is detected, the value will be:

```
NONE
```

---

## Mechanism Mapping

The script converts MAVISp predictions and annotations into broader mechanism categories.

- STABILITY
- PTM STABILITY 
- LONG RANGE 
- LOCAL INTERACTION 
- FUNCTIONAL SITE 
- EARLY FOLDING
- SS-BOND GAIN 
- SS-BOND LOSS 
- PHOSPHO LOSS 
- DENOVO PHOSPHO 

If the column **damaging conformations** contains a positive value, the following mechanism is added:
- STABILITY CONFORMATION DEPENDENT


Multiple mechanisms are concatenated using underscores (`_`).

Example:

```
LONG RANGE_STABILITY
```

---

## LoF/GoF Assignment

If **no specific mechanism** is detected, the script examines the **DeMaSk delta fitness** score.

Given threshold `d`:

- score ≤ `-d` → `LoF`
- score ≥ `d` → `GoF`
- otherwise → `NONE`

The default threshold is:

```
0.25
```

This can be changed using:

```bash
-d 0.5
```

---

## Logging

The script reports missing data using Python logging.

Typical warnings include:

- missing dataset files
- missing `HGVSp` column
- variants not found in a dataset
- proteins absent from all prediction modes

---

## Example Run 
cd example \
module load python \
python ../mech_ids.py -i test.csv -m 29052026_ALL 
