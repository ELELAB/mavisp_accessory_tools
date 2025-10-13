import itertools
import pandas as pd
import random

def random_mutation():
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    return f"{random.choice(aa)}{random.randint(1, 50)}{random.choice(aa)}"

def demask_value(label):
    return random.uniform(-0.5, -0.25) if label == "Damaging" else random.uniform(-0.249, 0.249)

def gemme_value(label):
    return random.uniform(-5, -3) if label == "Damaging" else random.uniform(-2.9, 2.9)

def dhfr_pca_value(label):
    if label == "Damaging":
        return random.uniform(-1, -0.161)
    elif label == "Neutral":
        return random.uniform(-0.161, 1)
    else:
        return None

# Colonne
columns = [
    "Stability classification, alphafold, (RaSP, FoldX)",           # binaria
    "Experimental data (DHFR-PCA yeast, Stability (MAVEdb))",      # numerica
    "Experimental data classification (DHFR-PCA yeast, Stability (MAVEdb))",  # binaria
    "DeMaSk delta fitness",                                         # binaria -> numerica
    "GEMME Score",                                                  # binaria -> numerica
    "Local Int. classification (PMS2_AFmulti)"                      # binaria
]

classes = ["Damaging", "Neutral"]

rows = []

# Generiamo tutte le combinazioni delle 5 colonne binarie (tranne la numerica)
for combo in itertools.product(classes, repeat=5):
    mutation = random_mutation()
    exp_class = combo[2]  # colonna di classificazione DHFR
    row = {
        "Mutation": mutation,
        columns[0]: combo[0],
        columns[2]: exp_class,
        columns[1]: dhfr_pca_value(exp_class),  # numero coerente con la classificazione
        columns[3]: demask_value(combo[3]),
        columns[4]: gemme_value(combo[4]),
        columns[5]: combo[4]  # Local Int. puoi usare lo stesso combo[4] o random
    }
    rows.append(row)

df = pd.DataFrame(rows)
df["Relative Side Chain Solvent Accessibility in wild-type"] = 20
df.to_csv("test-simple_mode.csv", index=False)
print(f"Generated CSV with {len(df)} simulated mutations.")

