#!/usr/bin/env python3

# Copyright (C) 2025 Eszter Toldi
# Technical University of Denmark, Danish Cancer Institute

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd

# Load fireprot dataset
df = pd.read_csv("fireprot_train.csv")

columns_keep = ["pdb_id_corrected", "uniprot_id", "sequence", "protein_name"]
df_small = df[columns_keep]
df_unique = df_small.drop_duplicates().reset_index(drop=True)

# Save to a new file
df_unique.to_csv("fireprot_clean.csv", index=False)

print(f"Saved {len(df_unique)} unique entries to fireprot_clean.csv")
