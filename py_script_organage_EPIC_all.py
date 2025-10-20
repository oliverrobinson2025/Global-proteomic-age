from organage import OrganAge
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy.interpolate import interp1d
from scipy import stats
import numpy as np

EPIC4ND = pd.read_csv("Input_prot_clocks_all_age_sex.csv")
EPIC4ND = EPIC4ND.set_index("ID")
print(EPIC4ND)

df_prot = pd.read_csv("Input_prot_clocks_all.csv")
df_prot = df_prot.set_index("ID")
print(df_prot)

data = OrganAge.CreateOrganAgeObject()
data.add_data(EPIC4ND, df_prot)
data.normalize(assay_version="v4.1")  #requires "v4" 5k, or "v4.1" 7k
res = data.estimate_organ_ages()
print(res)

res.to_csv("output/Results_EPIC_ALL_ageblood_organage.csv")
