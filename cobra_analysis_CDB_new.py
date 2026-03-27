import cobra
import pandas as pd
import os
from cobra.flux_analysis import pfba

metadata = pd.read_csv("/home/cjn40747/metabolome_map/model_metadata.csv")

# ----------------------------
# MEDIA (same as yours)
# ----------------------------
CDB_medium = {
    "EX_glc__D_e": 10,
    "EX_k_e": 1000,
    "EX_pi_e": 1000,
    "EX_mg2_e": 1000,
    "EX_so4_e": 1000,
    "EX_h2o_e": 1000,
    "EX_na1_e": 1000,
    "EX_cl_e": 1000,
    "EX_nh4_e": 1000,
    "EX_ca2_e": 1000,
    "EX_hco3_e": 1000,
    "EX_ala__L_e": 5,
    "EX_arg__L_e": 5,
    "EX_asn__L_e": 5,
    "EX_asp__L_e": 5,
    "EX_glu__L_e": 5,
    "EX_gly_e": 5,
    "EX_his__L_e": 5,
    "EX_ile__L_e": 5, 
    "EX_leu__L_e": 5,
    "EX_lys__L_e": 5,
    "EX_met__L_e": 5,
    "EX_phe__L_e": 5,
    "EX_pro__L_e": 5,
    "EX_ser__L_e": 5,
    "EX_thr__L_e": 5,
    "EX_trp__L_e": 5,
    "EX_tyr__L_e": 5,
    "EX_val__L_e": 5,
    "EX_thm_e": 0.1,
    "EX_ribflv_e": 0.1,
    "EX_btn_e": 0.01,
    "EX_fol_e": 0.01,
    "EX_ade_e": 1,
    "EX_gua_e": 1,
    "EX_uri_e": 1,
    "EX_o2_e": 1000,
    "EX_cobalt2_e": 1000,
    "EX_cu2_e": 1000,
    "EX_fe2_e": 1000,
    "EX_fe3_e": 1000,
    "EX_mn2_e": 1000,
    "EX_mobd_e": 1000,
    "EX_ni2_e": 1000,
    "EX_zn2_e": 1000
}

# ----------------------------
# MEDIUM SETTING (CORRECT)
# ----------------------------
def set_medium(model, include_glucose=True):

    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 0

    for rxn_id, value in CDB_medium.items():
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = -value
            rxn.upper_bound = 1000

    if not include_glucose and "EX_glc__D_e" in model.reactions:
        model.reactions.EX_glc__D_e.lower_bound = 0


# ----------------------------
# GROWTH CONSTRAINT
# ----------------------------
def constrain_growth(model):
    sol = model.optimize()
    growth = sol.objective_value

    if growth < 1e-6:
        return None

    biomass = list(model.objective.variables)[0].name
    biomass_rxn = model.reactions.get_by_id(biomass)

    biomass_rxn.lower_bound = growth * 0.9
    biomass_rxn.upper_bound = growth * 0.9

    return growth


# ----------------------------
# GET SECRETIONS
# ----------------------------
def get_secretions(model):
    sol = pfba(model)
    sec = {}

    for rxn in model.exchanges:
        flux = sol.fluxes[rxn.id]
        if flux > 1e-6:
            met = list(rxn.metabolites.keys())[0]
            sec[met.id] = flux

    return sec


# ----------------------------
# RUN MODEL
# ----------------------------
def run_model(model_path):

    model = cobra.io.read_sbml_model(model_path)

    with model:
        set_medium(model, False)
        if constrain_growth(model) is None:
            return None
        no_glc = get_secretions(model)

    with model:
        set_medium(model, True)
        if constrain_growth(model) is None:
            return None
        glc = get_secretions(model)

    df = pd.DataFrame({
        "no_glc": no_glc,
        "glc": glc
    }).fillna(0)

    df["delta"] = df["glc"] - df["no_glc"]

    return df["delta"]


# ----------------------------
# COLLECT ALL MODELS
# ----------------------------
all_data = []

for _, row in metadata.iterrows():
    print("Processing:", row["model_path"])

    delta = run_model(row["model_path"])

    if delta is None:
        print("Skipped (no growth)")
        continue

    df = pd.DataFrame(delta)
    df.columns = ["delta"]
    df["label"] = row["label"]

    all_data.append(df)

combined = pd.concat(all_data)
combined.reset_index(inplace=True)
combined.rename(columns={"index": "metabolite"}, inplace=True)

# ----------------------------
# GROUP ANALYSIS
# ----------------------------
summary = combined.groupby(["metabolite", "label"])["delta"].agg([
    "mean",
    "median",
    "count"
]).reset_index()

pivot_mean = summary.pivot(index="metabolite", columns="label", values="mean").fillna(0)
pivot_count = summary.pivot(index="metabolite", columns="label", values="count").fillna(0)

# prevalence = how many models secrete it
pivot_prev = pivot_count / metadata["label"].value_counts()

# ----------------------------
# FINAL MERGED TABLE
# ----------------------------
final = pd.DataFrame(index=pivot_mean.index)

final["dissolver_mean"] = pivot_mean.get("dissolver", 0)
final["non_dissolver_mean"] = pivot_mean.get("non_dissolver", 0)

final["dissolver_prevalence"] = pivot_prev.get("dissolver", 0)
final["non_dissolver_prevalence"] = pivot_prev.get("non_dissolver", 0)

final["difference"] = final["dissolver_mean"] - final["non_dissolver_mean"]

# ----------------------------
# FILTER: SIGNATURES
# ----------------------------
signatures = final[
    (final["dissolver_mean"] > 1) &
    (final["difference"] > 1) &
    (final["dissolver_prevalence"] > 0.7) &
    (final["non_dissolver_prevalence"] < 0.3)
]

final.to_csv("/work/lylab/cjn40747/metabolome/all_metabolites_comparison.csv")
signatures.to_csv("/work/lylab/cjn40747/metabolome/calcite_dissolution_signatures.csv")

print("\nTop candidate metabolites:")
print(signatures.sort_values("difference", ascending=False).head(20))