import cobra
import pandas as pd
import sys
import os
from cobra.flux_analysis import pfba

model_path = sys.argv[1]
model = cobra.io.read_sbml_model(model_path)
base = os.path.basename(model_path.replace(".xml", ""))


# ----------------------------
# MEDIA
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
    "EX_o2_e": 1000,
}


def set_medium(model, include_glucose=True):
    m = CDB_medium.copy()
    if not include_glucose:
        m["EX_glc__D_e"] = 0
    model.medium = m


# ----------------------------
# FORCE SAME GROWTH STATE
# ----------------------------
def constrain_growth(model, fraction=0.9):
    sol = model.optimize()

    biomass_rxns = [r for r in model.reactions if "biomass" in r.id.lower()]
    if not biomass_rxns:
        raise ValueError("No biomass reaction found")

    biomass = biomass_rxns[0]

    biomass.lower_bound = sol.objective_value * fraction
    biomass.upper_bound = sol.objective_value * fraction

    return sol.objective_value


# ----------------------------
# GET SECRETION (pFBA)
# ----------------------------
def get_secretions(model):
    sol = pfba(model)

    secretions = {}

    for rxn in model.exchanges:
        flux = sol.fluxes[rxn.id]

        # secretion convention check (COBRA: usually negative = export)
        if abs(flux) > 1e-6:
            met = list(rxn.metabolites.keys())[0]
            secretions[met.id] = flux

    return secretions


# ----------------------------
# RUN CONDITION
# ----------------------------
def run_condition(model, include_glc):
    with model:
        set_medium(model, include_glc)

        # optional gas constraints (keep consistent but minimal)
        model.reactions.EX_o2_e.lower_bound = -10
        model.reactions.EX_co2_e.lower_bound = -5
        model.reactions.EX_hco3_e.lower_bound = -5

        constrain_growth(model)

        return get_secretions(model)


# ----------------------------
# MAIN
# ----------------------------
with model:
    secretions_no_glc = run_condition(model, False)

with model:
    secretions_glc = run_condition(model, True)


# ----------------------------
# BUILD DATAFRAME
# ----------------------------
df = pd.DataFrame({
    "no_glucose": secretions_no_glc,
    "glucose": secretions_glc
}).fillna(0)

df["difference"] = df["glucose"] - df["no_glucose"]

df["fold_change"] = (df["glucose"] + 1e-9) / (df["no_glucose"] + 1e-9)

df = df.sort_values("difference", ascending=False)

# metabolite names
names = {m.id: m.name for m in model.metabolites}
df["metabolite_name"] = df.index.map(names)

df = df[[
    "metabolite_name",
    "no_glucose",
    "glucose",
    "difference",
    "fold_change"
]]

# ----------------------------
# FILTERING (more meaningful now)
# ----------------------------
df_filtered = df[
    (abs(df["difference"]) > 1e-3) &
    (df["fold_change"] > 2)
]

df.to_csv(f"/work/lylab/cjn40747/metabolome/{base}_PFBA_secretion.csv")
df_filtered.to_csv(f"/work/lylab/cjn40747/metabolome/{base}_PFBA_secretion_filtered.csv")