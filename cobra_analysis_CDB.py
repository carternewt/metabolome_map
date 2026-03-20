import cobra
import pandas as pd
import sys
import os
from cobra.flux_analysis import flux_variability_analysis

model_path = sys.argv[1]
model = cobra.io.read_sbml_model(model_path)
base=os.path.basename(model_path.replace(".xml", ""))

GG_medium = {
    "EX_glc__D_e": 10,
    "EX_glu__L_e": 5,
    "EX_mg2_e": 1000,
    "EX_so4_e": 1000,
    "EX_k_e": 1000,
    "EX_pi_e": 1000,
    "EX_na1_e": 1000,
    "EX_cl_e": 1000,
    "EX_ca2_e": 1000,
    "EX_h2o_e": 1000,
    "EX_co2_e": 1000,
    "EX_nh4_e": 1000,
    "EX_hco3_e": 1000,
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

def set_CDB_medium(model, include_glucose=True):
    medium = CDB_medium.copy()
    if not include_glucose:
        medium["EX_glc__D_e"] = 0
    model.medium = medium

def find_secretions(model):
    sol = model.optimize()
    growth = sol.objective_value
    print("Growth rate:", growth)
    fva = flux_variability_analysis(
        model,
        reaction_list=[rxn.id for rxn in model.exchanges],
        fraction_of_optimum=0.99
    )
    secretions = {}
    for rxn in model.exchanges:
        max_flux = fva.loc[rxn.id, "maximum"]
        if max_flux > 1e-6:
            met = list(rxn.metabolites.keys())[0]
            secretions[met.id] = max_flux
    return secretions

def find_secretions_fba(model):
    sol = model.optimize()
    growth = sol.objective_value
    print("FBA Growth rate:", growth)
    secretions = {}
    for rxn in model.exchanges:
        flux = sol.fluxes[rxn.id]
        if flux > 1e-3:
            met = list(rxn.metabolites.keys())[0]
            secretions[met.id] = flux
    return secretions

if __name__ == "__main__":
    model = cobra.io.read_sbml_model(model_path)
    with model:
        set_CDB_medium(model, include_glucose=False)
        model.reactions.EX_o2_e.lower_bound = -10
        model.reactions.EX_h_e.lower_bound = -10
        model.reactions.EX_h_e.upper_bound = 10
        model.reactions.EX_co2_e.lower_bound = -5
        model.reactions.EX_hco3_e.lower_bound = -5
        secretions_no_glc = find_secretions(model)

    with model:
        set_CDB_medium(model, include_glucose=True)
        model.reactions.EX_o2_e.lower_bound = -10
        model.reactions.EX_h_e.lower_bound = -10
        model.reactions.EX_h_e.upper_bound = 10
        model.reactions.EX_co2_e.lower_bound = -5
        model.reactions.EX_hco3_e.lower_bound = -5
        secretions_with_glc = find_secretions(model)
    
    with model:
        set_CDB_medium(model, include_glucose=False)
        model.reactions.EX_o2_e.lower_bound = -10
        model.reactions.EX_h_e.lower_bound = -10
        model.reactions.EX_h_e.upper_bound = 10
        model.reactions.EX_co2_e.lower_bound = -5
        model.reactions.EX_hco3_e.lower_bound = -5
        fba_no_glc = find_secretions_fba(model)

    with model:
        set_CDB_medium(model, include_glucose=True)
        model.reactions.EX_o2_e.lower_bound = -10
        model.reactions.EX_h_e.lower_bound = -10
        model.reactions.EX_h_e.upper_bound = 10
        model.reactions.EX_co2_e.lower_bound = -5
        model.reactions.EX_hco3_e.lower_bound = -5
        fba_with_glc = find_secretions_fba(model)

    df = pd.DataFrame({
        "CDB_max_secretion": secretions_no_glc,
        "CDB_glucose_max_secretion": secretions_with_glc
    }).fillna(0)

    df["difference"] = df["CDB_glucose_max_secretion"] - df["CDB_max_secretion"]

    df = df.sort_values("difference", ascending=False)

    names = {m.id: m.name for m in model.metabolites}
    df["metabolite_name"] = df.index.map(names)

    df = df[[
        "metabolite_name",
        "CDB_max_secretion",
        "CDB_glucose_max_secretion",
        "difference"
    ]]

    df["fold_change"] = (df["CDB_glucose_max_secretion"] + 1e-6) / (df["CDB_max_secretion"] + 1e-6)

    df_fva_filtered = df[(df["difference"] >= 4) & (df["CDB_glucose_max_secretion"] >= 2) & (df["fold_change"] >= 2)]

    df_fva_filtered.to_csv("/work/lylab/cjn40747/metabolome/{base}_filtered_CDB.csv")
    df.to_csv("/work/lylab/cjn40747/metabolome/{base}_CDB.csv")

    df_fba = pd.DataFrame({
    "CDB_FBA_secretion": fba_no_glc,
    "CDB_glucose_FBA_secretion": fba_with_glc
    }).fillna(0)

    df_fba["difference"] = df_fba["CDB_glucose_FBA_secretion"] - df_fba["CDB_FBA_secretion"]

    df_fba = df_fba.sort_values("difference", ascending=False)

    names = {m.id: m.name for m in model.metabolites}
    df_fba["metabolite_name"] = df_fba.index.map(names)

    df_fba = df_fba[[
        "metabolite_name",
        "CDB_FBA_secretion",
        "CDB_glucose_FBA_secretion",
        "difference"
    ]]

    df_fba_filtered = df_fba[
    (df_fba["difference"] >= 2) &
    (df_fba["CDB_glucose_FBA_secretion"] >= 1)]

    df_fba_filtered.to_csv("/work/lylab/cjn40747/metabolome/{base}_filtered_CDB.csv")
    df_fba.to_csv("/work/lylab/cjn40747/metabolome/{base}_CDB.csv")