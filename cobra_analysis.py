import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

model_path = "/work/lylab/cjn40747/metabolome/F7_5.xml"
model = cobra.io.read_sbml_model(model_path)

GG_medium = {
    "EX_glc__D_e": 55.5,
    "EX_glu__L_e": 34.0,
    "EX_mg2_e": 0.811,
    "EX_so4_e": 0.811,
    "EX_k_e": 3.67,
    "EX_pi_e": 3.67,
    "EX_na1_e": 3.42,
    "EX_cl_e": 3.42,
    "EX_ca2_e": 49.9,
    "EX_h_e": 1000
}

def set_GG_medium(model, include_glucose=True):
    medium = GG_medium.copy()
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
        fraction_of_optimum=0.9
    )
    secretions = {}
    for rxn in model.exchanges:
        max_flux = fva.loc[rxn.id, "maximum"]
        if max_flux > 1e-6:
            met = list(rxn.metabolites.keys())[0]
            secretions[met.id] = max_flux
    return secretions

if __name__ == "__main__":
    model = cobra.io.read_sbml_model(model_path)
    with model:
        set_GG_medium(model, include_glucose=False)
        secretions_no_glc = find_secretions(model)

    with model:
        set_GG_medium(model, include_glucose=True)
        secretions_with_glc = find_secretions(model)

    df = pd.DataFrame({
        "GG_max_secretion": secretions_no_glc,
        "GG_glucose_max_secretion": secretions_with_glc
    }).fillna(0)

    df["difference"] = df["GG_glucose_max_secretion"] - df["GG_max_secretion"]

    df = df.sort_values("difference", ascending=False)

    names = {m.id: m.name for m in model.metabolites}
    df["metabolite_name"] = df.index.map(names)

    df = df[[
        "metabolite_name",
        "GG_max_secretion",
        "GG_glucose_max_secretion",
        "difference"
    ]]

    print(df.head(30))

    df.to_csv("/work/lylab/cjn40747/metabolome/secreted_metabolites_FVA.csv")