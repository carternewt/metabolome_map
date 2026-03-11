import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

model_path = "/work/lylab/cjn40747/metabolome/F7_5.xml"
model = cobra.io.read_sbml_model(model_path)

def set_lb_medium(model, glucose=False):
    medium = model.medium.copy()
    if "EX_glc__D_e" in model.reactions:
        medium["EX_glc__D_e"] = 0
    if glucose:
        medium["EX_glc__D_e"] = 10
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
        set_lb_medium(model, glucose=False)
        lb_secretions = find_secretions(model)

    with model:
        set_lb_medium(model, glucose=True)
        glc_secretions = find_secretions(model)

    df = pd.DataFrame({
        "LB_max_secretion": lb_secretions,
        "LB_glucose_max_secretion": glc_secretions
    }).fillna(0)

    df["difference"] = df["LB_glucose_max_secretion"] - df["LB_max_secretion"]

    df = df.sort_values("difference", ascending=False)

    names = {m.id: m.name for m in model.metabolites}
    df["metabolite_name"] = df.index.map(names)

    df = df[[
        "metabolite_name",
        "LB_max_secretion",
        "LB_glucose_max_secretion",
        "difference"
    ]]

    print(df.head(30))

    df.to_csv("/work/lylab/cjn40747/metabolome/secreted_metabolites_FVA.csv")