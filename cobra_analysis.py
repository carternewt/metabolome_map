import cobra 
import pandas as pd

model_path = "work/lylab/cjn40747/metabolome/F7_5.tsv"
model = cobra.io.read_sbml_model(model_path)

def set_lb_medium(model, glucose=False):
    medium = model.medium.copy()
    if "EX_glc__D_e" in model.reactions:
        if glucose:
            medium["EX_glc__D_e"] = 10
        else:
            medium["EX_glc__D_e"] = 0
    model.medium = medium

def get_secreted_metabolites(model):
    solution = model.optimize()
    secretions = {}
    for rxn in model.exchanges:
        flux = solution.fluxes[rxn.id]
        if flux > 1e-6:
            met = list(rxn.metabolites.keys())[0]
            secretions[met.id] = flux
    return secretions

with model:
    set_lb_medium(model, glucose=False)
    lb_secretions = get_secreted_metabolites(model)

with model:
    set_lb_medium(model, glucose=True)
    glc_secretions = get_secreted_metabolites(model)

df = pd.DataFrame({
    "LB_flux": lb_secretions,
    "LB_glucose_flux": glc_secretions
}).fillna(0)

df["difference"] = df["LB_glucose_flux"] - df["LB_flux"]

df = df.sort_values("difference", ascending=False)

print(df.head(50))

df.to_csv("/work/lylab/cjn40747/metabolome/secreted_metabolites_comparison.csv")