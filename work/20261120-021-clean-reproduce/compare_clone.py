"""Compare the clean clone outputs (artifacts/clean_checkout/artifacts/) with the main checkout and the committed
014 fit_results.csv. Run from the repo root after the clean-clone reproduce."""
import numpy as np, json, csv
a=np.load("artifacts/templates/spt_bb_templates_v1.npz"); b=np.load("artifacts/clean_checkout/artifacts/templates/spt_bb_templates_v1.npz")
print("templates max rel diff", max(float(np.max(np.abs(a[k]-b[k])/(np.abs(a[k])+1e-300))) for k in a.files))
print("cb_scale", json.load(open("artifacts/templates/cb_scale.json"))["scale"], json.load(open("artifacts/clean_checkout/artifacts/templates/cb_scale.json"))["scale"])
A=list(csv.DictReader(open("work/20261120-014-data-reproduce-fit/fit_results.csv"))); B=list(csv.DictReader(open("artifacts/clean_checkout/artifacts/reproduce/fit_results.csv")))
print(len(A), len(B))
for x,y in zip(A,B):
    d=[(k,x[k],y[k]) for k in x if x[k]!=y[k]]
    if d: print(x["dataset"],x["model"],d)
