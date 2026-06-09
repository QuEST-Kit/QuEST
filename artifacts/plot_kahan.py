#!/usr/bin/env python3
"""
Visualise the cost/benefit of compensated (Kahan) vs naive accumulation in
cpu_statevec_anyCtrlAnyTargDenseMatr_sub(), per qcomp precision.

Input : results.csv (produced by bench_kahan)
Output: kahan_accuracy.png, kahan_runtime.png, kahan_costbenefit.png
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
with open("results.csv") as f:
    for r in csv.DictReader(f):
        rows.append(r)

precs = ["fp1", "fp2", "fp4"]
labels = {
    "fp1": "fp1  (float, 24-bit)",
    "fp2": "fp2  (double, 53-bit)",
    "fp4": "fp4* (long double = double here, 53-bit)",
}
colors = {"fp1": "#d62728", "fp2": "#1f77b4", "fp4": "#2ca02c"}

def series(p, key):
    xs, ys = [], []
    for r in rows:
        if r["precision"] == p:
            xs.append(int(r["numTargets"]))
            ys.append(float(r[key]))
    return xs, ys

# ---------------- accuracy ----------------
fig, ax = plt.subplots(figsize=(8, 5.5))
for p in precs:
    x, yn = series(p, "err_naive")
    _, yk = series(p, "err_kahan")
    ax.plot(x, yn, "--o", color=colors[p], label=f"{labels[p]} - naive", alpha=0.9)
    ax.plot(x, yk, "-s", color=colors[p], label=f"{labels[p]} - Kahan", alpha=0.9)
ax.set_yscale("log")
ax.set_xlabel("number of target qubits  (matrix is 2^n x 2^n)")
ax.set_ylabel("worst-case abs error vs __float128 (113-bit) reference")
ax.set_title("Accuracy: Kahan vs naive dense-matrix accumulation\n"
             "(adversarial ill-conditioned matrix, single-CPU)")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig("kahan_accuracy.png", dpi=130)

# ---------------- runtime ----------------
fig, ax = plt.subplots(figsize=(8, 5.5))
for p in precs:
    x, tn = series(p, "ms_naive")
    _, tk = series(p, "ms_kahan")
    ax.plot(x, tn, "--o", color=colors[p], label=f"{labels[p]} - naive")
    ax.plot(x, tk, "-s", color=colors[p], label=f"{labels[p]} - Kahan")
ax.set_yscale("log")
ax.set_xlabel("number of target qubits  (matrix is 2^n x 2^n)")
ax.set_ylabel("runtime per applyCompMatr  [ms]  (single CPU, -O3)")
ax.set_title("Runtime cost: Kahan vs naive dense-matrix accumulation")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig("kahan_runtime.png", dpi=130)

# ---------------- combined cost/benefit ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.5))
for p in precs:
    x, yn = series(p, "err_naive")
    _, yk = series(p, "err_kahan")
    # accuracy benefit = error reduction factor (naive/kahan)
    benefit = [(n / k) if k > 0 else 1.0 for n, k in zip(yn, yk)]
    a1.plot(x, benefit, "-o", color=colors[p], label=labels[p])
a1.axhline(1.0, color="k", lw=0.8, ls=":")
a1.set_yscale("log")
a1.set_xlabel("number of target qubits")
a1.set_ylabel("accuracy BENEFIT  =  err_naive / err_kahan  (>1 = Kahan better)")
a1.set_title("Benefit: error-reduction factor from Kahan")
a1.grid(True, which="both", alpha=0.3)
a1.legend(fontsize=8)

for p in precs:
    x, tn = series(p, "ms_naive")
    _, tk = series(p, "ms_kahan")
    cost = [k / n if n > 0 else 1.0 for n, k in zip(tn, tk)]
    a2.plot(x, cost, "-s", color=colors[p], label=labels[p])
a2.axhline(1.0, color="k", lw=0.8, ls=":")
a2.set_xlabel("number of target qubits")
a2.set_ylabel("runtime COST  =  ms_kahan / ms_naive  (>1 = Kahan slower)")
a2.set_title("Cost: runtime slowdown from Kahan")
a2.grid(True, alpha=0.3)
a2.legend(fontsize=8)
fig.suptitle("Cost / Benefit of Kahan summation in cpu_statevec_anyCtrlAnyTargDenseMatr_sub()",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("kahan_costbenefit.png", dpi=130)

print("wrote kahan_accuracy.png, kahan_runtime.png, kahan_costbenefit.png")
