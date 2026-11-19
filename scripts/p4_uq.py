"""Phase 4 (C2): tolerance uncertainty quantification with a verified surrogate.

Steps, in plain language
 1. Surrogate.  A Gaussian-process (GP) regression model learns log L(100 MHz) from the
    256 'train' FastHenry runs (11 inputs: 5 design + 6 tolerance variables).  A GP gives
    both a prediction and its own uncertainty (sigma).
 2. Calibrated error bars (split conformal prediction).  On 64 independent 'cal' runs we
    measure how wrong the GP is relative to its own sigma, s_i = |y_i - mu_i| / sigma_i,
    and take the ceil((n+1)(1-alpha))-th smallest s as the multiplier q.  The interval
    mu +- q sigma then contains the true solver value with probability >= 1-alpha for
    new inputs from the same distribution -- a guarantee that does not rely on the GP
    being right.  Checked on 64 further 'test' runs.
 3. Tolerance propagation.  At the nominal design, the six tolerance variables are
    sampled uniformly within their ranges (1e5 samples through the GP) and compared with
    40 direct FastHenry Monte Carlo runs ('mc_nominal').
 4. Sobol sensitivity indices (SALib, Saltelli sampling): the share of the variance of L
    caused by each tolerance alone (S1) and including interactions (ST), at the nominal
    (3313) and default (7628) stackups.
Fixed seeds.  Outputs in results/phase4/: gp_metrics.csv, conformal.csv, tolerance_mc.csv,
sobol.csv, gp_model.pkl, fig_p4_parity.(pdf|png), fig_p4_sobol.(pdf|png), fig_p4_mc.(pdf|png)
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

sys.path.insert(0, str(Path(__file__).parent))
from p4_campaign import HI, LO, NOMINAL, TOL, VARS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase4"


class Surrogate:
    """GP on inputs scaled to [0,1], output log(L/nH)."""

    def __init__(self, seed=0):
        k = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=np.ones(len(VARS)), length_scale_bounds=(1e-2, 1e3),
                                                     nu=2.5) + WhiteKernel(1e-6, (1e-10, 1e-2))
        self.gp = GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=8, random_state=seed)
        self.q = {}

    def _x(self, X):
        return (np.asarray(X, float) - LO) / (HI - LO)

    def fit(self, X, L):
        self.gp.fit(self._x(X), np.log(L))
        return self

    def predict(self, X):
        mu, sd = self.gp.predict(self._x(X), return_std=True)
        return mu, sd                                   # in log(L)

    def calibrate(self, X, L, alphas=(0.1, 0.05)):
        mu, sd = self.predict(X)
        s = np.abs(np.log(L) - mu) / sd
        n = len(s)
        for a in alphas:
            k = int(np.ceil((n + 1) * (1 - a)))
            self.q[a] = float(np.sort(s)[min(k, n) - 1])
        return s

    def interval(self, X, alpha):
        mu, sd = self.predict(X)
        return np.exp(mu - self.q[alpha] * sd), np.exp(mu), np.exp(mu + self.q[alpha] * sd)


def load():
    C = pd.read_csv(OUT / "campaign.csv")
    C = C[C.err.fillna("") == ""]
    return C


def sample_tol(n, rng):
    u = rng.random((n, len(TOL)))
    return LO[5:] + u * (HI[5:] - LO[5:])


def X_at(design: dict, tol_samples):
    d = np.array([design[v] for v in VARS[:5]], float)
    return np.hstack([np.tile(d, (len(tol_samples), 1)), tol_samples])


def main():
    C = load()
    tr, ca, te, mc = (C[C.set == s] for s in ("train", "cal", "test", "mc_nominal"))
    sur = Surrogate().fit(tr[VARS].values, tr.L_nH.values)
    s_cal = sur.calibrate(ca[VARS].values, ca.L_nH.values)
    rows = []
    for name, d in (("train", tr), ("cal", ca), ("test", te), ("mc_nominal", mc)):
        mu, sd = sur.predict(d[VARS].values)
        e = 100 * (np.exp(mu) / d.L_nH.values - 1)
        rows.append(dict(set=name, n=len(d), rmse_pct=float(np.sqrt(np.mean(e ** 2))), max_abs_pct=float(np.abs(e).max()),
                         mean_sigma_pct=float(100 * np.mean(sd))))
    M = pd.DataFrame(rows)
    M["kernel"] = str(sur.gp.kernel_)
    M.to_csv(OUT / "gp_metrics.csv", index=False)
    print(M.drop(columns="kernel").round(3).to_string(), "\n", sur.gp.kernel_, flush=True)
    crow = []
    for a in (0.1, 0.05):
        for name, d in (("test", te), ("mc_nominal", mc)):
            lo, mid, hi = sur.interval(d[VARS].values, a)
            cov = np.mean((d.L_nH.values >= lo) & (d.L_nH.values <= hi))
            crow.append(dict(alpha=a, target_coverage=1 - a, q=sur.q[a], set=name, n=len(d), coverage=float(cov),
                             mean_halfwidth_pct=float(100 * np.mean((hi - lo) / 2 / mid))))
    Cf = pd.DataFrame(crow)
    Cf.to_csv(OUT / "conformal.csv", index=False)
    print(Cf.round(4).to_string(), flush=True)
    with open(OUT / "gp_model.pkl", "wb") as fh:
        pickle.dump(sur, fh)

    # ---- tolerance Monte Carlo at the nominal design: surrogate vs direct FastHenry
    rng = np.random.default_rng(21)
    Xn = X_at(NOMINAL, sample_tol(100_000, rng))
    mu, _ = sur.predict(Xn)
    Ls = np.exp(mu)
    trow = [dict(method="GP surrogate (1e5)", n=len(Ls), mean_nH=Ls.mean(), std_nH=Ls.std(),
                 p2_5_nH=np.percentile(Ls, 2.5), p97_5_nH=np.percentile(Ls, 97.5)),
            dict(method="FastHenry direct MC", n=len(mc), mean_nH=mc.L_nH.mean(), std_nH=mc.L_nH.std(ddof=1),
                 p2_5_nH=np.percentile(mc.L_nH, 2.5), p97_5_nH=np.percentile(mc.L_nH, 97.5))]
    for d1n, name in ((0.076, "1080"), (0.210, "7628")):
        dz = dict(NOMINAL, d1=d1n)
        m2, _ = sur.predict(X_at(dz, sample_tol(100_000, rng)))
        L2 = np.exp(m2)
        trow.append(dict(method=f"GP surrogate, {name} stackup", n=len(L2), mean_nH=L2.mean(), std_nH=L2.std(),
                         p2_5_nH=np.percentile(L2, 2.5), p97_5_nH=np.percentile(L2, 97.5)))
    T = pd.DataFrame(trow)
    T["rel_halfwidth95_pct"] = 100 * (T.p97_5_nH - T.p2_5_nH) / 2 / T.mean_nH
    T.to_csv(OUT / "tolerance_mc.csv", index=False)
    print(T.round(4).to_string(), flush=True)

    # ---- Sobol indices over the tolerance variables
    from SALib.analyze import sobol
    from SALib.sample import sobol as sobol_sample
    prob = dict(num_vars=len(TOL), names=TOL, bounds=[[a, b] for a, b in zip(LO[5:], HI[5:])])
    srows = []
    for d1n, name in ((0.099, "3313 (nominal)"), (0.210, "7628 (default)")):
        Xs = sobol_sample.sample(prob, 4096, calc_second_order=False, seed=3)
        m, _ = sur.predict(X_at(dict(NOMINAL, d1=d1n), Xs))
        Si = sobol.analyze(prob, np.exp(m), calc_second_order=False, seed=3)
        for k, v in enumerate(TOL):
            srows.append(dict(stackup=name, variable=v, S1=Si["S1"][k], S1_conf=Si["S1_conf"][k], ST=Si["ST"][k],
                              ST_conf=Si["ST_conf"][k]))
    S = pd.DataFrame(srows)
    S.to_csv(OUT / "sobol.csv", index=False)
    print(S.round(4).to_string(), flush=True)
    figures(sur, te, mc, Ls, S)


def figures(sur, te, mc, Ls, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lo, mid, hi = sur.interval(te[VARS].values, 0.05)
    fig, ax = plt.subplots(figsize=(3.6, 3.2), constrained_layout=True)
    ax.errorbar(te.L_nH, mid, yerr=[mid - lo, hi - mid], fmt="o", ms=3, lw=0.8, capsize=0, label="test set, 95 % conformal")
    a = [te.L_nH.min() * 0.95, te.L_nH.max() * 1.05]
    ax.plot(a, a, "k--", lw=0.8)
    ax.set_xlabel("FastHenry $L$ (nH)")
    ax.set_ylabel("surrogate $L$ (nH)")
    ax.legend(fontsize=7)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p4_parity.{ext}", dpi=200)
    fig, ax = plt.subplots(figsize=(4.4, 3.0), constrained_layout=True)
    labels = {"t_d1": "prepreg $d_1$", "t_cu": "outer Cu", "t_cuin": "inner Cu", "hole": "via hole",
              "bump": "bump height", "dx_die": "die offset"}
    for k, (st, g) in enumerate(S.groupby("stackup", sort=False)):
        x = np.arange(len(g)) + (k - 0.5) * 0.38
        ax.bar(x, g.ST, 0.36, yerr=g.ST_conf, label=st, capsize=2)
    ax.set_xticks(range(len(TOL)), [labels[v] for v in TOL], rotation=30, fontsize=7)
    ax.set_ylabel("total Sobol index $S_T$")
    ax.legend(fontsize=7)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p4_sobol.{ext}", dpi=200)
    fig, ax = plt.subplots(figsize=(3.6, 2.8), constrained_layout=True)
    ax.hist(Ls, bins=60, density=True, alpha=0.6, label="GP, $10^5$ samples")
    ax.hist(mc.L_nH, bins=12, density=True, histtype="step", lw=1.5, color="k", label="FastHenry, 40 runs")
    ax.set_xlabel("$L$ at 100 MHz (nH)")
    ax.set_ylabel("density")
    ax.legend(fontsize=7)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p4_mc.{ext}", dpi=200)


if __name__ == "__main__":
    main()
