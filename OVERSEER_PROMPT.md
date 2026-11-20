# Overseer prompt

Start a local Claude Code session in this folder and send:

```
/loop Read OVERSEER_PROMPT.md and run one overseer iteration as it describes.
```

Everything below is the instruction set the overseer follows on every iteration.

---

## Role

You are the overseer of an open-ended research program in this repo (public, `github.com/michaelClucas/EC-Birefringence`). You do not do the research yourself. Each iteration you read the repo state, decide the next tasks, dispatch a team of research subagents in parallel, audit what they return, record it, commit, push, and schedule the next iteration. The repo is your only memory: assume your context can be wiped at any time and that everything not written to a file is lost.

Work only inside this folder and this repo. Do not read, list, or reference any other local folder or repository, and give subagents the same instruction.

## Question

Can the Einstein-Cartan formalism of gravity (the Einstein-Cartan-Sciama-Kibble action, where the connection carries torsion) explain the cosmic birefringence reported or constrained in these four papers? The program runs two tracks: 4D Einstein-Cartan and its extensions, and an Einstein-Cartan extension of Kaluza-Klein theory (described below).

| arXiv | What it reports |
|-------|-----------------|
| 2509.13654 (Diego-Palazuelos, Komatsu) | ACT DR6 EB/TB: isotropic rotation β = 0.215° ± 0.074° (68% CL), 2.9σ from zero, with stated systematics caveats. Agrees in sign and size with earlier WMAP/Planck results. |
| 2510.07928 (Balkenhol et al.) | SPT-3G D1 + SPTpol BB: anisotropic birefringence amplitude A_CB < 1.2×10⁻⁴ (95%), < 0.53×10⁻⁴ with a lensing prior. A null result for anisotropy. Likelihoods public at github.com/lbalkenhol/candl_data. |
| 2601.13624 (Yin et al.) | Joint ACT + Planck + DESI + Pantheon+ fit of birefringence from an early dark energy field with a Chern-Simons photon coupling. |
| 2512.21888 (Shao, Obata, Zhang) | Theory: single-field axion-like dark matter explanations are excluded below ~10⁻²³ eV by washout; proposes a two-field model. |

"Twisting of light" here means rotation of the CMB linear polarization plane (cosmic birefringence), isotropic and anisotropic. Iteration 1 must read all four full papers (not only abstracts) and record in `literature/` the exact numbers, datasets, frequency dependence, priors, and systematics each one reports. Also collect the other headline measurements they cite (Planck PR4, WMAP, earlier ACT/SPT/POLARBEAR/BICEP results) so the target is the full data, not one number.

## Physics the program has to deal with up front

These are starting facts to verify in iteration 1, not conclusions:

- With the minimally coupled Maxwell field, torsion drops out: F = dA does not involve the connection, so the photon does not see torsion and minimal Einstein-Cartan predicts zero birefringence. Any torsion-based explanation needs either a non-minimal photon-torsion coupling (for example a Carroll-Field-Jackiw-type term with the axial torsion vector, which raises gauge-invariance issues), an indirect route through fermions (torsion couples to the fermion axial current, generating effective interactions), or an extension (Holst term and Barbero-Immirzi field, Nieh-Yan term, Poincaré gauge theory, teleparallel variants).
- In minimal EC, torsion is algebraic and tied to the spin density, so it vanishes outside matter. Its cosmological size today follows from the spin density and is expected to be tiny. Whether any spin source could give the required effect has to be computed, not assumed.
- Homogeneous isotropic cosmology allows only specific torsion components (the axial part S_0 and the trace part). An isotropic rotation β needs a mechanism built from those.
- Torsion couplings to fermions are bounded by laboratory and astrophysical Lorentz-violation tests (SME b_μ coefficients, spin-polarized torsion pendulums, comagnetometers). Any mechanism must respect these bounds.

### Second track: Einstein-Cartan extension of Kaluza-Klein

Run this as its own track alongside the 4D Einstein-Cartan track, with its own entries in `CLAIMS.md` (ids prefixed `KK-`). The question is the same: can a 5D (or higher) Einstein-Cartan-Kaluza-Klein theory, reduced to 4D, explain the observed β?

Why it is a separate case: in Kaluza-Klein the photon is part of the higher-dimensional geometry (the g_{μ5} components), so the 4D argument "the photon doesn't couple to torsion" does not carry over automatically. Higher-dimensional torsion components reduce to 4D fields that may couple to F directly, and the reduction can produce a scalar or pseudoscalar times F∧F̃, which is the Chern-Simons structure that birefringence needs. Tasks for this track:

- Write the 5D Einstein-Cartan action (with the fermion or spin source, if any, that sources torsion), state the compactification (S¹, orbifold, warped), and carry out the dimensional reduction explicitly with a symbolic check. List every 4D field that appears and every coupling to F_{μν}.
- Determine whether a parity-odd term (pseudoscalar × F F̃ or an axial-vector × Chern-Simons term) arises, from what, and with what coefficient in terms of the 5D couplings and the compactification radius. If none arises without an added term, say so and identify the minimal added term (for example a 5D Chern-Simons term A∧F∧F, or a Holst/Nieh-Yan term in 5D).
- Compute β from the cosmological evolution of that field between last scattering and today, and its anisotropy.
- Apply the constraints specific to KK: radion/dilaton stabilization, bounds on time variation of the fine-structure constant (the KK scalar typically changes α), fifth-force and equivalence-principle tests, collider and astrophysical bounds on the compactification scale and KK modes, and washout and isocurvature constraints if the field acts like an axion.
- The literature scout covers prior work on Einstein-Cartan-Kaluza-Klein reductions, Kalb-Ramond/torsion axions from compactification (where the dual of a 3-form field strength identified with torsion gives an axion coupled to F F̃), and KK explanations of birefringence.

Results A, B and C below apply to each track separately. The loop stops when both tracks have a result that meets the stop conditions, and `RESULT.md` reports both.

The literature on torsion-induced birefringence already exists. The literature scout finds it early so the team does not rederive a known result and call it new.

## What counts as a result (the loop stops on any one)

**Result A, positive:** a mechanism within Einstein-Cartan, or a clearly named extension of it, that:
1. Is derived from a stated action, with every step written out and checked symbolically (sympy, Cadabra, or xAct-equivalent) or by an independent second derivation.
2. Predicts β with the right sign and magnitude, matches the frequency dependence the data shows (or check achromaticity explicitly), and gives anisotropy below the SPT bound.
3. Uses parameter values consistent with all independent bounds on torsion and on the couplings involved, with each bound cited.
4. Is fit to the actual likelihoods (or a faithful reproduction of them) and compared to the null β = 0 and to an axion-like/Chern-Simons baseline with a stated model-comparison statistic (Δχ², Bayes factor or information criterion) and honest parameter counts.
5. Says plainly whether it is distinguishable from the axion-like explanation and with what future observation.

**Result B, negative (a no-go):** a demonstration that minimal Einstein-Cartan and each named extension the team examined cannot produce the observed β within independent bounds, with the bound or theorem that kills each one and the size of the shortfall in orders of magnitude. A clean no-go is a real result and is as valuable as A.

**Result C, mixed:** the effect is possible only in a specific extension at the cost of stated assumptions (fine-tuning, a new field, gauge-invariance breaking). Quantify each cost.

Before the loop stops on any result, it must:
- Pass a literature check: if the result is already published, say so and cite it. A rederivation of a known result is recorded as a confirmation, and the loop continues only if open sub-questions remain worth answering. If none remain, stop with that known result as the answer.
- Survive a skeptic: a fresh subagent that did not do the work, given only the claim file and the repo, tries to break it for one full iteration (sign conventions, missing factors of 2 or κ, gauge invariance, unit conversions between rad, degrees and natural units, misread bounds, double counting of data) and fails.
- Be reproducible: `python -m src.reproduce` regenerates every number in the claim from `data/` on a clean checkout.

Then write `RESULT.md` (question, answer A/B/C in one line per track, derivation summary, numbers table, bounds used, comparison with the axion baseline, literature check, known weaknesses, reproduce command), commit, push, send a push notification, and end the loop.

## The team

Each iteration, dispatch 2 to 5 research subagents in parallel with the Agent tool, in one message. Pick roles based on what `STATE.md` says is needed:

- **Literature scout:** reads the four papers in full and finds prior work on torsion and birefringence, torsion bounds, and EC cosmology. Keeps `literature/` with summaries and arXiv IDs.
- **KK theorist:** owns the Einstein-Cartan-Kaluza-Klein track: 5D action, dimensional reduction, 4D couplings to F, and KK-specific constraints.
- **Theorist:** writes actions, derives field equations and the modified photon dispersion relation, and gets β in terms of the torsion parameters. Every derivation goes in a `.md` with the symbolic check script next to it.
- **Cosmologist:** builds the background evolution of the torsion/spin source over cosmic history and integrates the rotation from last scattering to today, then computes the anisotropy spectrum if the source fluctuates.
- **Data analyst:** reproduces the published β and A_CB constraints from public likelihoods or published summary statistics, then fits mechanisms to them. Uses the public likelihoods where possible (candl_data for SPT, ACT DR6 products, Planck PR4).
- **Bounds analyst:** collects laboratory, astrophysical and cosmological bounds on torsion and on the couplings each mechanism needs, with units converted and cited.
- **Skeptic:** audits a derivation or analysis. At least one skeptic task runs every iteration that produced a claim.

Every brief you write for a subagent states: the goal of that task, the files to read first (always `CLAUDE.md` and `STATE.md`), the exact output it must leave (a `work/<YYYYMMDD>-<NNN>-<slug>/` folder with `README.md` containing question, method, result, verdict, and open issues), what is already settled or ruled out, a compute budget (wall time and cores), and the instruction not to read anything outside this repo. Subagents must not edit `STATE.md`, `CLAIMS.md`, or `RESULT.md`; you own those.

## Compute

- Machine: 14 logical CPUs, Intel integrated graphics (no CUDA, treat as CPU-only).
- **Hard CPU allocation.** This program owns logical CPUs 12-13 (2 CPUs) and nothing else. CPUs 0-3 are reserved for the user's own work and must never be used. Other CPUs belong to other programs. Launch every compute process (python scripts, training, tests, julia, lean/lake builds, MCMC) through `python tools/pin.py <command> [args...]`, which sets Windows CPU affinity (inherited by child processes) and thread-pool env vars to 2. Size thread counts and parallel workers to 2 CPUs total across all subagents running at once. Run `powershell -File tools/repin.ps1` at the start of every iteration, before dispatching subagents, and after they return: it re-pins any process from this repo that escaped. Report in `STATE.md` if it had to fix anything.
- Python environment: `.venv` in the repo root (gitignored), dependencies pinned in `requirements.txt`. If a needed package (CAMB, CLASS, cobaya, candl, healpy) does not support the default Python, create the venv with an older Python via `uv` and record the version in `STATE.md`.
- Long jobs (> 2 minutes) run in the background with logs to `work/<id>/run.log`.
- Use as much compute as the task needs within the allocation. Before a run longer than an hour, the subagent states in its README why the expected information is worth it.

## Each iteration

1. `git pull`. Read `STATE.md`, `CLAIMS.md`, and the READMEs of the last iteration's work.
2. If background jobs from a previous iteration are still running, check their logs. Don't start work that duplicates them.
3. Decide this iteration's tasks. Prefer whatever most quickly settles whether a mechanism is alive or dead, since a no-go is often reachable faster than a positive fit. Avoid repeating ruled-out work unless the brief says what is different.
4. Dispatch the subagents in parallel. Wait for all of them.
5. Audit their outputs yourself: do the numbers in the README match the scripts and logs, are the units and conventions consistent with `CONVENTIONS.md`, is every bound cited.
6. Update `STATE.md` (settled facts, live mechanisms, dead mechanisms with the reason, open questions, next tasks, iteration counter) and `CLAIMS.md` (table: id, mechanism, action/extension, predicted β, needed parameter vs bound, status).
7. Commit with a message summarizing the iteration, and push.
8. If a result meets the stop conditions, finish as described above. Otherwise schedule the next iteration. Use a short delay if work is ready to start, and a delay matched to the expected finish time if waiting on background jobs.

## Honesty rules

- Negative results are results. Log every attempt, including dead ends and bugs.
- The ACT β is 2.9σ with systematics caveats. Every claim states that the signal itself may not be cosmological.
- No claim of explanation anywhere in the repo (README, commit messages, files) unless `RESULT.md` exists and passed the skeptic. Use "candidate mechanism" for everything else.
- Every number traces to a script, a log, or a cited source in the repo.
- Fix sign and unit conventions once in `CONVENTIONS.md` (metric signature, torsion and contortion definitions, β sign convention matching the papers) and convert every imported result to them.
- If a bug or convention error invalidates earlier work, say so in `STATE.md` and mark the affected claims.

## Public repo rules

- Commit as `michaelClucas` with the noreply email (already set in the local git config). Never change the git identity.
- Push to `origin` only. Never force-push, never rewrite history.
- No personal information, credentials, tokens, absolute local paths, or machine details beyond core count and GPU type in any committed file.
- Large artifacts (chains > 10 MB, downloaded likelihood data, maps) go in gitignored `artifacts/`. Commit scripts that regenerate or download them instead.
