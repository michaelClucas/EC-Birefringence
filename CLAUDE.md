# Repo rules (overseer and all subagents)

Research program: can Einstein-Cartan torsion explain the cosmic birefringence (CMB polarization rotation) reported in arXiv 2509.13654, 2510.07928, 2601.13624, 2512.21888? The full program, stop conditions, and team structure are in `OVERSEER_PROMPT.md`. Read it before doing anything.

- Work only inside this repo. Do not read, list, or reference any other local folder or repository.
- `STATE.md` is the current state of the program. Read it first. Only the overseer edits it, `CLAIMS.md`, and `RESULT.md`.
- Each task lives in `work/<YYYYMMDD>-<NNN>-<slug>/` with a `README.md` (question, method, result, verdict, open issues) and the scripts needed to rerun it.
- Every derivation gets a symbolic or independent second check. Use the conventions in `CONVENTIONS.md` and convert imported results to them.
- Minimal Einstein-Cartan with a minimally coupled photon gives zero birefringence. Any mechanism must name the coupling or extension that changes this.
- Every bound and measured number is cited (arXiv ID, table/equation) in `data/` or `literature/`. Don't hardcode measured numbers in scripts: load them from `data/`.
- Python lives in `.venv/`. Add dependencies to `requirements.txt` with pinned versions.
- 14 CPU cores shared with other work: this program uses at most 4 unless load is low. Intel iGPU, no CUDA. Long jobs run in the background and log to `work/<id>/run.log`.
- Repo is public. No personal info, credentials, or absolute local paths in committed files. Don't change the git identity, don't force-push.
- Use "candidate mechanism" for every finding unless `RESULT.md` exists.
