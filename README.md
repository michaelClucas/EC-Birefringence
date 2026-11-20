# EC-Birefringence

An autonomous research program testing whether Einstein-Cartan torsion can explain cosmic birefringence, the rotation of the CMB polarization plane reported in ACT DR6 (arXiv:2509.13654) and constrained in arXiv:2510.07928, 2601.13624, and 2512.21888.

Run by a Claude Code overseer agent that dispatches teams of research subagents in a loop. See `OVERSEER_PROMPT.md` for the program and its stop conditions, `STATE.md` for current status, and `CLAIMS.md` for candidate mechanisms.

Possible answers are positive, negative (a no-go), or mixed. Nothing here is a conclusion unless `RESULT.md` exists.

## Layout
- `literature/`: paper notes and bounds with citations
- `data/`: measured inputs with sources
- `src/`: shared code
- `work/`: one folder per task, including dead ends
