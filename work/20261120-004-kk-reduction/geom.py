"""Numeric Riemann-Cartan curvature at a point, CONVENTIONS.md index order.

nabla_mu V^lam = d_mu V^lam + Gam[lam, mu, nu] V^nu   (derivative index = 2nd slot)
R^rho_{sig mu nu} = d_mu Gam^rho_{nu sig} - d_nu Gam^rho_{mu sig}
                    + Gam^rho_{mu lam} Gam^lam_{nu sig} - Gam^rho_{nu lam} Gam^lam_{mu sig}
Ricci R_{sig nu} = R^mu_{sig mu nu},  R = g^{sig nu} R_{sig nu}.

Inputs are arrays at one point: g[i,j], dg[a,i,j] = d_a g_ij, ddg[a,b,i,j].
"""
import numpy as np


def christoffel(g, dg, ddg):
    gi = np.linalg.inv(g)
    dgi = -np.einsum("im,amn,nj->aij", gi, dg, gi)
    # L[s,m,n] = d_m g_sn + d_n g_sm - d_s g_mn
    L = np.einsum("msn->smn", dg) + np.einsum("nsm->smn", dg) - dg
    dL = (np.einsum("amsn->asmn", ddg) + np.einsum("ansm->asmn", ddg)
          - ddg)
    Gam = 0.5 * np.einsum("rs,smn->rmn", gi, L)
    dGam = 0.5 * (np.einsum("ars,smn->armn", dgi, L)
                  + np.einsum("rs,asmn->armn", gi, dL))
    return gi, Gam, dGam


def ricci_scalar(gi, Gam, dGam):
    """Gam[r,m,n] = Gam^r_{mn}; dGam[a,r,m,n] = d_a Gam^r_{mn}."""
    # R^r_{s m n}
    Rie = (np.einsum("mrns->rsmn", dGam) - np.einsum("nrms->rsmn", dGam)
           + np.einsum("rml,lns->rsmn", Gam, Gam)
           - np.einsum("rnl,lms->rsmn", Gam, Gam))
    Ric = np.einsum("msmn->sn", Rie)
    return np.einsum("sn,sn->", gi, Ric)
