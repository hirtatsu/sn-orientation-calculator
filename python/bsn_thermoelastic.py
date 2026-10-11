"""Thermoelastic response of a beta-Sn grain under in-plane constraint.

Reference implementation of the calculation in index.html.

Model
-----
Single beta-Sn crystal. In-plane total strain is fixed
(eps_xx = eps_yy = eps_xy = 0) and the surface normal ND (= z) is
traction-free (sig_zz = sig_xz = sig_yz = 0). The response is linear in
the temperature change dT; results are reported per kelvin.

Material data: Xian et al., Nat. Commun. 15 (2024),
doi:10.1038/s41467-024-48532-6, Supplementary Table 4.

Usage
-----
    python bsn_thermoelastic.py PHI1 PHI PHI2      # Bunge Euler angles in degrees
"""
import sys
import numpy as np

LAT_A, LAT_C = 5.831, 3.182  # lattice parameters (Angstrom)

# family, plane (hkl), direction [uvw], CRSS ratio
SLIP = [
    ("{100}<001>", (-1, 0, 0), (0, 0, -1), 1.0), ("{100}<001>", (0, -1, 0), (0, 0, -1), 1.0),
    ("{110}<001>", (-1, -1, 0), (0, 0, -1), 1.0), ("{110}<001>", (-1, 1, 0), (0, 0, -1), 1.0),
    ("{100}<010>", (-1, 0, 0), (0, -1, 0), 1.05), ("{100}<010>", (0, -1, 0), (-1, 0, 0), 1.05),
    ("{110}<1-11>/2", (-1, -1, 0), (-1, 1, -1), 1.1), ("{110}<1-11>/2", (-1, -1, 0), (-1, 1, 1), 1.1),
    ("{110}<1-11>/2", (-1, 1, 0), (-1, -1, -1), 1.1), ("{110}<1-11>/2", (-1, 1, 0), (-1, -1, 1), 1.1),
    ("{110}<1-10>", (-1, -1, 0), (-1, 1, 0), 1.2), ("{110}<1-10>", (-1, 1, 0), (-1, -1, 0), 1.2),
    ("{100}<011>", (-1, 0, 0), (0, -1, -1), 1.25), ("{100}<011>", (-1, 0, 0), (0, -1, 1), 1.25),
    ("{100}<011>", (0, -1, 0), (-1, 0, -1), 1.25), ("{100}<011>", (0, -1, 0), (-1, 0, 1), 1.25),
    ("{001}<010>", (0, 0, -1), (0, -1, 0), 1.3), ("{001}<010>", (0, 0, -1), (-1, 0, 0), 1.3),
    ("{001}<110>", (0, 0, -1), (-1, -1, 0), 1.4), ("{001}<110>", (0, 0, -1), (-1, 1, 0), 1.4),
    ("{011}<01-1>", (0, -1, -1), (0, -1, 1), 1.5), ("{011}<01-1>", (0, -1, 1), (0, -1, -1), 1.5),
    ("{011}<01-1>", (-1, 0, -1), (-1, 0, 1), 1.5), ("{011}<01-1>", (-1, 0, 1), (-1, 0, -1), 1.5),
    ("{211}<01-1>", (-2, -1, -1), (0, -1, 1), 1.5), ("{211}<01-1>", (-2, -1, 1), (0, -1, -1), 1.5),
    ("{211}<01-1>", (-2, 1, -1), (0, 1, 1), 1.5), ("{211}<01-1>", (-2, 1, 1), (0, 1, -1), 1.5),
    ("{211}<01-1>", (-1, -2, -1), (-1, 0, 1), 1.5), ("{211}<01-1>", (-1, -2, 1), (-1, 0, -1), 1.5),
    ("{211}<01-1>", (-1, 2, -1), (-1, 0, 1), 1.5), ("{211}<01-1>", (-1, 2, 1), (-1, 0, -1), 1.5),
]


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


SLIP_N = np.array([_unit((h / LAT_A, k / LAT_A, l / LAT_C)) for _, (h, k, l), _, _ in SLIP])
SLIP_S = np.array([_unit((u * LAT_A, v * LAT_A, w * LAT_C)) for _, _, (u, v, w), _ in SLIP])
SLIP_R = np.array([r for *_, r in SLIP])


def props(T=25.0):
    """Stiffness (MPa, Voigt 11,22,33,23,13,12, crystal frame) and CTEs (1/K) at T (deg C)."""
    E1, E3 = (25.3 - 0.083 * T) * 1e3, (68.3 - 0.072 * T) * 1e3
    G12, G13 = (24.4 - 0.020 * T) * 1e3, (22.5 - 0.026 * T) * 1e3
    n12, n13 = 0.46, 0.22
    S = np.zeros((6, 6))
    S[0, 0] = S[1, 1] = 1 / E1
    S[2, 2] = 1 / E3
    S[0, 1] = S[1, 0] = -n12 / E1
    S[0, 2] = S[2, 0] = S[1, 2] = S[2, 1] = -n13 / E1
    S[3, 3] = S[4, 4] = 1 / G13
    S[5, 5] = 1 / G12
    aa = 14.195e-6 + 25.53e-9 * T
    ac = 28.662e-6 + 56.79e-9 * T
    return np.linalg.inv(S), aa, ac


def bond(a):
    """6x6 Bond matrix for stiffness: C' = K C K^T for the rotation v' = a v."""
    P = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    K = np.zeros((6, 6))
    for I, (i, j) in enumerate(P):
        for J, (k, l) in enumerate(P):
            if J < 3:
                K[I, J] = a[i, k] * a[j, l]
            else:
                K[I, J] = a[i, k] * a[j, l] + a[i, l] * a[j, k]
    return K


def bunge_g(phi1, Phi, phi2):
    """Passive Bunge orientation matrix (crystal = g . sample), angles in degrees."""
    p1, P, p2 = np.radians([phi1, Phi, phi2])
    c1, s1, c, s, c2, s2 = np.cos(p1), np.sin(p1), np.cos(P), np.sin(P), np.cos(p2), np.sin(p2)
    return np.array([
        [c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
        [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
        [s1 * s, -c1 * s, c],
    ])


def solve(M, T=25.0, dT=-1.0):
    """M: rows = sample axes x, y, ND in crystal coordinates (v_sample = M v_crystal).

    Returns stress tensor in the sample frame (MPa), elastic strain energy density
    per dT^2 (J m^-3 K^-2), K_tau = max |tau|/r_CRSS per |dT| (MPa K^-1), and the
    resolved shear stress of every slip system (MPa).
    """
    C, aa, ac = props(T)
    K = bond(M)
    Cs = K @ C @ K.T
    As = M @ np.diag([aa, aa, ac]) @ M.T
    th = dT * np.array([As[0, 0], As[1, 1], As[2, 2], 2 * As[1, 2], 2 * As[0, 2], 2 * As[0, 1]])
    f = [2, 3, 4]                       # free components: eps_zz, gamma_yz, gamma_xz
    e = np.zeros(6)
    e[f] = np.linalg.solve(Cs[np.ix_(f, f)], Cs[f] @ th)
    el = e - th
    sv = Cs @ el
    sig = np.array([[sv[0], sv[5], sv[4]], [sv[5], sv[1], sv[3]], [sv[4], sv[3], sv[2]]])
    W = 0.5 * sv @ el * 1e6 / dT ** 2
    s = SLIP_S @ M.T
    n = SLIP_N @ M.T
    tau = np.einsum("ia,ab,ib->i", s, sig, n)
    K_tau = np.max(np.abs(tau) / SLIP_R) / abs(dT)
    return sig, W, K_tau, tau


def ipf_angles(nd):
    """theta_001, theta_100 and the azimuth phi (deg) after reduction to the 4/mmm standard triangle."""
    x, y, z = np.abs(nd)
    if y > x:
        x, y = y, x
    return (np.degrees(np.arccos(min(1.0, z))), np.degrees(np.arccos(min(1.0, x))),
            np.degrees(np.arctan2(y, x)))


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 1
    phi1, Phi, phi2 = map(float, argv)
    g = bunge_g(phi1, Phi, phi2)
    M = g.T
    nd = g[:, 2]
    sig, W, K_tau, tau = solve(M)
    t001, t100, phi = ipf_angles(nd)
    print(f"ND (crystal frame): {np.round(nd, 4)}")
    print(f"theta_001 = {t001:.2f} deg, theta_100 = {t100:.2f} deg, phi = {phi:.2f} deg")
    print(f"K_tau = {K_tau:.4f} MPa/K   W = {W:.3f} J m^-3 K^-2")
    print("stress per kelvin of cooling (MPa/K, sample frame):")
    print(np.array2string(sig, precision=4, suppress_small=True))
    order = np.argsort(-np.abs(tau) / SLIP_R)[:6]
    print("top slip systems:")
    for i in order:
        fam, hkl, uvw, r = SLIP[i]
        print(f"  {fam:15s} {str(hkl):13s} {str(uvw):13s} r={r:.2f}  tau={tau[i]: .4f}  |tau|/r={abs(tau[i]) / r:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
