# β-Sn Orientation Calculator

A browser-based calculator for the thermoelastic response of a β-Sn grain whose in-plane strain is constrained, as in a thin solder layer bonded to stiff parts.

Enter Bunge Euler angles (φ₁, Φ, φ₂). The calculator returns:

- the position of the sample normal (ND) on the inverse pole figure, with θ[001], θ⟨100⟩ and the azimuth φ;
- the maximum relative resolved shear stress coefficient *K*τ = max |τ| / *r*CRSS (MPa K⁻¹) and the slip systems ranked by |τ| / *r*CRSS;
- the stress tensor in the sample frame (MPa K⁻¹);
- the elastic strain energy density (J m⁻³ K⁻²).

Maps of *K*τ and the strain energy over the standard triangle are shown for reference.

## Use

Open `index.html` in a web browser. No installation or server is required.

A Python implementation of the same calculation is in `python/bsn_thermoelastic.py` (requires NumPy):

```
python python/bsn_thermoelastic.py 30 55 20
```

## Model

- Single β-Sn crystal (tetragonal, 4/mmm).
- In-plane total strain fixed: ε<sub>xx</sub> = ε<sub>yy</sub> = ε<sub>xy</sub> = 0.
- Traction-free along ND: σ<sub>zz</sub> = σ<sub>xz</sub> = σ<sub>yz</sub> = 0.
- Uniform temperature change Δ*T*. All quantities are linear in Δ*T* and reported per kelvin; *K*τ and the energy (per Δ*T*²) do not depend on the sign of Δ*T*. Stresses are given for cooling and reverse sign on heating.

The free strain components (ε<sub>zz</sub>, γ<sub>yz</sub>, γ<sub>xz</sub>) are obtained from the traction-free condition, and the stress follows from σ = **C**(ε − αΔ*T*) with the stiffness and thermal expansion tensors rotated to the sample frame. The strain energy density is *W* = ½ σ : ε<sup>el</sup>.

### Material data

Elastic constants, thermal expansion coefficients, the 32 slip systems and their CRSS ratios are taken from Supplementary Table 4 of Xian et al. (2024) and evaluated at 25 °C. Lattice parameters: *a* = 5.831 Å, *c* = 3.182 Å.

### Orientation convention

Bunge Euler angles, passive convention (crystal = **g** · sample), as exported by EDAX OIM. ND is the sample *z* axis. *K*τ and *W* depend only on the ND direction in the crystal; the in-plane stress components also depend on φ₁.

The inverse pole figure is a stereographic projection of the 001–100–110 standard triangle. θ[001] is the angle between ND and [001]; φ is the azimuth about [001], 0° on the 001–100 edge and 45° on the 001–110 edge.

## Verification

For [001] ∥ ND the problem is equi-biaxial and has the closed-form solution *W* = *M* α<sub>a</sub>², with *M* = *C*₁₁ + *C*₁₂ − 2*C*₁₃²/*C*₃₃. The page shows this comparison. The JavaScript and Python implementations agree to machine precision.

| ND ∥ | *K*τ (MPa K⁻¹) | *W* (J m⁻³ K⁻²) |
|---|---|---|
| [001] | 0.179 | 9.46 |
| [101] | 1.017 | 39.6 |
| [100] | 0.722 | 45.5 |
| [111] | 0.953 | 52.5 |
| [110] | 0.998 | 65.5 |

## Reference

J. W. Xian et al., "The role of microstructure in the thermal fatigue of solder joints," *Nature Communications* **15** (2024). https://doi.org/10.1038/s41467-024-48532-6

## License

MIT License. See `LICENSE`.
