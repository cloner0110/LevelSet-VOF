# Coupled Level-Set / VOF Solver for OpenFOAM

![OpenFOAM](https://img.shields.io/badge/OpenFOAM-v2312--v2606-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-green)

`interLevelSetFoam` simulates two incompressible, immiscible fluids. The volume
fraction is transported with VOF/MULES, which conserves mass. At every outer
corrector, a signed-distance level set is rebuilt from that volume fraction.
The level set then supplies the interface curvature and the smooth density
and viscosity that the momentum equation uses.

![rising bubble](Tutorial/bubbleColumn/result/gifLSVOF.gif)

## How it works

All of the level-set work lives in one class, `levelSetCoupledMixture`
(`src/levelSetCoupledMixture`). It extends OpenFOAM's
`immiscibleIncompressibleTwoPhaseMixture`, so the rest of the solver
(MULES, turbulence, `fvOptions`, dynamic meshes) is unchanged. With `h` the
reference cell size, every rebuild does the following:

| Step | Operation |
|---|---|
| Seed | `d0 = (2 alpha1 - 1) * seedCoeff * h` |
| Redistance | `d <- d + dTau * sgn(d0) * (1 - \|grad d\|)`, `dTau = pseudoTimeCoeff * h`, `nRedistanceIter` times |
| Heaviside | `H = (1 + x + sin(pi x)/pi) / 2`, `x = clamp(d/eps, -1, 1)`, `eps = widthCoeff * h` |
| Curvature | `kappa = -div(n_f)`, `n_f = (grad d)_f / \|(grad d)_f\|` |
| Surface tension | `f_f = sigma * kappa_f * snGrad(H)` |
| Properties | `rho = H rho1 + (1-H) rho2`, `mu = H mu1 + (1-H) mu2` |

The density and the mass flux `rhoPhi` both come from `H`. The momentum
equation also removes the remaining mass imbalance,
`ddt(rho) + div(rhoPhi)`. Together, these keep `ddt(rho U)` and
`div(rhoPhi U)` consistent with each other. The coupling follows the
simple coupled level-set/VOF approach of Albadawi et al. (2013).

## Settings

Every entry is optional. They all go in a `levelSet` sub-dictionary of
`constant/transportProperties`. The surface tension coefficient is the usual
`sigma` entry.

```
levelSet
{
    cellSize         0.00625;   // h [m]; default: smallest cell spacing
    widthCoeff       1.5;       // interface half-width eps / h
    seedCoeff        0.75;      // |seed distance| / h
    pseudoTimeCoeff  0.1;       // redistance pseudo-time step / h
    nRedistanceIter  15;        // default: ceil(widthCoeff/pseudoTimeCoeff)
    curvatureFrom    levelSet;  // levelSet | volumeFraction
    propertiesFrom   levelSet;  // levelSet | volumeFraction
}
```

On graded meshes, set `cellSize` explicitly to the cell size near the
interface. The solver writes three extra fields: `levelSet` (signed distance),
`heaviside` and `levelSetCurvature`.

## Build

```bash
git clone https://github.com/cloner0110/LevelSet-VOF.git
cd LevelSet-VOF/src && wmake
```

Tested with OpenFOAM.com v2306, v2312 (CI), v2512 (CI) and v2606.

## Tutorial: 2-D rising bubble

`Tutorial/bubbleColumn` is test case 2 of Hysing et al. (2009). A circular
bubble of radius 0.25 starts at (0.5, 0.5) in a closed 1 × 2 box, with
ρ1/ρ2 = 1000, μ1/μ2 = 100, σ = 1.96 and g = 0.98.
- The bubble is set with `setAlphaField`, which integrates the volume fraction
  in each cut cell.
- The top and bottom walls are no-slip, the side walls are free-slip, and
  p_rgh is `fixedFluxPressure` on all walls.
- The case writes the bubble's centre of mass and mean rise velocity at every
  time step, to `postProcessing/bubble/0/volFieldValue.dat`.

```bash
cd Tutorial/bubbleColumn && ./Allrun
```

The results below were run with OpenFOAM v2606. The 1/h = 160 runs are the
tutorial as shipped, on 4 MPI ranks. The 1/h = 80 runs use `(80 160 1)` cells
and `deltaT 7.8125e-4`. The reference is TP2D on its finest mesh, from Hysing
et al. (2009), Tables XIII and XVI.

| Run | V<sub>max,1</sub> (t) | V<sub>max,2</sub> (t) | y<sub>c</sub>(t = 3) |
|---|---|---|---|
| TP2D, 1/h = 640 (reference) | 0.2524 (0.733) | 0.2434 (2.071) | 1.1380 |
| FreeLIFE, 1/h = 160 | 0.2514 (0.728) | 0.2440 (1.984) | 1.1249 |
| MooNMD | 0.2502 (0.732) | 0.2393 (2.060) | 1.1376 |
| this solver, 1/h = 80 | 0.2491 (0.734) | 0.2287 (1.913) | 1.1140 |
| this solver, 1/h = 160 | 0.2472 (0.738) | 0.2315 (1.982) | 1.1174 |
| volume-fraction coupling, 1/h = 160 | 0.2471 (0.734) | 0.2309 (2.001) | 1.1140 |

The phase volume is conserved to all printed digits.
- With `curvatureFrom` and `propertiesFrom` set to `volumeFraction`, the solver
  reproduces `interFoam`.
- The level-set coupling follows the reference more closely than that setting
  after the first maximum.

Results published for this tutorial before the setup fix (PR #4) came from a
tilted, 3% oversized initial bubble and a top boundary that let fluid through.
They should not be used.

## References

- Albadawi, A., Donoghue, D.B., Robinson, A.J., Murray, D.B., Delauré, Y.M.C.
  (2013). Influence of surface tension implementation in Volume of Fluid and
  coupled Volume of Fluid with Level Set methods for bubble growth and
  detachment. *Int. J. Multiphase Flow* 53:11–28.
- Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S.,
  Tobiska, L. (2009). Quantitative benchmark computations of two-dimensional
  bubble dynamics. *Int. J. Numer. Meth. Fluids* 60:1259–1288.
- Taran, A. (2026). A consistent coupled level-set and volume-of-fluid solver
  for two-phase flow in OpenFOAM. *arXiv preprint* arXiv:2610.24225.
  [https://arxiv.org/abs/2610.24225v1](https://arxiv.org/abs/2610.24225v1)

## License

This solver is derived from OpenFOAM's `interFoam` and links against
OpenFOAM libraries. It is therefore distributed under the
[GNU General Public License v3](LICENSE), like OpenFOAM itself. Files that
derive from OpenFOAM keep the original OpenFOAM Foundation and OpenCFD
copyright lines.

OpenFOAM® is a registered trade mark of OpenCFD Ltd. This project is not
approved or endorsed by OpenCFD Ltd.

## Contact

- **Email**: [amirhossein.taran@ucdconnect.ie](mailto:amirhossein.taran@ucdconnect.ie)
- **Issues**: [GitHub Issues](https://github.com/cloner0110/LevelSet-VOF/issues)

## Citation

If you use this solver, or want to cite it, please cite the following paper:

> Taran, A. (2026). A consistent coupled level-set and volume-of-fluid solver
> for two-phase flow in OpenFOAM. *arXiv preprint* arXiv:2610.24225.
> [https://arxiv.org/abs/2610.24225v1](https://arxiv.org/abs/2610.24225v1)

```bibtex
@misc{taran2026interLevelSetFoam,
  author        = {Taran, Amirhossein},
  title         = {A consistent coupled level-set and volume-of-fluid solver
                   for two-phase flow in {OpenFOAM}},
  year          = {2026},
  eprint        = {2610.24225},
  archivePrefix = {arXiv},
  primaryClass  = {physics.flu-dyn},
  doi           = {10.48550/arXiv.2610.24225},
  url           = {https://arxiv.org/abs/2610.24225v1}
}
```
