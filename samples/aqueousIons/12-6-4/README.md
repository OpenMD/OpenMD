**12-6-4 Aqueous Ion Solutions**

The interactions present in aqueous ion solutions need a bit of care. Charge-charge and charge-dipole terms are generally covered by the direct Coulomb interactions between point charges on the ions and the water molecules. The 12-6-4 model includes contributions from ion-induced dipole interactions in addition to the standard 12-6 van der Waals interaction. The modified non-bonded interaction looks like:

$$V_{12-6-4} = \frac{C_{12}}{r^{12}} + \frac{C_6}{r^6}  + \frac{C_4}{r^4}$$

which is handled as an inverse power series in OpenMD. The $C_{12}$, $C_6$, and $C_4$ parameters between an ions and the relevant water model are part of the force field for 12-6-4 models.  However for ionic solutions, we need to generate **_all_** pairs of $C_4$ terms for the inter-ionic contributions. We use the following model:

$$C_{4}(i-j) = C_{4}(i-H2O)  \left(\frac{\alpha_j}{\alpha_{H2O}}\right)$$ 
which was suggested in the original 12-6-4 paper:

- "Systematic Parameterization of Monovalent Ions Employing the Nonbonded Model," Pengfei Li, Lin Frank Song, and Kenneth M. Merz Jr., *Journal of Chemical Theory and Computation*, **11** (4), 1645-1657, (2015) DOI: 10.1021/ct500918t 

A few pieces of information are required:

1. Atomic polarizability values $(\alpha_j)$. These can be taken either from Sangster & Atwood:

    - M J L Sangster and R M Atwood, *J. Phys. C: Solid State Phys*, **11**, 1541 (1978) DOI 10.1088/0022-3719/11/8/015

    or from Miller:

    - "Additivity methods in molecular polarizability," Kenneth J. Miller, *Journal of the American Chemical Society*, **112** (23), 8533-8542 (1990). DOI: 10.1021/ja00179a044

2. Water polarizability values $(\alpha_{H2O})$.  We are using $1.444 Å^3$, which was taken from Eisenberg and Kauzmann:

    - Eisenberg, D. S.; Kauzmann, W. *The Structure and Properties of Water* (Oxford University Press: Oxford, U.K., 1969)

3. $C_4$ values for ion $i$ with the water model. These are taken from the relevant force field for a specific water model, and are model-dependent. In the example in this directory, we have taken the $C_4$ values for Na+ and Cl- ions from the `LiMerzIons12-6-4.SPCE.frc` file in the `forceFields` directory.

Additionally, the $C_{12}$ and $C_6$ parameters are also required when we override Lennard-Jones interactions with explicit pair interactions. These can be generated as follows:

$$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}$$
$$\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}$$
$$C_{12} = 4 \epsilon_{ij} \sigma_{ij}^{12}$$
$$C_{6} = - 4 \epsilon_{ij} \sigma_{ij}^{6}$$

where $\sigma$ and $\epsilon$ values for each atom are given in the frc file. Once all relevant parameters have been calculated for a pair of atom types, they should be added to the `NonBondedInteractionTypes` section of the frc file:

```
atype1 atype2 InversePowerSeries  12  C12   6  C6   4  C4
```
where the `C12`, `C6`, and `C4` are replaced with the values calculated above.