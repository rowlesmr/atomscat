Creates a CIF file detailing the Inverse Mott-Bethe coefficients for the atomic scattering factors of all elements, most cations, and a couple of anions.

The coefficients can be used to calculate the atomic scattering factor for a given element/ion by

$$
f(s; Z_0) = Z_0 - 8 \pi r_B s^2 \left[ c + \sum_{i=1}^N a_i \exp\left(-b_i s^2 \right) \right]
$$

where $s=\dfrac{\sin\theta}{\lambda}$, $\theta$ is the Bragg angle, $\lambda$ is the X-ray wavelength in angstroms, $Z_0$ is the number of electrons, $r_B$ is the Bohr radius (0.529177210544 Å), $a_i$ is the list of Gaussian scalar coefficients, $b_i$ is the list of Gaussian exponent coefficients, and $c$ is a constant.


Coefficients come from 
- Thorkildsen, G. (2023). Acta Cryst. A79, 318-330 https://doi.org/10.1107/S2053273323003996
- Thorkildsen, G. (2023). Acta Cryst. A79, 318-330 https://doi.org/10.1107/S2053273323010550.

Original data come from 
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Acta Cryst. A79, 59-79 https://doi.org/10.1107/S2053273322010944
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Acta Cryst. A79, 229-245 https://doi.org/10.1107/S205327332300116X.

There is a new paper out, and once those data are reduced, this will be updated with those (Copeland, H., Watanabe, Y. & Volkov, A. (2026) Acta Cryst. A82, 200-217; https://doi.org/10.1107/S2053273326003554).
Raw data available: https://doi.org/10.5281/zenodo.20289472


