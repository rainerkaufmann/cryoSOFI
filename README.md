## cryo-SOFI
<h3>Swapped mean and absolute value calculation bug in version 1.0</h3>
<br><br>
The cryo-SOFI reconstruction code version 1.0 contained a bug that can lead to images with slightly reduced resolution compared to the correct implementation.
<br><br>
In the cryo-SOFI code version 1.0, the order of calculating mean and absolute value for the correlation factors had accidentally been swapped. This favours noise contribution (e.g. shot-noise from non-fluctuating signals or background) in the resulting SOFI reconstructions. The following simulation with 9 molecules on a pattern with 320 nm spacing and a 415 nm wide (FWHM) PSF shows, that this can lead to a loss in achievable resolution compared to the reconstruction with code that uses the correct order of calculating mean and absolute value (cryoSOFI v1.1).
<br><br>
<img src="https://github.com/rainerkaufmann/cryoSOFI/blob/master/9mol_2on_4off_10000fs_comparison.png" width="1024">
<br>
Simulated dataset of 9 molecules on a pattern with 320 nm spacing and a PSF with a 415 nm FWHM. Simulation includes a constant background and a background of non-fluctuating signals at the molecule positions. Third order cross correlation SOFI (XC3) images were calculated using code as originally deposited on GitHub (v1.0), code with correct order of mean and absolute value calculations (v1.1) and an improved implementation (v1.5).
<br><br>
While the simulation indicates that the cryo-SOFI code version 1.0 can impair the achievable resolution, we did not observe this for our experimental datasets. The example below shows the results of the different cryo-SOFI software versions for Dendra2-Lifeact-labelled actin in a vitrified cell (same dataset as in Fig. 1 of the PNAS article https://doi.org/10.1073/pnas.1810690116).
<br><br>
<img src="https://github.com/rainerkaufmann/cryoSOFI/blob/master/actin_comparison.png" width="1024">
<br>
Cryo-SOFI of Dendra2-Lifeact-labelled actin in a vitrified XC cell (rat Rous sarcoma cell line).
<br><br>
The cryo-SOFI reconstruction code version 1.5 contains several additional improvements such as automatic calculation of pixel combinations, calculation of every tau combination between 0 and 4, and GPU based calculations to achieve much faster performance.

