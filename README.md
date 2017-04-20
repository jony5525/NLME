# NLME
Objective functions that can be used for parameter estimation in Non Linear Mixed Effect Models

Examples based on code from Joakim Nyberg joakim.nyberg@farmbio.uu.se
The different OFV approximations are written by Joakim Nyberg but some of the code;
(EBE estimation, LinMatrixL, LinMatrixH, V) comes from the PopED sofware written by Andrew Hooker (andrew.hooker@farmbio.uu.se),
Joakim Nyberg (joakim.nyberg@farmbio.uu.se) and Sebastian Ueckert (sebastian.ueckert@farmbio.uu.se)

Likelihoods are assuming that parameters are described as point estimates from distributions (e.g. similar to the NONMEM software) as opposed to distributions with link functions (e.g. similar to the MONOLIX software).

References:
M Davidian, D Giltinan: Nonlinear Models for Repeated Measurements data, Chapman and Hall/CRC, 1989.
Wang, Y. J. Derivations of various NONMEM estimation methods Pharmacokinet Pharmacodyn (2007) 34: 575. doi:10.1007/s10928-007-9060-6.
M Lavielle, Mixed Effects Models for the Population Approach: Models, Tasks, Methods and Tools uly 14, 2014 by Chapman and Hall/CRC

