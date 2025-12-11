# SIESTA

The Sea Ice Ecosystem State (SIESTA) model assimilates relevant satellite sea ice observations in an effort to produce a state estimate of sea ice physics. SIESTA includes dynamic, thermodynamic, and hydrodynamic processes necessary to simulate the vertical distribution and evolution of sea ice characteristics that include salinity, porosity, brine flux, and spectral light transmission and scattering [Arrigo et al., 1993; Arrigo and Sullivan, 1994; Bitz et al., 2001; Holland et al., 2006; Saenz and Arrigo, 2012]. Sea ice, snow cover, and associated physics and tracers are represented in the vertical dimension by the sea ice ecosystem model of Saenz and Arrigo [2012] using a maximum of 42 horizontal ice layers and up to 26 snow layers. This ecosystem model includes macronutrient cycling (NO3, NH4, P, Si), and an ice alga species represented by empirical ice community parameterizations, including a spectral-dependent photon absorption and photoadaptation.

Note, this code was written a long time ago under pressure of limited compute resources and a 'just get it done' outlook of a PhD program. Apologies for the single-massive-subroutine structure - it optimized best on the computer I had.

For a clean, modularized version of the SIESTA FORTRAN 1-dimensional (vertical) snow and ice physics, look at the KEIpy (or KEI) model. The sea ice ecosystem code and sea ice dynamics code as still only found in this repository.

Please reach out to blsaenz if you have questions.

References:

Saenz BT, Arrigo KR (2012) Simulation of a sea ice ecosystem using a hybrid model for slush layer desalination. Journal of Geophysical Research 117:5007.

Saenz BT, Arrigo KR (2014) Annual primary production in Antarctic sea ice during 2005-2006 from a sea ice state estimate. Journal of Geophysical Research Oceans 119:3645–3678.

Tedesco et al. Intercomparison study reveals pathways for improving the representation of sea-ice biogeochemistry in models.  In Review.
