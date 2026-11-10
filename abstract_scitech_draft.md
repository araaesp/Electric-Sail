# Draft — awaiting comparison with the previous abstract

## Title

Solar Wind Variability and Electric Sail Performance: Numerical Analysis of Orbital Energy and Thrust Sensitivity

## Abstract

Electric sails provide continuous propulsion through the interaction between positively charged tethers and the solar wind. Their performance depends on plasma conditions and applied tether voltage, making solar wind variability relevant to interplanetary mission analysis. This work investigates how solar wind density, temperature, and speed affect electric sail thrust and the evolution of spacecraft orbital energy, and examines whether historical wind conditions reveal favorable launch periods.

A numerical simulator couples gravitational dynamics with an electric sail thrust model and integrates the equations of motion using an adaptive Runge–Kutta method. Historical solar wind measurements are screened and interpolated to construct time dependent reference conditions at one astronomical unit. An approximate radial thrust dependence proportional to heliocentric distance raised to the power of minus seven sixths is applied. Missions begin from the same approximately circular heliocentric state at one astronomical unit, with a fixed sail attitude, allowing comparisons under controlled initial conditions.

Three complementary experiments are considered. First, six year missions launched on different dates are compared through their final specific orbital energy. Second, four ten year scenarios combine independent maxima and minima of wind properties to explore limiting model responses; these combinations represent synthetic conditions rather than simultaneous historical observations. Third, a sensitivity study comprises 75 distinct ten year simulations. Each wind property is varied from its tenth to its ninetieth percentile in ten percentile increments, while the other two remain fixed at their tenth, fiftieth, or ninetieth percentiles. These simulations use a tether voltage of 4 kV. A separate thrust calculation examines eight voltages between 2 and 20 kV.

The historical launch comparisons produce both positive and negative energy differences, providing no simple ordering of launch dates within the sampled cases. The controlled experiments demonstrate that simultaneously maximizing wind properties does not necessarily maximize orbital energy gain. In the adopted thrust model, increasing wind speed can reduce thrust beyond a voltage dependent maximum, and increasing tether voltage shifts this maximum toward higher wind speeds. These findings motivate assessing launch conditions through their combined propulsive effect rather than through a single indicator of solar activity. The results characterize model sensitivity and support further investigation of launch opportunities, without establishing a universal optimal launch date.

## Evidence and pending work

- Source: current code and saved outputs of cells 16–19 in `src/analise_perfil_real_direto_invertido.ipynb`.
- Cell 16: 34 distinct missions and 45 pairs; each separation group includes energy differences of both signs.
- Cell 19: saved output records completion of 75 ten year missions, fixed attitude phi = -30 degrees and theta = 0 degrees, at 4000 V.
- Cell 18: saved thrust maxima shift from approximately 145 km/s at 2 kV to 462.5 km/s at 20 kV.
- Cell 20 contains the P10–P90 amplitude analysis but has no saved outputs; its numerical results are not claimed here.
- The processed dataset contains interpolated missing values. Percentiles therefore describe the processed series.
- Pending: receive and analyze the previous congress abstract, compare its scope and claims with current work, revise this draft, and enter the agreed base text in the open conference form without submitting.
