within PEMFCModel.Membrane.PEMFC;
model ElectroChem "PEMFC | Electrochemical model"
  
  extends .PEMFCModel.Membrane.Templates.PartialMembrane(
    final n_e = n_e_exch,
    final S_reac_an = Medium_an.stoichiometry({"H2"}, {-1}),
    final S_reac_cath = Medium_cath.stoichiometry({"H2O","O2"}, {1,-0.5}),
    final an_names = {"H2","H2O"},
    final cath_names = {"O2","H2O"},
    redeclare model ConcentrationLoss = PEMFCModel.Membrane.Losses.ConcentrationLoss.PEMFCElectroChem(
      alpha = alpha,
      j_lim = j_lim,
      final T_cell = T_cell),
    redeclare model OhmicLoss = PEMFCModel.Membrane.Losses.OhmicLoss.PEMFCElectroChem(
      final T_cell = T_cell, 
      final lambda_mean = waterContent.lambda,
      final z = z),
    redeclare model ActivationLoss = PEMFCModel.Membrane.Losses.ActivationLoss.PEMFCElectroChem(
      alpha = alpha,
      n_e = n_e,
      final z = z,
      j0_ref = j0_ref,
      Rc = Rc,
      gamma_ORR = gamma_ORR,
      Ea_ORR = Ea_ORR,
      Ea_k_H2 = Ea_k_H2,
      A_k_H2 = A_k_H2,
      final T_cell = T_cell,
      final pH2 = pH2,
      final pO2 = pO2),
    redeclare replaceable model WaterDiffusion = PEMFCModel.Membrane.MassTransport.WaterDiffusion.Dutta,
    redeclare replaceable model ElectroOsmoticDrag = .PEMFCModel.Membrane.MassTransport.ElectroOsmoticDrag.Springer);

  // Parameters of membrane and catalyst
  final parameter Real n_e_exch = 2 "Number of exchanged electrons";
  parameter .Modelica.Units.SI.CurrentDensity j0_ref = 5e-5 "Reference exchange current density for activation loss" annotation (Dialog(group="Characteristics and polarization"));
  parameter Real Rc = 100 "Roughness factor of the cathode, i.e. catalyst active area per unit cathode area" annotation (Dialog(group="Characteristics and polarization"));
  parameter Real gamma_ORR = 0.5 "O2 reduction reaction dependency factor on O2 partial pressure" annotation (Dialog(group="Characteristics and polarization"));
  parameter .FuelCell.Internal.Units.ActivationEnergy Ea_ORR = 66e3 "Activation energy for ORR on the catalyst" annotation (Dialog(group="Characteristics and polarization"));
  parameter .FuelCell.Internal.Units.ActivationEnergy Ea_k_H2 = 21.03e3 "Activation energy for H2 permeation" annotation (Dialog(group="Characteristics and polarization"));
  parameter Real A_k_H2 = 6.6e-11 "Pre-exponetial factor for H2 permeation" annotation (Dialog(group="Characteristics and polarization"));
  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation (Dialog(group="Characteristics and polarization"));
  parameter .Modelica.Units.SI.CurrentDensity j_lim = 20000 "Limiting exchange current density" annotation (Dialog(group="Characteristics and polarization"));

  .Modelica.Units.SI.MolarInternalEnergy g_reaction[N] "Gibbs free energy of reaction";

protected
    
  constant .Modelica.Units.SI.Temperature Tref_25C = 298.15 "Reference temperature at 25C";

  .Modelica.Units.SI.Pressure pH2O[N] = p_cath_partial[:, iCath[2]];
  .Modelica.Units.SI.Pressure pH2O_an[N] = p_an_partial[:, iAn[2]];
  .Modelica.Units.SI.Pressure pH2[N] = p_an_partial[:, iAn[1]];
  .Modelica.Units.SI.Pressure pO2[N] = p_cath_partial[:, iCath[1]];

equation
    
  assert(pin_n.i >= -1e-5,"Cell external current should be positive, check operating conditions.");

  //Nernst potential
  g_reaction = cathode.g_formation[iCath[2]] - anode.g_formation[iAn[1]] - 0.5*cathode.g_formation[iCath[1]];
  E0_cell = -1/(n_e*FuelCell.Internal.Units.F)*(g_reaction + .Modelica.Constants.R*T_cell .* .Modelica.Math.log(pH2O*.FuelCell.Internal.Units.pNTP
    ^0.5 ./ (pH2 .* pO2 .^ 0.5)));

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<h4>PEMFC ElectroChem</h4>
<p>This is a model of a proton exchange membrane fuel cell. The model has the following connectors:</p>
<ul>
<li>Component mass transfer connectors of reaction gas type (ideal gas) to the anode and cathode sides</li>
<li>Electrical pins (positive and negative)</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Heat port connectors for heat transfer between membrane and flow channel </li>
</ul>
<h4>Parametrization</h4>
<p>The model has parameters determining the electrochemical characteristics, which can be tuned by the user.</p>
<p>For example,the catalyst properties of the model is defined here, by the following parameters:</p>
<ul>
<li>The parameter <span style=\"font-family: Courier New;\">Rc</span>, which describes the catalyst roughness factor of the cathode surface, i.e. area of catalyst per unit area of the cathode</li>
<li>The parameter <span style=\"font-family: Courier New;\">j0_ref</span>, which describes the exchange current density of the metal catalyst.</li>
</ul>
<p>There is also a set of predefined parameters, which the user should not override, as they are common for all membranes of a stack and therefore should be defined for the entire stack.</p>
<p>For example,the discretization of the model is defined here, by the following parameters: </p>
<ul>
<li>The parameter <span style=\"font-family: Courier New;\">N</span>, which describes the number of discretizations in the flow directions (the membrane is divided in <span style=\"font-family: Courier New;\">N</span> segments) </li>
<li>The parameter <span style=\"font-family: Courier New;\">n_cell</span>, which describes the number of cells in one membrane; the membrane can be used to describe a single cell membrane, or a number of cells (for lumped stacks) </li>
</ul>
<h4>Assumptions</h4>
<h5>Cell Voltage Model</h5>
<p>The cell voltage V<sub>cell</sub> is obtained from</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-m7OIyWl5.png\" alt=\"V_cell = E0 - E_ohmic - E_act - E_conc\"/></p>
<p>where E<sub>ohmic</sub>, E<sub>act</sub> and E<sub>conc</sub> are the Ohmic loss, the activation loss and the concentration loss.</p>
<p>E0 is the Nernst potential given by</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-fSrt1Kqc.png\" alt=\"E0 = -deltaG/(2*F) -R*T/(2*F) *ln(p_H2O*p_ref^0.5/(p_H2*p_O2^0.5))\"/></p>
<p>where deltaG is Gibbs free energy and p<sub>ref</sub> is the standard pressure at 0.1MPa.</p>
<h5>Activation Loss</h5>
<p>Assuming overvoltage at the anode is negligible compared to that of the cathode, the exchange current density is parameterized by the activation energy <span style=\"font-family: Courier New;\">Ea_ORR</span> and oxygen partial pressure dependency <span style=\"font-family: Courier New;\">gamma_ORR</span> of the oxygen reduction reaction. Fuel crossover is also taken into account by a model describing hydrogen permeation through Nafion membrane.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-t0bxqVlt.png\" alt=\"E_act = RT/(alpha*nF)*ln((j+j_loss)/j_0)\"/></p>
<h5>Ohmic Loss</h5>
<p>Humidity effects at both the anode and cathode are included in the calculation of internal resistance (ASR) given as:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-AehaYlYY.png\" alt=\"ASR = z/sigma\"/></p>
<p>where z is the membrane thickness and &sigma; is the membrane conductivity, based on the water content &lambda; and cell temperature, fitted as follows:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-JZ9tzsiN.png\" alt=\"sigma = (0.5139*lambda-0.326)*exp(1268*(1/303 - 1/T))\"/></p>
<p>A simple average of water vapor activities at the anode and cathode is used to calculate the average water content of the membrane.</p>
<h5>Concentration Loss</h5>
<p>Concentration loss is parameterized by the maximum exchange current density <span style=\"font-family: Courier New;\">j_L</span> and charge transfer coefficient <span style=\"font-family: Courier New;\">alpha </span>and is given by:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-GO84oQpM.png\" alt=\"E_conc = (RT/(nF))*(1+ 1/alpha)*ln(j_L/(j_L-j))\"/></p>
<h5>Heat Transfer Model</h5>
<p>All electro-chemical heat is absorbed by the cell and then dissipated through heat connectors. </p>
<h5>Mass Transfer</h5>
<p>The mass transfers to the anode and cathode channels are calculated in the model using Faraday&apos;s law of electrolysis.</p>
<p>The following water and gas transport implementations are provided:
<ul>
<li>Water content</li>
<li>Water diffusion</li>
<li>Water permeation</li>
<li>Electro-osmotic drag</li>
<li>Gas diffusion</li>
</ul>
Details on the implemented correlations could be found under the subpackage <a href=\"modelica://FuelCell.Membranes.MassTransport\">Membranes.MassTransport</a>.</p>
<h5>Polarization curve</h5>
<p>The polarization curve of the default parametrization for single cell is shown below. </p>
<p><img src=\"modelica://FuelCell/Resources/images/membranes/PEMFCElectroChem.png\"/></p>
<h4>References</h4>
<dl><dt>Springer, T. E., Zawodzinski, T. A., and Gottesfeld, S.:</dt>
<p style=\"margin-left: 30px;\"><h4>Polymer electrolyte fuel cell model.</h4></p>
<p style=\"margin-left: 30px;\">Journal of the Electrochemical Society, 1991, 138(8), 2334-2342.</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1149/1.2085971\">10.1149/1.2085971</a></p>
<dl><dt>Kocha, S.S., Deliang Yang, J. and Yi, J.S.:</dt>
<p style=\"margin-left: 30px;\"><h4>Characterization of gas crossover and its implications in PEM fuel cells.</h4></p>
<p style=\"margin-left: 30px;\">AIChE Journal, 2006, 52(5), 1916-1925.</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1002/aic.10780\">10.1002/aic.10780</a></p>
<dl><dt>Pukrushpan, Jay Tawee:</dt>
<p style=\"margin-left: 30px;\"><h4>Modeling and control of fuel cell systems and fuel processors.</h4></p>
<p style=\"margin-left: 30px;\">PhD thesis, University of Michigan, 2003.</p>
<p style=\"margin-left: 30px;\"><a href=\"https://api.semanticscholar.org/CorpusID:116943188\">Link to reference</a></p>
<dl><dt>Buck, A.L.:</dt>
<p style=\"margin-left: 30px;\"><h4>Buck Research CR-1A User&apos;s Manual</h4></p>
<p style=\"margin-left: 30px;\">Appendix 1, pp. 21.</p>
<p style=\"margin-left: 30px;\"><a href=\"https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf\">Link to reference</a></p>
<dl><dt>O&apos;Hayre, R., Cha, S., Colella, W.G. and Prinz, F.B.:</dt>
<p style=\"margin-left: 30px;\"><h4>Fuel Cell Fundamentals</h4></p>
<p style=\"margin-left: 30px;\">Wiley, 2016, pp-203-236</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1002/9781119191766\">10.1002/9781119191766</a></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end ElectroChem;
