within PEMFCModel4par;
model Empirical "PEMFC | Empirical models"
  
  extends .PEMFCModel.Membrane.Templates.PartialMembrane(
    final n_e = n_e_exch,
    final S_reac_an = Medium_an.stoichiometry({"H2"}, {-1}),
    final S_reac_cath = Medium_cath.stoichiometry({"H2O","O2"}, {1,-0.5}),
    final an_names = {"H2","H2O"},
    final cath_names = {"O2","H2O"}, 
    redeclare model ConcentrationLoss = PEMFCModel4par.PEMFCExponentialFit(m = m, n = n), 
    redeclare model OhmicLoss = .PEMFCModel.Membrane.Losses.OhmicLoss.PEMFCElectroChem(final T_cell = T_cell, final lambda_mean = waterContent.lambda),
    redeclare model ActivationLoss = PEMFCModel4par.PEMFCEmpirical(
      alpha = alpha, 
      j_loss = j_loss,
      j_0 = j_0,
      final T_cell = T_cell),
    redeclare replaceable model ContaminantsLoss = .PEMFCModel.Membrane.Losses.ContaminantsLoss.ZeroLoss,redeclare replaceable model WaterDiffusion = .PEMFCModel.Membrane.MassTransport.WaterDiffusion.Dutta,redeclare replaceable model ElectroOsmoticDrag = .PEMFCModel.Membrane.MassTransport.ElectroOsmoticDrag.Springer,CF_N_eod = fill(CF_0_eod,N),CF_0_eod = 0.8);
  
  final parameter Real n_e_exch = 2 "Number of exchanged electrons";
  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  input Real c1(unit "V/K", start = 0.85e-3) "Voltage derivative by temperature" annotation(Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_0 (start = 1) "Exchange current density for activation loss" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_loss (start = 100) "Activation current density loss" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  input .Modelica.Units.SI.Voltage m (start = 3e-4) "Pre-exponential factor for concentration loss" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  input Real n (unit = "m2/A", start = 3.2e-4) "Exponential factor for concentration loss" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
    
protected
  
  constant .Modelica.Units.SI.Temperature Tref_25C = 298.15 "Reference temperature at 25°C";
  .Modelica.Units.SI.Pressure pH2O_cath[N] = p_cath_partial[:, iCath[2]];
  .Modelica.Units.SI.Pressure pH2O_an[N] = p_an_partial[:, iAn[2]];
  .Modelica.Units.SI.Pressure pH2[N] = p_an_partial[:, iAn[1]];
  .Modelica.Units.SI.Pressure pO2[N] = p_cath_partial[:, iCath[1]];
    
equation
  
  assert(pin_n.i>= -1e-5,"Stack external current should be positive, check operating conditions!");

  for i in 1:N loop
    
    E0_cell[i] = E0_ref - c1*(T_cell[i] - Tref_25C) + 4.3085e-5*T_cell[i]*(.Modelica.Math.log(max(.Modelica.Constants.eps, pH2[i]/.FuelCell.Internal.Units.pNTP)) + 0.5*.Modelica.Math.log(max(.Modelica.Constants.eps,pO2[i]/.FuelCell.Internal.Units.pNTP)));
  
  end for;

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Documentation(revisions="<html>Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.</html>", info="<html>
<h4>PEMFC Empirical</h4>
<p>This is a model of a proton exchange membrane fuel cell. The model has the following connectors:</p>
<ul>
<li>Component mass transfer connectors of reaction gas type (ideal gas) to the anode and cathode sides</li>
<li>Electrical pins (positive and negative)</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Heat port connectors for heat transfer between membrane and flow channel </li>
</ul>
<h4>Parametrization</h4>
<p>The membrane consists of several parameters whose default values are derived experimentally using the <a href=\"modelica://FuelCell.Membranes.PEMFC.ElectroChem\">ElectroChem membrane</a> as the reference.</p>
<p>There is also a set of predefined parameters, which the user should not override, as they are common for all membranes of a stack and therefore should be defined for the entire stack.</p>
<p>For example,the discretization of the model is defined here, by the following parameters: </p>
<ul>
<li>The parameter <span style=\"font-family: Courier New;\">N</span>, which describes the number of discretizations in the flow directions (the membrane is divided in <span style=\"font-family: Courier New;\">N</span> segments) </li>
<li>The parameter <span style=\"font-family: Courier New;\">n_cell</span>, which describes the number of cells in one membrane; the membrane can be used to describe a single cell membrane, or a number of cells (for lumped stacks) </li>
</ul>
<h4>Assumptions</h4>
<h5>Cell Voltage Model</h5>
<p>The cell voltage V_cell is obtained from</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-m7OIyWl5.png\" alt=\"V_cell = E0 - E_ohmic - E_act - E_conc\"/></p>
<p>where E<sub>ohmic</sub>, E<sub>act</sub> and E<sub>conc</sub> are the Ohmic loss, the activation loss and the concentration loss.</p>
<p>E0 is the Nernst potential given by</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-0fXrbWRp.png\" alt=\"E0 = 1.229 -0.85e-3*(T-298.15)+4.3085e-5*T*[ln*p_H2+1/2*ln*p_O2]\"/></p>
<h5>Activation Loss</h5>
<p>Assuming overvoltage at the anode is negligible compared to that of the cathode, the activation loss is parameterized by the constant exchange current density <span style=\"font-family: Courier New;\">j_0</span>, the constant current density loss <span style=\"font-family: Courier New;\">j_n</span> and charge transfer coefficient <span style=\"font-family: Courier New;\">alpha</span>.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-jXVIbpS8.png\" alt=\"E_act = RT/(alpha*nF)*ln((j+j_n)/j_0)\"/></p>
<h5>Ohmic Loss</h5>
<p>Humidity effects at both the anode and cathode are included in the calculation of internal resistance (ASR) given as:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-AehaYlYY.png\" alt=\"ASR = z/sigma\"/></p>
<p>where z is the membrane thickness and &sigma; is the membrane conductivity, based on the water content &lambda; and cell temperature, fitted as follows:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-JZ9tzsiN.png\" alt=\"sigma = (0.5139*lambda-0.326)*exp(1268*(1/303 - 1/T))\"/></p>
<p>A simple average of water vapor activities at the anode and cathode is used to calculate the average water content of the membrane.</p>
<h5>Concentration Loss</h5>
<p>Concentration loss is parameterized by the experimentally derived values of m and n coefficients.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-S8QU47jP.png\" alt=\"E_conc = m*exp(n*j)\"/></p>
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
<p><img src=\"modelica://FuelCell/Resources/images/membranes/PEMFCEmpirical.png\"/></p>
<h4>References</h4>
<dl><dt>Dicks, A.L, and Rand, D.A.:</dt>
<p style=\"margin-left: 30px;\"><h4>Fuel cell systems explained, 2nd edition</h4></p>
<p style=\"margin-left: 30px;\">John Wiley &amp; Sons, 2018.</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1002/9781118878330\">10.1002/9781118878330</a></p>
<dl><dt>Pukrushpan, Jay Tawee:</dt>
<p style=\"margin-left: 30px;\"><h4>Modeling and control of fuel cell systems and fuel processors.</h4></p>
<p style=\"margin-left: 30px;\">PhD thesis, University of Michigan, 2003.</p>
<p style=\"margin-left: 30px;\"><a href=\"https://api.semanticscholar.org/CorpusID:116943188\">Link to reference</a></p>
<dl><dt>O&apos;Hayre, R., Cha, S., Colella, W.G. and Prinz, F.B.:</dt>
<p style=\"margin-left: 30px;\"><h4>Fuel Cell Fundamentals</h4></p>
<p style=\"margin-left: 30px;\">Wiley, 2016, pp-203-236</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1002/9781119191766\">10.1002/9781119191766</a></p>
</html>"));

end Empirical;