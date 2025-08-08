within PEMFCModel.Membrane.PEMFC;
model Simplified "PEMFC | Using simplified Butler-Volmer characteristics"
  
  extends PEMFCModel.Membrane.Templates.PartialMembrane(
    final n_e = n_e_exch,
    final S_reac_an = Medium_an.stoichiometry({"H2"}, {-1}),
    final S_reac_cath = Medium_cath.stoichiometry({"H2O","O2"}, {1,-0.5}),
    final an_names = {"H2"},
    final cath_names = {"O2","H2O"},
    redeclare model ConcentrationLoss = PEMFCModel.Membrane.Losses.ConcentrationLoss.ZeroLoss,
    redeclare model OhmicLoss = PEMFCModel.Membrane.Losses.OhmicLoss.PEMFCSimplified(R = R),
    redeclare model ActivationLoss = PEMFCModel.Membrane.Losses.ActivationLoss.PEMFCSimplified(a = a, j0 = j0));

  final parameter Real n_e_exch = 2.0 "Number of exchanged electrons";
  parameter .FuelCell.Internal.Units.AreaSpecificResistance R = 1e-6 "Area specific resistance" annotation (Dialog(group = "Characteristics and polarization"));
  parameter .Modelica.Units.SI.Voltage a = 0.0302 "Tafel slope" annotation (Dialog(group = "Characteristics and polarization"));
  parameter .Modelica.Units.SI.CurrentDensity j0 = 1 "Exchange current density" annotation (Dialog(group = "Characteristics and polarization"));

initial equation
  
  assert(pin_n.i>=-1e-5,"Cell external current should be positive, check operating conditions.");
    
equation
  
  E0_cell = ones(N)*E0_ref;

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<h4>PEMFCSimplified</h4>
<p>This is a model of a proton exchange membrane fuel cell. The model has the following connectors:</p>
<ul>
<li>Component mass transfer connectors of reaction gas type (ideal gas) to the anode and cathode sides</li>
<li>Electrical pins (positive and negative)</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Heat port connectors for heat transfer between membrane and flow channel </li>
</ul>
<h4>Parametrization</h4>
<p>The model has parameters determining the electrochemical characteristics, which the user can set. There is also a set of predefined parameters, which the user should not override, as they are common for all membranes of a stack and therefore should be defined for the entire stack. The discretization of the model is defined here, by the following parameters: </p>
<ul>
<li>The parameter <span style=\"font-family: Courier New;\">N</span>, which describes the number of discretizations in the flow directions (the membrane is divided in <span style=\"font-family: Courier New;\">N</span> segments) </li>
<li>The parameter <span style=\"font-family: Courier New;\">n_cell</span>, which describes the number of cells in one membrane; the membrane can be used to describe a single cell membrane, or a number of cells (for lumped stacks)</li>
</ul>
<h4>Assumptions</h4>
<h5>Cell Voltage Model</h5>
<p>The cell voltage V<sub>cell</sub> is obtained from open circuit voltage <i>E0 </i>subtracted for activation losses and ohmic resistance losses as follows</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-4doUZEiv.png\" alt=\"V_cell = E0 - E_act - R*j\"/></p>
<p>where <i>R</i> is the area specific resistance (ohm*m2)</p>
<p>Activation losses are based on the Butler-Volmer expression simplified to the Tafel expression:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-aCZp5rzZ.png\" alt=\"Eact=a*ln(j/j0+1)\"/></p>
<p>where <i>a</i> represents the Tafel slope, <i>j</i> the current density and <i>j0</i> the so called exchange current density. </p>
<h5>Heat Transfer Model</h5>
<p>All electro-chemical heat is absorbed by the cell and then dissipated through heat connectors.</p>
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
<p><img src=\"modelica://FuelCell/Resources/images/membranes/PEMFCSimplified.png\"/></p>
<h4>References</h4>
<dl><dt>O&apos;Hayre, R., Cha, S., Colella, W.G. and Prinz, F.B.:</dt>
<p style=\"margin-left: 30px;\"><h4>Fuel Cell Fundamentals</h4></p>
<p style=\"margin-left: 30px;\">Wiley, 2016, pp-203-236</p>
<p style=\"margin-left: 30px;\">DOI: <a href=\"https://doi.org/10.1002/9781119191766\">10.1002/9781119191766</a></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end Simplified;
