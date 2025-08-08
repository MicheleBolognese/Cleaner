within PEMFCModel4par;

model CoolStack "PEMFC stack model with cooling"
  
  extends .PEMFCModel.Stacks.Templates.CoolStack(
    redeclare model SubStack = PEMFCModel4par.SubStack,
    subStack(
    c1 = c1,
    alpha = alpha,
    j_0 = j_0,
    j_loss = j_loss,
    m_conc = m_conc,
    n_conc = n_conc),
    coolingPipe(
      redeclare replaceable model Friction = Friction_cooling,
      redeclare replaceable model HeatTransfer = HeatTransfer_cooling,
      NA = 1,
      NB = 1));
  
  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation (Dialog(enable = enable_setting,tab = "Substack",group = "Characteristics and polarization"));
  input Real c1(unit = "V/K", start = 0.85e-3) "Voltage derivative by temperature" annotation (Dialog(enable = enable_setting,tab = "Substack",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_0 (start = 1) "Exchange current density for activation loss" annotation (Dialog(enable = enable_setting,tab = "Substack",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_loss (start = 100) "Activation current density loss" annotation (Dialog(enable = enable_setting,tab = "Substack",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.Voltage m_conc (start = 3e-4) "Pre-exponential factor for concentration loss" annotation (Dialog(enable = enable_setting,tab = "Substack",group = "Characteristics and polarization"));
  input Real n_conc (unit = "m2/A", start = 3.2e-4) " Exponential factor for concentration loss" annotation(Dialog(group = "Characteristics and polarization",tab = "Substack",enable = enable_setting));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>PEMFCStack</h4>
<p>Predefined proton exchange membrane fuel cell (PEMFC) stack. Extends from <a href=\"modelica://FuelCell.Stacks.Templates.CoolStack\">FuelCell.Stacks.Templates.CoolStack</a>.
The model contains the following components:</p>
<p><ul>
<li>TwoPhase water cooling pipe</li>
<li>Feed and drain connectors of reaction gas type (ideal gas) for both anode and cathode</li>
<li>Feed and drain connectors of twophase water type for the cooling channel</li>
<li>Electrical pins (positive and negative)</li>
<li>Proton exchange membrane fuel cell substack model</li>
</ul></p>
<p>The model supports both co-flow and counter-flow configurations. The default setup is co-flow, to obtain counter-flow the flow in the cathode channel is reversed.</p>
<h4>Parametrization</h4>
<p>The discretization of the model is defined in the following way:</p>
<p><ul>
<li>The parameter <code>N</code>, which describes the number of discretizations in the flow directions (the membrane is divided in <code>N</code> segments) </li>
<li>The parameter <code>n_cell</code>, which describes the number of cells in the substack; the substack can be used to describe a single cell, a lumped subset of a stack or a complete stack</li>
</ul></p>

<p><br/><b><font style=\"color: #d2492a; \">Substack model</font></b> </p>
<p>The stack is parameterized to use a (<a href=\"modelica://FuelCell.Stacks.PEMFC.SubStack\">FuelCell.Stacks.PEMFC.PEMFCSubStack</a>) substack.
</p>
<p><br/><b><font style=\"color: #d2492a; \">Cooling pipe</font></b></p>
<p>TwoPhase water cooling pipe discretized in <code>N</code> volumes in the fluid flow direction (<a href=\"modelica://Modelon.ThermoFluid.FlowChannels.DistributedTwoPhase\">Modelon.ThermoFluid.FlowChannels.DistributedTwoPhase</a>). Friction and heat transfer models in the pipe are partial and needs to be redeclared before use.</p>
</html>"));
    
end CoolStack;
