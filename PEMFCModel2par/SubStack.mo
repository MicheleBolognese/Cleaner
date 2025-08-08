within PEMFCModel2par;

model SubStack "PEMFC substack"
    
  extends .PEMFCModel.Stacks.Templates.SubStack(
    redeclare model Membrane = .PEMFCModel4par.Empirical, 
    cell(
      final h_conv_an = anode_channel.channel.alpha, 
      final h_conv_cath = cathode_channel.channel.alpha,
      final pstart = (p_start_in_anode + p_start_in_cathode)/2, 
      final Tstart = (T_start_in_anode + T_start_in_cathode)/2,
      final T_from_h = false,
      final c1 = c1,
      final alpha = alpha,
      final j_0 = j_0,
      final j_loss = j_loss,
      final m = m_conc,
      final n = n_conc),
      final addProxToAnode = false,
      redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.DetailedWallFriction,
      redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.DetailedWallFriction,
      redeclare replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
      redeclare replaceable model HeatTransfer_cathode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable);
  
  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation (Dialog(enable = enable_setting,tab = "Cell model",group = "Characteristics and polarization"));
  input Real c1(unit = "V/K", start = 0.85e-3) "Voltage derivative by temperature" annotation(Dialog(enable = enable_setting,tab = "Cell model",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_0 (start = 1) "Exchange current density for activation loss" annotation (Dialog(enable = enable_setting,tab = "Cell model",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.CurrentDensity j_loss (start = 100) "Activation current density loss" annotation (Dialog(enable = enable_setting,tab = "Cell model",group = "Characteristics and polarization"));
  input .Modelica.Units.SI.Voltage m_conc (start = 3e-4) "Pre-exponential factor for concentration loss" annotation (Dialog(enable = enable_setting,tab = "Cell model",group = "Characteristics and polarization"));
  input Real n_conc (start = 3.2e-4) " Exponential factor for concentration loss" annotation(Dialog(group = "Characteristics and polarization",tab = "Cell model",enable = enable_setting));
    
  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>PEMFCSubStack</h4>
<p>Predefined proton exchange membrane fuel cell (PEMFC) substack with Butler-Volmer characteristics. Extends from <a href=\"modelica://FuelCell.Stacks.Templates.SubStack\">FuelCell.Stacks.Templates.SubStack</a>.
Can be used to describe a single cell or arbitrary number of cells. The model consists of the following components:</p>
<p><ul>
<li>Feed and drain connectors of reaction gas type (ideal gas) for both anode and cathode</li>
<li>Electrical pins (positive and negative)</li>
<li>Proton exchange fuel cell membrane model</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Fluid flow channels for the cathode and the anode side</li>
</ul></p>
<p>The model supports both co-flow and counter-flow configurations. The default setup is co-flow, to obtain counter-flow the flow on the cathode side is reversed.</p>
<h4>Parametrization</h4>
<p>The discretization of the model is defined in the following way:</p>
<p><ul>
<li>The parameter <code>N</code>, which describes the number of discretizations in the flow directions (the membrane is divided in <code>N</code> segments) </li>
<li>The parameter <code>n_cell</code>, which describes the number of cells in the substack; the substack can be used to describe a single cell, a lumped subset of a stack or a complete stack</li>
</ul></p>
<h4>Membrane Model</h4>
<p>The membrane is an essential part of the fuel cell as it describes the fuel cell reactions and mass transfer between anode and cathode.
<br/>The PEMFCSubStack uses the <a href=\"modelica://FuelCell.Membranes.PEMFC.Simplified\">FuelCell.Membranes.PEMFC.PEMFCSimplified</a> membrane by default. </p>
<h4>Flow Channels</h4>
<p>The anode and cathode channels are FuelCell.Pipes.FlowChannel_mass pipes and consists of <code>N</code> discrete volumes in the fluid flow direction, each with unique pressure, temperature and mass content. Each volume has a mass flow port for mass flow to the membrane and a heat flow port for heat transfer between channel and membrane No reactions occur in either of the channels.</p>
<h4>Conditional Prox Reactor</h4>
<p>The conditional preferential oxidation (prox) reactor before the anode inlet is enabled or disabled using the parameter <code>addproxToAnode</code>. The prox reactor is only needed for cases where extra air is added just before the stack to oxidate any excess fraction of carbon monoxide (CO), before the stack inlet. This is done to prevent the degradation of PEMFC membranes due to too high fractions of carbon monoxide. A typical PEM fuel cell can handle up to 25 ppm CO.</p>
</html>"));
    
end SubStack;
