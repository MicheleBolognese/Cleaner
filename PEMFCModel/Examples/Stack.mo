within PEMFCModel.Examples;

model Stack "PEMFC stack model with cooling"
  extends .FuelCell.Stacks.Templates.CoolStack(
    redeclare replaceable model SubStack = .PEMFCModel.Examples.SubStack_1,
    L_cooling=00813,
    D_cooling=0.0375,
    N=4,
    n_cell=11,
    coolingPipe(
      redeclare model Friction = .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction,
      redeclare model HeatTransfer = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.CoefficientsPerPhase (k_2ph=10000,k_vap = 5000,k_liq = 5000),
      NA=1,
      NB=1),A_cell = 300e-4,M_stack = 10.1,positiveFlow_cooling = false,subStack(L_anode = 0.0813,D_anode = 0.023,L_cathode = 0.0813,D_cathode = 0.0445,redeclare replaceable model Membrane = .PEMFCModel.Examples.Empirical,addProxToAnode = false,redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 0.4,dp0 = 0.1e5,m_flow0 = 2.901e-4),redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 1.2,dp0 = 0.15e5,m_flow0 = 1.096e-3)));

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
end Stack;
