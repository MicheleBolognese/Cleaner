within SAMU;

model CoolStack_1 "Substack with cooling pipe"

  replaceable package Medium_an =
      .FuelCell.Media.PreDefined.IdealGases.NASAReformateShort
    constrainedby .FuelCell.Media.Templates.ReactionGas "Anode medium" annotation(choicesAllMatching);
  replaceable package Medium_cath =
      .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir                               constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathode medium"  annotation(choicesAllMatching);
  replaceable package Medium_cooling =
      .FuelCell.Media.PreDefined.TwoPhase.WaterIF97 "Water medium"
                    annotation(choicesAllMatching);

  .FuelCell.SubComponents.ComponentSummaries.CoolStackSummary summary(
    N=N,
    n_cell=n_cell,
    V_stack=subStack.summary.V_stack,
    P_stack=subStack.summary.P_stack,
    j_external=subStack.summary.j_external,
    T_stack=subStack.summary.T_stack,
    Q_stack=subStack.summary.Q_stack,
    anode_stoich=subStack.summary.anode_stoich,
    cathode_stoich=subStack.summary.cathode_stoich,
    m_flow_an=subStack.summary.m_flow_an,
    m_flow_cath=subStack.summary.m_flow_cath,
    T_cool=coolingPipe.T,
    m_flow_cool=-coolingPipe.portB[1].m_flow,
    dp_anode=(subStack.feed_anode.p - subStack.drain_anode.p),
    dp_cathode=(subStack.feed_cathode.p - subStack.drain_cathode.p),
    dp_cooling=(coolingPipe.portA[1].p - coolingPipe.portB[1].p))
    annotation (Placement(transformation(extent={{-76,74},{-56,94}})));

 /* Substack parameters */
  replaceable model SubStack = .FuelCell.Stacks.Templates.SubStack constrainedby .FuelCell.Stacks.Templates.SubStack
                                       "Substack model to use"                                                                             annotation(choicesAllMatching,Dialog(group="Stack configuration"));

  parameter Integer N(min = 1) = 1
    "Number of discretization nodes in the flow direction"                                annotation(Dialog(group="Stack configuration"));
  parameter Integer n_cell(min = 1) = 1 "number of cells in stack"               annotation(Dialog(group="Stack configuration"));

  parameter .Modelica.Units.SI.Area A_cell=180e-4 "active cell area" annotation (Dialog(group="Stack data"));
  parameter .Modelica.Units.SI.SpecificHeatCapacity Cp_cell=500 "specific heat capacity of cell material"
    annotation (Dialog(group="Stack data"));
  parameter .Modelica.Units.SI.CoefficientOfHeatTransfer kc_cell=250 "cell internal heat transfer coefficient"
    annotation (Dialog(group="Stack data"));

  parameter .Modelica.Units.SI.Mass M_stack=30 "total mass of stack" annotation (Dialog(group="Stack data"));

  /* Cooling pipe*/

  .Modelon.ThermoFluid.FlowChannels.DistributedTwoPhase coolingPipe(
    redeclare package Medium = Medium_cooling,
    n_channels=fill(n_channels_cooling, N),
    L=fill(L_cooling/N, N),
    Dhyd=fill(Dhyd_cooling, N),
    A=fill(A_cooling, N),
    levels=(1:N)*height_cooling/N,
    A_heat=fill(A_heat_cooling/N, N),
    includeStaticHead=includeStaticHead_cooling,
    g=g_cooling,
    initOpt=initOpt_cooling,
    p_start_in=p_start_in_cooling,
    p_start_out=p_start_out_cooling,
    p_start=p_start_cooling,
    initFromEnthalpy=initFromEnthalpy_cooling,
    h_start_in=h_start_in_cooling,
    h_start_out=h_start_out_cooling,
    h_start=h_start_cooling,
    T_start_in=T_start_in_cooling,
    T_start_out=T_start_out_cooling,
    T_start=T_start_cooling,
    X_start=X_start_cooling,
    positiveFlow=positiveFlow_cooling,
    from_dp=from_dp_cooling,
    dp_smooth=dp_smooth_cooling,
    mflow_smooth=mflow_smooth_cooling,
    generateEventForReversal=generateEventForReversal_cooling,
    CF_PressureLoss=CF_PressureLoss_cooling,
    CF_HeatTransfer=CF_HeatTransfer_cooling,
    m_flow_start=m_flow_start_cooling,
    NA=1,
    NB=1) "cooling pipe model"
    annotation (Placement(transformation(extent={{10,-60},{30,-40}})));

  extends .FuelCell.Stacks.Interfaces.LumpedCoolingGeometry;
  extends .FuelCell.Stacks.Interfaces.CoolingInitParameters(final ni_cooling=N, X_start_cooling=Medium_cooling.reference_X);
  extends .FuelCell.Stacks.Interfaces.AdvancedCoolingParameters;

  /* Initialization */
  extends .FuelCell.Stacks.Interfaces.DistributedAnodeInitialization(
     final n=N, X_start_anode=Medium_an.reference_X);
  extends .FuelCell.Stacks.Interfaces.DistributedCathodeInitialization(
     final n=N, X_start_cathode=Medium_cath.reference_X);

  /* Advanced parameters */
  extends .FuelCell.Stacks.Interfaces.AdvancedParameters(final
      useHeatTransfer_anode=false, final useHeatTransfer_cathode=false);

  SubStack subStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    N=N,
    A_cell=A_cell,
    n_cell=n_cell,
    Cp_cell=Cp_cell,
    kc=kc_cell,
    M_cell=M_stack,
    initOpt_anode=initOpt_anode,
    p_start_in_anode=p_start_in_anode,
    p_start_out_anode=p_start_out_anode,
    p_start_anode=p_start_anode,
    initFromEnthalpy_anode=initFromEnthalpy_anode,
    h_start_in_anode=h_start_in_anode,
    h_start_out_anode=h_start_out_anode,
    h_start_anode=h_start_anode,
    T_start_in_anode=T_start_in_anode,
    T_start_out_anode=T_start_out_anode,
    T_start_anode=T_start_anode,
    X_start_anode=X_start_anode,
    m_flow_start_anode=m_flow_start_anode,
    initOpt_cathode=initOpt_cathode,
    p_start_in_cathode=p_start_in_cathode,
    p_start_out_cathode=p_start_out_cathode,
    p_start_cathode=p_start_cathode,
    initFromEnthalpy_cathode=initFromEnthalpy_cathode,
    h_start_in_cathode=h_start_in_cathode,
    h_start_out_cathode=h_start_out_cathode,
    h_start_cathode=h_start_cathode,
    T_start_in_cathode=T_start_in_cathode,
    T_start_out_cathode=T_start_out_cathode,
    T_start_cathode=T_start_cathode,
    X_start_cathode=X_start_cathode,
    m_flow_start_cathode=m_flow_start_cathode)
                   annotation (Placement(transformation(extent={{-20,-22},{20,18}},
                  rotation=0)));

  .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p
    annotation (Placement(transformation(extent={{-10,86},{10,106}}),
        iconTransformation(extent={{30,56},{50,76}})));
  .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n
    annotation (Placement(transformation(extent={{-10,-106},{10,-86}}),
        iconTransformation(extent={{-50,56},{-30,76}})));

  .FuelCell.Interfaces.WaterVolumePort feed_cooling(redeclare package Medium =
        Medium_cooling)
    annotation (Placement(transformation(extent={{-110,-60},{-90,-40}}),
        iconTransformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-40,-66})));
  .FuelCell.Interfaces.WaterFlowPort drain_cooling(redeclare package Medium =
        Medium_cooling)
    annotation (Placement(transformation(extent={{90,-60},{110,-40}}),
        iconTransformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={40,-66})));
  .FuelCell.Interfaces.GasVolumePort
                      feed_an(redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-110,60},{-90,80}}),
        iconTransformation(extent={{-100,30},{-80,50}})));
  .FuelCell.Interfaces.GasVolumePort
                      feed_cath(redeclare replaceable package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{-110,-20},{-90,0}}),
        iconTransformation(extent={{-100,-50},{-80,-30}})));
  .FuelCell.Interfaces.GasFlowPort
                        drain_an(redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{90,60},{110,80}}),
        iconTransformation(extent={{80,30},{100,50}})));
  .FuelCell.Interfaces.GasFlowPort
                        drain_cath(redeclare replaceable package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{90,-20},{110,0}}),
        iconTransformation(extent={{80,-50},{100,-30}})));

  parameter Real CF_PressureLoss_cooling=1.0
    "Calibration factor for pressure drop in cooling pipe"
    annotation (Dialog(tab="Cooling", group="Calibration factors"));
  parameter Real CF_HeatTransfer_cooling=1.0
    "Calibration factor for heat transfer in cooling pipe"
    annotation (Dialog(tab="Cooling", group="Calibration factors"));
equation
  connect(feed_cath, subStack.feed_cathode) annotation (Line(
      points={{-100,-10},{-18,-10}},
      color={255,128,0},
      smooth=Smooth.None));

  connect(subStack.wall, coolingPipe.q) annotation (Line(
      points={{8,-2},{26,-2},{26,-36},{20,-36},{20,-45}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(feed_cooling, coolingPipe.portA[1]) annotation (Line(
      points={{-100,-50},{10,-50}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(coolingPipe.portB[1], drain_cooling) annotation (Line(
      points={{30,-50},{100,-50}},
      color={0,190,0},
      smooth=Smooth.None));
  connect(subStack.drain_cathode, drain_cath) annotation (Line(
      points={{18,-10},{100,-10}},
      color={255,128,0},
      smooth=Smooth.None));

  connect(subStack.drain_anode, drain_an) annotation (Line(
      points={{18,6},{80,6},{80,70},{100,70}},
      color={255,128,0},
      smooth=Smooth.None));

  connect(feed_an, subStack.feed_anode) annotation (Line(
      points={{-100,70},{-80,70},{-80,6},{-18,6}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(pin_p, subStack.pin_p) annotation (Line(
      points={{0,96},{0,11.2}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(pin_n, subStack.pin_n) annotation (Line(
      points={{0,-96},{0,-15.2}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-80},{100,80}}),
                                      graphics={
        Rectangle(
          extent={{-72,24},{72,-24}},
          lineColor={215,215,215},
          fillColor={175,175,175},
          fillPattern=FillPattern.CrossDiag),
        Rectangle(
          extent={{-80,56},{80,24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={215,215,215}),
        Rectangle(
          extent={{-80,-24},{80,-56}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={215,215,215})}),
    Documentation(info="<html>
<h4>CoolStack</h4>
<p>This template is suitable for easy setup of complete low-temperature fuel cell stacks with water cooling. The model contains the following components:</p>
<p><ul>
<li>TwoPhase water cooling pipe</li>
<li>Feed and drain connectors of reaction gas type (ideal gas) for both anode and cathode</li>
<li>Feed and drain connectors of twophase water type for the cooling channel</li>
<li>Electrical pins (positive and negative)</li>
<li>Replaceable substack model (used as placeholder) - needs to be replaced before simulation</li>
</ul></p>
<p>The model supports both co-flow and counter-flow configurations. The default setup is co-flow, to obtain counter-flow the flow in the cathode channel is reversed.</p>

<h4>Parametrization</h4>
<p>The discretization of the model is defined in the following way:</p>
<p><ul>
<li>The parameter <code>N</code>, which describes the number of discretizations in the flow directions (the membrane is divided in <code>N</code> segments) </li>
<li>The parameter <code>n_cell</code>, which describes the number of cells in the substack; the substack can be used to describe a single cell, a lumped subset of a stack or a complete stack</li>
</ul></p>

<h4>Substack Model</h4>
<p>The Stack template uses a replaceable model (<a href=\"modelica://FuelCell.Stacks.Templates.SubStack\">FuelCell.Stacks.Templates.SubStack</a>) that needs to be replaced with a finalized and parameterized substack model before use. Any preparameterized substack model based on the SubStack template may be used, see for instance the predefined <a href=\"modelica://FuelCell.Stacks.PEMFC\">PEMFC</a> substack. 
The number of cells modelled by the substack is described by the parameter <code>n_cell</code> and the number discrete segments in the fluid flow direction of each substack is the described by the parameter <code>N</code>.</p>
<h4>Cooling Pipe</h4>
<p>TwoPhase water cooling pipe discretized in <code>N</code> volumes in the fluid flow direction (<a href=\"modelica://Modelon.ThermoFluid.FlowChannels.DistributedTwoPhase\">Modelon.ThermoFluid.FlowChannels.DistributedTwoPhase</a>). Friction and heat transfer models in the pipe are partial and needs to be redeclared before use.</p>
<p></p>
</html>", revisions="<html>
Copyright &copy; 2004-2026, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end CoolStack_1;
