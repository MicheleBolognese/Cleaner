within SAMU.HollowFiberTest2_EB_inline;

model FlowChannel
  "Pipe model with distributed volumes and pressure drop, without reaction model"
   extends .FuelCell.Icons.Pipe;

//    replaceable package Medium =
//        FuelCell.Media.PreDefined.IdealGases.NASAMoistAir
//              constrainedby FuelCell.Media.Templates.ReactionGas   "Medium" annotation(choicesAllMatching);
  replaceable package Medium =
      .FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir                          constrainedby
    .Modelon.Media.Templates.IdealGasWithCondensingLiquid                                                                                                                                                                                                         "Medium" annotation(choicesAllMatching);

  .Modelon.ThermoFluid.FlowChannels.Records.PipeSummary
                                           summary(
    dp=channel.summary.dp,
    m_flow=channel.summary.m_flow,
    T_in=channel.summary.T_in,
    T_out=channel.summary.T_out,
    h_in=channel.summary.h_in,
    h_out=channel.summary.h_out,
    d_in=channel.summary.d_in,
    M=channel.summary.M,
    V=channel.summary.V)
    annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
  .FuelCell.Interfaces.CondensingGasVolumePort        portA(redeclare package Medium =
        Medium)
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}},
          rotation=0), iconTransformation(extent={{-100,-10},{-80,10}})));
  .FuelCell.Interfaces.CondensingGasFlowPort        portB(redeclare package Medium =
        Medium)
    annotation (Placement(transformation(extent={{80,-10},{100,10}},
          rotation=0), iconTransformation(extent={{80,-10},{100,10}})));

  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.DistributedGeometry;

  replaceable model Friction =
      .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss
  constrainedby
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance
    "Friction model"   annotation(choicesAllMatching, Dialog(group="Characteristic"));

  /* Heat Transfer*/
  parameter Boolean useHeatTransfer=false
    "Consider heat transfer effects (not used when connected to membrane)"                                         annotation(Evaluate=true, Dialog(group="Heat transfer"));
  replaceable model HeatTransfer =
      .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient
    constrainedby
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase
    "Definition of heat transfer coefficient" annotation(choicesAllMatching,Dialog(enable = useHeatTransfer,group="Heat transfer"));

  /*Additional parameters*/
  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.DistributedInitParameters(
    final ni=n,
    X_start = Medium.reference_X,
    C_start = fill(0.0, Medium.nC));
  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.AdvancedParameters;
  parameter Real CF_PressureLoss=1.0 "Calibration factor for pressure drop"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  parameter Real CF_HeatTransfer=1.0 "Calibration factor for heat transfer"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));

  parameter Boolean includeStaticHead = false "Consider static head"
    annotation(Evaluate=true,Dialog(group="Static head"));
  parameter .Modelica.Units.SI.Length height=0 "Height(outlet) - height(inlet)"
    annotation (Dialog(enable=includeStaticHead, group="Static head"));

  parameter .Modelica.Units.SI.Acceleration g=.Modelica.Constants.g_n "Acceleration of gravity"
    annotation (Dialog(enable=includeStaticHead, group="Static head"));

  .Modelica.Units.SI.MassFlowRate mX_flow[n,Medium.nS]=zeros(n, Medium.nS) "Mass transfer rate"
    annotation (Dialog(group="Mass transfer"));

  .Modelica.Units.SI.HeatFlowRate H_flow[n]=zeros(n) "Diffusion heat flow rate";

  .Modelica.Units.SI.Power checkEnergyBalance
    "Energy balance check for the channel, should be zero";

  .SAMU.HollowFiberTest2_EB_inline.DistributedChannel channel(
    redeclare .FuelCell.Interfaces.CondensingGasVolumePort portA,
    redeclare .FuelCell.Interfaces.CondensingGasFlowPort portB,
    redeclare package Medium = Medium,
    n_channels=n_channels,
    L=L,
    Dhyd=Dhyd,
    A=A,
    A_heat=A_heat,
    initOpt=initOpt,
    p_start_in=p_start_in,
    p_start_out=p_start_out,
    p_start=p_start,
    initFromEnthalpy=initFromEnthalpy,
    h_start_in=h_start_in,
    h_start_out=h_start_out,
    h_start=h_start,
    T_start_in=T_start_in,
    T_start_out=T_start_out,
    T_start=T_start,
    X_start=X_start,
    C_start=C_start,
    m_flow_start=m_flow_start,
    positiveFlow=positiveFlow,
    from_dp=from_dp,
    dp_smooth=dp_smooth,
    mflow_smooth=mflow_smooth,
    generateEventForReversal=generateEventForReversal,
    redeclare model Friction = Friction,
    redeclare model HeatTransfer = HeatTransfer,
    g=g,
    levels=(1:n)*height/n,
    includeStaticHead=includeStaticHead,
    CF_PressureLoss=CF_PressureLoss,
    CF_HeatTransfer=CF_HeatTransfer,
    rMX = mX_flow,
    Q_wall=if useHeatTransfer then (H_flow + wall.Q_flow) else H_flow,
    NA=1,
    NB=1)
    annotation (Placement(transformation(extent={{-20,-20},{20,20}})));

  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] wall(each T(nominal=500))
                                           annotation (Placement(transformation(
          extent={{20,40},{60,80}}, rotation=0), iconTransformation(extent={{30,40},
            {50,60}})));


equation
   assert(sum(L) >= abs(height), "The parameter L must be >= abs(height).");

  if useHeatTransfer then
    connect(channel.q, wall) annotation (Line(
      points={{0,10},{0,60},{40,60}},
      color={191,0,0},
      thickness=0.5));
  else
    connect(channel.q_fluid,wall);

    checkEnergyBalance = sum(channel.checkEnergyBalance);
  end if;

  connect(portA, channel.portA[1])
    annotation (Line(points={{-90,0},{-20,0}}, color={209,60,0}));
  connect(channel.portB[1], portB)
    annotation (Line(points={{20,0},{90,0}}, color={209,60,0}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),  Icon(graphics={
        Text(
          extent={{-100,-62},{100,-82}},
          lineColor={0,0,0},
          textString="%name"),     Line(
          points={{-64,-52},{56,-52}},
          color={0,0,0},
          thickness=0.5), Polygon(
          points={{68,-52},{56,-48},{56,-56},{68,-52}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<html>
<h4>FlowChannel</h4>
<p>Model of a pipe for discretized humidifier with <code>n</code> discrete volumes arranged in the fluid flow direction.  All control volumes have a unique pressure, temperature and mass content. 
Mass transfer and  diffusion heat flow rate are also accounted for using the variables <code>mX_flow</code> and <code>H_flow</code>.</p>
<h4>Connection Principles</h4>
<p>The first control volume is connected to <code>portA</code>, a flow model is connected to <code>portB</code>. To assure a numerically sound model of alternating volume and flow models, connect a flow model to <code>portA</code> and a volume model to <code>portB</code>.</p>
<h4>Assumptions</h4>
<p><ul>
<li>The total pipe length must be specified. All control volumes are assumed to have the same length, i.e. <code>L/n</code>.</li>
<li>The cross sectional area is the same for all control volumes. The volume is calculated from length and cross sectional area.</li>
<li>The mass flow rate into the first control volume is given by the connector <code>portA.</code> The mass flow rate between each control volume and from the last control volume to <code>portB</code> is calculated in a replaceable model <code>friction</code>, from the pressure drop between the two components.</li>
<li>Heat flow rate into each control volume from the heat connector is calculated. The heat transfer coefficient is calculated in the replaceable model <code>htcoeff</code>.</li>
</ul></p>
</html>", revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end FlowChannel;
