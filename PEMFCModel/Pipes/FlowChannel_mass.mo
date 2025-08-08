within PEMFCModel.Pipes;
model FlowChannel_mass "Pipe model with distributed volumes and pressure drop"
  
  extends .FuelCell.Icons.Pipe;

  replaceable package Medium = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir
    constrainedby .FuelCell.Media.Templates.ReactionGas "Medium" annotation(choicesAllMatching);

  .Modelon.ThermoFluid.FlowChannels.Records.PipeSummary summary(
    dp = channel.summary.dp,
    m_flow = channel.summary.m_flow,
    T_in = channel.summary.T_in,
    T_out = channel.summary.T_out,
    h_in = channel.summary.h_in,
    h_out = channel.summary.h_out,
    d_in = channel.summary.d_in,
    M = channel.summary.M,
    V = channel.summary.V) annotation (Placement(transformation(extent={{-9.608678980103516,-59.60867898010352},{10.391321019896484,-39.60867898010352}},rotation = 0.0,origin = {0.0,0.0})));
    
  .FuelCell.Interfaces.GasVolumePort portA(redeclare package Medium = Medium) annotation (Placement(transformation(extent={{-100,-10},{-80,10}},
          rotation=0), iconTransformation(extent={{-100,-10},{-80,10}})));
  .FuelCell.Interfaces.GasFlowPort portB(redeclare package Medium = Medium) annotation (Placement(transformation(extent={{80,-10},{100,10}},
          rotation=0), iconTransformation(extent={{80,-10},{100,10}})));

  parameter Integer n(min = 1) = 1 "Number of control volumes" annotation (Evaluate=true,Dialog(group="Geometry"));
  parameter Boolean reaction_occurrence = false "If true, reactions occur in the channel: choose your reaction model" annotation(Dialog(group="Reaction"));

  extends PEMFCModel.Pipes.Interfaces.LumpedGeometry;

  /* Friction */
  replaceable model Friction = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss
    constrainedby .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Friction model" annotation (choicesAllMatching, Dialog(group="Flow resistance"));

  /* Heat Transfer */
  parameter Boolean useHeatTransfer = false "If true, heat transfer between gas and channel wall (channel.q) considered (disabled when connected to a membrane)" annotation (Evaluate, Dialog(group="Heat transfer"));
  replaceable model HeatTransfer = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient(alpha = 0) constrainedby 
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Definition of heat transfer coefficient" annotation (choicesAllMatching,Dialog(enable = true, group="Heat transfer"));

  /* Reaction */
  replaceable model Reaction = PEMFCModel.Reaction.Templates.ZeroReaction constrainedby PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for channel" annotation (choicesAllMatching,Dialog(group="Reaction",enable = reaction_occurrence));
 

  /* Additional parameters */
  extends PEMFCModel.Pipes.Interfaces.DistributedChannelInitialization(
    final ni = n,
    X_start = Medium.reference_X,
    C_start = fill(0.0, Medium.nC));
    
  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.AdvancedParameters;
  parameter Real CF_PressureLoss = 1.0 "Calibration factor for pressure drop" annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  parameter Real CF_HeatTransfer = 1.0 "Calibration factor for heat transfer"annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
    
  input .Modelica.Units.SI.SpecificEnthalpy h_inflow "Specific enthalpy at inlet";

  //package WaterMedium = .Modelon.Media.PreDefined.TwoPhase.WaterIF97 "Water medium model";
  parameter Real X_weight(min = 0, max = 1) = 0.05 "Weight factor for inlet mass-based composition in lumped case";

  //.Modelica.Units.SI.MoleFraction y_sat "Saturated water mole fraction at exit";

  parameter Boolean includeStaticHead = false " = true to inclue static head" annotation (Evaluate,Dialog(group="Static head"));
  parameter .Modelica.Units.SI.Length height = 0 "Height(outlet) - height(inlet)" annotation (Dialog(enable=includeStaticHead, group="Static head"));

  .PEMFCModel.Pipes.ReactionChannel channel(
    redeclare package Medium = Medium,
    reaction_occurrence = reaction_occurrence,
    n_channels = fill(n_channels, n),
    L = fill(L/n, n),
    Dhyd = fill(Dhyd, n),
    A = fill(A, n),
    A_heat = fill(A_heat/n, n),
    initOpt = initOpt,
    p_start_in = p_start_in,
    p_start_out = p_start_out,
    p_start = p_start,
    initFromEnthalpy = initFromEnthalpy,
    h_start_in = h_start_in,
    h_start_out = h_start_out,
    h_start = h_start,
    T_start_in = T_start_in,
    T_start_out = T_start_out,
    T_start = T_start,
    X_start = X_start,
    m_flow_start = m_flow_start,
    positiveFlow = positiveFlow,
    from_dp = from_dp,
    dp_smooth = dp_smooth,
    mflow_smooth = mflow_smooth,
    generateEventForReversal = generateEventForReversal,
    useHeatTransfer = useHeatTransfer,
    redeclare model Friction = Friction,
    redeclare model HeatTransfer = HeatTransfer,
    redeclare model Reaction = Reaction,
    levels = (1:n)*height/n,
    includeStaticHead = includeStaticHead,
    CF_PressureLoss = CF_PressureLoss,
    CF_HeatTransfer = CF_HeatTransfer,
    rMX = port.mX_flow,
    Q_extra = port.H_flow, //in place of Q_wall = wall.Q_flow + port.H_flow
    NA = 1,
    NB = 1)
    annotation (Placement(transformation(extent={{-20.0,-20.0},{20.0,20.0}},rotation = 0.0,origin = {0.0,0.0})));
        
  .FuelCell.Interfaces.MassPort_a port[n](redeclare package Medium = Medium) annotation (Placement(transformation(extent={{-60.253138083394525,28.0},{-20.253138083394525,68.0}}, rotation=0.0,origin = {0.0,0.0}),
        iconTransformation(extent={{-40,34},{-20,54}})));
  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a wall[n](each T(nominal = 500)) annotation (Placement(transformation(
          extent={{20.0,40.0},{60.0,80.0}}, rotation=0.0,origin = {0.0,0.0}), iconTransformation(extent={{30,40},
            {50,60}})));

  .Modelica.Units.SI.MassFlowRate checkMassBalance = abs(channel.reaction.dmdt -(portA.m_flow + portB.m_flow + sum(port.m_flow))); 
  .Modelica.Units.SI.Power checkEnergyBalance = abs(channel.dUdt - (portA.m_flow*h_inflow +  portB.m_flow*portB.h_outflow + sum(port.H_flow)) - sum(wall.Q_flow));
        
equation

  assert(L >= abs(height), "The parameter L must be >= abs(height).");

  // Mass port connections
  port.p = channel.p;
  port.h = channel.gas.h;
  port.X = if n == 1 then {X_weight*actualStream(portA.X_outflow) + (1-X_weight)*channel.gas[1].X} else channel.gas.X;
  //y_sat = min(1, WaterMedium.saturationPressure_TX(channel.T[end])/portB.p);
  
  if useHeatTransfer then
    
    connect(channel.q, wall) annotation (Line(
      points={{0,10},{0,60},{40,60}},
      color={191,0,0},
      thickness=0.5,
      smooth=Smooth.None));
        
  else
    
    connect(channel.q_fluid, wall);
        
  end if;

  connect(portA, channel.portA[1]) annotation (Line(
      points={{-90,0},{-20,0}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(channel.portB[1], portB) annotation (Line(
      points={{20,0},{90,0}},
      color={255,128,0},
      smooth=Smooth.None));

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
<h4>FlowChannel_mass</h4>
<p>Model of a pipe with <code>n</code> discrete volumes arranged in the fluid flow direction.  All control volumes have a unique pressure, temperature and mass content. The model also includes <code>n</code> component mass flow ports for mass transfer from pipe to cell membrane model.</p>
<h4>Reaction Model</h4>
<p>Any dynamic reaction model that inherits from the <a href=\"modelica://FuelCell.Reactions.Templates.DynamicReaction\">FuelCell.Reactions.Templates.DynamicReaction</a> can be specified by the user to give dynamic reactions along the flow direction on both the primary and the secondary side. The DynamicReaction model is only used as a constraining class and does not represent any reaction by itself.</p>
<p>A typical use case where reactions takes place in a flow channel is internal reforming on the anode channel of a solid oxide fuel cell. </p>
<h4>Connection Principles</h4>
<p>The first control volume is connected to <code>portA</code>, a flow model is connected to <code>portB</code>. To assure a numerically sound model of alternating volume and flow models, connect a flow model to <code>portA</code> and a volume model to <code>portB</code>.</p>
<h4>Assumptions</h4>
<p><ul>
<li>The total pipe length must be specified. All control volumes are assumed to have the same length, i.e. <code>L/n</code>.</li>
<li>The cross sectional area is the same for all control volumes. The volume is calculated from length and cross sectional area.</li>
<li>The mass flow rate into the first control volume is given by the connector <code>portA.</code> The mass flow rate between each control volume and from the last control volume to <code>portB</code> is calculated in a replaceable model <code>friction</code>, from the pressure drop between the two components.</li>
<li>Heat flow rate into each control volume from the heat connector is calculated. The heat transfer coefficient is calculated in the replaceable model <code>htcoeff</code>.</li>
<li>Dynamic reactions possible in the flow direction.</li>
</ul></p>
<p></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end FlowChannel_mass;
