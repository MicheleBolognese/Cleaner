within PEMFCModel.Stacks.Templates;
model SubStack
  
  "Substack template with membrane and single flow paths on anode and cathode sides"
  extends PEMFCModel.Stacks.Interfaces.PartialSubStack(
    redeclare replaceable model Membrane = PEMFCModel.Membrane.Templates.PartialMembrane, 
    cell(
      final h_conv_an = anode_channel.channel.alpha, 
      final h_conv_cath = cathode_channel.channel.alpha,
      final pstart = (p_start_in_anode + p_start_in_cathode)/2, 
      final Tstart = (T_start_in_anode + T_start_in_cathode)/2));

  /* Geometry parameters */
  extends PEMFCModel.Stacks.Interfaces.LumpedAnodeGeometry(final enable_setting_anode_geometry = enable_setting);
  extends PEMFCModel.Stacks.Interfaces.LumpedCathodeGeometry(final enable_setting_cathode_geometry = enable_setting);

  /* Initialization */
  extends PEMFCModel.Stacks.Interfaces.DistributedAnodeInitialization(final n = N, final enable_setting_anode_init = enable_setting, X_start_anode = Medium_an.reference_X);
  extends PEMFCModel.Stacks.Interfaces.DistributedCathodeInitialization(final n = N, final enable_setting_cathode_init = enable_setting, X_start_cathode = Medium_cath.reference_X);

  /* Advanced parameters */
  extends PEMFCModel.Stacks.Interfaces.AdvancedStackParameters(final enable_setting_advanced = enable_setting);
  parameter .Modelica.Units.SI.Density d0_prox = 400 "Nominal density in prox loss (only used if addProxToAnode = true)" annotation (Dialog(enable = addProxToAnode, tab = "Advanced", group = "Reactions"));
  parameter .Modelica.Units.SI.Pressure dp0_prox = 100 "Nominal pressure drop in prox loss (only used if addProxToAnode = true)" annotation (Dialog(enable = addProxToAnode, tab = "Advanced",group = "Reactions"));
  parameter .Modelica.Units.SI.MassFlowRate m_flow_nom_prox = 0.1 "Nominal mass flow rate in prox loss (only used if addProxToAnode = true)" annotation (Dialog(enable = addProxToAnode, tab = "Advanced", group = "Reactions"));

  replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss constrainedby
  .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Flow definition through anodic channel " annotation (choicesAllMatching,Dialog(tab="Correlations",group="Flow resistance",enable = enable_setting));
  
  replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient (alpha = 250) constrainedby
  .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Definition of heat transfer coefficient at anodic side" annotation(choicesAllMatching, Dialog(tab = "Correlations", group = "Heat transfer",enable = enable_setting));
  
  parameter Boolean reaction_occurrence_an = false "If true, reactions occur in the anode channel" annotation (Dialog(tab = "Advanced", group = "Reactions", enable = enable_setting));  
  replaceable model Reaction_anode = PEMFCModel.Reaction.Templates.ZeroReaction constrainedby
    PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for anode channel" annotation(choicesAllMatching, Dialog(enable = reaction_occurrence_an, tab = "Advanced", group = "Reactions"));
    
  replaceable model Friction_cathode =
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss constrainedby
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Flow definition through cathodic channel" annotation (choicesAllMatching, Dialog(tab = "Correlations", group = "Flow resistance",enable = enable_setting));
    
  replaceable model HeatTransfer_cathode =
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient (alpha = 250) constrainedby
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Definition of heat transfer coefficient at cathodic side" annotation(choicesAllMatching, Dialog(enable = true, tab = "Correlations", group = "Heat transfer",enable = enable_setting));
  
  parameter Boolean reaction_occurrence_cath = false "If true, reactions occur in the cathode channel" annotation (Dialog(tab = "Advanced", group = "Reactions", enable = enable_setting));  
  replaceable model Reaction_cathode = PEMFCModel.Reaction.Templates.ZeroReaction constrainedby
    PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for cathode channel" annotation(choicesAllMatching, Dialog(enable = reaction_occurrence_cath, tab = "Advanced", group = "Reactions"));

  parameter Real CF_AnodeSidePressureLoss = 1.0 "Calibration factor for pressure drop at anodic side" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_AnodeSideHeatTransfer = 1.0 "Calibration factor for heat transfer at anodic side" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CathodeSidePressureLoss = 1.0 "Calibration factor for pressure drop at cathodic side" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CathodeSideHeatTransfer = 1.0 "Calibration factor for heat transfer at cathodic side" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));

  input .Modelica.Units.SI.SpecificEnthalpy h_inflow_an "Specific enthalpy at anode inlet";
  input .Modelica.Units.SI.SpecificEnthalpy h_inflow_cath "Specific enthalpy at cathode inlet";
  .Modelica.Units.SI.MassFlowRate dmdt = anode_channel.channel.reaction.dmdt + cathode_channel.channel.reaction.dmdt;
  .Modelica.Units.SI.Power dUdt = anode_channel.channel.dUdt + cathode_channel.channel.dUdt + cell.dUdt; 
    
  .Modelica.Units.SI.MassFlowRate checkMassBalance = abs(dmdt - (feed_anode.m_flow + feed_cathode.m_flow + drain_anode.m_flow + drain_cathode.m_flow));
  .Modelica.Units.SI.Power checkEnergyBalance = abs(dUdt - (feed_anode.m_flow*h_inflow_an + feed_cathode.m_flow*h_inflow_cath + drain_anode.m_flow*drain_anode.h_outflow + drain_cathode.m_flow*drain_cathode.h_outflow) + summary.P_stack - sum(wall.Q_flow));
     
  .Modelica.Units.SI.MoleFraction H2_utilization "H2 utilization";
  
  .FuelCell.SubComponents.ComponentSummaries.SubstackSummary summary(N = N,
  n_cell = n_cell,
  V_stack = V_stack,
  P_stack = P_stack,
  j_external = cell.j,
  T_stack = T_stack,
  Q_stack = cell.Q_stack,
  anode_stoich = anode_stoich,
  cathode_stoich = cathode_stoich,
  m_flow_an = -anode_channel.portB.m_flow,
  m_flow_cath = -cathode_channel.portB.m_flow,
  dp_anode = feed_anode.p-drain_anode.p,
  dp_cathode = feed_cathode.p-drain_cathode.p) annotation (Placement(transformation(extent={{57.934607626155255,92.81678993997663},{77.93460762615526,112.81678993997663}},rotation = 0.0,origin = {0.0,0.0})));

  .PEMFCModel.Pipes.FlowChannel_mass anode_channel(
    redeclare package Medium = Medium_an,
    A = A_anode,
    A_heat = A_heat_anode*n_cell,
    initOpt = initOpt_anode,
    p_start_in = p_start_in_anode,
    p_start_out = p_start_out_anode,
    p_start = p_start_anode,
    initFromEnthalpy = initFromEnthalpy_anode,
    h_start_in = h_start_in_anode,
    h_start_out = h_start_out_anode,
    h_start = h_start_anode,
    T_start_in = T_start_in_anode,
    T_start_out = T_start_out_anode,
    T_start = T_start_anode,
    X_start = X_start_anode,
    m_flow_start = m_flow_start_anode,
    L = L_anode,
    Dhyd = Dhyd_anode,
    V = V_anode*n_cell,
    positiveFlow = positiveFlow_anode,
    from_dp = from_dp_anode,
    dp_smooth = dp_smooth,
    mflow_smooth = mflow_smooth,
    n_channels = n_channels_anode,
    final useHeatTransfer = false,
    generateEventForReversal = generateEventForReversal_anode,
    redeclare model Friction = Friction_anode,
    redeclare model HeatTransfer = HeatTransfer_anode,
    reaction_occurrence = reaction_occurrence_an,
    redeclare model Reaction =  Reaction_anode,
    CF_PressureLoss = CF_AnodeSidePressureLoss,
    CF_HeatTransfer = CF_AnodeSideHeatTransfer,
    n = N,
    D = D_anode,
    C = C_anode,
    h_inflow = h_inflow_an) annotation (Placement(transformation(extent={{-20.0,58.230125523012546},{20.0,18.230125523012546}},
          rotation=0.0,origin = {0.0,0.0})));

  .PEMFCModel.Pipes.FlowChannel_mass cathode_channel(
    redeclare package Medium = Medium_cath,
    A = A_cathode,
    A_heat = A_heat_cathode*n_cell,
    initOpt = initOpt_cathode,
    p_start_out = p_start_out_cathode,
    p_start = p_start_cathode,
    initFromEnthalpy = initFromEnthalpy_cathode,
    h_start_in = h_start_in_cathode,
    h_start_out = h_start_out_cathode,
    h_start = h_start_cathode,
    T_start_in = T_start_in_cathode,
    T_start_out = T_start_out_cathode,
    T_start = T_start_cathode,
    X_start = X_start_cathode,
    m_flow_start = m_flow_start_cathode,
    p_start_in = p_start_in_cathode,
    L = L_cathode,
    Dhyd = Dhyd_cathode,
    V = V_cathode*n_cell,
    positiveFlow = positiveFlow_cathode,
    from_dp = from_dp_cathode,
    dp_smooth = dp_smooth,
    mflow_smooth = mflow_smooth,
    n_channels = n_channels_cathode,
    final useHeatTransfer = false,
    generateEventForReversal = generateEventForReversal_cathode,
    redeclare model Friction = Friction_cathode,
    redeclare model HeatTransfer = HeatTransfer_cathode,
    reaction_occurrence = reaction_occurrence_cath,
    redeclare model Reaction =  Reaction_cathode,
    CF_PressureLoss = CF_CathodeSidePressureLoss,
    CF_HeatTransfer = CF_CathodeSideHeatTransfer,
    n = N,
    D = D_cathode,
    C = C_cathode,
    h_inflow = h_inflow_cath) annotation (Placement(transformation(extent={{-20.0,-60.0},{20.0,-20.0}},
          rotation=0.0,origin = {0.0,0.0})));

   .FuelCell.FlowResistances.TurbulentLoss proxLoss(
      redeclare package Medium = Medium_an,
      dp0 = dp0_prox,
      m_flow0 = m_flow_nom_prox,
      d0 = d0_prox) if addProxToAnode annotation (Placement(transformation(extent={{-52,34},{-40,46}})));
    
    .FuelCell.Reactors.PrOx.InvariantProx prox(
    redeclare package Medium = Medium_an,
    V_tot(displayUnit = "l") = 0.001,
    pstart = p_start_in_anode,
    Tstart = T_start_in_anode,
    Xstart = X_start_anode,
    m_flow_start = m_flow_nom_prox) if addProxToAnode annotation (Placement(transformation(extent={{-76,32},{-60,48}})));

protected
    
  parameter Integer[1] i_H2 = Medium_an.substanceIndexVector({"H2"}) "Index for H2 in the anodic medium" annotation (Evaluate = true);
  parameter Integer[1] i_O2 = Medium_cath.substanceIndexVector({"O2"}) "Index for O2 in the cathodic medium" annotation(Evaluate = true);
  .Modelica.Units.SI.MolarFlowRate anode_stoich_num "H2 molar flow rate at anode inlet";
  .Modelica.Units.SI.MolarFlowRate anode_stoich_den "Molar flow rate of reacting H2";
  .Modelica.Units.SI.MolarFlowRate cathode_stoich_num "O2 molar flow rate at cathode inlet";
  .Modelica.Units.SI.MolarFlowRate cathode_stoich_den "Molar flow rate of reacting O2";

initial equation
    
  assert(i_H2[1] > 0, "Substack: the specified anodic medium does not contain the required substance: H2", level = AssertionLevel.error);
  assert(i_O2[1] > 0, "Substack: The specified cathodic medium does not contain the required substance: O2", level = AssertionLevel.error);
    
equation
    
  T_stack = cell.T_cell;
  V_stack = cell.V_stack;
  P_stack = cell.P_stack;

  // Stoichiometry at anode and cathode side
  anode_stoich_num = max(0, anode_channel.portA.m_flow*X_feed_an[i_H2[1]]/Medium_an.MMX[i_H2[1]]);
  anode_stoich_den = max(.Modelica.Constants.eps, (anode_channel.portA.m_flow*X_feed_an[i_H2[1]]+anode_channel.portB.m_flow*anode_channel.portB.X_outflow[i_H2[1]])/Medium_an.MMX[i_H2[1]]);
  anode_stoich = if anode_stoich_num/anode_stoich_den <100 then anode_stoich_num/anode_stoich_den else 100;

  cathode_stoich_num = max(0, cathode_channel.portA.m_flow*X_feed_cath[i_O2[1]]/Medium_cath.MMX[i_O2[1]]);
  cathode_stoich_den = max(.Modelica.Constants.eps, (cathode_channel.portA.m_flow*X_feed_cath[i_O2[1]]+cathode_channel.portB.m_flow*cathode_channel.portB.X_outflow[i_O2[1]])/Medium_cath.MMX[i_O2[1]]);
  cathode_stoich = if cathode_stoich_num/cathode_stoich_den < 100 then cathode_stoich_num/cathode_stoich_den else 100;

  H2_utilization = min(1, pin_n.i/(2*.FuelCell.Internal.Units.F)*n_cell/anode_stoich_num);

  connect(cell.pin_p, pin_p) annotation (Line(
      points={{-18,10},{-30,10},{-30,66},{0,66}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(cell.pin_n, pin_n) annotation (Line(
      points={{-18,-10},{-30,-10},{-30,-70.5564850647429},{0,-70.5564850647429},{0,-66}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(anode_channel.port, cell.port_an) annotation (Line(
      points={{-6,29.430125523012546},{-6,8.8}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(cathode_channel.port, cell.port_cath) annotation (Line(
      points={{-6,-31.2},{-6,-23.2},{-6,-23.2},{-6,-8.8}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(anode_channel.wall, cell.wall_an) annotation (Line(
      points={{8,28.230125523012546},{8,10}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(cathode_channel.wall, cell.wall_cath) annotation (Line(
      points={{8,-30},{8,-24},{8,-24},{8,-10}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(cell.wall, wall) annotation (Line(
      points={{22,0},{40,0}},
      color={191,0,0},
      smooth=Smooth.None));
  if addProxToAnode then
    
    connect(feed_anode, prox.feed) annotation (Line(points={{-90,40},{-75.2,40}},
      color={255,128,0},
      smooth=Smooth.None));
    connect(prox.drain, proxLoss.feed) annotation (Line(
      points={{-60.8,40},{-51.4,40}},
      color={255,128,0},
      smooth=Smooth.None));
    connect(proxLoss.drain, anode_channel.portA)  annotation (Line(
      points={{-40.6,40},{-18,40}},
      color={255,128,0},
      smooth=Smooth.None));

  else
    
    connect(feed_anode,anode_channel.portA);
    
  end if;
  connect(anode_channel.portB, drain_anode) annotation (Line(
      points={{18,38.230125523012546},{90,38.230125523012546},{90,40}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(feed_cathode, cathode_channel.portA) annotation (Line(
      points={{-88,-40},{-18,-40}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(cathode_channel.portB, drain_cathode) annotation (Line(
      points={{18,-40},{90,-40}},
      color={255,128,0},
      smooth=Smooth.None));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
                       Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
    Documentation(info="<html>
<h4>SubStack</h4>
<p>SubStack template suitable for easy setup of PEMFC, SOFC or user defined substack configurations. The model extends from FuelCell.Stacks.Interfaces.PartialSubstack.
and can be used to describe either a single cell or an arbitrary number of cells. The model contains the following components:</p>
<p><ul>
<li>Feed and drain connectors of reaction gas type (ideal gas) for both anode and cathode</li>
<li>Electrical pins (positive and negative)</li>
<li>Replaceable partial cell membrane model (used as placeholder) - needs to be replaced before simulation</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Fluid flow channels for the cathode and the anode side</li>
<li>Conditional prox reactor before inlet to the anode channel for oxidation of excess fractions of carbon monoxide</li>
</ul></p>
<p>The model supports both co-flow and counter-flow configurations. The default setup is co-flow, to obtain counter-flow the flow on the cathode side is reversed.</p>

<h4>Parametrization</h4>
<p>The discretization of the model is defined in the following way:</p>
<p><ul>
<li>The parameter <code>N</code>, which describes the number of discretizations in the flow directions (the membrane is divided in <code>N</code> segments) </li>
<li>The parameter <code>n_cell</code>, which describes the number of cells in the substack; the substack can be used to describe a single cell, a lumped subset of a stack or a complete stack</li>
</ul></p>

<h4>Membrane model</h4>
<p>The membrane is an essential part of the fuel cell as it describes the fuel cell reactions and mass transfer between anode and cathode.
The SubStack template uses a replaceable partial cell membrane model (<a hi_H2=\"modelica://FuelCell.Membranes.Templates.PartialCell\">FuelCell.Membranes.Templates.PartialCell</a>) that needs to be replaced with a finalized and parameterized membrane model before use. Any membrane model based on the CellPartial membrane may be used, see for instance the predefined <a hi_H2=\"modelica://FuelCell.Membranes.PEMFC\">PEMFC</a> and <a hi_H2=\"modelica://FuelCell.Membranes.SOFC\">SOFC</a> membrane models.</p>
<h4>Flow channels</h4>
<p>The anode and cathode channels are FuelCell.Pipes.FlowChannel_mass pipes and consist of <code>N</code> discrete volumes in the fluid flow direction, each with unique pressure, temperature and mass content. Each volume has a mass flow port for mass flow to the membrane and a heat flow port for heat transfer between channel and membrane. Some modeling options for the channels are:</p>
<p><ul>
<li>Possibility to add intenal reactions to any of the pipes (internal i_H2orming in the case of SOFC)</li>
<li>Possibility to change the friction model that is used to calculate the pressure drop over a pipe</li>
<li>Possibility to change or disable the heat transfer model used for heat transfer between channel and membrane</li>
</ul></p>
<h4>Conditional prox reactor</h4>
<p>The conditional pi_H2erential oxidation (prox) reactor before the anode inlet is enabled or disabled using the parameter <code>addproxToAnode</code>. The prox reactor is only needed for cases where extra i_O2 is added just before the stack to oxidate any excess fraction of carbon monoxide (CO), before the stack inlet. This is done to prevent the degradation of PEMFC membranes due to too high fractions of carbon monoxide. A typical PEM fuel cell can handle up to 25 ppm CO.</p>
<p></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end SubStack;
