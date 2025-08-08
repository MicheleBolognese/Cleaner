within PEMFCModel.Stacks.Templates;
model CoolStack "Substack with cooling pipe"

  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformateShort constrainedby 
    .FuelCell.Media.Templates.ReactionGas "Anodic medium" annotation (choicesAllMatching,Dialog(enable = enable_setting));
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium" annotation (choicesAllMatching,Dialog(enable = enable_setting));
  replaceable package Medium_cooling = .FuelCell.Media.PreDefined.TwoPhase.WaterIF97 "Cooling medium" annotation (choicesAllMatching,Dialog(enable = enable_setting));
  
  parameter Boolean enable_setting = true "If true, selecatbility of parameters enabled";  
  
  .FuelCell.SubComponents.ComponentSummaries.CoolStackSummary summary(
    N = N,
    n_cell = n_cell,
    V_stack = subStack.summary.V_stack,
    P_stack = subStack.summary.P_stack,
    j_external = subStack.summary.j_external,
    T_stack = subStack.summary.T_stack,
    Q_stack = subStack.summary.Q_stack,
    anode_stoich = subStack.summary.anode_stoich,
    cathode_stoich = subStack.summary.cathode_stoich,
    m_flow_an = subStack.summary.m_flow_an,
    m_flow_cath = subStack.summary.m_flow_cath,
    T_cool = coolingPipe.T,
    m_flow_cool = -coolingPipe.portB[1].m_flow,
    dp_anode = subStack.feed_anode.p-subStack.drain_anode.p,
    dp_cathode = subStack.feed_cathode.p-subStack.drain_cathode.p,
    dp_cooling = coolingPipe.portA[1].p-coolingPipe.portB[1].p) annotation (Placement(transformation(extent={{-76,74},{-56,94}})));
  
  replaceable model SubStack = PEMFCModel.Stacks.Templates.SubStack(enable_setting = true) constrainedby PEMFCModel.Stacks.Templates.SubStack(enable_setting = false) "Substack model to use" annotation (choicesAllMatching,Dialog(tab = "Substack"));
  
  /* Substack parameters */  
  parameter Integer N(min = 1) = 4 "Number of along-the-channel discretization nodes" annotation (Dialog(tab = "Substack",enable = enable_setting));
  parameter Integer n_cell(min = 1) = 50 "Number of series-connected cells in the stack" annotation(Dialog(tab = "Substack",enable = enable_setting));
  parameter .Modelica.Units.SI.Area A_cell = 180e-4 "Cell active area" annotation (Dialog(tab = "Substack",enable = enable_setting));
  parameter .Modelica.Units.SI.SpecificHeatCapacity c_stack = 500 "Specific heat capacity of stack material" annotation (Dialog(tab = "Substack",enable = enable_setting));
  parameter .Modelica.Units.SI.Mass M_stack = 30 "Stack mass" annotation (Dialog(tab = "Substack",enable = enable_setting));
  parameter Real z = 20e-6 "Membrane thickness" annotation (Dialog(tab = "Substack", group = "Membrane characteristics", enable = enable_setting));
  parameter Real rho_dry_m = 1980 "Dry membrane density" annotation (Dialog(tab = "Substack", group = "Membrane characteristics", enable = enable_setting));
  parameter Real EW_m = 1.1 "Membrane equivalent weight" annotation (Dialog(tab = "Substack", group = "Membrane characteristics", enable = enable_setting));
  parameter String diffusiveSpecies[:] = fill("",0) "Species diffusing through membrane" annotation (Dialog(tab = "Substack", group = "Membrane characteristics", enable = enable_setting));
  parameter Real E0_ref = 1.229 "Single cell Nernst's potential at standard conditions" annotation (Dialog(tab = "Substack", group = "Characteristics and polarization", enable = enable_setting));
  parameter String contaminantsSpecies[:] = fill("",0) "Species affecting voltage loss due to contaminants" annotation (Dialog(tab = "Substack", group = "Characteristics and polarization", enable = enable_setting));       
  parameter Boolean includeCellConduction = false "If true, along-the-channel thermal conduction included" annotation(Dialog(tab = "Substack",group = "Along-the-channel thermal conduction",enable = enable_setting));
  parameter .Modelica.Units.SI.ThermalConductivity lambda_cell = 20 "Cell internal thermal conductivity" annotation (Dialog(enable = includeCellConduction, group = "Along-the-channel thermal conduction", tab = "Substack"));
  parameter .Modelica.Units.SI.Area A_crosssection_cell = 18e-3*2e-3 "Cross section area of a single cell" annotation (Dialog(enable = includeCellConduction, group = "Along-the-channel thermal conduction", tab = "Substack"));
  parameter .Modelica.Units.SI.Length length_cell = 0.15 "Cell length" annotation (Dialog(enable = includeCellConduction, tab = "Substack",group = "Along-the-channel thermal conduction"));
      
  /* Geometry */
  extends PEMFCModel.Stacks.Interfaces.LumpedAnodeGeometry(final enable_setting_anode_geometry = enable_setting);
  extends PEMFCModel.Stacks.Interfaces.LumpedCathodeGeometry(final enable_setting_cathode_geometry = enable_setting);
  extends PEMFCModel.Stacks.Interfaces.LumpedCoolingGeometry(final enable_setting_cooling_geometry = enable_setting);
  
  /* Initialization */
  extends PEMFCModel.Stacks.Interfaces.DistributedAnodeInitialization(final enable_setting_anode_init = enable_setting, final n = N, X_start_anode = Medium_an.reference_X);
  extends PEMFCModel.Stacks.Interfaces.DistributedCathodeInitialization(final enable_setting_cathode_init = enable_setting, final n = N, X_start_cathode = Medium_cath.reference_X);
  extends PEMFCModel.Stacks.Interfaces.CoolingInitialization(final enable_setting_cooling_init = enable_setting, final ni_cooling = N, X_start_cooling = Medium_cooling.reference_X);

  /* Advanced parameters */
  extends PEMFCModel.Stacks.Interfaces.AdvancedStackParameters(final enable_setting_advanced = enable_setting);
  extends PEMFCModel.Stacks.Interfaces.AdvancedCoolingParameters(final enable_setting_cooling_advanced = enable_setting);
  
  parameter Boolean reaction_occurrence_an = false "If true, reactions occur in the anode channel" annotation (Dialog(tab = "Substack", group = "Reactions", enable = enable_setting));  
  replaceable model Reaction_anode = PEMFCModel.Reaction.Templates.ZeroReaction constrainedby
    PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for anode channel" annotation(choicesAllMatching, Dialog(enable = reaction_occurrence_an, tab = "Substack", group = "Reactions"));
  parameter Boolean reaction_occurrence_cath = false "If true, reactions occur in the cathode channel" annotation (Dialog(tab = "Substack", group = "Reactions", enable = enable_setting));  
  replaceable model Reaction_cathode = PEMFCModel.Reaction.Templates.ZeroReaction constrainedby
    PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for cathode channel" annotation(choicesAllMatching, Dialog(enable = reaction_occurrence_cath, tab = "Substack", group = "Reactions"));
   
  replaceable model Friction_anode = Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss constrainedby 
    Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Friction model for anode channel flow" annotation(Dialog(tab = "Correlations",group = "Flow resistance"));
  replaceable model Friction_cathode = Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss constrainedby 
    Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Friction model for cathode channel flow" annotation(Dialog(tab = "Correlations",group = "Flow resistance"));
  replaceable model Friction_cooling = Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction constrainedby
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.Interfaces.TwoPhaseResistance "Friction model for cooling channel flow" annotation(Dialog(tab = "Correlations",group = "Flow resistance")); 
    
  replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient constrainedby 
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Heat transfer coefficient for anode channel flow" annotation(Dialog(tab = "Correlations",group = "Heat transfer"));
  replaceable model HeatTransfer_cathode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient constrainedby
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Heat transfer coefficient for cathode channel flow" annotation(Dialog(tab = "Correlations",group = "Heat transfer"));
  replaceable model HeatTransfer_cooling = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.CoefficientsPerPhase constrainedby
     .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.Interfaces.TwoPhaseHeatTransfer "Heat transfer coefficient for cooling channel flow" annotation(Dialog(tab = "Correlations",group = "Heat transfer"));

  parameter Real CF_AnodeSidePressureLoss = 1.0 "Calibration factor for pressure drop in anode channel" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_AnodeSideHeatTransfer = 1.0 "Calibration factor for heat transfer in anode channel" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CathodeSidePressureLoss = 1.0 "Calibration factor for pressure drop in cathode channel" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CathodeSideHeatTransfer = 1.0 "Calibration factor for heat transfer in cathode channel" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CoolingPressureLoss = 1.0 "Calibration factor for pressure drop in cooling pipe" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
  parameter Real CF_CoolingHeatTransfer = 1.0 "Calibration factor for heat transfer in cooling pipe" annotation (Dialog(tab = "Correlations", group = "Calibration factors",enable = enable_setting));
      
  input .Modelica.Units.SI.MassFraction X_feed_an[Medium_an.nS] "Mass-based composition at anode inlet";
  input .Modelica.Units.SI.MassFraction X_feed_cath[Medium_cath.nS] "Mass-based composition at cathode inlet";
  input .Modelica.Units.SI.SpecificEnthalpy h_inflow_an "Specific enthalpy at anode inlet";
  input .Modelica.Units.SI.SpecificEnthalpy h_inflow_cath "Specific enthalpy at cathode inlet";
  input .Modelica.Units.SI.SpecificEnthalpy h_inflow_cooling "Specific enthalpy at cooling pipe inlet";

  /* Cooling pipe*/
  PEMFCModel.Pipes.DistributedTwoPhase coolingPipe(
    redeclare package Medium = Medium_cooling,
    n_channels = fill(n_channels_cooling, N),
    L = fill(L_cooling/N, N),
    Dhyd = fill(Dhyd_cooling, N),
    A = fill(A_cooling, N),
    levels = (1:N)*height_cooling/N,
    A_heat = fill(A_heat_cooling/N, N),
    includeStaticHead = includeStaticHead_cooling,
    g = g_cooling,
    initOpt = initOpt_cooling,
    p_start_in = p_start_in_cooling,
    p_start_out = p_start_out_cooling,
    p_start = p_start_cooling,
    initFromEnthalpy = initFromEnthalpy_cooling,
    h_start_in = h_start_in_cooling,
    h_start_out = h_start_out_cooling,
    h_start = h_start_cooling,
    T_start_in = T_start_in_cooling,
    T_start_out = T_start_out_cooling,
    T_start = T_start_cooling,
    X_start = X_start_cooling,
    positiveFlow = positiveFlow_cooling,
    from_dp = from_dp_cooling,
    dp_smooth = dp_smooth_cooling,
    mflow_smooth = mflow_smooth_cooling,
    generateEventForReversal = generateEventForReversal_cooling,
    CF_PressureLoss = CF_CoolingPressureLoss,
    CF_HeatTransfer = CF_CoolingHeatTransfer,
    m_flow_start = m_flow_start_cooling,
    NA = 1,
    NB = 1,
    redeclare model Friction = Friction_cooling,
    redeclare model HeatTransfer = HeatTransfer_cooling,
    CF_PressureLoss = CF_CoolingPressureLoss,
    CF_HeatTransfer = CF_CoolingHeatTransfer) "Cooling pipe model" annotation (Placement(transformation(extent={{16.485484525575046,-59.97103326538163},{36.48548452557504,-39.97103326538163}},rotation = 0.0,origin = {0.0,0.0})));

  SubStack subStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    N = N,
    A_cell = A_cell,
    n_cell = n_cell,
    c_stack = c_stack,
    M_stack = M_stack,
    initOpt_anode = initOpt_anode,
    p_start_in_anode = p_start_in_anode,
    p_start_out_anode = p_start_out_anode,
    p_start_anode = p_start_anode,
    initFromEnthalpy_anode = initFromEnthalpy_anode,
    h_start_in_anode = h_start_in_anode,
    h_start_out_anode = h_start_out_anode,
    h_start_anode = h_start_anode,
    T_start_in_anode = T_start_in_anode,
    T_start_out_anode = T_start_out_anode,
    T_start_anode = T_start_anode,
    X_start_anode = X_start_anode,
    m_flow_start_anode = m_flow_start_anode,
    initOpt_cathode = initOpt_cathode,
    p_start_in_cathode = p_start_in_cathode,
    p_start_out_cathode = p_start_out_cathode,
    p_start_cathode = p_start_cathode,
    initFromEnthalpy_cathode = initFromEnthalpy_cathode,
    h_start_in_cathode = h_start_in_cathode,
    h_start_out_cathode = h_start_out_cathode,
    h_start_cathode = h_start_cathode,
    T_start_in_cathode = T_start_in_cathode,
    T_start_out_cathode = T_start_out_cathode,
    T_start_cathode = T_start_cathode,
    X_start_cathode = X_start_cathode,
    m_flow_start_cathode = m_flow_start_cathode,
    h_inflow_an = h_inflow_an,
    h_inflow_cath = h_inflow_cath,
    X_feed_an = X_feed_an,
    X_feed_cath = X_feed_cath,
    n_channels_anode = n_channels_anode,
    L_anode = L_anode,
    D_anode = D_anode,
    n_channels_cathode = n_channels_cathode,
    L_cathode = L_cathode,
    D_cathode = D_cathode,
    dp_smooth = dp_smooth,
    mflow_smooth = mflow_smooth,
    from_dp_anode = from_dp_anode,
    positiveFlow_anode = positiveFlow_anode,
    generateEventForReversal_anode = generateEventForReversal_anode,
    from_dp_cathode = from_dp_cathode,
    positiveFlow_cathode = positiveFlow_cathode,
    generateEventForReversal_cathode = generateEventForReversal_cathode,
    redeclare model Friction_anode = Friction_anode,
    redeclare model Friction_cathode = Friction_cathode,
    redeclare model HeatTransfer_anode =  HeatTransfer_anode,
    redeclare model HeatTransfer_cathode =  HeatTransfer_cathode,
    reaction_occurrence_an = reaction_occurrence_an,
    redeclare model Reaction_anode = Reaction_anode,
    reaction_occurrence_cath = reaction_occurrence_cath,
    redeclare model Reaction_cathode = Reaction_cathode,
    CF_AnodeSidePressureLoss = CF_AnodeSidePressureLoss,
    CF_AnodeSideHeatTransfer = CF_AnodeSideHeatTransfer,
    CF_CathodeSidePressureLoss = CF_CathodeSidePressureLoss,
    CF_CathodeSideHeatTransfer = CF_CathodeSideHeatTransfer,
    z = z,
    rho_dry_m = rho_dry_m,
    EW_m = EW_m,
    diffusiveSpecies = diffusiveSpecies,
    E0_ref = E0_ref,
    contaminantsSpecies = contaminantsSpecies,includeCellConduction = includeCellConduction,lambda_cell = lambda_cell,A_crosssection_cell = A_crosssection_cell,length_cell = length_cell) annotation (Placement(transformation(extent={{-22.867725444177402,-17.773439003303366},{17.132274555822598,22.226560996696634}},
                  rotation=0.0,origin = {0.0,0.0})));


  .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-13.220921470205052,-102.7605202795218},{6.779078529794947,-82.7605202795218}},rotation = 0.0,origin = {0.0,0.0}),
        iconTransformation(extent={{-50,56},{-30,76}})));
    
  .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-12.87308097585588,85.98960418496083},{7.1269190241441205,105.98960418496083}},rotation = 0.0,origin = {0.0,0.0}),
        iconTransformation(extent={{30,56},{50,76}})));

  .FuelCell.Interfaces.WaterVolumePort feed_cooling(redeclare package Medium = Medium_cooling) annotation (Placement(transformation(extent={{-110,-60},{-90,-40}}),
        iconTransformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-40,-66})));
    
  .FuelCell.Interfaces.WaterFlowPort drain_cooling(redeclare package Medium = Medium_cooling) annotation (Placement(transformation(extent={{90,-60},{110,-40}}),
        iconTransformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={40,-66})));
    
  .FuelCell.Interfaces.GasVolumePort feed_an(redeclare replaceable package Medium = Medium_an) annotation (Placement(transformation(extent={{-110,60},{-90,80}}),
        iconTransformation(extent={{-100,30},{-80,50}})));
    
  .FuelCell.Interfaces.GasVolumePort feed_cath(redeclare replaceable package Medium = Medium_cath) annotation (Placement(transformation(extent={{-110,-20},{-90,0}}),
        iconTransformation(extent={{-100,-50},{-80,-30}})));
    
  .FuelCell.Interfaces.GasFlowPort drain_an(redeclare replaceable package Medium = Medium_an) annotation (Placement(transformation(extent={{90,60},{110,80}}),
        iconTransformation(extent={{80,30},{100,50}})));
    
  .FuelCell.Interfaces.GasFlowPort drain_cath(redeclare replaceable package Medium = Medium_cath) annotation (Placement(transformation(extent={{90,-20},{110,0}}),
        iconTransformation(extent={{80,-50},{100,-30}})));

equation
  
  connect(feed_an, subStack.feed_anode) annotation (Line(
      points={{-100,70},{-80,70},{-80,10.226560996696634},{-20.867725444177402,10.226560996696634}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(subStack.drain_anode, drain_an) annotation (Line(
      points={{15.1322745558226,10.226560996696634},{80,10.226560996696634},{80,70},{100,70}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(feed_cath, subStack.feed_cathode) annotation (Line(
      points={{-100,-10},{-20.867725444177402,-10},{-20.867725444177402,-5.773439003303366}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(subStack.drain_cathode, drain_cath) annotation (Line(
      points={{15.1322745558226,-5.773439003303366},{15.1322745558226,-10},{100,-10}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(subStack.wall, coolingPipe.q) annotation (Line(
      points={{5.132274555822599,2.2265609966966338},{26.485484525575043,2.2265609966966338},{26.485484525575043,-44.97103326538163}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(feed_cooling, coolingPipe.portA[1]) annotation (Line(
      points={{-100,-50},{16.485484525575046,-50},{16.485484525575046,-49.97103326538163}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(coolingPipe.portB[1], drain_cooling) annotation (Line(
      points={{36.48548452557504,-49.97103326538163},{36.48548452557504,-50},{100,-50}},
      color={0,190,0},
      smooth=Smooth.None));
  connect(pin_n,subStack.pin_n) annotation(Line(points = {{-2.8730809758558795,95.98960418496083},{-2.8730809758558795,39.466630798178485},{-2.8677254441774007,39.466630798178485},{-2.8677254441774007,15.426560996696635}},color = {0,0,255}));
  connect(subStack.pin_p,pin_p) annotation(Line(points = {{-2.8677254441774007,-10.973439003303367},{-2.8677254441774007,-92.7605202795218},{-3.2209214702050533,-92.7605202795218}},color = {0,0,255}));
    
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
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end CoolStack;