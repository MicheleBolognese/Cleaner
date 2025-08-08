within PEMFCModel.Pipes;
model ReactionChannel
  
  "Flow channel for components including reactions using the finite volume method"

  extends .FuelCell.Icons.Pipe;
  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.Channel(
    redeclare replaceable package Medium = .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong constrainedby .FuelCell.Media.Templates.ReactionGas,
    X_start = Medium.reference_X,
    summary(
      dp = pA-pB,
      m_flow = m_flow_mean,
      T_in = TA,
      d_in = Medium.density(state[1]),
      T_out = gas[end].T,
      h_in = gas[1].h,
      h_out = gas[end].h,
      M = sum(M),
      V = sum(V)),
      redeclare .FuelCell.Interfaces.GasVolumePort portA,
      redeclare .FuelCell.Interfaces.GasFlowPort portB,
      T(each stateSelect = StateSelect.prefer));
 
  parameter Boolean reaction_occurrence = false "If true, reactions occur in the channel: choose your reaction model" annotation(Dialog(group="Reaction"));
   
  /* Friction */
  replaceable model Friction = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss
    constrainedby .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance "Friction model" annotation(choicesAllMatching, Dialog(group="Flow resistance"));

  /* Heat Transfer */
  parameter Boolean useHeatTransfer = false " = true if heat transfer between gas and channel wall is enabled (disabled with a membrane)" annotation(Dialog(group = "Heat transfer"));
  replaceable model HeatTransfer = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient(alpha = 0)
    constrainedby .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase "Definition of heat transfer coefficient" annotation(choicesAllMatching, Dialog(group = "Heat transfer", enable = true));

  /* Reaction */
  replaceable model Reaction = .PEMFCModel.Reaction.Templates.ZeroReaction constrainedby 
    PEMFCModel.Reaction.Templates.DynamicReaction "Reaction model for channel" annotation (choicesAllMatching,Dialog(group="Reaction", enable = reaction_occurrence));

  /* Advanced parameters */
  parameter Boolean dp_asState = false "Use pressure difference as state variable" annotation (Evaluate=true, Dialog(tab="Advanced",group="Numerics"));
  parameter .Modelon.ThermoFluid.Choices.FrictionDistribution frictionDistribution = .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric "Discretization scheme for friction and control volumes (see info)" annotation (Evaluate=true, Dialog(tab="Advanced",group="Numerics"));
  final parameter Integer n_fric = if (frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
    frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol) then n
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then n+1
    else n-1 annotation (Evaluate=true);

  parameter Real CF_PressureLoss = 1.0 "Calibration factor for pressure drop" annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  parameter Real CF_HeatTransfer = 1.0 "Calibration factor for heat transfer" annotation (Dialog(tab="Advanced", group="Calibration factors"));

  parameter Boolean includeStaticHead = false " = true to include static head" annotation (Evaluate,Dialog(group="Static head"));
  parameter .Modelica.Units.SI.Length levels[n] = fill(0, n) "Relative levels of control volume outlets (inlet port at 0)" annotation (Dialog(enable=includeStaticHead, group="Static head"));
  parameter .Modelica.Units.SI.Acceleration g = .Modelica.Constants.g_n "Gravitational acceleration" annotation (Dialog(enable=includeStaticHead, group="Static head"));

  Friction friction[n_fric](
    each m_flow(start = m_flow_start),
    redeclare each package Medium = Medium,
    dp = dp_fric,
    stateA_inflow = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then state[2:n+1]
                    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then state[1:n]
                    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then state[1:n+1]
                    else state[2:n],
    stateB_inflow = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then state[3:n+2]
                    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then state[2:n+1]
                    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then state[2:n+2]
                    else state[3:n+1],
    final A = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(1, {sum(A)/n}, A)
              elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then A[1:n-1]
              else A,
    each final L = L_total,
    final Dhyd = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(1, {sum(Dhyd)/n}, Dhyd)
                 elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then Dhyd[1:n-1]
                 else Dhyd,
    final n_channels = n_channels_fric,
    final lengthFraction = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(1, {1/(n+1)}, L/sum(L)*(n/(n+1)))
                           elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then L[1:n-1]/sum(L[1:n-1])
                           else L/sum(L),
    each final from_dp = from_dp,
    each final positiveFlow = positiveFlow,
    final dp_smooth = dp_smooth * (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(1, {1/(n+1)}, L/sum(L)*(n/(n+1)))
                                   elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then L[1:n-1]/sum(L[1:n-1])
                                   else L/sum(L)),
    final mflow_smooth = mflow_smooth./n_channels_fric,
    each final F_user = CF_PressureLoss);

  HeatTransfer htcoeff[n](
    redeclare each package Medium = Medium,
    final A = A,
    final Dhyd = Dhyd,
    final L = L,
    each final F_user = CF_HeatTransfer,
    m_flow = m_flow[1:n]./n_channels,
    stateA = state[2:n+1],
    stateB = state[2:n+1]);

  parameter String volName = getInstanceName() "Volume-name for better diagnosis";

  Reaction reaction(
    N = n,
    redeclare package Medium = Medium,
    volName = volName,
    T = T,
    pstart = p_start[1],
    Tstart = T_start[1],
    Xout_start = X_start,
    V = V_tot,
    mX_flow = dMX + rMX,
    X_in = Xi_i[1:n,:],
    p = p) "Dynamic reaction model in channel"  annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
//     X_in=(if n == 1 then {Xi_i} else [transpose([Xi_i]); gas[1:n - 1].X]),

  .Modelon.ThermoFluid.Interfaces.TempHeatPort q_fluid[n] "Connector exposing fluid temperature" annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
        iconTransformation(extent={{-10,-10},{10,10}})));
    
  Medium.ReactionProperties gas[n](
    p(start = p_start[1:n]) = p,
    T(start = T_start[2:n+1]) = T,
    X(start = fill(X_start[1:nS], n)) = reaction.X_out_real[1:n, :]) "Gas nodal properties";

  .Modelon.Media.Units.AbsolutePressure p[n](start = p_start[1:n], each stateSelect = StateSelect.prefer) "Gas pressure in each control volume";
  .Modelica.Units.SI.Mass M[n] "Total mass in each control volume";
  .Modelica.Units.SI.Energy U[n] "Internal energy of each control volume";
  .Modelica.Units.SI.Mass MX[n, nS] "Mass of each component in each control volume";

  .Modelica.Units.SI.Pressure dp[n + 1] "Pressure difference between portA, control volumes and portB";
  .Modelica.Units.SI.Pressure sh[n] "Static head between control volumes";

  .Modelica.Units.SI.AbsolutePressure pA "Pressure in portA";
  .Modelica.Units.SI.AbsolutePressure pB "Pressure in portB";

  .Modelica.Units.SI.MassFlowRate m_flow_mean "Average mass flow rate in pipe, positive from portA to portB";

  Medium.ThermodynamicState state[n+2] "Thermodynamic state in portA, control volumes and portB";

  .Modelica.Units.SI.Temperature TA "Upstream temperature when flow A -> B";
  .Modelica.Units.SI.Temperature TB "Upstream temperature when flow B -> A";

  .Modelica.Units.SI.HeatFlowRate Q[n] "Heat flow rate into each control volume";
  .Modelica.Units.SI.HeatFlowRate Q_tot "Total heat flow rate into pipe";

  .Modelica.Units.SI.HeatFlowRate Q_wall[n] = if useHeatTransfer then q.Q_flow else fill(0, n) "Heat flow rate from connector q into each control volume";
  .Modelica.Units.SI.HeatFlowRate Q_fluid[n] "Heat flow rate from connector q_fluid into each control volume";
  .Modelica.Units.SI.HeatFlowRate Q_extra[n] = fill(0, n) "Extra energy per unit time into each control volume";
    
  .Modelica.Units.SI.CoefficientOfHeatTransfer alpha[n] "Heat transfer coefficient";

  .Modelica.Units.SI.MassFlowRate rMX[n, nS] = zeros(n, nS) "Mass residual of each component in each control volume for balance";
  .Modelica.Units.SI.Power dUdt = sum(dUdt_cv) "Time-derivative of channel internal energy";  

protected
    
  .Modelica.Units.SI.MassFlowRate mflow_A_in "Total inlet mass flow rate at portA";
  .Modelica.Units.SI.MassFlowRate mflow_B_in "Total inlet mass flow rate at portB";
  .Modelica.Units.SI.MassFlowRate m_flow_fric[n_fric] = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
     then m_flow[2:n + 1] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
     then m_flow[1:n] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
     then m_flow else m_flow[2:n] "Mass flow rates through control volume boundaries, positive from portA to portB";
  .Modelica.Units.SI.SpecificEnthalpy hA "Mixture specific enthalpy of inlet flows at portA";
  .Modelica.Units.SI.SpecificEnthalpy hB "Mixture specific enthalpy of inlet flows at portB";
  .Modelica.Units.SI.MassFraction XA[Medium.nS] "Mass-based composition of inlet flows at portA";
  .Modelica.Units.SI.MassFraction XB[Medium.nS] "Mass-based composition of inlet flows at portB";
  .Modelica.Units.SI.Power dUdt_cv[n] "Time-derivative of internal energy of each control volume";

  final parameter Real n_channels_fric[n_fric] = (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(1, {n_channels[1]}, n_channels)
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then n_channels[1:n-1]
    else n_channels);
    
  .Modelica.Units.SI.Pressure dp_fric[n_fric](each nominal = 0.1, stateSelect = if dp_asState then StateSelect.always else StateSelect.default) "Pressure drops from friction models, vector length depending on discretization option";

  .Modelica.Units.SI.Temperature T_wall[n] "Wall temperature, if heat transfer enabled";

  constant Integer nS = Medium.nS "Number of medium components";
  .Modelon.Media.Units.MassFlowRate m_flow[n + 1](each start = m_flow_start, each fixed = false) "Mass flow rates of fluid across control volume boundaries";
  .Modelon.Media.Units.MassFlowRate mX_flow[n + 1, Medium.nS] "Mass flow rates of each component across control volume boundaries";
  .Modelon.Media.Units.EnthalpyFlowRate H_flow[n + 1] "Enthalpy flow rates of fluid across control volume boundaries";
  .Modelon.Media.Units.MassFraction Xi_i[n + 1, Medium.nS] = [transpose([XA]); gas.X] "Mass-based composition at inlet of each control volume";
  .Modelica.Units.SI.EnergyFlowRate dH[n] "Net enthalpy flow rate through each control volume";
  .Modelica.Units.SI.MassFlowRate dMX[n, nS] "Net mass flow rate of each component through each control volume";
  .Modelica.Units.SI.MolarFlowRate dZX[n, Medium.nS] = reaction.Zx_flow + reaction.rZ "Net molar flow rate of each component through each control volume = Time-derivative of amount of moles of each component in each control volume";
  .Modelica.Units.SI.HeatCapacity dUTZ[n] "Derivative of internal energy by temperature at constant moles in each control volume";
  .Modelica.Units.SI.MoleFraction Z[n, Medium.nS] = reaction.y_out "Moler-based composition in each control volume";

initial equation
    
  T = if initFromEnthalpy then Medium.temperature(Medium.setState_phX(p_start[1:n], h_start[2:n+1], X_start)) else T_start[2:n+1];
    
equation

  M = {Medium.MMX*reaction.Zx[i, :] for i in 1:n} "Total mass in each control volume";
  gas.d = M./V;
  U = {M[i]*gas[i].u for i in 1:n} "Internal energy of each control volume";
  MX = {M[i]*gas[i].X for i in 1:n} "Mass of each component in each control volume";

  for i in 1:n loop
    
    dUTZ[i] = Medium.specificHeatCapacityCp(gas[i].state)*M[i] - reaction.Ztot[i]*.Modelica.Constants.R;
    dUdt_cv[i] = der(T[i])*dUTZ[i] + gas[i].dUZT[:]*dZX[i, :];
    dUdt_cv[i] = dH[i] + Q[i] "Energy balance equation";
    
  end for;
  
  // Balances
  H_flow[1] = sum({portA[i].m_flow*actualStream(portA[i].h_outflow) for i in 1:NA});
  H_flow[end] = -sum({portB[i].m_flow*actualStream(portB[i].h_outflow) for i in 1:NB});
  mX_flow[1,:] = {sum({portA[i].m_flow*actualStream(portA[i].X_outflow[k]) for i in 1:NA}) for k in 1:Medium.nS};
  mX_flow[end,:] = {sum(-{portB[i].m_flow*actualStream(portB[i].X_outflow[k]) for i in 1:NB}) for k in 1:Medium.nS};

  for i in 2:n loop

    if generateEventForReversal then

      H_flow[i] = smooth(0, if m_flow[i] >= 0 then gas[i-1].h*m_flow[i] else gas[i].h*m_flow[i]);
      mX_flow[i,:] = smooth(0, if m_flow[i] >= 0 then gas[i-1].X*m_flow[i] else gas[i].X*m_flow[i]);

    else

      H_flow[i] = noEvent(if m_flow[i] >= 0 then gas[i-1].h*m_flow[i] else gas[i].h*m_flow[i]);
      mX_flow[i,:] = noEvent(if m_flow[i] >= 0 then gas[i-1].X*m_flow[i] else gas[i].X*m_flow[i]);

    end if;

  end for;

  // Convection transport terms
  dH = H_flow[1:n]-H_flow[2:n+1];
  dMX = mX_flow[1:n,:]-mX_flow[2:n+1,:];

  sh = cat(1, {levels[1]}, levels[2:n] - levels[1:(n-1)]) .* gas.d * g;

  // Boundary conditions
  portA.p = fill(pA, NA);
  portB.p = fill(pB, NB);
  portA.h_outflow = fill(gas[1].h, NA);
  portB.h_outflow = fill(gas[end].h, NB);
  portA.X_outflow = fill(gas[1].X, NA);
  portB.X_outflow = fill(gas[end].X, NB); 
  q.T = T_wall;
  q_fluid.Q_flow = Q_fluid;
  q_fluid.T = T;
  Q_tot = sum(Q);
  m_flow[1] = sum(portA.m_flow);
  m_flow[end] = -sum(portB.m_flow);
  m_flow_mean = sum(m_flow)/(n+1);

  // Definition of dp vector as differences between ports and control volumes
  dp = if includeStaticHead then

    cat(1,{pA-gas[1].p}, gas[1:n-1].p-gas[2:n].p-sh[1:n-1], {gas[n].p-pB-sh[n]})

  else

    cat(1,{pA-gas[1].p}, gas[1:n-1].p-gas[2:n].p, {gas[n].p-pB});

  // Coupling of dp vector to dp from friction of 0 depending on discretization option
  dp = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then cat(1, {0}, dp_fric)
       elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then cat(1, dp_fric, {0})
       elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then dp_fric
       else cat(1, {0}, dp_fric, {0});

  // Mass flow from friction model
  friction.m_flow = m_flow_fric ./ n_channels_fric;

  mflow_A_in = sum({max(0, portA[i].m_flow) for i in 1:NA});
  mflow_B_in = sum({max(0, portB[i].m_flow) for i in 1:NB});
  if NA == 1 then

    hA = inStream(portA[1].h_outflow);
    XA = inStream(portA[1].X_outflow);

  else

    hA = sum(inStream(portA.h_outflow).*{max(0, portA[i].m_flow) for i in 1:NA})/max(1e-10, mflow_A_in);
    for k in 1:Medium.nS-1 loop

      XA[k] = sum({inStream(portA[i].X_outflow[k])*max(0, portA[i].m_flow) for i in 1:NA})/max(1e-10, mflow_A_in);

    end for;
    sum(XA) = 1;

  end if;
  if NB == 1 then

    hB = inStream(portB[1].h_outflow);
    XB = inStream(portB[1].X_outflow);

  else

    hB = sum(inStream(portB.h_outflow).*{max(0, portB[i].m_flow) for i in 1:NB})/max(1e-10, mflow_B_in);
    for k in 1:Medium.nS-1 loop

      XB[k] = sum({inStream(portB[i].X_outflow[k])*max(0, portB[i].m_flow) for i in 1:NB})/max(1e-10, mflow_B_in);

    end for;
    sum(XB) = 1;

  end if;

  /*Thermodynamic states in cv*/
  state[2:n+1] = gas.state;
  if Medium.analyticInverseTfromh then

    state[1] = Medium.setState_phX(pA, hA, XA);
    TA = Medium.temperature(state[1]);
    state[n+2] = Medium.setState_phX(pB, hB, XB);
    TB = Medium.temperature(state[n+2]);

  else

    state[1] = Medium.setState_pTX(pA, TA, XA);
    hA = Medium.specificEnthalpy_pTX(pA, TA, XA);
    state[n+2] = Medium.setState_pTX(pB, TB, XB);
    hB = Medium.specificEnthalpy_pTX(pB, TB, XB);

  end if;

  /* Heat transfer */
  Q = Q_wall + Q_fluid + reaction.Q + Q_extra;

  for i in 1:n loop

    alpha[i] = htcoeff[i].alphaA;
    if useHeatTransfer then 

      Q_wall[i] = alpha[i] * n_channels[i] * A_heat[i] * (T_wall[i] - gas[i].T);

    else
      
      T_wall[i] = gas[i].T;
   
    end if;

  end for;

  annotation (Documentation(info="<html>
<h4>ReactionChannel</h4>
<p>Model of a pipe with <code>n</code> discrete volumes arranged in the fluid flow direction.  All control volumes have a unique pressure, temperature and mass content.</p>
<h4>Reaction Model</h4>
<p>Any dynamic reaction model that inherits from the <a href=\"modelica://FuelCell.Reactions.Templates.DynamicReaction\">FuelCell.Reactions.Templates.DynamicReaction</a> can be specified by the user to give dynamic reactions along the flow direction on both the primary and the secondary side. The DynamicReaction model is only used as a constraining class and does not represent any reaction by itself.</p>
<h4><span style=\"color:#d2492a\">Connection Principles</span></h4>
<p>By default, the first control volume is connected to <code>portA</code> and a flow model is connected to <code>portB</code>. To assure a numerically sound model of alternating volume and flow models, connect a flow model to <code>portA</code> and a volume model to <code>portB</code>. It is possible to reconfigure the model to have either a volume or friction model at any of the ports, please see the section on flow segmentation below for more information.</p>
<p>The thermal connector <code>q</code> is assumed to be connected to a wall or other external object with a temperature. The selected heat transfer correlation computes heat transfer between the fluid and object connected to <code>q</code>. The thermal connector <code>q_fluid</code> exposes the control volume temperatures and can be used to connect external heat transfer models to the channel, or prescribe heat flow as boundary condition.</p>
<h4>Assumptions</h4>
<ul>
<li>The cross sectional area must be specified for each control volume. The volume is calculated from length and cross sectional area.</li>
<li>The mass flow rate into the first control volume is given by the connector <code>portA.</code> The mass flow rate between each control volume and from the last control volume to <code>portB</code> is calculated in a replaceable model <code>friction</code>, from the pressure drop between the two components.</li>
<li>Heat flow rate into each control volume from the heat connector is calculated. The heat transfer coefficient is calculated in the replaceable model <code>htcoeff</code>.</li>
<li>Dynamic reactions possible in the flow direction.</li>
</ul>
<h4>Parametrization</h4>
<h5>General Tab</h5>
<p><u>Medium</u></p>
<ul>
<li><code>Medium</code>: Defines the package defining medium properties.</li>
</ul>
<p><u>Geometry</u></p>
<p>Each parameter in this group is provided per channel segment. The number of elements in the parameter vector <code>L</code> determines the number of segments <code>n</code>. All other parameters must be provided as vectors with this length.</p>
<ul>
<li><code>n_channels</code>: Allows the component to model multiple parallel channels under identical operating conditions. The model operates by computing all properties for a single channel and then multiplying the mass flow rate and transferred heat with the number of channels. Another way to model multiple parallel channels is to set <code>n_channels</code> = 1 and provide the total cross sectional and heat transfer areas of all channels.</li>
<li><code>L</code>: Length of each channel segment.</li>
<li><code>Dhyd</code>: Hydraulic diameter of the channel.</li>
<li><code>A</code>: Cross-sectional area.</li>
<li><code>V</code>: Channel volume (total in case of multiple channels). Automatically computed for most channel geometries from area and length.</li>
</ul>
<p><u>Static head</u></p>
<ul>
<li><code>includeStaticHead</code>: Set to true to account for pressure difference between the inlet and oulet port due to different elevations.</li>
<li><code>levels</code>: Relative levels of the segments, 0 at inlet port. Only needed if static head is included.</li>
<li><code>g</code>: Acceleration of gravity. The standard gravitational acceleration on earth is the default value. Only needed if static head is included.</li>
</ul>
<p><u>Heat transfer</u></p>
<ul>
<li><code>A_heat</code>: Area transferring heat (orthogonal to heat flow direction) for each segment.</li>
<li><code>HeatTransfer</code>: Replaceable model to characterize the heat transfer coefficient. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase\">this interface class</a>.</li>
</ul>
<p><u>Flow resistance</u></p>
<ul>
<li><code>Friction</code>: Replaceable model to characterize the flow resistance, i.e. a correlation between mass flow rate and pressure drop. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance\">this interface class</a>.</li>
</ul>
<p><br><u>Reaction</u></p>
<ul>
<li><code>Reaction</code>: Replaceable model to characterize reactions in the channel. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://FuelCell.Reactions.Templates.DynamicReaction\">this interface class</a>.</li>
</ul>
<h5>Initialization Tab</h5>
<p><u>Initialization option</u></p>
<ul>
<li><code>initOpt</code>: Determines if initialization is done from fixed values, or if steady-state initialization should be used. Options for steady pressure or steady temperature/enthalpy only is also available.</li>
</ul>
<p><u>Pressure</u></p>
<ul>
<li><code>p_start_in</code>: Initial inlet pressure (in inlet port and first control volume) if the default expression for p_start is not modified.</li>
<li><code>p_start_out</code>: Initial pressure in the outlet port if the default expression for p_start is not modified.</li>
<li><code>p_start</code>: Vector of initial pressures in all control volumes, plus the pressure in the outlet port. Note that the last element is not used for initialization, this element is introduced so that when connecting multiple components in series, the same initial inlet pressure can be used in the second pipe as the initial outlet pressure of the first pipe and still maintaining an initial pressure gradient between all control volumes.</li>
</ul>
<p><u>Enthalpy</u></p>
<ul>
<li><code>initFromEnthlpy</code>: Set to true to initialize by specifying initial specific enthalpies. If false, initial temperatures have to be specified instead.</li>
<li><code>h_start_in</code>: Initial inflowing specific enthalpy if the default expression for h_start is not modified.</li>
<li><code>h_start_out</code>: Initial specific enthalpy in the last control volume if the default expression for h_start is not modified.</li>
<li><code>h_start</code>: Initial specific enthalpies. The first element is assumed for the inlet flow, the following values are initial values for all control volumes. Note that the first element is not used for initialization, but present so that the same initial inlet enthalpy can be specified as the outlet enthalpy of the upstream component.</li>
</ul>
<p><u>Temperature</u></p>
<p>The parameters in this group are used if initFromEnthlpy = false.</p>
<ul>
<li><code>T_start_in</code>: Initial inflowing temperature if the default expression for T_start is not modified.</li>
<li><code>T_start_out</code>: Initial temperature in the last control volume if the default expression for T_start is not modified.</li>
<li><code>T_start</code>: Initial control volume temperatures. The first element is assumed for the inlet flow, the following values are initial values for all control volumes. Note that the first element is not used for initialization, but present so that the same initial inlet temperature can be specified as the outlet temperature of the upstream component.</li>
</ul>
<p><u>Mass fractions</u></p>
<ul>
<li><code>X_start</code>: Initial mass fraction for all fluid species. Identical initial compositions are prescribed to all control volumes.</li>
</ul>
<p><u>Mass flow start</u></p>
<ul>
<li><code>m_flow_start</code>: Initial guess value for mass flow rates in cases when the mass flow rate must be computed by solving non-linear systems of equations.</li>
</ul>
<h5>Advanced Tab</h5>
<p><u>Numerics</u></p>
<ul>
<li><code>positiveFlow</code>: Set to true if only positive flow rate can be assumed. In this case upstream fluid properties are independent of the direction of the flow, which generally reduce the model complexity. Note that enabling this parameter does not prevent reversing flow.</li>
<li><code>from_dp</code>: This parameter defines for which causality the flow resistance model is expressed explicitly. Set to true if mass flow rate is computed from pressure drop. Note that all correlations cannot be explicitly expressed for both causalities, and in such cases this flag is ignored. Also note that this parameter does not control the causality of the experiment, only the selected boundary conditions and state selection does that.</li>
<li><code>dp_smooth</code>: Pressure drop interval used for regularization, used when from_dp = true. It is up to the flow resistance model to use this parameter. Most correlations will introduce regularization in this interval to increase numerical robustness.</li>
<li><code>mflow_smooth</code>: Mass flow rate interval used for regularization, used when from_dp = false. It is up to the flow resistance model to use this parameter. Most correlations will introduce regularization in this interval to increase numerical robustness.</li>
<li><code>dp_asState</code>: Setting this parameter to true enforces the pressure drops to be used as numerical states instead of the absolute pressures. This can be beneficial when there is low pressure drop between the segments.</li>
<li><code>frictionDistribution</code>: This parameter controls if flow resistance (friction) or control volume models are exposed in the fluid connectors. See further description below.</li>
</ul>
<p><u>Event handling</u></p>
<ul>
<li><code>generateEventForReversal</code>: Set to true to enable model events when the flow rate changes sign. This will give more accurate results but slow down simulation speed when there are several zero crossings.</li>
</ul>
<p><u>Calibration factors</u></p>
<ul>
<li><code>CF_PressureLoss</code>: The resulting pressure drop due to friction is multiplied with this factor.</li>
<li><code>CF_HeatTransfer</code>: The resulting heat transfer is multiplied with this factor.</li>
</ul>
<h4>Flow Segmentation</h4>
<p>The parameter <code>frictionDistribution</code> can be used to control the type of model exposed in the fluid connector. This option is available to facilitate creation of numerically sound system models with alternating control volume and flow resistance models. The available options are:</p>
<ul>
<li><code>Volume&nbsp;-&nbsp;...&nbsp;-&nbsp;Friction</code> (default). Control volume model at portA and friction model at portB.</li>
<li><code>Friction&nbsp;-&nbsp;...&nbsp;-&nbsp;Volume</code>. Friction model at portA and volume model at portB.</li>
<li><code>Friction&nbsp;-&nbsp;...&nbsp;-&nbsp;Friction</code>. Friction models at both portA and portB.</li>
<li><code>Volume&nbsp;-&nbsp;...&nbsp;-&nbsp;Volume</code>. Volume models at both portA and portB.</li>
</ul>
<p><br>The number of control volumes do not change with the different options above and neither does the number of parameters that have to be set. In the two first options, there are equal numbers of flow resistance and control volume instances in the model, so the geometric parameters are propagated directly from the segment to the flow resistance models. In the two last options there is one more or one less friction model than the number of control volumes.</p>
<p>In the case there is one more friction model, the first friction model will use the average area and hydraulic diameter and account for a factor 1/(n + 1) of the length, whereas the remaining friction models will account for the same relative amount of the remaining length as they do according to the paraemters.</p>
<p>In the case there is one less friction model, the entered parmeter value for area and hydraulic diameter for the last segment are ignored. The length entered for the last segment is accounted for by the remaining segments in such a way that they retain the relative length per segment.</p>
<p><b>Note:</b> Due to the handling of geometry parameters described above, option 3 and 4 should be avoided if the cross section area is not invariate along the pipe, the segments are not of equal lengths, or the segments do not include equal number of channels.</p>
<p></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={Line(
          points={{-64,-54},{56,-54}},
          color={0,0,0},
          thickness=0.5), Polygon(
          points={{68,-54},{56,-50},{56,-58},{68,-54}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-58,8},{-42,-8}},
          lineColor={255,128,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,128,0}),
        Ellipse(
          extent={{42,8},{58,-8}},
          lineColor={255,128,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,128,0})}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
            graphics));

end ReactionChannel;