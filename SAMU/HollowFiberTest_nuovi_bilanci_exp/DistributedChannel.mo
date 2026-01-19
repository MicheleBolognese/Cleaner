within SAMU.HollowFiberTest_nuovi_bilanci_exp;

model DistributedChannel
  "Flow channel (mass transfer) model with distributed volumes and pressure drop"

  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.Channel(summary(
    dp=pA - pB,
    M=Mtot,
    V=Vtot,
    m_flow=m_flow_mean,
    T_out=volume[end].T,
    h_out=volume[end].h,
    T_in=TA,
    h_in=hA,
    d_in=Medium.density(state[1])),
    redeclare .Modelon.ThermoFluid.Interfaces.VolumePort portA,
    redeclare .Modelon.ThermoFluid.Interfaces.FlowPort portB);

  .Modelon.ThermoFluid.Interfaces.TempHeatPort[n] q_fluid
    "Connector exposing fluid temperature"
    annotation (Placement(transformation(extent={{-10.0,-10.0},{10.0,10.0}},rotation = 0.0,origin = {0.0,0.0})));

  replaceable model Friction =
      .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss
                                                                                               constrainedby
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Interfaces.PartialPipeResistance
    "Friction model"   annotation(choicesAllMatching, Dialog(group="Flow resistance"));

  /* Heat Transfer*/
  replaceable model HeatTransfer =
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient(alpha=0)
    constrainedby
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase
    "Heat transfer coefficient model" annotation(choicesAllMatching,Dialog(group="Heat transfer"));

  parameter Boolean includeStaticHead = false "Consider static head"
    annotation(Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Momentum Balance"));
  parameter .Modelica.Units.SI.Length[n] levels=fill(0, n)
    "Relative levels of control volume outlets (inlet port at 0)" annotation (Dialog(
      enable=includeStaticHead,
      tab="Balance Equations",
      group="Momentum Balance"));
  parameter .Modelica.Units.SI.Acceleration g=.Modelica.Constants.g_n "Acceleration of gravity" annotation (Dialog(
      enable=includeStaticHead,
      tab="Balance Equations",
      group="Momentum Balance"));
  parameter Boolean includeAcceleration=false
    "If true, includes acceleration term in the momentum balance"
    annotation (Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Momentum Balance"));
  parameter Boolean kineticEnergyInBalance=false
    "If true, includes kinetic energy in the energy balance"
    annotation (Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Energy Balance"));

  /* Advanced parameters */
  parameter Boolean dp_asState=false
    "Use pressure difference as state variable"                                  annotation(Evaluate=true, Dialog(tab="Advanced",group="Numerics"));
  parameter .Modelon.ThermoFluid.Choices.FrictionDistribution frictionDistribution=
    .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
    "Discretization scheme for friction and control volumes (see info)" annotation(Evaluate=true, Dialog(tab="Advanced",group="Numerics"));
  final parameter Integer n_fric = if (frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
    frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol) then n
    elseif frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then n+1
    else n-1 annotation(Evaluate=true);

  input Real CF_PressureLoss=1.0 "Calibration factor for pressure drop"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  input Real CF_HeatTransfer=1.0 "Calibration factor for heat transfer"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));

  parameter Real flowFraction0 = 1 "Fraction of nominal mass flow rate through this component"
    annotation(Dialog(tab="Advanced", group="Interpretation of nominal parameters in pressure drop correlations"));
  parameter Real dpFraction0 = 1 "Fraction of nominal pressure drop over this component"
    annotation(Dialog(tab="Advanced", group="Interpretation of nominal parameters in pressure drop correlations"));
  parameter Boolean useMeanTempDrivenQ=false
    "If true, heat flow in each control volume is driven by the average temperature" annotation(Dialog(tab="Advanced", group="Heat transfer"));

  /* Length calibration */
  input Real CF_length=1.0 "Calibration factor for pipe length"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  .Modelica.Units.SI.Length L_total_internal=L_total*CF_length "Modified pipe length for design calculation";
  .Modelica.Units.SI.Volume[n] V_internal=V*CF_length "Modified pipe volume for design calculation";
  .Modelica.Units.SI.Area[n] A_heat_internal=A_heat*CF_length
    "Modified heat transfer area for design calculation";

  Friction[n_fric] friction(
    each m_flow(start=m_flow_start),
    redeclare each package Medium = Medium,
    dp=dp,
    stateA_inflow=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then
        state[2:n+1]
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then
        state[1:n]
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        state[1:n+1]
      else
        state[2:n]),
    stateB_inflow=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then
        state[3:n+2]
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then
        state[2:n+1]
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        state[2:n+2]
      else
        state[3:n+1]),
    final A = (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {sum(A)/n}, A)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        A[1:n-1]
      else A),
    each final L=L_total_internal,
    final Dhyd=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {sum(Dhyd)/n}, Dhyd)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        Dhyd[1:n-1]
      else Dhyd),
    final n_channels=n_channels_fric,
    final lengthFraction = (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {1/(n+1)}, L/sum(L)*(n/(n+1)))
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        L[1:n-1]/sum(L[1:n-1])
      else L/sum(L)),
    each final from_dp = from_dp,
    each final positiveFlow = positiveFlow,
    final dp_smooth = dp_smooth *
      (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {1/(n+1)}, L/sum(L)*(n/(n+1)))
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        L[1:n-1]/sum(L[1:n-1])
      else L/sum(L)),
    final mflow_smooth = mflow_smooth./n_channels_fric,
    each final F_user=CF_PressureLoss,
    each final flowFraction0 = flowFraction0,
    each final dpFraction0 = dpFraction0,
    final outletInstance = cat(1, fill(false, n_fric-1), {true})) "Flow resistance model";

  HeatTransfer[n] htcoeff(
    redeclare each package Medium = Medium,
    final A = A,
    final Dhyd = Dhyd,
    each final L = L_total,
    each final F_user = CF_HeatTransfer,
    m_flow = m_flow[1:n]./n_channels,
    stateA = state[2:n+1],
    stateB = state[2:n+1],
    each final CF_length=CF_length) "Heat transfer model";

  Medium.VolumeDynamics[n] volume(
    each final thermalDynamics=true,
    each initOpt= initOpt,
    p_start = p_start[1:n],
    T_start = if initFromEnthalpy then Medium.temperature(Medium.setState_phX(p_start[1:n], h_start[2:n+1], Medium.reference_X)) else T_start[2:n+1],
    h_start = if initFromEnthalpy then h_start[2:n+1] else Medium.specificEnthalpy(Medium.setState_pTX(p_start[1:n],T_start[2:n+1],Medium.reference_X)),
    each C_start = C_start,
    each X_start = X_start,
    dE = dE,
    dMX = transpose(dMX) + rMX,
    dMC = transpose(dMC),
    final V_tot = V_internal,
    final T_in=T) "Control volume models";

  .Modelica.Units.SI.Pressure[n_fric] dp(start=if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
         then p_start[1:n] - p_start[2:n + 1] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n - 1] - p_start[2:n]) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n] - p_start[2:n + 1]) else p_start[1:n - 1] - p_start[2:n], each stateSelect=if dp_asState
         then StateSelect.always else StateSelect.default)
    "Pressure drops from friction models, vector length depending on discretization option";

  /* Fluid properties in control volumes */
  .Modelica.Units.SI.Density[n] d;
  .Modelica.Units.SI.SpecificEnthalpy[n] h "Control volume specific enthalpies";
  .Modelica.Units.SI.Pressure[n] p(start=p_start[1:n]) "Control volume pressures";
  .Modelica.Units.SI.Mass[n] M "Control volume masses";
  .Modelica.Units.SI.Mass Mtot "Total mass";
  .Modelica.Units.SI.Volume Vtot "Total volume";

  .Modelica.Units.SI.Pressure[n] sh "Static head between control volumes";
  .Modelica.Units.SI.MassFlowRate m_flow_mean "Average mass flow rate in pipe, positive from portA to portB";
  .Modelica.Units.SI.MassFlowRate[n + 1] m_flow(each start=m_flow_start)
    "Mass flow over cv boundaries, positive from portA towards portB";
  .Modelica.Units.SI.AbsolutePressure pA "Pressure in portA";
  .Modelica.Units.SI.AbsolutePressure pB "Pressure in portB";
  .Modelon.Media.Units.Temperature TA "Upstream temperature when flow A -> B";
  .Modelon.Media.Units.Temperature TB "Upstream temperature when flow B -> A";
  .Modelon.Media.Units.Temperature TA_out "Outlet temperature when flow B -> A";
  .Modelon.Media.Units.Temperature TB_out "Outlet temperature when flow A -> B";
  .Modelica.Units.SI.Density dA "Inlet density, nominal direction, portA";
  .Modelica.Units.SI.Density dB "Inlet density at portB";
  Medium.ThermodynamicState[n+2] state
    "Thermodynamic states in control volumes and instream from ports";
  Medium.ThermodynamicState stateA_out
    "Thermodynamic state flowing out of portA (in case of back flow)";
  Medium.ThermodynamicState stateB_out
    "Thermodynamic state flowing out of portB (in case of forward flow)";
  .Modelica.Units.SI.HeatFlowRate[n] Q "Heat flow rate into control volumes";
  .Modelica.Units.SI.HeatFlowRate Q_tot "Total heat flow rate into pipe";
  .Modelica.Units.SI.HeatFlowRate[n] Q_wall=q.Q_flow "Heat flow rate from connector q";

  .Modelica.Units.SI.HeatFlowRate[n] dE "Internal energy derivative";
  .Modelica.Units.SI.EnthalpyFlowRate H_flow[n + 1] "Enthalpy flow over control volume boundaries";
  .Modelica.Units.SI.MassFlowRate MX_flow[Medium.nS,n + 1] "Species mass flow rate over cv boundaries";
  Real MC_flow[Medium.nC,n+1] "Trace components flow rate over cv boundaries";

  .Modelica.Units.SI.CoefficientOfHeatTransfer[n] alpha "Heat transfer coefficient";

  .Modelica.Units.SI.Temperature T_fluid[n] "Temperature for heat transfer modeling";

  //Variables for deatiled balance formulations
  Real[n+1] drhodx(each final quantity="Density", each final unit="kg/m3")
    "Density differences between adjacent control volumes";
  .Modelica.Units.SI.PressureDifference[n + 1] dpaccel
    "Resulting pressure drop from acceleration due to specific volume change";

  .Modelica.Units.SI.Velocity v[n + 1] "Fluid velocities at control volume boundaries";
  .Modelica.Units.SI.Power Ek_flow[n + 1] "Kinetic energy flow m_flow*v^2/2 over control volume boundaries";
  .Modelica.Units.SI.Power dEk_flow[n] "Difference of kinetic energy flow over boundaries, per control volume";
  .Modelica.Units.SI.MassFlowRate[n,Medium.nS] rMX=zeros(n, Medium.nS) "Mass residual for balance from diffusion"
    annotation (Dialog(group="Mass transfer"));

protected
  .Modelica.Units.SI.MassFlowRate mflow_A_in "Total inlet flow rate at portA";
  .Modelica.Units.SI.MassFlowRate mflow_B_in "Total inlet flow rate at portB";
  .Modelica.Units.SI.MassFlowRate[n_fric] m_flow_fric=if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
       then m_flow[2:n + 1] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
       then m_flow[1:n] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
       then m_flow else m_flow[2:n] "Mass flow over cv boundaries, positive from portA towards portB";
public
  .Modelica.Units.SI.Pressure[n + 1] dp_internal(start=if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
         then cat(
        1,
        {0},
        p_start[1:n] - p_start[2:n + 1]) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n - 1] - p_start[2:n],
        {0}) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n] - p_start[2:n + 1]) else cat(
        1,
        {0},
        p_start[1:n - 1] - p_start[2:n],
        {0}) "Pressure drop between control volumes and ports");
//protected
  .Modelica.Units.SI.SpecificEnthalpy hA "Mixture enthalpy of inlet flows at portA";
  .Modelica.Units.SI.SpecificEnthalpy hB "Mixture enthalpy of inlet flows at portB";
  .Modelica.Units.SI.MassFraction[Medium.nS] XA "Mixture mass fractions of inlet flows at portA";
  .Modelica.Units.SI.MassFraction[Medium.nS] XB "Mixture mass fractions of inlet flows at portB";

  final parameter Real[n_fric] n_channels_fric=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {n_channels[1]}, n_channels)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        n_channels[1:n-1]
      else n_channels);
  .Modelica.Units.SI.Temperature[n] T_wall "Wall temperature";
  .Modelica.Units.SI.HeatFlowRate[n] Q_fluid "Heat flow rate from connector q_fluid";
  .Modelica.Units.SI.MassFlowRate[Medium.nS,n] dMX "Mass derivative of substances";
  Real[Medium.nC,n] dMC "Trace component derivatives";


equation
  assert( not (frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol and n == 1),
    "Volumes on both ports is not supported for n = 1 (one control volume)");

  // Definition of dp vector as differences between ports and control volumes
  dp_internal = if includeStaticHead then
    cat(1,{pA-volume[1].p+dpaccel[1]}, volume[1:n-1].p-volume[2:n].p-sh[1:n-1]+dpaccel[2:n], {volume[n].p-pB-sh[n]+dpaccel[n]})
    else
    cat(1,{pA-volume[1].p+dpaccel[1]}, volume[1:n-1].p-volume[2:n].p+dpaccel[2:n], {volume[n].p-pB+dpaccel[n]});

  // Coupling of dp_internal vector to dp from friction of 0 depending on discretization option
  dp_internal = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then
      cat(1, {0}, dp)
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then
      cat(1, dp, {0})
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      dp
    else
      cat(1, {0}, dp, {0});

  drhodx = cat(1, {dA - d[1]}, {d[i - 1] - d[i] for i in 2:n}, {d[n] - dB});

  dpaccel = if includeAcceleration then
    cat(1,{-(m_flow[1]*m_flow[1]/max(1e-6, (A[1]*A[1]*d[1]*dA)))*drhodx[1]},
    {-(m_flow[i]*m_flow[i]/max(1e-6, (A[i]*A[i]*d[i - 1]*d[i])))*drhodx[i] for i in 2:n},
    {-(m_flow[n + 1]*m_flow[n + 1]/max(1e-6, (A[n]*A[n]*d[n]*dB)))*
      drhodx[n+1]})
    else
      fill(0, n + 1);

  /* Fluid properties in control volumes */
  M = d .* V_internal;
  p = volume.p;
  h = volume.h;
  d = volume.d;
  Vtot=sum(V_internal);
  Mtot=sum(M);

  // Mass flow from friction model
  friction.m_flow = m_flow_fric ./ n_channels_fric;

  sh = cat(1, {levels[1]}, levels[2:n] - levels[1:(n-1)]) .* volume.d * g;
  m_flow[1] = sum(portA.m_flow);
  m_flow[n+1] = - sum(portB.m_flow);
  portA.p = fill(pA, NA);
  portB.p = fill(pB, NB);
  portA.h_outflow = fill(volume[1].h, NA);
  portB.h_outflow = fill(volume[n].h, NB);
  portA.X_outflow = fill(volume[1].X, NA);
  portB.X_outflow = fill(volume[n].X, NB);
  portA.C_outflow = fill(volume[1].C, NA);
  portB.C_outflow = fill(volume[n].C, NB);
  q.T = T_wall;
  q_fluid.Q_flow = Q_fluid;
  q_fluid.T = T;
  T = volume.T;
  Q_tot = sum(Q);

  m_flow_mean = sum(m_flow)/(n+1);
  mflow_A_in = sum({max(0, portA[i].m_flow) for i in 1:NA});
  mflow_B_in = sum({max(0, portB[i].m_flow) for i in 1:NB});
  if NA == 1 then
    hA = inStream(portA[1].h_outflow);
    XA = inStream(portA[1].X_outflow);
  else
    hA = sum(inStream(portA.h_outflow).*{max(1e-10, portA[i].m_flow) for i in 1:NA})/
      max(NA*1e-10, mflow_A_in);
      for k in 1:Medium.nS-1 loop
       XA[k] = sum({inStream(portA[i].X_outflow[k])*max(1e-10, portA[i].m_flow) for i in 1:NA})/
        max(NA*1e-10, mflow_A_in);
      end for;
    sum(XA) = 1;
  end if;
  if NB == 1 then
    hB = inStream(portB[1].h_outflow);
    XB = inStream(portB[1].X_outflow);
  else
    hB = sum(inStream(portB.h_outflow).*{max(1e-10, portB[i].m_flow) for i in 1:NB})/
     max(NB*1e-10, mflow_B_in);
    for k in 1:Medium.nS-1 loop
       XB[k] = sum({inStream(portB[i].X_outflow[k])*max(1e-10, portB[i].m_flow) for i in 1:NB})/
        max(NB*1e-10, mflow_B_in);
    end for;
    sum(XB) = 1;
  end if;

  T_fluid = if useMeanTempDrivenQ then
    (if n>1 then
      cat(1, {.Modelon.Math.Smoothing.spliceFunction((TA+volume[1].T)/2, (volume[1].T+volume[2].T)/2, m_flow[1], mflow_smooth)},
        {.Modelon.Math.Smoothing.spliceFunction((volume[i-1].T+volume[i].T)/2, (volume[i].T+volume[i+1].T)/2, m_flow[i], mflow_smooth) for i in 2:n-1},
        {.Modelon.Math.Smoothing.spliceFunction((volume[n-1].T+volume[n].T)/2, (volume[n].T+TB)/2, m_flow[n], mflow_smooth)})
      else
        fill(.Modelon.Math.Smoothing.spliceFunction((TA+volume[1].T)/2, (volume[1].T+TB)/2, m_flow[1], mflow_smooth),n))
    else
      volume.T;

  dA = Medium.density(state[1]);
  dB = Medium.density(state[n+2]);

  /* Propreties flow over cv boundaries */
  H_flow[1] = sum({portA[i].m_flow*actualStream(portA[i].h_outflow) for i in 1:NA});
  H_flow[n+1] = -sum({portB[i].m_flow*actualStream(portB[i].h_outflow) for i in 1:NB});
  MX_flow[:,1] = {sum({portA[i].m_flow*actualStream(portA[i].X_outflow[k]) for i in 1:NA}) for k in 1:Medium.nS};
  MX_flow[:,n+1] = {-sum({portB[i].m_flow*actualStream(portB[i].X_outflow[k]) for i in 1:NB}) for k in 1:Medium.nS};
  MC_flow[:,1] = {sum({portA[i].m_flow*actualStream(portA[i].C_outflow[k]) for i in 1:NA}) for k in 1:Medium.nC};
  MC_flow[:,n+1] = {-sum({portB[i].m_flow*actualStream(portB[i].C_outflow[k]) for i in 1:NB}) for k in 1:Medium.nC};
  for i in 2:n loop
    if generateEventForReversal then
      H_flow[i] = smooth(0, if m_flow[i] >=0 then m_flow[i]*volume[i-1].h else m_flow[i]*volume[i].h);
      MX_flow[:,i] = smooth(0, if m_flow[i] >=0 then m_flow[i]*volume[i-1].X else m_flow[i]*volume[i].X);
      MC_flow[:,i] = smooth(0, if m_flow[i] >=0 then m_flow[i]*volume[i-1].C else m_flow[i]*volume[i].C);
    else
      H_flow[i] = noEvent(if m_flow[i] >=0 then m_flow[i]*volume[i-1].h else m_flow[i]*volume[i].h);
      MX_flow[:,i] = noEvent(if m_flow[i] >=0 then m_flow[i]*volume[i-1].X else m_flow[i]*volume[i].X);
      MC_flow[:,i] = noEvent(if m_flow[i] >=0 then m_flow[i]*volume[i-1].C else m_flow[i]*volume[i].C);
    end if;
  end for;

  /*Thermodynamic states in cv*/
  if Medium.analyticInverseTfromh then
    state[1] = Medium.setState_phX(pA, hA, XA);
    TA = Medium.temperature(state[1]);
    state[n+2] = Medium.setState_phX(pB, hB, XB);
    Medium.temperature(state[n+2]) = TB;
    stateA_out =Medium.setState_phX(pA, h[1], volume[1].X);
    TA_out = Medium.temperature(stateA_out);
    stateB_out =Medium.setState_phX(pB, h[n], volume[n].X);
    TB_out = Medium.temperature(stateB_out);
  else
    state[1] = Medium.setState_pTX(pA, TA, XA);
    hA = Medium.specificEnthalpy(state[1]);
    state[n+2] = Medium.setState_pTX(pB, TB, XB);
    Medium.specificEnthalpy(state[n+2]) = hB;
    stateA_out = Medium.setState_pTX(pA, TA_out, volume[1].X);
    h[1] = Medium.specificEnthalpy(stateA_out);
    stateB_out = Medium.setState_pTX(pB, TB_out, volume[n].X);
    h[n] = Medium.specificEnthalpy(stateB_out);
  end if;
  state[2:n+1] = volume.state;

  /* Velocity & kinetic energy */
  v = cat(1, {noEvent(if sum({portA[i].m_flow for i in 1:NA}) >= 0 then sum(
       {portA[i].m_flow for i in 1:NA})/max(1e-6, (dA*A[1])) else sum({portA[i].m_flow
      for i in 1:NA})/max(1e-6, (d[1]*A[1])))},
        {if generateEventForReversal then
          smooth(0, if m_flow[i] >= 0 then m_flow[i]/max(1e-6, (d[
          i - 1]*A[i])) else m_flow[i]/max(1e-6, (d[i]*A[i])))
          else
          noEvent(if m_flow[i] >= 0 then m_flow[i]/max(1e-6, (d[i -
          1]*A[i])) else m_flow[i]/max(1e-6, (d[i]*A[i]))) for i in 2:n},
        {noEvent(if -sum({portB[i].m_flow for i in 1:NB}) >= 0 then
          -sum({portB[i].m_flow for i in 1:NB})/max(1e-6, (d[n]*A[n])) else
          -sum({portB[i].m_flow for i in 1:NB})/max(1e-6, (dB*A[n])))});

  Ek_flow = if kineticEnergyInBalance then
    cat(1, {(sum({portA[i].m_flow for i in 1:NA})*
        .Modelon.Math.Smoothing.regExp(v[1], 0.0001, 2))/2},
      {(m_flow[i]*.Modelon.Math.Smoothing.regExp(v[i], 0.0001, 2))/2 for i in 2:n},
      {(-sum({portB[i].m_flow for i in 1:NB})*
        .Modelon.Math.Smoothing.regExp(v[n + 1], 0.0001, 2))/2})
      else
        fill(0,n+1);

  dEk_flow = Ek_flow[1:n] - Ek_flow[2:n + 1];

  /* Net mass and enthalpy flow into control volumes */
  dE = H_flow[1:n] - H_flow[2:n+1] + Q_fluid + dEk_flow;
  dMX = MX_flow[:,1:n] - MX_flow[:,2:n+1];
  dMC = MC_flow[:,1:n] - MC_flow[:,2:n+1];

  /* Heat transfer */
  Q = Q_wall + Q_fluid;

  for i in 1:n loop
    alpha[i] = htcoeff[i].alphaA;
    Q_wall[i] = alpha[i] * n_channels[i] * A_heat_internal[i] * (T_wall[i] - T_fluid[i]);
  end for;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                      graphics), Icon(graphics={
        Rectangle(
          extent={{-80,40},{80,-40}},
          lineColor={0,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.HorizontalCylinder),
        Ellipse(
          extent={{-58,8},{-42,-8}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={0,0,0}),
        Ellipse(
          extent={{42,8},{58,-8}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={0,0,0})}),
    Documentation(info="<html>
<p>A pipe modeled as <code><span style=\"font-family: Courier New,courier;\">n</span></code> discrete control volumes and distributed flow resistance. Each control volume introduces an independent fluid state. The numerical state variables that are actually used depends on the selected medium model.</p>
<p>In this component, geometrical parameters are provided per control volume segment to allow for a wide range of configurations. The length of the parameter <code><span style=\"font-family: Courier New,courier;\">L</span></code> determines the number of segments, and all other geometric parameters must be assigned as vectors with the same length. If all segments are of the same length and cross-section, a model with simplified parameterization is available <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.DistributedPipe\">here</a>.</p>
<h4>Connection to other components</h4>
<p>By default, the first control volume is connected to <code><span style=\"font-family: Courier New,courier;\">portA</span></code> and a flow model is connected to <code><span style=\"font-family: Courier New,courier;\">portB</span></code>. To assure a numerically sound model of alternating volume and flow models, connect a flow model to <code><span style=\"font-family: Courier New,courier;\">portA</span></code> and a volume model to <code><span style=\"font-family: Courier New,courier;\">portB</span></code>. It is possible to reconfigure the model to have either a volume or friction model at any of the ports, please see the section on flow segmentation below for more information.</p>
<p>The thermal connector <code><span style=\"font-family: Courier New,courier;\">q</span></code> is assumed to be connected to a wall or other external object with a temperature. The selected heat transfer correlation computes heat transfer between the fluid and object connected to <code><span style=\"font-family: Courier New,courier;\">q</span></code>. The thermal connector <code><span style=\"font-family: Courier New,courier;\">q_fluid</span></code> exposes the control volume temperatures and can be used to connect external heat transfer models to the channel, or prescribe heat flow as boundary condition.</p>
<p>This model accounts for mass transfer through the variable <code>rMX</code>.</p>
<h4>Model assumptions</h4>
<ul>
<li>The pipe is discretized in n segments from inlet to outlet.</li>
<li>The mass flow rate into the first control volume is given by the connector <code><span style=\"font-family: Courier New,courier;\">portA.</span></code> The mass flow rate between each control volume and from the last control volume to <code><span style=\"font-family: Courier New,courier;\">portB </span></code>is calculated in a replaceable model <code><span style=\"font-family: Courier New,courier;\">Friction, </span></code>from the pressure drop between the two adjacent control volumes.</li>
<li>Heat transfer to each control volume from the heat connector is modeled. The heat transfer coefficient is modeled by the replaceable model <code><span style=\"font-family: Courier New,courier;\">HeatTransfer</span></code>.</li>
</ul>
<h4>Parameterization</h4>
<p><i>General</i> tab:</p>
<p><u>Medium</u></p>
<ul>
<li>Medium: Defines the package defining medium properties.</li>
</ul>
<p><u>Geometry</u></p>
<p>Each parameter in this group is provided per channel segment. The number of elements provided for the parameter L determines the number of segments n. All other parameters must be provided as vectors with this length.</p>
<ul>
<li>n_channels: Allows the component to model multiple parallel channels under identical operating conditions. The model operates by computing all properties for a single channel and then multiplying the mass flow rate and transferred heat with the number of channels. Another way to model multiple parallel channels is to set n_channels = 1 and provide the total cross sectional and heat transfer areas of all channels.</li>
<li>L: Length of each channel segment.</li>
<li>Dhyd: Hydraulic diameter of the channel.</li>
<li>A: Cross-sectional area.</li>
<li>V: Channel volume (total in case of multiple channels). Automatically computed for most channel geometries from area and length.</li>
</ul>
<p><u>Static head</u></p>
<ul>
<li>includeStaticHead: Set to true to account for pressure difference between the inlet and oulet port due to different elevations.</li>
<li>levels: Relative levels of the segments, 0 at inlet port. Only needed if static head is included.</li>
<li>g: Acceleration of gravity. The standard gravitational acceleration on earth is the default value. Only needed if static head is included.</li>
</ul>
<p><u>Heat transfer</u></p>
<ul>
<li>A_heat: Area transferring heat (orthogonal to heat flow direction) for each segment.</li>
<li>HeatTransfer: Replaceable model to characterize the heat transfer coefficient. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.Interfaces.HTCoefficientsBase\">this interface class</a>.</li>
<li>useMeanTempDrivenQ: If true, this would configure the heat flow to be driven by the mean temperature in each control volume. This would reduce the number of segments needed in the case where temperature difference between the fluids is small compared to the temperature gradients.</li>
</ul>
<p><u>Flow resistance</u></p>
<ul>
<li>Friction: Replaceable model to characterize the flow resistance, i.e. a correlation between mass flow rate and pressure drop. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.Templates.PartialPipeResistance\">this interface class</a>.</li>
</ul>
<p><i>Initialization</i> tab:</p>
<p><u>Initialization option</u></p>
<ul>
<li>initOpt: Determines if initialization is done from fixed values, or if steady-state initialization should be used. Options for steady pressure or steady temperature/enthalpy only is also available.</li>
</ul>
<p><u>Pressure</u></p>
<ul>
<li>p_start_in: Initial inlet pressure (in inlet port and first control volume) if the default expression for p_start is not modified.</li>
<li>p_start_out: Initial pressure in the outlet port if the default expression for p_start is not modified.</li>
<li>p_start: Vector of initial pressures in all control volumes, plus the pressure in the outlet port. Note that the last element is not used for initialization, this element is introduced so that when connecting multiple components in series, the same initial inlet pressure can be used in the second pipe as the initial outlet pressure of the first pipe and still maintaining an initial pressure gradient between all control volumes.</li>
</ul>
<p><u>Enthalpy</u></p>
<ul>
<li>initFromEnthlpy: Set to true to initialize by specifying initial specific enthalpies. If false, initial temperatures have to be specified instead.</li>
<li>h_start_in: Initial inflowing specific enthalpy if the default expression for h_start is not modified.</li>
<li>h_start_out: Initial specific enthalpy in the last control volume if the default expression for h_start is not modified.</li>
<li>h_start: Initial specific enthalpies. The first element is assumed for the inlet flow, the following values are initial values for all control volumes. Note that the first element is not used for initialization, but present so that the same initial inlet enthalpy can be specified as the outlet enthalpy of the upstream component.</li>
</ul>
<p><u>Temperature</u></p>
<p>The parameters in this group are used if initFromEnthlpy = false.</p>
<ul>
<li>T_start_in: Initial inflowing temperature if the default expression for T_start is not modified.</li>
<li>T_start_out: Initial temperature in the last control volume if the default expression for T_start is not modified.</li>
<li>T_start: Initial control volume temperatures. The first element is assumed for the inlet flow, the following values are initial values for all control volumes. Note that the first element is not used for initialization, but present so that the same initial inlet temperature can be specified as the outlet temperature of the upstream component.</li>
</ul>
<p><u>Mass fractions</u></p>
<ul>
<li>X_start: Initial mass fraction for all fluid species. Identical initial compositions are prescribed to all control volumes.</li>
</ul>
<p><u>Mass flow start</u></p>
<ul>
<li>m_flow_start: Initial guess value for mass flow rates in cases when the mass flow rate must be computed by solving non-linear systems of equations.</li>
</ul>
<p><i>Advanced</i> tab:</p>
<p><u>Numerics</u></p>
<ul>
<li>positiveFlow: Set to true if only positive flow rate can be assumed. In this case upstream fluid properties are independent of the direction of the flow, which generally reduce the model complexity. Note that enabling this parameter does not prevent reversing flow.</li>
<li>from_dp: This parameter defines for which causality the flow resistance model is expressed explicitly. Set to true if mass flow rate is computed from pressure drop. Note that all correlations cannot be explicitly expressed for both causalities, and in such cases this flag is ignored. Also note that this parameter does not control the causality of the experiment, only the selected boundary conditions and state selection does that.</li>
<li>dp_smooth: Pressure drop interval used for regularization, used when from_dp = true. It is up to the flow resistance model to use this parameter. Most correlations will introduce regularization in this interval to increase numerical robustness.</li>
<li>mflow_smooth: Mass flow rate interval used for regularization, used when from_dp = false. It is up to the flow resistance model to use this parameter. Most correlations will introduce regularization in this interval to increase numerical robustness.</li>
<li>dp_asState: Setting this parameter to true enforces the pressure drops to be used as numerical states instead of the absolute pressures. This can be beneficial when there is low pressure drop between the segments.</li>
<li>frictionDistribution: This parameter controls if flow resistance (friction) or control volume models are exposed in the fluid connectors. See further description below.</li>
</ul>
<p><u>Event handling</u></p>
<ul>
<li>generateEventForReversal: Set to true to enable model events when the flow rate changes sign. This will give more accurate results but slow down simulation speed when there are several zero crossings.</li>
</ul>
<p><u>Calibration factors</u></p>
<ul>
<li>CF_PressureLoss: The resulting pressure drop due to friction is multiplied with this factor.</li>
<li>CF_HeatTransfer: The resulting heat transfer is multiplied with this factor.</li>
</ul>
<h4>Distribution of flow resistance</h4>
<p>The parameter <code><span style=\"font-family: Courier New,courier;\">frictionDistribution</span></code> can be used to control the type of model exposed in the fluid connector. This option is available to facilitate creation of numerically sound system models with alternating control volume and flow resistance models. The available options are:</p>
<ul>
<li><code><span style=\"font-family: Courier New,courier; color: #006400;\">&quot;Volume&nbsp;-&nbsp;...&nbsp;-&nbsp;Friction&quot; </span></code>(default). Control volume model at portA and friction model at portB.</li>
<li><code><span style=\"font-family: Courier New,courier; color: #006400;\">&quot;Friction&nbsp;-&nbsp;...&nbsp;-&nbsp;Volume&quot;</span></code>. Friction model at portA and volume model at portB.</li>
<li><code><span style=\"font-family: Courier New,courier; color: #006400;\">&quot;Friction&nbsp;-&nbsp;...&nbsp;-&nbsp;Friction&quot;</span></code>. Friction models at both portA and portB.</li>
<li><code><span style=\"font-family: Courier New,courier; color: #006400;\">&quot;Volume&nbsp;-&nbsp;...&nbsp;-&nbsp;Volume&quot;</span></code>. Volume models at both portA and portB.</li>
</ul>
<p><br>The number of control volumes do not change with the different options above and neither does the number of parameters that have to be set. In the two first options, there are equal numbers of flow resistance and control volume instances in the model, so the geometric parameters are propagated directly from the segment to the flow resistance models. In the two last options there is one more or one less friction model than the number of control volumes.</p>
<p>In the case there is one more friction model, the first friction model will use the average area and hydraulic diameter and account for a factor 1/(n + 1) of the length, whereas the remaining friction models will account for the same relative amount of the remaining length as they do according to the parameters.</p>
<p>In the case there is one less friction model, the entered parameter value for area and hydraulic diameter for the last segment are ignored. The length entered for the last segment is accounted for by the remaining segments in such a way that they retain the relative length per segment.</p>
<p><b>Note:</b> Due to the handling of geometry parameters described above, option 3 and 4 should be avoided if the cross section area is not invariate along the pipe, the segments are not of equal lengths, or the segments to not include equal number of channels.</p>
</html>", revisions="<html>Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.</html>"));
end DistributedChannel;
