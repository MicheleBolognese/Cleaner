within PEMFCModel.Pipes;

model DistributedTwoPhase "Distributed channel for two phase media"

  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.Channel(
  summary(
    dp=pA - pB,
    m_flow=m_flow_mean,
    d_in=Medium.density(state[1]),
    T_in=TA,
    T_out=Medium.temperature_phX(pB,volume[end].h,volume[end].X),
    h_in=hA,
    h_out=volume[end].h,
    M=Mtot,
    V=Vtot),
    redeclare
      .Modelon.ThermoFluid.Interfaces.ApplicationSpecific.TwoPhaseVolumePort portA,
    redeclare
      .Modelon.ThermoFluid.Interfaces.ApplicationSpecific.TwoPhaseFlowPort portB,
    redeclare replaceable package Medium =
        .Modelon.Media.PreDefined.TwoPhase.WaterIF97 constrainedby
      .Modelon.Media.Interfaces.TwoPhaseMedium);
  extends .Modelon.ThermoFluid.FlowChannels.Interfaces.DistributedGeometrySizing(
    L_start = sum(L),
    A_heat_start = sum(A_heat));

  .Modelon.ThermoFluid.Interfaces.TempHeatPort[n] q_fluid
    "Connector exposing fluid temperature"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  replaceable model Friction =
      .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction constrainedby
    .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.Interfaces.TwoPhaseResistance
    "Friction model"   annotation(choicesAllMatching=true, Dialog(group="Flow resistance"));

  input Real CF_PressureLoss=1.0 "Calibration factor for pressure drop"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));

  /* Heat Transfer*/
  parameter Boolean adiabatic = false "Disable wall connector heat transfer if true"
    annotation(Evaluate = true, Dialog(group="Heat transfer"));
  replaceable model HeatTransfer =
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.CoefficientsPerPhase
    constrainedby
    .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.Interfaces.TwoPhaseHeatTransfer
    "Definition of heat transfer coefficient" annotation(choicesAllMatching=true,Dialog(group="Heat transfer", enable = not adiabatic));

  input Real CF_HeatTransfer=1.0 "Calibration factor for heat transfer"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));

  parameter Boolean includeStaticHead = false "Consider static head"
    annotation(Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Momentum Balance"));
  parameter .Modelica.Units.SI.Length[n] levels=fill(0, n)
    "Relative levels of control volume outlets (inlet port at 0)" annotation (
      Dialog(
      enable=includeStaticHead,
      tab="Balance Equations",
      group="Momentum Balance"));
  parameter .Modelica.Units.SI.Acceleration g=.Modelica.Constants.g_n
    "Acceleration of gravity" annotation (Dialog(
      enable=includeStaticHead,
      tab="Balance Equations",
      group="Momentum Balance"));

  /* Advanced parameters */
  parameter .Modelica.Units.SI.MassFlowRate mflow_zero=0
    "Region around zero flow where upstream properties are considered constant"
    annotation (Dialog(tab="Advanced", group="Numerics"));

  parameter .Modelica.Units.SI.TemperatureDifference dT_smooth_frac=2
    "Temperature smoothing interval to include multi-phase HTC"
    annotation (Dialog(tab="Advanced", group="Numerics"));

  parameter Boolean flow_based_fractions_for_dp=false
    "If true use flow rate to determine fractions of phases per segment in dp calculation"
    annotation (Dialog(tab="Advanced", group="Numerics"));

  parameter .Modelon.ThermoFluid.Choices.FrictionDistribution frictionDistribution=
    .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric
    "Discretization scheme for friction and control volumes (see info)" annotation(Evaluate=true, Dialog(tab="Advanced",group="Numerics"));
  final parameter .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal upstreamhAtReversal = .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont
    "Smoothness of upstream enthalpy at flow reversal" annotation(Evaluate=true, Dialog(tab="Advanced",group="Numerics", enable = not positiveFlow));
  final parameter .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal upstreamxAtReversal = .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont
    "Smoothness of upstream vapor quality at flow reversal" annotation(Evaluate=true, Dialog(tab="Advanced",group="Numerics", enable = not positiveFlow));

  final parameter Integer n_fric = if (frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
    frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol) then n
    elseif frictionDistribution ==
    .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then n+1
  else n-1 annotation(Evaluate=true);

  parameter Real flowFraction0 = 1 "Fraction of nominal mass flow rate through this component"
    annotation(Dialog(tab="Advanced", group="Interpretation of nominal parameters in pressure drop correlations"));
  parameter Real pos_rel_in = 0 "Relative position of channel inlet in flow direction (for correlations parameterized with inlet and outlet conditions)"
    annotation(Dialog(tab="Advanced", group="Interpretation of nominal parameters in pressure drop correlations"));
  parameter Real pos_rel_out = 1 "Relative position of channel outlet in flow direction (for correlations parameterized with inlet and outlet conditions)"
    annotation(Dialog(tab="Advanced", group="Interpretation of nominal parameters in pressure drop correlations"));
  parameter Boolean useMeanTempDrivenQ=false
    "If true, heat flow in each control volume is driven by the average temperature" annotation(Evaluate=true,Dialog(tab="Advanced", group="Heat transfer"));
  parameter Boolean includeAcceleration=false
    "If true, includes the acceleration due to specific volume in momentum balance equation"
    annotation (Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Momentum Balance"));
  parameter Boolean kineticEnergyInBalance=false
    "If true, includes the kinetic energy in energy balance"
    annotation (Evaluate=true,Dialog(tab="Balance Equations",
                                       group="Energy Balance"));

  parameter .Modelon.Types.ThermoStates stateChoice=.Modelon.Types.ThermoStates.ph_States
    "State selection"
    annotation (Dialog(tab="Advanced", group="Numerics"), Evaluate=true);
  parameter Boolean dpAsState = false
      "Prefer pressure differences as states over absolute pressures (can help low flow / low dp simulation)"
    annotation (Dialog(tab="Advanced", group="Numerics"), Evaluate=true);

  /* Length calibration */
  input Real CF_length=1.0 "Calibration factor for pipe length"
    annotation (Dialog(tab="Advanced",
  group="Calibration factors"));
  .Modelica.Units.SI.Length[n] L_internal=if not settings_TF.usePbS then L*CF_length elseif not sizing or length_fixed then fill(L_par/n,n)*CF_length else fill(L_set/n,n)
    "Modified pipe length for design calculation";
  .Modelica.Units.SI.Volume[n] V_internal=V*CF_length
    "Modified pipe volume for design calculation";
  .Modelica.Units.SI.Area[n] A_heat_internal=if not settings_TF.usePbS then A_heat*CF_length elseif not sizing or A_heat_fixed then fill(A_heat_par/n,n)*CF_length else fill(A_heat_set/n,n)
    "Modified heat transfer area for design calculation";

  Friction friction(
    redeclare package Medium = Medium,
    final n=n_fric,
    final n_channels=n_channels_fric,
    final F_user = CF_PressureLoss,
    final A=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {sum(A)/n}, A)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        A[1:n-1]
      else A),
    final L=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {sum(L_internal)/(n+1)}, L_internal*(n/(n+1)))
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        sum(L_internal) * L_internal[1:n-1]/sum(L_internal[1:n-1])
      else L_internal),
    final Dhyd=(if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, {sum(Dhyd)/n}, Dhyd)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        Dhyd[1:n-1]
      else Dhyd),
    final m_flow(start=m_flow_start*ones(n_fric)./n_channels_fric)=
      m_flow_fric./n_channels_fric,
    final d = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, d, {dB}, {dB})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
      d else cat(1, d, {dB}),
    final d_up_one = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, d_up_one, {d_up_one[n]})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
      d_up_one[1:n_fric] else d_up_one,
    final eta = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, eta, {etaB}, {etaB})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
      eta else cat(1, eta, {etaB}),
    final eta_up_one =  if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, eta_up_one, {eta_up_one[n]})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
      eta_up_one[1:n_fric] else eta_up_one,
    final onePhaseOutFraction = onePhaseOutFraction_fric,
    final twoPhaseFraction = twoPhaseFraction_fric,
    final sigma = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, sigma,{sigmaB}, {sigmaB})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
      sigma else cat(1, sigma,{sigmaB}),
    final sat = if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      cat(1, sat, {satB}, {satB})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        sat else cat(1, sat,{satB}),
    final from_dp = if not settings_TF.usePbS then from_dp else false,
    final positiveFlow = if not settings_TF.usePbS then positiveFlow else true,
    final dp_smooth = dp_smooth,
    final mflow_smooth = mflow_smooth,
    final flowFraction0 = flowFraction0,
    final pos_rel_in = pos_rel_in,
    final pos_rel_out = pos_rel_out,
    dp(start=if frictionDistribution ==
        .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then p_start[1
        :n] - p_start[2:n + 1] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n - 1] - p_start[2:n]) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n] - p_start[2:n + 1]) else p_start[1:n - 1] - p_start[2:n],
        stateSelect = if dpAsState then StateSelect.prefer else StateSelect.default,
        each nominal=1e5)) "Flow resistance model";

  HeatTransfer htcoeff(
    final n=n,
    redeclare package Medium = Medium,
    final F_user = CF_HeatTransfer,
    final A=A,
    final A_heat=if not settings_TF.usePbS then A_heat else fill(A_heat_set/n,n),
    final Dhyd=Dhyd,
    pcrit=Medium.criticalPressure,
    final L=if not settings_TF.usePbS then sum(L) elseif not sizing or length_fixed then L_par else L_set,
    final m_flow=0.5*(m_flow[1:n]+m_flow[2:n+1])./n_channels,
    final T=T,
    final T_wall=q.T,
    final p=p,
    final cp=cp,
    final Re=Re,
    final eta=eta,
    final lam=lam,
    final cp_up_one=  cp_up_one,
    final Re_up_one= Re_up_one,
    final eta_up_one=  eta_up_one,
    final lam_up_one= lam_up_one,
    final Fr=Fr,
    final Re_liq=Re_liq,
    final sat=sat,
    final onePhaseOutFraction=onePhaseOutFraction,
    final twoPhaseFraction=twoPhaseFraction,
    final CF_length=CF_length) "Heat transfer model";

  Medium.VolumeDynamics[n] volume(
    each final quasiStatic = settings_TF.usePbS,
    each thermalDynamics = true,
    each initOpt = initOpt,
    each stateChoice = stateChoice,
    p_start = p_start[1:n],
    T_start = if initFromEnthalpy then Medium.temperature_phX(p_start[1:n], h_start[2:n+1], Medium.reference_X) else T_start[2:n+1],
    h_start = if initFromEnthalpy then h_start[2:n+1] else Medium.specificEnthalpy_pTX(p_start[1:n],T_start[2:n+1],Medium.reference_X),
    each X_start = X_start,
    each C_start = C_start,
    dE = dE,
    dMX = transpose(dMX),
    dMC = transpose(dMC),
    dM_bulk = m_flow[1:n] - m_flow[2:n+1],
    final V_tot = V_internal,
    final T_in=T) "Control volume models";

  .Modelica.Units.SI.Mass Mtot "Total mass";
  .Modelica.Units.SI.Volume Vtot "Total volume";
  .Modelica.Units.SI.HeatFlowRate Q_tot "Total heat flow rate into pipe";

  .Modelica.Units.SI.AbsolutePressure pA "Pressure in portA";
  .Modelica.Units.SI.AbsolutePressure pB "Pressure in portB";
  .Modelon.Media.Units.Temperature TA "Upstream temperature when flow A -> B";
  .Modelon.Media.Units.Temperature TB "Upstream temperature when flow B -> A";
  .Modelon.Media.Units.Temperature TA_out "Outlet temperature when flow B -> A";
  .Modelon.Media.Units.Temperature TB_out "Outlet temperature when flow A -> B";
  Medium.ThermodynamicState[n+2] state
    "Thermodynamic states in control volumes and instream from ports";
  Medium.ThermodynamicState stateA_out
    "Thermodynamic state flowing out of portA (in case of back flow)";
  Medium.ThermodynamicState stateB_out
    "Thermodynamic state flowing out of portB (in case of forward flow)";

  .Modelica.Units.SI.Mass[n] M "Control volume masses";

  .Modelica.Units.SI.Pressure[n] sh "Static head between control volumes";
  .Modelica.Units.SI.AbsolutePressure[n] p(start=p_start[1:n])
    "Control volume pressures";
  .Modelica.Units.SI.Temperature[n] T "Control volume temperatures";
  input Real[n] x(each unit="1") = {noEvent(if p_red[i] < 1.0 then (h[i] - h_liq[i])/max(h_vap[i] - h_liq[i], 1e-6) else 1.0) for i in 1:n}
    "Equilibrium thermodynamic quality (<0 for subcooled, >1 for superheated)";
  .Modelica.Units.SI.MassFraction[n] quality "Steam quality (0 <= quality <= 1)";
  .Modelica.Units.SI.SpecificEnthalpy[n] h(start=volume.h_start) "Control volume specific enthalpies";
  .Modelica.Units.SI.MassFlowRate[n + 1] m_flow(each start=m_flow_start)
    "Mass flow over cv boundaries, positive from portA towards portB";
  .Modelica.Units.SI.MassFlowRate m_flow_mean
    "Average mass flow rate in pipe, positive from portA to portB";
  .Modelica.Units.SI.Pressure[n_fric] dp(start=if frictionDistribution ==
        .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then p_start[1
        :n] - p_start[2:n + 1] elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n - 1] - p_start[2:n]) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n] - p_start[2:n + 1]) else p_start[1:n - 1] - p_start[2:n],
        stateSelect = if dpAsState then StateSelect.prefer else StateSelect.default,
        each nominal=1e5)
    "Pressure drops from friction models, vector length depending on discretization option";

  input .Modelica.Units.SI.Density[n] d=volume.d;

  .Modelica.Units.SI.HeatFlowRate[n] Q "Heat flow rate into control volumes";
  .Modelica.Units.SI.CoefficientOfHeatTransfer[n] alpha
    "Heat transfer coefficient";

  .Modelica.Units.SI.Temperature T_fluid[n]
    "Temperature for heat transfer modeling";

  .Modelica.Units.SI.MassFlowRate mflow_A_in "Total inlet flow rate at portA";
  .Modelica.Units.SI.MassFlowRate mflow_B_in "Total inlet flow rate at portB";

protected
  .Modelica.Units.SI.MassFlowRate[n_fric] m_flow_fric=
    if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then
      m_flow[2:n + 1]
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then
      m_flow[1:n]
    elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
      m_flow
    else m_flow[2:n]
    "Mass flow over cv boundaries, positive from portA towards portB";

   final parameter Real[n_fric] n_channels_fric=
     if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
       cat(1, {n_channels[1]}, n_channels)
     elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
       n_channels[1:n-1]
     else n_channels;

public
  .Modelica.Units.SI.Pressure[n + 1] dp_internal(start=if frictionDistribution
         == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then cat(
        1,
        {0},
        p_start[1:n] - p_start[2:n + 1]) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n - 1] - p_start[2:n],
        {0}) elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric
         then cat(
        1,
        {(p_start[1] - p_start[n + 1])/n_fric},
        p_start[1:n] - p_start[2:n + 1]) else cat(
        1,
        {0},
        p_start[1:n - 1] - p_start[2:n],
        {0})) "Pressure drop after each control volume";

  .Modelica.Units.SI.SpecificEnthalpy hA
    "Mixture enthalpy of inlet flows at portA";
  .Modelica.Units.SI.SpecificEnthalpy hB
    "Mixture enthalpy of inlet flows at portB";

  input .Modelica.Units.SI.VolumeFraction[n] onePhaseOutFractionZero=
    {.Modelon.ThermoFluid.Functions.onePhaseOutFraction(
        {h[i], h[i]},
        {h_liq[i], h_liq[i]},
        {h_vap[i], h_vap[i]},
        p_red[i], 10)  for i in 1:n}
    "Fraction of control volume in single phase region at zero flow";

  .Modelica.Units.SI.VolumeFraction[n] onePhaseOutFraction
    "Fraction of control volume in single phase region, outlet side";

  input .Modelica.Units.SI.VolumeFraction[n] twoPhaseFractionZero=
      {max(0, 1-onePhaseOutFractionZero[i]) for i in 1:n}
      "Fraction of control volume in two phase region at zero flow";

  .Modelica.Units.SI.VolumeFraction[n] twoPhaseFraction
    "Fraction of control volume in two-phase region";
  .Modelica.Units.SI.VolumeFraction[n] onePhaseInFraction
    "Fraction of control volume in single phase region, inflow";

protected

  .Modelica.Units.SI.VolumeFraction[n_fric] onePhaseOutFraction_fric=
    if flow_based_fractions_for_dp then
      (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, onePhaseOutFraction[1:n-1],
          {max(0, 2*(onePhaseOutFraction[n]-0.5))},
          {min(1, 2*onePhaseOutFraction[n])})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        cat(1, onePhaseOutFraction[1:n-2],
          {0.5*(onePhaseOutFraction[n-1]+onePhaseOutFraction[n])})
      else
        onePhaseOutFraction)
      else
      (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, onePhaseOutFractionZero[1:n-1],
          {max(0, 2*(onePhaseOutFractionZero[n]-0.5))},
          {min(1, 2*onePhaseOutFractionZero[n])})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        cat(1, onePhaseOutFractionZero[1:n-2],
          {0.5*(onePhaseOutFractionZero[n-1]+onePhaseOutFractionZero[n])})
      else
        onePhaseOutFractionZero);

   .Modelica.Units.SI.VolumeFraction[n_fric] twoPhaseFraction_fric=
    if flow_based_fractions_for_dp then
      (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, twoPhaseFraction[1:n-1],
          { 2*(max(0,twoPhaseFraction[n]+onePhaseOutFraction[n] - 0.5)
           - max(0, onePhaseOutFraction[n]-0.5))},
          { 2*max(0, min(0.5, twoPhaseFraction[n]) - onePhaseOutFraction[n])})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        cat(1, twoPhaseFraction[1:n-2],
          {0.5*(twoPhaseFraction[n-1]+twoPhaseFraction[n])})
      else
        twoPhaseFraction)
      else
      (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        cat(1, twoPhaseFractionZero[1:n-1],
          { 2*(max(0,twoPhaseFractionZero[n]+onePhaseOutFractionZero[n] - 0.5)
           - max(0, onePhaseOutFractionZero[n]-0.5))},
          { 2*max(0, min(0.5, twoPhaseFractionZero[n]) - onePhaseOutFractionZero[n])})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        cat(1, twoPhaseFractionZero[1:n-2],
          {0.5*(twoPhaseFractionZero[n-1]+twoPhaseFractionZero[n])})
      else
        twoPhaseFractionZero);

  public
  .Modelica.Units.SI.MassFraction[Medium.nS] XA
    "Mixture mass fractions of inlet flows at portA";
  .Modelica.Units.SI.MassFraction[Medium.nS] XB
    "Mixture mass fractions of inlet flows at portB";

  .Modelica.Units.SI.Density dA "Inlet density, nominal direction, portA";
  .Modelica.Units.SI.Density dB "Inlet density, nominal direction, portB";
  .Modelica.Units.SI.Density dA_out "Outlet density, nominal direction, portA";
  .Modelica.Units.SI.Density dB_out "Outlet density at portB";
  .Modelica.Units.SI.DynamicViscosity etaA
    "Inlet viscosity, nominal direction, portA";
  .Modelica.Units.SI.DynamicViscosity etaB
    "Inlet viscosity, nominal direction, portB";
  .Modelica.Units.SI.SurfaceTension sigmaA
    "Inlet surface tension, nominal direction, portA";
  .Modelica.Units.SI.SurfaceTension sigmaB
    "Inlet surface tension, nominal direction, portB";
  Medium.PhaseBoundaryProps satA
    "Inlet saturation properties, nominal direction, portA";
  Medium.PhaseBoundaryProps satB
    "Inlet saturation properties, nominal direction, portB";
  Real p_redA "Reduced pressure: 1 at critical point, portA inlet";
  Real p_redB "Reduced pressure: 1 at critical point, portB inlet";
  .Modelica.Units.SI.AbsolutePressure p_satA(start=p_start[1])
    "Saturation pressure, portA inlet";
  .Modelica.Units.SI.AbsolutePressure p_satB(start=p_start[end])
    "Saturation pressure, portB inlet";
  .Modelon.Media.Units.Temperature T_vapA "Dew line temperature, portA inlet";
  .Modelon.Media.Units.Temperature T_vapB "Dew line temperature, portB inlet";
  .Modelon.Media.Units.Temperature T_liqA "Bubble line temperature, portA inlet";
  .Modelon.Media.Units.Temperature T_liqB "Bubble line temperature, portB inlet";
  .Modelica.Units.SI.SpecificEnthalpy h_liqA(
    start=Medium.bubbleEnthalpy_pX(min(p_start[1], Medium.criticalPressure), Medium.reference_X)) "Liquid enthalpy, portA inlet";
  .Modelica.Units.SI.SpecificEnthalpy h_vapA(
    start=Medium.dewEnthalpy_pX(min(p_start[1], Medium.criticalPressure), Medium.reference_X)) "Vapour enthalpy, portA inlet";
  .Modelica.Units.SI.SpecificEnthalpy h_liqB(
    start=Medium.bubbleEnthalpy_pX(min(p_start[n+1], Medium.criticalPressure), Medium.reference_X)) "Liquid enthalpy, portB inlet";
  .Modelica.Units.SI.SpecificEnthalpy h_vapB(
    start=Medium.dewEnthalpy_pX(min(p_start[n+1], Medium.criticalPressure), Medium.reference_X)) "Vapour enthalpy, portB inlet";

  .Modelica.Units.SI.Temperature[n] T_wall
    "Wall temperature, if heat transfer enabled";
  .Modelica.Units.SI.HeatFlowRate[n] Q_wall "Heat flow rate from connector q";
  .Modelica.Units.SI.HeatFlowRate[n] Q_fluid
    "Heat flow rate from connector q_fluid";
  .Modelica.Units.SI.HeatFlowRate[n-1] Q_cond_int
    "Internal thermal conduction";
  Real[n] p_red(start=p_start[1:n]./Medium.criticalPressure) "Reduced pressure: 1 at critical point";
  .Modelica.Units.SI.AbsolutePressure[n] p_sat(start=p_start[1:n]) "Saturation pressures";
  .Modelica.Units.SI.Temperature[n] T_vap(
    start={Medium.dewTemperature_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n}) "Dew line temperature";
  .Modelica.Units.SI.Temperature[n] T_liq(
    start={Medium.bubbleTemperature_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n}) "Bubble line temperature";
  .Modelica.Units.SI.SpecificEnthalpy[n + 1] h_liq(
    start={Medium.bubbleEnthalpy_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n+1}) "Liquid enthalpy";
  .Modelica.Units.SI.SpecificEnthalpy[n + 1] h_vap(
    start={Medium.dewEnthalpy_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n+1}) "Vapour enthalpy";
  .Modelon.Media.Units.Temperature T_vap_out "Dew line temperature";
  .Modelon.Media.Units.Temperature T_liq_out "Bubble line temperature";
  .Modelon.Media.Units.MassFraction[Medium.nS,n] X "Mass fractions";
  .Modelica.Units.SI.SpecificHeatCapacity[n] cp
    "Heat capacity at constant pressure";
  .Modelica.Units.SI.SpecificHeatCapacity[n] cv
    "Heat capacity at constant volume";
  .Modelica.Units.SI.ReynoldsNumber Re[n] "Reynolds number";
  input .Modelica.Units.SI.DynamicViscosity[n] eta(
    start={Medium.dynamicViscosity_dTX(
      Medium.density_phX(
        volume[i].p_start,
        volume[i].h_start,
        volume[i].X_start),
      volume[i].T_start,
      volume[i].X_start) for i in 1:n})
    ={Medium.dynamicViscosity_dTX(
      d[i],
      T[i],
      X[:, i]) for i in 1:n};
  .Modelica.Units.SI.ThermalConductivity lam[n] "Thermal conductivity";
  .Modelica.Units.SI.FroudeNumber Fr[n] "Froude number";
  .Modelica.Units.SI.ReynoldsNumber Re_liq[n] "Reynolds number, boiling curve";
  .Modelica.Units.SI.SurfaceTension sigma[n]
    "Surface tension in two phase region";
  Medium.PhaseBoundaryProps[n] sat(
    d_liq(start={Medium.bubbleDensity_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n}),
    d_vap(start={Medium.dewDensity_pX(min(p_start[i], Medium.criticalPressure), Medium.reference_X) for i in 1:n})) "Saturation properties";

  //Variables for detailed balance formulations
  .Modelica.Units.SI.Velocity v[n + 1]
    "Fluid velocities at control volume boundaries";

  .Modelica.Units.SI.Power Ek_flow[n + 1]
    "Kinetic energy flow m_flow*v^2/2 over control volume boundaries";

  .Modelica.Units.SI.Power dEk_flow[n]
    "Difference of kinetic energy flow over boundaries, per control volume";

  .Modelica.Units.SI.PressureDifference[n+1] dpaccel
    "Resulting pressure drop from acceleration due to specific volume change";


protected
  .Modelica.Units.SI.MassFlowRate m_flow_help[n] =
   { (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
      frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then m_flow[i+1] else m_flow[i])  for i in 1:n}
      "Flow rate used for determining upstream direction per CV";

  .Modelica.Units.SI.Length L2in_prel[n] =
    { abs(m_flow[i]+m_flow[i+1])/2* L[i] / max(1e-15, htcoeff.alpha2ph[i] * n_channels[i]*A_heat_internal[i])
     *(if noEvent(T_wall[i] >= T_vap[i]) then
        noEvent(if m_flow_help[i] >= 0 then
            max(0, (if i == n then h_vapB else h_vap[i+1]) - max((if i==1 then hA else h[i-1]), h_liq[i])) else
            max(0, (if i == 1 then h_vapA else h_vap[i-1]) - max((if i==n then hB else h[i+1]), h_liq[i])))
          / max(1e-6, T_wall[i] - T_vap[i])
    elseif noEvent(T_wall[i] <= T_liq[i]) then
        noEvent(if m_flow_help[i] >= 0 then
            min(0, (if i==n then h_liqB else h_liq[i+1]) - min((if i==1 then hA else h[i-1]), h_vap[i])) else
            min(0, (if i==1 then h_liqA else h_liq[i-1]) - min((if i==n then hB else h[i+1]), h_vap[i])))
          / min(-1e-6, T_wall[i] - T_liq[i])
      else 1e4)  for i in 1:n}
     "Length of twophase inlet flow region into single phase CV";

   .Modelica.Units.SI.Length L2in[n] =
     {if noEvent(h_up[i] > (if m_flow_help[i] >= 0 then (if i==n then h_vapB else h_vap[i+1]) else (if i==1 then h_vapA else h_vap[i-1])))
      then .Modelon.Math.Smoothing.spliceFunction(
         L2in_prel[i], 0,
         min(T[i], T_liq[i]) - T_wall[i] - dT_smooth_frac, dT_smooth_frac)
      elseif noEvent(h_up[i] < (if m_flow_help[i] >= 0 then (if i==n then h_liqB else h_liq[i+1]) else (if i==1 then h_liqA else h_liq[i-1])))
        then .Modelon.Math.Smoothing.spliceFunction(
         L2in_prel[i], 0,
         T_wall[i] - max(T[i], T_vap[i]) - dT_smooth_frac, dT_smooth_frac)
      else
         .Modelon.Math.Smoothing.spliceFunction(
         L2in_prel[i], 0,
         if noEvent(h[i]>0.5*(h_liq[i]+h_vap[i]))
           then T_wall[i] - max(T[i], T_vap[i]) - dT_smooth_frac
           else min(T[i], T_liq[i]) - T_wall[i] - dT_smooth_frac, dT_smooth_frac)
     for i in 1:n} "Smooth activation wrt Twall of phase lengths";

  .Modelica.Units.SI.Length L1in_prel[n] =
    { -abs(m_flow[i]+m_flow[i+1])/2*cp_up_one[i]* L[i] / max(1e-15, htcoeff.alpha1phIn[i] * n_channels[i]*A_heat_internal[i])
     *noEvent(if T_up_one[i] >= T_wall[i] then
     log(min(1, max(1e-4, T_vap[i] - T_wall[i])/max(1e-4, T_up_one[i] - T_wall[i])))
    else
      log(min(1, min(-1e-4, T_liq[i] - T_wall[i])/min(-1e-4, T_up_one[i] - T_wall[i]))))
        for i in 1:n}
     "Length of single phase inlet flow region into two-phase or opposite phase CV";

  .Modelica.Units.SI.Length L1in[n] =
     {if noEvent(h_up[i] > (if m_flow_help[i] >= 0 then (if i==n then h_vapB else h_vap[i+1]) else (if i==1 then h_vapA else h_vap[i-1])))
      then .Modelon.Math.Smoothing.spliceFunction(
         L1in_prel[i], 0,
         min(T[i],T_liq[i]) - T_wall[i] - dT_smooth_frac, dT_smooth_frac)
      elseif noEvent(h_up[i] < (if m_flow_help[i] >= 0 then (if i==n then h_liqB else h_liq[i+1]) else (if i==1 then h_liqA else h_liq[i-1])))
        then .Modelon.Math.Smoothing.spliceFunction(
         L1in_prel[i], 0,
         T_wall[i] - max(T[i], T_vap[i]) - dT_smooth_frac, dT_smooth_frac)
      else
         .Modelon.Math.Smoothing.spliceFunction(
         L1in_prel[i], 0,
         if noEvent(h[i]>0.5*(h_liq[i]+h_vap[i])) then T_wall[i] - max(T[i],T_vap[i]) - dT_smooth_frac
         else min(T[i],T_liq[i]) - T_wall[i] - dT_smooth_frac, dT_smooth_frac)
     for i in 1:n} "Smooth activation wrt Twall of phase lengths";

  .Modelica.Units.SI.HeatFlowRate[n] Q_sens_zero "Zero flow heat flow rate from connector q";
  .Modelica.Units.SI.HeatFlowRate[n] Q_sens "Heat flow rate from connector q";
  .Modelica.Units.SI.HeatFlowRate[n] Q_latRateLim "Max transferable latent heat based on area, HTC, dT";
  //Modelica.Units.SI.HeatFlowRate[n] Q_latFluidLim "Max transferable latent heat based on area, HTC, dT";
  .Modelica.Units.SI.HeatFlowRate[n] Q_lat "Actual latent heat";

  .Modelica.Units.SI.MassFlowRate[Medium.nS,n] dMX
    "Mass derivative of substances";
  Real[Medium.nC,n] dMC "Trace component derivatives";
  .Modelica.Units.SI.HeatFlowRate[n] dE "Internal energy derivative";
  .Modelica.Units.SI.EnthalpyFlowRate H_flow[n + 1]
    "Enthalpy flow over control volume boundaries";
  .Modelica.Units.SI.MassFlowRate MX_flow[Medium.nS,n + 1]
    "Species mass flow rate over cv boundaries";
  Real MC_flow[Medium.nC,n + 1] "Trace components flow rate over cv boundaries";

  //Real[n] x_up(each unit="1") "Upstream vapor quality, flow direction dependent";
  .Modelica.Units.SI.SpecificEnthalpy[n] h_up "Flow dir dependent upstream specific enthalpies";
  .Modelica.Units.SI.Temperature[n] T_up_one "Temperature of upstream single phase fraction";
  .Modelica.Units.SI.Density[n] d_up_one "Density of upstream single phase fraction";
  .Modelica.Units.SI.DynamicViscosity[n] eta_up_one "Dynamic viscosity of upstream single phase fraction";
    .Modelica.Units.SI.SpecificHeatCapacity[n] cp_up_one
    "Heat capacity at constant pressure of upstream single phase fraction";
  .Modelica.Units.SI.ReynoldsNumber Re_up_one[n] "Reynolds number of upstream single phase fraction";
  .Modelica.Units.SI.ThermalConductivity lam_up_one[n] "thermal conductivity, upstream single phase";

  parameter Real par=0 "Auxiliary parameter" annotation (Evaluate=false);
  // Constants for division to bring to dimensionless ratios in held equation
  constant .Modelica.Units.SI.Pressure  dp_1 = 1 "Unity mass flow rate";
  constant .Modelica.Units.SI.Length L_1 = 1 "Unity pressure";
  constant .Modelica.Units.SI.HeatFlowRate Q_flow_1 = 1 "Unity heat flow rate";
  constant .Modelica.Units.SI.Area A_heat_1 = 1 "Unity heat transfer area";
  constant .Modelon.Media.Units.Temperature T_out_1 = 1 "Unity temperature";
initial equation

  // Non-PbS equations: Prescribe L_par, the total length if not settings_TF.usePbS
  if not settings_TF.usePbS then
    // Define length via prescribed pressure drop or directly
    if sizing and not length_fixed then
      sum(dp_internal) = dp0 + par*L_par/L_1;
    else
      L_par = L_total;
    end if;
    if sizing and not A_heat_fixed then
      if thermalOpt == .Modelon.ThermoFluid.Choices.ChannelThermalOpt.Q_flow then
        Q_tot = Q_flow0 + par*A_heat_par/A_heat_1;
      else
        T[end] = T_out0 + par*A_heat_par/A_heat_1;
      end if;
    else
      A_heat_par = sum(A_heat);
    end if;
  else
    // Define directly as dummy value
    L_par = L_total;
    A_heat_par = sum(A_heat);
  end if;
equation
  /* Aggregate quantities*/
  Mtot = sum(M) "Total internal fluid mass";
  Vtot=sum(V_internal) "Total volume";
  Q_tot = sum(Q) "Total heat transfer into channel";
  m_flow_mean = sum(m_flow)/(n+1) "Mass flow rate averaged over channel length";

  /* Boundary conditions */
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

  if adiabatic then
    T_wall = T;
  else
    q.Q_flow = Q_wall;
  end if;

  /* Inlet flow properties in vectorized ports */
  mflow_A_in = sum({max(0, portA[i].m_flow) for i in 1:NA})
    "Sum of inlet flows, portA";
  mflow_B_in = sum({max(0, portB[i].m_flow) for i in 1:NB})
    "Sum of inlet flows, portB";

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

  /* Populate vector of thermodynamic states in ports and control volumes */
  //State records from port intet flows
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
  //State records from control volume dynamic states
  state[2:n+1] = volume.state;

  // Definition of dp_internal vector as differences between ports and control volumes

  dp_internal = if includeStaticHead then
    cat(1,{pA-volume[1].p}, volume[1:n-1].p-volume[2:n].p-sh[1:n-1], {volume[n].p-pB-sh[n]})
    else
    cat(1,{pA-volume[1].p}, volume[1:n-1].p-volume[2:n].p, {volume[n].p-pB});

  sh = cat(1, {levels[1]}, levels[2:n] - levels[1:(n-1)]) .* d * g;

    // Definition of drhodx  as differences between ports and control volumes density
    // Definition of dpaccel as pressure drop due to the change in specific volume

  dpaccel = if includeAcceleration then
    cat(1,
    {d[1]*abs(.Modelon.Math.Smoothing.regExp(m_flow[1]/(d[1]*A[1]*n_channels[1]), 0.0001, 2)) -
     noEvent(if m_flow[1] >= 0 then dA else dA_out)*
       abs(.Modelon.Math.Smoothing.regExp(m_flow[1]/(noEvent(if m_flow[1] >= 0 then dA else dA_out)*A[1]*n_channels[1]), 0.0001, 2))},
    {d[i]*abs(.Modelon.Math.Smoothing.regExp(m_flow[i]/(d[i]*A[i]*n_channels[i]), 0.0001, 2)) -
      d[i-1]*abs(.Modelon.Math.Smoothing.regExp(m_flow[i]/(d[i-1]*A[i]*n_channels[i]), 0.0001, 2))  for i in 2:n},
    {noEvent(if m_flow[n+1] >= 0 then dB_out else dB)
      *abs(.Modelon.Math.Smoothing.regExp(m_flow[n+1]/(noEvent(if m_flow[n+1] >= 0 then dB_out else dB)*A[n]*n_channels[n]), 0.0001, 2)) -
     d[n]*abs(.Modelon.Math.Smoothing.regExp(m_flow[n+1]/(d[n]*A[n]*n_channels[n]), 0.0001, 2))})
    else fill(0, n+1);

  /* Pressure drop to friction model */
  dp = friction.dp;
  dp_internal = (if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric then
        cat(1, {0}, dp)
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricVol then
        cat(1, dp, {0})
      elseif frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.FricFric then
        dp else cat(1, {0}, dp, {0}))
    + dpaccel;
  /* Fluid properties in control volumes */
  M = d .* V_internal;
  p = volume.p;
  h = volume.h;
  T = volume.T;
  p_red = p/Medium.criticalPressure;

  dA = Medium.density(state[1]);
  dB = Medium.density(state[n+2]);
  dA_out = Medium.density(stateA_out);
  dB_out = Medium.density(stateB_out);

  etaA = Medium.dynamicViscosity_dTX(dA, TA, XA);
  etaB = Medium.dynamicViscosity_dTX(dB, TB, XB);
  p_redA = pA/Medium.criticalPressure;
  p_redB = pB/Medium.criticalPressure;
  p_satA = min(pA, Medium.criticalPressure);
  p_satB = min(pB, Medium.criticalPressure);
  T_vapA = Medium.dewTemperature_pX(p_satA, XA);
  T_vapB = Medium.dewTemperature_pX(p_satB, XB);
  T_liqA = Medium.bubbleTemperature_pX(p_satA, XA);
  T_liqB = Medium.bubbleTemperature_pX(p_satB, XB);
  h_liqA = Medium.bubbleEnthalpy_pX(p_satA, XA);
  h_vapA = Medium.dewEnthalpy_pX(p_satA, XA);
  h_liqB = Medium.bubbleEnthalpy_pX(p_satB, XB);
  h_vapB = Medium.dewEnthalpy_pX(p_satB, XB);
  h_vap[n+1] = Medium.dewEnthalpy_pX(pB, X[:,n]);
  h_liq[n+1] = Medium.bubbleEnthalpy_pX(pB, X[:,n]);
  T_liq_out = Medium.bubbleTemperature_pX(pB, X[:,n]);
  T_vap_out = Medium.dewTemperature_pX(pB, X[:,n]);

  sigmaA = Medium.surfaceTension_pX(p_satA, XA);
  sigmaB = Medium.surfaceTension_pX(p_satB, XB);

  satA.d_vap = Medium.dewDensity_pX(p_satA, XA);
  satA.d_liq = Medium.bubbleDensity_pX(p_satA, XA);
  satA.eta_vap = Medium.dynamicViscosity_dTX(satA.d_vap, T_vapA, XA);
  satA.eta_liq = Medium.dynamicViscosity_dTX(satA.d_liq, T_liqA, XA);
  satA.cp_vap = min(max(Medium.specificHeatCapacityCp_phX(p_satA, h_vapA, XA), 500),5e5);
  satA.cp_liq = min(max(Medium.specificHeatCapacityCp_phX(p_satA, h_liqA, XA), 500),5e5);
  satA.lam_vap = min(max(Medium.thermalConductivity_dTX(satA.d_vap, T_vapA, XA), 1.0e-4),2.0);
  satA.lam_liq = min(max(Medium.thermalConductivity_dTX(satA.d_liq, T_liqA, XA), 1.0e-4),2.0);
  satA.x = noEvent(if p_redA < 1.0 then (hA - h_liqA)/max(h_vapA - h_liqA, 1e-6) else 1.0);

  satB.d_vap = Medium.dewDensity_pX(p_satB, XB);
  satB.d_liq = Medium.bubbleDensity_pX(p_satB, XB);
  satB.eta_vap = Medium.dynamicViscosity_dTX(satB.d_vap, T_vapB, XB);
  satB.eta_liq = Medium.dynamicViscosity_dTX(satB.d_liq, T_liqB, XB);
  satB.cp_vap = min(max(Medium.specificHeatCapacityCp_phX(p_satB, h_vapB, XB), 500),5e5);
  satB.cp_liq = min(max(Medium.specificHeatCapacityCp_phX(p_satB, h_liqB, XB), 500),5e5);
  satB.lam_vap = min(max(Medium.thermalConductivity_dTX(satB.d_vap, T_vapB, XB), 1.0e-4),2.0);
  satB.lam_liq = min(max(Medium.thermalConductivity_dTX(satB.d_liq, T_liqB, XB), 1.0e-4),2.0);
  satB.x = noEvent(if p_redB < 1.0 then (hB - h_liqB)/max(h_vapB - h_liqB, 1e-6) else 1.0);

  cp = Medium.specificHeatCapacityCp(volume.state);
  cv = Medium.specificHeatCapacityCv(volume.state);
  for i in 1:n loop
    X[:,i] = volume[i].X;
    p_sat[i] = min(p[i], Medium.criticalPressure);
    sigma[i] = Medium.surfaceTension_pX(p_sat[i],X[:,i]);
    lam[i] = min(max(smooth(0, noEvent(if p_red[i] < 1.0 then
           .Modelon.Math.Smoothing.spliceFunction((1.0 - sat[i].x)*sat[i].lam_liq + sat[i].x*sat[i].lam_vap,
                                             Medium.thermalConductivity_dTX(d[i], T[i],X[:,i]),
                                             min(sat[i].x+0.02,1.0-sat[i].x),0.02) else
                            Medium.thermalConductivity_dTX(d[i], T[i],X[:,i]))),1e-4),2.0);
    T_vap[i] = Medium.dewTemperature_pX(p_sat[i],X[:,i]);
    T_liq[i] = Medium.bubbleTemperature_pX(p_sat[i],X[:,i]);
    h_vap[i] = Medium.dewEnthalpy_pX(p_sat[i], X[:,i]);
    h_liq[i] = Medium.bubbleEnthalpy_pX(p_sat[i], X[:,i]);
    sat[i].d_vap = Medium.dewDensity_pX(p_sat[i], X[:,i]);
    sat[i].d_liq = Medium.bubbleDensity_pX(p_sat[i], X[:,i]);
    sat[i].eta_vap = Medium.dynamicViscosity_dTX(sat[i].d_vap, T_vap[i], X[:,i]);
    sat[i].eta_liq = Medium.dynamicViscosity_dTX(sat[i].d_liq, T_liq[i], X[:,i]);
    // cp and lambda go haywire at the critical point: safeguards needed
    sat[i].cp_vap = min(max(Medium.specificHeatCapacityCp_phX(p_sat[i], h_vap[i], X[:,i]), 500),5e5);
    sat[i].cp_liq = min(max(Medium.specificHeatCapacityCp_phX(p_sat[i], h_liq[i], X[:,i]), 500),5e5);
    sat[i].lam_vap = min(max(Medium.thermalConductivity_dTX(sat[i].d_vap, T_vap[i], X[:,i]), 1.0e-4),2.0);
    sat[i].lam_liq = min(max(Medium.thermalConductivity_dTX(sat[i].d_liq, T_liq[i], X[:,i]), 1.0e-4),2.0);

    sat[i].x = x[i]
      "For correlation use";
    quality[i] = max(0.0, min(1.0, x[i]));

    /* Flow regime */
    Re[i] = .Modelon.Math.Smoothing.above(
        .Modelon.ThermoFluid.Functions.CharacteristicNumbers.ReynoldsNumber(
        0.5*(m_flow[i]+m_flow[i+1])/n_channels[i], Dhyd[i], A[i], eta[i]), 1, 10);
    Re_liq[i] = .Modelon.Math.Smoothing.above(
        .Modelon.ThermoFluid.Functions.CharacteristicNumbers.ReynoldsNumber(
        0.5*(m_flow[i]+m_flow[i+1])/n_channels[i], Dhyd[i], A[i], sat[i].eta_liq), 1, 10);
    Fr[i] = .Modelon.Math.Smoothing.above(
        .Modelon.ThermoFluid.Functions.CharacteristicNumbers.FroudeNumber(
        .Modelon.Math.Smoothing.smoothMax(abs(0.5*(m_flow[i]+m_flow[i+1])), mflow_smooth, mflow_smooth)/n_channels[i], Dhyd[i], A[i], sat[i].d_liq), 1, 5);
  end for;

  T_fluid = if useMeanTempDrivenQ then
    (if n>1 then
      cat(1, {.Modelon.Math.Smoothing.spliceFunction((TA+volume[1].T)/2, (volume[1].T+volume[2].T)/2, m_flow[1], mflow_smooth)},
        {.Modelon.Math.Smoothing.spliceFunction((volume[i-1].T+volume[i].T)/2, (volume[i].T+volume[i+1].T)/2, m_flow[i], mflow_smooth) for i in 2:n-1},
        {.Modelon.Math.Smoothing.spliceFunction((volume[n-1].T+volume[n].T)/2, (volume[n].T+TB)/2, m_flow[n], mflow_smooth)})
      else
        fill(.Modelon.Math.Smoothing.spliceFunction((TA+volume[1].T)/2, (volume[1].T+TB)/2, m_flow[1], mflow_smooth),n))
    else
      volume.T;

  /* Properties flow over control volume boundaries */
  H_flow[1] = sum({portA[i].m_flow*actualStream(portA[i].h_outflow) for i in 1:NA});
  H_flow[n+1] = -sum({portB[i].m_flow*actualStream(portB[i].h_outflow) for i in 1:NB});
  MX_flow[:,1] = {sum({portA[i].m_flow*actualStream(portA[i].X_outflow[k]) for i in 1:NA}) for k in 1:Medium.nS};
  if not settings_TF.usePbS then
    MX_flow[:,n+1] = {-sum({portB[i].m_flow*actualStream(portB[i].X_outflow[k]) for i in 1:NB}) for k in 1:Medium.nS};
  else
    MX_flow[:,n+1] = {-sum({portB[i].m_flow*portB[i].X_outflow[k] for i in 1:NB}) for k in 1:Medium.nS};
  end if;
  MC_flow[:,1] = {sum({portA[i].m_flow*actualStream(portA[i].C_outflow[k]) for i in 1:NA}) for k in 1:Medium.nC};
  MC_flow[:,n+1] = {-sum({portB[i].m_flow*actualStream(portB[i].C_outflow[k]) for i in 1:NB}) for k in 1:Medium.nC};

  // Calculating the kinetic energy
  v = cat(1, {noEvent(if m_flow[1] >= 0 then m_flow[1]/(dA*A[1]*n_channels[1]) else
    m_flow[1]/(dA_out*A[1]*n_channels[1]))},
       {m_flow[i]/(A[i]*n_channels[i]*noEvent(if m_flow[i] >= 0 then d[i - 1] else d[i])) for i in 2:n},
       {noEvent(if m_flow[n+1] >= 0 then
          m_flow[n+1]/(dB_out*A[n]*n_channels[n]) else
          m_flow[n+1]/(dB*A[n]*n_channels[n]))});

  Ek_flow = if kineticEnergyInBalance then
    cat(1, {(abs(m_flow[1])*.Modelon.Math.Smoothing.regExp(v[1], 0.0001, 2))/2},
      {(abs(m_flow[i])*.Modelon.Math.Smoothing.regExp(v[i], 0.0001, 2))/2 for i in 2:n},
      {(abs(m_flow[n+1])*
        .Modelon.Math.Smoothing.regExp(v[n + 1], 0.0001, 2))/2})
      else
        fill(0,n+1);

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

  dEk_flow = Ek_flow[1:n] - Ek_flow[2:n + 1];

  /* Net mass and enthalpy flow into control volumes */
  if not settings_TF.usePbS then
    dE =H_flow[1:n] - H_flow[2:n + 1] + Q + dEk_flow
     - cat(1, Q_cond_int, {0}) + cat(1, {0}, Q_cond_int);
  else
    volume[1].h = sum({portA[i].m_flow*inStream(portA[i].h_outflow) for i in 1:NA})/m_flow[1]
                     + (Q[1]+ volume[1].dE)/m_flow[1]
        annotation(__Modelon(ResidualEquation(iterationVariable=h[1])));
    for i in 2:n loop
      volume[i].h = volume[i-1].h + (Q[i] + volume[i].dE)/m_flow[i]
        annotation(__Modelon(ResidualEquation(iterationVariable=h[i])));
    end for;
  end if;
  dMX = MX_flow[:,1:n] - MX_flow[:,2:n+1];
  dMC = MC_flow[:,1:n] - MC_flow[:,2:n+1];

  /* Heat transfer */
  Q = Q_wall + Q_fluid;
  alpha = htcoeff.alpha;
  for i in 1:n loop
    Q_wall[i] = Q_sens[i] + Q_lat[i];
    Q_sens[i] = onePhaseOutFraction[i]*Q_sens_zero[i];
    Q_sens_zero[i] = if adiabatic then 0 else
      n_channels[i]*A_heat_internal[i]*htcoeff.alpha[i]*(T_wall[i] - T_fluid[i]);

    Q_lat[i] = Q_latRateLim[i];

    Q_latRateLim[i] = if adiabatic then 0 else
      n_channels[i] * A_heat_internal[i] *
      (twoPhaseFraction[i]*htcoeff.alpha2ph[i]*(T_wall[i] - min(T_vap[i], max(T_liq[i], T_fluid[i])))
       +max(0,(1-onePhaseOutFraction[i]-twoPhaseFraction[i]))
         *htcoeff.alpha1phIn[i]*(T_wall[i] - min(T_vap[i], max(T_liq[i], T_fluid[i]))));

  end for;

  for i in 1:n-1 loop
    Q_cond_int[i]=(T[i] - T[i+1])/(
      L[i]/(2*A[i]*n_channels[i]*lam[i]) + L[i+1]/(2*A[i+1]*n_channels[i+1]*lam[i+1]));
  end for;

  //Up-stream properties for CVs partially in two-phase

  if positiveFlow then
    h_up = cat(1, {hA}, h[1:n-1]);
    //x_up = cat(1, {satA.x}, x[1:n-1]);
    else

    if n == 1 then
      if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
        frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        if upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont then
          h_up[1] = noEvent(if m_flow[2] >= 0 then hA else hB);
        elseif upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Linear then
          h_up[1] = .Modelon.Math.Smoothing.linearTransition(hA,
            .Modelon.Math.Smoothing.linearTransition(h[1], hB, m_flow[2] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[2] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
        else
          h_up[1] = .Modelon.Math.Smoothing.spliceFunction(hA,
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[1], hA, hB),hB, m_flow[2] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[2] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
        end if;
      else
        //if frictionDistribution == FrictionDistribution.FricVol then
        if upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont then
          h_up[1] = noEvent(if m_flow[1] >= 0 then hA else hB);
        elseif upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Linear then
          h_up[1] = .Modelon.Math.Smoothing.linearTransition(hA,
            .Modelon.Math.Smoothing.linearTransition(h[1], hB, m_flow[1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
        else
          h_up[1] = .Modelon.Math.Smoothing.spliceFunction(hA,
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[1], hA, hB), hB, m_flow[1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
        end if;
      end if;
    else
      if frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolFric or
        frictionDistribution == .Modelon.ThermoFluid.Choices.FrictionDistribution.VolVol then
        if upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont then
          h_up[1] = noEvent(if m_flow[2] >= 0 then hA else h[2]);
          h_up[n] = noEvent(if m_flow[n+1] < 0 then hB else h[n-1]);
          for i in 2:(n-1) loop
            h_up[i] = noEvent(if m_flow[i+1] >= 0 then h[i-1] else h[i+1]);
          end for;
        elseif upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Linear then
          h_up[1] = .Modelon.Math.Smoothing.linearTransition(hA,
            .Modelon.Math.Smoothing.linearTransition(h[1], h[2], m_flow[2] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[2] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          h_up[n_fric] = .Modelon.Math.Smoothing.linearTransition(h[n-1],
            .Modelon.Math.Smoothing.linearTransition(h[n], hB, m_flow[n+1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[n+1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          for i in 2:(n_fric-1) loop
            h_up[i] = .Modelon.Math.Smoothing.linearTransition(h[i-1],
              .Modelon.Math.Smoothing.linearTransition(h[i], h[i+1], m_flow[i+1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
              m_flow[i+1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          end for;
        else
          h_up[1] = .Modelon.Math.Smoothing.spliceFunction(hA,
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[1], hA, h[2]), h[2], m_flow[2] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[2] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          h_up[n_fric] = .Modelon.Math.Smoothing.spliceFunction(h[n-1],
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[n], h[n-1], hB), hB, m_flow[n+1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[n+1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          for i in 2:(n_fric-1) loop
            h_up[i] = .Modelon.Math.Smoothing.spliceFunction(h[i-1],
              .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[i], h[i-1], h[i+1]), h[i+1], m_flow[i+1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
              m_flow[i+1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          end for;
        end if;
      else
        //if frictionDistribution == FrictionDistribution.FricVol then
        if upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Discont then
          h_up[1] = noEvent(if m_flow[1] >= 0 then hA else h[2]);
          h_up[n] = noEvent(if m_flow[n] < 0 then hB else h[n-1]);
          for i in 2:(n-1) loop
            h_up[i] = noEvent(if m_flow[i] >= 0 then h[i-1] else h[i+1]);
          end for;
        elseif upstreamhAtReversal == .Modelon.ThermoFluid.Choices.UpstreamPropsAtReversal.Linear then
          h_up[1] = .Modelon.Math.Smoothing.linearTransition(hA,
            .Modelon.Math.Smoothing.linearTransition(h[1], h[2], m_flow[1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          h_up[n_fric] = .Modelon.Math.Smoothing.linearTransition(h[n-1],
            .Modelon.Math.Smoothing.linearTransition(h[n], hB, m_flow[n] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[n] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          for i in 2:(n_fric-1) loop
            h_up[i] = .Modelon.Math.Smoothing.linearTransition(h[i-1],
              .Modelon.Math.Smoothing.linearTransition(h[i], h[i+1], m_flow[i] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
              m_flow[i] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          end for;
        else
          h_up[1] = .Modelon.Math.Smoothing.spliceFunction(hA,
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[1], hA, h[2]), h[2], m_flow[1] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[1] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          h_up[n_fric] = .Modelon.Math.Smoothing.spliceFunction(h[n-1],
            .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[n], h[n-1], hB), hB, m_flow[n] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
            m_flow[n] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          for i in 2:(n_fric-1) loop
            h_up[i] = .Modelon.Math.Smoothing.spliceFunction(h[i-1],
              .Modelon.Math.Smoothing.spliceFunction(.Modelon.Math.RealScalar.between(h[i], h[i-1], h[i+1]), h[i+1], m_flow[i] + 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero)),
              m_flow[i] - 0.5*(mflow_smooth+mflow_zero), 0.5*(mflow_smooth-mflow_zero));
          end for;
        end if;
      end if;
    end if;

  end if;

  for i in 1:n loop
    T_up_one[i] = Medium.temperature_phX(p[i], h_up[i], volume[i].X);
    d_up_one[i] = noEvent(if h_up[i] > (h_liq[i]+h_vap[i])*0.5 then sat[i].d_vap else sat[i].d_liq);
    eta_up_one[i] = noEvent(if h_up[i] > (h_liq[i]+h_vap[i])*0.5 then sat[i].eta_vap else sat[i].eta_liq);
    cp_up_one[i] = noEvent(if h_up[i] > (h_liq[i]+h_vap[i])*0.5 then sat[i].cp_vap else sat[i].cp_liq);
    lam_up_one[i] = noEvent(if h_up[i] > (h_liq[i]+h_vap[i])*0.5 then sat[i].lam_vap else sat[i].lam_liq);
    Re_up_one[i] = .Modelon.Math.Smoothing.above(
      .Modelon.ThermoFluid.Functions.CharacteristicNumbers.ReynoldsNumber(
      abs(0.5*(m_flow[i]+m_flow[i+1]))/n_channels[i], Dhyd[i], A[i], eta_up_one[min(i,n_fric)]), 1, 10);
  end for;

  //Phase fraction, flow rate dependency calculation
  for i in 1:n loop

  onePhaseInFraction[i] = 
    .Modelon.Math.Smoothing.spliceFunction(
      0,
      if noEvent(
        (h_up[i] <= (if m_flow_help[i] >= 0 then (if i == n then h_liqB else h_liq[i+1]) else
           (if i == 1 then h_liqA else h_liq[i-1])) and
          h[i] <= (h_vap[i] + h_liq[i])*0.5)
        or (h_up[i] >= (if m_flow_help[i] >= 0 then (if i == n then h_vapB else h_vap[i+1]) else
           (if i == 1 then h_vapA else h_vap[i-1])) and
          h[i] >= (h_vap[i] + h_liq[i])*0.5)
        or (h[i] <=h_vap[i] and h[i] >= h_liq[i])) then
            max(0, min(L1in[i]/L[i], twoPhaseFractionZero[i]))
        else max(0, min(L1in[i]/L[i], 1)),
      p_red[i]-0.9999,
      1e-4);

  twoPhaseFraction[i] = if noEvent(
    (h_up[i] <= (if m_flow_help[i] >= 0 then (if i == n then h_liqB else h_liq[i+1]) else
       (if i == 1 then h_liqA else h_liq[i-1])) and
      h[i] <= (h_vap[i] + h_liq[i])*0.5)
    or (h_up[i] >= (if m_flow_help[i] >= 0 then (if i == n then h_vapB else h_vap[i+1]) else
      (if i == 1 then h_vapA else h_vap[i-1])) and
      h[i] >= (h_vap[i] + h_liq[i])*0.5)
    or (h[i] <=h_vap[i] and h[i] >= h_liq[i])) then
        max(0, twoPhaseFractionZero[i] - onePhaseInFraction[i])
    else max(0, min(1, twoPhaseFractionZero[i] + 
        .Modelon.Math.Smoothing.spliceFunction(
          0,
          max(0, L2in[i])/L[i],
        p_red[i]-0.9999,
        1e-4)) 
        - onePhaseInFraction[i]);

    onePhaseOutFraction[i] = 1 - twoPhaseFraction[i] - onePhaseInFraction[i];
  end for;


  /* Sizing */
  // PbS equations: Prescribe L_set, the total length if settings_TF.usePbS
  if not settings_TF.usePbS or length_fixed then
    // Use intermediate parameter variable L_par defined in initial equation permanently (L_set otherwise unused)
    L_set = L_par;
  else
    // Solve for L_set from prescribed pressure drop (in case of the latter the IV is held)
    par*L_set/L_1 = if sizing then (sum(dp_internal) - dp0)/dp_1 else (L_set - L_total)/L_1
      annotation (__Modelon(ResidualEquation(iterationVariable(enabled=settings_TF.usePbS, hold=not sizing)=L_set, enabled=settings_TF.usePbS, hold=not sizing)));
  end if;

  // PbS equations: Prescribe A_heat_set, the heat transfer area if settings_TF.usePbS
  if not settings_TF.usePbS or A_heat_fixed then
    // Use intermediate parameter variable A_heat_par defined in initial equation permanently (A_heat_set otherwise unused)
    A_heat_set = A_heat_par;
  else
    // Solve for A_heat_set from prescribed heat flow rate (in case of the latter the IV is held)
    if sizing then
      par*A_heat_set/A_heat_1 = if thermalOpt == .Modelon.ThermoFluid.Choices.ChannelThermalOpt.Q_flow then (Q_tot - Q_flow0)/Q_flow_1 else (T[end]-T_out0)/T_out_1
        annotation (__Modelon(ResidualEquation(iterationVariable(enabled=settings_TF.usePbS)=A_heat_set, enabled=settings_TF.usePbS)));
    else
      par*A_heat_set/A_heat_1 = (A_heat_set - sum(A_heat))/A_heat_1;
    end if;
  end if;

  annotation (defaultComponentName="pipe",
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}),
                                      graphics={
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
<p>Model of a pipe with <span style=\"font-family: Courier New;\">n</span> discrete control volumes and distributed flow resistance. Each control volume introduces an independent fluid state. The numerical state variables that are actually used depends on the selected medium model.</p>
<p>In this component, geometrical parameters are provided per control volume to allow for a wide range of configurations.</p>
<h4>Connection to other components</h4>
<p>The first control volume is connected to <span style=\"font-family: Courier New;\">portA</span>, a flow model is connected to <span style=\"font-family: Courier New;\">portB</span>. To assure a numerically sound model of alternating volume and flow models, connect a flow model to <span style=\"font-family: Courier New;\">portA</span> and a volume model to <span style=\"font-family: Courier New;\">portB</span>.</p>
<p>Heat transfer is computed between the fluid and heat transfer connector <span style=\"font-family: Courier New;\">q</span>. A wall temperature is assumed to be connected to <span style=\"font-family: Courier New;\">q</span>. The connector <span style=\"font-family: Courier New;\">q_fluid</span> exposes the internal fluid temperature and it is possible to add heat from an additional source through this port.</p>
<h4>Model assumptions</h4>
<ul>
<li>The pipe is discretized in n segments from inlet to outlet.</li>
<li>The mass flow rate into the first control volume is given by the connector <span style=\"font-family: Courier New;\">portA.</span> The mass flow rate between each control volume and from the last control volume to <span style=\"font-family: Courier New;\">portB </span>is calculated in a replaceable model <span style=\"font-family: Courier New;\">Friction, </span>from the pressure drop between the two components.</li>
<li>Heat flow rate into each control volume from the heat connector is calculated. The heat transfer coefficient is calculated in the replaceable model <span style=\"font-family: Courier New;\">HeatTransfer</span>.</li>
</ul>
<h4>Parameterization</h4>
<p><i>General</i> tab:</p>
<p><u>Medium</u></p>
<ul>
<li>Medium: Defines the package defining medium properties.</li>
</ul>
<p><u>Geometry</u></p>
<p>Each parameter in this group is provided per channel segment. The number of elements provided for the parameter L determines the number of segments n. All other paraemters must be provided as vectors with this length.</p>
<ul>
<li>n_channels: Allows the component to model multiple parallel channels under identical operating conditions. The model operates by computing all properties for a single channel and then multiplying the mass flow rate and transferred heat with the number of channels. Another way to model multiple parallel channels is to set n_channels = 1 and provide the total cross sectional and heat transfer areas of all channels.</li>
<li>L: Length of each channel segment.</li>
<li>Dhyd: Hydraulic diameter of the channel.</li>
<li>A: Cross-sectional area.</li>
<li>V: Channel volume (total in case of multiple channels). Automatically computed for most channel geometries from area and length. </li>
</ul>
<p><u>Heat transfer</u></p>
<ul>
<li>A_heat: Area transferring heat (orthogonal to heat flow direction) for each segment.</li>
<li>HeatTransfer: Replaceable model to characterize the heat transfer coefficient. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.Interfaces.TwoPhaseHeatTransfer\">this interface class</a>.</li>
<li>useMeanTempDrivenQ: If true, this would configure the heat flow to be driven by the mean temperature in each control volume. This would reduce the number of segments needed in the case where temperature difference between the fluids is small compared to the temperature gradients.</li>
</ul>
<p><u>Flow resistance</u></p>
<ul>
<li>Friction: Replaceable model to characterize the flow resistance, i.e. a correlation between mass flow rate and pressure drop. Pre-defined models are available in a drop-down menu. Custom models must extend <a href=\"modelica://Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.Interfaces.TwoPhaseResistance\">this interface class</a>.</li>
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
<li>initFromEnthalpy: Set to true to initialize by specifying initial specific enthalpies. If false, initial temperatures have to be specified instead.</li>
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
<p><br><br><i>Balance Equations tab:</i></p>
<p><u>Static head</u></p>
<ul>
<li>includeStaticHead: Set to true to account for pressure difference between the inlet and oulet port due to different elevations.</li>
<li>levels: Relative levels of the segments, 0 at inlet port. Only needed if static head is included.</li>
<li>g: Acceleration of gravity. The standard gravitational acceleration on earth is the default value. Only needed if static head is included.</li>
</ul>
<p><u>includeAcceleration: </u></p>
<p>This option is used to enable the acceleration due to specific volume which causes the pressure drop. </p>
<p>includeAcceleration - False: </p>
<p>This option considers only the frictional pressure drop.</p>
<p><img src=\"modelica://Modelon/Resources/Images/equations/useAcceldp_1.png\"/></p>
<p>includeAcceleration -True: </p>
<p>This option considers both frictional pressure drop and pressure drop due to the specific volume change</p>
<p><img src=\"modelica://Modelon/Resources/Images/equations/useAcceldp_2.png\"/> </p>
<p>Where, </p>
<p>G Mass flux. <img src=\"modelica://Modelon/Resources/Images/equations/massflux.png\"/> </p>
<p>v specific volume </p>
<p><u>kineticEnergyInBalance: </u></p>
<p>This option is used to enable the varying the energy due to change in velocity. </p>
<p>kineticEnergyInBalance - False: </p>
<p><img src=\"modelica://Modelon/Resources/Images/equations/useMechPower_1.png\"/> </p>
<p>The energy balance only considered the balancing of enthalpy change and additional heat flow transfer by wall </p>
<p>kineticEnergyInBalance True: </p>
<p><img src=\"modelica://Modelon/Resources/Images/equations/useMechPower_2.png\"/> </p>
<p>Where, </p>
<p>changeInKineticEnergy - kineticEnergy[1:n] - kineticEnergy[2:n+1] </p>
<p>kineticEnergy - m_flow * velocity^2 / 2.</p>
<p>The energy balance is also considering the change in energy due to the velocity change. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end DistributedTwoPhase;
