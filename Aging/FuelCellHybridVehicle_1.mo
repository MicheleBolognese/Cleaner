within Aging;
model FuelCellHybridVehicle_1 "Fuel cell for a hybrid vehicle"
  extends .Modelon.Icons.Experiment;

  package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformateShort
    "Medium at anode side";
  package Medium_cath =
      .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir "Medium at cathode side";
  package CondMedium_cath =
      .FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir
    "Condensing medium at cathode side";
  package MediumAir = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir
    "Air medium";
  package MediumWater = .Modelon.Media.PreDefined.TwoPhase.WaterIF97
    "Water medium";

  constant Integer[1] cath=Medium_cath.substanceIndexVector({"O2"})
    "index of required substances on the cathode side" annotation (Evaluate=true);
  .Modelica.Units.SI.MassFraction[Medium_cath.nS] X_air={0,0,X_H2O_ratio[1],X_H2O_ratio[2]*0.767,X_H2O_ratio[2]*
      0.233} "Inlet air composition: Ar, CO2, H2O, N2, O2";
  final parameter .Modelica.Units.SI.MassFraction X_H2O_ratio[2](each fixed=false, start={0.001,0.999})
    "Mass fraction of water and air";

  // Ambient conditions
  parameter .Modelica.Units.SI.AbsolutePressure p_ambient=101300 "Ambient pressure"
    annotation (Dialog(group="Ambient conditions"));
  parameter .Modelica.Units.SI.Temperature T_ambient=298.15 "Ambient temperature"
    annotation (Dialog(group="Ambient conditions"));
  parameter Real RH_ambient=0.40 "Relative Humidity of inlet air" annotation(Dialog(group="Ambient conditions"));

  // Fuel (anode) boundary and initial conditions
  parameter .Modelica.Units.SI.AbsolutePressure p_in_fuel=150000 "Inlet fuel pressure"
    annotation (Dialog(group="Fuel conditions"));
  parameter .Modelica.Units.SI.Temperature Tin_fuel=T_ambient "Inlet fuel temperature"
    annotation (Dialog(group="Fuel conditions"));
  parameter .Modelica.Units.SI.MoleFraction[Medium_an.nS] X_fuel=Medium_an.moleToMassFractions({1,0,0,0,0,0},
      Medium_an.MMX) "Inlet fuel composition  - NasaReformate: H2, CO, CO2,  H2O, N2, O2"
    annotation (Dialog(group="Fuel conditions"));
  final parameter .Modelica.Units.SI.MoleFraction[Medium_an.nS] X_start_an_condencing={X_fuel[4],X_fuel[2],X_fuel[
      3],X_fuel[1],X_fuel[5],X_fuel[6]} "Stack anode side starts filled with same medium as anode side";

  // Initialization
  parameter .Modelica.Units.SI.Temperature T_stack_start(displayUnit="K")=300
                                                               "Stack initial temperature"
    annotation (Dialog(tab="Initialization"));
  parameter .Modelica.Units.SI.Temperature T_humidifier_start=max(330.15, T_stack_start - 5)
    "Start value for humidifier temperature, must not be below 57 C" annotation (Dialog(tab="Initialization"));
  parameter .Modelica.Units.SI.Pressure pstart_an=p_in_fuel "Start pressure (anode side)"
    annotation (Dialog(tab="Initialization"));
  parameter .Modelica.Units.SI.Pressure pstart_cath=250000 "Start pressure (cathode side)"
    annotation (Dialog(tab="Initialization"));
  parameter .Modelica.Units.SI.MassFlowRate m_fuel_start=1e-3 "Start value for fuel flow rate"
    annotation (Dialog(tab="Initialization"));
  parameter Real waterFlow_init = 0.0001
    "Initialization value for the water flow" annotation(Dialog(tab="Initialization"));

  // Controls
  parameter Real RH_stack_setpoint=0.9 "Desired RH of air in stack" annotation(Dialog(group="Controls"));
  parameter Real cath_stoich=1.9
    "Cathode stoichiometry, expected to be about 2" annotation(Dialog(group="Controls"));
  parameter Real comp_speed=1230*2*.Modelica.Constants.pi/60
    "Compressor idle speed and start value" annotation(Dialog(group="Controls"));

  // Medium index for fluid convertion
  parameter Integer[Medium_an.nS] index_an={4,2,3,1,5,6}
    "Index of anode medium for fluid convertion" annotation(Dialog(tab="Advanced"));
  parameter Integer[Medium_cath.nS] index_cath={3,2,1,4,5}
    "Index of cathode medium for fluid convertion"  annotation(Dialog(tab="Advanced"));

  .Modelica.Units.SI.Power P_stack_net=coolStack.subStack.P_stack - compressor.compressor.shaftPower -
      pumpCooling.W_tot - condenser_pump.pumpShaft.W_tot "Net power of stack";
  .Modelica.Units.SI.Efficiency eta_stack=max(0, P_stack_net/max(coolStack.subStack.P_stack, 1)*coolStack.subStack.cell.eff_cell)
    "Efficiency of the stack";
  Real RH_stack_measured=CondMedium_cath.relativeHumidity_pTX(
      coolStack.drain_cath.p,
      coolStack.subStack.T_stack[1],
      coolStack.drain_cath.X_outflow) "Relative Humidity in stack";
  .Modelica.Units.SI.Temperature T_pressureHandler=heatExchanger.summary.T_out
    "Temperature pressure handler component";

  .FuelCell.Sources.GasPressureBoundary sourceHydrogen(
    fluidPort(m_flow(start=m_fuel_start)),
    X=X_fuel,
    T=Tin_fuel,
    use_Th_in=false,
    redeclare package Medium = Medium_an,
    pressureUnit=.Modelon.ThermoFluid.Choices.RealPressureUnit.bar,
    use_p_in=false,
    p=p_in_fuel) annotation (Placement(transformation(extent={{-94,64},{-78,80}},
          rotation=0)));
  .FuelCell.Stacks.PEMFC.CondStack coolStack(
    n_cell=370,
    use_wall=true,
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    pstart_cath=pstart_cath,
    redeclare package CondensingMedium_an = .FuelCell.Media.PreDefined.CondensingGases.CondensingReformateShort,
    redeclare package CondensingMedium_cath = .FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir,
    includeStaticHead_cooling=false,
    redeclare package Medium_cooling = MediumWater,
    coolingPipe(
      NA=1,
      NB=1,
      redeclare model Friction = .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction (
          h0_in=1.5e5,
          h0_out=2e5,
          mflow0=0.4)),
    Xstart_an=X_fuel,
    subStack(
      redeclare replaceable model Membrane = .FuelCell.Membranes.PEMFC.Empirical,
      addProxToAnode=false,
      rhonom_cath=2.5,
      m_flow_nom_cath=1e-1,
      anode_channel(Xstart=X_start_an_condencing, water_flow(start=waterFlow_init)),
      cathode_channel(water_flow(start=waterFlow_init)),
      m_flow_nom_an=1e-3,
      rhonom_an=0.1,
      dp_nom_an=2000,
      dp_nom_cath=150000),
    pstart_an=pstart_an - 400,
    Tstart_an=T_stack_start,
    Tstart_cath=T_stack_start) annotation (Placement(transformation(extent={{12,-22},{66,22}})));

  .FuelCell.Sources.CondensingGasPressureBoundary sinkCathode(redeclare package
      Medium = CondMedium_cath)           annotation (Placement(transformation(
        extent={{8,8},{-8,-8}},
        rotation=180,
        origin={-60,-60})));
  .Modelon.Utilities.Deprecated.Compressors.DynamicCompressor compressor(
    redeclare package Medium = Medium_cath,
    redeclare model CharMap =
        .Modelon.Utilities.Deprecated.Compressors.Characteristics.Dynamic.TableBasedFlowFromPR (
        flowMap=[0,1,1.5,2,2.5; 1500,9e-2,5e-2,3e-2,1e-2; 2000,13e-2,9e-2,6e-2,3e-2;
            2500,17e-2,13e-2,11e-2,9e-2; 3000,21e-2,17e-2,15e-2,13e-2],
        normalizedCorrection=true,
        effMap=[0,1.5,2,2.5,3; 500,0.45,0.4,0.4,0.4; 1000,0.65,0.6,0.55,0.55; 1500,
            0.7,0.75,0.7,0.7; 2000,0.7,0.8,0.8,0.85; 2500,0.7,0.85,0.85,0.8]),
    eta_mech=0.75,
    T_start=T_ambient,
    positiveFlow=false)
    annotation (Placement(transformation(extent={{-70,-20},{-50,0}})));
  .Modelica.Mechanics.Rotational.Sources.Speed motorSpeed annotation (Placement(
        transformation(
        extent={{-6.0,-6.0},{6.0,6.0}},
        rotation=-90.0,
        origin={-58.72135631647848,20.319660920880377})));
  .FuelCell.Sources.CondensingGasPressureBoundary sourceCathode(
    T=T_ambient,
    p=p_ambient,
    use_X_in=true,
    redeclare package Medium = Medium_cath) annotation (Placement(
        transformation(
        extent={{8,8},{-8,-8}},
        rotation=180,
        origin={-84,-10})));
  .Modelica.Blocks.Sources.RealExpression InAir[Medium_cath.nS](y=X_air)
    annotation (Placement(transformation(extent={{-100,2},{-84,18}})));

  .FuelCell.Interfaces.FluidConverter fluidConverter(
    redeclare package Medium_a =
        .FuelCell.Media.PreDefined.CondensingGases.CondensingReformateShort,
    index=index_an,
    redeclare package Medium_b = Medium_an,
    Tstart=T_stack_start)
                   annotation (Placement(transformation(
        extent={{-6,6},{6,-6}},
        rotation=180,
        origin={70,88})));
  .Modelon.ThermoFluid.Pumps.PumpShaft pumpCooling(
    redeclare package Medium = MediumWater,
    mflow_smooth=1e-5,
    m_flow_start=0.0015,
    redeclare function efficiencyCharacteristic =
        .Modelon.ThermoFluid.Pumps.PumpCharacteristics.constantEfficiency,
    effParam=.Modelon.ThermoFluid.Choices.EfficiencyParameterization.Power,
    redeclare function powerCharacteristic =
        .Modelon.ThermoFluid.Pumps.PumpCharacteristics.linearPower (q_nom={0.01e-3,
            6e-3}, W_nom={0.5,300}),
    T_start=T_ambient,
    redeclare replaceable function flowCharacteristic =
        .Modelon.ThermoFluid.Pumps.PumpCharacteristics.linearFlow (head_nom={6,0},
          q_nom={0,0.10}),
    N_nom=5000) annotation (Placement(transformation(
        extent={{-8,8},{8,-8}},
        rotation=90,
        origin={28,-52})));
  .FuelCell.Sources.CondensingGasPressureBoundary coolingAirOut(redeclare package
              Medium = MediumAir,
    p=p_ambient,
    T=T_ambient + 35)
    annotation (Placement(transformation(
        extent={{8,8},{-8,-8}},
        rotation=180,
        origin={-60,-86})));
  .FuelCell.Sources.GasFlowBoundary coolingAirIn(
    redeclare package Medium =
        MediumAir,
    T=T_ambient,
    m_flow=11)
              annotation (Placement(transformation(extent={{48,-94},{32,-78}})));
  .FuelCell.Examples.PEMFC.Components.PressureHandler pressureHandler(redeclare package
              Medium = MediumWater, T=T_ambient)
    annotation (Placement(transformation(extent={{48,-80},{28,-60}})));
  .FuelCell.HeatExchangers.GasGas.EpsNTU heatExchanger(
    dp_prim_start=10000,
    T_prim_start=298.15,
    dp_sec_start=10000,
    T_sec_start=298.15,
    redeclare function epsFun =
        .Modelon.ThermoFluid.HeatExchangers.Functions.counterFlowEps,
    redeclare model Friction_sec =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.DarcyWeisbach (
         zeta_A=20),
    redeclare model Friction_prim =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.DarcyWeisbach (
         zeta_A=10),
    redeclare model HeatTransfer_prim =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    redeclare model HeatTransfer_sec =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    heatFlowDynamics=false,
    redeclare package PrimaryMedium = MediumWater,
    redeclare package SecondaryMedium = MediumAir,
    Dhyd_prim=2e-3,
    A_prim=7e-5,
    Dhyd_sec=2e-3,
    A_sec=7e-5,
    n_channels_prim=600,
    n_channels_sec=300,
    L_prim=0.5,
    A_heat_prim=0.7,
    L_sec=0.6,
    A_heat_sec=0.084,
    A_heat=0.084*300*100 + 200*0.7)
    annotation (Placement(transformation(extent={{0,-90},{-20,-70}})));

  .FuelCell.Interfaces.FluidConverter fluidConverter1(
    redeclare package Medium_a = Medium_cath,
    redeclare package Medium_b = CondMedium_cath,
    index=index_cath,
    Tstart=T_ambient)
                   annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=0,
        origin={-40,-10})));
  .FuelCell.Examples.PEMFC.Components.Controls.CoolingSystemControl coolingSystemControl(T_stack_ref=
        343.15, k=-2000)
                       annotation (Placement(transformation(rotation=0, extent={{62,-60},
            {78,-44}})));

  .Modelon.ThermoFluid.Volumes.TwoportVolume volume(redeclare package Medium =
        Medium_an,
    V=0.016,
    p_start(displayUnit="Pa") = 130000,
    T_start=T_stack_start,
                   X_start=X_fuel)
    annotation (Placement(transformation(extent={{34,82},{46,94}})));
  .Modelica.Blocks.Sources.RealExpression RH_stack(y=RH_stack_measured) annotation (
      Placement(transformation(
        extent={{8,-8},{-8,8}},
        rotation=180,
        origin={-84,-44})));
  .Modelon.Blocks.Control.LimPID PID(
    yMax=1500,
    yMin=0,
    Ti=9,
    k=300,
    Td=0.5,
    controllerType=.Modelica.Blocks.Types.SimpleController.PI,
    initType=.Modelica.Blocks.Types.Init.InitialOutput,
    limitsAtInit=false,
    y_start=(PID.yMax + PID.yMin)/2) annotation (Placement(transformation(
        extent={{6,6},{-6,-6}},
        rotation=180,
        origin={-60,-30})));
  .Modelica.Blocks.Sources.RealExpression RH_ref(y=RH_stack_setpoint)
    annotation (Placement(transformation(extent={{-92,-38},{-76,-22}})));
  .Modelon.ThermoFluid.Sensors.MassflowSensor massflowSensor1(redeclare package
      Medium = Medium_an)
    annotation (Placement(transformation(extent={{-9.644134498036415,18.45079286829087},{2.3558655019635846,30.45079286829087}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelon.ThermoFluid.Sensors.MassflowSensor massflowSensor2(redeclare package
      Medium =
      .FuelCell.Media.PreDefined.CondensingGases.CondensingReformateShort)
    annotation (Placement(transformation(extent={{74,6},{86,18}})));
  .FuelCell.Examples.PEMFC.Components.Controls.CompressorControl compressorControl(compSpeed=comp_speed)
                       annotation (Placement(transformation(rotation=0.0, extent={{-41.80470291740819,43.518389374181815},{-25.804702917408193,59.518389374181815}},origin = {0.0,0.0})));
  .FuelCell.Examples.PEMFC.Components.PumpSystemCondensing condenser_pump(
    redeclare package Medium_cath = CondMedium_cath,
    redeclare package WaterMedium = MediumWater,
    p_system=pstart_cath,
    T_separator_start=T_humidifier_start,
    X_cath={0.115,0,0,0.75,0.135},
    separator(water_flow(start=waterFlow_init)),
    pumpShaft(m_flow_start=0.1))
    annotation (Placement(transformation(extent={{-20,-22},{0,-42}})));

  .FuelCell.Volumes.IdealHumidifier humidifier(
    redeclare package Medium_Gas = CondMedium_cath,
    redeclare package WaterMedium = MediumWater,
    p_start=pstart_cath,
    T_start=T_humidifier_start,
    X_start={0.115,0,0,0.75,0.135},
    heatOn=true,
    simpleMixer(pstart=230000),
    preheat_pipe(A={7.85e-5}))
    annotation (Placement(transformation(extent={{-20,-18},{0,2}})));
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    duration=100,
    offset=0,
    startTime=1,
    I=500)
    annotation (Placement(transformation(
        origin={39,33},
        extent={{-11,-11},{11,11}},
        rotation=180)));
  .Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{10,24},{22,36}},     rotation=
           0)));
  .Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
        T_ambient) annotation (Placement(transformation(
        extent={{-4,-4},{4,4}},
        rotation=90,
        origin={-10,-50})));
  .Modelon.ThermoFluid.Valves.ValveCompressible valve(
    redeclare package Medium = Medium_an,
    T_start=Tin_fuel,
    mflow_start=m_fuel_start,
    mflow_smooth=1e-7,
    flowCoeff=.Modelon.ThermoFluid.Choices.CvTypes.Kv,
    Kv=15,
    Fxt_full=0.8)
    annotation (Placement(transformation(extent={{-72,64},{-58,80}})));
  .Modelon.ThermoFluid.Volumes.TwoportVolume volume1(
    redeclare package Medium = Medium_an,
    V(displayUnit="m3") = 0.016,
    initFromEnthalpy=false,
    p_start=pstart_an,
    h_start=5e5,
    X_start=X_fuel)
    annotation (Placement(transformation(extent={{-51.72183265916388,66.0},{-39.72183265916388,78.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelon.Blocks.Control.LimPID PID_valve(
    yMax=1,
    yMin=0,
    Ti=0.3,
    k=2,
    Td=0.5,
    controllerType=.Modelica.Blocks.Types.SimpleController.PI,
    initType=.Modelica.Blocks.Types.Init.InitialOutput,
    limitsAtInit=false,
    y_start=0.1)
               annotation (Placement(transformation(
        extent={{-6.0,-6.0},{6.0,6.0}},
        rotation=180.0,
        origin={-53.07192272804079,123.18146453884918})));
  .Modelica.Blocks.Sources.RealExpression m_flow_suc_sp(y=0.00008)
    annotation (Placement(transformation(extent={{9.662654997465413,115.18146453884916},{-6.337345002534587,131.18146453884916}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelica.Blocks.Sources.RealExpression m_flow_suc(y=massflowSensor2.y)
    annotation (Placement(transformation(extent={{10.541078179053159,126.0102336300655},{-5.458921820946841,142.0102336300655}},rotation = 0.0,origin = {0.0,0.0})));
    .Modelon.Utilities.Deprecated.Compressors.DynamicCompressor compressor2(redeclare replaceable package Medium = Medium_an,redeclare replaceable model CharMap = .Modelon.Utilities.Deprecated.Compressors.Characteristics.Dynamic.TableBasedSAE) annotation(Placement(transformation(extent = {{-32.810508279613686,62.28338117474988},{-12.810508279613689,82.28338117474988}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 100) annotation(Placement(transformation(extent = {{-49.399316905157804,91.12906315789718},{-34.22231756057346,106.30606250248152}},origin = {0.0,0.0},rotation = 0.0)));
initial equation
  RH_ambient =
    .Modelon.Media.PreDefined.CondensingGases.MoistAir.relativeHumidity_pTX(
    p_ambient,
    T_ambient,
    X_H2O_ratio) "Air composition calculation ";
  X_H2O_ratio[1] + X_H2O_ratio[2] = 1 "Air composition constraint ";

equation
  connect(sourceCathode.fluidPort, compressor.portA)
    annotation (Line(points={{-76.8,-10},{-70,-10}},
                                               color={209,60,0}));
  connect(coolingAirIn.fluidPort, heatExchanger.feed_sec)
    annotation (Line(points={{32.8,-86},{-1,-86}}, color={255,128,0}));
  connect(heatExchanger.drain_sec, coolingAirOut.fluidPort)
    annotation (Line(points={{-19,-86},{-52.8,-86}}, color={255,128,0}));
  connect(InAir.y, sourceCathode.X_in)
    annotation (Line(points={{-83.2,10},{-79.2,10},{-79.2,-2}},
                                                            color={0,0,127}));
  connect(volume.portB, fluidConverter.drain)
    annotation (Line(points={{46,88},{65.8,88}}, color={85,170,255}));
  connect(coolStack.feed_an, massflowSensor1.portB)
    annotation (Line(points={{14.7,11},{14.7,24.45079286829087},{2.3558655019635846,24.45079286829087}},
                                            color={255,128,0}));
  connect(coolStack.drain_an, massflowSensor2.portA)
    annotation (Line(points={{63.3,11},{63.3,12},{74,12}},
                                             color={209,60,0}));
  connect(massflowSensor2.portB, fluidConverter.feed) annotation (Line(points={{86,12},
          {90,12},{90,88},{74.2,88}},      color={85,170,255}));
  connect(fluidConverter1.feed, compressor.portB)
    annotation (Line(points={{-44.2,-10},{-50,-10}},
                                                 color={255,128,0}));
  connect(compressor.flange, motorSpeed.flange)
    annotation (Line(points={{-60,-3},{-60,14.31966092088038},{-58.72135631647848,14.31966092088038}},color={0,0,0}));
  connect(RH_ref.y, PID.u_s)
    annotation (Line(points={{-75.2,-30},{-67.2,-30}}, color={0,0,127}));
  connect(RH_stack.y, PID.u_m) annotation (Line(points={{-75.2,-44},{-60,-44},{-60,
          -37.2}}, color={0,0,127}));

  connect(current.n, coolStack.pin_n)
    annotation (Line(points={{28,33},{28,18.15},{28.2,18.15}},
                                                    color={0,0,255}));
  connect(current.p, coolStack.pin_p)
    annotation (Line(points={{50,33},{50,18.15},{49.8,18.15}},
                                                    color={0,0,255}));
  connect(fixedTemperature.port, condenser_pump.cooler[1])
    annotation (Line(points={{-10,-46},{-10,-42}}, color={191,0,0}));
  connect(pumpCooling.shaft, coolingSystemControl.flange) annotation (Line(
        points={{35.2,-52},{62,-52}},                         color={0,0,0}));
  connect(coolingSystemControl.port, coolStack.wall1[1]) annotation (Line(
        points={{78.16,-52},{90,-52},{90,0},{66,0}}, color={191,0,0}));
  connect(ground.p, coolStack.pin_n) annotation (Line(points={{16,36},{28.2,36},
          {28.2,18.15}}, color={0,0,255}));
  connect(massflowSensor2.y, compressorControl.u1)
    annotation (Line(points={{80,17.4},{80,51.518389374181815},{-25.804702917408193,51.518389374181815}}, color={0,0,127}));
  connect(massflowSensor1.y, compressorControl.u2)
    annotation (Line(points={{-3.6441344980364154,29.850792868290867},{-3.6441344980364154,36.8353839895272},{-33.80470291740819,36.8353839895272},{-33.80470291740819,44.31838937418181}}, color={0,0,127}));
  connect(compressorControl.y, motorSpeed.w_ref)
    annotation (Line(points={{-42.60470291740819,51.518389374181815},{-58.72135631647848,51.518389374181815},{-58.72135631647848,27.519660920880373}}, color={0,0,127}));
  connect(humidifier.coolingFeed, coolStack.drain_cooling) annotation (Line(
        points={{0,-2},{8,-2},{8,-30},{49.8,-30},{49.8,-18.15}}, color={0,190,0}));
  connect(humidifier.humidifierDrain, coolStack.feed_cath) annotation (Line(
        points={{0,-10},{14,-10},{14,-11},{14.7,-11}}, color={209,60,0}));
  connect(humidifier.humidifierFeed, fluidConverter1.drain)
    annotation (Line(points={{-20,-10},{-35.8,-10}}, color={209,60,0}));
  connect(humidifier.simpleMixerPressure, condenser_pump.p_tank) annotation (
      Line(points={{0,-14},{4,-14},{4,-36},{-1,-36}}, color={0,0,127}));
  connect(PID.y, condenser_pump.pump_speed)
    annotation (Line(points={{-53.4,-30},{-19,-30}}, color={0,0,127}));
  connect(condenser_pump.separatorFeed, coolStack.drain_cath) annotation (Line(
        points={{0,-39},{80,-39},{80,-11},{63.3,-11}}, color={209,60,0}));
  connect(condenser_pump.separatorDrain, sinkCathode.fluidPort) annotation (
      Line(points={{-20,-39},{-40,-39},{-40,-60},{-52.8,-60}}, color={209,60,0}));
  connect(humidifier.coolingDrain, heatExchanger.feed_prim) annotation (Line(
        points={{-20,-2},{-30,-2},{-30,-74},{-19,-74}}, color={0,190,0}));
  connect(humidifier.flowPort, condenser_pump.portB)
    annotation (Line(points={{-10,-17.8},{-10,-22}}, color={255,128,0}));
  connect(pressureHandler.port, heatExchanger.drain_prim)
    annotation (Line(points={{28,-72},{28,-74},{-1,-74}}, color={0,190,0}));
  connect(pressureHandler.port1, pumpCooling.portA)
    annotation (Line(points={{28,-67},{28,-60}}, color={0,190,0}));
  connect(pumpCooling.portB, coolStack.feed_cooling) annotation (Line(points={{28,
          -44},{28,-18.15},{28.2,-18.15}}, color={0,190,0}));
  connect(volume1.portA,valve. portB) annotation (Line(points={{-51.72183265916388,72},{-58,72}},
                                   color={85,170,255}));
  connect(m_flow_suc_sp.y, PID_valve.u_s)
    annotation (Line(points={{-7.1373450025345875,123.18146453884916},{-45.871922728040786,123.18146453884916},{-45.871922728040786,123.18146453884918}}, color={0,0,127}));
  connect(sourceHydrogen.fluidPort, valve.portA)
    annotation (Line(points={{-78.8,72},{-72,72}}, color={85,170,255}));
  connect(PID_valve.u_m, m_flow_suc.y) annotation (Line(points={{-53.07192272804079,130.38146453884917},{-53.07192272804079,134.0102336300655},{-6.258921820946842,134.0102336300655}}, color={0,0,127}));
  connect(PID_valve.y, valve.opening) annotation (Line(
      points={{-59.67192272804079,123.18146453884918},{-65,123.18146453884918},{-65,79.2}},
      color={0,0,127},
      pattern=LinePattern.Dash));
    connect(volume1.portB,compressor2.portA) annotation(Line(points = {{-39.72183265916388,72},{-32.810508279613686,72},{-32.810508279613686,72.28338117474988}},color = {255,128,0}));
    connect(constantSpeed.flange,compressor2.flange) annotation(Line(points = {{-34.22231756057346,98.71756283018935},{-22.81050827961369,98.71756283018935},{-22.81050827961369,79.28338117474988}},color = {0,0,0}));
    connect(volume.portA,compressor2.portA) annotation(Line(points = {{34,88},{-32.810508279613686,88},{-32.810508279613686,72.28338117474988}},color = {255,128,0}));
    connect(compressor2.portB,massflowSensor1.portA) annotation(Line(points = {{-12.810508279613689,72.28338117474988},{-12.810508279613689,24.45079286829087},{-9.644134498036415,24.45079286829087}},color = {255,128,0}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,120}})),
    experiment(
      StopTime=500,
      Tolerance=1e-08),
    Documentation(info="<html>
<p><b>FuelCellHybridVehicle</b></p>
<p>This fuel cell system example consists of a proton exchange membrane fuel cell (PEMFC) stack, a humidifier, a heat exchanger 
and several fluid machines for air, hydrogen and coolant flow. In order to improve efficiency and circumvent the use of after-burner, 
the excessive hydrogen on the anode side is recirculated in the system through an ejector model. 
A proportional controller is used to control the air flow on the cathode side. 
The controller is based on hydrogen consumption on the anode side and stoichiometry. 
A PID controller controls the cooling water pump and uses stack temperature as input. 
As it is important to keep the cells hydrated, the model also contains a humidifier that humidifies incoming air to the fuel cell to a desired level. 
The humidification is performed by using the water produced in the cathode.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"),
    Icon(coordinateSystem(extent={{-100,-100},{100,120}})));
end FuelCellHybridVehicle_1;
