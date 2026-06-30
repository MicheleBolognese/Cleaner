within SAMU.Humidifier_funzionante;
// model Test_1

//   .SAMU.Humidifier_funzionante.DiscretizedGasGasHumidifier discretizedGasGasHumidifier(
//     from_dp_sec = false,
//     positiveFlow_prim = false,
//     positiveFlow_sec = false,
//     WallSteadyStateInit = true,
//     n = 5,
//     from_dp_prim = false,
//     dp_smooth = 10,
//     redeclare replaceable model Friction_prim = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss,
//     redeclare replaceable model Friction_sec = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss
//   ) annotation(Placement(transformation(extent = {{12,-2},{32,18}})));

//   .FuelCell.Sources.CondensingGasPressureBoundary wetSink(p = 131722)
//     annotation(Placement(transformation(extent = {{-82,26},{-62,46}})));
//   .FuelCell.Sources.CondensingGasPressureBoundary drySink(p = 131722)
//     annotation(Placement(transformation(extent = {{128,24.38},{108,44.38}})));

//   // Sensori 
//   .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor
//     annotation(Placement(transformation(extent = {{-38,-14},{-18,-34}})));
//   .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot(displayUnits = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{-19.7,-69.7},{15.7,-34.3}})));
//   .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay(displayMassUnit = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{-48,-43.5},{-28,-60.5}})));
//   .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot2(displayUnits = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{-23.33,46.67},{15.33,85.33}})));
//   .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay2(displayMassUnit = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{-50,56.5},{-30,73.5}})));
//   .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor2
//     annotation(Placement(transformation(extent = {{-16,26},{-36,46}})));
//   .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot3(displayUnits = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{39.31,47.31},{76.69,84.69}})));
//   .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay3(displayMassUnit = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{82,57.5},{102,74.5}})));
//   .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor3
//     annotation(Placement(transformation(extent = {{68,24.38},{88,44.38}})));
//   .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot4(displayUnits = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{41.72,-68.28},{74.28,-35.72}})));
//   .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay4(displayMassUnit = true,
//       redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
//     annotation(Placement(transformation(extent = {{88,-41.5},{108,-58.5}})));
//   .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor4
//     annotation(Placement(transformation(extent = {{90,-13.62},{70,-33.62}})));
//     .FuelCell.Sources.CondensingGasFlowBoundary gasFlowSource(use_flow_in = true,T = 323.15,X = {0.005,0.001,0.013,0.752,0.23}) annotation(Placement(transformation(extent = {{-82,-34},{-62,-14}},origin = {0,0},rotation = 0)));
//     .FuelCell.Sources.CondensingGasFlowBoundary gasFlowSource2(use_flow_in = true,T = 343.15,X = {0.162,0.001,0.011,0.633,0.194}) annotation(Placement(transformation(extent = {{132.0,-34.0},{112.0,-14.0}},origin = {0.0,0.0},rotation = 0.0)));
//     .Modelica.Blocks.Sources.CombiTimeTable dryTable (    
//       table = [0, 0.001;
//              250, 0.0014;
//              500, 0.0015;
//              750, 0.002;
//              1000, 0.0025;
//              1250, 0.003;
//              1500, 0.0035;
//              1750, 0.004;
//              2000, 0.004;
//              2250, 0.0025;
//              2500, 0.0015],
//      timeScale = 1,
//      startTime = 0
//      ) annotation(Placement(transformation(extent = {{-122.0,-48.0},{-102.0,-28.0}},origin = {0.0,0.0},rotation = 0.0)));
//    .Modelica.Blocks.Sources.CombiTimeTable wetTable(startTime = 0,timeScale = 1,table = [0, 0.001;
//              250, 0.0014;
//              500, 0.0015;
//              750, 0.002;
//              1000, 0.0025;
//              1250, 0.003;
//              1500, 0.0035;
//              1750, 0.004;
//              2000, 0.004;
//              2250, 0.0025;
//              2500, 0.0015]) annotation(Placement(transformation(extent = {{176.0,-52.0},{156.0,-32.0}},origin = {0.0,0.0},rotation = 0.0)));

// equation

//   // Tutte le altre connessioni rimangono identiche
//   connect(airCompositionDisplay.data,condensingGasMultiDisplaySensor.u)
//     annotation(Line(points={{-38,-43.5},{-38,-33.77},{-27.95,-33.77},{-27.95,-24.05}},color={0,0,0}));
//   connect(multiDisplay_phTmdot.y,condensingGasMultiDisplaySensor.u)
//     annotation(Line(points={{-2,-52},{-2,-34},{-27.95,-34},{-27.95,-24.05}},color={0,0,0}));
//   connect(airCompositionDisplay2.data,condensingGasMultiDisplaySensor2.u)
//     annotation(Line(points={{-40,56.5},{-40,46.275},{-26.05,46.275},{-26.05,36.05}},color={0,0,0}));
//   connect(multiDisplay_phTmdot2.y,condensingGasMultiDisplaySensor2.u)
//     annotation(Line(points={{-4,66},{-4,46},{-26.05,46},{-26.05,36.05}},color={0,0,0}));
//   connect(condensingGasMultiDisplaySensor2.portA,discretizedGasGasHumidifier.portB_prim)
//     annotation(Line(points={{-20,36},{-9,36},{-9,12},{12,12}},color={209,60,0}));
//   connect(condensingGasMultiDisplaySensor.portB,discretizedGasGasHumidifier.portA_sec)
//     annotation(Line(points={{-22,-24},{-10,-24},{-10,8},{12,8}},color={209,60,0}));
//   connect(discretizedGasGasHumidifier.portB_sec,condensingGasMultiDisplaySensor3.portA)
//     annotation(Line(points={{32,8},{50,8},{50,34.38},{72,34.38}},color={209,60,0}));
//   connect(airCompositionDisplay3.data,condensingGasMultiDisplaySensor3.u)
//     annotation(Line(points={{92,57.5},{92,46.17},{78.05,46.17},{78.05,34.42}},color={0,0,0}));
//   connect(condensingGasMultiDisplaySensor3.u,multiDisplay_phTmdot3.y)
//     annotation(Line(points={{78.05,34.42},{78.05,46.39},{58,46.39},{58,66}},color={0,0,0}));
//   connect(condensingGasMultiDisplaySensor4.portB,discretizedGasGasHumidifier.portA_prim)
//     annotation(Line(points={{74,-23.62},{51,-23.62},{51,4},{32,4}},color={209,60,0}));
//   connect(airCompositionDisplay4.data,condensingGasMultiDisplaySensor4.u)
//     annotation(Line(points={{98,-41.5},{98,-33.61},{79.95,-33.61},{79.95,-23.66}},color={0,0,0}));
//   connect(multiDisplay_phTmdot4.y,condensingGasMultiDisplaySensor4.u)
//     annotation(Line(points={{58,-52},{58,-33.61},{79.95,-33.61},{79.95,-23.66}},color={0,0,0}));
//   connect(drySink.fluidPort,condensingGasMultiDisplaySensor3.portB)
//     annotation(Line(points={{109,34.39},{84,34.39}},color={209,60,0}));
//   connect(wetSink.fluidPort,condensingGasMultiDisplaySensor2.portB)
//     annotation(Line(points={{-63,36},{-32,36}},color={209,60,0}));
//     connect(gasFlowSource.fluidPort,condensingGasMultiDisplaySensor.portA) annotation(Line(points = {{-63,-24},{-34,-24}},color = {209,60,0}));
//     connect(gasFlowSource2.fluidPort,condensingGasMultiDisplaySensor4.portA) annotation(Line(points = {{113,-24},{100,-24},{100,-23.62},{86,-23.62}},color = {209,60,0}));
//     connect(dryTable.y[1],gasFlowSource.m_flow_in) annotation(Line(points = {{-101,-38},{-95,-38},{-95,-9},{-76,-9},{-76,-15}},color = {0,0,127}));
//     connect(wetTable.y[1],gasFlowSource2.m_flow_in) annotation(Line(points = {{155,-42},{149,-42},{149,-9},{126,-9},{126,-15}},color = {0,0,127}));

//   annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
//                   graphics={Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100,-100},{100,100}}),
//                              Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
// end Test_1;

model Test_1

  .SAMU.Humidifier_funzionante.DiscretizedGasGasHumidifier discretizedGasGasHumidifier(
    from_dp_sec = false,
    positiveFlow_prim = false,
    positiveFlow_sec = false,
    WallSteadyStateInit = true,
    n = 5,
    from_dp_prim = false,
    dp_smooth = 10,
    redeclare replaceable model Friction_prim = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss,
    redeclare replaceable model Friction_sec = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss
  ) annotation(Placement(transformation(extent = {{12.0,-2.0},{32.0,18.0}},rotation = 0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.CondensingGasPressureBoundary wetSink(p = 131722)
    annotation(Placement(transformation(extent = {{-82,26},{-62,46}})));
  .FuelCell.Sources.CondensingGasPressureBoundary drySink(p = 131722)
    annotation(Placement(transformation(extent = {{128,24.38},{108,44.38}})));

  // Sensori 
  .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor
    annotation(Placement(transformation(extent = {{-38,-14},{-18,-34}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot(displayUnits = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{-19.7,-69.7},{15.7,-34.3}})));
  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay(displayMassUnit = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{-48,-43.5},{-28,-60.5}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot2(displayUnits = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{-23.33,46.67},{15.33,85.33}})));
  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay2(displayMassUnit = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{-50,56.5},{-30,73.5}})));
  .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor2
    annotation(Placement(transformation(extent = {{-16,26},{-36,46}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot3(displayUnits = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{39.31,47.31},{76.69,84.69}})));
  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay3(displayMassUnit = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{82,57.5},{102,74.5}})));
  .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor3
    annotation(Placement(transformation(extent = {{68,24.38},{88,44.38}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot4(displayUnits = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{41.72,-68.28},{74.28,-35.72}})));
  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay4(displayMassUnit = true,
      redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingMoistAir)
    annotation(Placement(transformation(extent = {{88,-41.5},{108,-58.5}})));
  .FuelCell.Sensors.CondensingGasMultiDisplaySensor condensingGasMultiDisplaySensor4
    annotation(Placement(transformation(extent = {{90,-13.62},{70,-33.62}})));
    .FuelCell.Sources.CondensingGasFlowBoundary gasFlowSource(use_flow_in = true,T = 298.15,X = {0.005,0.001,0.013,0.752,0.23}) annotation(Placement(transformation(extent = {{-82.0,-34.0},{-62.0,-14.0}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Sources.CondensingGasFlowBoundary gasFlowSource2(use_flow_in = true,T = 343.15,X = {0.162,0.001,0.011,0.633,0.194}) annotation(Placement(transformation(extent = {{132.0,-34.0},{112.0,-14.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.CombiTimeTable dryTable (    
      table = [0, 0.001; 125,0.001;
             250, 0.0016; 375,0.0016;
             500, 0.0022; 625,0.0022;
             750, 0.0028; 875,0.0028;
             1000, 0.0034; 1125,0.0034;
             1250, 0.0038; 1375,0.00383;
             1500, 0.0044; 1625,0.0044;
             1750, 0.005; 1875,0.005],
     timeScale = 1,
     startTime = 0
     ) annotation(Placement(transformation(extent = {{-122.0,-48.0},{-102.0,-28.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.CombiTimeTable wetTable(startTime = 0,timeScale = 1,table = [0, 0.001; 125,0.001;
             250, 0.0016; 375,0.0016;
             500, 0.0022; 625,0.0022;
             750, 0.0028; 875,0.0028;
             1000, 0.0034; 1125,0.0034;
             1250, 0.0038; 1375,0.00383;
             1500, 0.0044; 1625,0.0044;
             1750, 0.005; 1875,0.005]) annotation(Placement(transformation(extent = {{176.0,-52.0},{156.0,-32.0}},origin = {0.0,0.0},rotation = 0.0)));
 
equation

  // Tutte le altre connessioni rimangono identiche
  connect(airCompositionDisplay.data,condensingGasMultiDisplaySensor.u)
    annotation(Line(points={{-38,-43.5},{-38,-33.77},{-27.95,-33.77},{-27.95,-24.05}},color={0,0,0}));
  connect(multiDisplay_phTmdot.y,condensingGasMultiDisplaySensor.u)
    annotation(Line(points={{-2,-52},{-2,-34},{-27.95,-34},{-27.95,-24.05}},color={0,0,0}));
  connect(airCompositionDisplay2.data,condensingGasMultiDisplaySensor2.u)
    annotation(Line(points={{-40,56.5},{-40,46.275},{-26.05,46.275},{-26.05,36.05}},color={0,0,0}));
  connect(multiDisplay_phTmdot2.y,condensingGasMultiDisplaySensor2.u)
    annotation(Line(points={{-4,66},{-4,46},{-26.05,46},{-26.05,36.05}},color={0,0,0}));
  connect(condensingGasMultiDisplaySensor2.portA,discretizedGasGasHumidifier.portB_prim)
    annotation(Line(points={{-20,36},{-9,36},{-9,12},{12,12}},color={209,60,0}));
  connect(condensingGasMultiDisplaySensor.portB,discretizedGasGasHumidifier.portA_sec)
    annotation(Line(points={{-22,-24},{-10,-24},{-10,8},{12,8}},color={209,60,0}));
  connect(discretizedGasGasHumidifier.portB_sec,condensingGasMultiDisplaySensor3.portA)
    annotation(Line(points={{32,8},{50,8},{50,34.38},{72,34.38}},color={209,60,0}));
  connect(airCompositionDisplay3.data,condensingGasMultiDisplaySensor3.u)
    annotation(Line(points={{92,57.5},{92,46.17},{78.05,46.17},{78.05,34.42}},color={0,0,0}));
  connect(condensingGasMultiDisplaySensor3.u,multiDisplay_phTmdot3.y)
    annotation(Line(points={{78.05,34.42},{78.05,46.39},{58,46.39},{58,66}},color={0,0,0}));
  connect(condensingGasMultiDisplaySensor4.portB,discretizedGasGasHumidifier.portA_prim)
    annotation(Line(points={{74,-23.62},{51,-23.62},{51,4},{32,4}},color={209,60,0}));
  connect(airCompositionDisplay4.data,condensingGasMultiDisplaySensor4.u)
    annotation(Line(points={{98,-41.5},{98,-33.61},{79.95,-33.61},{79.95,-23.66}},color={0,0,0}));
  connect(multiDisplay_phTmdot4.y,condensingGasMultiDisplaySensor4.u)
    annotation(Line(points={{58,-52},{58,-33.61},{79.95,-33.61},{79.95,-23.66}},color={0,0,0}));
  connect(drySink.fluidPort,condensingGasMultiDisplaySensor3.portB)
    annotation(Line(points={{109,34.39},{84,34.39}},color={209,60,0}));
  connect(wetSink.fluidPort,condensingGasMultiDisplaySensor2.portB)
    annotation(Line(points={{-63,36},{-32,36}},color={209,60,0}));
    connect(gasFlowSource.fluidPort,condensingGasMultiDisplaySensor.portA) annotation(Line(points = {{-63,-24},{-34,-24}},color = {209,60,0}));
    connect(gasFlowSource2.fluidPort,condensingGasMultiDisplaySensor4.portA) annotation(Line(points = {{113,-24},{100,-24},{100,-23.62},{86,-23.62}},color = {209,60,0}));
    connect(dryTable.y[1],gasFlowSource.m_flow_in) annotation(Line(points = {{-101,-38},{-95,-38},{-95,-9},{-76,-9},{-76,-15}},color = {0,0,127}));
    connect(wetTable.y[1],gasFlowSource2.m_flow_in) annotation(Line(points = {{155,-42},{149,-42},{149,-9},{126,-9},{126,-15}},color = {0,0,127}));

  annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                  graphics={Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100,-100},{100,100}}),
                             Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end Test_1;