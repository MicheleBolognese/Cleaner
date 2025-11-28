within SAMU.Ema_ShellAndTube;
model test_shellAndTube
    .FuelCell.Sources.CondensingGasFlowBoundary sourceSecondary(m_flow = 0.5,redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingReformateLong) annotation(Placement(transformation(extent = {{-103.9,-13.06},{-83.9,6.94}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Sources.CondensingGasPressureBoundary sinkSecondary(redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingReformateLong) annotation(Placement(transformation(extent = {{2.51,9.01},{-17.49,29.01}},origin = {0.0,0.0},rotation = 0.0)));
    parameter Integer n = 10;
    .FuelCell.Sources.CondensingGasFlowBoundary sourcePrimary(m_flow = 0.5,redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingReformateLong) annotation(Placement(transformation(extent = {{3.24,-14.05},{-16.76,5.95}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Sources.CondensingGasPressureBoundary sinkPrimary(redeclare replaceable package Medium = FuelCell.Media.PreDefined.CondensingGases.CondensingReformateLong) annotation(Placement(transformation(extent = {{-103.97,9.96},{-83.97,29.96}},origin = {0.0,0.0},rotation = 0.0)));
    .SAMU.Ema_ShellAndTube.DiscretizedGasGasHumidifier discretizedGasGasHumidifier_Ema(n = n,useHeatTransfer = false) annotation(Placement(transformation(extent = {{-66.0,6.0},{-46.0,26.0}},origin = {0.0,0.0},rotation = 0.0)));
equation
    connect(sourceSecondary.fluidPort,discretizedGasGasHumidifier_Ema.portA_sec) annotation(Line(points = {{-84.9,-3.06},{-84.9,16},{-66,16}},color = {209,60,0}));
    connect(discretizedGasGasHumidifier_Ema.portB_sec,sinkSecondary.fluidPort) annotation(Line(points = {{-46,16},{-16.49,16},{-16.49,19.01}},color = {209,60,0}));
    connect(sourcePrimary.fluidPort,discretizedGasGasHumidifier_Ema.portA_prim) annotation(Line(points = {{-15.76,-4.05},{-28.5,-4.05},{-28.5,12},{-46,12}},color = {209,60,0}));
    connect(discretizedGasGasHumidifier_Ema.portB_prim,sinkPrimary.fluidPort) annotation(Line(points = {{-66,20},{-66,19.96},{-84.97,19.96}},color = {209,60,0}));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end test_shellAndTube;
