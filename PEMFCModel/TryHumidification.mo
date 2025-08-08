within PEMFCModel;
model TryHumidification
    .PEMFCModel.Humidification humidification(T_in = 70 + 273.15,y_dry_in = {0.7,0,0.3,0},redeclare replaceable package Medium = .PEMFCModel4par.Medium.NASASteamH2N2O2) annotation(Placement(transformation(extent = {{-42.0,22.0},{-22.0,42.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression realExpression(y = 1.9e5) annotation(Placement(transformation(extent = {{-108,46},{-88,66}},origin = {0,0},rotation = 0)));
    .Modelica.Blocks.Sources.RealExpression realExpression2(y = 29.59) annotation(Placement(transformation(extent = {{-108.0,10.0},{-88.0,30.0}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Sources.GasFlowBoundary gasFlowSource(use_flow_in = true,use_X_in = true,redeclare replaceable package Medium = .PEMFCModel4par.Medium.NASASteamH2N2O2,T = 70 + 273.15) annotation(Placement(transformation(extent = {{11.342461695456645,5.342461695456645},{64.65753830454335,58.65753830454335}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Sources.GasPressureBoundary gasPressureSource(p = 1.9e5,T = 70 + 237.15,use_X_in = true,redeclare replaceable package Medium = .PEMFCModel4par.Medium.NASASteamH2N2O2,use_p_in = true) annotation(Placement(transformation(extent = {{294.0141864847823,5.469064323478285},{243.4816665171869,56.001584291073605}},origin = {0.0,0.0},rotation = 0.0)));
equation
    connect(realExpression.y,humidification.p_in) annotation(Line(points = {{-87,56},{-65.5,56},{-65.5,38},{-44,38}},color = {0,0,127}));
    connect(realExpression2.y,humidification.V_flow_dry_in) annotation(Line(points = {{-87,20},{-65.5,20},{-65.5,26},{-44,26}},color = {0,0,127}));
    connect(humidification.m_flow_wet_in,gasFlowSource.m_flow_in) annotation(Line(points = {{-19.975,36.1},{-13.975000000000001,36.1},{-13.975000000000001,64.65753830454335},{22.00547701727399,64.65753830454335},{22.00547701727399,58.65753830454335}},color = {0,0,127}));
    connect(humidification.x_wet_in,gasFlowSource.X_in) annotation(Line(points = {{-19.625,28.7},{-13.625,28.7},{-13.625,64.65753830454335},{53.99452298272601,64.65753830454335},{53.99452298272601,58.65753830454335}},color = {0,0,127}));
    connect(gasFlowSource.fluidPort,gasPressureSource.fluidPort) annotation(Line(points = {{61.991784474089016,32},{155.24723057598382,32},{155.24723057598382,30.735324307275945},{246.00829251556667,30.735324307275945}},color = {255,128,0}));
    connect(realExpression.y,gasPressureSource.p_in) annotation(Line(points = {{-87,56},{-81,56},{-81,106.20382352178511},{283.9076824912632,106.20382352178511},{283.9076824912632,56.001584291073605}},color = {0,0,127}));
    connect(humidification.x_wet_in,gasPressureSource.X_in) annotation(Line(points = {{-19.625,28.7},{-13.625,28.7},{-13.625,89.29008591486323},{253.588170510706,89.29008591486323},{253.588170510706,56.001584291073605}},color = {0,0,127}));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end TryHumidification;
