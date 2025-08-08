within Examples;
model RC_circuit
    .Modelica.Electrical.Analog.Basic.Resistor resistor(useHeatPort = false,R = 1) annotation(Placement(transformation(extent = {{-58.0,10.0},{-38.0,30.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Basic.Ground ground annotation(Placement(transformation(extent = {{50,0},{70,20}},origin = {0,0},rotation = 0)));
    .Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 1) annotation(Placement(transformation(extent = {{26.0,56.0},{46.0,76.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 0) annotation(Placement(transformation(extent = {{-38.0,34.0},{-58.0,54.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch(Ron = 0,Goff = 0) annotation(Placement(transformation(extent = {{-10.0,-10.0},{10.0,10.0}},origin = {-82.0,44.0},rotation = -90.0)));
    .Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = 1) annotation(Placement(transformation(extent = {{-42.0,56.0},{-62.0,76.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Ideal.IdealOpeningSwitch switch2(Ron = 0,Goff = 0) annotation(Placement(transformation(extent = {{-10.0,-10.0},{10.0,10.0}},origin = {60.0,46.0},rotation = -90.0)));
    .Modelica.Blocks.Sources.BooleanStep booleanStep2(startTime = 5) annotation(Placement(transformation(extent = {{126.0,36.0},{106.0,56.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch3(Ron = 0,Goff = 0) annotation(Placement(transformation(extent = {{-10.0,-10.0},{10.0,10.0}},origin = {-6.0,44.0},rotation = -90.0)));
    .Modelica.Blocks.Sources.BooleanStep booleanStep3(startTime = 5) annotation(Placement(transformation(extent = {{40.0,34.0},{20.0,54.0}},origin = {0.0,0.0},rotation = 0.0)));
equation
    connect(ground.p,resistor.n) annotation(Line(points = {{60,20},{-38,20}},color = {0,0,255}));
    connect(switch.n,resistor.p) annotation(Line(points = {{-82,34},{-82,20},{-58,20}},color = {0,0,255}));
    connect(booleanStep.y,switch.control) annotation(Line(points = {{-59,44},{-70,44}},color = {255,0,255}));
    connect(constantVoltage.p,capacitor.p) annotation(Line(points = {{26,66},{-42,66}},color = {0,0,255}));
    connect(capacitor.n,switch.p) annotation(Line(points = {{-62,66},{-82,66},{-82,54}},color = {0,0,255}));
    connect(booleanStep2.y,switch2.control) annotation(Line(points = {{105,46},{72,46}},color = {255,0,255}));
    connect(constantVoltage.n,switch2.p) annotation(Line(points = {{46,66},{60,66},{60,56}},color = {0,0,255}));
    connect(switch2.n,ground.p) annotation(Line(points = {{60,36},{60,20}},color = {0,0,255}));
    connect(constantVoltage.p,switch3.p) annotation(Line(points = {{26,66},{-6.000000000000002,66},{-6.000000000000002,54}},color = {0,0,255}));
    connect(switch3.n,resistor.n) annotation(Line(points = {{-5.999999999999998,34},{-5.999999999999998,20},{-38,20}},color = {0,0,255}));
    connect(booleanStep3.y,switch3.control) annotation(Line(points = {{19,44},{6,44}},color = {255,0,255}));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end RC_circuit;
