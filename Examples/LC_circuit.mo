within Examples;

model LC_circuit
    .Modelica.Electrical.Analog.Basic.Ground ground annotation(Placement(transformation(extent = {{46.0,-28.0},{66.0,-8.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 1) annotation(Placement(transformation(extent = {{-42.0,4.0},{-22.0,24.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = 1) annotation(Placement(transformation(extent = {{16.0,72.0},{-4.0,92.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Basic.Inductor inductor(L = 1) annotation(Placement(transformation(extent = {{-32.0,72.0},{-52.0,92.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-32,116},{-52,96}},origin = {0,0},rotation = 0)));
    .Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor2 annotation(Placement(transformation(extent = {{14,116},{-6,96}},origin = {0,0},rotation = 0)));
    .Modelica.Blocks.Math.Add add annotation(Placement(transformation(extent = {{50.0,110.0},{70.0,130.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Interaction.Show.RealValue realValue annotation(Placement(transformation(extent = {{92,110},{112,130}},origin = {0,0},rotation = 0)));
equation
    connect(ground.p,constantVoltage.n) annotation(Line(points = {{56,-8},{56,14},{-22,14}},color = {0,0,255}));
    connect(constantVoltage.p,inductor.n) annotation(Line(points = {{-42,14},{-52,14},{-52,82}},color = {0,0,255}));
    connect(inductor.p,capacitor.n) annotation(Line(points = {{-32,82},{-4,82}},color = {0,0,255}));
    connect(capacitor.p,constantVoltage.n) annotation(Line(points = {{16,82},{22,82},{22,14},{-22,14}},color = {0,0,255}));
    connect(voltageSensor.n,inductor.n) annotation(Line(points = {{-52,106},{-58,106},{-58,82},{-52,82}},color = {0,0,255}));
    connect(voltageSensor.p,inductor.p) annotation(Line(points = {{-32,106},{-26,106},{-26,82},{-32,82}},color = {0,0,255}));
    connect(voltageSensor2.n,capacitor.n) annotation(Line(points = {{-6,106},{-12,106},{-12,82},{-4,82}},color = {0,0,255}));
    connect(voltageSensor2.p,capacitor.p) annotation(Line(points = {{14,106},{22,106},{22,82},{16,82}},color = {0,0,255}));
    connect(voltageSensor.v,add.u1) annotation(Line(points = {{-42,117},{-42,123},{48,123},{48,126}},color = {0,0,127}));
    connect(voltageSensor2.v,add.u2) annotation(Line(points = {{4,117},{4,123},{42,123},{42,114},{48,114}},color = {0,0,127}));
    connect(add.y,realValue.numberPort) annotation(Line(points = {{71,120},{90.5,120}},color = {0,0,127}));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end LC_circuit;
