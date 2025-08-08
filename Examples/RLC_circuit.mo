within Examples;

model RLC_circuit
    .Modelica.Electrical.Analog.Basic.Resistor resistor(useHeatPort = false,R = 1) annotation(Placement(transformation(extent = {{-38.0,10.0},{-58.0,30.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Basic.Ground ground annotation(Placement(transformation(extent = {{50,0},{70,20}},origin = {0,0},rotation = 0)));
    .Modelica.Electrical.Analog.Basic.Inductor inductor(L = 1) annotation(Placement(transformation(extent = {{-60.0,56.0},{-40.0,76.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = 1) annotation(Placement(transformation(extent = {{-10.0,-10.0},{10.0,10.0}},origin = {-66.0,48.0},rotation = 90.0)));
    .Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 1) annotation(Placement(transformation(extent = {{46.0,56.0},{26.0,76.0}},origin = {0.0,0.0},rotation = 0.0)));
equation
    connect(inductor.p,capacitor.n) annotation(Line(points = {{-60,66},{-66,66},{-66,58}},color = {0,0,255}));
    connect(capacitor.p,resistor.n) annotation(Line(points = {{-66,38},{-66,20},{-58,20}},color = {0,0,255}));
    connect(resistor.p,ground.p) annotation(Line(points = {{-38,20},{60,20}},color = {0,0,255}));
    connect(constantVoltage.p,ground.p) annotation(Line(points = {{46,66},{60,66},{60,20}},color = {0,0,255}));
    connect(inductor.n,constantVoltage.n) annotation(Line(points = {{-40,66},{26,66}},color = {0,0,255}));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end RLC_circuit;
