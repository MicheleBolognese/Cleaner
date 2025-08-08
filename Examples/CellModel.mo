within Examples;
model CellModel
    .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(Placement(transformation(extent = {{26.0,56.0},{46.0,76.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(Placement(transformation(extent = {{26,-74},{46,-54}},origin = {0,0},rotation = 0)));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end CellModel;
