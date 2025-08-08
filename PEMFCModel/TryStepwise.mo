within PEMFCModel;
model TryStepwise
    .PEMFCModel.Stepwise stepwise(N = 4,duration = {1125,450,450,1575},value = {1.4,1.5,1.9,2.2}) annotation(Placement(transformation(extent = {{-70,14},{-50,34}},origin = {0,0},rotation = 0)));
    annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
end TryStepwise;
