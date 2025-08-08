within PEMFCModel;
block Stepwise "Generate stepwise signal of type Real"
  
  extends .Modelica.Blocks.Interfaces.SO;
  parameter Integer N(min = 0) = 1 "Number of steps (including offset)";
  parameter Real value[N] = fill(1,N) "Signal value at each step" annotation(Dialog(groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/Step.png"));
  parameter Real duration[N] = fill(1,N) "Duration of each step";
  Integer i_value "Index corresponding to the signal value at a certain time instant";

algorithm

  for i in 1:N loop
  
    if time >= sum(duration[1:i-1]) and time < sum(duration[1:i]) then
    
      i_value := i;
  
    end if;

  end for;

equation
   
    y = value[i_value];

  annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
          Polygon(
            points={{-80,90},{-88,68},{-72,68},{-80,90}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
          Polygon(
            points={{90,-70},{68,-62},{68,-78},{90,-70}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-70},{0,-70},{0,50},{80,50}}),
          Text(
            extent={{-150,-150},{150,-110}},
            textString="")}),
      Documentation(info="<html>
<p>
The Real output y is a step signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/Step.png\"
     alt=\"Step.png\">
</p>

</html>"));

end Stepwise;