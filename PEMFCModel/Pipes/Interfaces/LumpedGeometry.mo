within PEMFCModel.Pipes.Interfaces;
partial model LumpedGeometry "Geometry parameters for lumped circular flow channel or pipe"
  
  constant Real pi = .Modelica.Constants.pi;

  parameter Real n_channels(min = 1.0) = 1.0 "Number of parallel channels" annotation(Dialog(group="Geometry"));
  parameter .Modelica.Units.SI.Length L "Channel length" annotation (Dialog(group="Geometry"));
  parameter .Modelica.Units.SI.Length D "Diameter of single channel" annotation (Dialog(group="Geometry"));

  final parameter .Modelica.Units.SI.Length Dhyd = .Modelon.ThermoFluid.Functions.CharacteristicNumbers.HydraulicDiameter(A, C) "Hydraulic diameter of single channel" annotation (Dialog(group="Geometry"));
  final parameter .Modelica.Units.SI.Area A = pi*D*D/4 "Cross sectional area of single channel" annotation (Dialog(group="Geometry"));
  final parameter .Modelica.Units.SI.Length C = pi*D "Circumference of single channel" annotation (Dialog(group="Geometry"));
  final parameter .Modelica.Units.SI.Volume V = n_channels*A*L "Total volume of all channels" annotation (Dialog(group="Geometry"));

  final parameter .Modelica.Units.SI.Area A_heat = C*L "Heat transfer area of single channel" annotation (Dialog(group="Heat transfer"));

end LumpedGeometry;