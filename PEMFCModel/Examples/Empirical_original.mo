within PEMFCModel.Examples;

model Empirical_original "Test PEMFC Empirical membrane"
  extends .FuelCell.Membranes.Experiments.Templates.PartialTestMembrane(
    T_anode=350,
    p_anode=1.013e5,
    p_cathode=1.013e5,
    X_anode=Medium_an.moleToMassFractions({0.78,0,0.01,0.01,0.2,0,0},Medium_an.MMX),
    X_cathode=Medium_cath.moleToMassFractions({0,0,0.135,0.645,0.22}, Medium_cath.MMX),
    current(I=720, duration=3600,offset =0 ,startTime = 0),
    redeclare .FuelCell.Membranes.PEMFC.Empirical cellMembrane(N = 3,
      redeclare model WaterDiffusion =
          .FuelCell.Membranes.MassTransport.WaterDiffusion.Gandomi,
      redeclare model WaterPermeation =
          .FuelCell.Membranes.MassTransport.WaterPermeation.Generic,
      redeclare model ElectroOsmoticDrag =
          .FuelCell.Membranes.MassTransport.ElectroOsmoticDrag.Generic,redeclare replaceable model GasDiffusion = .FuelCell.Membranes.MassTransport.GasDiffusion.ZeroFlow));

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=3600,
      __Dymola_NumberOfIntervals=600,
      __Dymola_Algorithm="Dassl"),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end Empirical_original;
