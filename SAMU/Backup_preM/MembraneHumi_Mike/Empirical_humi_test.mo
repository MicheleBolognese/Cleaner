within SAMU.Backup_preM.MembraneHumi_Mike;

model Empirical_humi_test "Test PEMFC Empirical membrane"
  extends .SAMU.Backup_preM.MembraneHumi_Mike.PartialTestumi(
    T_anode=350,
    p_anode=300000,
    p_cathode=250000,
    X_anode=Medium_an.moleToMassFractions({0.65,0,0,0,0.35,0,0}, Medium_an.MMX),
    X_cathode=Medium_cath.moleToMassFractions({0,0,0.135,0.645,0.22}, Medium_cath.MMX),
    redeclare .SAMU.Backup_preM.MembraneHumi_Mike.Empirical_humidificatore_def cellMembrane(
      redeclare model WaterDiffusion =
          .FuelCell.Membranes.MassTransport.WaterDiffusion.Gandomi,
      redeclare model WaterPermeation =
          .FuelCell.Membranes.MassTransport.WaterPermeation.Generic,
                redeclare
        replaceable model                                                                                    GasDiffusion =
              .FuelCell.Membranes.MassTransport.GasDiffusion.Generic (diffusiveSpecies={"N2"})));

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=3600,
      __Dymola_NumberOfIntervals=600,
      __Dymola_Algorithm="Dassl"),
    Documentation(revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end Empirical_humi_test;
