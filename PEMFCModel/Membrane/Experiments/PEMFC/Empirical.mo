within PEMFCModel.Membrane.Experiments.PEMFC;
model Empirical "Test PEMFC Empirical membrane"
    
  extends PEMFCModel.Membrane.Experiments.Templates.PartialTestMembrane(
    T_anode = 70 + 273.15,
    p_anode = 1.2e5,
    p_cathode = 1e5,
    X_anode = Medium_an.moleToMassFractions({0.56,0,0,0.2,0.24,0}, Medium_an.MMX),
    X_cathode = Medium_cath.moleToMassFractions({0,0,0.2,0.632,0.168}, Medium_cath.MMX),
    current(I = 0, duration = 0,offset = 450,startTime = 3600),
    redeclare PEMFCModel.Membrane.PEMFC.Empirical cellMembrane(N = 1, C_cell = 50, h_conv_an = fill(20,cellMembrane.N), h_conv_cath = fill(20,cellMembrane.N), pstart = p_anode,n_cell = 1,A_cell = 300e-4,enable_setting = true,z = 50e-6),redeclare replaceable package Medium_an = .PEMFCModel.Medium.NASALowQualityHydrogen);

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
    
end Empirical;
