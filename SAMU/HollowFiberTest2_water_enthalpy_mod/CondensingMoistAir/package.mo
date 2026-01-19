within SAMU.HollowFiberTest2_water_enthalpy_mod;

package CondensingMoistAir
  "Moist air with condensed water fraction: {H2O, CO2, Ar, N2, O2}"
  extends .Modelon.Icons.Fluid;
  extends .Modelon.Media.Templates.IdealGasWithCondensingLiquid(
    redeclare package EquationsOfState =
        .FuelCell.Media.EquationsOfState.Examples.GasAndLiquid.CondensingMoistAir,
    redeclare package TransportProperties =
      .Modelon.Media.TransportProperties.Templates.Nasa.IdealGasNasa7TransportProperties
        (data = {.Modelon.Media.DataDefinitions.NasaData.H2O, .Modelon.Media.DataDefinitions.NasaData.CO2, .Modelon.Media.DataDefinitions.NasaData.Ar, .Modelon.Media.DataDefinitions.NasaData.N2, .Modelon.Media.DataDefinitions.NasaData.O2},
         substanceNames=substanceNames,
         fluidConstants = fluidConstants),
    fluidConstants = {.Modelon.Media.DataDefinitions.FluidData.H2O, .Modelon.Media.DataDefinitions.FluidData.CO2, .Modelon.Media.DataDefinitions.FluidData.Ar, .Modelon.Media.DataDefinitions.FluidData.N2, .Modelon.Media.DataDefinitions.FluidData.O2},
    reference_X={0.01,0.005,0.005,0.75,0.23},
    Tfreezing = 273.15,
    substanceNames={.Modelon.Media.DataDefinitions.NasaData.H2O.name, .Modelon.Media.DataDefinitions.NasaData.CO2.name, .Modelon.Media.DataDefinitions.NasaData.Ar.name, .Modelon.Media.DataDefinitions.NasaData.N2.name, .Modelon.Media.DataDefinitions.NasaData.O2.name},
    MMX=EquationsOfState.Gas.MMX,
    R=EquationsOfState.Gas.R,
    hasCriticalData=true,
    final limits = EquationsOfState.limits,
    CompressibilityType = .Modelon.Media.Interfaces.Types.Compressibility.FullyCompressible,
    analyticInverseTfromh = EquationsOfState.analyticInverseTfromh,
    ThermoStates = .Modelon.Media.Interfaces.Types.IndependentVariables.pTX,
    FixedComposition=false);

  annotation (Protection(access=Access.packageDuplicate), Documentation(revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end CondensingMoistAir_1;
