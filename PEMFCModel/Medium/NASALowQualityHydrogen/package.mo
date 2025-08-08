within PEMFCModel.Medium;
package NASALowQualityHydrogen "NASA industrial quality hydrogen gas: H2, CO, CO2, H2O, N2, O2"
    
    extends .FuelCell.Media.PreDefined.IdealGases.NASAReformateShort(reference_X = moleToMassFractions({0.48,0.006,0.07,0.162,0,0},MMX));
    
end NASALowQualityHydrogen;
