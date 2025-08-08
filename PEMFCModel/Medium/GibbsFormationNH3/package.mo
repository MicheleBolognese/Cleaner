within PEMFCModel.Medium;
package GibbsFormationNH3
    
    extends .FuelCell.Media.DataDefinitions.GibbsFormation;
    
 constant GibbsData NH3(
     TMIN = 0,
     TMAX = 6000,
     coeff = {0,0,0}) "Ammonia";  
      
end GibbsFormationNH3;
