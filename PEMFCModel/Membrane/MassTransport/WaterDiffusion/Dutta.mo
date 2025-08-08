within PEMFCModel.Membrane.MassTransport.WaterDiffusion;
 model Dutta "Correlation by Dutta for water diffusion"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterDiffusion;
   
protected
  
    .Modelica.Units.SI.DiffusionCoefficient D_lambda[N];

equation
  
  c_an_H2O = rho_dry_m/EW_m*lambda_an;
  c_cath_H2O = rho_dry_m/EW_m*lambda_cath;
  for i in 1:N loop
     
    D_diff[i] = CF_N[i]*D_lambda[i]*exp(2416.0*(1/303.15-1/T_cell[i]));
    D_lambda[i] = Modelon.Math.Smoothing.spliceFunction(
      1.25e-10,
      Modelon.Math.Smoothing.spliceFunction(
        1e-10*(3.0-1.67*(lambda[i]-3.0)),
        Modelon.Math.Smoothing.spliceFunction(
          1e-10*(1.0+2.0*(lambda[i]-2.0)),
          1e-10,
          lambda[i]-2.0,
          0.01),
        lambda[i]-3.0,
        0.01),
      lambda[i]-4.5,
      0.01);

  end for;

end Dutta;
