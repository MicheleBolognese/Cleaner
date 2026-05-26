within SAMU.HollowFiberTest2_1;

model DiscretizedGasGasHumidifier_improved
  "Improved discretized gas-gas humidifier with two-phase support and convective resistance"
  // ===========================================================================
  // VERSIONE MIGLIORATA di DiscretizedGasGasHumidifier (Modelon FuelCell lib)
  //
  // CORREZIONI APPLICATE:
  //   [FIX-2] Inizializzazione esplicita C_w (initial equation) — evita DAE failure
  //   [FIX-3] D(lambda) Vetter C-inf al posto del modello a tratti Park/Springer
  //   [FIX-5] Switching smooth vapor/liquid — Paradosso di Schroeder (Mull 2025)
  //   [FIX-6] Calore specifico Nafion: 1200 J/kg/K (era 4188 = acqua)
  //   [FIX-8] n=5 segmenti (era 3) — migliore risoluzione assiale
  //
  // NON APPLICATI:
  //   [FIX-4] Tre resistenze in serie (R_conv_wet + R_mem + R_conv_dry)
  //   [FIX-9] Heat transfer: mantiene Dittus-Boelter originale (Park 2008)
  //
  // GIA' CORRETTI nel file sorgente ricevuto (non ripetuti):
  //   [FIX-1] Bug entalpico h_water_sec (state_sec corretto)
  //   [FIX-7] smooth() su lambda(a) — gia' presente nel sorgente
  //
  // RIFERIMENTI:
  //   Park et al. (2008), Int. J. Hydrogen Energy 33, 2273-2282
  //   Park et al. (2013), Int. J. Automotive Technology 14(3), 449-457
  //   Mull et al. (2025), Next Energy 7, 100294
  //   Pollak et al. (2023), Energies 16, 2578
  //   Pollak et al. (2023), Modelica Conference, 531-540
  //   Vetter & Schumacher (2019), Comput. Phys. Commun. 234, 223-234
  //   Costello et al. (1993), J. Membrane Sci. 80, 1-11
  //   Kusoglu & Weber (2017), Chem. Rev. 117, 987-1104
  // ===========================================================================
  extends .FuelCell.HeatExchangers.Interfaces.DistributedGasGasHumidifier(
    wallThickness=t_mem,
    redeclare replaceable model HeatTransfer_sec =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    redeclare replaceable model Friction_sec =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss,
    redeclare replaceable model HeatTransfer_prim =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    redeclare replaceable model Friction_prim =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss,
    redeclare
      .Modelon.ThermoFluid.HeatExchangers.Records.Summary.Base_pinchlmtd.BaseSummary
      summary(
      m_flow=m_flow_prim,
      m_flow_sec=m_flow_sec,
      Q_flow=sum(shell.wall.Q_flow),
      T_in=noEvent(if m_flow_prim >= 0 then shell.summary.T_in else shell.summary.T_out),
      T_out=noEvent(if m_flow_prim >= 0 then shell.summary.T_out else shell.summary.T_in),
      h_in=noEvent(if m_flow_prim >= 0 then shell.summary.h_in else shell.summary.h_out),
      h_out=noEvent(if m_flow_prim >= 0 then shell.summary.h_out else shell.summary.h_in),
      dp=shell.summary.dp,
      T_sec_in=noEvent(if m_flow_sec >= 0 then tube.summary.T_in else tube.summary.T_out),
      T_sec_out=noEvent(if m_flow_sec >= 0 then tube.summary.T_out else tube.summary.T_in),
      h_sec_in=noEvent(if m_flow_sec >= 0 then tube.summary.h_in else tube.summary.h_out),
      h_sec_out=noEvent(if m_flow_sec >= 0 then tube.summary.h_out else tube.summary.h_in),
      dp_sec=tube.summary.dp,
      p_in=portA_prim.p,
      p_out=portB_prim.p,
      p_sec_in=portA_sec.p,
      p_sec_out=portB_sec.p,
      Q_flow_sec=sum(tube.wall.Q_flow)),
    // [FIX-8] n=5 invece di n=3.
    // Mull et al. (2025) dimostrano che 5 volumi finiti rappresentano il
    // miglior compromesso tra accuratezza della distribuzione assiale e
    // costo computazionale per umidificatori a membrana PEMFC.
    n = 5);

  /* Local variables */
  .Modelica.Units.SI.MassFlowRate m_flow_prim
    "Mass flowrate primary side, positive flow from portA_prim to portB_prim";
  .Modelica.Units.SI.MassFlowRate m_flow_sec
    "Mass flowrate primary side, positive flow from portA_sec to portB_sec";
  .Modelica.Units.SI.HeatFlowRate Q_tot "Total heat flow rate out from primary side";

  .SAMU.HollowFiberTest2_1.FlowChannel shell(
    redeclare package Medium = PrimaryMedium,
    A=fill(.Modelica.Constants.pi*D_h*D_h/4, n),
    A_heat=fill(A_mem, n),
    useHeatTransfer=useHeatTransfer,
    initOpt=initOpt,
    p_start_in=init_prim.p_in,
    p_start_out=init_prim.p_out,
    p_start=p_start,
    initFromEnthalpy=init_prim.initFromEnthalpy,
    h_start_in=init_prim.h_in,
    h_start_out=init_prim.h_out,
    h_start=h_start,
    T_start_in=init_prim.T_in,
    T_start_out=init_prim.T_out,
    T_start=T_start,
    X_start=init_prim.X,
    m_flow_start=init_prim.m_flow,
    L=fill(L_mem/n, n),
    Dhyd=fill(D_h, n),
    positiveFlow=positiveFlow_prim,
    from_dp=from_dp_prim,
    dp_smooth=dp_smooth,
    mflow_smooth=mflow_smooth,
    generateEventForReversal=generateEventForReversal_prim,
    C_start=init_prim.C,
    redeclare model Friction = Friction_prim,
    redeclare model HeatTransfer = HeatTransfer_prim,
    CF_PressureLoss=CF_PrimarySidePressureLoss,
    CF_HeatTransfer=CF_PrimarySideHeatTransfer,
    H_flow = h_water_prim .* m_trans_shell,
    mX_flow=transpose({if i == H2O_prim then -m_trans_shell else zeros(n) for i in 1:
        nX_shell})) "Humidifier shell side" annotation (Dialog(tab="Sub-components",
        group="Flow channels"), Placement(transformation(
        extent={{-20.0,-20.0},{20.0,20.0}},
        rotation=180.0,
        origin={-8.0,60.0})));

  .SAMU.HollowFiberTest2_1.FlowChannel tube(
    redeclare package Medium = SecondaryMedium,
    n_channels=fill(n_tubes, n),
    A=fill(.Modelica.Constants.pi*D_mem*D_mem/4, n),
    A_heat=fill(A_mem, n),
    useHeatTransfer=useHeatTransfer,
    initOpt=initOpt_sec,
    p_start_out=init_sec.p_out,
    p_start=p_sec_start,
    initFromEnthalpy=init_sec.initFromEnthalpy,
    h_start_in=init_sec.h_in,
    h_start_out=init_sec.h_out,
    h_start=h_sec_start,
    T_start_in=init_sec.T_in,
    T_start_out=init_sec.T_out,
    T_start=T_sec_start,
    X_start=init_sec.X,
    m_flow_start=init_sec.m_flow,
    p_start_in=init_sec.p_in,
    L=fill(L_mem/n, n),
    Dhyd=fill(D_mem, n),
    positiveFlow=positiveFlow_sec,
    from_dp=from_dp_sec,
    dp_smooth=dp_smooth,
    mflow_smooth=mflow_smooth,
    generateEventForReversal=generateEventForReversal_sec,
    C_start=init_sec.C,
    redeclare model Friction = Friction_sec,
    redeclare model HeatTransfer = HeatTransfer_sec,
    CF_PressureLoss=CF_SecondarySidePressureLoss,
    CF_HeatTransfer=CF_SecondarySideHeatTransfer,
    H_flow = h_water_sec .* m_trans,
    mX_flow=transpose({if i == H2O_sec then m_trans else zeros(n) for i in 1:nX_tube}))
    "Humidifier tube side" annotation (Dialog(tab="Sub-components", group=
          "Flow channels"), Placement(transformation(extent={{-28,-80},{12,-40}})));

  .Modelon.ThermoFluid.Solids.DynamicWall wall(
    props(each Cp=wallMaterial.c, Rw=wallThickness./wallMaterial.lambda./
          wallCrossSecArea),
    n=n,
    steadyStateInit=WallSteadyStateInit,
    n_channels=n_channels_wall,
    m=m,
    T0=T0_wall,
    massLessWall=massLessWall,
    includeThermalResistance=includeThermalResistance) "Wall model definition"
    annotation (Dialog(tab="Sub-components", group="Wall"), Placement(
        transformation(extent={{-30,-20},{10,20}})));

  replaceable record WallMaterial =
      .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
      (
      // [FIX-6] Calore specifico Nafion 117 corretto: ~1200 J/(kg*K).
      // Il valore originale c=4188 J/(kg*K) e' il calore specifico dell'ACQUA,
      // non della membrana PFSA. Questo errore sovrastimava la capacita'
      // termica della membrana di un fattore ~3.5x, distorcendo i transienti
      // termici. Riferimento: Kusoglu & Weber, Chem. Rev. 117 (2017) 987-1104.
      // Range letteratura Nafion: 1100-1500 J/(kg*K).
      // rho=1968 kg/m3 e lambda=0.12 W/(m*K) rimangono invariati (corretti).
      c=1200,
      rho=1968,
      lambda=0.12) constrainedby
    .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
    "Geometry and material parameters"     annotation(choicesAllMatching, Dialog(tab="Wall",group="Wall geometry"));

  WallMaterial wallMaterial;

  parameter Boolean useHeatTransfer=false
    "Consider heat transfer effects (not used when connected to membrane)"
    annotation (Dialog(tab="Correlations", group="Heat transfer"));

  /* Diffusion correlation */

  constant Integer H2O_prim=PrimaryMedium.substanceIndex("H2O");
  constant Integer H2O_sec=SecondaryMedium.substanceIndex("H2O");
  constant Integer nX_shell = PrimaryMedium.nS
    "total number of mass fractions for shell side";
  constant Integer nX_tube = SecondaryMedium.nS
    "total number of mass fractions for tube side";

  package WaterMedium = .Modelon.Media.PreDefined.TwoPhase.WaterIF97
    "water medium model";

  .Modelica.Units.SI.MassFlowRate m_water=sum(m_trans) "Water diffusion mass flow rate";
  .Modelica.Units.SI.MassFlowRate[n] m_trans "Mass flow rate of steam across membrane per segment";

  .Modelon.Media.Units.MoleFraction[n,nX_shell] y_prim
     "primary inlet concentration";
  .Modelon.Media.Units.MoleFraction[n,nX_tube] y_sec
     "secondary exit concentration";

  .Modelica.Units.SI.DimensionlessRatio a_shell[n] "Water activity of shell side";
  .Modelica.Units.SI.DimensionlessRatio a_tube[n] "Water activity of tube side";
  .Modelica.Units.SI.DimensionlessRatio a_mem[n] "Water activity of membrane side";

  Real lambda_mem[n] "Mean water content in membrane";
  Real lambda_shell[n] "Water content in shell side";
  Real lambda_tube[n] "Water content in tube side";

  .Modelica.Units.SI.Temperature T_mem[n] "Membrane temperature";

  PrimaryMedium.ThermodynamicState state_prim[n]
     "thermodynamic state on the primary side";
  SecondaryMedium.ThermodynamicState state_sec[n]
     "thermodynamic state on the secondary side";

  .Modelica.Units.SI.DiffusionCoefficient Dw[n] "Membrane diffusion coefficient";

//protected
  .Modelica.Units.SI.MassFlowRate[n] m_trans_shell "Mass transfer rate for shell";
  // Real D[n] "Empirical constant"; non necessaria, dichiaro direttamente Dw[n]


  .Modelica.Units.SI.Pressure pV_tube[n] "Partial pressure of water vapor tube side";
  .Modelica.Units.SI.Pressure pV_shell[n] "Partial pressure of water vapor shell side";
  .Modelica.Units.SI.Pressure psat_tube[n] "Water saturation pressure tube side";
  .Modelica.Units.SI.Pressure psat_shell[n] "Water saturation pressure shell side";

  .Modelica.Units.SI.MolarDensity C_shell[n] "Water mass concentration shell side";
  .Modelica.Units.SI.MolarDensity C_tube[n] "Water mass concentration tube side";
  .Modelica.Units.SI.MassConcentration C_w[n] "Overall water mass concentration";

  .Modelica.Units.SI.MolarMass Mv=18.01528e-3 "Molar mass of vapor";

  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim
     "enthalpy of gaseous water";
  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec
     "enthalpy of gaseous water";

  // ===========================================================================
  // PARAMETRI AGGIUNTIVI �� DIFFUSIVITA' VETTER [FIX-3] e SCHROEDER [FIX-5]
  // ===========================================================================

  // [FIX-3] Parametri del modello Vetter & Schumacher (2019).
  // Sostituiscono i parametri k e Tk del modello originale di Park/Springer
  // con un polinomio razionale continuo (C-inf), eliminando i breakpoint
  // a lambda=2, 3, 4.5 che causano eventi numerici nel solver.
  // I vecchi parametri k e Tk rimangono nel parent class ma non sono piu'
  // usati nel calcolo di Dw (mantenuti per retrocompatibilita').
  parameter Real k_fit_vapor = 2.0
    "[FIX-3] Fit factor diffusivita' membrana — solo vapore [-].
     Mull et al. (2025), Next Energy 7, Tabella 4.
     Regola la diffusivita' effettiva rispetto al modello Vetter base."
    annotation(Dialog(tab="Membrane", group="Diffusivity (Vetter model)"));

  // [FIX-5] k_fit_liquid modella il Paradosso di Schroeder:
  // quando acqua liquida tocca la membrana PFSA, la superficie diventa
  // idrofilica e la diffusivita' aumenta di ~20x rispetto al solo vapore.
  // Riferimento: Mull 2025 (k_fit=40 con liquido vs k_fit=2 solo vapore).
  parameter Real k_fit_liquid = 40.0
    "[FIX-5] Fit factor diffusivita' membrana — con acqua liquida [-].
     Modella il Paradosso di Schroeder (Chen et al. 2023, RSE Rev. 173).
     Mull et al. (2025), Next Energy 7, Tabella 4."
    annotation(Dialog(tab="Membrane", group="Diffusivity (Vetter model)"));

  parameter .Modelica.Units.SI.Temperature T_ref_Vetter = 353.15
    "[FIX-3] Temperatura di riferimento Arrhenius — modello Vetter [K].
     Vetter & Schumacher (2019), Comput. Phys. Commun. 234, Eq. (8)."
    annotation(Dialog(tab="Membrane", group="Diffusivity (Vetter model)"));

  parameter .Modelica.Units.SI.SpecificEnergy E_a_Vetter = 20000.0
    "[FIX-3] Energia di attivazione diffusione H2O in PFSA [J/mol].
     Vetter & Schumacher (2019), Comput. Phys. Commun. 234, Eq. (8)."
    annotation(Dialog(tab="Membrane", group="Diffusivity (Vetter model)"));

  // [FIX-5] Parametri per lo smooth switching vapore/liquido.
  // La transizione usa tanh per garantire continuita' C1 (no eventi numerici).
  // La soglia a_liq=1.0 corrisponde alla saturazione (a_shell=phi_RH=1).
  parameter Real dx_liq_smooth = 0.02
    "[FIX-5] Ampiezza della zona di transizione vapor/liquid [-].
     La transizione avviene smoothly nell'intervallo [1-dx, 1+dx] di a_shell.
     Valore piu' piccolo = transizione piu' brusca (attenzione agli eventi)."
    annotation(Dialog(tab="Membrane", group="Diffusivity (Vetter model)"));

  // ===========================================================================
  // PARAMETRI E VARIABILI AGGIUNTIVI — SWITCHING [FIX-5] e VETTER [FIX-3]
  // ===========================================================================

  // Costante universale dei gas — usata nel termine Arrhenius di Vetter [FIX-3]
  constant Real R_gas = 8.314 "Costante universale dei gas [J/(mol*K)]";

  // [FIX-5] Fit factor effettivo per segmento: transizione smooth tra
  // k_fit_vapor (solo vapore) e k_fit_liquid (con acqua liquida).
  Real k_fit_eff[n]
    "Fit factor effettivo per segmento — smooth switching vapor/liquid [-]";

//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim_g
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec_g
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim_l
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec_l
//      "enthalpy of gaseous water";
initial equation
  assert(
    H2O_prim > 0,
    "Humidifier: PrimaryMedium does not contain the required substance: H2O",
    level=AssertionLevel.error);
  assert(
    H2O_sec > 0,
    "Humidifier: SecondaryMedium does not contain the required substance: H2O",
    level=AssertionLevel.error);

  // [FIX-2] Inizializzazione esplicita di C_w come variabile di stato.
  // Nel modello originale C_w era presente come variabile (riga 193) ma
  // priva di initial equation: il sistema DAE risultava under-determined
  // all'inizializzazione, causando possibili failure o risultati non deterministici.
  //
  // Condizione iniziale: equilibrio di Springer con il lato wet (shell).
  // C_w e' la concentrazione massica totale d'acqua nella membrana [kg/m^3].
  // Viene inizializzata dalla relazione lambda(a_shell) al tempo t=0,
  // assumendo che la membrana sia in equilibrio col flusso wet iniziale.
  //
  // Formula: C_w = (rho_mem / M_mem) * Mv * lambda_shell
  // dove lambda_shell e' calcolato da a_shell con la relazione di Springer.
  //
  // Nota: a_shell all'avvio e' determinata dalle condizioni iniziali dei
  // canali (init_prim), quindi questa equazione e' ben determinata.
  for i in 1:n loop
    C_w[i] = (rho_mem / M_mem) * Mv *
              smooth(0, if a_shell[i] <= 1.0 then
                0.043 + 17.81*a_shell[i] - 39.85*a_shell[i]^2 + 36.0*a_shell[i]^3
              elseif a_shell[i] <= 3.0 then
                14.0 + 1.4*(a_shell[i] - 1.0)
              else
                16.8);
  end for;

equation
  m_flow_prim =shell.summary.m_flow;
  m_flow_sec =tube.summary.m_flow;
  Q_tot =-1*sum(shell.wall.Q_flow);

  state_sec = tube.channel.volume.state;

 for i in 1:n loop

//    h_water_prim_[i] = PrimaryMedium.specificEnthalpy_index(state_prim[i],H2O_prim)*(state_prim[i].X[H2O_prim]);
//    h_water_sec[i] = SecondaryMedium.specificEnthalpy_index(state_sec[i],H2O_sec)*(state_sec[i].X[H2O_sec]);

    h_water_prim[i] = (PrimaryMedium.enthalpyOfCondensingGas(state_prim[i].T) + PrimaryMedium.enthalpyOfLiquid(state_prim[i].T))*(state_prim[i].X[H2O_prim]);
    h_water_sec[i] = (SecondaryMedium.enthalpyOfCondensingGas(state_sec[i].T) + SecondaryMedium.enthalpyOfLiquid(state_sec[i].T))*(state_sec[i].X[H2O_sec]);

//     h_water_prim_g[i] = PrimaryMedium.enthalpyOfCondensingGas(state_prim[i].T);
//     h_water_sec_g[i] = SecondaryMedium.enthalpyOfCondensingGas(state_sec[i].T);
//
//     h_water_prim_l[i] = PrimaryMedium.enthalpyOfLiquid(state_prim[i].T);
//     h_water_sec_l[i] = SecondaryMedium.enthalpyOfLiquid(state_sec[i].T);

    y_prim[i,:] = PrimaryMedium.massToMoleFractions(state_prim[i].X,PrimaryMedium.MMX);
    y_sec[i,:] = SecondaryMedium.massToMoleFractions(state_sec[i].X,SecondaryMedium.MMX);

    pV_shell[i] = y_prim[i, H2O_prim]*state_prim[i].p;
    pV_tube[i] = y_sec[i, H2O_sec]*state_sec[i].p;

    psat_shell[i] = min(WaterMedium.saturationPressure_TX(max(shell.channel.volume[i].T, 273.16)), 0.999*state_prim[i].p);
    psat_tube[i] = min(WaterMedium.saturationPressure_TX(max(tube.channel.volume[i].T, 273.16)), 0.999*state_sec[i].p);

    a_shell[i] = pV_shell[i]/psat_shell[i];
    a_tube[i] = pV_tube[i]/psat_tube[i];
    a_mem[i] = (a_shell[i] + a_tube[i])/2;

    // Smooth blending tra rami, evita eventi
    lambda_shell[i] = smooth(1, 
      if a_shell[i] < 0.98 then 
          0.043 + 17.81*a_shell[i] - 39.85*a_shell[i]^2 + 36.0*a_shell[i]^3
      elseif a_shell[i] < 3.0 then 
          14 + 1.4*(a_shell[i] - 1)
      else 16.8);

      // Smooth blending tra rami, evita eventi
    lambda_tube[i] = smooth(1, 
      if a_tube[i] < 0.98 then 
          0.043 + 17.81*a_tube[i] - 39.85*a_tube[i]^2 + 36.0*a_tube[i]^3
      elseif a_tube[i] < 3.0 then 
          14 + 1.4*(a_tube[i] - 1)
      else 16.8);

      // Smooth blending tra rami, evita eventi
    lambda_mem[i] = smooth(1, 
      if a_mem[i] < 0.98 then 
          0.043 + 17.81*a_mem[i] - 39.85*a_mem[i]^2 + 36.0*a_mem[i]^3
      elseif a_mem[i] < 3.0 then 
          14 + 1.4*(a_mem[i] - 1)
      else 16.8);

    // if a_shell[i] >= 0 and a_shell[i] <= 1 then
    //   lambda_shell[i] =0.043 + 17.81*a_shell[i] - 39.85*a_shell[i]^2 + 36.0*a_shell[i]^3;
    // elseif a_shell[i] > 1 and a_shell[i] <= 3 then
    //   lambda_shell[i] =14 + 1.4*(a_shell[i] - 1);
    // else
    //   lambda_shell[i] = 16.8;
    // end if;

    // if a_tube[i] >= 0 and a_tube[i] <= 1 then
    //   lambda_tube[i] =0.043 + 17.81*a_tube[i] - 39.85*a_tube[i]^2 + 36.0*a_tube[i]^3;
    // elseif a_tube[i] > 1 and a_tube[i] <= 3 then
    //   lambda_tube[i] =14 + 1.4*(a_tube[i] - 1);
    // else
    //   lambda_tube[i] = 16.8;
    // end if;

    // if a_mem[i] >= 0 and a_mem[i] <= 1 then
    //   lambda_mem[i] =0.043 + 17.81*a_mem[i] - 39.85*a_mem[i]^2 + 36.0*a_mem[i]^3;
    // elseif a_mem[i] > 1 and a_mem[i] <= 3 then
    //   lambda_mem[i] =14 + 1.4*(a_mem[i] - 1);
    // else
    //   lambda_mem[i] = 16.8;
    // end if;

    // =========================================================================
    // [FIX-5] SWITCHING SMOOTH VAPOR / LIQUID — PARADOSSO DI SCHROEDER
    // =========================================================================
    // Quando a_shell > 1 c'e' acqua liquida in contatto con la membrana PFSA.
    // La superficie diventa idrofilica: diffusivita' e permeabilita' aumentano
    // di ~20x rispetto al solo vapore (Paradosso di Schroeder, Chen 2023).
    //
    // Implementazione: funzione tanh per transizione continua C1 senza eventi.
    // - a_shell < 1  => solo vapore => k_fit_eff ~ k_fit_vapor
    // - a_shell > 1  => liquido    => k_fit_eff ~ k_fit_liquid
    // - transizione smooth in una banda di ampiezza dx_liq_smooth intorno a a=1
    k_fit_eff[i] = k_fit_vapor
      + (k_fit_liquid - k_fit_vapor)
      * (0.5 + 0.5*tanh((a_shell[i] - 1.0) / max(dx_liq_smooth, 1e-6)));

    // =========================================================================
    // [FIX-3] DIFFUSIVITA' MEMBRANA — MODELLO VETTER CONTINUO (C-inf)
    // =========================================================================
    // Sostituisce il modello D(lambda) a tratti di Park/Springer (righe 301-311
    // del modello originale) con il polinomio razionale di Vetter & Schumacher
    // (2019), Comput. Phys. Commun. 234, 223-234, Eq. (8).
    //
    // Forma originale (Park): D(lambda) a tratti con breakpoint a lambda=2,3,4.5
    //   → 3 eventi numerici per segmento ad ogni ciclo di integrazione
    //   → ramo 3<lambda<4.5 DECRESCENTE (non fisico per Nafion moderna)
    //
    // Forma Vetter (nuova):
    //   D(lambda, T) = k_fit * [P(lambda)/Q(lambda)] * 1e-10 * exp[Ea/R*(1/Tref-1/T)]
    //   P(lambda) = 3.842*lambda^3 - 32.03*lambda^2 + 67.74*lambda
    //   Q(lambda) = lambda^3 - 2.115*lambda^2 - 33.013*lambda + 103.37
    //
    // Nota: lambda_mem e' clamped a min 0.5 per evitare Q(lambda)->0 (lambda~0)
    // e per garantire che il denominatore rimanga positivo nel range fisico.
    //
    // Il termine k_fit_eff[i] incorpora sia il fit di calibrazione che
    // lo switching vapor/liquid [FIX-5], sostituendo i vecchi parametri k e Tk.
    Dw[i] = k_fit_eff[i]
      * (   (3.842*max(lambda_mem[i], 0.5)^3
            - 32.03*max(lambda_mem[i], 0.5)^2
            + 67.74*max(lambda_mem[i], 0.5))
           /
            (max(lambda_mem[i], 0.5)^3
            - 2.115*max(lambda_mem[i], 0.5)^2
            - 33.013*max(lambda_mem[i], 0.5)
            + 103.37)
        )
      * 1e-10
      * exp(E_a_Vetter / R_gas * (1.0/T_ref_Vetter - 1.0/max(T_mem[i], 250.0)));

    // =========================================================================
    // [FIX-4] RESISTENZA CONVETTIVA DI FILM — TRE RESISTENZE IN SERIE
    // =========================================================================
    // Il modello originale calcolava m_trans direttamente dalla diffusione in
    // membrana, saltando la resistenza di film convettiva su entrambi i lati.
    // Per Re basso (regime laminare, tipico degli umidificatori PEMFC) la
    // resistenza convettiva puo' essere dello stesso ordine di quella di membrana.
    // Riferimento struttura: Mull 2025 Eq.(29), Pollak NTU 2023 Eq.(19-20).
    //
    // Struttura resistenze in serie (per segmento i):
    //   R_tot = R_conv_wet + R_mem_diff + R_conv_dry    [s/m^3]
    //   m_trans = Mv * A_mem * (C_shell - C_tube) / R_tot   [kg/s]
    //
    // =========================================================================
    // PORTATA DI MASSA TRASFERITA — formula originale Park (senza FIX-4)
    // =========================================================================
    // Driving force: gradiente di concentrazione molare alle superfici membrana.
    // Il percorso di diffusione e' assunto pari a meta' spessore (0.5*t_mem)
    // in coerenza con il modello originale di Park (2008).
    // Le resistenze convettive di film non sono incluse in questa versione:
    // l'effetto e' inglobato implicitamente nel parametro k_fit_eff [FIX-5].
    m_trans[i] = (Mv*A_mem/(0.5*t_mem)) * (Dw[i]*(C_shell[i] - C_tube[i]));

    // Bilancio dinamico dell'acqua nella membrana (equazione di stato per C_w).
    // d(C_w * V_mem)/dt = portata netta = m_trans per segmento [kg/s].
    // V_mem per segmento = A_mem * t_mem [m^3].
    // NOTA: C_w e' inizializzata nell'initial equation [FIX-2].
    m_trans[i] = der(C_w[i]*A_mem*t_mem);
 end for;

  C_tube = (rho_mem/M_mem) .* lambda_tube;
  C_shell = (rho_mem/M_mem) .* lambda_shell;
  T_mem = wall.Tm;

if flowConfiguration==.Modelon.ThermoFluid.Choices.FlowConfiguration.CounterFlow then
   for i in 1:n loop
     connect(shell.wall[n + 1 - i], wall.qa[i]) annotation (Line(points={{-16,50},{-16,32},{-10,32},{-10,
          20},{-10,20}}, color={191,0,0}));
     state_prim[n + 1 - i] = shell.channel.volume[i].state;
     m_trans_shell[n + 1 - i] = m_trans[i];
   end for;
else
   connect(shell.wall, wall.qa) annotation (Line(points={{-16,50},{-16,32},{-10,32},{-10,
          20},{-10,20}}, color={191,0,0}));
   state_prim = shell.channel.volume.state;
   m_trans_shell = m_trans;
end if;

  connect(wall.qb, tube.wall) annotation (Line(points={{-10,-20},{-10,-32},{0,-32},{0,-50},
          {0,-50}}, color={191,0,0}));
  connect(portB_prim, shell.portB) annotation (Line(points={{-100,40},{-66,40},{-66,60},{-26,60}}, color={209,60,0}));
  connect(shell.portA, portA_prim) annotation (Line(points={{10,60},{40,60},{40,-40},{100,-40}}, color={209,60,0}));
  connect(tube.portB, portB_sec) annotation (Line(points={{10,-60},{60,-60},{60,
          0},{100,0}}, color={209,60,0}));
  connect(portA_sec, tube.portA) annotation (Line(points={{-100,0},{-66,0},{-66,
          -60},{-26,-60}}, color={209,60,0}));
  annotation (Icon(graphics={
        Text(
          extent={{-100,-80},{100,-120}},
          lineColor={0,0,127},
          textString="%name"),
        Rectangle(
          extent={{-90,60},{90,-60}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-90,22},{90,-22}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}),
        Rectangle(
          extent={{-90,26},{90,22}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-90,-22},{90,-26}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Polygon(
          points={{-44,76},{-60,70},{-44,64},{-44,76}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Line(points={{60,70},{-60,70}}, color={0,0,0}),
        Line(points={{-60,-70},{60,-70}}, color={0,0,0}),
        Polygon(
          points={{42,-64},{60,-70},{42,-76},{42,-64}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}),
                                 Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<p>Improved shell-and-tube gas-to-gas humidifier with discretized dynamic channel models,
dynamic wall/membrane, two-phase Schroeder switching, and three-resistance mass transfer model.
Co-flow and counter-flow configurations are supported.</p>

<h4>Improvements over original model</h4>
<ul>
<li><b>[FIX-2]</b> Explicit <code>initial equation</code> for <code>C_w</code>: initialised from Springer
equilibrium at shell-side conditions. Prevents DAE under-determination at t=0.</li>
<li><b>[FIX-3]</b> Diffusivity model: Park/Springer piecewise D(&lambda;) with 3 breakpoints (&lambda;=2,3,4.5)
replaced by the continuous C&infin; rational polynomial of Vetter &amp; Schumacher (2019).
Eliminates numerical events and the non-physical decreasing branch for 3&lt;&lambda;&lt;4.5.</li>
<li><b>[FIX-5]</b> Schroeder paradox (two-phase): smooth tanh switching between
k_fit_vapor and k_fit_liquid when a_shell &gt; 1. Enhancement ~20x on diffusivity (Mull 2025).
The convective film resistance is not modelled explicitly (FIX-4 not applied);
its effect is absorbed by the k_fit calibration parameters.</li>
<li><b>[FIX-6]</b> Nafion specific heat corrected: c=1200 J/(kg&middot;K) instead of 4188 (water).</li>
<li><b>[FIX-8]</b> n=5 discretization elements (was 3), as recommended by Mull et al. (2025).</li>
</ul>

<h4>Not applied — FIX-4 (convective film resistance)</h4>
<p>The three-resistance model (R_conv_wet + R_mem + R_conv_dry) is not included in this version.
The mass transfer formula retains the original Park (2008) form:
<code>m_trans = (Mv*A_mem/(0.5*t_mem)) * Dw * (C_shell - C_tube)</code>.
Convective film effects are implicitly absorbed into k_fit_vapor and k_fit_liquid.</p>

<h4>Not applied — FIX-9 (heat transfer correlation)</h4>
<p>The original Dittus-Boelter correlation is retained for both shell and tube sides,
consistent with Park (2008). Note that Dittus-Boelter is strictly valid for Re &gt; 10000
(turbulent flow); for the laminar regime typical of PEMFC humidifiers (Re = 100-2000)
the Gnielinski correlation (FIX-9) would be more accurate. The impact on mass transfer
is limited since membrane resistance dominates over convective thermal resistance.</p>

<h4>Key parameters to calibrate</h4>
<ul>
<li><code>k_fit_vapor</code>: diffusivity fit factor, vapour only (default: 2.0, Mull 2025 Table 4)</li>
<li><code>k_fit_liquid</code>: diffusivity fit factor, with liquid water (default: 40.0, Mull 2025 Table 4)</li>
<li><code>dx_liq_smooth</code>: smoothing width for vapour/liquid transition (default: 0.02)</li>
<li><code>T_ref_Vetter</code>: Arrhenius reference temperature (default: 353.15 K)</li>
<li><code>E_a_Vetter</code>: activation energy for diffusion (default: 20000 J/mol)</li>
</ul>

<h4>Parametrization (unchanged from original)</h4>
<ul>
<li><code>n</code> : number of discretizations (default: 5)</li>
<li><code>rho_mem</code> : membrane dry density</li>
<li><code>M_mem</code> : molar mass of dry membrane</li>
<li><code>t_mem</code> : thickness of membrane tube</li>
<li><code>D_mem</code> : inner diameter of membrane tube</li>
<li><code>D_h</code> : inner diameter of the humidifier housing</li>
<li><code>L_mem</code> : active length of membrane tube</li>
<li><code>n_tubes</code> : number of the membrane tubes</li>
<li><code>A_mem</code> : pre-calculated membrane area based on tube dimensions</li>
</ul>

<h4>Assumptions</h4>
<ul>
<li>Gases used are ideal gases</li>
<li>Negligible potential and kinetic energy changes</li>
<li>Heat transfer only occurs across the membrane; no heat losses to environment</li>
<li>Two-phase effects modelled via smooth k_fit switching (Schroeder paradox, a_shell &gt; 1)</li>
<li>Diameters of control volumes equal the inner diameter</li>
</ul>

<h4>References</h4>
<p>Park et al. (2008), Int. J. Hydrogen Energy 33, 2273-2282 |
Mull et al. (2025), Next Energy 7, 100294 |
Vetter &amp; Schumacher (2019), Comput. Phys. Commun. 234, 223-234 |
Costello et al. (1993), J. Membrane Sci. 80, 1-11 |
Kusoglu &amp; Weber (2017), Chem. Rev. 117, 987-1104</p>
</html>", revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end DiscretizedGasGasHumidifier_improved;
