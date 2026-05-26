model HollowFiberHumidifier_Enhanced
  "Shell-and-tube hollow fiber humidifier - physically enhanced model.
   Corrects diffusion path, enthalpy of transported water, and adds
   convective mass transfer resistance on both sides.
   Reference: Park et al. (2008), Huizing et al. (2012)"

  extends .FuelCell.HeatExchangers.Interfaces.DistributedGasGasHumidifier(
    wallThickness = t_mem,
    // Abilita heat transfer: accoppiato al mass transfer tramite T_mem
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
      m_flow     = m_flow_prim,
      m_flow_sec = m_flow_sec,
      Q_flow     = sum(shell.wall.Q_flow),
      T_in    = noEvent(if m_flow_prim >= 0 then shell.summary.T_in  else shell.summary.T_out),
      T_out   = noEvent(if m_flow_prim >= 0 then shell.summary.T_out else shell.summary.T_in),
      h_in    = noEvent(if m_flow_prim >= 0 then shell.summary.h_in  else shell.summary.h_out),
      h_out   = noEvent(if m_flow_prim >= 0 then shell.summary.h_out else shell.summary.h_in),
      dp      = shell.summary.dp,
      T_sec_in  = noEvent(if m_flow_sec >= 0 then tube.summary.T_in  else tube.summary.T_out),
      T_sec_out = noEvent(if m_flow_sec >= 0 then tube.summary.T_out else tube.summary.T_in),
      h_sec_in  = noEvent(if m_flow_sec >= 0 then tube.summary.h_in  else tube.summary.h_out),
      h_sec_out = noEvent(if m_flow_sec >= 0 then tube.summary.h_out else tube.summary.h_in),
      dp_sec    = tube.summary.dp,
      p_in      = portA_prim.p,
      p_out     = portB_prim.p,
      p_sec_in  = portA_sec.p,
      p_sec_out = portB_sec.p,
      Q_flow_sec = sum(tube.wall.Q_flow)));

  // =========================================================================
  // PARAMETRI GEOMETRICI — Hollow Fiber Bundle
  // =========================================================================
  parameter Integer    n_fibers = 500
    "Numero di fibre cavi nel bundle"
    annotation(Dialog(tab="Geometry", group="Hollow Fiber"));
  parameter .Modelica.Units.SI.Length d_i = 0.30e-3
    "Diametro interno della fibra [m]"
    annotation(Dialog(tab="Geometry", group="Hollow Fiber"));
  parameter .Modelica.Units.SI.Length t_mem = 0.05e-3
    "Spessore della parete della fibra (membrana) [m]"
    annotation(Dialog(tab="Geometry", group="Hollow Fiber"));
  parameter .Modelica.Units.SI.Length L_mem = 0.20
    "Lunghezza attiva delle fibre [m]"
    annotation(Dialog(tab="Geometry", group="Hollow Fiber"));
  parameter .Modelica.Units.SI.Length D_shell = 0.10
    "Diametro interno del contenitore (shell) [m]"
    annotation(Dialog(tab="Geometry", group="Shell"));

  // Parametri derivati — geometria
  final parameter .Modelica.Units.SI.Length d_o = d_i + 2*t_mem
    "Diametro esterno della fibra [m]";

  // Diametro idraulico lato shell (bundle di fibre nel contenitore):
  // D_h,shell = (D_shell^2 - N_f*d_o^2) / (D_shell + N_f*d_o)
  // NOTA: questa è l'approssimazione di Happel per bundle circolari
  final parameter .Modelica.Units.SI.Length D_h_shell =
    (D_shell^2 - n_fibers*d_o^2) / (D_shell + n_fibers*d_o)
    "Diametro idraulico lato shell [m]";

  // Area di trasferimento basata sul diametro MEDIO della fibra (log-mean area)
  // Più accurata di usare solo d_i o d_o
  final parameter .Modelica.Units.SI.Length d_lm =
    (d_o - d_i) / Modelica.Math.log(d_o/d_i)
    "Diametro medio logaritmico fibra [m]";
  final parameter .Modelica.Units.SI.Area A_mem_tot =
    Modelica.Constants.pi * d_lm * L_mem * n_fibers
    "Area totale membrana (log-mean) [m²]";
  final parameter .Modelica.Units.SI.Area A_mem = A_mem_tot / n
    "Area membrana per segmento [m²]";

  // =========================================================================
  // PARAMETRI MEMBRANA (Nafion tubolare / PFSA)
  // =========================================================================
  parameter Real rho_mem = 1970
    "Densità secca della membrana [kg/m³] — Nafion: ~1970"
    annotation(Dialog(tab="Membrane"));
  parameter Real M_mem = 1.1
    "Massa equivalente della membrana [kg/mol EW] — Nafion: ~1.1"
    annotation(Dialog(tab="Membrane"));
  parameter Real k = 2416
    "Costante di Arrhenius per Dw [K] — tipico Nafion: 2416"
    annotation(Dialog(tab="Membrane"));
  parameter Real Tk = 303.15
    "Temperatura di riferimento Arrhenius [K]"
    annotation(Dialog(tab="Membrane"));

  // =========================================================================
  // PARAMETRI MASS TRANSFER CONVETTIVO
  // =========================================================================
  parameter Boolean useConvectiveMassTransfer = true
    "Considera resistenza convettiva al mass transfer su entrambi i lati"
    annotation(Dialog(tab="Correlations", group="Mass Transfer"));
  parameter Real Sc_shell = 0.6
    "Numero di Schmidt lato shell (aria umida ~273K, 3bar)"
    annotation(Dialog(tab="Correlations", group="Mass Transfer",
               enable=useConvectiveMassTransfer));
  parameter Real Sc_tube = 0.6
    "Numero di Schmidt lato tube (aria secca)"
    annotation(Dialog(tab="Correlations", group="Mass Transfer",
               enable=useConvectiveMassTransfer));
  // Diffusività vapore acqua in aria [m²/s] — Chapman-Enskog a T_ref
  parameter .Modelica.Units.SI.DiffusionCoefficient D_va_ref = 2.56e-5
    "Diffusività H2O in aria a 298 K, 1 atm [m²/s]"
    annotation(Dialog(tab="Correlations", group="Mass Transfer"));

  // =========================================================================
  // VARIABILI LOCALI
  // =========================================================================
  .Modelica.Units.SI.MassFlowRate m_flow_prim
    "Portata massica lato primario (shell)";
  .Modelica.Units.SI.MassFlowRate m_flow_sec
    "Portata massica lato secondario (tube)";
  .Modelica.Units.SI.HeatFlowRate Q_tot
    "Flusso termico totale uscente dal lato primario";

  .Modelica.Units.SI.MassFlowRate[n] m_trans
    "Portata massica vapore attraverso la membrana per segmento [kg/s]";
  .Modelica.Units.SI.MassFlowRate m_water = sum(m_trans)
    "Portata massica totale di vapore trasferito";

  // Attività acqua
  .Modelica.Units.SI.DimensionlessRatio a_shell[n] "Attività acqua lato shell";
  .Modelica.Units.SI.DimensionlessRatio a_tube[n]  "Attività acqua lato tube";
  .Modelica.Units.SI.DimensionlessRatio a_mem[n]   "Attività acqua in membrana";

  // Contenuto d'acqua nella membrana (lambda = mol H2O / mol SO3-)
  Real lambda_shell[n] "Contenuto d'acqua interfaccia shell-membrana";
  Real lambda_tube[n]  "Contenuto d'acqua interfaccia tube-membrana";
  Real lambda_mem[n]   "Contenuto d'acqua medio in membrana";

  // Concentrazioni molari di acqua nella membrana [mol/m³]
  .Modelica.Units.SI.MolarDensity C_shell[n]
    "Concentrazione H2O lato shell della membrana";
  .Modelica.Units.SI.MolarDensity C_tube[n]
    "Concentrazione H2O lato tube della membrana";

  // Coefficienti di diffusione
  Real D_emp[n]
    "Coefficiente empirico D(lambda_mem) — Motupally et al.";
  .Modelica.Units.SI.DiffusionCoefficient Dw[n]
    "Coefficiente di diffusione effettivo con correzione Arrhenius [m²/s]";

  // FIX: resistenza convettiva al mass transfer
  .Modelica.Units.SI.CoefficientOfHeatTransfer[n] km_shell
    "Coeff. convettivo mass transfer lato shell [m/s] — calcolato da Re, Sc";
  .Modelica.Units.SI.CoefficientOfHeatTransfer[n] km_tube
    "Coeff. convettivo mass transfer lato tube [m/s]";
  Real[n] R_mem   "Resistenza membrana [s/m]";
  Real[n] R_shell "Resistenza convettiva lato shell [s/m]";
  Real[n] R_tube  "Resistenza convettiva lato tube [s/m]";
  Real[n] R_tot   "Resistenza totale al trasferimento di massa [s/m]";

  // Pressioni parziali e di saturazione
  .Modelica.Units.SI.Pressure pV_shell[n]   "Pressione parziale vapore shell";
  .Modelica.Units.SI.Pressure pV_tube[n]    "Pressione parziale vapore tube";
  .Modelica.Units.SI.Pressure psat_shell[n] "Pressione di saturazione shell";
  .Modelica.Units.SI.Pressure psat_tube[n]  "Pressione di saturazione tube";

  // Frazioni molare
  .Modelon.Media.Units.MoleFraction[n, nX_shell] y_prim;
  .Modelon.Media.Units.MoleFraction[n, nX_tube]  y_sec;

  // Entalpia trasportata (CORRETTA: solo h_vapore, non h_liq + h_vap)
  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim
    "Entalpia specifica del vapore trasportato — lato shell";
  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec
    "Entalpia specifica del vapore trasportato — lato tube";

  // Temperatura membrana (dalla DynamicWall)
  .Modelica.Units.SI.Temperature T_mem[n];

  // Stati termodinamici
  PrimaryMedium.ThermodynamicState   state_prim[n];
  SecondaryMedium.ThermodynamicState state_sec[n];

  // Indici sostanze
  constant Integer H2O_prim  = PrimaryMedium.substanceIndex("H2O");
  constant Integer H2O_sec   = SecondaryMedium.substanceIndex("H2O");
  constant Integer nX_shell  = PrimaryMedium.nS;
  constant Integer nX_tube   = SecondaryMedium.nS;
  constant .Modelica.Units.SI.MolarMass Mv = 18.01528e-3 "Massa molare H2O";

  package WaterMedium = .Modelon.Media.PreDefined.TwoPhase.WaterIF97;

  // Variabile ausiliaria per counter-flow
  .Modelica.Units.SI.MassFlowRate[n] m_trans_shell;

  // =========================================================================
  // CANALI E PARETE
  // =========================================================================
  .FuelCell.Pipes.FlowChannel shell(
    redeclare package Medium = PrimaryMedium,
    // CORRETTO: area sezione trasversale del lato shell (contenitore - fibre)
    A      = fill(Modelica.Constants.pi/4 * (D_shell^2 - n_fibers*d_o^2), n),
    A_heat = fill(A_mem, n),
    useHeatTransfer = true,   // ABILITATO: necessario per T_mem corretta
    n      = n,
    L      = fill(L_mem/n, n),
    Dhyd   = fill(D_h_shell, n),  // CORRETTO: D_h del bundle
    H_flow   = h_water_prim .* m_trans_shell,
    mX_flow  = transpose({if i == H2O_prim then -m_trans_shell
                          else zeros(n) for i in 1:nX_shell}),
    // ... altri parametri di inizializzazione come nell'originale
    positiveFlow = positiveFlow_prim,
    initOpt     = initOpt,
    p_start_in  = init_prim.p_in,
    p_start_out = init_prim.p_out,
    h_start_in  = init_prim.h_in_used,
    h_start_out = init_prim.h_out_used,
    X_start     = init_prim.X,
    m_flow_start = init_prim.m_flow,
    redeclare model Friction      = Friction_prim,
    redeclare model HeatTransfer  = HeatTransfer_prim,
    CF_PressureLoss  = CF_PrimarySidePressureLoss,
    CF_HeatTransfer  = CF_PrimarySideHeatTransfer)
    "Lato shell — gas umido (es. scarico catodico)";

  .FuelCell.Pipes.FlowChannel tube(
    redeclare package Medium = SecondaryMedium,
    n_channels = fill(n_fibers, n),
    // CORRETTO: area sezione interna totale delle fibre
    A      = fill(n_fibers * Modelica.Constants.pi/4 * d_i^2, n),
    A_heat = fill(A_mem, n),
    useHeatTransfer = true,
    n      = n,
    L      = fill(L_mem/n, n),
    Dhyd   = fill(d_i, n),    // Diametro idraulico fibra = d_i (tubolare)
    H_flow   = h_water_sec .* m_trans,
    mX_flow  = transpose({if i == H2O_sec then m_trans
                          else zeros(n) for i in 1:nX_tube}),
    positiveFlow = positiveFlow_sec,
    initOpt     = initOpt_sec,
    p_start_in  = init_sec.p_in_used,
    p_start_out = init_sec.p_out_used,
    h_start_in  = init_sec.h_in_used,
    h_start_out = init_sec.h_out_used,
    X_start     = init_sec.X,
    m_flow_start = init_sec.m_flow,
    redeclare model Friction      = Friction_sec,
    redeclare model HeatTransfer  = HeatTransfer_sec,
    CF_PressureLoss  = CF_SecondarySidePressureLoss,
    CF_HeatTransfer  = CF_SecondarySideHeatTransfer)
    "Lato tube — gas secco (aria alimentazione catodo)";

  // Proprietà materiale membrana CORRETTE per Nafion/PFSA
  replaceable record WallMaterial =
    .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
    (
      c      = 1150,   // CORRETTO: Cp Nafion secco ~1100-1250 J/(kg·K)
      rho    = 1970,   // kg/m³ — Nafion secco
      lambda = 0.13)   // W/(m·K) — Nafion umido
    constrainedby
    .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
    "Proprietà termofisiche membrana (default: Nafion)"
    annotation(choicesAllMatching, Dialog(tab="Wall", group="Wall geometry"));

  WallMaterial wallMaterial;

  .Modelon.ThermoFluid.Solids.DynamicWall wall(
    props(each Cp = wallMaterial.c,
          Rw = wallThickness ./ wallMaterial.lambda ./ wallCrossSecArea),
    n = n,
    steadyStateInit = WallSteadyStateInit,
    n_channels = n_channels_wall,
    m  = m,
    T0 = T0_wall,
    massLessWall = massLessWall,
    includeThermalResistance = includeThermalResistance,
    surfaceState = false)
    "Parete dinamica = membrana hollow fiber";

initial equation
  assert(H2O_prim > 0,
    "HollowFiberHumidifier: PrimaryMedium non contiene H2O",
    level=AssertionLevel.error);
  assert(H2O_sec > 0,
    "HollowFiberHumidifier: SecondaryMedium non contiene H2O",
    level=AssertionLevel.error);
  assert(d_i > 0 and t_mem > 0,
    "HollowFiberHumidifier: geometria fibra non valida",
    level=AssertionLevel.error);
  assert(n_fibers * Modelica.Constants.pi/4 * d_o^2
         < Modelica.Constants.pi/4 * D_shell^2,
    "HollowFiberHumidifier: le fibre non entrano nel contenitore! Ridurre n_fibers o d_o",
    level=AssertionLevel.error);

equation
  m_flow_prim = shell.summary.m_flow;
  m_flow_sec  = tube.summary.m_flow;
  Q_tot       = -1*sum(shell.wall.Q_flow);
  state_sec   = tube.channel.volume.state;
  T_mem       = wall.Tm;

  for i in 1:n loop

    // ------------------------------------------------------------------
    // 1. STATI TERMODINAMICI e FRAZIONI MOLARE
    // ------------------------------------------------------------------
    y_prim[i,:] = PrimaryMedium.massToMoleFractions(
                    state_prim[i].X, PrimaryMedium.MMX);
    y_sec[i,:]  = SecondaryMedium.massToMoleFractions(
                    state_sec[i].X, SecondaryMedium.MMX);

    // Pressioni parziali vapore
    pV_shell[i] = y_prim[i, H2O_prim] * state_prim[i].p;
    pV_tube[i]  = y_sec[i,  H2O_sec]  * state_sec[i].p;

    // Pressione di saturazione (guard: min con 0.999*p evita a > 1 per discontinuità)
    psat_shell[i] = min(
      WaterMedium.saturationPressure_TX(max(shell.channel.volume[i].T, 273.16)),
      0.999 * state_prim[i].p);
    psat_tube[i]  = min(
      WaterMedium.saturationPressure_TX(max(tube.channel.volume[i].T, 273.16)),
      0.999 * state_sec[i].p);

    // Attività acqua (clamp a [0, 3] per stabilità numerica)
    a_shell[i] = smooth(1, max(0, min(3, pV_shell[i] / psat_shell[i])));
    a_tube[i]  = smooth(1, max(0, min(3, pV_tube[i]  / psat_tube[i])));

    // Attività membrana: media aritmetica (Springer 1991)
    a_mem[i] = (a_shell[i] + a_tube[i]) / 2;

    // ------------------------------------------------------------------
    // 2. CORRELAZIONE λ(a) — Springer et al. (1991), Nafion 117
    //    NOTA: sostituibile con correlazioni per altri ionomeri
    // ------------------------------------------------------------------
    lambda_shell[i] = if a_shell[i] <= 1 then
        0.043 + 17.81*a_shell[i] - 39.85*a_shell[i]^2 + 36.0*a_shell[i]^3
      elseif a_shell[i] <= 3 then
        14.0 + 1.4*(a_shell[i] - 1)
      else 16.8;

    lambda_tube[i] = if a_tube[i] <= 1 then
        0.043 + 17.81*a_tube[i] - 39.85*a_tube[i]^2 + 36.0*a_tube[i]^3
      elseif a_tube[i] <= 3 then
        14.0 + 1.4*(a_tube[i] - 1)
      else 16.8;

    lambda_mem[i] = if a_mem[i] <= 1 then
        0.043 + 17.81*a_mem[i] - 39.85*a_mem[i]^2 + 36.0*a_mem[i]^3
      elseif a_mem[i] <= 3 then
        14.0 + 1.4*(a_mem[i] - 1)
      else 16.8;

    // ------------------------------------------------------------------
    // 3. CONCENTRAZIONI MOLARI IN MEMBRANA [mol/m³]
    //    C = (rho_mem / M_mem) * lambda   [EW normalizzato]
    // ------------------------------------------------------------------
    C_shell[i] = (rho_mem / M_mem) * lambda_shell[i];
    C_tube[i]  = (rho_mem / M_mem) * lambda_tube[i];

    // ------------------------------------------------------------------
    // 4. COEFFICIENTE DI DIFFUSIONE — Motupally et al. (2000) + Arrhenius
    // ------------------------------------------------------------------
    D_emp[i] = if lambda_mem[i] < 2 then
        1.0e-6
      elseif lambda_mem[i] < 3 then
        1.0e-6 * (1.0 + 2.0*(lambda_mem[i] - 2.0))
      elseif lambda_mem[i] < 4.5 then
        1.0e-6 * (3.0 - (5.0/3.0)*(lambda_mem[i] - 3.0))
      else
        1.25e-6;

    // Correzione Arrhenius: Dw(T) = D_emp * exp(k*(1/Tk - 1/T_mem))
    Dw[i] = D_emp[i] * Modelica.Math.exp(k * (1.0/Tk - 1.0/T_mem[i]));

    // ------------------------------------------------------------------
    // 5. RESISTENZE AL TRASFERIMENTO DI MASSA
    //
    //    Resistenza membrana:   R_mem   = t_mem / Dw      [s/m]
    //    Resistenza convettiva: R_shell = 1 / km_shell    [s/m]
    //                           R_tube  = 1 / km_tube
    //
    //    Per gas in moto laminare in condotti (Re < 2300):
    //    Sh = 3.66 (Nusselt termico analogo, condizione T uniforme)
    //    km = Sh * D_va / D_h
    //
    //    D_va(T,p) = D_va_ref * (T/298)^1.75 / (p/101325)  [Chapman-Enskog]
    // ------------------------------------------------------------------
    R_mem[i] = t_mem / Dw[i];  // CORRETTO: path length = t_mem (non 0.5*t_mem)

    if useConvectiveMassTransfer then
      // Diffusività H2O in aria corretta per T e p (Chapman-Enskog)
      // Shell side
      km_shell[i] = 3.66 *
        (D_va_ref * (shell.channel.volume[i].T / 298.0)^1.75
          / (state_prim[i].p / 101325.0)) / D_h_shell;
      // Tube side
      km_tube[i]  = 3.66 *
        (D_va_ref * (tube.channel.volume[i].T / 298.0)^1.75
          / (state_sec[i].p / 101325.0)) / d_i;
      R_shell[i] = 1.0 / km_shell[i];
      R_tube[i]  = 1.0 / km_tube[i];
    else
      km_shell[i] = 1e10;  // Resistenza convettiva nulla → R ≈ 0
      km_tube[i]  = 1e10;
      R_shell[i]  = 0.0;
      R_tube[i]   = 0.0;
    end if;

    R_tot[i] = R_shell[i] + R_mem[i] + R_tube[i];

    // ------------------------------------------------------------------
    // 6. FLUSSO DI MASSA ATTRAVERSO LA MEMBRANA [kg/s]
    //
    //    m_trans = Mv * A_mem * ΔC / R_tot
    //    FIX: path length CORRETTA = t_mem (non 0.5*t_mem come nell'originale)
    //    FIX: rimossa la variabile di stato C_w disconnessa
    // ------------------------------------------------------------------
    m_trans[i] = Mv * A_mem * (C_shell[i] - C_tube[i]) / R_tot[i];

    // ------------------------------------------------------------------
    // 7. ENTALPIA DEL VAPORE TRASPORTATO
    //    FIX: h_vapore = enthalpyOfCondensingGas(T) [J/kg]
    //    L'originale sommava erroneamente h_gas + h_liq
    // ------------------------------------------------------------------
    h_water_prim[i] = PrimaryMedium.enthalpyOfCondensingGas(state_prim[i].T);
    h_water_sec[i]  = SecondaryMedium.enthalpyOfCondensingGas(state_sec[i].T);

  end for;

  // ------------------------------------------------------------------
  // 8. CONFIGURAZIONE FLUSSO (CO-FLOW / COUNTER-FLOW)
  // ------------------------------------------------------------------
  if flowConfiguration == .Modelon.ThermoFluid.Choices.FlowConfiguration.CounterFlow then
    for i in 1:n loop
      connect(shell.wall[n + 1 - i], wall.qa[i]);
      state_prim[n + 1 - i]  = shell.channel.volume[i].state;
      m_trans_shell[n + 1 - i] = m_trans[i];
    end for;
  else
    connect(shell.wall, wall.qa);
    state_prim    = shell.channel.volume.state;
    m_trans_shell = m_trans;
  end if;

  connect(wall.qb,     tube.wall);
  connect(portB_prim,  shell.portB);
  connect(shell.portA, portA_prim);
  connect(tube.portB,  portB_sec);
  connect(portA_sec,   tube.portA);

  annotation (Documentation(info="<html>
<h4>Umidificatore Shell-and-Tube a Fibre Cave — Modello Fisico Avanzato</h4>
<p>Miglioramenti rispetto a DiscretizedGasGasHumidifier:</p>
<ul>
  <li><b>Geometria hollow fiber</b>: D_h_shell calcolato dal bundle (Happel),
      area di trasferimento basata sul diametro medio logaritmico.</li>
  <li><b>Path di diffusione corretto</b>: usa t_mem (non 0.5*t_mem).</li>
  <li><b>Resistenza convettiva</b>: Sh=3.66 (laminare) su entrambi i lati,
      con D_va(T,p) Chapman-Enskog.</li>
  <li><b>Entalpia vapore corretta</b>: enthalpyOfCondensingGas(T).</li>
  <li><b>Rimossa C_w disconnessa</b>: eliminato stato ODE non fisico.</li>
  <li><b>Cp Nafion corretto</b>: 1150 J/(kg·K) invece di 4188.</li>
  <li><b>useHeatTransfer=true</b>: abilitato di default.</li>
</ul>
<h4>Riferimenti</h4>
<ul>
  <li>Park et al., Int. J. Hydrogen Energy 33 (2008) 2273-2282</li>
  <li>Springer et al., J. Electrochem. Soc. 138 (1991) 2334</li>
  <li>Motupally et al., J. Electrochem. Soc. 147 (2000) 3171</li>
  <li>Huizing et al., J. Membr. Sci. 369 (2011) 127-135</li>
</ul>
</html>"));

end HollowFiberHumidifier_Enhanced;