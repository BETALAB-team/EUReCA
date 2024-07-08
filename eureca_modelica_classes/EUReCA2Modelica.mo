within ;
package EUReCA2Modelica

  model Zone_EUReCA "Thermal zone based on 5R1C network"
    parameter Modelica.Units.SI.ThermalConductance Htr_w "Heat transfer through glazed elements [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_is "Coupling conductance betwee air and surface nodes [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_ms "Heat transfer through opaque elements [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_em "Heat transfer through opaque elements [W/K]";
    parameter Modelica.Units.SI.HeatCapacity Cm "Zone thermal capacity [J/K]";
    parameter Modelica.Units.SI.HeatFlowRate cooling_design_flow_rate "Cooling maximum capacity [W]";
    parameter Modelica.Units.SI.HeatFlowRate heating_design_flow_rate "Heating maximum capacity [W]";
    parameter String path_to_sched_file "Path to schedule file";
    parameter String path_to_mos_weather_file "Path to weather file";

    Modelica.Blocks.Interfaces.RealOutput TAir(final unit = "K", displayUnit = "degC") "Room air temperature" annotation (
      Placement(transformation(origin={-82,18},  extent = {{140, 70}, {160, 90}}), transformation(extent = {{140, 70}, {160, 90}})));
    Modelica.Blocks.Interfaces.RealOutput TSur(final unit = "K", displayUnit = "degC")
      "Average inside surface temperature"                                                                                  annotation (
      Placement(transformation(origin={-82,36},  extent = {{140, -10}, {160, 10}}), transformation(extent = {{140, -10}, {160, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus "Weather data bus" annotation (
      Placement(transformation(origin={-300,-34},  extent = {{68, 78}, {132, 142}}), transformation(origin = {-94, 20}, extent = {{68, 78}, {132, 142}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heaPorAir
      "Heat port to air node"                                                             annotation (
      Placement(transformation(extent={{-54,66},{-34,86}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heaPorSur "Heat port to surface temperatures" annotation (
      Placement(transformation(extent={{-54,-14},{-34,6}})));

     Modelica.Thermal.HeatTransfer.Components.ThermalConductor HTra(G=Htr_em) "Heat transfer through opaque elements" annotation (
      Placement(transformation(origin={-75,-85},  extent = {{-11, -11}, {11, 11}})));
     Modelica.Thermal.HeatTransfer.Components.ThermalConductor HWin(G=Htr_w) "Heat transfer through glazed elements" annotation (
      Placement(transformation(origin={-86,-4},  extent = {{0, -10}, {20, 10}})));
     Modelica.Thermal.HeatTransfer.Components.ThermalConductor HThe(G=Htr_is) "Coupling conductance betwee air and surface nodes" annotation (
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin={-44,36})));
     Modelica.Thermal.HeatTransfer.Components.ThermalConductor HMas(G=Htr_ms) "Coupling conductance between surface and mass nodes" annotation (
      Placement(transformation(origin={-44,-50},   extent = {{-10, -10}, {10, 10}}, rotation = 90)));
     Modelica.Thermal.HeatTransfer.Components.HeatCapacitor capMas(C = Cm, T(displayUnit = "degC", fixed = true, start = 293.15)) "Zone thermal capacity" annotation (
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin={-44,-108})));
    Buildings.HeatTransfer.Sources.PrescribedTemperature TExt "External air temperature" annotation (
      Placement(transformation(extent={{-148,-28},{-128,-8}})));
    Buildings.HeatTransfer.Sources.PrescribedTemperature TVen "Supply air temperature" annotation (
      Placement(transformation(extent={{-150,66},{-130,86}})));
    Buildings.Controls.Continuous.LimPID conCooPID(
      Ti=300,
      k=0.1,
      reverseActing=false,
      strict=true)                                                                                          annotation (
      Placement(transformation(origin={138,-38},    extent = {{-6, 6}, {6, -6}})));
    Buildings.Controls.Continuous.LimPID conHeaPID(
      Ti=300,
      k=0.1,
      reverseActing=true,
      strict=true)                                                                                         annotation (
      Placement(transformation(origin={130,0},      extent = {{-6, 6}, {6, -6}}, rotation = -0)));
    Modelica.Blocks.Math.Gain gaiCoo(k=cooling_design_flow_rate*2)     annotation (
      Placement(transformation(origin={176,-38},   extent = {{-6, -6}, {6, 6}})));
    Modelica.Blocks.Math.Gain gaiHea(k=heating_design_flow_rate*2)     annotation (
      Placement(transformation(origin={174,0},     extent = {{-6, -6}, {6, 6}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow preHeaCoo annotation (
      Placement(transformation(origin={180,-72},    extent = {{58, 54}, {70, 66}})));
    Modelica.Blocks.Math.Sum sumHeaCoo(nin=2)   annotation (
      Placement(transformation(origin={180,-72},    extent = {{44, 56}, {52, 64}})));
    Modelica.Blocks.Routing.Multiplex2 multiplex2 annotation (
      Placement(transformation(origin={180,-72},    extent = {{30, 56}, {38, 64}})));
    Modelica.Blocks.Math.Max p_design_cool(u1=cooling_design_flow_rate)   annotation (
      Placement(transformation(origin={196,-36},   extent = {{-4, -4}, {4, 4}})));
    Modelica.Blocks.Math.Min p_design_heat(u2=heating_design_flow_rate)
      annotation (Placement(transformation(extent={{192,-6},{200,2}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
      tableOnFile=true,
      smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
      tableName="tab1",
      fileName= path_to_sched_file,
      timeScale=1,
      columns={2,3,4,5,6,7,8,9,10,11,12,13,14})                                                                                                                                    annotation (
      Placement(transformation(origin={79,-117},    extent={{-11,-11},{11,11}})));
      // "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/data_Zone_1.txt",

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam=path_to_mos_weather_file)                                  "weather data" annotation (
      Placement(transformation(origin={-162,56},   extent = {{-80, 10}, {-60, 30}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b targetHeatCoolDemand
      annotation (Placement(transformation(extent={{104,114},{124,134}})));
    InfFlowThermalConductor infFlowThermalConductor
      annotation (Placement(transformation(extent={{-86,66},{-66,86}})));
    Modelica.Blocks.Interfaces.RealOutput targetDHWDemand(final unit="W",
        displayUnit="W") annotation (Placement(transformation(origin={0,-90},
            extent={{140,-10},{160,10}}), transformation(extent={{140,-10},{160,10}})));
  protected
    parameter Real ratSur = 4.5 "Ratio between the internal surfaces area and the floor area";
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor senTAir "Air temperature sensor" annotation (
      Placement(transformation(origin={37,97},      extent={{-13,-13},{13,13}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor senTSur "Surface temperature sensor" annotation (
      Placement(transformation(origin={114,125},   extent={{-88,-99},{-66,-77}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heaAir annotation (
      Placement(transformation(origin={0,76},     extent = {{10, -10}, {-10, 10}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heaSur annotation (
      Placement(transformation(                  extent={{6,-14},{-14,6}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heaMas annotation (
      Placement(transformation(origin={-8,-84},    extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  equation
    connect(heaPorSur, HWin.port_b) annotation (
      Line(points={{-44,-4},{-66,-4}},  color = {191, 0, 0}));
    connect(heaPorAir, HThe.port_b) annotation (
      Line(points={{-44,76},{-44,46}},    color = {191, 0, 0}));
    connect(heaPorAir, heaAir.port) annotation (
      Line(points={{-44,76},{-10,76}},    color = {191, 0, 0}));
    connect(heaSur.port, heaPorSur) annotation (
      Line(points={{-14,-4},{-44,-4}},  color = {191, 0, 0}));
    connect(heaPorSur, HMas.port_b) annotation (
      Line(points={{-44,-4},{-44,-40}},   color = {191, 0, 0}));
    connect(weaBus.TDryBul, TExt.T) annotation (
      Line(points={{-200,76},{-200,-18},{-150,-18}},                  color = {255, 204, 51}));
    connect(heaPorSur, HThe.port_a) annotation (
      Line(points={{-44,-4},{-44,26}},   color = {191, 0, 0}));
    connect(heaPorSur, heaPorSur) annotation (
      Line(points={{-44,-4},{-44,-4}},  color = {191, 0, 0}));
    connect(HWin.port_a, TExt.port) annotation (
      Line(points={{-86,-4},{-122,-4},{-122,-18},{-128,-18}},
                                         color = {191, 0, 0}));
    connect(HTra.port_a, TExt.port) annotation (
      Line(points={{-86,-85},{-86,-86},{-122,-86},{-122,-18},{-128,-18}},
                                                                 color = {191, 0, 0}));
    connect(heaMas.port, HTra.port_b) annotation (
      Line(points={{-18,-84},{-45,-84},{-45,-85},{-64,-85}},      color = {191, 0, 0}));
    connect(capMas.port, HTra.port_b) annotation (
      Line(points={{-44,-98},{-44,-85},{-64,-85}},     color = {191, 0, 0}));
    connect(HMas.port_a, HTra.port_b) annotation (
      Line(points={{-44,-60},{-44,-85},{-64,-85}},     color = {191, 0, 0}));
    connect(weaBus.TDryBul, TVen.T) annotation (
      Line(points={{-200,76},{-152,76}},                   color = {255, 204, 51}));
    connect(senTSur.port, heaPorSur) annotation (
      Line(points={{26,37},{-24,37},{-24,-4},{-44,-4}},      color = {191, 0, 0}));
    connect(senTAir.port, heaPorAir) annotation (
      Line(points={{24,97},{-24,97},{-24,76},{-44,76}},          color = {191, 0, 0}));
    connect(senTAir.T, TAir) annotation (
      Line(points={{51.3,97},{51.3,98},{68,98}},
                                             color = {0, 0, 127}));
    connect(senTSur.T, TSur) annotation (
      Line(points={{49.1,37},{54,37},{54,36},{68,36}},
                                           color = {0, 0, 127}));
  connect(conHeaPID.y,gaiHea. u) annotation (
      Line(points={{136.6,0},{166.8,0}},     color = {0, 0, 127}));
  connect(conCooPID.y,gaiCoo. u) annotation (
      Line(points={{144.6,-38},{168.8,-38}},                    color = {0, 0, 127}));
    connect(multiplex2.y,sumHeaCoo. u)
      annotation (Line(points={{218.4,-12},{223.2,-12}},
                                                       color={0,0,127}));
    connect(sumHeaCoo.y,preHeaCoo. Q_flow)
      annotation (Line(points={{232.4,-12},{238,-12}},
                                                     color={0,0,127}));
    connect(weaDat.weaBus, weaBus) annotation (Line(
        points={{-222,76},{-200,76}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(TAir, conHeaPID.u_m) annotation (Line(points={{68,98},{130,98},{130,7.2}},
                    color={0,0,127}));
    connect(TAir, conCooPID.u_m) annotation (Line(points={{68,98},{130,98},{130,14},
            {138,14},{138,-30.8}},     color={0,0,127}));
    connect(combiTimeTable.y[4], heaMas.Q_flow) annotation (Line(points={{91.1,-117},
            {91.1,-118},{96,-118},{96,-84},{2,-84}},           color={0,0,127}));
    connect(combiTimeTable.y[3], heaSur.Q_flow) annotation (Line(points={{91.1,-117},
            {91.1,-118},{96,-118},{96,-4},{6,-4}},
          color={0,0,127}));
    connect(combiTimeTable.y[2], heaAir.Q_flow) annotation (Line(points={{91.1,-117},
            {91.1,-118},{96,-118},{96,-4},{86,-4},{86,76},{10,76}},
                              color={0,0,127}));
    connect(combiTimeTable.y[9], conHeaPID.u_s) annotation (Line(points={{91.1,-117},
            {91.1,-118},{96,-118},{96,0},{122.8,0}},          color={0,0,127}));
    connect(combiTimeTable.y[10], conCooPID.u_s) annotation (Line(points={{91.1,-117},
            {91.1,-118},{96,-118},{96,-38},{130.8,-38}},                color={
            0,0,127}));
    connect(gaiHea.y, p_design_heat.u1) annotation (Line(points={{180.6,0},{180.6,
            0.4},{191.2,0.4}},       color={0,0,127}));
    connect(p_design_heat.y, multiplex2.u1[1]) annotation (Line(points={{200.4,-2},
            {206,-2},{206,-9.6},{209.2,-9.6}},   color={0,0,127}));
    connect(gaiCoo.y, p_design_cool.u2) annotation (Line(points={{182.6,-38},{182.6,
            -38.4},{191.2,-38.4}},       color={0,0,127}));
    connect(p_design_cool.y, multiplex2.u2[1]) annotation (Line(points={{200.4,-36},
            {204,-36},{204,-14.4},{209.2,-14.4}},      color={0,0,127}));
    connect(preHeaCoo.port, targetHeatCoolDemand) annotation (Line(points={{250,-12},
            {254,-12},{254,124},{114,124}}, color={191,0,0}));
    connect(TVen.port, infFlowThermalConductor.port_a)
      annotation (Line(points={{-130,76},{-86,76}}, color={191,0,0}));
    connect(infFlowThermalConductor.port_b, heaPorAir)
      annotation (Line(points={{-66,76},{-44,76}}, color={191,0,0}));
    connect(combiTimeTable.y[5], infFlowThermalConductor.mass_flow_rate)
      annotation (Line(points={{91.1,-117},{91.1,-118},{96,-118},{96,-4},{86,-4},{
            86,76},{18,76},{18,60},{-86.6,60},{-86.6,68}},
                  color={0,0,127}));
    connect(targetDHWDemand, combiTimeTable.y[13]) annotation (Line(points={{150,-90},
            {144,-90},{144,-120},{91.1,-120},{91.1,-117}}, color={0,0,127}));
    annotation (
      defaultComponentName = "zon",
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -140}, {140, 140}}), graphics={  Rectangle(fillColor = {95, 95, 95}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-140, 140}, {140, -140}}), Rectangle(lineColor = {117, 148, 176}, fillColor = DynamicSelect({170, 213, 255}, {((min(1, max(0, (1 - ((heaPorAir.T - 295.15)/10))))*28) + (min(1, max(0, ((heaPorAir.T - 295.15)/10)))*255)), (min(1, max(0, (1 - ((heaPorAir.T - 295.15)/10))))*108), (min(1, max(0, (1 - ((heaPorAir.T - 295.15)/10))))*200)}), fillPattern = FillPattern.Solid, extent = {{-120, 122}, {118, -122}}), Text(textColor = {0, 0, 255}, extent = {{-104, 174}, {118, 142}}, textString = "%name"), Text(extent = {{60, -82}, {118, -126}}, textString = "ISO"), Text(textColor = {0, 0, 88}, extent = {{88, 94}, {136, 62}}, textString = "TAir"), Text(textColor = {0, 0, 88}, extent = {{82, 18}, {130, -14}}, textString = "TSur"), Text(textColor = {255, 255, 255}, extent = {{48, 108}, {-72, 58}}, textString = DynamicSelect("", String((heaPorAir.T - 273.15), ".1f")))}),
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -140}, {140, 140}})),
      Documentation(info = "<html>
<p>
This is a lumped-capacity simplified building model based  on the 5R1C
network presented in the ISO 13790:2008 Standard. The simplified 5R1C model uses
five thermal resistances and one thermal capacity to reproduce the
transient thermal behaviour of buildings. The thermal zone is modeled with three
temperature nodes, the indoor air temperature <code>TAir</code>, the envelope internal
surface temperature <code>TSur</code> and the zone's mass temperature <code>TMas</code>
(the heat port is not shown in the figure), and two boundary
condition nodes, supply air temperature <code>TSup</code> and the external air temperature
<code>TExt</code>. The five resistances are related to heat transfer by ventilation <code>HVen</code>,
windows <code>HWin</code>, opaque components (split between <code>HTra</code> and <code>HMas</code>) and heat
transfer between the internal surfaces of walls and the air temperature <code>HThe</code>.
The thermal capacity <code>Cm</code> includes the thermal capacity of the entire zone. The heating and/or
cooling demand is found by calculating the heating and/or cooling power &Phi;HC that
needs to be supplied to, or extracted from, the internal air node to maintain a
certain set-point. Internal, &Phi;int , and solar, &Phi;sol, heat gains are input values,
which are split in three components.
</p>
<br>
<p align=\"center\">
<img src=\"modelica://Buildings/Resources/Images/ThermalZones/ISO13790/Zone/5R1CNetwork.png\" alt=\"image\"/>
</p>
<br>
The ventilation heat transfer coefficient <i>H<sub>ven</sub></i> is calculated using
<p align=\"center\" style=\"font-style:italic;\">
H<sub>ven</sub> = &rho;<sub>a</sub> c<sub>a</sub> &sum;<sub>k</sub>V&#775;<sub>k</sub>,
</p>
where <i>&rho;<sub>a</sub></i> is the density of air, <i>c<sub>a</sub></i> is the specific
heat capacity of air and <i>V&#775;<sub>k</sub></i> is the k-th volumetric external air
flow rate.
The coupling conductance <i>H<sub>the</sub></i> is given by
<p align=\"center\" style=\"font-style:italic;\">
H<sub>the</sub> = h<sub>as</sub> A<sub>tot</sub>,
</p>
where <i>h<sub>as</sub></i> is the heat transfer coefficient between the air
node the surface node, with a fixed value of <i>3.45 W/m<sup>2</sup>K</i>, and
<i>A<sub>tot</sub></i> is the area of all surfaces facing the building zone.
The thermal transmission coefficient of windows <i>H<sub>win</sub></i> is calculated using
<p align=\"center\" style=\"font-style:italic;\">
H<sub>win</sub> =
&sum;<sub>k</sub>U<sub>win,k</sub>A<sub>win,k</sub>,
</p>
where <i>U<sub>win,k</sub></i> is the thermal transmittance of window element
k of the building envelope and <i>A<sub>k</sub></i>  is the area of the window
element k of the building envelope. The coupling conductance <i>H<sub>mas</sub></i> is given by
<p align=\"center\" style=\"font-style:italic;\">
H<sub>mas</sub> =h<sub>ms</sub> f<sub>ms</sub> A<sub>f</sub>,
</p>
where <i>h<sub>ms</sub></i> is the heat transfer coefficient between the mass
node and the surface node, with fixed value of <i>9.1 W/m<sup>2</sup>K</i>,
<i>f<sub>ms</sub></i> is a correction factor, and <i>A<sub>f</sub></i>
is the floor area. The correction factor <i>f<sub>ms</sub></i> can be assumed as
<i>2.5</i> for light and medium building constructions, and <i>3</i> for heavy constructions.
The coupling conductance <i>H<sub>tra</sub></i> is calculated using
<p align=\"center\" style=\"font-style:italic;\">
H<sub>tra</sub> =
1 &frasl; (1 &frasl; H<sub>op</sub> - 1 &frasl; H<sub>mas</sub>),
</p>
where <i>H<sub>op</sub></i> is the thermal transmission coefficient of opaque elements.
The three heat gains components are calculated using
<p align=\"center\" style=\"font-style:italic;\">
&Phi;<sub>air</sub> = 0.5 &Phi;<sub>int</sub>,
</p>
<p align=\"center\" style=\"font-style:italic;\">
&Phi;<sub>sur</sub> = (1-f<sub>ms</sub> A<sub>f</sub> &frasl; A<sub>tot</sub>
-H<sub>win</sub> &frasl; h<sub>ms</sub> A<sub>tot</sub>)(0.5 &Phi;<sub>int</sub>+
&Phi;<sub>sol</sub>),
</p>
<p align=\"center\" style=\"font-style:italic;\">
&Phi;<sub>mas</sub> = f<sub>ms</sub> A<sub>f</sub> &frasl; A<sub>tot</sub> (0.5&Phi;<sub>int</sub> +
&Phi;<sub>sol</sub>).
</p>
<h4>Tips for parametrization</h4>
<ul>
<li>
The parameters <code>AWin</code>, <code>AWal</code>, <code>surTil</code> and <code>surAzi</code>
must have the same dimension of <code>nOrientations</code> .
</li>
<li>
The areas in <code>AWal</code> must account only for the opaque parts of the walls (excluding windows).
The floor and roof area is entered through <code>AFlo</code> and <code>ARoo</code>
and must not be entered as part of <code>AWal</code>.
</li>
<li>
If a wall contains only opaque parts, the corresponding window area must be set to <i>0</i>.
</li>
</ul>
</html>", revisions = "<html>
<ul>
<li>
Mar 16, 2022, by Alessandro Maccarini:<br/>
First implementation.
</li>
</ul>
</html>"));
  end Zone_EUReCA;

  model Building_EUReCA "Illustrates the use of the 5R1C thermal zone in free-floating conditions"
    parameter Modelica.Units.SI.ThermalConductance Htr_w "Heat transfer through glazed elements [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_is "Coupling conductance betwee air and surface nodes [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_ms "Heat transfer through opaque elements [W/K]";
    parameter Modelica.Units.SI.ThermalConductance Htr_em "Heat transfer through opaque elements [W/K]";
    parameter Modelica.Units.SI.HeatCapacity Cm "Zone thermal capacity [J/K]";
    parameter Modelica.Units.SI.HeatFlowRate cooling_design_flow_rate "Cooling maximum capacity [W]";
    parameter Modelica.Units.SI.HeatFlowRate heating_design_flow_rate "Heating maximum capacity [W]";
    parameter String path_to_sched_file "Path to schedule file";
    parameter String path_to_mos_weather_file "Path to weather file";


    Zone_EUReCA zona_EUReCA(
      Htr_w=Htr_w,
      Htr_is=Htr_is,
      Htr_ms=Htr_ms,
      Htr_em=Htr_em,
      Cm=Cm,
      cooling_design_flow_rate=cooling_design_flow_rate,
      heating_design_flow_rate=heating_design_flow_rate,
      path_to_sched_file=
          path_to_sched_file,
      path_to_mos_weather_file=path_to_mos_weather_file)
      annotation (Placement(transformation(extent={{-8,-14},{20,14}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b targetHeatCoolDemand
      "Heat port to air node"
      annotation (Placement(transformation(extent={{-10,92},{10,112}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    PrescribedHeatFlowDiverter prescribedHeatFlowDiverter(portA_percent=1.)
      annotation (Placement(transformation(extent={{-76,-10},{-56,10}})));
    Modelica.Blocks.Interfaces.RealOutput targetDHWDemand(final unit="W",
        displayUnit="W") annotation (Placement(transformation(origin={-44,0},
            extent={{140,-10},{160,10}}), transformation(extent={{140,-10},{160,
              10}})));
  equation
    connect(prescribedHeatFlowDiverter.portA, zona_EUReCA.heaPorAir) annotation (
        Line(points={{-56,4},{-24,4},{-24,20},{1.6,20},{1.6,7.6}},   color={191,0,
            0}));
    connect(prescribedHeatFlowDiverter.portB, zona_EUReCA.heaPorSur) annotation (
        Line(points={{-56,-4},{-14,-4},{-14,-0.4},{1.6,-0.4}},  color={191,0,0}));
    connect(port_a, prescribedHeatFlowDiverter.port_in)
      annotation (Line(points={{-100,0},{-76,0}}, color={191,0,0}));
    connect(zona_EUReCA.targetHeatCoolDemand, targetHeatCoolDemand) annotation (
       Line(points={{17.4,12.4},{18,12.4},{18,86},{0,86},{0,102}}, color={191,0,
            0}));
    connect(zona_EUReCA.targetDHWDemand, targetDHWDemand) annotation (Line(
          points={{21,-9},{92,-9},{92,0},{106,0}},   color={0,0,127}));
    annotation (
      experiment(
        StopTime=31536000,
        __Dymola_NumberOfIntervals=17520,
        Tolerance=1e-06,
        __Dymola_Algorithm="Dassl"),
      __Dymola_Commands(file = "modelica://Buildings/Resources/Scripts/Dymola/ThermalZones/ISO13790/Examples/FreeFloating.mos" "Simulate and plot"),
      Documentation(info = "<html>
  <p>
  This model illustrates the use of <a href=\"modelica://Buildings.ThermalZones.ISO13790.Zone5R1C.Zone\">
  Buildings.ThermalZones.ISO13790.Zone5R1C.Zone</a> in a free-floating case 
  (i.e. no heating or cooling)
  </p>
  </html>", revisions = "<html>
  <ul>
  <li>
  Mar 16, 2022, by Alessandro Maccarini:<br/>
  First implementation.
  </li>
  </ul>
  </html>"),
      Icon(graphics={
          Rectangle(extent={{-100,40},{100,-100}}, lineColor={255,128,0}),
          Polygon(
            points={{-120,40},{122,40},{0,120},{-120,40}},
            lineColor={28,108,200},
            fillColor={255,170,85},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-76,-2},{70,-40}},
            textColor={255,128,0},
            textString="BETA_Lab
EUReCA")}));
  end Building_EUReCA;

  model InfFlowThermalConductor
    "Lumped thermal element transporting heat without storing it"
    extends Modelica.Thermal.HeatTransfer.Interfaces.Element1D;

    Modelica.Blocks.Interfaces.RealInput mass_flow_rate "mass flow rate for infiltration input [kg/s]"
      annotation (Placement(transformation(extent={{-126,-100},{-86,-60}})));
  equation
    Q_flow = mass_flow_rate*1005*dT;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={
          Rectangle(
            extent={{-90,70},{90,-70}},
            pattern=LinePattern.None,
            fillColor={192,192,192},
            fillPattern=FillPattern.Backward),
          Line(
            points={{-90,70},{-90,-70}},
            thickness=0.5),
          Line(
            points={{90,70},{90,-70}},
            thickness=0.5),
          Text(
            extent={{-150,120},{150,80}},
            textString="%name",
            textColor={0,0,255}),
          Text(
            extent={{-150,-80},{150,-110}},
            textString="G=%G")}),
      Documentation(info="<html>
<p>
This is a model for transport of heat without storing it; see also:
<a href=\"modelica://Modelica.Thermal.HeatTransfer.Components.ThermalResistor\">ThermalResistor</a>.
It may be used for complicated geometries where
the thermal conductance G (= inverse of thermal resistance)
is determined by measurements and is assumed to be constant
over the range of operations. If the component consists mainly of
one type of material and a regular geometry, it may be calculated,
e.g., with one of the following equations:
</p>
<ul>
<li><p>
    Conductance for a <strong>box</strong> geometry under the assumption
    that heat flows along the box length:</p>
    <blockquote><pre>
G = k*A/L
k: Thermal conductivity (material constant)
A: Area of box
L: Length of box
    </pre></blockquote>
    </li>
<li><p>
    Conductance for a <strong>cylindrical</strong> geometry under the assumption
    that heat flows from the inside to the outside radius
    of the cylinder:</p>
    <blockquote><pre>
G = 2*pi*k*L/log(r_out/r_in)
pi   : Modelica.Constants.pi
k    : Thermal conductivity (material constant)
L    : Length of cylinder
log  : Modelica.Math.log;
r_out: Outer radius of cylinder
r_in : Inner radius of cylinder
    </pre></blockquote>
    </li>
</ul>
<blockquote><pre>
Typical values for k at 20 degC in W/(m.K):
  aluminium   220
  concrete      1
  copper      384
  iron         74
  silver      407
  steel        45 .. 15 (V2A)
  wood         0.1 ... 0.2
</pre></blockquote>
</html>"));
  end InfFlowThermalConductor;

  model PrescribedHeatFlowDiverter
    "Prescribed heat flow boundary condition"
    parameter Real portA_percent(min = 0, max=1) "Flux percentage to port A";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_in
          annotation (Placement(transformation(
          origin={-100,0},
          extent={{20,-20},{-20,20}},
          rotation=180)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b portA
      annotation (Placement(transformation(extent={{90,30},{110,50}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b portB
      annotation (Placement(transformation(extent={{90,-50},{110,-30}})));
  equation
    portA.Q_flow = -port_in.Q_flow*portA_percent;
    portB.Q_flow = -port_in.Q_flow*(1-portA_percent);
    portB.T = port_in.T;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={
          Line(
            points={{-60,-20},{40,-20}},
            color={191,0,0},
            thickness=0.5),
          Line(
            points={{-60,20},{40,20}},
            color={191,0,0},
            thickness=0.5),
          Line(
            points={{-80,0},{-60,-20}},
            color={191,0,0},
            thickness=0.5),
          Line(
            points={{-80,0},{-60,20}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{40,0},{40,40},{70,20},{40,0}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{40,-40},{40,0},{70,-20},{40,-40}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{70,40},{90,-40}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,100},{150,60}},
            textString="%name",
            textColor={0,0,255})}),
      Documentation(info="<html>
<p>
This model allows a specified amount of heat flow rate to be \"injected\"
into a thermal system at a given port.  The amount of heat
is given by the input signal Q_flow into the model. The heat flows into the
component to which the component PrescribedHeatFlow is connected,
if the input signal is positive.
</p>
<p>
If parameter alpha is &lt;&gt; 0, the heat flow is multiplied by (1 + alpha*(port.T - T_ref))
in order to simulate temperature dependent losses (which are given with respect to reference temperature T_ref).
</p>
</html>"));
  end PrescribedHeatFlowDiverter;

  annotation (
    uses(Modelica(version = "4.0.0"), Buildings(version = "10.0.0")));
end EUReCA2Modelica;
