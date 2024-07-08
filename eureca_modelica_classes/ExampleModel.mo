within ;
model ExampleModel
  EUReCA2Modelica.Building_EUReCA building_EUReCA(
    Htr_w=500,
    Htr_is=3000,
    Htr_ms=7000,
    Htr_em=2000,
    Cm=60000000,
    cooling_design_flow_rate=-10000,
    heating_design_flow_rate=50000,
    path_to_sched_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/data_Zone_1.txt",
    path_to_mos_weather_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/ITA_Venezia-Tessera.161050_IGDG.mos")
    annotation (Placement(transformation(extent={{30,52},{50,72}})));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=31536000,
      Interval=3600,
      __Dymola_Algorithm="Dassl"));
end ExampleModel;
