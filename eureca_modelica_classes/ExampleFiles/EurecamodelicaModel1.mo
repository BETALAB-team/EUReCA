within ;
model EurecaModelicaModel1

  EUReCA2Modelica.Building_EUReCA bd_EUReCA_Bd_1(
    Htr_w=307.08,
    Htr_is=13768.02,
    Htr_ms=32552.29,
    Htr_em=2135.74,
    Cm=345521229.59,
    cooling_design_flow_rate=-22779,
    heating_design_flow_rate=87113,
    path_to_sched_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/data_Zone_1.txt",
    path_to_mos_weather_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/ITA_Venezia-Tessera.161050_IGDG.mos")
    annotation (Placement(transformation(extent={{50,50},{70,70}})));
  


  EUReCA2Modelica.Building_EUReCA bd_EUReCA_Bd_2(
    Htr_w=153.54,
    Htr_is=12903.17,
    Htr_ms=30803.47,
    Htr_em=1787.98,
    Cm=330728909.55,
    cooling_design_flow_rate=-20587,
    heating_design_flow_rate=71906,
    path_to_sched_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/data_Zone_1.txt",
    path_to_mos_weather_file=
        "C:/Users/pratenr82256/Desktop/ModelicaEUReCA/ClassiDaEUReCA/ITA_Venezia-Tessera.161050_IGDG.mos")
    annotation (Placement(transformation(extent={{50,-70},{70,-50}})));
  

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=31536000,
      Interval=3600,
      __Dymola_Algorithm="Dassl"))
end EurecaModelicaModel1;
;