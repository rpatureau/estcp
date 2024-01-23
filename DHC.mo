within estcp;
package DHC "Models for district heating and cooling systems"
  extends Modelica.Icons.VariantsPackage;

  package Loads "Models for computing thermal loads served by the DES"
    extends Modelica.Icons.VariantsPackage;

    package Combined "Package with models for loads"
      extends Modelica.Icons.VariantsPackage;

      model BuildingTimeSeriesWithETS_chiller
        "Model of a building with loads provided as time series, connected to an ETS"
        extends
          estcp.DHC.Loads.Combined.BaseClasses.PartialBuildingWithETS_chiller(
            redeclare
            Buildings.Experimental.DHC.Loads.BaseClasses.Examples.BaseClasses.BuildingTimeSeries bui(
            final filNam=filNam,
            T_aHeaWat_nominal=THeaWatSup_nominal,
            T_bHeaWat_nominal=THeaWatSup_nominal - 5,
            T_aChiWat_nominal=TChiWatSup_nominal,
            T_bChiWat_nominal=TChiWatSup_nominal + 5), ets(
            QChiWat_flow_nominal=QCoo_flow_nominal,
            QHeaWat_flow_nominal=QHea_flow_nominal,
            QHotWat_flow_nominal=QHot_flow_nominal));
        parameter String filNam
          "Library path of the file with thermal loads as time series";
        final parameter Modelica.Units.SI.HeatFlowRate QCoo_flow_nominal(
          max=-Modelica.Constants.eps)=
          bui.facMul * bui.QCoo_flow_nominal
          "Space cooling design load (<=0)"
          annotation (Dialog(group="Design parameter"));
        final parameter Modelica.Units.SI.HeatFlowRate QHea_flow_nominal(
          min=Modelica.Constants.eps)=
          bui.facMul * bui.QHea_flow_nominal
          "Space heating design load (>=0)"
          annotation (Dialog(group="Design parameter"));
        final parameter Modelica.Units.SI.HeatFlowRate QHot_flow_nominal(
          min=Modelica.Constants.eps)=
          bui.facMul * Buildings.Experimental.DHC.Loads.BaseClasses.getPeakLoad(
            string="#Peak water heating load",
            filNam=Modelica.Utilities.Files.loadResource(filNam))
          "Hot water design load (>=0)"
          annotation (Dialog(group="Design parameter"));
        Buildings.Controls.OBC.CDL.Reals.MultiplyByParameter loaHeaNor(
          k=1/QHea_flow_nominal) "Normalized heating load"
          annotation (Placement(transformation(extent={{-200,-110},{-180,-90}})));
        Buildings.Controls.OBC.CDL.Reals.GreaterThreshold enaHeaCoo[2](each t=1e-4)
          "Threshold comparison to enable heating and cooling"
          annotation (Placement(transformation(extent={{-110,-130},{-90,-110}})));
        Buildings.Controls.OBC.CDL.Reals.MultiplyByParameter loaCooNor(k=1/
              QCoo_flow_nominal) "Normalized cooling load"
          annotation (Placement(transformation(extent={{-200,-150},{-180,-130}})));
      equation
        connect(enaHeaCoo[1].y, ets.uHea) annotation (Line(points={{-88,-120},{-40,-120},
                {-40,-48},{-34,-48}},       color={255,0,255}));
        connect(enaHeaCoo[2].y, ets.uCoo) annotation (Line(points={{-88,-120},{-40,-120},
                {-40,-54},{-34,-54}},       color={255,0,255}));
        connect(loaHeaNor.y, enaHeaCoo[1].u) annotation (Line(points={{-178,-100},{
                -120,-100},{-120,-120},{-112,-120}}, color={0,0,127}));
        connect(loaCooNor.y, enaHeaCoo[2].u) annotation (Line(points={{-178,-140},{
                -120,-140},{-120,-120},{-112,-120}}, color={0,0,127}));
        connect(bui.QReqHea_flow, loaHeaNor.u) annotation (Line(points={{20,4},{20,-6},
                {-218,-6},{-218,-100},{-202,-100}}, color={0,0,127}));
        connect(bui.QReqCoo_flow, loaCooNor.u) annotation (Line(points={{24,4},{24,-4},
                {-220,-4},{-220,-140},{-202,-140}}, color={0,0,127}));
        connect(loaHeaNor.y, resTHeaWatSup.u) annotation (Line(points={{-178,-100},{
                -120,-100},{-120,-40},{-112,-40}},  color={0,0,127}));
        annotation (
          Documentation(info="<html>
<p>
This model is composed of a heat pump based energy transfer station model
<a href=\"modelica://Buildings.Experimental.DHC.EnergyTransferStations.Combined.HeatPumpHeatExchanger\">
Buildings.Experimental.DHC.EnergyTransferStations.Combined.HeatPumpHeatExchanger</a>
connected to a simplified building model where the space heating, cooling
and hot water loads are provided as time series.
</p>
<h4>Scaling</h4>
<p>
The parameter <code>bui.facMul</code> is the multiplier factor
applied to the building loads that are provided as time series.
It is used to represent multiple identical buildings served by
a unique energy transfer station.
The parameter <code>facMul</code> is the multiplier factor
applied to the whole system composed of the building(s) and the
energy transfer station.
It is used to represent multiple identical ETSs served by
the DHC system.
So, if for instance the overall heating and cooling efficiency is
equal to <i>1</i>, then the load on the district loop
is the load provided as time series multiplied by <i>facMul * bui.facMul</i>.
</p>
<p>
Note that the parameters <code>QCoo_flow_nominal</code>, <code>QHea_flow_nominal</code>
and <code>QHot_flow_nominal</code> are the <i>ETS</i> design values. They include
the building loads multiplier factor <code>bui.facMul</code> but not
the building and ETS multiplier factor <code>facMul</code>.
</p>
</html>",       revisions="<html>
<ul>
<li>
November 21, 2022, by David Blum:<br/>
Change <code>bui.facMulHea</code> and <code>bui.facMulCoo</code> to be default.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/2302\">
issue 2302</a>.
</li>
<li>
February 23, 2021, by Antoine Gautier:<br/>
First implementation.
</li>
</ul>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-300,-300},{
                  300,300}})));
      end BuildingTimeSeriesWithETS_chiller;

      package Examples "This package contains example models"
        extends Modelica.Icons.ExamplesPackage;

        model BuildingTimeSeriesWithETS_chiller
          "Example model of a building with loads provided as time series for heat pump heating and free cooling in an ambient district network"
          extends Modelica.Icons.Example;
          package Medium=Buildings.Media.Water
            "Medium model";
          Buildings.Fluid.Sources.Boundary_pT supAmbWat(
            redeclare package Medium = Medium,
            p(displayUnit="bar"),
            use_T_in=true,
            T=280.15,
            nPorts=1) "Ambient water supply"
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},rotation=0,origin={-50,-10})));
          Buildings.Fluid.Sources.Boundary_pT sinAmbWat(
            redeclare package Medium = Medium,
            p(displayUnit="bar"),
            nPorts=1) "Sink for ambient water"
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},rotation=0,origin={-50,-70})));
          Buildings.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium = Medium)
            "Mass flow rate sensor"
            annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
          Modelica.Blocks.Sources.Constant TDisSup(k(
              unit="K",
              displayUnit="degC") = 288.15)
            "District supply temperature"
            annotation (Placement(transformation(extent={{-92,-16},{-72,4}})));
          estcp.DHC.Loads.Combined.BuildingTimeSeriesWithETS_chiller bui(
            redeclare package MediumSer = Medium,
            redeclare package MediumBui = Medium,
            bui(facMul=10),
            allowFlowReversalSer=true,
            filNam=
                "modelica://Buildings/Resources/Data/Experimental/DHC/Loads/Examples/SwissOffice_20190916.mos")
            annotation (Placement(transformation(extent={{40,-20},{60,0}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant TChiWatSupSet(k=bui.TChiWatSup_nominal)
            "Chilled water supply temperature set point"
            annotation (Placement(transformation(extent={{-90,60},{-70,80}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant THeaWatSupMaxSet(k=bui.THeaWatSup_nominal)
            "Heating water supply temperature set point - Maximum value"
            annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant THeaWatSupMinSet(
           k(final unit="K",
             displayUnit="degC") = 301.15)
           "Heating water supply temperature set point - Minimum value"
            annotation (Placement(transformation(extent={{0,60},{20,80}})));
        equation
          connect(supAmbWat.ports[1], senMasFlo.port_a)
            annotation (Line(points={{-40,-10},{-20,-10}},
                                                         color={0,127,255}));
          connect(TDisSup.y,supAmbWat. T_in)
            annotation (Line(points={{-71,-6},{-62,-6}}, color={0,0,127}));
          connect(senMasFlo.port_b, bui.port_aSerAmb) annotation (Line(points={{0,-10},
                  {40,-10}},              color={0,127,255}));
          connect(sinAmbWat.ports[1], bui.port_bSerAmb) annotation (Line(points={{-40,-70},
                  {70,-70},{70,-10},{60,-10}}, color={0,127,255}));
          connect(THeaWatSupMinSet.y, bui.THeaWatSupMinSet) annotation (Line(points={{22,70},
                  {34,70},{34,0},{38,0},{38,-1}},                 color={0,0,127}));
          connect(THeaWatSupMaxSet.y, bui.THeaWatSupMaxSet) annotation (Line(points={{-18,70},
                  {-10,70},{-10,34},{32,34},{32,-3},{38,-3}},     color={0,0,127}));
          connect(TChiWatSupSet.y, bui.TChiWatSupSet) annotation (Line(points={{-68,70},
                  {-52,70},{-52,32},{30,32},{30,-5},{38,-5}},
                                                        color={0,0,127}));
          annotation (
            Icon(
              coordinateSystem(
                preserveAspectRatio=false)),
            Diagram(
                coordinateSystem(
                preserveAspectRatio=false)),
            __Dymola_Commands(
              file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Loads/Combined/Examples/BuildingTimeSeriesWithETS.mos" "Simulate and plot"),
            experiment(
              StopTime=864000,
              Tolerance=1e-06),
            Documentation(info="<html>
<p>
Example model of a building with loads provided as time series for heat
pump space heating, heat pump domestic hot water heating,
and free cooling in an ambient district network.
</p>
</html>",         revisions="<html>
<ul>
<li>
May 3, 2023, by David Blum:<br/>
First implementation.
</li>
</ul>
</html>"));
        end BuildingTimeSeriesWithETS_chiller;
        annotation (Documentation(info="<html>
<p>
This package contains an example illustrating the use of the model in
<a href=\"modelica://Buildings.Experimental.DHC.Loads.Combined\">
Buildings.Experimental.DHC.Loads.Combined</a>.
</p>
</html>"));
      end Examples;

      package BaseClasses "Package with base classes that are used by multiple models"
        extends Modelica.Icons.BasesPackage;

        model PartialBuildingWithETS_chiller
          "Partial model with ETS model and partial building model"
          extends Buildings.Experimental.DHC.Loads.BaseClasses.PartialBuildingWithPartialETS(
            nPorts_heaWat=1,
            nPorts_chiWat=1,
            redeclare Buildings.Experimental.DHC.EnergyTransferStations.Combined.ChillerBorefield ets(
            hex(show_T=true),
              WSE(show_T=true),
              conCon=Buildings.Experimental.DHC.EnergyTransferStations.Types.ConnectionConfiguration.Pump,
              dp1Hex_nominal=20E3,
              dp2Hex_nominal=20E3,
              QHex_flow_nominal=abs(QChiWat_flow_nominal),
              T_a1Hex_nominal=282.15,
              T_b1Hex_nominal=278.15,
              T_a2Hex_nominal=276.15,
              T_b2Hex_nominal=280.15,
              have_WSE=true,
              QWSE_flow_nominal=QChiWat_flow_nominal,
              dpCon_nominal=15E3,
              dpEva_nominal=15E3,
              final datChi=datChi,
              T_a1WSE_nominal=281.15,
              T_b1WSE_nominal=286.15,
              T_a2WSE_nominal=288.15,
              T_b2WSE_nominal=283.15));
          parameter Modelica.Units.SI.TemperatureDifference dT_nominal(min=0) = 4
            "Water temperature drop/increase accross load and source-side HX (always positive)"
            annotation (Dialog(group="ETS model parameters"));
          parameter Modelica.Units.SI.Temperature TChiWatSup_nominal=18 + 273.15
            "Chilled water supply temperature"
            annotation (Dialog(group="ETS model parameters"));
          parameter Modelica.Units.SI.Temperature THeaWatSup_nominal=38 + 273.15
            "Heating water supply temperature"
            annotation (Dialog(group="ETS model parameters"));
          parameter Modelica.Units.SI.Pressure dp_nominal=50000
            "Pressure difference at nominal flow rate (for each flow leg)"
            annotation (Dialog(group="ETS model parameters"));
          parameter Real COPHeaWat_nominal(final unit="1") = 4.0
            "COP of heat pump for heating water production"
            annotation (Dialog(group="ETS model parameters"));
          parameter Real COPHotWat_nominal(final unit="1") = 2.3
            "COP of heat pump for hot water production"
            annotation (Dialog(group="ETS model parameters", enable=have_hotWat));
          // IO CONNECTORS
          Buildings.Controls.OBC.CDL.Interfaces.RealInput TChiWatSupSet(
            final unit="K",
            displayUnit="degC")
            "Chilled water supply temperature set point"
            annotation (Placement(
                transformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-320,80}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-120,50})));
          Buildings.Controls.OBC.CDL.Interfaces.RealInput THeaWatSupMaxSet(
            final unit="K",
            displayUnit="degC")
            "Heating water supply temperature set point - Maximum value"
            annotation (
              Placement(transformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-320,120}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-120,70})));
          Buildings.Controls.OBC.CDL.Interfaces.RealInput THeaWatSupMinSet(
            final unit="K",
            displayUnit="degC")
            "Heating water supply temperature set point - Minimum value"
            annotation (
              Placement(transformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-320,160}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=0,
                origin={-120,90})));
          // COMPONENTS
          Buildings.Controls.OBC.CDL.Reals.Line resTHeaWatSup
            "HW supply temperature reset"
            annotation (Placement(transformation(extent={{-110,-50},{-90,-30}})));
          Buildings.Controls.OBC.CDL.Reals.Sources.Constant zer(k=0)
            "Zero"
            annotation (Placement(transformation(extent={{-180,-30},{-160,-10}})));
          Buildings.Controls.OBC.CDL.Reals.Sources.Constant one(k=1)
            "One"
            annotation (Placement(transformation(extent={{-180,-70},{-160,-50}})));
          Buildings.Controls.OBC.CDL.Reals.MultiplyByParameter mulPPumETS(u(final
                unit="W"), final k=facMul) if have_pum "Scaling"
            annotation (Placement(transformation(extent={{270,-10},{290,10}})));
          Buildings.Controls.OBC.CDL.Interfaces.RealOutput PPumETS(
            final unit="W") if have_pum
            "ETS pump power"
            annotation (Placement(
                transformation(extent={{300,-20},{340,20}}),iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=90,
                origin={70,120})));
          parameter Buildings.Fluid.Chillers.Data.ElectricEIR.Generic datChi(
            QEva_flow_nominal=QChiWat_flow_nominal,
            COP_nominal=3.8,
            PLRMax=1,
            PLRMinUnl=0.3,
            PLRMin=0.3,
            etaMotor=1,
            mEva_flow_nominal=abs(QChiWat_flow_nominal)/4186/4,
            mCon_flow_nominal=abs(QChiWat_flow_nominal)*(1+1/datChi.COP_nominal)/4186/8,
            TEvaLvg_nominal=276.15,
            capFunT={1.72,0.02,0,-0.02,0,0},
            EIRFunT={0.28,-0.02,0,0.02,0,0},
            EIRFunPLR={0.1,0.9,0},
            TEvaLvgMin=276.15,
            TEvaLvgMax=288.15,
            TConEnt_nominal=315.15,
            TConEntMin=291.15,
            TConEntMax=328.15)
            "Chiller performance data"
            annotation (Placement(transformation(extent={{-250,260},{-230,280}})));
        equation
          connect(TChiWatSupSet, ets.TChiWatSupSet) annotation (Line(points={{-320,80},{
                  -132,80},{-132,-66},{-34,-66}},color={0,0,127}));
          connect(resTHeaWatSup.y, ets.THeaWatSupSet) annotation (Line(points={{-88,-40},
                  {-60,-40},{-60,-60},{-34,-60}}, color={0,0,127}));
          connect(THeaWatSupMaxSet, resTHeaWatSup.f2) annotation (Line(points={{-320,120},
                  {-280,120},{-280,-48},{-112,-48}}, color={0,0,127}));
          connect(THeaWatSupMinSet, resTHeaWatSup.f1) annotation (Line(points={{-320,160},
                  {-276,160},{-276,-36},{-112,-36}}, color={0,0,127}));
          connect(one.y, resTHeaWatSup.x2) annotation (Line(points={{-158,-60},{-126,-60},
                  {-126,-44},{-112,-44}}, color={0,0,127}));
          connect(zer.y, resTHeaWatSup.x1) annotation (Line(points={{-158,-20},{-116,-20},
                  {-116,-32},{-112,-32}}, color={0,0,127}));
          connect(mulPPumETS.y, PPumETS)
            annotation (Line(points={{292,0},{320,0}},   color={0,0,127}));
          connect(ets.PPum, mulPPumETS.u) annotation (Line(points={{34,-60},{240,-60},{
                  240,0},{268,0}},   color={0,0,127}));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)),
            Documentation(info="<html>
<p>
This model is composed of a heat pump based energy transfer station model 
<a href=\"modelica://Buildings.Experimental.DHC.EnergyTransferStations.Combined.HeatPumpHeatExchanger\">
Buildings.Experimental.DHC.EnergyTransferStations.Combined.HeatPumpHeatExchanger</a>
connected to a repleacable building load model. 
</p>
</html>",         revisions="<html>
<ul>
<li>
February 23, 2021, by Antoine Gautier:<br/>
First implementation.
</li>
</ul>
</html>"));
        end PartialBuildingWithETS_chiller;
      annotation (Documentation(info="<html>
<p>
This package contains base classes that are used to construct the classes in
<a href=\"modelica://Buildings.Experimental.DHC.Loads.Combined\">
Buildings.Experimental.DHC.Loads.Combined</a>.
</p>
</html>"));
      end BaseClasses;
    annotation (Documentation(info="<html>
<p>
This package contains models of building loads that are used to
build example models of DHC systems.
</p>
</html>"));
    end Combined;

    annotation (
      preferredView="info",
      Documentation(
        info="<html>
<p>
This package contains models for the thermal and domestic hot water demand
prediction in buildings.
</p>
</html>"));
  end Loads;

  package Plants "Package of models for central plants"
    extends Modelica.Icons.VariantsPackage;

    package Reservoir "Package of models for district-scale thermal reservoirs"
      extends Modelica.Icons.Package;
      model BoreField "Geothermal borefield model"
        extends Buildings.Fluid.Geothermal.Borefields.TwoUTubes(
          final energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          final tLoaAgg(displayUnit="h") = 3600,
          final nSeg=5,
          TExt0_start=282.55,
          final z0=10,
          final dT_dz=0.02,
          final dynFil=true,
          borFieDat(
            filDat=Buildings.Fluid.Geothermal.Borefields.Data.Filling.Bentonite(
                kFil=2.0,
                cFil=3040,
                dFil=1450),
            soiDat=Buildings.Fluid.Geothermal.Borefields.Data.Soil.SandStone(
                kSoi=2.3,
                cSoi=1000,
                dSoi=2600),
            conDat=Buildings.Fluid.Geothermal.Borefields.Data.Configuration.Example(
              borCon=Buildings.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.DoubleUTubeParallel,
              dp_nominal=35000,
              hBor=300,
              rBor=0.095,
              nBor=350,
              cooBor=cooBor,
              dBor=1,
              rTub=0.02,
              kTub=0.5,
              eTub=0.0037,
              xC=0.05)),
          show_T=true);
        /*
  Some parameters (such as nBor) cannot be propagated down to
  borFieDat.conDat otherwise Dymola fails to expand.
  We assign them literally within borFieDat.conDat and propagate them up here
  to compute dependent parameters.
  */
        parameter Integer nBor = borFieDat.conDat.nBor
          "Number of boreholes"
          annotation(Evaluate=true);
        parameter Real dxyBor = 10
          "Distance between boreholes";
        final parameter Modelica.Units.SI.Length cooBor[nBor,2]={dxyBor*{mod(i - 1,
            10),floor((i - 1)/10)} for i in 1:nBor}
          "Cartesian coordinates of the boreholes in meters";
        Buildings.Controls.OBC.CDL.Interfaces.RealOutput Q_flow(final unit="W")
          "Rate at which heat is extracted from soil"
          annotation (Placement(transformation(extent={{100,-50},{120,-30}}),
            iconTransformation(extent={{100,-64},{140,-24}})));
      equation
        connect(gaiQ_flow.y, Q_flow) annotation (Line(points={{1,80},{14,80},{14,54},{
                96,54},{96,-40},{110,-40}},
                                          color={0,0,127}));
        annotation (Documentation(info="<html>
<p>
This model represents a borefield composed of 350 boreholes,
with the following main assumptions.
</p>
<ul>
<li>
The soil is made of sandstone.
</li>
<li>
The boreholes are filled with a bentonite grout.
</li>
<li>
The boreholes have a height of 300 m and a diameter of 190 mm.
They are discretized vertically in five segments.
</li>
<li>
A distance of 10 m between each borehole is considered.
</li>
<li>
HDPE pipes with a diameter of 40 mm are considered, in a
double U-tube parallel configuration.
</li>
</ul>
</html>",       revisions="<html>
<ul>
<li>
May 31, 2023, by Michael Wetter:<br/>
Removed <code>final</code> modifier for <code>borFieDat</code> to allow record to be replaced
in models that extend this model.
</li>
<li>
February 23, 2021, by Antoine Gautier:<br/>
Updated documentation.
</li>
<li>
January 12, 2020, by Michael Wetter:<br/>
Added documentation.
</li>
</ul>
</html>"));
      end BoreField;
    annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains models for district scale thermal reservoirs.
</p>
</html>"));
    end Reservoir;

    annotation (
      preferredView="info",
      Documentation(
        info="<html>
<p>
This package contains models for central plants.
</p>
</html>"));
  end Plants;

  package Examples "Collection of cases study"
    extends Modelica.Icons.ExamplesPackage;

    package Combined "Package of example models for DHC systems"
      extends Modelica.Icons.VariantsPackage;

      model SeriesConstantFlow_chiller
        "Example of series connection with constant district water mass flow rate"
        extends estcp.DHC.Examples.Combined.BaseClasses.PartialSeries_chiller(
            redeclare estcp.DHC.Loads.Combined.BuildingTimeSeriesWithETS_chiller bui[
            nBui](final filNam=filNam), datDes(
            mPumDis_flow_nominal=95,
            mPipDis_flow_nominal=95,
            dp_length_nominal=250,
            epsPla=0.935));
        parameter String filNam[nBui]={
          "modelica://Buildings/Resources/Data/Experimental/DHC/Loads/Examples/SwissOffice_20190916.mos",
          "modelica://Buildings/Resources/Data/Experimental/DHC/Loads/Examples/SwissResidential_20190916.mos",
          "modelica://Buildings/Resources/Data/Experimental/DHC/Loads/Examples/SwissHospital_20190916.mos"}
          "Library paths of the files with thermal loads as time series";
        Modelica.Blocks.Sources.Constant masFloMaiPum(
          k=datDes.mPumDis_flow_nominal)
          "Distribution pump mass flow rate"
          annotation (Placement(transformation(extent={{-280,-70},{-260,-50}})));
        Modelica.Blocks.Sources.Constant masFloDisPla(
          k=datDes.mPla_flow_nominal)
          "District water flow rate to plant"
          annotation (Placement(transformation(extent={{-250,10},{-230,30}})));
      equation
        connect(masFloMaiPum.y, pumDis.m_flow_in) annotation (Line(points={{-259,-60},
                {60,-60},{60,-60},{68,-60}}, color={0,0,127}));
        connect(pumSto.m_flow_in, masFloMaiPum.y) annotation (Line(points={{-180,-68},
                {-180,-60},{-259,-60}}, color={0,0,127}));
        connect(masFloDisPla.y, pla.mPum_flow) annotation (Line(points={{-229,20},
                {-184,20},{-184,4.66667},{-161.333,4.66667}},
                                        color={0,0,127}));
        annotation (
        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-360,-260},{360,260}})),
          __Dymola_Commands(
        file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Examples/Combined/SeriesConstantFlow.mos"
        "Simulate and plot"),
        experiment(
            StopTime=604800,
            Tolerance=1e-06),
          Documentation(info="<html>
<p>
This is a model of a so-called \"reservoir network\" (Sommer 2020), i.e., a fifth
generation district system with unidirectional mass flow rate in the
district loop, and energy transfer stations connected in series.
In this model, the temperature of the district loop is stabilized through
the operation of the plant and the borefield.
The main circulation pump has a constant mass flow rate.
Each substation takes water from the main district loop and feeds its return water back
into the main district loop downstream from the intake.
The pipes of the main loop are designed for a pressure drop of
<code>dpDis_length_nominal=250</code> Pa/m at the design flow rate.
</p>
<h4>References</h4>
<p>
Sommer T., Sulzer M., Wetter M., Sotnikov A., Mennel S., Stettler C.
<i>The reservoir network: A new network topology for district heating
and cooling.</i>
Energy, Volume 199, 15 May 2020, 117418.
</p>
</html>",       revisions="<html>
<ul>
<li>
February 23, 2021, by Antoine Gautier:<br/>
Refactored with base classes from the <code>DHC</code> package.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1769\">
issue 1769</a>.
</li>
<li>
January 12, 2020, by Michael Wetter:<br/>
Added documentation.
</li>
</ul>
</html>"));
      end SeriesConstantFlow_chiller;

      package BaseClasses "Package with base classes that are used by multiple models"
        extends Modelica.Icons.BasesPackage;

        partial model PartialSeries_chiller "Partial model for series network"
          extends Modelica.Icons.Example;
          package Medium = Buildings.Media.Water "Medium model";
          constant Real facMul = 10
            "Building loads multiplier factor";
          parameter Real dpDis_length_nominal(final unit="Pa/m") = 250
            "Pressure drop per pipe length at nominal flow rate - Distribution line";
          parameter Real dpCon_length_nominal(final unit="Pa/m") = 250
            "Pressure drop per pipe length at nominal flow rate - Connection line";
          parameter Boolean allowFlowReversalSer = true
            "Set to true to allow flow reversal in the service lines"
            annotation(Dialog(tab="Assumptions"), Evaluate=true);
          parameter Boolean allowFlowReversalBui = false
            "Set to true to allow flow reversal for in-building systems"
            annotation(Dialog(tab="Assumptions"), Evaluate=true);
          parameter Integer nBui = datDes.nBui
            "Number of buildings connected to DHC system"
            annotation (Evaluate=true);
          inner parameter
            Buildings.Experimental.DHC.Examples.Combined.BaseClasses.DesignDataSeries datDes(final
              mCon_flow_nominal=bui.ets.hex.m1_flow_nominal)  "Design data"
            annotation (Placement(transformation(extent={{-340,220},{-320,240}})));
          // COMPONENTS
          estcp.DHC.Plants.Reservoir.BoreField borFie(redeclare final package
              Medium = Medium) "Bore field" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-130,-80})));
          Buildings.Experimental.DHC.EnergyTransferStations.BaseClasses.Pump_m_flow pumDis(
            redeclare final package Medium = Medium,
            final m_flow_nominal=datDes.mPumDis_flow_nominal,
            final allowFlowReversal=allowFlowReversalSer,
            dp_nominal=150E3) "Distribution pump" annotation (Placement(
                transformation(
                extent={{10,-10},{-10,10}},
                rotation=90,
                origin={80,-60})));
          Buildings.Fluid.Sources.Boundary_pT bou(
            redeclare final package Medium=Medium,
            final nPorts=1)
            "Boundary pressure condition representing the expansion vessel"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={112,-20})));
          Buildings.Experimental.DHC.EnergyTransferStations.BaseClasses.Pump_m_flow pumSto(
              redeclare final package Medium = Medium, m_flow_nominal=datDes.mSto_flow_nominal)
            "Bore field pump" annotation (Placement(transformation(
                extent={{10,10},{-10,-10}},
                rotation=180,
                origin={-180,-80})));
          Buildings.Experimental.DHC.Networks.Combined.BaseClasses.ConnectionSeriesStandard conPla(
            redeclare final package Medium = Medium,
            final mDis_flow_nominal=datDes.mPipDis_flow_nominal,
            final mCon_flow_nominal=datDes.mPla_flow_nominal,
            lDis=0,
            lCon=0,
            final dhDis=0.2,
            final dhCon=0.2,
            final allowFlowReversal=allowFlowReversalSer)
            "Connection to the plant (pressure drop lumped in plant and network model)"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={-80,-10})));
          Buildings.Experimental.DHC.Networks.Combined.BaseClasses.ConnectionSeriesStandard conSto(
            redeclare final package Medium = Medium,
            final mDis_flow_nominal=datDes.mPipDis_flow_nominal,
            final mCon_flow_nominal=datDes.mSto_flow_nominal,
            lDis=0,
            lCon=0,
            final dhDis=0.2,
            final dhCon=0.2,
            final allowFlowReversal=allowFlowReversalSer)
            "Connection to the bore field (pressure drop lumped in plant and network model)"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={-80,-90})));
          Buildings.Experimental.DHC.Plants.Heating.SewageHeatRecovery pla(
            redeclare final package Medium = Medium,
            final mSew_flow_nominal=datDes.mPla_flow_nominal,
            final mDis_flow_nominal=datDes.mPla_flow_nominal,
            final dpSew_nominal=datDes.dpPla_nominal,
            final dpDis_nominal=datDes.dpPla_nominal,
            final epsHex=datDes.epsPla) "Sewage heat recovery plant"
            annotation (Placement(transformation(extent={{-160,-10},{-140,10}})));
          Buildings.Experimental.DHC.Networks.Combined.UnidirectionalSeries dis(
            redeclare final package Medium = Medium,
            final nCon=nBui,
            show_TOut=true,
            final mDis_flow_nominal=datDes.mPipDis_flow_nominal,
            final mCon_flow_nominal=datDes.mCon_flow_nominal,
            final dp_length_nominal=datDes.dp_length_nominal,
            final lDis=datDes.lDis,
            final lCon=datDes.lCon,
            final lEnd=datDes.lEnd,
            final allowFlowReversal=allowFlowReversalSer)
            "Distribution network"
            annotation (Placement(transformation(extent={{-20,130},{20,150}})));
          replaceable
            Loads.Combined.BaseClasses.PartialBuildingWithETS_chiller   bui[nBui]
            constrainedby
            Buildings.Experimental.DHC.Loads.Combined.BaseClasses.PartialBuildingWithETS(
            bui(each final facMul=facMul),
            redeclare each final package MediumBui = Medium,
            redeclare each final package MediumSer = Medium,
            each final allowFlowReversalBui=allowFlowReversalBui,
            each final allowFlowReversalSer=allowFlowReversalSer)
            "Building and ETS"
            annotation (Placement(transformation(extent={{-10,170},{10,190}})));
          Buildings.Fluid.Sensors.TemperatureTwoPort TDisWatSup(redeclare
              final package Medium = Medium, final m_flow_nominal=datDes.mPumDis_flow_nominal)
            "District water supply temperature" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={-80,20})));
          Buildings.Fluid.Sensors.TemperatureTwoPort TDisWatRet(redeclare
              final package Medium = Medium, final m_flow_nominal=datDes.mPumDis_flow_nominal)
            "District water return temperature" annotation (Placement(
                transformation(
                extent={{10,-10},{-10,10}},
                rotation=90,
                origin={80,0})));
          Buildings.Fluid.Sensors.TemperatureTwoPort TDisWatBorLvg(redeclare
              final package Medium = Medium, final m_flow_nominal=datDes.mPumDis_flow_nominal)
            "District water borefield leaving temperature" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={-80,-40})));
          Modelica.Blocks.Sources.Constant TSewWat(k=273.15 + 17)
            "Sewage water temperature"
            annotation (Placement(transformation(extent={{-280,30},{-260,50}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant THeaWatSupMaxSet[nBui](
            k=bui.THeaWatSup_nominal)
            "Heating water supply temperature set point - Maximum value"
            annotation (Placement(transformation(extent={{-250,210},{-230,230}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant TChiWatSupSet[nBui](
            k=bui.TChiWatSup_nominal)
            "Chilled water supply temperature set point"
            annotation (Placement(transformation(extent={{-220,190},{-200,210}})));
         Buildings.Controls.OBC.CDL.Reals.Sources.Constant THeaWatSupMinSet[nBui](
            each k=28 + 273.15)
            "Heating water supply temperature set point - Minimum value"
            annotation (Placement(transformation(extent={{-280,230},{-260,250}})));
         Buildings.Controls.OBC.CDL.Reals.MultiSum PPumETS(
            final nin=nBui)
            "ETS pump power"
            annotation (Placement(transformation(extent={{140,190},{160,210}})));
          Modelica.Blocks.Continuous.Integrator EPumETS(
            initType=Modelica.Blocks.Types.Init.InitialState)
            "ETS pump electric energy"
            annotation (Placement(transformation(extent={{220,190},{240,210}})));
          Modelica.Blocks.Continuous.Integrator EPumDis(
            initType=Modelica.Blocks.Types.Init.InitialState)
            "Distribution pump electric energy"
            annotation (Placement(transformation(extent={{220,-90},{240,-70}})));
          Modelica.Blocks.Continuous.Integrator EPumSto(
            initType=Modelica.Blocks.Types.Init.InitialState)
            "Storage pump electric energy"
            annotation (Placement(transformation(extent={{220,-150},{240,-130}})));
          Modelica.Blocks.Continuous.Integrator EPumPla(initType=Modelica.Blocks.Types.Init.InitialState)
            "Plant pump electric energy"
            annotation (Placement(transformation(extent={{220,30},{240,50}})));
         Buildings.Controls.OBC.CDL.Reals.MultiSum EPum(nin=4)
            "Total pump electric energy"
            annotation (Placement(transformation(extent={{280,110},{300,130}})));
         Buildings.Controls.OBC.CDL.Reals.MultiSum PHeaPump(nin=nBui)
            "Heat pump power"
            annotation (Placement(transformation(extent={{140,150},{160,170}})));
          Modelica.Blocks.Continuous.Integrator EHeaPum(
            initType=Modelica.Blocks.Types.Init.InitialState)
            "Heat pump electric energy"
            annotation (Placement(transformation(extent={{220,150},{240,170}})));
         Buildings.Controls.OBC.CDL.Reals.MultiSum ETot(nin=2) "Total electric energy"
            annotation (Placement(transformation(extent={{320,150},{340,170}})));
          Buildings.Experimental.DHC.Loads.BaseClasses.ConstraintViolation conVio(
            final uMin(
              final unit="K",
              displayUnit="degC") = datDes.TLooMin,
            final uMax(
              final unit="K",
              displayUnit="degC") = datDes.TLooMax,
            final nu=3 + nBui,
            u(each final unit="K", each displayUnit="degC"))
            "Check if loop temperatures are within given range"
            annotation (Placement(transformation(extent={{320,10},{340,30}})));
        equation
          connect(dis.TOut, conVio.u[4:4+nBui-1]);
          connect(bou.ports[1], pumDis.port_a)
            annotation (Line(points={{102,-20},{80,-20},{80,-50}}, color={0,127,255}));
          connect(borFie.port_b, conSto.port_aCon) annotation (Line(points={{-120,-80},
                  {-100,-80},{-100,-84},{-90,-84}}, color={0,127,255}));
          connect(pumDis.port_b, conSto.port_aDis) annotation (Line(points={{80,-70},{
                  80,-120},{-80,-120},{-80,-100}}, color={0,127,255}));
          connect(borFie.port_a, pumSto.port_b)
            annotation (Line(points={{-140,-80},{-170,-80}}, color={0,127,255}));
          connect(conSto.port_bCon, pumSto.port_a) annotation (Line(points={{-90,-90},{
                  -100,-90},{-100,-100},{-200,-100},{-200,-80},{-190,-80}}, color={0,
                  127,255}));
          connect(conPla.port_bDis, TDisWatSup.port_a)
            annotation (Line(points={{-80,0},{-80,10}}, color={0,127,255}));
          connect(TDisWatSup.port_b, dis.port_aDisSup) annotation (Line(points={{-80,30},
                  {-80,140},{-20,140}}, color={0,127,255}));
          connect(dis.port_bDisSup, TDisWatRet.port_a)
            annotation (Line(points={{20,140},{80,140},{80,10}}, color={0,127,255}));
          connect(TDisWatRet.port_b, pumDis.port_a)
            annotation (Line(points={{80,-10},{80,-50}}, color={0,127,255}));
          connect(conSto.port_bDis, TDisWatBorLvg.port_a)
            annotation (Line(points={{-80,-80},{-80,-50}}, color={0,127,255}));
          connect(TDisWatBorLvg.port_b, conPla.port_aDis)
            annotation (Line(points={{-80,-30},{-80,-20}}, color={0,127,255}));
          connect(bui.port_bSerAmb, dis.ports_aCon) annotation (Line(points={{10,180},{20,
                  180},{20,160},{12,160},{12,150}}, color={0,127,255}));
          connect(dis.ports_bCon, bui.port_aSerAmb) annotation (Line(points={{-12,150},{
                  -12,160},{-20,160},{-20,180},{-10,180}}, color={0,127,255}));
          connect(TSewWat.y, pla.TSewWat) annotation (Line(points={{-259,40},{
                  -180,40},{-180,7.33333},{-161.333,7.33333}},
                                      color={0,0,127}));
          connect(THeaWatSupMaxSet.y, bui.THeaWatSupMaxSet) annotation (Line(points={{-228,
                  220},{-20,220},{-20,187},{-12,187}}, color={0,0,127}));
          connect(TChiWatSupSet.y, bui.TChiWatSupSet) annotation (Line(points={{-198,200},
                  {-24,200},{-24,185},{-12,185}},      color={0,0,127}));
          connect(pla.port_bSerAmb, conPla.port_aCon) annotation (Line(points={{-140,1.33333},
                  {-100,1.33333},{-100,-4},{-90,-4}}, color={0,127,255}));
          connect(conPla.port_bCon, pla.port_aSerAmb) annotation (Line(points={{-90,-10},
                  {-100,-10},{-100,-20},{-200,-20},{-200,1.33333},{-160,1.33333}},
                color={0,127,255}));
          connect(THeaWatSupMinSet.y, bui.THeaWatSupMinSet) annotation (Line(points={{-258,
                  240},{-16,240},{-16,189},{-12,189}}, color={0,0,127}));
          connect(bui.PPumETS, PPumETS.u)
            annotation (Line(points={{7,192},{7,200},{138,200}}, color={0,0,127}));
          connect(PPumETS.y, EPumETS.u)
            annotation (Line(points={{162,200},{218,200}}, color={0,0,127}));
          connect(pumDis.P, EPumDis.u)
            annotation (Line(points={{71,-71},{71,-80},{218,-80}}, color={0,0,127}));
          connect(pumSto.P, EPumSto.u) annotation (Line(points={{-169,-71},{-160,-71},{-160,
                  -140},{218,-140}}, color={0,0,127}));
          connect(pla.PPum, EPumPla.u) annotation (Line(points={{-138.667,
                  5.33333},{-120,5.33333},{-120,40},{218,40}},
                                                     color={0,0,127}));
          connect(EPumETS.y, EPum.u[1]) annotation (Line(points={{241,200},{260,200},{260,
                  119.25},{278,119.25}},
                                       color={0,0,127}));
          connect(EPumPla.y, EPum.u[2]) annotation (Line(points={{241,40},{260,40},{260,
                  119.75},{278,119.75}},
                                       color={0,0,127}));
          connect(EPumDis.y, EPum.u[3]) annotation (Line(points={{241,-80},{262,-80},{262,
                  120.25},{278,120.25}},
                                       color={0,0,127}));
          connect(EPumSto.y, EPum.u[4]) annotation (Line(points={{241,-140},{264,-140},{
                  264,120.75},{278,120.75}},
                                           color={0,0,127}));
          connect(PHeaPump.y, EHeaPum.u)
            annotation (Line(points={{162,160},{218,160}}, color={0,0,127}));
          connect(EHeaPum.y, ETot.u[1]) annotation (Line(points={{241,160},{300,160},{300,
                  159.5},{318,159.5}}, color={0,0,127}));
          connect(EPum.y, ETot.u[2]) annotation (Line(points={{302,120},{310,120},{310,160.5},
                  {318,160.5}},    color={0,0,127}));
          connect(TDisWatSup.T, conVio.u[1]) annotation (Line(points={{-91,20},{-100,20},
                  {-100,12},{-60,12},{-60,20},{318,20}},           color={0,0,127}));
          connect(TDisWatBorLvg.T, conVio.u[2]) annotation (Line(points={{-91,-40},{-100,
                  -40},{-100,-30},{-60,-30},{-60,-40},{300,-40},{300,20},{318,20}},
                                                                color={0,0,127}));
          connect(TDisWatRet.T, conVio.u[3]) annotation (Line(points={{69,6.66134e-16},{
                  60,6.66134e-16},{60,20},{318,20}},           color={0,0,127}));
          connect(bui.PCoo, PHeaPump.u) annotation (Line(points={{12,187},{132,
                  187},{132,160},{138,160}}, color={0,0,127}));
          annotation (Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-360,-260},{360,260}})),
              Documentation(revisions="<html>
<ul>
<li>
June 2, 2023, by Michael Wetter:<br/>
Added units to <code>conVio</code>.
</li>
<li>
November 16, 2022, by Michael Wetter:<br/>
Set correct nominal pressure for distribution pump.
</li>
<li>
February 23, 2021, by Antoine Gautier:<br/>
Refactored with base classes from the <code>DHC</code> package.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1769\">
issue 1769</a>.
</li>
<li>
January 16, 2020, by Michael Wetter:<br/>
Added documentation.
</li>
</ul>
</html>",         info="<html>
<p>
Partial model that is used by the reservoir network models.
The reservoir network models extend this model, add controls,
and configure some component sizes.
</p>
</html>"));
        end PartialSeries_chiller;
      annotation (Documentation(info="<html>
<p>
This package contains base classes that are used to construct the classes in
<a href=\"modelica://Buildings.Experimental.DHC.Examples.Combined\">
Buildings.Experimental.DHC.Examples.Combined</a>.
</p>
</html>"));
      end BaseClasses;
      annotation (
        preferredView="info",
        Documentation(
          info="<html>
<p>
This package contains example models for
district heating and cooling systems.
</p>
</html>"));
    end Combined;

    annotation (
      preferredView="info",
      Documentation(
        info="<html>
<p>
This package contains district heating and cooling case studies to show how the 
developed models can be used for design and operation.
</p>
</html>"));
  end Examples;

  annotation (
    preferredView="info",
    Documentation(
      info="<html>
<p>
This package contains models for district heating and cooling (DHC) systems.
</p>
</html>"));
end DHC;
