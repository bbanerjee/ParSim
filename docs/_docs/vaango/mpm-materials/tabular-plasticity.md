---
title: Tabular plasticity 
permalink: /docs/vaango/mpm-materials/tabular-plasticity/
---

* Contents
{:toc}

The tabular plasticity model is designed to simulate hypoelastic-plasticity using 
tabulated data for elastic moduli and the yield function for $$J_2$$ and $$J_2-I_1$$ models
of perfect plasticity. 

The model has been designed with the high-rate deformation of soils in mind but can 
also be used for non-hardening metals. Sample input files for linear elastic materials
and von Mises plasticity can be found in the [repository](https://github.com/bbanerjee/ParSim/tree/master/Vaango/src/StandAlone/inputs/MPM/TabularModels/TabularPlasticity).

Since the tabular plasticity model was designed for materials that have almost no tensile strength,
the inputs are expected in the `compression positive`  convention.  Note that the general
convention used in the Vaango code is that `tension is positive and compression is negative` .
Conversions are done internally in the code to make sure that signs are consistent.

### Shear modulus model
  The shear modulus is either a constant ($$G_0$$) or determined using a Poisson's ratio ($$\nu$$)
  from the tabular bulk modulus, $$K(p)$$, assuming a linear elastic model:
  <div>
  $$
    G = \frac{3K(1-2\nu)}{2(1+\nu)}
  $$
  </div>
  This relation is activated if $$\nu \in [-1.0, 0.5)$$, otherwise the constant shear modulus
  is used.

### Bulk modulus model
![Bulk modulus model]({{site.url}}/assets/img/mpm/material-models/TabularHydrostat.png){:class="img-responsive" height="250px" border="5px double red" background-color="white"}

  The bulk modulus is determined from a table of unloading curves.  Each unloading 
  curve is associated with a Hencky plastic volumetric strain ($$\bar{\varepsilon_v^p}$$).
  These strains are associated in the JSON file with the key `PlasticStrainVol`
  and added as the first independent variable in the Vaango `ups` input file.

  For each plastic strain value, the JSON file has to contain an associated data
  set of `TotalStrainVol` (the total Hencky volumetric strain, $$\bar{\varepsilon_v}$$)
  and the `Pressure` (the mean stress, $$\bar{p}$$).

### Yield function
![Yield function]({{site.url}}/assets/img/mpm/material-models/TabularYieldFn.png){:class="img-responsive" height="250px" border="5px double red" background-color="white"}

  The tabular yield condition is
  <div>
  $$
     f = \sqrt{J_2} - g(\bar{p}) = 0
  $$
  </div>
  The function $$g(\bar{p})$$ is provided in tabular form.  Only one such function is
  allowed.  The independent variable in the associated JSON file will have the name
  `Pressure` while the dependent variable will have the name `SqrtJ2`.
### Input file formats
The inputs for the tabular plasticity model are typically specified as follows.
{% highlight xml %}
<?xml version='1.0' encoding='ISO-8859-1' ?>
<Uintah_specification>
  <Meta>
      <title>Tabular_Verification_Test_05_Uniaxial_Strain_Compression_DP_with_LoadUnload</title>
  </Meta>
  <SimulationComponent type="mpm" />
  <Time>
      <maxTime> 8.0 </maxTime>
      <initTime> 0.0 </initTime>
      <delt_min> 1.0e-8 </delt_min>
      <delt_max> 0.01 </delt_max>
      <timestep_multiplier> 0.3 </timestep_multiplier>
  </Time>
  <DataArchiver>
      <filebase>TabularTest_05_UniaxialStrainLoadUnloadNonLinDPNonLin.uda</filebase>
      <outputInterval>1.0e-3</outputInterval>
      <outputInitTimestep/>
      <save label = "p.x"/>
      <save label = "p.color"/>
      <save label = "p.temperature"/>
      <save label = "p.velocity"/>
      <save label = "p.particleID"/>
      <save label = "p.stress"/>
      <save label = "g.mass"/>
      <save label = "p.deformationGradient"/>
      <save label = "g.acceleration"/>
      <save label = "p.plasticVolStrain"/>
      <save label = "p.elasticVolStrain"/>
      <checkpoint cycle = "2" timestepInterval = "2000"/>
  </DataArchiver>
  <MPM>
    <time_integrator>              explicit   </time_integrator>
    <interpolator>                 linear     </interpolator>
    <use_load_curves>              false      </use_load_curves>
    <minimum_particle_mass>        1.0e-15    </minimum_particle_mass>
    <minimum_mass_for_acc>         1.0e-15    </minimum_mass_for_acc>
    <maximum_particle_velocity>    1.0e5      </maximum_particle_velocity>
    <artificial_damping_coeff>     0.0        </artificial_damping_coeff>
    <artificial_viscosity>         true       </artificial_viscosity>
    <artificial_viscosity_heating> false      </artificial_viscosity_heating>
    <do_contact_friction_heating>  false      </do_contact_friction_heating>
    <create_new_particles>         false      </create_new_particles>
    <use_momentum_form>            false      </use_momentum_form>
    <with_color>                   true       </with_color>
    <use_prescribed_deformation>   true       </use_prescribed_deformation>
    <prescribed_deformation_file>  TabularTest_05_PrescribedDeformation.inp   </prescribed_deformation_file>
    <erosion algorithm = "none"/>
  </MPM>
  <PhysicalConstants>
      <gravity>[0,0,0]</gravity>
  </PhysicalConstants>
  <MaterialProperties>
    <MPM>
      <material name="TabularPlastic">
        <density>1050</density>
        <melt_temp>3695.0</melt_temp>
        <room_temp>294.0</room_temp>
        <thermal_conductivity>174.0e-7</thermal_conductivity>
        <specific_heat>134.0e-8</specific_heat>
        <constitutive_model type="tabular_plasticity">
          <elastic_moduli_model type="tabular">
            <filename>TabularTest_05_Elastic.json</filename>
            <independent_variables>PlasticStrainVol, TotalStrainVol</independent_variables>
            <dependent_variables>Pressure</dependent_variables>
            <interpolation type="linear"/>
            <G0>3500</G0>
            <nu>0.35</nu>
          </elastic_moduli_model>
          <plastic_yield_condition type="tabular">
            <filename>TabularTest_05_Yield.json</filename>
            <independent_variables>Pressure</independent_variables>
            <dependent_variables>SqrtJ2</dependent_variables>
            <interpolation type="linear"/>
          </plastic_yield_condition>
        </constitutive_model> 
        <geom_object>
          <box label = "Plate1">
            <min>[0.0,0.0,0.0]</min>
            <max>[1.0,1.0,1.0]</max>
          </box>
          <res>[1,1,1]</res>
          <velocity>[0.0,0.0,0.0]</velocity>
          <temperature>294</temperature>
          <color>0</color>
        </geom_object>
      </material>
      <contact>
        <type>null</type>
        <materials>[0]</materials>
        <mu>0.1</mu>
      </contact>
    </MPM>
  </MaterialProperties>
  <Grid>
      <BoundaryConditions>                      
      </BoundaryConditions>
      <Level>
        <Box label = "1">
            <lower>[-2.0, -2.0, -2.0]</lower>
            <upper>[3.0, 3.0, 3.0]</upper>
            <resolution>[5,5,5]</resolution>
            <extraCells>[0,0,0]</extraCells>
            <patches>[1,1,1]</patches>
        </Box>
      </Level>
  </Grid>
</Uintah_specification>
{% endhighlight %}

#### The prescribed deformation file
The input file listed above is for a single particle driven by a prescribed deformation gradient.  The
deformation file is `TabularTest_05_PrescribedDeformation.inp` which contains the following:
{% highlight xml %}
0  1.0  0  0  0  1  0  0  0  1  0  0  0  0
1  3.0  0  0  0  1  0  0  0  1  0  0  0  0  
3  1.0  0  0  0  1  0  0  0  1  0  0  0  0
5  0.5  0  0  0  1  0  0  0  1  0  0  0  0
7  1.0  0  0  0  1  0  0  0  1  0  0  0  0
8  1.02 0  0  0  1  0  0  0  1  0  0  0  0 
{% endhighlight %}
The first column is the time (seconds), the second column is deformation gradient $$F_{xx}$$.  This
deformation represents uniaxial strain in the $$x$$-direction.

#### The bulk-modulus model
The bulk modulus model is non-linear ($$\tanh$$-form) and is input as a JSON file containing data
that covers the range expected during a simulation.  A sample file  
(`TabularTest_05_Elastic.json`) is shown below:
{% highlight json %}
{"Vaango_tabular_data": {
  "Meta" : {
    "title" : "Nonlinear elastic data"
  },
  "Data" : {
    "PlasticStrainVol" : [-2.5, 2.5],
    "Data" : [{
      "TotalStrainVol" : [-15.000000, -14.310345, -13.620690, -12.931034, -12.241379, -11.551724, -10.862069, -10.172414, -9.482759, -8.793103, -8.103448, -7.413793, -6.724138, -6.034483, -5.344828, -4.655172, -3.965517, -3.275862, -2.586207, -1.896552, -1.206897, -0.517241, 0.172414, 0.862069, 1.551724, 2.241379, 2.931034, 3.620690, 4.310345, 5.000000],
      "Pressure" : [-2959.842894, -2943.464146, -2920.494168, -2888.366995, -2843.600634, -2781.548344, -2696.156593, -2579.811517, -2423.426002, -2217.004505, -1950.975669, -1618.496214, -1218.556435, -759.014896, -257.981931, 257.981931, 759.014896, 1218.556435, 1618.496214, 1950.975669, 2217.004505, 2423.426002, 2579.811517, 2696.156593, 2781.548344, 2843.600634, 2888.366995, 2920.494168, 2943.464146, 2959.842894]
    }, {
      "TotalStrainVol" : [-5.000000, -4.310345, -3.620690, -2.931034, -2.241379, -1.551724, -0.862069, -0.172414, 0.517241, 1.206897, 1.896552, 2.586207, 3.275862, 3.965517, 4.655172, 5.344828, 6.034483, 6.724138, 7.413793, 8.103448, 8.793103, 9.482759, 10.172414, 10.862069, 11.551724, 12.241379, 12.931034, 13.620690, 14.310345, 15.000000],
      "Pressure" : [-2959.842894, -2943.464146, -2920.494168, -2888.366995, -2843.600634, -2781.548344, -2696.156593, -2579.811517, -2423.426002, -2217.004505, -1950.975669, -1618.496214, -1218.556435, -759.014896, -257.981931, 257.981931, 759.014896, 1218.556435, 1618.496214, 1950.975669, 2217.004505, 2423.426002, 2579.811517, 2696.156593, 2781.548344, 2843.600634, 2888.366995, 2920.494168, 2943.464146, 2959.842894]
    }]
  }
}}
{% endhighlight %}

#### The yield function
The yield function in this example is a nonlinear Drucker-Prager model.  The input file
(`TabularTest_05_Yield.json`) is shown below:
{% highlight json %}
{"Vaango_tabular_data": {
  "Meta" : {
    "title" : "Quadratic Drucker-Prager yield data"
  },
  "Data" : {
    "Pressure" : [-333.333333, -298.850575, -264.367816, -229.885057, -195.402299, -160.919540, -126.436782, -91.954023, -57.471264, -22.988506, 11.494253, 45.977011, 80.459770, 114.942529, 149.425287, 183.908046, 218.390805, 252.873563, 287.356322, 321.839080, 356.321839, 390.804598, 425.287356, 459.770115, 494.252874, 528.735632, 563.218391, 597.701149, 632.183908, 666.666667],
    "SqrtJ2" :   [0.000000, 45.485883, 64.326752, 78.783860, 90.971765, 101.709526, 111.417203, 120.344334, 128.653504, 136.457648, 143.838990, 150.859606, 157.567719, 164.001682, 170.192589, 176.166066,181.943530, 187.543098, 192.980256, 198.268366, 203.419051, 208.442500, 213.347701, 218.142630, 222.834406, 227.429413, 231.933403, 236.351579, 240.688667, 244.948974]
  }
}}
{% endhighlight %}

