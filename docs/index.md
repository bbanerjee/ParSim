---
title: Parresia Engineering Simulations
layout: splash
author_profile: true
header:
  overlay_color: "#000"
  overlay_filter: "0.5"
  overlay_image: /assets/img/mpm/foam_expansion.png
  actions:
    - label: "Download Brochure"
      url: "assets/brochure/ParresiaBrochure.pdf"
    - label: "Download ParSim"
      url: "https://github.com/bbanerjee/ParSim"
  caption: "Photo: Generation of a metal foam with Vaango &copy; Biswajit Banerjee, 2012"
excerpt: "We use open-source tools to bring parallel multiphysics engineering simulations to small and medium businesses.  Here you will find technical papers, expository blogs, and source code used in our work."

intro:
  - excerpt: '<strong>Technical reports and papers</strong>:  You can download draft reports of some of our work on various aspects of computational engineering and simulation from the links on this page.  More information on the types of analysis that we do can be found in our brochure. <strong>Contact</strong>: Dr. Biswajit Banerjee, b.banerjee.nz@gmail.com'

feature_row1:
  - image_path: assets/img/reports/mpm_simulations_jan_2021.png
    alt: "comparison of uxo penetration simulation"
    title: "Material point method simulations of penetration of sand modeled with tabular, support-vector, and neural network models of elastoplasticity"
    excerpt: "Penetration of high-speed projectiles into soils and the effect of explosions in soils on objects above the surface have interested civil and military engineers for decades. With the advent of faster computers and better numerical algorithms, many intractable problems in these domains have become possible to solve. In this paper, we use the Material Point Method (MPM), and a tabular elastoplastic model for soils, to simulate penetration in dry sand. ..."
    url: "assets/tech_reports/MPMSimulations_Jan_2021.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row2:
  - image_path: assets/img/reports/model_compare_dec_2020.png
    alt: "comparison of singe particle tests"
    title: "Comparison of tabular, support-vector, and neural network models for granular elastoplasticity"
    excerpt: "Experiments and micromechanical simulations generate tabular data for material behavior.  Typically, models are fit to these material data before engineering simulations can be performed.  In this paper, we compare the response of elastic moduli models for a dry, poorly-graded, sand that use linear interpolation, support vector regression fits, and multilayer perceptron neural networks, respectively. ..."
    url: "assets/tech_reports/ModelCompare_Dec_2020.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row3:
  - image_path: assets/img/reports/tabular_model_dec_2020.png
    alt: "tabular elastoplasticity"
    title: "Tabular models for high pressure and high strain-rate granular plasticity"
    excerpt: "Experiments and micromechanical simulations generate tabular data for material behavior.  Typically, models are fit to these material data before engineering simulations can be performed.  It is frequently discovered that existing models cannot express the experimental data adequately, and new models have to be developed and fit. This process is undesirable and a preferable approach is to directly use the tabular data without model building. In this work, we discuss such a tabular model and associated numerical algorithms in the context of the elastoplastic behavior of a poorly graded concrete sand ..."
    url: "assets/tech_reports/TabularModel_Dec_2020.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row4:
  - image_path: assets/img/reports/tabular_interp_nov_2020.png
    alt: "tabular interpolation"
    title: "Interpolating tabular data for granular material models"
    excerpt: "Experimental and microscale simulation data are used to design material models for granular materials. These data are collected in tabular form, often as the outcome of a design-of-experiments process when multiple independent variables are expected to affect the result of an experiment. Tabular data are collected densely for one independent variable, typically the strain. Data are obtained sparsely in the other dimensions. In this paper, we discuss possible approaches to using tabular data directly in material models without attempting to design and fit closed-form expressions...."
    url: "assets/tech_reports/TabularInterp_Nov_2020.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row5:
  - image_path: assets/img/reports/neural_net_sep_2020.png
    alt: "neural network model"
    title: "Multilayer perceptron neural networks as multi-variable material models"
    excerpt: "Material models that depend on multiple independent variables are often necessary for accurate numerical simulations, particularly for applications that involve large stresses and deformations. It is rare that purely physics-based models are used in simulations because of the attendant computational cost. Instead, experimental and microscale simulation data are expressed as phenomenological models and fed into simulations. As the number of independent variables increases, such models are not only difficult to design but also need exponentially larger amounts of data to parameterize accurately. In this paper we examine an alternative procedure for developing multi-variable phenomenological models via multi-layer perceptron neural networks..."
    url: "assets/tech_reports/Neural_net_Sep_2020.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row6:
  - image_path: assets/img/reports/support_vector_sep_2020.png
    alt: "support vetor model"
    title: "Support vector regression for fitting multi-variable material models"
    excerpt: "Analytical expressions for phenomenological material models that depend on multiple independent variables are notoriously difficult to design. Parameter determination is also intimately tied with the model design process. Soils that exhibit elastic-plastic coupling are particularly prone to the design problem. It is not uncommon to have to redesign models for every new soil that is characterized experimentally. An unstated assumption in soil mechanics is that small inaccuracies in material models do not affect the predictive capabilities of those models significantly. First, we demonstrate that such an assumption in not warranted, particularly in the large deformation, non-monotonic loading, regime. We then proceed to explore support vector regression to replace analytical models..."
    url: "assets/tech_reports/Support_vector_Sep_2020.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"

feature_row7:
  - image_path: assets/img/reports/arena_model_sep_2017.png
    alt: "ARENA material model"
    title: "Theory, verification, and validation of the ARENA constitutive model"
    excerpt: "The Arena continuum-scale model for sand and/or clay under high-rate loading conditions is presented. Our scope is limited to adiabatic load/unload conditions in order to focus on model features that most crucial for simulations of buried explosives and similar phenomena that involve shock compression followed by free expansion (possibly with re-compression when ejecta impacts an object).  Evidence is provided that such conditions fall in a realm for which there is no substantial difference between additive or multiplicative inelasticity approaches. The Arena model is implemented in a Material Point Method (MPM) code and details of the implementation and algorithms are discussed...."
    url: "assets/tech_reports/ArenaModel_Sep_2017.pdf"
    btn_label: "Download report"
    btn_class: "btn--primary"
---

{% include feature_row id="intro" type="left" %}

{% include feature_row id="feature_row1" type="left" %}

{% include feature_row id="feature_row2" type="left" %}

{% include feature_row id="feature_row3" type="left" %}

{% include feature_row id="feature_row4" type="left" %}

{% include feature_row id="feature_row5" type="left" %}
 
{% include feature_row id="feature_row6" type="left" %}

{% include feature_row id="feature_row7" type="left" %}






