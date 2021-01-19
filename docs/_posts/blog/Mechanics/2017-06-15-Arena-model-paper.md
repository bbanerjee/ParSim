---
layout: posts
title:  "The ARENA model for partially saturated soils"
subheadline: "Biswajit Banerjee"
description: "A constitutive model for high-rate loading of soils"
date:  2017-06-15 10:30:00
categories:
    - Mechanics
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
I developed the ARENA model for partially saturated soils last year (2016).  Before that
we had tried using a typical Drucker-Prager with cap model and then a modified Cam-Clay
model but could not get these models to either represent observed experimental data or
produce robust solutions under large deformations.  The new model, though not perfect,
produced unexpectedly good predictions for partially saturated soils after being
calibrated using only dry soil data.

The model and its predictions have finally been made publicly available (with a large
amount of detail).  You can download the entire 15Mb report at [ResearchGate](https://www.researchgate.net/publication/317578167_Theory_verification_and_validation_of_the_ARENA_constitutive_model_for_applications_to_high-rate_loading_of_fully_or_partially_saturated_granular_media).

The citation for the report is:

> B. Banerjee and R. M. Brannon, 2017, "Theory, verification, and validation of the ARENA constitutive model for applications to high-rate loading of fully or partially saturated granular media", Technical Report #PAR-10021867-1516.v1, Parresia Research Limited, DOI: 10.13140/RG.2.2.10671.53922.

We hope to publish a much abridged version of this report as a book chapter soon.

#### Some predictions from ARENA ####
One of the striking features of ARENA is that it can predict the behavior of *partially saturated*
soils under high-rate compression loading even though the model parameters are calibrated using *dry* soil properties.  You can see a predicted curve compared with experimental data in the figure below.  The material is Colorado Mason sand containing 18% water by weight and the test is a split-Hopkinson pressure bar experiment that produces an average strain-rate of around 350/s.

<img style="width:500px" alt="ARENA_prediction" src="{{site.url}}/assets/blogimg/MasonSandUniaxialStrainSHPB-081612-003_Sig_t.png"/>

A frame from a typical low-resolution MPM simulation of an explosion in ARENA soil can be seen
in the image below.

<img style="width:600px" alt="ARENA_prediction" src="{{site.url}}/assets/blogimg/Centrifuge_VHull_BoulderClay_20g_13ww_midPBC_00120.png"/>


#### Remarks ####
Though the ARENA model produces excellent predictions in compression, its behavior in tension (particularly under disaggregation conditions) is not very accurate.  Our report contains a discussion section that identifies numerous research questions that arose during our work.  I hope some of you will find those interesting enough to improve upon our model.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

