---
layout: docs
title:  "Setting up the Three.js renderer"
subheadline: "Biswajit Banerjee"
description: "Javascript scientific visualization - Part 4"
date:  2017-02-27 10:30:00
categories:
    - Javascript
    - Typescript
    - Threejs
    - Vue
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

##### Introduction #####
The previous articles in this series were about:

1. [Part 1: Reading VTK format particles with Javascript in a browser](http://www.parresianz.com/javascript/typescript/vue/vuex/xml/reading-vtk-particles/)
2. [Part 2: Saving the read-in particle data in a Vuex store](http://www.parresianz.com/javascript/typescript/vue/vuex/vuex-store/)
3. [Part 3: Initialization of a store and the user interface](http://www.parresianz.com/javascript/typescript/threejs/vue/vtk-threejs/) before plotting the particles with Three.js.

In this article we will discuss the actual process of visualization of the particles.  
The particles we are interested in are [oriented ellipsoids](https://en.wikipedia.org/wiki/Ellipsoid).

##### The graphics panel #####
In the `main-panel` Vue component introduced in Part 3, we had declared a graphics panel
component called `three-graphics-panel`.  A simple Vue template for this component can be
written as:

{% highlight html %}
<template>
  <div id='three-graphics-container'>
    <div class="uk-card uk-card-default uk-card-large">
      <div class="uk-card-body">
        <!-- The renderer Vue component -->
        <three-renderer v-bind:size="{w:500, h:500}">
          <!-- The scene Vue component -->
          <three-scene v-bind:size="size"> 
            <!-- The camera Vue component -->
            <three-camera  v-bind:size="size"
                           v-bind:position="{x: 100,  z: 15 }">
            </three-camera>
            <!-- The particles Vue component -->
            <three-ellipsoid-particles> </three-ellipsoid-particles>
          </three-scene>
        </three-renderer>
      </div>
    </div>
  </div>
</template>
<script src="./ThreeGraphicsPanel.ts"> </script>
{% endhighlight %}

As before, the CSS classes used in this example are from Uikit3.  The graphics panel uses
four child Vue components, `three-renderer`, `three-scene`, `three-camera`, and
`three-ellipsoid-particles`.  The implementation of `ThreeGraphicsPanel.ts` is also
straightforward in this case:

{% highlight js %}
import * as Vue from "vue";
import {Component} from "av-ts";
@Component({
  name: 'three-graphics-panel'
})
export default class ThreeGraphicsPanel extends Vue {
  name: string = "ThreeGraphicsPanel";
  size: Object = {
    w: 200,
    h: 200
  };
}
{% endhighlight %}

The `size` property is passed on to the children of three-graphics-panel and determines the
default size of the window in which the particles will be displayed.  This property can be
accessed by the children by using the
["passing data with props approach."](https://vuejs.org/v2/guide/components.html#Passing-Data-with-Props)
{:.notice}

##### The renderer template #####
The first component instantiated in `three-graphics-panel` is the `three-renderer` component which
is passed the `size` argument using the `v-bind` directive.  Let use first examine
the `ThreeRenderer.vue` code:

{% highlight html %}
<template>
  <div id="three-renderer-div">
    <slot></slot>
    <div ref="three-graphics-container"
         id="three-renderer-graphics-div">
    </div>
  </div>
</template>
<script src="./ThreeRenderer.ts"> </script>
{% endhighlight %}

The first new aspect to notice is the use of the
[`slot` tag](https://vuejs.org/v2/guide/components.html#Content-Distribution-with-Slots).  This
is needed to make sure that the `three-scene`, `three-camera`, and `three-ellipsoid-particles`
from the `three-graphics-panel` parent components are correctly instantiated and interpolated
into the DOM.
{:.notice}

The second new aspect is the use of the [`ref` attribute](https://vuejs.org/v2/api/#ref)
in the `div` element after the `slot`.  This element is used to create the HTML canvas as
a child of the renderer component.
{:.notice}

We could add event handlers to this template to handle interactions with the camera.  That
aspect of the user interface will be discussed in a future article.
{:.notice--info}

##### The renderer implementation #####
The renderer component is implemented in `ThreeRenderer.ts` which contains the following
code:

{% highlight js %}
import * as Vue from "vue";
import THREE = require('three');
import { Component, Lifecycle, Prop, p } from 'av-ts';
import Store from '../vuex/Store';
interface ScreenSize {
  left : number, top : number, width : number, height : number
}
@Component({
  name: 'ThreeRenderer'
})
export default class ThreeRenderer extends Vue {
  // These are the "prop" variables used by Vue
  @Prop
  size = p({
    type: Object, // { w, h }
    required: true
  });
  @Prop
  renderer = p({
    type: THREE.WebGLRenderer
  });
  // These are the "data" variables used by Vue
  public  d_renderer: THREE.WebGLRenderer;
  private d_size: ScreenSize;
  // These are the computed properties
  get getScene() {
    return Store.getters.scene;
  }
  get getCamera() {
    return Store.getters.camera;
  }
  // These are the lifecycle hooks
  @Lifecycle
  created() {
    this.d_renderer = this.renderer;
    if (!(this.d_renderer instanceof THREE.WebGLRenderer)) {
      this.d_renderer = new THREE.WebGLRenderer({ antialias: true });
    }
    let w = (<any>this.size).w;
    let h = (<any>this.size).h;
    this.d_size = {left: 0, top: 0, width: w, height: w};
    this.d_renderer.setSize(w, h);
    this.d_renderer.setClearColor(0x52576e)
  }
  @Lifecycle
  mounted() {
    this.d_renderer.domElement.setAttribute("id", "three-graphics-canvas");
    if ((this.$refs)["three-graphics-container"]) {
      let el  = (this.$refs)["three-graphics-container"];
      (<any>el).appendChild(this.d_renderer.domElement);

      // Update based on actual screen size
      this.setActualScreenSize(this.d_renderer.domElement);
      this.d_renderer.setSize(this.d_size.width, this.d_size.height);
      let camera = Store.getters.camera;
      camera.aspect = this.d_size.width / this.d_size.height;
      camera.updateProjectionMatrix();
      Store.commit('SET_CAMERA', camera);
    }
    this.animate();
  }
  @Lifecycle
  beforeDestroy() {
    Store.commit('DELETE_SCENE');
    Store.commit('DELETE_CAMERA');
  }
  // If there is a resize, compute the new screen size
  private setActualScreenSize(domElement: HTMLCanvasElement) : void {
    let box = domElement.getBoundingClientRect();
    let d = domElement.ownerDocument.docsumentElement;
    this.d_size.left = box.left + window.docsXOffset - d.clientLeft;
    this.d_size.top = box.top + window.docsYOffset - d.clientTop;
    this.d_size.width = box.width;
    this.d_size.height = box.height;
  }
  // Start the animation loop
  private animate() {
    requestAnimationFrame( this.animate );
    this.render();
  }
  // Do the actual rendering using the THREE.js render method
  private render() {
    this.d_renderer.render(Store.getters.scene, Store.getters.camera);
  }
}
{% endhighlight %}

The main thing to note here is that we are using the THREE.WebGLRenderer, and 
the creation of the `canvas` element as a child of the renderer div element.  
The canvas element is created in the `mount` stage of the instance lifecycle and
all the plotting is done on that element.
{:.notice}

#### Remarks ####
In the next part of this article we'll discuss more details of setting up the scene and camera.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

