---
title:  "Setting up the Three.js scene and camera"
subheadline: "Biswajit Banerjee"
description: "Javascript scientific visualization - Part 5"
date:  2017-02-28 10:30:00
categories:
    - Javascript
    - Typescript
    - Threejs
    - Vue
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---


##### Introduction #####
In the previous articles in this series we talked about:

1. [Part 1: Reading VTK format particles with Javascript in a browser]({{ site.baseurl }}/javascript/typescript/vue/vuex/xml/reading-vtk-particles/)
2. [Part 2: Saving the read-in particle data in a Vuex store]({{ site.baseurl }}/javascript/typescript/vue/vuex/vuex-store/)
3. [Part 3: Initialization of a store and the user interface]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs/).
4. [Part 4: Setting up the Three.js renderer]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs-renderer/)
<!--more-->

Let us now discuss how we can set up the scene and the camera.

##### Flash-back: The graphics panel #####
Recall that the `three-graphics-panel` template introduced the components `three-scene`
and `three-camera` in the file `ThreeGraphicsPanel.vue`.

{% highlight html %}
<template>
  <div id='three-graphics-container'>
    <div class="uk-card uk-card-default uk-card-large">
      <div class="uk-card-body">
        <three-renderer v-bind:size="{w:500, h:500}">
          <three-scene v-bind:size="size"> 
            <three-camera v-bind:size="size" v-bind:position="{x: 100,  z: 15 }">
            </three-camera>
            <three-ellipsoid-particles> </three-ellipsoid-particles>
          </three-scene>
        </three-renderer>
      </div>
    </div>
  </div>
</template>
<script src="./ThreeGraphicsPanel.ts"> </script>
{% endhighlight %}

The `three-scene` component is a child of `three-renderer` and has two children
of its own: `three-camera` and `three-ellipsoid-particles`.
{:.notice}

The `size` parameter is passed down from `three-renderer` to its children.
{:.notice}

The `three-camera` component takes an additional `position` property as input.
{:.notice}

##### The `three-scene` component #####
The `three-scene` component has a slot for the `three-camera` and
`three-ellipsoid-particles` and the corresponding 'ThreeScene.vue` file contains:

{% highlight html %}
<template>
    <div id="three-scene-div">
        <slot></slot>
    </div>
</template>
<script src="./ThreeScene.ts"> </script>
{% endhighlight %}

Note that we could have included the scene in `three-renderer` instead of creating a
new component for it.  That choice depends on the use-case, e.g., if there are
multiple scenes assigned to a single renderer and we would like to switch between them
then it makes sense to make the scene into a component.
{:.notice--info}

The scene component is implemented in `ThreeScene.ts`:

{% highlight js %}
import * as Vue from "vue";
import THREE = require("three");
import { Component, Lifecycle, Prop, p } from 'av-ts';
import Store from './Store';
@Component({
  name: 'ThreeScene'
})
export default class ThreeRenderer extends Vue {
  @Prop
  size = p({
    type: Object, // { w, h }
    required: true
  });
  @Lifecycle
  created() {
    // Get the scene from the store
    let scene = Store.getters.scene;
    // If a scene doesn't exist, create  one
    if (!scene) {
      // Create the scene
      scene = new THREE.Scene();
      // Add the lights
      let dirLight = new THREE.DirectionalLight( 0xffffff );
      dirLight.position.set( 1, 1, 1 );
      scene.add( dirLight );
      dirLight = new THREE.DirectionalLight( 0x002288 );
      dirLight.position.set( -1, -1, -1 );
      scene.add( dirLight );
      let ambLight = new THREE.AmbientLight( 0x222222 );
      scene.add( ambLight );
      // Update the scene
      Store.commit('SET_SCENE', scene);
    }
  }
  @Lifecycle
  mounted() {
    // Get the scene and the camera
    let scene = Store.getters.scene;
    let camera = Store.getters.camera;
    if (camera) {
      // Add the camera to the scene
      scene.add(camera);
      // Update the scene
      Store.commit('SET_SCENE', scene);
    } else {
      console.log("Camera has not been created yet");
    }
  }
}
{% endhighlight %}

During the creation phase, we add lights to the scene but not the camera(s).  This is
because the camera component has not been instantiated yet.  The camera is added to
the scene during the mount stage.
{:.notice}

We don't use the `size` property here, but keep it around in case it's needed by
child components.
{:.notice--info}

##### The `three-camera` component #####
Now we are ready to create a camera that can be used to view the scene.  The template
file `ThreeCamera.vue` is essentially empty:

{% highlight html %}
<template>
</template>
<script src="./ThreeCamera.ts"> </script>
{% endhighlight %}

The empty `template` tags look redundant, but if we don't include them in `ThreeCamera.vue`, the `ts-loader` fails with a compilation error.
{:.notice--warning}

The implementation of the camera component in `ThreeCamera.ts` contains:

{% highlight js %}
import * as Vue from "vue";
import THREE = require("three");
import { Component, Lifecycle, Prop, p } from 'av-ts';
import Store from "./Store";
import assign from './util';
@Component({
  name: 'Camera',
})
export default class Camera extends Vue {
  // "props"
  @Prop
  size = p({
    type: Object, // { w, h }
  });
  @Prop
  position = p({
    type: Object  // { x, y, z }
  })  
  // "data"
  d_size: any;
  d_position: any;
  @Lifecycle
  created() {
    // Keep local copies of the properties in case we need to modify them
    this.d_size = this.size;
    this.d_position = this.position;
    let camera = Store.getters.camera;
    if (!camera) {
      let w = (<any>this.size).w;
      let h = (<any>this.size).h;
      camera = new THREE.PerspectiveCamera(75, w / h, 0.1, 1000);
      // Update the position if the template specifies that
      if (this.position) {
        this.updateCameraPosition(camera, this.position);
        camera.updateProjectionMatrix(); // Don't forget this update
      }
      Store.commit('SET_CAMERA', camera);
    }
  }
  private updateCameraPosition(camera: THREE.Camera, pos: any) : void {
    console.log(pos);
    assign(camera.position, pos);
    camera.lookAt(new THREE.Vector3(0.0, 0.0, 0.0));
  }
}
{% endhighlight %}

We use a [`THREE.PerspectiveCamera`](https://threejs.org/docs/api/cameras/PerspectiveCamera.html)
in this example.  The four parameters are the field of view angle in degrees, the
aspect ratio, the near plane distance, and the far plane distance of the camera frustum.
The camera position is determined by the input parameters and not by the positions of
the objects in the scene.
{: .notice}

The aspect ratio is determined at this stage by the input parameters and not by the actual
window dimensions.  Also, I recommended that you use the `lookAt` method to make sure
that the camera is oriented as you would expect.
{: .notice--info}

Don't forget to update the projection matrix of the `PerspectiveCamera` after you change
its position.
{: .notice--warning}

#### Remarks ####
Now that the scene and the camera have been set up, we can see how ellipsoid are created and displayed in the next part of this series.


