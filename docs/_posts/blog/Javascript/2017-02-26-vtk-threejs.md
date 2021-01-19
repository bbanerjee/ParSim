---
layout: posts
title:  "Plotting VTK particles with Three.js"
subheadline: "Biswajit Banerjee"
description: "Javascript scientific visualization - Part 3"
date:  2017-02-26 10:30:00
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

##### Introduction #####
In [Part 2](http://www.parresianz.com/javascript/typescript/vue/vuex/vuex-store/)
of this series, we showed how the particle data are saved into a Vuex store.  Now we are ready
to visualize the data.  In this article, we will discuss how [three.js](https://threejs.org/)
can be used to display the particles. Interaction with the particles will be discussed in
a future article.

##### A Vuex store for graphics #####
To make our manipulations of the data more convenient, we use a store for the graphics state too.
The setup is similar to that discussed in Part 2, and we have a state, getters, mutations, and
actions.  In this case, we implement the state in a file called `ThreeGraphicsState.ts`:

{% highlight js %}
import THREE = require('three');
export class ThreeGraphicsState {
  scene      : THREE.Scene;
  camera     : THREE.Camera;
  threeObjectsCreated : boolean;
  constructor(w: number, h: number) {
    this.scene = <THREE.Scene>null;
    this.camera = <THREE.Camera>null;
    this.threeObjectsCreated = false;
  }
};
{% endhighlight %}

The state contains the scene, one camera (multiple cameras can be added if necessary), and
a flag to indicate whether the store contains objects or not.  The associated getters are
straightforward, but there are a few more mutation functions in this case:

{% highlight js %}
import THREE = require('three');
import {Mutation, MutationTree} from 'vuex';
import {ThreeGraphicsState} from './ThreeGraphicsState';
export function SET_SCENE(state: ThreeGraphicsState, scene: THREE.Scene) {
  state.scene = scene;
  console.log("set scene");
}
export function SET_CAMERA(state: ThreeGraphicsState, camera: THREE.Camera) {
  state.camera = camera;
  console.log("set camera");
}
export function DELETE_SCENE(state: ThreeGraphicsState, message: string) {
  state.scene = null;
  console.log("deleted scene" + message);
}
export function DELETE_CAMERA(state: ThreeGraphicsState, message: string) {
  state.camera = null;
  console.log("deleted camera" + message);
}
export function ADD_THREE_OBJECT(state: ThreeGraphicsState, object: any) {
  state.scene.add(object);
}
export function THREE_OBJECTS_CREATED(state: ThreeGraphicsState, value: boolean) {
  state.threeObjectsCreated = value;
}
export default <MutationTree<ThreeGraphicsState>> {
  SET_SCENE,
  SET_CAMERA,
  DELETE_SCENE,
  DELETE_CAMERA,
  ADD_THREE_OBJECT,
  THREE_OBJECTS_CREATED,
}
{% endhighlight %}

We can also add functions that can delete objects if needed.

##### Interacting with Vue #####
[Vue2's](https://vuejs.org/) capabilities simplify the development of the user interface, but a
significant penalty has to be paid in terms of performance.  In this section we will
take a glimpse at how Vue makes development of the user interface easier.

Vue2 was designed for business applications and is not the optimal choice for graphics-heavy
applications.
{:.notice}

In the `Main.ts` file that acts as the entry point, we declare the following:

{% highlight js %}
import * as Vue from "vue";
// Include the .vue files that contains the templates, scripts, and styles for each Vue component
var MainPanel               = require("./MainPanel.vue").default;
var ThreeGraphicsPanel      = require("./ThreeGraphicsPanel.vue").default;
var ThreeRenderer           = require("./ThreeRenderer.vue").default;
var ThreeScene              = require("./ThreeScene.vue").default;
var ThreeCamera             = require("./ThreeCamera.vue").default;
var ThreeEllipsoidParticles = require("./ThreeEllipsoidParticles.vue").default;
// Register the three.js graphics Vue components
Vue.component('three-graphics-panel',      ThreeGraphicsPanel);
Vue.component('three-renderer',            ThreeRenderer);
Vue.component('three-scene',               ThreeScene);
Vue.component('three-camera',              ThreeCamera);
Vue.component('three-ellipsoid-particles', ThreeEllipsoidParticles);
class Main {
  public vm : Vue;
  constructor() {
    this.vm = new Vue ({
        store: Store,
        el: '#main-panel',
        components: {
            'main-panel' : MainPanel
        },
        render: h => {
            return h('main-panel');
        }
    });
  }
}
{% endhighlight %}

In the `MainPanel.vue` file, we set up the basic components used in the user interface; in
this case just a single window for display the particle data.  I have used
[Uikit3](https://getuikit.com/) for some of the CSS styling that I refer to in this series, mainly
because it is _relatively_ lightweight.

{% highlight html %}
<template>
  <div id='vaango-main-panel' 
       class="uk-container uk-container-center uk-margin-remove-top uk-margin-remove-left 
              uk-margin-remove-right uk-margin-large-bottom"> 
    <div id="container" class="uk-flex uk-flex-around">
      <!-- The THREE graphics Vue component -->
      <three-graphics-panel class="uk-width-xxlarge">
      </three-graphics-panel>
    </div>
  </div>
</template>
<script src="./MainPanel.ts"> </script>
<style>
 .uk-container {
   padding-left: 5px;
   padding-right: 5px;
 }
</style>
{% endhighlight %}

Note that we can specify the source code and the styling in the same file.  This has the potential
of greatly reducing the complexity of the code.
{:.notice}

The `MainPanel.ts` file is just an entry panel that acts as a container.  Therefore, it contains
little content, but is useful to illustrate the layout of a typical Vue2 component using
[av-ts](https://github.com/HerringtonDarkholme/av-ts) Typescript decorators.  Here's the code

{% highlight js %}
import * as Vue from "vue";
import {Component} from "av-ts";
@Component({
  name: 'main-panel',
})
export default class MainPanel extends Vue {
  name : string = "MainPanel";
}
{% endhighlight %}

#### Conclusion ####
We are now ready to get into the details of the actual graphics components.  I'll discuss how
I implemented that in the next part of this series.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

