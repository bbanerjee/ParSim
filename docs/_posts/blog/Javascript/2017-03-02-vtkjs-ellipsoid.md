---
layout: posts
title:  "Plotting particles with vtk.js"
subheadline: "Biswajit Banerjee"
description: "Javascript scientific visualization - Part 7"
date:  2017-03-02 10:30:00
categories:
    - Javascript
    - Typescript
    - vtkjs 
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
In the previous articles in this series we talked about:

1. [Part 1: Reading VTK format particles with Javascript in a browser]({{ site.baseurl }}/javascript/typescript/vue/vuex/xml/reading-vtk-particles/)
2. [Part 2: Saving the read-in particle data in a Vuex store]({{ site.baseurl }}/javascript/typescript/vue/vuex/vuex-store/)
3. [Part 3: Initialization of a store and the user interface]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs/).
4. [Part 4: Setting up the Three.js renderer]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs-renderer/)
5. [Part 5: Setting up the Three.js scene and camera]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs-camera/)
5. [Part 6: Creating and plotting Three.js ellipsoids]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs-ellipsoid/)

Let us now explore the [vtk.js](https://kitware.github.io/vtk-js/) library to see if the process
of displaying particles can be simplified given that we are using VTK format input files.

##### Registering VTK components #####
As discussed in [Part 3]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs/), we first
need to register components that are specific to VTK.  In our case, we need

{% highlight js %}
var VtkRenderer = require("./graphics/VtkRenderer.vue").default;
var VtkEllipsoidParticles = require("./graphics/VtkEllipsoidParticles.vue").default;
Vue.component('vtk-renderer', VtkRenderer);
Vue.component('vtk-ellipsoid-particles', VtkEllipsoidParticles);
{% endhighlight %}

##### Creating Vuex stores for VTK data #####
Once again, we repeat the process discussed in
[Part 3]({{ site.baseurl }}/javascript/typescript/threejs/vue/vtk-threejs/).  We have discussed
the `Store` in Part3, and in this case we create the files

* `VTKGraphicsModule.ts` - similar to `ThreeGraphicsModule.ts`
* `VTKGraphicsState.ts` - see below
* `VTKGraphicsGetters.ts` - similar to the getters for the three.js Vuex store
* `VTKGraphicsMutations.ts` - see below
* `VTKGraphicsActions.ts` - similar to what we have seen before

The `VTKGraphicsState.ts` file now contains a list of actors, sources, and mappers.  There is also
a flag indicating that VTK actors have been created:

{% highlight js %}
export class VTKGraphicsState {
  actors     : any [];
  sources    : any [];
  mappers    : any [];
  areVTKActorsCreated  : boolean;
  constructor() {
    this.actors  = [];
    this.sources = [];
    this.mappers = [];
    this.areVTKActorsCreated = false;
  }
};
{% endhighlight %}

The `VTKGraphicsMutations.ts` file contains the following:

{% highlight js %}
import {Mutation, MutationTree} from 'vuex';
import {VTKGraphicsState} from './VTKGraphicsState';
export function ADD_VTK_ACTOR(state: VTKGraphicsState, actor: any) {
  state.actors.push(actor);
}
export function ADD_VTK_SOURCE(state: VTKGraphicsState, source: any) {
  state.sources.push(source);
}
export function ADD_VTK_MAPPER(state: VTKGraphicsState, mapper: any) {
  state.mappers.push(mapper);
}
export function VTK_ACTORS_CREATED(state: VTKGraphicsState, value: boolean) {
  state.areVTKActorsCreated = value;
}
export default <MutationTree<VTKGraphicsState>> {
  ADD_VTK_ACTOR,
  ADD_VTK_SOURCE,
  ADD_VTK_MAPPER,
  VTK_ACTORS_CREATED
}
{% endhighlight %}

##### The VTK graphics panel #####
The VTK graphics panel is simpler in this case.

{% highlight html %}
<template>
  <div id='vtk-graphics-container'>
    <div class="uk-card uk-card-default uk-card-large">
      <div class="uk-card-body">
        <vtk-renderer v-bind:size="size">
          <vtk-ellipsoid-particles> </vtk-ellipsoid-particles>
        </vtk-renderer>
      </div>
    </div>
  </div>
</template>
<script src="./VtkGraphicsPanel.ts"> </script>
{% endhighlight %}

The graphics panel template is simpler mainly because each actor is associated with a single
mapper and a single source in our implementation.  Better options will, one hopes,
become available after vtk.js reaches version 1.
{:.notice-warning}

##### The VTK renderer #####
We use a named slot for the template in the file `VtkRenderer.vue`:

{% highlight html %}
<template>
  <div>
    <slot name="vtk">
    <div ref="vtk-graphics-container"></div>
    <slot>
  </div>
</template>
<script src="./VtkRenderer.ts"> </script>
{% endhighlight %}

The code for the renderer is slight more complicated in this case and we add
an interactor to interact with the plot.  The `VtkRenderer.ts` code is listed below.

{% highlight js %}
import * as Vue from "vue";
import { Component, Lifecycle, Watch, Prop, p } from 'av-ts';
import Store from './Store';
import * as vtkRenderWindow           from 'vtk.js/Sources/Rendering/Core/RenderWindow';
import * as vtkRenderer               from 'vtk.js/Sources/Rendering/Core/Renderer';
import * as vtkOpenGLRenderWindow     from 'vtk.js/Sources/Rendering/OpenGL/RenderWindow';
import * as vtkRenderWindowInteractor from 'vtk.js/Sources/Rendering/Core/RenderWindowInteractor';
import * as vtkTexture                from 'vtk.js/Sources/Rendering/Core/Texture';
@Component({
  name: 'VtkRenderer'
})
export default class VtkRenderer extends Vue {
  @Prop
  size = p({
    type: Object, // { w, h }
    required: true
  });
  @Prop
  renderWindow = p({
    type: vtkRenderWindow
  });
  @Prop
  renderer = p({
    type: vtkRenderer
  });
  @Prop
  openGLRenderWindow = p({
    type: vtkOpenGLRenderWindow
  });
  @Prop
  interactor = p({
    type: vtkRenderWindowInteractor
  });
  private _renderWindow: any;
  private _renderer: any;
  private _openGLRenderWindow: any;
  private _interactor: any;
  @Lifecycle
  created() {
    this._renderWindow = this.renderWindow;
    this._renderer = this.renderer;
    this._openGLRenderWindow = this.openGLRenderWindow;
    this._interactor = this.interactor;
    // Create VTK render window and renderer
    this._renderWindow = vtkRenderWindow.newInstance();
    this._renderer = vtkRenderer.newInstance();
    this._renderWindow.addRenderer(this._renderer);
    this._renderer.setBackground(0.32, 0.34, 0.43);
    // Create OpenGL renderwindow
    this._openGLRenderWindow = vtkOpenGLRenderWindow.newInstance();
    this._renderWindow.addView(this._openGLRenderWindow);
    // Create interactor
    this._interactor = vtkRenderWindowInteractor.newInstance();
    // Add watch to check for data updates
    var self = this;
    Store.watch(
      function() { return Store.getters.areVTKActorsCreated; },
      function() {
        if (Store.getters.areVTKActorsCreated) {
          self.addActors();
          self._renderer.resetCamera();
          self._renderWindow.render();
        }
      }
    );
  }
  @Lifecycle
  mounted() {
    if ((this.$refs)["vtk-graphics-container"]) {
      let el  = (this.$refs)["vtk-graphics-container"];
      this._openGLRenderWindow.setContainer(el);
      // Set the size of the window
      let w = (<any>this.size).w;
      let h = (<any>this.size).h;
      this._openGLRenderWindow.setSize(w, h);
      // Add the actors from the store
      this.addActors();
      // Interactor
      this._interactor.setView(this._openGLRenderWindow);
      this._interactor.initialize();
      this._interactor.bindEvents(el);
      this._renderWindow.render();
    }
  }
  // Get actors from the store and add to the scene
  private addActors() {
    var self = this;
    let actors  = Store.getters.actors;
    actors.map(function(actor : any, index : number){
      self._renderer.addActor(actor);
    });
  }
}
{% endhighlight %}

The procedure used in this code is identical to that explained in the standard book on VTK.
You will be able to access that book [here](http://www.vtk.org/vtk-textbook/).
{:.notice--info}

The interactor allows only rotation and pan.  For zoom capabilities, you will need to you
a more capable interactor from the vtk.js library.
{:.notice--warning}

##### The VTK particles component #####
To draw the particles, we use the following code

{% highlight js %}
import * as Vue from "vue";
import {Data, Component, Lifecycle, Watch, Prop, p } from 'av-ts';
import * as vtkActor        from 'vtk.js/Sources/Rendering/Core/Actor';
import * as vtkSphereSource from 'vtk.js/Sources/Filters/Sources/SphereSource';
import * as vtkMapper       from 'vtk.js/Sources/Rendering/Core/Mapper';
import Store from "./Store";
@Component({
  name: 'VtkEllipsoidParticles'
})
export default class VtkEllipsoidParticles extends Vue {
  @Lifecycle
  public created() {
    var self = this;
    Store.watch(function() {
                  return Store.getters.isParticleReadComplete;
                },
                function() {
                  if (Store.getters.isParticleReadComplete)
                    self.createVTKParticles();
                });
  }
  private createVTKParticles() {
    // See below
  } 
{% endhighlight %}

This part of the code is identical to that used for the three.js particles.

##### Creating the particles #####
The particle axis data are in the form of angles between the ellipsoid axes and the world coordinate axes.
These are converted directly into the appropriate rotation matrix.  Sphere objects are then created at
the origin, rotated, scaled, and translated to their actual positions.  The sphere objects are then
transformed into `SphereBufferGeometry` objects to make their manipulation slightly more efficient.
Finally, a "material" shading model is added to make sure that the image displayed isn't flat and
a triangulated mesh is generated for each object.

{% highlight js %}
  private createVTKParticles() {
    // Get the particle data
    let particles = Store.getters.particleData;
    // Extract the radius and center
    let radii = particles["Radius"];
    let centers = particles["Position"];
    // Loop through particles
    radii.map(function(radius : any, index : number){
      // Create the _mapper
      const mapper = vtkMapper.newInstance();
      // Create the actor
      const actor = vtkActor.newInstance();
      actor.getProperty().setEdgeVisibility(true);
      actor.getProperty().setEdgeColor(1.0, 0.5, 0.5);
      // Create the source
      const sphere = vtkSphereSource.newInstance();
      sphere.setPhiResolution(10);
      sphere.setThetaResolution(10);
      // Get the radius ratios
      let ratio = [1.0, radius[1]/radius[0], radius[2]/radius[0]];
      sphere.setRadius(radius[0]);
      let center = centers[index];
      sphere.setCenter(center[0], center[1], center[2]);
      // Set up the connections
      mapper.setInputConnection(sphere.getOutputPort());
      actor.setMapper(mapper);
      // Save the data
      Store.commit('ADD_VTK_ACTOR',  actor);
      Store.commit('ADD_VTK_SOURCE', sphere);
      Store.commit('ADD_VTK_MAPPER', mapper);
    });
    Store.commit('VTK_ACTORS_CREATED', true);
  }
{% endhighlight %}

Note that we did not rotate or scale the particles.  To do that you will have to write
your own code or wait until a more full-featured version of VTK is implemented by vtk.js.
{:.notice--warning}

#### Remarks ####
I ran into several issues while trying to use `vtk.js`.  The main ones were:

* The `vtk.js` library needs to be transpiled using `babel` to avoid complaints about `import` statements
in the `webpack` bundle.
* All events are captured by the interactors used by `vtk.js` and it is difficult to use Vue's
event handling capabilities with this library.
* Geometry transformations have not been implemented in `vtk.js` yet.
* Unstructured grid VTK XML files cannot yet be read in by `vtk.js` readers.
* The rendering process was much slower than that of `three.js`.
{:.notice}

The main advantages of `vtk.js` were:

* The familiar architecture and API which is similar to what we have done in scientific visualization
for more than 20 years.
* The easy application of interactors (which was not the case with Typescript, Vue, and three.js, and
an interactor had to be written).  We will discuss the three.js interactor in a future article.
{:.notice}

A plot of the spheres produced by our code can be seen below.

![Plot produced by Three.js]({{site.baseurl}}/assets/blogimg/VTKGraphicsPanel.jpg){:class="img-responsive center-image" height="450px" border="5px double red"}

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

