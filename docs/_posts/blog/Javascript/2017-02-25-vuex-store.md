---
layout: posts
title:  "Saving particle data in a Vuex store"
subheadline: "Biswajit Banerjee"
description: "Javascript scientific visualization - Part 2"
date:  2017-02-25 10:30:00
categories:
    - Javascript
    - Typescript
    - Vue
    - Vuex
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}

##### Introduction #####
In the [first part](http://www.parresianz.com/javascript/typescript/vue/vuex/xml/reading-vtk-particles/)
of this series on scientific visualization in Javascript, we saw how data in VTK unstructured grid
format produced by C++ simulation codes can be read in to a Javascript frontend.

Recall the `parseAndSaveVTKXML` function in Typescript syntax:

{% highlight js %}
public parseAndSaveVTKXML(xmlDoc : any) {
  let $xml = $(xmlDoc);
  ..............
  var pointData:any = {};
  ..............
  pointData["Position"] = arrayOfVec;
  // Save the data to a Vuex store
  Store.commit('SET_PARTICLE_DATA', pointData);
}
{% endhighlight %}

The `Store` in this code snippet, is a [Vuex store](https://vuex.vuejs.org/en/intro.html) that
we use to save and access data that will be needed by various components of the visualizer.  In
this article, I will discuss my implementation of the store.

##### Typescript typings for Vue and Vuex #####
I've found the typings provided by [av-ts](https://github.com/HerringtonDarkholme/av-ts) to
be more convenient than the ["official"](https://vuejs.org/v2/guide/typescript.html) version
for standard Vue-2, but for Vuex I use the official 
[version](https://github.com/vuejs/vuex/tree/Devi/types).  These type header files need
the "experimental decorators" option to be included in the `tsconfig.json` file used by
[Visual Studio Code](https://code.visualstudio.com/).

The compiler options in my `tsconfig.json` file are

{% highlight json %}
{
    "compilerOptions": {
        "outDir": "./dist",
        "sourceMap": false,
        "noImplicitAny": true,
        "target": "ES6",
        "lib": ["dom", "es6", "es2015.promise"],
        "module": "commonjs",
        "allowJs": true,
        "rootDirs": [
           "vaango_ui"
        ],
        "experimentalDecorators": true
    },
    ...........
}
{% endhighlight %}

##### The store #####
I create the Store in a file called `Store.ts`.  Because I am comparing the performance of
[Three.js](https://threejs.org/) and [vtk.js](https://kitware.github.io/vtk-js/), I include
two modules that correspond to states needed by these.  I also include a module for storing
and retrieving the particle data.  My `Store.ts` file is listed below:

{% highlight js %}
import * as Vue   from 'vue';
import * as Vuex  from 'vuex';
import {ThreeGraphicsModule} from './ThreeGraphicsModule';
import {VTKGraphicsModule}   from './VTKGraphicsModule';
import {ParticleModule}      from './ParticleModule';

Vue.use(Vuex);
export default new Vuex.Store({
  modules: {
             threeGraphics : new ThreeGraphicsModule(),
             vtkGraphics   : new VTKGraphicsModule(),
             particles     : new ParticleModule()  
           }
})
{% endhighlight %}

To make sure that the store is available to all components, we add the `store: Store` key-value pair
to the constructor in the file that acts as the entry point for the program, `main.ts`:

{% highlight js %}
import * as Vue from "vue";
import Store from "./Store";
class main {
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

##### The particle module #####
The particle module creates the particle state and defines actions, getters, and mutations on
this state.  My implementation is in the file `ParticleModule.ts` listed below:

{% highlight js %}
import { Module } from 'vuex';
import { ParticleState }    from './ParticleState';
import ParticleGetters      from './ParticleGetters';
import ParticleActions      from './ParticleActions';
import ParticleMutations    from './ParticleMutations';

export class ParticleModule implements Module<ParticleState, any> {
  public state : ParticleState;
  public getters   = ParticleGetters;
  public mutations = ParticleMutations;
  public actions   = ParticleActions;
  constructor() {
    this.state = new ParticleState();
  }
}
{% endhighlight %}

##### The particle state #####
The particle state contains the particle data that we want to store.  Our implementation of
`ParticleState.ts` is:

{% highlight js %}
export class ParticleState {
  public particleData: any;
  public particleReadComplete: boolean;
  constructor() {
    this.particleData = {};
    this.particleReadComplete = false;
  }
};
{% endhighlight %}

##### The particle getters #####
The getter functions are needed to access the particle store data from Vue components.  My
implementation of `ParticleGetters.ts` is shown below:

{% highlight js %}
import {Getter, GetterTree} from 'vuex';
import {ParticleState} from "./ParticleState";
export function particleData(state: ParticleState) {
  return state.particleData;
}
export function isParticleReadComplete(state: ParticleState) {
  return state.isParticleReadComplete;
}
export default <GetterTree<ParticleState, any>> {
  particleData,
  particleReadComplete
}
{% endhighlight %}

##### The particle mutations #####
To mutate the state of the store, we have to explicitly use mutation functions.  These are
implemented in `ParticleMutations.ts` as follows:

{% highlight js %}
import {Mutation, MutationTree} from 'vuex';
import {ParticleState} from './ParticleState';
export function SET_PARTICLE_DATA(state: ParticleState, particleData: any) {
  state.particleData = particleData;
  state.isParticleReadComplete = true;
  console.log("Set particleData");
}
export function DELETE_PARTICLE_DATA(state: ParticleState, message: string) {
  state.particleData = null;
  state.isParticleReadComplete = false;
  console.log("Deleted particle data" + message);
}
export default <MutationTree<ParticleState>> {
  SET_PARTICLE_DATA,
  DELETE_PARTICLE_DATA
}
{% endhighlight %}

The names of the functions are capitalized by convention.

The `SET_PARTICLE_DATA` mutation is how we save the particle data when we call
  `Store.commit('SET_PARTICLE_DATA', pointData);` in our function `parseAndSaveVTKXML`.
{: .notice}

##### The particle actions #####
Actions are similar to mutations but can perform asynchronous operations.  In our case, we don't
want to use actions at this stage but implement the functionality in `ParticleActions.ts`:

{% highlight js %}
import {Store, ActionTree, ActionContext} from 'vuex';
import {ParticleState} from "./ParticleState";
export function setParticleData(store: ActionContext<ParticleState, any>,
                                particleData: any)  {
  store.commit('SET_PARTICLE_DATA', particleData);
};
export function deleteParticleData(store: ActionContext<ParticleState, any>)  {
  store.commit('DELETE_PARTICLE_DATA', 'Deleting all particles');
};
export default <ActionTree<ParticleState, any>> {
  setParticleData,
  deleteParticleData
}
{% endhighlight %}

#### Remarks ####
Now that the particle data have been read in and stored, we can proceed with creating
a three-dimensional view that we can interact with.  We will explore a Three.js implementation 
in part 3 of this series.

Notice that we have used Typescript in our code.  This choice has the advantage that we can use
strong typing to catch errors during compile time.  However, we a limited to using libraries
that have been carefully assigned types by third-parties.  This is a serious limitation and
significantly reduces the attractiveness of Typescript as a development platform.
{: .notice}

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

