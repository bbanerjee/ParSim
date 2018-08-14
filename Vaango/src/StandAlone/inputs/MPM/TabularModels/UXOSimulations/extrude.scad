module linear_extrude_fs(height = 1, isteps = 20, twist = 0) {

  union() {

    for (i = [0 : 1: isteps-1]) {

       translate([0, 0, i*height/isteps]) 
       linear_extrude(
         height = height/isteps,
         twist = 0,
         scale = f_lefs((i+1)/isteps) / f_lefs(i/isteps)
       )
       scale(f_lefs(i/isteps))
       obj2D_lefs();
    }
  }
}

function f_lefs(x) = 
  let(span = 150, start = 20, normpos = 45)
  sin(x*span + start) / sin(normpos);

module obj2D_lefs() {
  translate([-4,-3])
  square([9,12]);
}

color([1,0,0]) linear_extrude_fs(height=20);