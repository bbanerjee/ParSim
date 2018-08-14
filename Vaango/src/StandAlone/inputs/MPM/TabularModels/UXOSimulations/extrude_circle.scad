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

stem_dia = 5;
stem_rad = stem_dia / 2;
stem_height = 20;
stem_faces = 12;

wing_thick = 2*PI*stem_rad/stem_faces;
echo(wing_thick);
num_wings = 6;
wing_length = 5;
wing_height = 10;
wing_points = [[0, 0], [wing_length, 0], [0, wing_height]];
wing_angle = 360/num_wings;
stem_face_angle = 360/stem_faces;

module obj2D_lefs() {
  translate([0,0])
  circle(stem_dia, $fn=stem_faces);
}

module uxo_stem(diameter = 5, height = 10) {
  cylinder(d = diameter, h = height, $fn=stem_faces);
}

module uxo_wing() {
  in_by = 1;
  translate([stem_rad - in_by, 0, -stem_height])
  rotate([90, 0, 0])
  linear_extrude(height = wing_thick, center = true)
  polygon(wing_points);
}



union() {
  color([1,0,0]) 
  linear_extrude_fs(height=30);

  translate([0, 0, -stem_height])
  uxo_stem(diameter = stem_dia, height = stem_height);

  for (i = [0 : num_wings-1]) {
    rotate([0, 0, 360/num_wings*i + stem_face_angle/2])
    uxo_wing();
  }
}