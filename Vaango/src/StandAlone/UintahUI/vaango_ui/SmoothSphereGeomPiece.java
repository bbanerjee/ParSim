package vaango_ui;

//**************************************************************************
//Class   : SmoothSphereGeomPiece
//Purpose : Creates a smooth sphere
//**************************************************************************
import java.lang.System;
import java.io.PrintWriter;

public class SmoothSphereGeomPiece extends GeomPiece {

  // Common data
  private Point d_center = null;
  private double d_outerRadius = 0.0;
  private double d_innerRadius = 0.0;
  private int d_numRadial = 0;

  public SmoothSphereGeomPiece() {
    d_name = new String("SmoothSphere");
    d_center = new Point();
    d_outerRadius = 0.0;
    d_innerRadius = 0.0;
    d_numRadial = 0;
  }

  public SmoothSphereGeomPiece(String name, Point center, double outerRadius, 
      double innerRadius) {
    d_name = new String(name);
    d_outerRadius = outerRadius;
    d_innerRadius = innerRadius;
    d_numRadial = 0;
    d_center = new Point(center.getX(), center.getY(), center.getZ());
  }

  public SmoothSphereGeomPiece(String name, Point center, double outerRadius, 
      double innerRadius,  int numRadial) {
    d_name = new String(name);
    d_outerRadius = outerRadius;
    d_innerRadius = innerRadius;
    d_numRadial = numRadial;
    d_center = new Point(center.getX(), center.getY(), center.getZ());
  }

  public String getName() {
    return d_name;
  }

  // Common Methods
  public void writeUintah(PrintWriter pw, String tab){

    String tab1 = new String(tab+"  ");
    pw.println(tab+"<smooth_sphere label=\""+d_name+"\">");
    pw.println(tab1+"<center> ["+d_center.getX()+", "+d_center.getY()+", "+
        d_center.getZ()+"] </center>");
    pw.println(tab1+"<outer_radius> "+d_outerRadius+" </outer_radius>");
    pw.println(tab1+"<inner_radius> "+d_innerRadius+" </inner_radius>");
    pw.println(tab1+"<num_radial_pts> "+d_numRadial+" </num_radial_pts>");
    pw.println(tab1+"<algorithm> equal_area </algorithm>");
    pw.println(tab+"</smooth_sphere>");
  }

  public void print(){
    String tab1 = new String("  ");
    System.out.println("<smooth_spherel label=\""+d_name+"\">");
    System.out.println(tab1+"<center> ["+d_center.getX()+", "+
        d_center.getY()+", "+ d_center.getZ()+"] </center>");
    System.out.println(tab1+"<outer_radius> "+d_outerRadius+" </outer_radius>");
    System.out.println(tab1+"<inner_radius> "+d_innerRadius+" </inner_radius>");
    System.out.println(tab1+"<num_radial_pts> "+d_numRadial+" </num_radial_pts>");
    System.out.println(tab1+"<algorithm> equal_area </algorithm>");
    System.out.println("</smooth_sphere>");
  }


}
