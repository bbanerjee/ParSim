package vaango_ui;

//**************************************************************************
//Class   : SphereGeomPiece
//Purpose : Creates a cylinder
//**************************************************************************
import java.lang.System;
import java.io.PrintWriter;


public class SphereGeomPiece extends GeomPiece {

  // Common data
  private Point d_center = null;
  private double d_radius = 0.0;

  public SphereGeomPiece() {
    d_name = new String("Sphere");
    d_center = new Point();
    d_radius = 0.0;
  }

  public SphereGeomPiece(String name, Point center, double radius) { 
    d_name = new String(name);
    d_radius = radius;
    d_center = new Point(center);
  }

  public String getName() {
    return d_name;
  }

  // Common Methods
  public void writeUintah(PrintWriter pw, String tab){

    String tab1 = new String(tab+"  ");
    pw.println(tab+"<sphere label=\""+d_name+"\">");
    pw.println(tab1+"<origin> ["+d_center.getX()+", "+d_center.getY()+", "+
        d_center.getZ()+"] </origin>");
    pw.println(tab1+"<radius> "+d_radius+" </radius>");
    pw.println(tab+"</sphere>");
  }

  public void print(){
    String tab1 = new String("  ");
    System.out.println("<sphere label=\""+d_name+"\">");
    System.out.println(tab1+"<origin> ["+d_center.getX()+", "+
        d_center.getY()+", "+ d_center.getZ()+"] </origin>");
    System.out.println(tab1+"<radius> "+d_radius+" </radius>");
    System.out.println("</sphere>");
  }


}
