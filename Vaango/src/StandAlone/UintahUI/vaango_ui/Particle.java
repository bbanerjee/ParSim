package vaango_ui;

import java.io.PrintWriter;
//**************************************************************************
// Program : Particle.java
// Purpose : Define the Particle objects
// Author  : Biswajit Banerjee
// Date    : 03/24/1999
// Mods    :
//**************************************************************************

//**************************************************************************
// Class   : Particle
// Purpose : Creates a Particle object
//**************************************************************************
public class Particle extends Object {

  // Data
  private int d_type = 1;
  private double d_radius = 0.0;
  private double d_length = 1.0;
  private double d_thickness = 0.0;
  private double d_rotation = 0.0;
  private Point d_center = null;
  private int d_matCode = 0;

  // Public static data
  static final int CIRCLE = 1;
  static final int SPHERE = 2;

  // Constructor
  public Particle() {
    d_type = CIRCLE;
    d_radius = 0.0;
    d_length = 1.0;
    d_thickness = 0.0;
    d_rotation = 0.0;
    d_center = new Point(0.0,0.0,0.0);
    d_matCode = 0;
  }

  public Particle(int type) {
    d_type = type;
    d_radius = 0.0;
    d_length = 1.0;
    d_thickness = 0.0;
    d_rotation = 0.0;
    d_center = new Point(0.0,0.0,0.0);
    d_matCode = 0;
  }

  public Particle(double radius, double length, Point center, int matCode) {
    d_type = CIRCLE;
    d_radius = radius;
    d_length = length;
    d_thickness = 0.0;
    d_rotation = 0.0;
    d_center = new Point(center);
    d_matCode = matCode;
  }

  public Particle(double radius, double length, double thickness,
                  Point center, int matCode) {
    d_type = CIRCLE;
    d_radius = radius;
    d_length = length;
    d_thickness = thickness;
    d_rotation = 0.0;
    d_center = new Point(center);
    d_matCode = matCode;
  }

  public Particle(int type, double radius, double length, double thickness,
                  Point center, int matCode) {
    d_type = type;
    d_radius = radius;
    d_length = length;
    d_thickness = thickness;
    d_rotation = 0.0;
    d_center = new Point(center);
    d_matCode = matCode;
  }

  public Particle(int type, double radius, double rotation, Point center,
		  int matCode) {
    d_type = type;
    d_radius = radius;
    d_length = 0.0;
    d_thickness = 0.0;
    d_rotation = rotation;
    d_center = new Point(center);
    d_matCode = matCode;
  }

  public Particle(int type, double radius, double rotation, Point center, 
                  int matCode, double thickness) {
    d_type = type;
    d_radius = radius;
    d_length = 0.0;
    d_thickness = thickness;
    d_rotation = rotation;
    d_center = new Point(center);
    d_matCode = matCode;
  }

  // Get Particle data
  public int getType() {return d_type;}
  public int getMatCode() {return d_matCode;}
  public double getRadius() {return d_radius;}
  public double getLength() {return d_length;}
  public double getThickness() {return d_thickness;}
  public double getRotation() {return d_rotation;}
  public Point getCenter() {return d_center;}
  public double getVolume() {
    if (d_type == CIRCLE) {
      return Math.PI*d_radius*d_radius*d_length; 
    } else {
      return Math.PI*d_radius*d_radius*d_radius*(4.0/3.0); 
    }
  }

  // Set Particles data
  public void setType(int type) {d_type = type;}
  public void setMatCode(int matCode) { d_matCode = matCode;}
  public void setRadius(double radius) { d_radius = radius;}
  public void setLength(double length) { d_length = length;}
  public void setThickness(double thickness) { d_thickness = thickness;}
  public void setRotation(double rotation) {d_rotation = rotation;}
  public void setCenter(Point center) { d_center = center;}

  // Print the particle data
  public void print() 
  {
    System.out.println("Material Code = "+d_matCode+" Type = "+d_type+
                       " Rad = "+d_radius+" Length = "+d_length+
                       " Thick = "+d_thickness+
                       " Rotation = "+d_rotation+ " Center = ["+
                       d_center.getX()+", "+d_center.getY()+", "+
                       d_center.getZ()+"]"); 
  }
  
  // Print the particle data to a file
  public void print(PrintWriter pw, String tab)
  {
    String tab1 = new String(tab+"  ");
    switch (d_type) 
    {
    case CIRCLE:
      pw.println(tab+"<cylinder label=\""+d_matCode+"\">");
      pw.println(tab1+"<bottom> ["+d_center.getX()+", "+d_center.getY()+", "+
                 d_center.getZ()+"] </bottom>");
      double zcoord = d_center.getZ() + d_length;
      pw.println(tab1+"<top> ["+d_center.getX()+", "+d_center.getY()+", "+
                 zcoord+"] </top>");
      pw.println(tab1+"<radius> "+d_radius+" </radius>");
      pw.println(tab1+"<thickness> "+d_thickness+" </thickness>");
      pw.println(tab+"</cylinder>");
      break;
      
    case SPHERE:
      pw.println(tab+"<sphere label=\""+d_matCode+"\">");
      pw.println(tab1+"<center> ["+d_center.getX()+", "+d_center.getY()+", "+
                 d_center.getZ()+"] </center>");
      pw.println(tab1+"<radius> "+d_radius+" </radius>");
      pw.println(tab1+"<thickness> "+d_thickness+" </thickness>");
      pw.println(tab+"</sphere>");
      break;
      
    default:
      System.out.println("Not output method for particle type "+d_type+" implemented yet.");
    }
  }
}
