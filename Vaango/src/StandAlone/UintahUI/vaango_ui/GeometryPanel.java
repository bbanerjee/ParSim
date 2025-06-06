package vaango_ui;
//**************************************************************************
// Program : GeometryPanel.java
// Purpose : Create geometry
// Author  : Biswajit Banerjee
// Date    : 05/12/2006
// Mods    :
//**************************************************************************

//************ IMPORTS **************

import java.awt.*;
import javax.swing.*;
import java.util.Vector;

//**************************************************************************
// Class   : GeometryPanel
// Purpose : Generate geometry (use particle distribution if available).
//**************************************************************************
public class GeometryPanel extends JPanel {

  /**
	 * 
	 */
	private static final long serialVersionUID = 6999161601032330600L;

// Data
  private double d_domainSize;

  private UintahInputPanel d_parentPanel = null;
  private Vector<GeomObject> d_geomObj = null;
  private Vector<GeomPiece> d_geomPiece = null;

  private InputGeometryPanel inputPanel = null;
  private DisplayGeometryFrame displayFrame = null;
  //private DisplayParticle3DFrame display3DFrame = null;

  // Constructor
  public GeometryPanel(ParticleList partList, 
                       Vector<GeomObject> geomObj,
                       UintahInputPanel parentPanel) {

    // Copy the arguments
    d_parentPanel = parentPanel;
    d_domainSize = 100.0;
    d_geomObj = geomObj;
    d_geomPiece = new Vector<GeomPiece>();

    // Create and add the relevant panels
    inputPanel = new InputGeometryPanel(partList, d_geomObj, d_geomPiece, this);
    displayFrame = new DisplayGeometryFrame(partList, d_geomPiece, this);
    displayFrame.setVisible(false);
    displayFrame.pack();
    //display3DFrame = new DisplayParticle3DFrame(partList, d_geomPiece);
    //display3DFrame.setVisible(true);
    //display3DFrame.pack();
 
    // Create a grid bag
    GridBagLayout gb = new GridBagLayout();
    GridBagConstraints gbc = new GridBagConstraints();
    setLayout(gb);

    // Set the constraints for the label
    UintahGui.setConstraints(gbc, GridBagConstraints.CENTER,
                             1.0, 1.0, 0, 0, 1, 1, 5);
    gb.setConstraints(inputPanel,gbc);
    add(inputPanel);
  }

  //-----------------------------------------------------------------------
  // Refresh
  //-----------------------------------------------------------------------
  public void refresh() {
    inputPanel.refresh();
  }

  //-------------------------------------------------------------------------
  // Update the components
  //-------------------------------------------------------------------------
  public void updatePanels() {
    validate();
    inputPanel.validate();
    displayFrame.validate();
    //display3DFrame.validate();
    d_parentPanel.updatePanels();
  }

  public UintahInputPanel getSuper() {
    return d_parentPanel;
  }

  public String getSimComponent() {
    return d_parentPanel.getSimComponent();
  }

  public void refreshDisplayGeometryFrame() {
    displayFrame.refresh();
    //display3DFrame.refresh();
  }

  public void setVisibleDisplayGeometryFrame(boolean visible) {
    displayFrame.setVisible(visible);
    //display3DFrame.setVisible(visible);
  }

  public void setDomainSize(double domainSize) {
    d_domainSize = domainSize;
  }

  public double getDomainSize() {
    return d_domainSize;
  }

  //-------------------------------------------------------------------------
  // Update geometry objects
  //-------------------------------------------------------------------------
  public void createPartListGeomObjects() {
    inputPanel.createPartListGeomObjects();
  }
}
