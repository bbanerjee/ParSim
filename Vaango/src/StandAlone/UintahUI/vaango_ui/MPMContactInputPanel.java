package vaango_ui;
//**************************************************************************
// Class   : MPMContactInputPanel.java
// Purpose : Create a panel that contains widgets to take inputs for
//           contact between MPM materials
// Author  : Biswajit Banerjee
// Date    : 05/17/2006
// Mods    :
//**************************************************************************


import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.Vector;

public class MPMContactInputPanel extends JPanel 
                                  implements ItemListener,
                                             ActionListener {

  /**
	 * 
	 */
	private static final long serialVersionUID = 8045934468645835561L;
// Data and components
  private String d_contactType = null;
  private double[] d_contactDir = null;
  private Vector<String> d_mpmMat = null;
  private int[] d_selMat = null;

  private JComboBox<String> contactTypeComB = null;
  private JList<String> contactMatList = null;
  private DefaultListModel<String> contactMatListModel = null;
  private JScrollPane contactMatSP = null;
  private JLabel frictionLabel = null;
  private DecimalField frictionEntry = null;
  private JLabel contactDirLabel = null;
  private JComboBox<String> contactDirComB = null;
  private JButton updateButton = null;

  public MPMContactInputPanel(Vector<String> mpmMat) {

    // Initialize
    d_contactType = new String("null");
    d_contactDir = new double[3];
    d_contactDir[0] = 0.0;
    d_contactDir[1] = 0.0;
    d_contactDir[2] = 0.0;
    d_mpmMat = mpmMat;

    // Create a gridbaglayout and constraints
    GridBagLayout gb = new GridBagLayout();
    GridBagConstraints gbc = new GridBagConstraints();
    setLayout(gb);

    JLabel contactTypeLabel = new JLabel("Contact Type");
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 0, 0);
    gb.setConstraints(contactTypeLabel, gbc);
    add(contactTypeLabel); 

    contactTypeComB = new JComboBox<String>();
    contactTypeComB.addItem("No Contact");
    contactTypeComB.addItem("Rigid Contact");
    contactTypeComB.addItem("Specified Velocity Contact");
    contactTypeComB.addItem("Single Velocity Contact");
    contactTypeComB.addItem("Approach Contact");
    contactTypeComB.addItem("Frictional Contact");
    contactTypeComB.addItemListener(this);
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 1, 0);
    gb.setConstraints(contactTypeComB, gbc);
    add(contactTypeComB); 
    
    JLabel contactMatLabel = new JLabel("Contact Materials");
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 0, 1);
    gb.setConstraints(contactMatLabel, gbc);
    add(contactMatLabel); 

    contactMatListModel = new DefaultListModel<String>();
    for (int ii = 0; ii < d_mpmMat.size(); ++ii) {
      contactMatListModel.addElement((String) d_mpmMat.elementAt(ii));
    }
    contactMatList = new JList<String>(contactMatListModel);
    contactMatList.setVisibleRowCount(5);
    contactMatSP = new JScrollPane(contactMatList);
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 1, 1);
    gb.setConstraints(contactMatSP, gbc);
    add(contactMatSP); 
    
    frictionLabel = new JLabel("Friction Coefficient");
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 0, 2);
    gb.setConstraints(frictionLabel, gbc);
    add(frictionLabel); 
    frictionLabel.setEnabled(false);

    frictionEntry = new DecimalField(0.0, 7);
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 1, 2);
    gb.setConstraints(frictionEntry, gbc);
    add(frictionEntry); 
    frictionEntry.setEnabled(false);

    contactDirLabel = new JLabel("Initial Contact Direction");
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 0, 3);
    gb.setConstraints(contactDirLabel, gbc);
    add(contactDirLabel); 
    contactDirLabel.setEnabled(false);

    contactDirComB = new JComboBox<String>();
    contactDirComB.addItem("X-Direction");
    contactDirComB.addItem("Y-Direction");
    contactDirComB.addItem("Z-Direction");
    contactDirComB.addItemListener(this);
    UintahGui.setConstraints(gbc, GridBagConstraints.NONE, 1, 3);
    gb.setConstraints(contactDirComB, gbc);
    add(contactDirComB); 
    contactDirComB.setEnabled(false);

    updateButton = new JButton("Update");
    UintahGui.setConstraints(gbc, 0, 4);
    gb.setConstraints(updateButton, gbc);
    add(updateButton); 
    updateButton.addActionListener(this);

  }

  @Override
public void actionPerformed(ActionEvent e) {
    d_selMat = contactMatList.getSelectedIndices();
    contactMatList.setSelectedIndices(d_selMat);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Listens for item picked in combo box and takes action as required.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  @Override
public void itemStateChanged(ItemEvent e) {
        
    // Get the object that has been selected
    Object source = e.getItemSelectable();

    // Get the item that has been selected
    String item = String.valueOf(e.getItem());

    if (source.equals(contactTypeComB)) {
      frictionLabel.setEnabled(false);
      frictionEntry.setEnabled(false);
      contactDirLabel.setEnabled(false);
      contactDirComB.setEnabled(false);
      if (item.equals("No Contact")) {
        d_contactType = "null";
      } else if (item.equals("Rigid Contact")) {
        d_contactType = "rigid";
        contactDirLabel.setEnabled(true);
        contactDirComB.setEnabled(true);
      } else if (item.equals("Specified Velocity Contact")) {
        d_contactType = "specified";
      } else if (item.equals("Single Velocity Contact")) {
        d_contactType = "single_velocity";
      } else if (item.equals("Approach Contact")) {
        d_contactType = "approach";
        frictionLabel.setEnabled(true);
        frictionEntry.setEnabled(true);
      } else if (item.equals("Frictional Contact")) {
        d_contactType = "friction";
        frictionLabel.setEnabled(true);
        frictionEntry.setEnabled(true);
      }
    } else if (source.equals(contactDirComB)) {
      if (item.equals("X-Direction")) {
        d_contactDir[0] = 1.0;
        d_contactDir[1] = 0.0;
        d_contactDir[2] = 0.0;
      } else if (item.equals("Y-Direction")) {
        d_contactDir[0] = 0.0;
        d_contactDir[1] = 1.0;
        d_contactDir[2] = 0.0;
      } else if (item.equals("Z-Direction")) {
        d_contactDir[0] = 0.0;
        d_contactDir[1] = 0.0;
        d_contactDir[2] = 1.0;
      }

    }
  }

  //--------------------------------------------------------------------
  // Refresh the panel
  //--------------------------------------------------------------------
  public void refresh() {
    contactMatListModel.removeAllElements();
    for (int ii = 0; ii < d_mpmMat.size(); ++ii) {
      contactMatListModel.addElement((String) d_mpmMat.elementAt(ii));
    }
    if (d_selMat != null) {
      contactMatList.setSelectedIndices(d_selMat);
    }
    validate();
  }

  //--------------------------------------------------------------------
  /** Write the contents out in Uintah format */
  //--------------------------------------------------------------------
  public void writeUintah(PrintWriter pw, String tab) {
      
    if (pw == null) return;

    if (d_contactType.equals("null")) return;
    
    String tab1 = new String(tab+"  ");
    try {
    pw.println(tab+"<contact>");
    pw.println(tab1+"<type> "+d_contactType+" </type>");
    if (d_contactType.equals("rigid")) {
      pw.println(tab1+"<direction> ["+d_contactDir[0]+", "+
                      d_contactDir[1]+", "+d_contactDir[2]+"] </direction>");
    }
    } catch (Exception e) {
      System.out.println("failed contact part 1");
    }
    try {
    if (d_selMat.length > 0) {
      pw.print(tab1+"<materials> [");
      for (int ii = 0; ii < d_selMat.length; ++ii) {
        if (ii == d_selMat.length - 1) {
          pw.print(d_selMat[ii]);
        } else {
          pw.print(d_selMat[ii]+",");
        }
      }
    }
    } catch (Exception e) {
      System.out.println("failed contact part 2");
    }
    try {
    pw.println("] </materials>");
    if (d_contactType.equals("friction") ||
        d_contactType.equals("approach")) {
      pw.println(tab1+"<mu> "+frictionEntry.getValue()+" </mu>");
    }
    pw.println(tab1+"<stop_time> 9999999.9 </stop_time>");
    pw.println(tab+"</contact>");
    } catch (Exception e) {
      System.out.println("failed contact part 3");
    }

  }
}
