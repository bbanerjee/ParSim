package vaango_ui;
//**************************************************************************
// Program : ViscoSCRAMMaterialPanel.java
// Purpose : Create a panel that contains widgets to take inputs for
//           ViscoSCRAM materials
// Author  : Biswajit Banerjee
// Date    : 05/04/2006
// Mods    :
//**************************************************************************

//import java.awt.event.*;
import java.io.*;
import javax.swing.*;

public class ViscoSCRAMMaterialPanel extends JPanel {

  /**
	 * 
	 */
	private static final long serialVersionUID = -7468356791990080074L;

// Data and components
  public ViscoSCRAMMaterialPanel() {

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Class   : ComboBoxListener
  // Purpose : Listens for item picked in combo box and takes action as
  //           required.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*private class ComboBoxListener implements ItemListener {
    @Override
	public void itemStateChanged(ItemEvent e) {
        
    }
  }*/

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Class   : CheckBoxListener
  // Purpose : Listens for item seleceted in check box and takes action as
  //           required.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*private class CheckBoxListener implements ItemListener {
    @Override
	public void itemStateChanged(ItemEvent e) {
        
    }
  }*/

  //--------------------------------------------------------------------
  /** Write the contents out in Uintah format */
  //--------------------------------------------------------------------
  public void writeUintah(PrintWriter pw) {
      
    if (pw == null) return;

    // Write the data
  }
}
