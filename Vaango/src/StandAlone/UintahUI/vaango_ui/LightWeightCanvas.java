package vaango_ui;
//**************************************************************************
// Program : LightWeightCanvas.java
// Purpose : A canvas that is light weight
// Author  : Biswajit Banerjee
// Date    : 12/7/1998
// Mods    :
//**************************************************************************

//************ IMPORTS **************
import java.awt.*;
import javax.swing.*;

//**************************************************************************
// Class   : LightWeightCanvas
// Purpose : Creates a light weight canvas
//**************************************************************************
public class LightWeightCanvas extends JComponent {

  // Data

  // Data that may be needed later

  /**
	 * 
	 */
	private static final long serialVersionUID = 4745691118012860762L;

// Constructor
  public LightWeightCanvas(int width, int height) {

    // Set the size of the component
    setSize(width, height);

    // Set the preferrred size of the component
    setPreferredSize(new Dimension(width, height));
  }

  // Paint the component
  @Override
public void paintComponent(Graphics g) {

    Dimension d = getSize();
    g.drawRect(0,0,d.width,d.height);
  }
}
