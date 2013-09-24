package vaango_ui;
//**************************************************************************
// Class   : DisplayParticle3DFrame
// Purpose : Draws a Canvas3D for displaying the 3D view of the particle list.
//**************************************************************************

//************ IMPORTS **************
import java.awt.*;
import java.util.Vector;
import java.util.Enumeration;
import javax.media.j3d.*;
import javax.vecmath.*;
import javax.swing.*;
//import com.sun.j3d.utils.geometry.Box;
import com.sun.j3d.utils.geometry.ColorCube;
import com.sun.j3d.utils.geometry.Sphere;
import com.sun.j3d.utils.geometry.Cylinder;
import com.sun.j3d.utils.geometry.Cone;
import com.sun.j3d.utils.universe.SimpleUniverse;

public class DisplayParticle3DFrame extends JFrame {

  // Data
  private static final long serialVersionUID = -2117067419917398344L;
	
  private DisplayParticle3DPanel d_panel = null;

  private static final int WIDTH = 500;	
  private static final int HEIGHT = 500;	

  //-------------------------------------------------------------------------
  // Constructor
  //-------------------------------------------------------------------------
  public DisplayParticle3DFrame(ParticleList partList)
  {
    // Set position and title
    setLocation(200, 100);
    setTitle("3D Geometry View");
      
    // Create the 3D view panel
    d_panel = new DisplayParticle3DPanel(partList);

    // Add the display panel to the frame
    getContentPane().add(d_panel);
    
    // Resize the frame and make visible
    setSize(WIDTH, HEIGHT);
    setVisible(true);
  }

  public void refresh() {
    setVisible(true);
	d_panel.refresh();
  }

  private class DisplayParticle3DPanel extends JPanel {

    /**
	 * 
	 */
	private static final long serialVersionUID = 3545544188004609798L;

    private static final int CIRCLE = 1;
    private static final int SPHERE = 2;

    private double d_domainSize;
    private ParticleList d_partList = null;
    private Vector<GeomPiece> d_geomPiece;

    private SimpleUniverse d_universe = null;
    private Locale d_locale = null;
    
	private BranchGroup d_baseGroup = null;
	private BranchGroup d_sceneGroup = null;
	private BranchGroup d_viewGroup = null;
	private View d_view = null;
	private Canvas3D d_canvas = null;
	
    public DisplayParticle3DPanel(ParticleList partList) {
	  
      d_partList = partList;
      if (d_partList != null) {
        d_domainSize = d_partList.getRVESize();
      } else {
    	d_domainSize = 100.0;
      }
      
      setLayout(new BorderLayout());
      initialize();
      //testInitialize();
      //testInitialize2();
    }

    private void initialize()
    {
      // Create canvas
	  d_canvas = createCanvas3D(false);
      add("Center", d_canvas);

      // Create universe
      d_universe = new SimpleUniverse(d_canvas);

      // Create base and scene branch group
      d_baseGroup = new BranchGroup();
      d_sceneGroup = new BranchGroup();

      // Create background
      Background background = createBackground();
      if (background != null) d_baseGroup.addChild(background);

      // Create lighting
	  Color3f lightColor = new Color3f(0.7f, 0.7f, 0.7f);
	  Vector3f lightDir = new Vector3f(-1.0f, -1.0f, -1.0f);
	  Color3f ambientColor = new Color3f(0.2f, 0.2f, 0.2f);
	  BoundingSphere bound_sph = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
	  AmbientLight ambientLight = new AmbientLight(ambientColor);
	  ambientLight.setInfluencingBounds(bound_sph);
	  DirectionalLight light1 = new DirectionalLight(lightColor, lightDir);
	  light1.setInfluencingBounds(bound_sph);

	  d_baseGroup.addChild(ambientLight);
	  d_baseGroup.addChild(light1);

	  // Create lighting
	  /*
	  Color3f light1Color = new Color3f(.1f, 1.4f, .1f); // green light
      BoundingSphere bounds =
         new BoundingSphere(new Point3d(0.0,0.0,0.0), 100.0);
      Vector3f light1Direction = new Vector3f(4.0f, -7.0f, -12.0f);
      DirectionalLight light1
         = new DirectionalLight(light1Color, light1Direction);
      light1.setInfluencingBounds(bounds);
      d_baseGroup.addChild(light1);
      */

      // Create the bounding box
	  Appearance appear = new Appearance();
	  //appear.setTransparencyAttributes(
	  //		  new TransparencyAttributes(TransparencyAttributes.FASTEST, 0.8f));
	  appear.setPolygonAttributes(new PolygonAttributes( PolygonAttributes.POLYGON_LINE, 
			  PolygonAttributes.CULL_NONE, 0) );
	  Color3f color = new Color3f(1.0f, 0.7f, 0.8f);
	  Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
	  appear.setMaterial(new Material(color, black, color, black, 80.0f));
      //d_baseGroup.addChild(new Box(1.0f, 1.0f, 1.0f, appear));
	  ColorCube cube = new ColorCube();
	  cube.setAppearance(appear);
      d_baseGroup.addChild(cube);

      // Create the particles
      //System.out.println("Calling create particles for 3D display.");
      //createParticles();
	  //d_baseGroup.addChild(createCylinder(20.0, 50.0, 50.0));

	  // Change the view direction
      Transform3D locator = new Transform3D();
      locator.lookAt(new Point3d(1.5, 1.5, 4.0), new Point3d(-0.25, -0.25, 0.0), 
                     new Vector3d(0.0, 1.0, 0.0));
      locator.invert();
      d_universe.getViewingPlatform().getViewPlatformTransform().setTransform(locator);
      //universe.getViewer().getView().setProjectionPolicy(View.PARALLEL_PROJECTION);
      
      // Add the scene to the universe
      d_universe.addBranchGraph(d_baseGroup); 
      
    }

	public void refresh() {
	
	  System.out.println("Refreshing DisplayGeometry3DFrame");
	  removeShape();

	  createParticles();
	  //d_sceneGroup.addChild(createSphere(10.0, 10.0, 20.0, 20.0));
	  //d_sceneGroup.addChild(createSphere(20.0, 20.0, 50.0, 70.0));
      d_universe.addBranchGraph(d_sceneGroup); 
	}

	// Remove a branchgroup from the scene and redraw
	private void removeShape() {
	  try {
		Enumeration<SceneGraphObject> en = d_sceneGroup.getAllChildren();  
		int index = 0;
		
		while (en.hasMoreElements() != false) {
		  d_sceneGroup.removeChild(index);
		  index++;
		}
	  } catch (Exception e) {
		System.out.println("Scenegraph not synchronized."); 
	  }
	}
    private void testInitialize()
    {
      GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
      Canvas3D canvas = new Canvas3D(config);
      add("Center", canvas);

      BranchGroup group = new BranchGroup();
      for (float x = -1.0f; x <= 1.0f; x = x + 0.1f)
      {
        Sphere sphere = new Sphere(0.05f);
        TransformGroup tg = new TransformGroup();
        Transform3D transform = new Transform3D();
        Vector3f vector = new Vector3f( x, .0f, .0f);
        transform.setTranslation(vector);
        tg.setTransform(transform);
        tg.addChild(sphere);
        group.addChild(tg);
      }

      for (float y = -1.0f; y <= 1.0f; y = y + 0.1f)
      {
        TransformGroup tg = new TransformGroup();
        Transform3D transform = new Transform3D();
        Cone cone = new Cone(0.05f, 0.1f);
        Vector3f vector = new Vector3f(.0f, y, .0f);
        transform.setTranslation(vector);
        tg.setTransform(transform);
        tg.addChild(cone);
        group.addChild(tg);
      }

      for (float z = -1.0f; z <= 1.0f; z = z+ 0.1f)
      {
        TransformGroup tg = new TransformGroup();
        Transform3D transform = new Transform3D();
        Cylinder cylinder = new Cylinder(0.05f, 0.1f);
        Vector3f vector = new Vector3f(.0f, .0f, z);
        transform.setTranslation(vector);
        tg.setTransform(transform);
        tg.addChild(cylinder);
        group.addChild(tg);
      }

      Color3f light1Color = new Color3f(.1f, 1.4f, .1f); // green light
      BoundingSphere bounds =
         new BoundingSphere(new Point3d(0.0,0.0,0.0), 100.0);
      Vector3f light1Direction = new Vector3f(4.0f, -7.0f, -12.0f);
      DirectionalLight light1
         = new DirectionalLight(light1Color, light1Direction);
      light1.setInfluencingBounds(bounds);
      group.addChild(light1);
      
      SimpleUniverse universe = new SimpleUniverse(canvas);
      universe.getViewingPlatform().setNominalViewingTransform();
      universe.addBranchGraph(group); 
    }

    //-------------------------------------------------------------------------
    // initialize
    //-------------------------------------------------------------------------
    private void testInitialize2() {

      // Create universe
      d_universe = new SimpleUniverse();

      // Create local
      d_locale = new Locale(d_universe);

	  // Create scene branch group
      d_sceneGroup = createSceneBranchGroup();

      // Create background and add to scene branch group 
      Background bg = createBackground();
      if (bg != null) d_sceneGroup.addChild(bg);
      
      // Create the view platform
      ViewPlatform vp = createViewPlatform();
      
      // Create view branch group
      d_viewGroup = createViewBranchGroup(vp);
      
      // Add the scene branch group to the universe
      d_locale.addBranchGraph(d_sceneGroup);
      
      // Add the view branch group to the universe
	  d_locale.addBranchGraph(d_viewGroup);	

	  // Create view
      d_view = createView(vp);

	  // Create canvas and add to view
	  d_canvas = createCanvas3D(false);
	  d_view.addCanvas3D(d_canvas);
	  add("Center", d_canvas);
	  
    }

    // Create the scene branch group
	private BranchGroup createSceneBranchGroup() {
	  
	  BranchGroup group = new BranchGroup();
	  group.setCapability(Group.ALLOW_CHILDREN_EXTEND);
	  group.setCapability(Group.ALLOW_CHILDREN_READ);
	  group.setCapability(Group.ALLOW_CHILDREN_WRITE);

	  // Create particles and add to branches
	  group.addChild(createSphere(0.5, 10.0, 20.0, 20.0));
	  //createParticles();

	  // Create colors and lights
	  Color3f lightColor = new Color3f(0.7f, 0.7f, 0.7f);
	  Vector3f lightDir = new Vector3f(-1.0f, -1.0f, -1.0f);
	  Color3f ambientColor = new Color3f(0.2f, 0.2f, 0.2f);
	  BoundingSphere bound_sph = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
	  AmbientLight ambientLight = new AmbientLight(ambientColor);
	  ambientLight.setInfluencingBounds(bound_sph);
	  DirectionalLight light1 = new DirectionalLight(lightColor, lightDir);
	  light1.setInfluencingBounds(bound_sph);

	  group.addChild(ambientLight);
	  group.addChild(light1);
	  
	  return group;
	}

    // Create the scene background 
	private Background createBackground() {
	  Background back = new Background(new Color3f(0.9f,0.9f,0.9f));
	  BoundingSphere bound_sph = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
	  back.setApplicationBounds(bound_sph);
	  return back;
	}

    // Create the view platform
	private ViewPlatform createViewPlatform() {
	  ViewPlatform vp = new ViewPlatform();
	  vp.setViewAttachPolicy(View.RELATIVE_TO_FIELD_OF_VIEW);
	  vp.setActivationRadius(100.0f);
	  return vp;
	}

	// Create view branch group
	private BranchGroup createViewBranchGroup(ViewPlatform vp) {
	  BranchGroup viewGroup = new BranchGroup();
	  viewGroup.addChild(vp);	  
	  return viewGroup;
	}

	// Create the canvas3D and view
	private View createView(ViewPlatform vp) {
	  View view = new View();	
	  
	  // Create physical body and environment
	  PhysicalBody body =  new PhysicalBody();
	  PhysicalEnvironment env = new PhysicalEnvironment();
	  view.setPhysicalEnvironment(env);
	  view.setPhysicalBody(body);
	  if (vp != null) {
		view.attachViewPlatform(vp);
	  }
	  view.setBackClipDistance(100.0);
	  view.setFrontClipDistance(1.0);
	  
	  return view;
	}

	// Create Canvas3D
	private Canvas3D createCanvas3D(boolean offscreen) {
	  GraphicsConfigTemplate3D gc3D = new GraphicsConfigTemplate3D();
	  gc3D.setSceneAntialiasing(GraphicsConfigTemplate.PREFERRED);
	  GraphicsDevice gd[] = GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices();
	  
	  Canvas3D canvas = new Canvas3D(gd[0].getBestConfiguration(gc3D), offscreen);
      canvas.setSize(WIDTH, HEIGHT);

	  return canvas;
	}

	private void createParticles()
	{
      int size = d_partList.size();
      if (!(size > 0)) return;

      // Update domain size
      d_domainSize = d_partList.getRVESize();

      // Find particle type
      Particle part = d_partList.getParticle(0);
      int type = part.getType();

      // Draw the particles
      if (type == CIRCLE) 
        createCylinders(size);
      else if (type == SPHERE)
        createSpheres(size);
	}

    //-------------------------------------------------------------------------
    // Check if the cylinder is outside the RVE
    //-------------------------------------------------------------------------
	private void createCylinders(int size) {

      for (int ii = 0; ii < size; ii++) {

        // Get the particle data
        Particle part = d_partList.getParticle(ii);
        double radius = part.getRadius();
        Point center = part.getCenter();
        double xCent = center.getX();
        double yCent = center.getY();

        if (cylinderIsOutsideRVE(radius, xCent, yCent)) {
          continue;
        }

	    d_sceneGroup.addChild(createCylinder(radius, xCent, yCent));
      }
		
	}

    //-------------------------------------------------------------------------
    // Check if the cylinder is outside the RVE
    //-------------------------------------------------------------------------
	private boolean cylinderIsOutsideRVE(double radius, double xCent, double yCent) {
      double distXPlus = d_domainSize - (xCent-radius);
      double distYPlus = d_domainSize - (yCent-radius);
      double distXMinus = xCent+radius;
      double distYMinus = yCent+radius;
      if (distXPlus <= 0.0) return true;
      if (distYPlus <= 0.0) return true;
      if (distXMinus <= 0.0) return true;
      if (distYMinus <= 0.0) return true;
      return false;
    }

    // Actually create cylinder
	private BranchGroup createCylinder(double radius, double xcent, double ycent) {

	  // Scale radius to lie between 0 and 2 and coordinates to lie between -1 and 1
	  System.out.println("Before:  radius = "+radius+" xcent = "+xcent+" ycent = " +ycent);
	  radius = 2.0*radius/d_domainSize;
	  xcent = 2.0*xcent/d_domainSize-1.0;
	  ycent = 2.0*ycent/d_domainSize-1.0;
	  double height = 2.0;
	  //double height = d_partList.getRVESize();

	  System.out.println("After: radius = "+radius+" xcent = "+xcent+" ycent = " +ycent);
	  // Create the object root
	  BranchGroup objRoot = new BranchGroup();
	  objRoot.setCapability(BranchGroup.ALLOW_DETACH);
	  
	  // The cylinder is initially aligned with the Y axis.  Make it aligned with Z by
	  // rotating around x by 90 degrees.
	  Transform3D rotate = new Transform3D();
	  rotate.rotX(Math.PI/2.0d);
	  
	  // Next translate it the correct (x, y) position
	  Transform3D translate = new Transform3D();
	  Vector3d vector = new Vector3d(xcent, ycent, 0.0);
	  translate.setTranslation(vector);
	  
	  // Create two transform groups
	  TransformGroup objRotate = new TransformGroup(rotate);
	  TransformGroup objTranslate = new TransformGroup(translate);
	  
	  // Combine the transformations
	  objRoot.addChild(objTranslate);
	  objTranslate.addChild(objRotate);
	  
	  // Add the cylinder to translations TransformGroup
	  Appearance appear = new Appearance();
	  Color3f color = new Color3f(1.0f, 0.7f, 0.8f);
	  Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
	  //appear.setTransparencyAttributes(
	  //		  new TransparencyAttributes(TransparencyAttributes.FASTEST, 0.3f));
	  appear.setMaterial(new Material(color, black, color, black, 80.0f));
	  Cylinder cylinder = new Cylinder((float) radius, (float) height, appear);
	  objRotate.addChild(cylinder);
	  objRoot.setUserData("Cylinder");

	  return objRoot;
	}

	private void createSpheres(int size) {
		
      for (int ii = 0; ii < size; ii++) {

        // Get the particle data
        Particle part = d_partList.getParticle(ii);
        double radius = part.getRadius();
        Point center = part.getCenter();
        double xCent = center.getX();
        double yCent = center.getY();
        double zCent = center.getZ();

        if (cylinderIsOutsideRVE(radius, xCent, yCent) ||
        	cylinderIsOutsideRVE(radius, yCent, zCent) ||
        	cylinderIsOutsideRVE(radius, zCent, xCent)) {
          continue;
        }

	    d_sceneGroup.addChild(createSphere(radius, xCent, yCent, zCent));
      }
	}

	private BranchGroup createSphere(double radius, double xcent, double ycent, double zcent) {
		
	  // Scale radius and coordinates to lie between 0 and 1
	  radius = 2.0*radius/d_domainSize;
	  xcent = 2.0*xcent/d_domainSize -1.0;
	  ycent = 2.0*ycent/d_domainSize -1.0;
	  zcent = 2.0*zcent/d_domainSize -1.0;

	  BranchGroup bg = new BranchGroup();
	  bg.setCapability(BranchGroup.ALLOW_DETACH);
	  
	  Appearance appear = new Appearance();
	  Color3f color = new Color3f(1.0f, 0.7f, 0.8f);
	  Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
	  //appear.setTransparencyAttributes(
	  //		  new TransparencyAttributes(TransparencyAttributes.FASTEST, 0.3f));
	  appear.setMaterial(new Material(color, black, color, black, 80.0f));

	  Sphere sphere = new Sphere((float) radius, appear);
	  TransformGroup tg = new TransformGroup();
	  Transform3D transform = new Transform3D();
	  Vector3d vector = new Vector3d(xcent, ycent, zcent);
	  transform.setTranslation(vector);
	  tg.setTransform(transform);
	  tg.addChild(sphere);
	  bg.addChild(tg);
	  //bg.addChild(new Sphere((float) radius, appear));
	  bg.setUserData("Sphere");

	  return bg;
	}

	/*
	private TransformGroup[] getViewTransformGroupArray() {
	  TransformGroup[] array = new TransformGroup[1];
	  array[0] = new TransformGroup();
	  
	  Transform3D tMatrix = new Transform3D();
	  tMatrix.setScale(getScale());
	  tMatrix.setTranslation(new Vector3d(0.0, 0.0, -20.0));
	  tMatrix.invert();
	  array[0].setTransform(tMatrix);
	  return array;
	}

	private double getScale() {
	  return 3.0;
	}
	*/


  }

}