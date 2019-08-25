//███╗   ███╗ █████╗ ██╗███╗   ██╗
//████╗ ████║██╔══██╗██║████╗  ██║
//██╔████╔██║███████║██║██╔██╗ ██║
//██║╚██╔╝██║██╔══██║██║██║╚██╗██║
//██║ ╚═╝ ██║██║  ██║██║██║ ╚████║
//╚═╝     ╚═╝╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝

Polygon c,p;
float mx,my,mz;      //mouse x,y and fake "z"
PVector mouse;
PVector g;           //gravitational force
PVector forcePoint;  //point on which the spring force acts
float timeInterval;  //
float pokeForce;     //
PVector centerOfMass;//
Boolean toggleCuboid, toggleDampening; 

void setup() {
  fullScreen(P3D);
  //size(1000, 600, P3D);
  noStroke();
  
  //make the buttons
  setupButtons();
  textSize(24);
  
  //initial settings
  toggleCuboid = false;
  toggleDampening = false;
  
  //set up the polygon with desired properties
  p = new Polygon();
  p.makeObj("teapot.obj");
  p.setupObj();
  
  //setup the constants
  mz = 0;
  g = new PVector(0,100,0);
  timeInterval = 1.0/60;   
  pokeForce = 10000;  
  
  /*
  c = new Polygon();
  c.minCorner = new PVector(-width/2, -height/2, -200);
  c.maxCorner = new PVector( width/2,  height/2,  200);
  c.mass = 99999999;
  c.makeCuboid();
  */
}    



//██████████ Carry out each frame
void draw() { 
  // define the mouse vectors
  mx = mouseX - width/2;
  my = mouseY - height/2;
  mouse = new PVector(mx, my, mz);
  lights();
  background(0);
  perspective();
  
  //dampen if dampening key is pressed
  if (toggleDampening) {
    p.dampen(0.01);
    translate( width/2 - 60, 0);
    text("DAMPENING", 10, 30 );
    translate(-width/2 + 60, 0);
  }
  
  //show the numerical values we get
  //printVector(0, p.position, "Position"); //print a vector
  showMatrix(p.I, "Tensor of Inertia", 0, 30);
  showMatrix(p.ICuboid, "Tensor of Inertia Approximation (Cuboid)", 0, height-130);
  translate(width/2, height/2, 0);
  
  //choose the point on which the spring force will act
  forcePoint = p.sRotated.getChild(0).getVertex(0);
  
  //scale(1);
  p.doFrame(timeInterval);    //process the forces on the object
  if (key != 'h'&&key != 'H') //allow hiding of the shape
  p.display();                //show the object
  
  /*
  PVector newPosition = p.maxCorner.copy();
  newPosition.sub(p.minCorner);
  newPosition.div(2);
  newPosition.add(p.minCorner);
  newPosition.sub(p.position);
  newPosition.mult(-1);  
  */
  
  //show the bounding box if key pressed
  if (key == 'q' || key == 'Q'){
    p.sBoundingRotated = p.constructCuboidPShape(p.verticesRotated);
    shape(p.sBoundingRotated);
  }
  
  //draw line of spring from mouse to the vertex it acts on
  drawVector(mouse, forcePoint);
  
  //move "mouse" in or out
  if ( ((key == 'S') || (key == 's')) && keyPressed) mz+= 10;
  if ( ((key == 'W' )||( key == 'w')) && keyPressed)  mz-= 10;
  
  if (key == ESC ) exit();
}
