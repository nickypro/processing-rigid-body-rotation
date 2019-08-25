import g4p_controls.*;

GButton importObjButton, useCuboidalButton, useIntegratedButton,
        sizeHalf, sizeDouble, pokeForceHalf, pokeForceDouble,
        massHalf, massDouble, springConstantHalf, springConstantDouble,
        resetButton;
        
void setupButtons() {
  G4P.setCtrlMode(GControlMode.CENTER);
  importObjButton      = new GButton(this, width-70 , 30,  128, 42, "Import Object");
  useCuboidalButton    = new GButton(this, width-70 , 80,  128, 42, "Use Cuboid Tensor");
  useIntegratedButton  = new GButton(this, width-70 , 130, 128, 42, "Use Integrated Tensor");
  sizeHalf             = new GButton(this, width-90 , 175, 85, 42, "size x0.5");
  sizeDouble           = new GButton(this, width-25 , 175, 40, 42, "x2");
  pokeForceHalf        = new GButton(this, width-90 , 220, 85, 42, "poke    x0.5");
  pokeForceDouble      = new GButton(this, width-25 , 220, 40, 42, "x2");
  massHalf             = new GButton(this, width-90 , 265, 85, 42, "mass   x0.5");
  massDouble           = new GButton(this, width-25 , 265, 40, 42, "x2");
  springConstantHalf   = new GButton(this, width-90 , 310, 85, 42, "springiness   x0.5");
  springConstantDouble = new GButton(this, width-25 , 310, 40, 42, "x2");
  resetButton          = new GButton(this, width-70 , 355, 128, 42, "RESET (R)");
}


//██████████ 
void keyPressed() {
  if (key == ' ') toggleDampening = !toggleDampening;
  
  if (key == '0' ) {
    mz = 0;
    timeInterval = 1.0/60;
  }
  //change time interval between frames to be more or less slow
  if ((key == 'A') || (key == 'a')) timeInterval *= 2;
  if ((key == 'D' )||( key == 'd')) timeInterval /= 2;
  
  //reset the position/velocity/orientaation of the body
  if ((key == 'r')||(key == 'R'))   p.resetPosition();
} 

void mousePressed() {
  //poke the object if close to it
  PVector mousePos = new PVector(mx,my,p.position.z);
  mousePos.sub(p.position);
  if(mousePos.mag()<p.maxCorner.mag())
    p.applyForce(new PVector(mx,my,p.position.z), new PVector(0,0, -pokeForce), 1.0/60);
}

//██████████ Button clicks
public void handleButtonEvents(GButton button, GEvent event) {
if (event == GEvent.CLICKED) { 
  
  //button click prompts import of a new object
  if(button == importObjButton) selectInput("Select an OBJ or SVG file", "fileSelected");
  
  //button click causes approximate cuboid tensor to be used
  if(button == useCuboidalButton){
    p.Iinv = matrixClone(p.ICuboid);
    p.Iinv.invert();
  }
  
  //button click causes integrated tensor to be used
  if(button == useIntegratedButton){
    p.Iinv = matrixClone(p.I);
    p.Iinv.invert();
  }
  
  if(button == sizeHalf) p.scaleShape(0.5);
  if(button == sizeDouble) p.scaleShape(2);
  
  if(button == pokeForceHalf) pokeForce /=2;
  if(button == pokeForceDouble) pokeForce *= 2;
  
  if(button == massHalf) p.scaleMass(0.5);
  if(button == massDouble) p.scaleMass(2);
  
  if(button == springConstantHalf) p.springConstant *= 0.5;
  if(button == springConstantDouble) p.springConstant *= 2;
  
  if(button == resetButton){
    p.resetPosition();
  }
}
} 
  
//██████████
void fileSelected(File selection) { 
  //load up the shape (stop existing to stop rendering while object not loaded)
  p.exists = false; 
  
  p.s = loadShape(selection.getAbsolutePath());
  p.sRotated = loadShape(selection.getAbsolutePath());
  p.setupObj();
  
  p.exists = true;
} 
