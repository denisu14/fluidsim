
final int R = 360;
final int iter = 1;
final int uP = 4;

Fluid_D fluid;


void settings()
{
  size(R,R);
}
void setup()
{

  fluid = new Fluid_D(0.1,0,0);
}

void draw() 
{
  background(0);
  fluid.t_step();
  fluid.renderD();
}

void mouseDragged()
{
  fluid.createDyeDensity(mouseX/uP, mouseY/uP, random(100,300));
  float mvx = mouseX - pmouseX;
  float mvy = mouseY - pmouseY;
  fluid.addVel(mouseX/uP, mouseY/uP, mvx, mvy);

}
  
