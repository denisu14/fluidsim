int get1DIndex(int x, int y)
{
  x = constrain(x, 0, R-1);
  y = constrain(y, 0, R-1);
  return x + (y*R);
}

class Fluid_D
{
  int size;
  float ts;
  float diff;
  float visc;
  
  float[] s;
  float[] dnsty;
  
  float[] Vx;
  float[] Vy;
  
  float[] Vxp;
  float[] Vyp;
  
  
  Fluid_D(float timestamp, float diffusion, float viscosity)
  {
    this.size = R;
    this.ts = timestamp;
    this.diff = diffusion;
    this.visc = viscosity;
    
    this.s = new float[R*R];
    this.dnsty = new float[R*R];
    
    this.Vx = new float[R*R];
    this.Vy = new float[R*R];
    
    this.Vxp = new float[R*R];
    this.Vyp = new float[R*R];
  }
  
  void createDyeDensity(int x, int y, float amount)
  {
    int index = get1DIndex(x,y);
    this.dnsty[index] += amount;
  }
  
  void addVel(int x, int y, float amountX, float amountY)
  {
    int index = get1DIndex(x,y);
    this.Vx[index] += amountX;
    this.Vy[index] += amountY;
  }
  
  void t_step() {
    float visc     = this.visc;
    float diff     = this.diff;
    float dt       = this.ts;
    float[] Vx      = this.Vx;
    float[] Vy      = this.Vy;
    float[] Vxp     = this.Vxp;
    float[] Vyp     = this.Vyp;
    float[] s       = this.s;
    float[] density = this.dnsty;

    diffuse(1, Vxp, Vx, visc, dt);
    diffuse(2, Vyp, Vy, visc, dt);

    project(Vxp, Vyp, Vx, Vy);

    advect(1, Vx, Vxp, Vxp, Vyp, dt);
    advect(2, Vy, Vyp, Vxp, Vyp, dt);

    project(Vx, Vy, Vxp, Vyp);

    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
  }
  void renderD()
  {
    for(int i = 0; i < R; i++){
      for(int j = 0; j < R; j++){
        float x = i*uP;
        float y = j*uP;
        float d = this.dnsty[get1DIndex(i,j)];
        fill(255,d);
        noStroke();
        square(x,y,uP);
    }
  }
 }
}



void diffuse (int b, float[] x, float[] xp, float diff, float ts)
{
    float a = ts * diff * (R - 2) * (R - 2);
    lin_solve(b, x, xp, a, 1 + 6 * a);
}

void lin_solve(int b, float[] x, float[] xp, float a, float c)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < R - 1; j++) {
                for (int i = 1; i < R - 1; i++) {
                    x[get1DIndex(i, j)] =
                        (xp[get1DIndex(i, j)]
                            + a*(    x[get1DIndex(i+1, j)]
                                    +x[get1DIndex(i-1, j)]
                                    +x[get1DIndex(i  , j+1)]
                                    +x[get1DIndex(i  , j-1)]
                                    +x[get1DIndex(i  , j)]
                                    +x[get1DIndex(i  , j)]
                           )) * cRecip;
                }
            }
        set_bnd(b, x);
    }
}

void project(float[] velocX, float[] velocY, float []p, float []div)
{
        for (int j = 1; j < R - 1; j++) {
            for (int i = 1; i < R - 1; i++) {
                div[get1DIndex(i, j)] = -0.5f*(
                         velocX[get1DIndex(i+1, j)]
                        -velocX[get1DIndex(i-1, j)]
                        +velocY[get1DIndex(i  , j+1)]
                        -velocY[get1DIndex(i  , j-1)]
                    )/R;
                p[get1DIndex(i, j)] = 0;
            }
    }
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    
        for (int j = 1; j < R - 1; j++) {
            for (int i = 1; i < R - 1; i++) {
                velocX[get1DIndex(i, j)] -= 0.5f * (  p[get1DIndex(i+1, j)]
                                                -p[get1DIndex(i-1, j)]) * R;
                velocY[get1DIndex(i, j)] -= 0.5f * (  p[get1DIndex(i, j+1)]
                                                -p[get1DIndex(i, j-1)]) * R;
            }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);

}

void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY, float dt)
{
    float i0, i1, j0, j1;
    
    float dtx = dt * (R - 2);
    float dty = dt  * (R - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Rfloat = R;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < R - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < R - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[get1DIndex(i, j)];
                tmp2 = dty * velocY[get1DIndex(i, j)];
                
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Rfloat + 0.5f) x = Rfloat + 0.5f; 
                i0 = floor(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Rfloat + 0.5f) y = Rfloat + 0.5f; 
                j0 = floor(y);
                j1 = j0 + 1.0f; 
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                
                int i0i = int(i0);
                int i1i = int(i1);
                int j0i = int(j0);
                int j1i = int(j1);
                
                d[get1DIndex(i, j)] = 
        s0 * (t0 * d0[get1DIndex(i0i, j0i)] + t1 * d0[get1DIndex(i0i, j1i)]) + s1 * (t0 * d0[get1DIndex(i1i, j0i)] + t1 * d0[get1DIndex(i1i, j1i)]);

        }
    }
    set_bnd(b, d);
}

void set_bnd(int b, float[] x)
{

        for(int i = 1; i < R - 1; i++) {
            x[get1DIndex(i, 0)] = b == 2 ? -x[get1DIndex(i, 1)] : x[get1DIndex(i, 1)];
            x[get1DIndex(i, R-1)] = b == 2 ? -x[get1DIndex(i, R-2)] : x[get1DIndex(i, R-2)];
        }

        for(int j = 1; j < R - 1; j++) {
            x[get1DIndex(0  , j)] = b == 1 ? -x[get1DIndex(1  , j)] : x[get1DIndex(1  , j)];
            x[get1DIndex(R-1, j)] = b == 1 ? -x[get1DIndex(R-2, j)] : x[get1DIndex(R-2, j)];
        }
    
  x[get1DIndex(0, 0)] = 0.5f * (x[get1DIndex(1, 0)] + x[get1DIndex(0, 1)]);
  x[get1DIndex(0, R-1)] = 0.5f * (x[get1DIndex(1, R-1)] + x[get1DIndex(0, R-2)]);
  x[get1DIndex(R-1, 0)] = 0.5f * (x[get1DIndex(R-2, 0)] + x[get1DIndex(R-1, 1)]);
  x[get1DIndex(R-1, R-1)] = 0.5f * (x[get1DIndex(R-2, R-1)] + x[get1DIndex(R-1, R-2)]);
}
