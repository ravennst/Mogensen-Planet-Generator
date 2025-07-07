typedef struct Vertex
{
  double h; /* altitude */
  double s; /* seed */
  double x,y,z; /* coordinates */
  double shadow; /* approximate rain shadow */
} vertex;

/* distance squared between vertices */
double dist2(vertex a, vertex b)
{
  double abx, aby, abz;
  abx = a.x-b.x; aby = a.y-b.y; abz = a.z-b.z;
  return abx*abx+aby*aby+abz*abz;
}

/* For the vertices of the tetrahedron */
vertex tetra[4];

  /* initialize vertices to slightly irregular tetrahedron */
  tetra[0].x = -sqrt(3.0)-0.20;
  tetra[0].y = -sqrt(3.0)-0.22;
  tetra[0].z = -sqrt(3.0)-0.23;

  tetra[1].x = -sqrt(3.0)-0.19;
  tetra[1].y = sqrt(3.0)+0.18;
  tetra[1].z = sqrt(3.0)+0.17;

  tetra[2].x = sqrt(3.0)+0.21;
  tetra[2].y = -sqrt(3.0)-0.24;
  tetra[2].z = sqrt(3.0)+0.15;

  tetra[3].x = sqrt(3.0)+0.24;
  tetra[3].y = sqrt(3.0)+0.22;
  tetra[3].z = -sqrt(3.0)-0.25;

  r1 = rseed;

  r1 = rand2(r1,r1);
  r2 = rand2(r1,r1);
  r3 = rand2(r1,r2);
  r4 = rand2(r2,r3);

  tetra[0].s = r1;
  tetra[1].s = r2;
  tetra[2].s = r3;
  tetra[3].s = r4;

  tetra[0].h = M;
  tetra[1].h = M;
  tetra[2].h = M;
  tetra[3].h = M;

  tetra[0].shadow = 0.0;
  tetra[1].shadow = 0.0;
  tetra[2].shadow = 0.0;
  tetra[3].shadow = 0.0;

double planet(a,b,c,d, x,y,z, level)
vertex a,b,c,d;             /* tetrahedron vertices */
double x,y,z;               /* goal point */
int level;                  /* levels to go */
{
  vertex e;
  double lab, lac, lad, lbc, lbd, lcd, maxlength;
  double es1, es2, es3;
  double eax,eay,eaz, epx,epy,epz;
  double ecx,ecy,ecz, edx,edy,edz;
  double x1,y1,z1,x2,y2,z2,l1,tmp;

  if (level>0) {

    /* make sure ab is longest edge */
    lab = dist2(a,b);
    lac = dist2(a,c);
    lad = dist2(a,d);
    lbc = dist2(b,c);
    lbd = dist2(b,d);
    lcd = dist2(c,d);

    maxlength = lab;
    if (lac > maxlength) maxlength = lac;
    if (lad > maxlength) maxlength = lad;
    if (lbc > maxlength) maxlength = lbc;
    if (lbd > maxlength) maxlength = lbd;
    if (lcd > maxlength) maxlength = lcd;

    if (lac == maxlength) return(planet(a,c,b,d, x,y,z, level));
    if (lad == maxlength) return(planet(a,d,b,c, x,y,z, level));
    if (lbc == maxlength) return(planet(b,c,a,d, x,y,z, level));
    if (lbd == maxlength) return(planet(b,d,a,c, x,y,z, level));
    if (lcd == maxlength) return(planet(c,d,a,b, x,y,z, level));

    if (level == 11) { /* save tetrahedron for caching */
      ssa = a; ssb = b; ssc = c; ssd = d;
    }

    /* ab is longest, so cut ab */
      e.s = rand2(a.s,b.s);
      es1 = rand2(e.s,e.s);
      es2 = 0.5+0.1*rand2(es1,es1);  /* find cut point */
      es3 = 1.0-es2;

      if (a.s<b.s) {
        e.x = es2*a.x+es3*b.x; e.y = es2*a.y+es3*b.y; e.z = es2*a.z+es3*b.z;
      } else if (a.s>b.s) {
        e.x = es3*a.x+es2*b.x; e.y = es3*a.y+es2*b.y; e.z = es3*a.z+es2*b.z;
      } else { /* as==bs, very unlikely to ever happen */
        e.x = 0.5*a.x+0.5*b.x; e.y = 0.5*a.y+0.5*b.y; e.z = 0.5*a.z+0.5*b.z;
      }

      /* new altitude is: */
      if (matchMap && lab > matchSize) { /* use map height */
        double l, xx, yy;
        l = sqrt(e.x*e.x+e.y*e.y+e.z*e.z);
        yy = asin(e.y/l)*23/PI+11.5;
        xx = atan2(e.x,e.z)*23.5/PI+23.5;
        e.h = cl0[(int)(xx+0.5)][(int)(yy+0.5)]*0.1/8.0;
      } else {
        if (lab>1.0) lab = pow(lab,0.5);
        /* decrease contribution for very long distances */
        e.h = 0.5*(a.h+b.h) /* average of end points */
          + e.s*dd1*pow(fabs(a.h-b.h),POWA)
          /* plus contribution for altitude diff */
          + es1*dd2*pow(lab,POW); /* plus contribution for distance */
      }

      /* calculate approximate rain shadow for new point */
      if (e.h <= 0.0 ) e.shadow = 0.0;
      else {
      x1 = 0.5*(a.x+b.x);
      x1 = a.h*(x1-a.x)+b.h*(x1-b.x);
      y1 = 0.5*(a.y+b.y);
      y1 = a.h*(y1-a.y)+b.h*(y1-b.y);
      z1 = 0.5*(a.z+b.z);
      z1 = a.h*(z1-a.z)+b.h*(z1-b.z);
      l1 = sqrt(x1*x1+y1*y1+z1*z1);
      if (l1==0.0) l1 = 1.0;
      tmp = sqrt(1.0-y*y);
      if (tmp<0.0001) tmp = 0.0001;
      x2 = x*x1+y*y1+z*z1;
      z2 = -z/tmp*x1+x/tmp*z1;
      if (lab > 0.04)
	e.shadow = (a.shadow + b.shadow- cos(PI*shade_angle/180.0)*z2/l1)/3.0;
      else
	e.shadow = (a.shadow + b.shadow)/2.0;
      }
      


      /* find out in which new tetrahedron target point is */
      eax = a.x-e.x; eay = a.y-e.y; eaz = a.z-e.z;
      ecx = c.x-e.x; ecy = c.y-e.y; ecz = c.z-e.z;
      edx = d.x-e.x; edy = d.y-e.y; edz = d.z-e.z;
      epx =   x-e.x; epy =   y-e.y; epz =   z-e.z;
      if ((eax*ecy*edz+eay*ecz*edx+eaz*ecx*edy
           -eaz*ecy*edx-eay*ecx*edz-eax*ecz*edy)*
          (epx*ecy*edz+epy*ecz*edx+epz*ecx*edy
           -epz*ecy*edx-epy*ecx*edz-epx*ecz*edy)>0.0) {
        /* point is inside acde */
        return(planet(c,d,a,e, x,y,z, level-1));
      } else {
        /* point is inside bcde */
        return(planet(c,d,b,e, x,y,z, level-1));
      }
  }

  
  else { /* level == 0 */    
    rainShadow  = 0.25*(a.shadow+b.shadow+c.shadow+d.shadow);
    return 0.25*(a.h+b.h+c.h+d.h);
  }
}