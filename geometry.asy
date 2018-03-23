settings.outformat = "pdf";
settings.prc = false;
settings.render = 4;
size(16cm,16cm);
import three;
import graph3;
usepackage("mathpazo");

currentprojection=orthographic((100,30,6),up=Z);
//currentlight=light(-0.15,-3.5,1.5);
real myopacity=1;//0.801

real a = 7.16, Lx=30, Ly=30, radius=4, height=12, bottomthickness=2, totalheight = bottomthickness*3+height;

axes3("$x$","$y$","$z$",min=(-1,-1,-1),max=(Lx+3,Ly+3,height+3*bottomthickness+3),arrow=Arrow3());

string lengthtext(real L){
    return format("%.2g a = ",L) +  format("%.1f\ \mathrm{\AA}",a*L);
}

pen graypen = gray + opacity(0.2);

void boxandcube(triple origin=O, real scalex, real scaley, real scalez, pen p){
    draw(shift(origin)*scale(scalex, scaley, scalez)*unitcube,p);
    draw(shift(origin)*scale(scalex, scaley, scalez)*unitbox,black);
}

draw((Lx,Ly*0.75,3*bottomthickness/2),L=Label("Rigid $\mathrm{SiO_2}$"));
draw((Lx,Ly*0.75,bottomthickness/2),L=Label("Moving $\mathrm{SiO_2}$"));
draw((Lx,Ly*0.75,5*bottomthickness/2),L=Label("Moving $\mathrm{SiO_2}$"));

boxandcube(Lx,Ly,3*bottomthickness,graypen);
draw(shift(0,0,bottomthickness)*scale(Lx, Ly, bottomthickness)*unitbox,black);

draw(shift(Lx/2,Ly/2,3*bottomthickness)*scale(radius,radius,height)*unitcylinder,graypen);
draw(shift(Lx/2,Ly/2,3*bottomthickness)*scale(radius,radius,height)*unitdisk,graypen);
draw(shift(Lx/2,Ly/2,3*bottomthickness+height)*scale(radius,radius,height)*unitdisk,graypen);
draw(shift(Lx/2,Ly/2,3*bottomthickness)*scale(radius,radius,height)*unitcircle3,black);
draw(shift(Lx/2,Ly/2,3*bottomthickness+height)*scale(radius,radius,height)*unitcircle3,black);

draw((0,0,3*bottomthickness)--(0,0,totalheight),black+dashed);
draw((Lx,0,3*bottomthickness)--(Lx,0,totalheight),black+dashed,L=Label("$h = "+lengthtext(height)+"$",align=W,position=0.5));
draw((0,Ly,3*bottomthickness)--(0,Ly,totalheight),black+dashed);
draw((Lx,Ly,3*bottomthickness)--(Lx,Ly,totalheight),black+dashed);

draw((0,0,totalheight)--(Lx,0,totalheight)--(Lx,Ly,totalheight)--(0,Ly,totalheight)--cycle,black+dashed);

draw((Lx/2,Lx/2,totalheight)--(Lx/2,Ly/2+radius,totalheight),black+dashed,L=Label("$R="+lengthtext(radius)+"$",align=N));
draw((Lx,0,0.5*bottomthickness),black,L=Label("$d = "+lengthtext(bottomthickness)+"$",align=W,position=0.5));
draw((Lx,0,1.5*bottomthickness),black,L=Label("$d = "+lengthtext(bottomthickness)+"$",align=W,position=0.5));
draw((Lx,0,2.5*bottomthickness),black,L=Label("$d = "+lengthtext(bottomthickness)+"$",align=W,position=0.5));

draw((Lx,Ly/2,0),black,L=Label("$L_x=L_y = "+lengthtext(Lx)+"$",align=S));

draw((Lx/2,Ly,1.25totalheight),black,L=Label(format("$a=%.2f\ \mathrm{\AA}$ (lattice parameter of $\beta$-cristobalite)",a),align=N));
