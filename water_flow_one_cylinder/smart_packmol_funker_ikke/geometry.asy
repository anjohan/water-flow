settings.outformat = "pdf";
settings.prc = false;
settings.render = 16;
size(16cm,16cm);
import three;
import graph3;
usepackage("mathpazo");

currentprojection=orthographic((100,30,6),up=Z);
//currentlight=light(-0.15,-3.5,1.5);
real myopacity=1;//0.801

real a = 7.16, Lx=27, Ly=27, radius=Lx/9, height=Lx/3, movingbottom=2, rigidbottom=1, totalheight = 2*movingbottom+rigidbottom+height, slabthickness=2*movingbottom+rigidbottom;

axes3("$x$","$y$","$z$",min=(-1,-1,-1),max=(Lx+3,Ly+3,totalheight+3),arrow=Arrow3());

string lengthtext(real L){
    return format("%.2g a = ",L) +  format("%.1f\ \mathrm{\AA}",a*L);
}

pen graypen = gray + opacity(0.2);

void boxandcube(triple origin=O, real scalex, real scaley, real scalez, pen p){
    draw(shift(origin)*scale(scalex, scaley, scalez)*unitcube,p);
    draw(shift(origin)*scale(scalex, scaley, scalez)*unitbox,black);
}

draw((Lx,Ly*0.75,movingbottom+0.5*rigidbottom),L=Label("Rigid $\mathrm{SiO_2}$"));
draw((Lx,Ly*0.75,0.5*movingbottom),L=Label("Moving $\mathrm{SiO_2}$"));
draw((Lx,Ly*0.75,1.5*movingbottom+rigidbottom),L=Label("Moving $\mathrm{SiO_2}$"));

boxandcube(Lx,Ly,slabthickness,graypen);
draw(shift(0,0,movingbottom)*scale(Lx, Ly, movingbottom)*unitbox,black);

draw(shift(Lx/2,Ly/2,slabthickness)*scale(radius,radius,height)*unitcylinder,graypen);
draw(shift(Lx/2,Ly/2,slabthickness)*scale(radius,radius,height)*unitdisk,graypen);
draw(shift(Lx/2,Ly/2,slabthickness+height)*scale(radius,radius,height)*unitdisk,graypen);
draw(shift(Lx/2,Ly/2,slabthickness)*scale(radius,radius,height)*unitcircle3,black);
draw(shift(Lx/2,Ly/2,slabthickness+height)*scale(radius,radius,height)*unitcircle3,black);

draw((0,0,slabthickness)--(0,0,totalheight),black+dashed);
draw((Lx,0,slabthickness)--(Lx,0,totalheight),black+dashed,L=Label("$h = "+lengthtext(height)+"$",align=W,position=0.5));
draw((0,Ly,slabthickness)--(0,Ly,totalheight),black+dashed);
draw((Lx,Ly,slabthickness)--(Lx,Ly,totalheight),black+dashed);

draw((0,0,totalheight)--(Lx,0,totalheight)--(Lx,Ly,totalheight)--(0,Ly,totalheight)--cycle,black+dashed);

draw((Lx/2,Lx/2,totalheight)--(Lx/2,Ly/2+radius,totalheight),black+dashed,L=Label("$R="+lengthtext(radius)+"$",align=2N));
draw((Lx,0,0.5*movingbottom),black,L=Label("$d = "+lengthtext(movingbottom)+"$",align=W,position=0.5));
draw((Lx,0,movingbottom+0.5*rigidbottom),black,L=Label("$d = "+lengthtext(rigidbottom)+"$",align=W,position=0.5));
draw((Lx,0,rigidbottom+1.5*movingbottom),black,L=Label("$d = "+lengthtext(movingbottom)+"$",align=W,position=0.5));

draw((Lx,Ly/2,0),black,L=Label("$L_x=L_y = "+lengthtext(Lx)+"$",align=S));

draw((Lx/2,Ly,1.25totalheight),black,L=Label(format("$a=%.2f\ \mathrm{\AA}$ (lattice parameter of $\beta$-cristobalite)",a),align=N));
