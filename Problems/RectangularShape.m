function elements = RectangularShape(sizeVec, mainradius)

r1=mainradius;
r2=((r1*sqrt(3)-r1)*1.0001);

x1=0:r1*2:sizeVec(1);
y1=0:r1*2:sizeVec(2);
z1=0:r1*2:sizeVec(3);

x2=r1:r1*2:sizeVec(1)-r1;
y2=r1:r1*2:sizeVec(2)-r1;
z2=r1:r1*2:sizeVec(3)-r1;

xy1=combvec(x1,y1);
xyz1=combvec(xy1,z1);

xy2=combvec(x2,y2);
xyz2=combvec(xy2,z2);

elements.node=[xyz1 xyz2];

elements.prop=[ones(1,size(xyz1,2))*1 ones(1,size(xyz2,2))*2];

elements.propMat=[ones(1,size(xyz1,2))*1 ones(1,size(xyz2,2))*1];


elements.radius=[r1*ones(1,size(xyz1,2)) r2*ones(1,size(xyz2,2))];

end