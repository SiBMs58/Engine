//
// Created by Siebe Mees on 02/04/2023.
//

#ifndef ENGINE_WIREFRAME3D_H
#define ENGINE_WIREFRAME3D_H

#include "3DFigure/Figure.h"

class Wireframe3D {
public:
    Figure createCube();
    Figure createTetrahedron();
    Figure createIcosahedron();
    Figure createOctahedron();
    Figure createDodecahedron();
    Figure createCone(const int n, const double h);
    Figure createCylinder(const int n, const double h);
    Figure createSphere(const double radius, const int n);
    Figure createTorus(const double r,const double R, const int n, const int m);

    Figure createBuckeyBall();
    Figure createMengerSponge();
};


#endif //ENGINE_WIREFRAME3D_H
