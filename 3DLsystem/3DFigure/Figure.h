//
// Created by Siebe Mees on 13/03/2023.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include "../../2DLsystem/2DLine/Color.h"
#include "Face.h"
#include "../vector3d.h"
#include <cmath>

using namespace std;

class Figure {
public:
    vector<Vector3D> points;
    vector<Face> faces;
    Color color;

    Matrix scaleFigure(const double);
    Matrix rotateX(const double);
    Matrix rotateY(const double);
    Matrix rotateZ(const double);
    Matrix translate(const Vector3D&);
};

typedef vector<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
