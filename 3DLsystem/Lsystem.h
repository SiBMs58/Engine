//
// Created by Siebe Mees on 17/03/2023.
//

#ifndef ENGINE_LSYSTEM_H
#define ENGINE_LSYSTEM_H

#include "../easy_image.h"
#include "3DFigure/Figure.h"
#include "../ini_configuration.h"
#include "../2DLsystem/2DLine/Line2D.h"
#include "../l_parser.h"
#include "../2DLsystem/ZBuffer.h"
#include "../2DLsystem/Lsystem2D.h"
#include <stack>

using Lines2D = vector<Line2D>;

class Lsystem {
public:
    void applyDrawFunctions(Figure &figure, const double scale, const double rotateX, const double rotateY, const double rotateZ, const vector<double> vectorCenter);
    vector<Figure> generateFigures(const ini::Configuration &configuration);
    void applyTransformation(Figure &fig, const Matrix &m);
    void applyTransformation(Figures3D &figs, const Matrix &m);
    Matrix eyePointTrans(Vector3D &eyepoint);
    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    Lines2D doProjection(const Figures3D &figures);
    Point2D generatePoint2D(const Vector3D &point, const double d);

    Figure drawLSystem(const LParser::LSystem3D &l_system);

    vector<Face> triangulate(const Face& face);
    img::EasyImage drawZbufTriangles(vector<Figure>& figures, const int size, const vector<double> &backgroundColor);

    void generateFractal(Figure& fig, Figures3D& fractal, const int nr_iterations, const double scale);

private:
    stack<Vector3D> positionStack;
    stack<vector<Vector3D>> angleStack;
};


#endif //ENGINE_LSYSTEM_H