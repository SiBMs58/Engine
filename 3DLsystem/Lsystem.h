//
// Created by Siebe Mees on 17/03/2023.
//

#ifndef ENGINE_LSYSTEM_H
#define ENGINE_LSYSTEM_H

#include "../easy_image.h"
#include "3DFigure/Figure.h"
#include "../ini_configuration.h"
#include "../2DLsystem/2DLine/Line2D.h"

using Lines2D = vector<Line2D>;

class Lsystem {
public:
    vector<Figure> generateFigures(const ini::Configuration &configuration);
    void applyTransformation(Figure &fig, const Matrix &m);
    void applyTransformation(Figures3D &figs, const Matrix &m);
    Matrix eyePointTrans(Vector3D &eyepoint);
    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    Lines2D doProjection(const Figures3D &figures);
    Point2D generatePoint2D(const Vector3D &point,const double d);
};


#endif //ENGINE_LSYSTEM_H
