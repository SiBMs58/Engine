//
// Created by Siebe Mees on 17/03/2023.
//

#ifndef ENGINE_LSYSTEM_H
#define ENGINE_LSYSTEM_H

#include "../easy_image.h"
#include "3DFigure/Figure.h"
#include "../ini_configuration.h"


class Lsystem {
public:
    vector<Figure> generateFigures(const ini::Configuration &configuration);
    void applyTransformation(Figure &fig, const Matrix &m);
    Matrix eyePointTrans(const Vector3D &eyepoint);
    void toPolar(const Vector3D &point);
};


#endif //ENGINE_LSYSTEM_H
