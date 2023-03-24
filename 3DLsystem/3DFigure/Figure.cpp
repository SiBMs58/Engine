//
// Created by Siebe Mees on 13/03/2023.
//

#include "Figure.h"


Matrix Figure::scaleFigure(const double scale){
    Matrix scaleMatrix;
    scaleMatrix(1,1) = scale;
    scaleMatrix(2,2) = scale;
    scaleMatrix(3,3) = scale;
    scaleMatrix(4,4) = 1;
    return scaleMatrix;
}

Matrix Figure::rotateX(double angle){
    // Convert to radians
    angle = angle*M_PI/180;
    Matrix rotateXMatrix;
    rotateXMatrix(1,1) = 1;
    rotateXMatrix(2,2) = cos(angle);
    rotateXMatrix(2,3) = sin(angle);
    rotateXMatrix(3,2) = -sin(angle);
    rotateXMatrix(3,3) = cos(angle);
    rotateXMatrix(4,4) = 1;
    return rotateXMatrix;
}

Matrix Figure::rotateY(double angle){
    // Convert to radians
    angle = angle*M_PI/180;
    Matrix rotateYMatrix;
    rotateYMatrix(1,1) = cos(angle);
    rotateYMatrix(1,3) = -sin(angle);
    rotateYMatrix(2,2) = 1;
    rotateYMatrix(1,3) = sin(angle);
    rotateYMatrix(3,3) = cos(angle);
    rotateYMatrix(4,4) = 1;
    return rotateYMatrix;
}

Matrix Figure::rotateZ(double angle){
    // Convert to radians
    angle = angle*M_PI/180;
    Matrix rotateZMatrix;
    rotateZMatrix(1, 1) = cos(angle);
    rotateZMatrix(1, 2) = sin(angle);
    rotateZMatrix(2, 1) = -sin(angle);
    rotateZMatrix(2, 2) = cos(angle);
    rotateZMatrix(3, 3) = 1;
    rotateZMatrix(4, 4) = 1;
    return rotateZMatrix;
}

Matrix Figure::translate(const Vector3D &vector) {
    Matrix translatieMatrix;
    translatieMatrix(1, 1) = 1;
    translatieMatrix(2, 2) = 1;
    translatieMatrix(3, 3) = 1;
    translatieMatrix(4, 4) = 1;
    translatieMatrix(4, 1) = vector.x;
    translatieMatrix(4, 2) = vector.y;
    translatieMatrix(4, 3) = vector.z;
    return translatieMatrix;
}