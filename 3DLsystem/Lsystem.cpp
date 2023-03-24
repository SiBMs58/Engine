//
// Created by Siebe Mees on 17/03/2023.
//

#include "Lsystem.h"

vector<Figure> Lsystem::generateFigures(const ini::Configuration &configuration) {
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    // Maak figures aan
    vector<Figure> figures;
    for (int i = 0; i < nrFigures; ++i) {
        // Maak een figuur aan
        Figure figure;
        int nrPoints = configuration["Figure"+to_string(i)]["nrPoints"].as_int_or_die();
        int nrLines = configuration["Figure"+to_string(i)]["nrLines"].as_int_or_die();

        // Point
        for (int k = 0; k < nrPoints; ++k) {
            // Maak een punt aan
            Vector3D points;
            vector<double> vectorPoints = configuration["Figure"+to_string(i)]["point" + to_string(k)].as_double_tuple_or_die();
            points.x = vectorPoints[0];
            points.y = vectorPoints[1];
            points.z = vectorPoints[2];
            figure.points.push_back(points);
        }

        // Face
        for (int j = 0; j < nrLines; ++j) {
            // Maak een Face aan
            Face face;
            face.point_indexes = configuration["Figure"+to_string(i)]["line" + to_string(j)].as_int_tuple_or_die();
            figure.faces.push_back(face);
        }

        // Color
        Color colorFigure;
        vector<double> vectorColorFigure = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_die();
        colorFigure.red = vectorColorFigure[0]*255;
        colorFigure.green = vectorColorFigure[1]*255;
        colorFigure.blue = vectorColorFigure[2]*255;

        figure.color = colorFigure;

        // Schalen
        double scale = configuration["Figure"+to_string(i)]["scale"].as_double_or_die();
        Matrix scaleMatrix = figure.scaleFigure(scale);
        // RotateX
        double rotateX = configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die();
        Matrix rotateXMatrix = figure.rotateX(rotateX);
        // RotateT
        double rotateY = configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die();
        Matrix rotateYMatrix = figure.rotateY(rotateY);
        // RotateZ
        double rotateZ = configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die();
        Matrix rotateZMatrix = figure.rotateZ(rotateZ);
        // Translatie over vector
        vector<double> vectorCenter = configuration["Figure"+to_string(i)]["center"].as_double_tuple_or_die();
        Vector3D translatieVector;
        translatieVector.x = vectorCenter[0];
        translatieVector.y = vectorCenter[1];
        translatieVector.z = vectorCenter[2];
        Matrix translationMatrix = figure.translate(translatieVector);

        Matrix transformationMatrix = scaleMatrix*rotateXMatrix*rotateYMatrix*rotateZMatrix*translationMatrix;
        applyTransformation(figure, transformationMatrix);

        // Color
        vector<double> figureColor = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_die();
        figure.color.red = figureColor[0]*255;
        figure.color.green = figureColor[1]*255;
        figure.color.blue = figureColor[2]*255;

        figures.push_back(figure);
    }
    vector<double> vectorEye = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyePoint;
    eyePoint.x = vectorEye[0];
    eyePoint.y = vectorEye[1];
    eyePoint.z = vectorEye[2];
    Matrix eyePointTransformationMatrix = eyePointTrans(eyePoint);
    applyTransformation(figures, eyePointTransformationMatrix);
    return figures;
}

void Lsystem::applyTransformation(Figure &fig, const Matrix &m) {
    for (int i = 0; i < fig.points.size(); ++i) {
        fig.points[i] = fig.points[i]*m;
    }
}

void Lsystem::applyTransformation(Figures3D &figs, const Matrix &m) {
    for (int i = 0; i < figs.size(); ++i) {
        applyTransformation(figs[i], m);
    }
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos((point.z/r));
}

Matrix Lsystem::eyePointTrans(Vector3D &eyepoint){
    double r = 0;
    double theta = 0;
    double phi =0;
    // 1. Omzetten van carthesische naar poolcoördinaten
    ::toPolar(eyepoint, theta, phi, r);
    // 2. eyePointTransformationMatrix opstellen:
    Matrix eyePointTransformationMatrix;
    eyePointTransformationMatrix(1, 1) = -sin(theta);
    eyePointTransformationMatrix(1, 2) = -cos(theta)*cos(phi);
    eyePointTransformationMatrix(1, 3) = cos(theta)*sin(phi);
    eyePointTransformationMatrix(2, 1) = cos(theta);
    eyePointTransformationMatrix(2, 2) = -sin(theta)*cos(phi);
    eyePointTransformationMatrix(2, 3) = sin(theta)*sin(phi);
    eyePointTransformationMatrix(3,2) = sin(phi);
    eyePointTransformationMatrix(3,3) = cos(phi);
    eyePointTransformationMatrix(4,3) = -r;
    eyePointTransformationMatrix(4,4) = 1;
    return eyePointTransformationMatrix;
}

Lines2D Lsystem::doProjection(const Figures3D &figures){
    // Maak lines2D aan
    Lines2D lines;
    // Loop trough all figures
    for (int i = 0; i < figures.size(); ++i) {
        // Loop trough all the faces
        for (int j = 0; j < figures[i].faces.size(); ++j) {
            for (int k = 0; k < figures[i].faces[j].point_indexes.size()-1; ++k) {
                int point_index1 = figures[i].faces[j].point_indexes[k];
                int point_index2 = figures[i].faces[j].point_indexes[k+1];
                Vector3D vectorPoint1;
                vectorPoint1.x = figures[i].points[point_index1].x;
                vectorPoint1.y = figures[i].points[point_index1].y;
                vectorPoint1.z = figures[i].points[point_index1].z;
                Vector3D vectorPoint2;
                vectorPoint2.x = figures[i].points[point_index2].x;
                vectorPoint2.y = figures[i].points[point_index2].y;
                vectorPoint2.z = figures[i].points[point_index2].z;
                Point2D point1 = generatePoint2D(vectorPoint1, 1);
                Point2D point2 = generatePoint2D(vectorPoint2, 1);
                Line2D line = Line2D(point1, point2, figures[i].color);
                lines.push_back(line);
            }
            // Add line from last point_index to first one
            int point_index1 = figures[i].faces[j].point_indexes[figures[i].faces[j].point_indexes.size()-1];
            int point_index2 = figures[i].faces[j].point_indexes[0];
            Vector3D vectorPoint1;
            vectorPoint1.x = figures[i].points[point_index1].x;
            vectorPoint1.y = figures[i].points[point_index1].y;
            vectorPoint1.z = figures[i].points[point_index1].z;
            Vector3D vectorPoint2;
            vectorPoint2.x = figures[i].points[point_index2].x;
            vectorPoint2.y = figures[i].points[point_index2].y;
            vectorPoint2.z = figures[i].points[point_index2].z;
            Point2D point1 = generatePoint2D(vectorPoint1, 1);
            Point2D point2 = generatePoint2D(vectorPoint2, 1);
            Line2D line = Line2D(point1, point2, figures[i].color);
            lines.push_back(line);
        }
    }
    return lines;
}

Point2D Lsystem::generatePoint2D(const Vector3D &vectorPoint,const double d) {
    Point2D point;
    point.x = (d*vectorPoint.x)/-vectorPoint.z;
    point.y = (d*vectorPoint.y)/-vectorPoint.z;
    return point;
}