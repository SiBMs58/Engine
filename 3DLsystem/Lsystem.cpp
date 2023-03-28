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
        string figureType = configuration["Figure"+to_string(i)]["type"].as_string_or_die();
        if (figureType == "Cube"){
            figure = createCube();
        } else if (figureType == "3DLSystem") {
            figure = createLSystem();
        } else if (figureType == "Tetrahedron") {
            figure = createTetrahedron();
        } else if (figureType == "Icosahedron"){
            figure = createIcosahedron();
        } else if (figureType == "Octahedron"){
            figure = createOctahedron();
        } else if (figureType == "Dodecahedron"){
            figure = createDodecahedron();
        } else if (figureType == "Cone") {
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure"+to_string(i)]["height"].as_double_or_die();
            figure = createCone(n, height);
        } else if (figureType == "Cylinder"){
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure"+to_string(i)]["height"].as_double_or_die();
            figure = createCylinder(n, height);
        } else if (figureType == "Sphere"){
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            figure = createSphere( ,n)
        }
        else if (figureType == "LineDrawing") {
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

void Lsystem::toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos((point.z/r));
}

Matrix Lsystem::eyePointTrans(Vector3D &eyepoint){
    double r = 0;
    double theta = 0;
    double phi =0;
    // 1. Omzetten van carthesische naar poolcoördinaten
    Lsystem::toPolar(eyepoint, theta, phi, r);
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

Figure Lsystem::createCube() {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D point1;
    point1.x = 1;
    point1.y = -1;
    point1.z = -1;
    figure.points.push_back(point1);
    Vector3D point2;
    point2.x = -1;
    point2.y = 1;
    point2.z = -1;
    figure.points.push_back(point2);
    Vector3D point3;
    point3.x = 1;
    point3.y = 1;
    point3.z = 1;
    figure.points.push_back(point3);
    Vector3D point4;
    point4.x = -1;
    point4.y = -1;
    point4.z = 1;
    figure.points.push_back(point4);
    Vector3D point5;
    point5.x = 1;
    point5.y = 1;
    point5.z = -1;
    figure.points.push_back(point5);
    Vector3D point6;
    point6.x = -1;
    point6.y = -1;
    point6.z = -1;
    figure.points.push_back(point6);
    Vector3D point7;
    point7.x = 1;
    point7.y = -1;
    point7.z = 1;
    figure.points.push_back(point7);
    Vector3D point8;
    point8.x = -1;
    point8.y = 1;
    point8.z = 1;
    figure.points.push_back(point8);
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(5);
    face1.point_indexes.push_back(3);
    face1.point_indexes.push_back(7);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(5);
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(8);
    face2.point_indexes.push_back(3);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face2);
    Face face3;
    face3.point_indexes.push_back(2);
    face3.point_indexes.push_back(6);
    face3.point_indexes.push_back(4);
    face3.point_indexes.push_back(8);
    for (int i = 0; i < face3.point_indexes.size(); ++i) {
        face3.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face3);
    Face face4;
    face4.point_indexes.push_back(6);
    face4.point_indexes.push_back(1);
    face4.point_indexes.push_back(7);
    face4.point_indexes.push_back(4);
    for (int i = 0; i < face4.point_indexes.size(); ++i) {
        face4.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face4);
    Face face5;
    face5.point_indexes.push_back(7);
    face5.point_indexes.push_back(3);
    face5.point_indexes.push_back(8);
    face5.point_indexes.push_back(4);
    for (int i = 0; i < face5.point_indexes.size(); ++i) {
        face5.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face5);
    Face face6;
    face6.point_indexes.push_back(1);
    face6.point_indexes.push_back(6);
    face6.point_indexes.push_back(2);
    face6.point_indexes.push_back(5);
    for (int i = 0; i < face6.point_indexes.size(); ++i) {
        face6.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face6);

    return figure;
}

Figure Lsystem::createTetrahedron() {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D point1;
    point1.x = 1;
    point1.y = -1;
    point1.z = -1;
    figure.points.push_back(point1);
    Vector3D point2;
    point2.x = -1;
    point2.y = 1;
    point2.z = -1;
    figure.points.push_back(point2);
    Vector3D point3;
    point3.x = 1;
    point3.y = 1;
    point3.z = 1;
    figure.points.push_back(point3);
    Vector3D point4;
    point4.x = -1;
    point4.y = -1;
    point4.z = 1;
    figure.points.push_back(point4);
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(3);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(4);
    face2.point_indexes.push_back(3);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face2);
    Face face3;
    face3.point_indexes.push_back(1);
    face3.point_indexes.push_back(4);
    face3.point_indexes.push_back(2);
    for (int i = 0; i < face3.point_indexes.size(); ++i) {
        face3.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face3);
    Face face4;
    face4.point_indexes.push_back(1);
    face4.point_indexes.push_back(3);
    face4.point_indexes.push_back(4);
    for (int i = 0; i < face4.point_indexes.size(); ++i) {
        face4.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face4);

    return figure;
}

Figure Lsystem::createIcosahedron() {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D point1;
    point1.x = 0;
    point1.y = 0;
    point1.z = sqrt(5)/2;
    figure.points.push_back(point1);
    for (int i = 2; i < 7; ++i) {
        Vector3D point;
        point.x = cos((i-2)*2*M_PI/5);
        point.y = sin((i-2)*2*M_PI/5);
        point.z = 0.5;
        figure.points.push_back(point);
    }
    for (int i = 7; i < 12; ++i) {
        Vector3D point;
        point.x = cos(M_PI/5+(i-7)*2*M_PI/5);
        point.y = sin(M_PI/5+(i-7)*2*M_PI/5);
        point.z = -0.5;
        figure.points.push_back(point);
    }
    Vector3D point12;
    point12.x = 0;
    point12.y = 0;
    point12.z = -sqrt(5)/2;
    figure.points.push_back(point12);
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(3);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(1);
    face2.point_indexes.push_back(3);
    face2.point_indexes.push_back(4);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face2);
    Face face3;
    face3.point_indexes.push_back(1);
    face3.point_indexes.push_back(4);
    face3.point_indexes.push_back(5);
    for (int i = 0; i < face3.point_indexes.size(); ++i) {
        face3.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face3);
    Face face4;
    face4.point_indexes.push_back(1);
    face4.point_indexes.push_back(5);
    face4.point_indexes.push_back(6);
    for (int i = 0; i < face4.point_indexes.size(); ++i) {
        face4.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face4);
    Face face5;
    face5.point_indexes.push_back(1);
    face5.point_indexes.push_back(6);
    face5.point_indexes.push_back(2);
    for (int i = 0; i < face5.point_indexes.size(); ++i) {
        face5.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face5);
    Face face6;
    face6.point_indexes.push_back(2);
    face6.point_indexes.push_back(7);
    face6.point_indexes.push_back(3);
    for (int i = 0; i < face6.point_indexes.size(); ++i) {
        face6.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face6);
    Face face7;
    face7.point_indexes.push_back(3);
    face7.point_indexes.push_back(7);
    face7.point_indexes.push_back(8);
    for (int i = 0; i < face7.point_indexes.size(); ++i) {
        face7.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face7);
    Face face8;
    face8.point_indexes.push_back(3);
    face8.point_indexes.push_back(8);
    face8.point_indexes.push_back(4);
    for (int i = 0; i < face8.point_indexes.size(); ++i) {
        face8.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face8);
    Face face9;
    face9.point_indexes.push_back(4);
    face9.point_indexes.push_back(8);
    face9.point_indexes.push_back(9);
    for (int i = 0; i < face9.point_indexes.size(); ++i) {
        face9.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face9);
    Face face10;
    face10.point_indexes.push_back(4);
    face10.point_indexes.push_back(9);
    face10.point_indexes.push_back(5);
    for (int i = 0; i < face10.point_indexes.size(); ++i) {
        face10.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face10);
    Face face11;
    face11.point_indexes.push_back(5);
    face11.point_indexes.push_back(9);
    face11.point_indexes.push_back(10);
    for (int i = 0; i < face11.point_indexes.size(); ++i) {
        face11.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face11);
    Face face12;
    face12.point_indexes.push_back(5);
    face12.point_indexes.push_back(10);
    face12.point_indexes.push_back(6);
    for (int i = 0; i < face12.point_indexes.size(); ++i) {
        face12.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face12);
    Face face13;
    face13.point_indexes.push_back(6);
    face13.point_indexes.push_back(10);
    face13.point_indexes.push_back(11);
    for (int i = 0; i < face13.point_indexes.size(); ++i) {
        face13.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face13);
    Face face14;
    face14.point_indexes.push_back(6);
    face14.point_indexes.push_back(11);
    face14.point_indexes.push_back(2);
    for (int i = 0; i < face14.point_indexes.size(); ++i) {
        face14.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face14);
    Face face15;
    face15.point_indexes.push_back(2);
    face15.point_indexes.push_back(11);
    face15.point_indexes.push_back(7);
    for (int i = 0; i < face15.point_indexes.size(); ++i) {
        face15.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face15);
    Face face16;
    face16.point_indexes.push_back(12);
    face16.point_indexes.push_back(8);
    face16.point_indexes.push_back(7);
    for (int i = 0; i < face16.point_indexes.size(); ++i) {
        face16.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face16);
    Face face17;
    face17.point_indexes.push_back(12);
    face17.point_indexes.push_back(9);
    face17.point_indexes.push_back(8);
    for (int i = 0; i < face17.point_indexes.size(); ++i) {
        face17.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face17);
    Face face18;
    face18.point_indexes.push_back(12);
    face18.point_indexes.push_back(10);
    face18.point_indexes.push_back(9);
    for (int i = 0; i < face18.point_indexes.size(); ++i) {
        face18.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face18);
    Face face19;
    face19.point_indexes.push_back(12);
    face19.point_indexes.push_back(11);
    face19.point_indexes.push_back(10);
    for (int i = 0; i < face19.point_indexes.size(); ++i) {
        face19.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face19);
    Face face20;
    face20.point_indexes.push_back(12);
    face20.point_indexes.push_back(7);
    face20.point_indexes.push_back(11);
    for (int i = 0; i < face20.point_indexes.size(); ++i) {
        face20.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face20);

    return figure;
}

Figure Lsystem::createOctahedron() {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D point1;
    point1.x = 1;
    point1.y = 0;
    point1.z = 0;
    figure.points.push_back(point1);
    Vector3D point2;
    point2.x = 0;
    point2.y = 1;
    point2.z = 0;
    figure.points.push_back(point2);
    Vector3D point3;
    point3.x = -1;
    point3.y = 0;
    point3.z = 0;
    figure.points.push_back(point3);
    Vector3D point4;
    point4.x = 0;
    point4.y = -1;
    point4.z = 0;
    figure.points.push_back(point4);
    Vector3D point5;
    point5.x = 0;
    point5.y = 0;
    point5.z = -1;
    figure.points.push_back(point5);
    Vector3D point6;
    point6.x = 0;
    point6.y = 0;
    point6.z = 1;
    figure.points.push_back(point6);
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(6);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(3);
    face2.point_indexes.push_back(6);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face2);
    Face face3;
    face3.point_indexes.push_back(3);
    face3.point_indexes.push_back(4);
    face3.point_indexes.push_back(6);
    for (int i = 0; i < face3.point_indexes.size(); ++i) {
        face3.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face3);
    Face face4;
    face4.point_indexes.push_back(4);
    face4.point_indexes.push_back(1);
    face4.point_indexes.push_back(6);
    for (int i = 0; i < face4.point_indexes.size(); ++i) {
        face4.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face4);
    Face face5;
    face5.point_indexes.push_back(2);
    face5.point_indexes.push_back(1);
    face5.point_indexes.push_back(5);
    for (int i = 0; i < face5.point_indexes.size(); ++i) {
        face5.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face5);
    Face face6;
    face6.point_indexes.push_back(3);
    face6.point_indexes.push_back(2);
    face6.point_indexes.push_back(5);
    for (int i = 0; i < face6.point_indexes.size(); ++i) {
        face6.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face6);
    Face face7;
    face7.point_indexes.push_back(4);
    face7.point_indexes.push_back(3);
    face7.point_indexes.push_back(5);
    for (int i = 0; i < face7.point_indexes.size(); ++i) {
        face7.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face7);
    Face face8;
    face8.point_indexes.push_back(1);
    face8.point_indexes.push_back(4);
    face8.point_indexes.push_back(5);
    for (int i = 0; i < face8.point_indexes.size(); ++i) {
        face8.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face8);

    return figure;
}

Figure Lsystem::createDodecahedron() {
    // Maak een figuur aan
    Figure figure;
    Figure icosahedron = createIcosahedron();
    // Punten
    for (int i = 0; i < icosahedron.faces.size(); ++i) {
        Vector3D point;
        vector<int> pointIndexes = icosahedron.faces[i].point_indexes;
        point.x = (icosahedron.points[pointIndexes[0]].x+icosahedron.points[pointIndexes[1]].x+icosahedron.points[pointIndexes[2]].x)/pointIndexes.size();
        point.y = (icosahedron.points[pointIndexes[0]].y+icosahedron.points[pointIndexes[1]].y+icosahedron.points[pointIndexes[2]].y)/pointIndexes.size();
        point.z = (icosahedron.points[pointIndexes[0]].z+icosahedron.points[pointIndexes[1]].z+icosahedron.points[pointIndexes[2]].z)/pointIndexes.size();
        figure.points.push_back(point);
    }
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(3);
    face1.point_indexes.push_back(4);
    face1.point_indexes.push_back(5);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(1);
    face2.point_indexes.push_back(6);
    face2.point_indexes.push_back(7);
    face2.point_indexes.push_back(8);
    face2.point_indexes.push_back(2);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face2);
    Face face3;
    face3.point_indexes.push_back(2);
    face3.point_indexes.push_back(8);
    face3.point_indexes.push_back(9);
    face3.point_indexes.push_back(10);
    face3.point_indexes.push_back(3);
    for (int i = 0; i < face3.point_indexes.size(); ++i) {
        face3.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face3);
    Face face4;
    face4.point_indexes.push_back(3);
    face4.point_indexes.push_back(10);
    face4.point_indexes.push_back(11);
    face4.point_indexes.push_back(12);
    face4.point_indexes.push_back(4);
    for (int i = 0; i < face4.point_indexes.size(); ++i) {
        face4.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face4);
    Face face5;
    face5.point_indexes.push_back(4);
    face5.point_indexes.push_back(12);
    face5.point_indexes.push_back(13);
    face5.point_indexes.push_back(14);
    face5.point_indexes.push_back(5);
    for (int i = 0; i < face5.point_indexes.size(); ++i) {
        face5.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face5);
    Face face6;
    face6.point_indexes.push_back(5);
    face6.point_indexes.push_back(14);
    face6.point_indexes.push_back(15);
    face6.point_indexes.push_back(6);
    face6.point_indexes.push_back(1);
    for (int i = 0; i < face6.point_indexes.size(); ++i) {
        face6.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face6);
    Face face7;
    face7.point_indexes.push_back(20);
    face7.point_indexes.push_back(19);
    face7.point_indexes.push_back(18);
    face7.point_indexes.push_back(17);
    face7.point_indexes.push_back(16);
    for (int i = 0; i < face7.point_indexes.size(); ++i) {
        face7.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face7);
    Face face8;
    face8.point_indexes.push_back(20);
    face8.point_indexes.push_back(15);
    face8.point_indexes.push_back(14);
    face8.point_indexes.push_back(13);
    face8.point_indexes.push_back(19);
    for (int i = 0; i < face8.point_indexes.size(); ++i) {
        face8.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face8);
    Face face9;
    face9.point_indexes.push_back(19);
    face9.point_indexes.push_back(13);
    face9.point_indexes.push_back(12);
    face9.point_indexes.push_back(11);
    face9.point_indexes.push_back(18);
    for (int i = 0; i < face9.point_indexes.size(); ++i) {
        face9.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face9);
    Face face10;
    face10.point_indexes.push_back(18);
    face10.point_indexes.push_back(11);
    face10.point_indexes.push_back(10);
    face10.point_indexes.push_back(9);
    face10.point_indexes.push_back(17);
    for (int i = 0; i < face10.point_indexes.size(); ++i) {
        face10.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face10);
    Face face11;
    face11.point_indexes.push_back(17);
    face11.point_indexes.push_back(9);
    face11.point_indexes.push_back(8);
    face11.point_indexes.push_back(7);
    face11.point_indexes.push_back(16);
    for (int i = 0; i < face11.point_indexes.size(); ++i) {
        face11.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face11);
    Face face12;
    face12.point_indexes.push_back(16);
    face12.point_indexes.push_back(7);
    face12.point_indexes.push_back(6);
    face12.point_indexes.push_back(15);
    face12.point_indexes.push_back(20);
    for (int i = 0; i < face12.point_indexes.size(); ++i) {
        face12.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face12);

    return figure;
}

Figure Lsystem::createCone(const int n, const double h) {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D pointN;
    pointN.x = 0;
    pointN.y = 0;
    pointN.z = h;
    figure.points.push_back(pointN);
    for (int i = 0; i < n+1; ++i) {
        Vector3D point;
        point.x = cos((2*i*M_PI)/n);
        point.y = sin((2*i*M_PI)/n);
        point.z = 0;
        figure.points.push_back(point);
    }
    // Faces
    for (int i = 0; i < n+1; ++i) {
        Face face;
        face.point_indexes.push_back(1);
        face.point_indexes.push_back(i+1);
        face.point_indexes.push_back(i+2);
        for (int i = 0; i < face.point_indexes.size(); ++i) {
            face.point_indexes[i] -= 1;
        }
        figure.faces.push_back(face);
    }

    return figure;
}

Figure Lsystem::createCylinder(const int n, const double h) {
    // Maak een figuur aan
    Figure figure;
    // Punten
    for (int i = 0; i < n+1; ++i) {
        Vector3D pointBottom;
        pointBottom.x = cos((2*i*M_PI)/n);
        pointBottom.y = sin((2*i*M_PI)/n);
        pointBottom.z = 0;
        figure.points.push_back(pointBottom);
    }
    for (int i = 0; i < n+1; ++i) {
        Vector3D pointTop;
        pointTop.x = cos((2*i*M_PI)/n);
        pointTop.y = sin((2*i*M_PI)/n);
        pointTop.z = h;
        figure.points.push_back(pointTop);
    }
    // Faces
    for (int i = 0; i < n; ++i) {
        Face face;
        face.point_indexes.push_back(i);
        face.point_indexes.push_back(i+1);
        face.point_indexes.push_back(n+i+2);
        face.point_indexes.push_back(n+i+1);
        figure.faces.push_back(face);
    }

    return figure;
}

Figure Lsystem::createSphere(const double radius, const int n) {

}

Figure Lsystem::createLSystem() {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Vector3D point1;
    point1.x = 1;
    point1.y = 0;
    point1.z = 0;
    figure.points.push_back(point1);
    Vector3D point2;
    point2.x = 0;
    point2.y = 1;
    point2.z = 0;
    figure.points.push_back(point2);
    Vector3D point3;
    point3.x = -1;
    point3.y = 0;
    point3.z = 0;
    figure.points.push_back(point3);
    // Faces
    Face face1;
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(3);
    for (int i = 0; i < face1.point_indexes.size(); ++i) {
        face1.point_indexes[i] -= 1;
    }
    figure.faces.push_back(face1);
    Face face2;
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(3);
    face2.point_indexes.push_back(1);
    for (int i = 0; i < face2.point_indexes.size(); ++i) {
        face2.point_indexes[i] -= 1;
    }

    return figure;
}