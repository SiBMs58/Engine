//
// Created by Siebe Mees on 17/03/2023.
//

#include <fstream>
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
        } else if (figureType == "3DLSystem"){
            // Maak een parser aan
            LParser::LSystem3D LParser3D;
            // 3DLSystem
            string inputFile = configuration["Figure" + to_string(i)]["inputfile"].as_string_or_die();
            // Creëer if stream
            ifstream file(inputFile);
            // Steek in LParser
            file >> LParser3D;
            figure = drawLSystem(LParser3D);
        }
        else if (figureType == "Tetrahedron") {
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
            figure = createSphere(1 ,n);
        }
        else if (figureType == "LineDrawing") {
            int nrPoints = configuration["Figure"+to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure"+to_string(i)]["nrLines"].as_int_or_die();
            // Point
            for (int k = 0; k < nrPoints; ++k) {
                vector<double> vectorPoints = configuration["Figure"+to_string(i)]["point" + to_string(k)].as_double_tuple_or_die();
                // Maak een punt aan
                Vector3D points = Vector3D::point(vectorPoints[0], vectorPoints[1], vectorPoints[2]);
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
        Vector3D translatieVector = Vector3D::point(vectorCenter[0], vectorCenter[1], vectorCenter[2]);
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
    Vector3D eyePoint = Vector3D::point(vectorEye[0], vectorEye[1], vectorEye[2]);
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
                Vector3D vectorPoint1 = Vector3D::point(figures[i].points[point_index1].x, figures[i].points[point_index1].y, figures[i].points[point_index1].z);
                Vector3D vectorPoint2 = Vector3D::point(figures[i].points[point_index2].x, figures[i].points[point_index2].y, figures[i].points[point_index2].z);
                Point2D point1 = generatePoint2D(vectorPoint1, 1);
                Point2D point2 = generatePoint2D(vectorPoint2, 1);
                Line2D line = Line2D(point1, point2, figures[i].color);
                lines.push_back(line);
            }
            // Add line from last point_index to first one
            int point_index1 = figures[i].faces[j].point_indexes[figures[i].faces[j].point_indexes.size()-1];
            int point_index2 = figures[i].faces[j].point_indexes[0];
            Vector3D vectorPoint1 = Vector3D::point(figures[i].points[point_index1].x, figures[i].points[point_index1].y, figures[i].points[point_index1].z);
            Vector3D vectorPoint2 = Vector3D::point(figures[i].points[point_index2].x, figures[i].points[point_index2].y, figures[i].points[point_index2].z);
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
    Vector3D point1 = Vector3D::point(1, -1, -1);
    figure.points.push_back(point1);
    Vector3D point2 = Vector3D::point(-1, 1, -1);
    figure.points.push_back(point2);
    Vector3D point3 = Vector3D::point(1, 1, 1);
    figure.points.push_back(point3);
    Vector3D point4 = Vector3D::point(-1, -1, 1);
    figure.points.push_back(point4);
    Vector3D point5 = Vector3D::point(1, 1, -1);
    figure.points.push_back(point5);
    Vector3D point6 = Vector3D::point(-1, -1, -1);
    figure.points.push_back(point6);
    Vector3D point7 = Vector3D::point(1, -1, 1);
    figure.points.push_back(point7);
    Vector3D point8 = Vector3D::point(-1, 1, 1);
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
    Vector3D point1 = Vector3D::point(1, -1, -1);
    figure.points.push_back(point1);
    Vector3D point2 = Vector3D::point(-1, 1, -1);
    figure.points.push_back(point2);
    Vector3D point3 = Vector3D::point(1, 1, 1);
    figure.points.push_back(point3);
    Vector3D point4 = Vector3D::point(-1, -1, 1);
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
    Vector3D point1 = Vector3D::point(0, 0, sqrt(5)/2);
    figure.points.push_back(point1);
    for (int i = 2; i < 7; ++i) {
        Vector3D point = Vector3D::point(cos((i-2)*2*M_PI/5), sin((i-2)*2*M_PI/5), 0.5);
        figure.points.push_back(point);
    }
    for (int i = 7; i < 12; ++i) {
        Vector3D point = Vector3D::point(cos(M_PI/5+(i-7)*2*M_PI/5), sin(M_PI/5+(i-7)*2*M_PI/5), -0.5);
        figure.points.push_back(point);
    }
    Vector3D point12 = Vector3D::point(0, 0, -sqrt(5)/2);
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
    Vector3D point1 = Vector3D::point(1,0,0);
    figure.points.push_back(point1);
    Vector3D point2 = Vector3D::point(0,1,0);
    figure.points.push_back(point2);
    Vector3D point3 = Vector3D::point(-1,0,0);
    figure.points.push_back(point3);
    Vector3D point4 = Vector3D::point(0,-1,0);
    figure.points.push_back(point4);
    Vector3D point5 = Vector3D::point(0,0,-1);
    figure.points.push_back(point5);
    Vector3D point6 = Vector3D::point(0,0,1);
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
        vector<int> pointIndexes = icosahedron.faces[i].point_indexes;
        Vector3D point = Vector3D::point(
                (icosahedron.points[pointIndexes[0]].x+icosahedron.points[pointIndexes[1]].x+icosahedron.points[pointIndexes[2]].x)/pointIndexes.size(),
                (icosahedron.points[pointIndexes[0]].y+icosahedron.points[pointIndexes[1]].y+icosahedron.points[pointIndexes[2]].y)/pointIndexes.size(),
                (icosahedron.points[pointIndexes[0]].z+icosahedron.points[pointIndexes[1]].z+icosahedron.points[pointIndexes[2]].z)/pointIndexes.size()
                );
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
    Vector3D pointN = Vector3D::point(0,0,h);
    figure.points.push_back(pointN);
    for (int i = 0; i < n+1; ++i) {
        Vector3D point = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), 0);
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
        Vector3D pointBottom = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), 0);
        figure.points.push_back(pointBottom);
    }
    for (int i = 0; i < n+1; ++i) {
        Vector3D pointTop = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), h);
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
    // TODO
    // Maak een figuur aan
    Figure figure;
    // 1. Genereer een icosahedron
    figure = createIcosahedron();
    // 2. Deel elke driehoek op in kleinere driehoeken
    int counter = 0;
    // 3. Herhaal stap 2 in totaal ‘n’ keer
    while (counter < n) {
        for (int j = 0; j < figure.faces.size(); ++j) {
            vector<int> pointIndexes = figure.faces[j].point_indexes;
            // Punten
            Vector3D pointA;
            pointA = figure.points[pointIndexes[0]];
            figure.points.push_back(pointA);
            Vector3D pointB;
            pointB = figure.points[pointIndexes[1]];
            figure.points.push_back(pointB);
            Vector3D pointC;
            pointC = figure.points[pointIndexes[2]];
            figure.points.push_back(pointC);
            Vector3D pointD;
            pointD = (pointA+pointB)/2;
            figure.points.push_back(pointD);
            Vector3D pointE;
            pointE = (pointA+pointC)/2;
            figure.points.push_back(pointE);
            Vector3D pointF;
            pointF = (pointB+pointC)/2;
            figure.points.push_back(pointF);
            // Faces
            Face face1; // ADE
            face1.point_indexes.push_back(figure.faces.size()+1);
            face1.point_indexes.push_back(figure.faces.size()+4);
            face1.point_indexes.push_back(figure.faces.size()+5);
            for (int i = 0; i < face1.point_indexes.size(); ++i) {
                face1.point_indexes[i] -= 1;
            }
            figure.faces.push_back(face1);
            Face face2; // BFD
            face2.point_indexes.push_back(2+figure.points.size());
            face2.point_indexes.push_back(6+figure.points.size());
            face2.point_indexes.push_back(4+figure.points.size());
            for (int i = 0; i < face2.point_indexes.size(); ++i) {
                face2.point_indexes[i] -= 1;
            }
            figure.faces.push_back(face2);
            Face face3; // CEF
            face3.point_indexes.push_back(3+figure.points.size());
            face3.point_indexes.push_back(5+figure.points.size());
            face3.point_indexes.push_back(6+figure.points.size());
            for (int i = 0; i < face3.point_indexes.size(); ++i) {
                face3.point_indexes[i] -= 1;
            }
            figure.faces.push_back(face3);
            Face face4; // DFE
            face4.point_indexes.push_back(4+figure.points.size());
            face4.point_indexes.push_back(6+figure.points.size());
            face4.point_indexes.push_back(5+figure.points.size());
            for (int i = 0; i < face4.point_indexes.size(); ++i) {
                face4.point_indexes[i] -= 1;
            }
            figure.faces.push_back(face4);

        }
        counter++;
    }
    // 4. Herschaal alle punten
    for (int i = 0; i < figure.points.size(); ++i) {
        double r = sqrt(pow(figure.points[i].x,2)+pow(figure.points[i].y,2)+pow(figure.points[i].z,2));
        figure.points[i].x = figure.points[i].x/r;
        figure.points[i].y = figure.points[i].y/r;
        figure.points[i].z = figure.points[i].z/r;
    }

    return figure;

}

string getReplacementRule(const LParser::LSystem3D &Lsystem) {
    string initiator = Lsystem.get_initiator();
    unsigned int nrIterations = Lsystem.get_nr_iterations();
    string replacementRule="";
    if (nrIterations > 0){
        // Loop nr of iterations
        for (unsigned int i = 0; i < nrIterations; ++i) {
            replacementRule ="";
            for (unsigned int i = 0; i < initiator.size(); ++i) {
                if (initiator[i] == '-' || initiator[i] == '+' || initiator[i] == ')' || initiator[i] == '(' || initiator[i] == '^' || initiator[i] == '&' || initiator[i] == 92 || initiator[i] == '/' || initiator[i] == '|'){
                    replacementRule += initiator[i];
                } else {
                    replacementRule += Lsystem.get_replacement(initiator[i]);
                }
            }
            initiator = replacementRule;
        }
    } else {
        for (unsigned int i = 0; i < initiator.size(); ++i) {
            if (initiator[i] == '-' || initiator[i] == '+' || initiator[i] == ')' || initiator[i] == '(' || initiator[i] == '^' || initiator[i] == '&' || initiator[i] == 92 || initiator[i] == '/' || initiator[i] == '|'){
                replacementRule += initiator[i];
            } else {
                replacementRule += Lsystem.get_replacement(initiator[i]);
            }
        }
    }
    return replacementRule;
}

Figure Lsystem::drawLSystem(const LParser::LSystem3D &l_system) {
    // Maak een figuur aan
    Figure figure;
    // Convert to radians
    set<char> alphabet = l_system.get_alphabet();
    double angle = l_system.get_angle();
    // Convert to radians
    angle = angle*M_PI/180;
    // Get total sting based on initiator
    string replacementRule = getReplacementRule(l_system);
    // 1. We starten in een willekeurige positie, typisch positie (x, y, z) = (0, 0, 0) en in de richting ↵ = ↵0 radialen.
    Vector3D startingPoint = Vector3D::point(0,0,0);
    Vector3D currentPosition = startingPoint;
    vector<Vector3D> startingAngle;
    Vector3D H = Vector3D::point(1,0,0);
    startingAngle.push_back(H);
    Vector3D L = Vector3D::point(0,1,0);
    startingAngle.push_back(L);
    Vector3D U = Vector3D::point(0,0,1);
    startingAngle.push_back(U);
    vector<Vector3D> currentAngle = startingAngle;
    // 2. We overlopen de string S = s1s2 . . . sk van links naar rechts en voeren afhankelijk van het huidige symbool si volgende actie uit:
    for (int i = 0; i < replacementRule.size(); ++i) {
        // Si = +
        if (replacementRule[i] == '+') {
            H = currentAngle[0];
            L = currentAngle[1];
            // Hnew = H cos(angle) + L sin(angle)
            currentAngle[0] = H*cos(angle) + L*sin(angle);
            // Lnew = - H sin(angle) - L cos(angle)
            currentAngle[1] = -H*sin(angle) + L*cos(angle);
            // Unew = U
        }
        // Si = -
        else if (replacementRule[i] == '-') {
            H = currentAngle[0];
            L = currentAngle[1];
            // Hnew = H cos(-angle) + L sin(-angle)
            currentAngle[0] = H*cos(-angle) + L*sin(-angle);
            // Lnew = - H sin(-angle) - L cos(angle)
            currentAngle[1] = -H*sin(-angle) + L*cos(-angle);
            // Unew = U
        }
        // Si = (
        else if (replacementRule[i] == '('){
            positionStack.push(currentPosition);
            angleStack.push(currentAngle);
        }
        // Si == )
        else if (replacementRule[i] == ')'){
            currentPosition = positionStack.top();
            positionStack.pop();
            currentAngle = angleStack.top();
            angleStack.pop();
        }
        // Si == ^
        else if (replacementRule[i] == '^'){
            H = currentAngle[0];
            U = currentAngle[2];
            // Hnew = H cos(angle) + U sin(angle)
            currentAngle[0] = H*cos(angle) + U*sin(angle);
            // Lnew = L
            // Unew = -H sin(angle) + U cos(angle)
            currentAngle[2] = -H*sin(angle) + U*cos(angle);
        }
        // Si == &
        else if (replacementRule[i] == '&'){
            H = currentAngle[0];
            U = currentAngle[2];
            // Hnew = H cos(-angle) + U sin(-angle)
            currentAngle[0] = H*cos(-angle) + U*sin(-angle);
            // Lnew = L
            // Unew = -H sin(angle) + U cos(angle)
            currentAngle[2] = -H*sin(-angle) + U*cos(-angle);
        }
        // Si == '\'
        else if (replacementRule[i] == '\\') {
            L = currentAngle[1];
            U = currentAngle[2];
            // Hnew = H
            // Lnew = L cos(angle) - U sin(angle)
            currentAngle[1] = L*cos(angle) - U*sin(angle);
            // Unew = L sin(angle) + U cos(angle)
            currentAngle[2] = L*sin(angle) + U*cos(angle);
        }
        // Si == /
        else if (replacementRule[i] == '/'){
            L = currentAngle[1];
            U = currentAngle[2];
            // Hnew = H
            // Lnew = L cos(-angle) - U sin(-angle)
            currentAngle[1] = L*cos(-angle) - U*sin(-angle);
            // Unew = L sin(angle) + U cos(angle)
            currentAngle[2] = L*sin(-angle) + U*cos(-angle);
        }
        // Si == |
        else if (replacementRule[i] == '|'){
            H = currentAngle[0];
            L = currentAngle[1];
            // Hnew = -H
            currentAngle[0] = -H;
            // Lnew = -L
            currentAngle[1] = -L;
            // Unew = U
        }
        else{
            // Si e A1
            if (l_system.draw(replacementRule[i])) {
                Vector3D point;
                figure.points.push_back(currentPosition);
                point = currentPosition+currentAngle[0];
                figure.points.push_back(point);
                Face face;
                face.point_indexes.push_back(figure.points.size()-1);
                face.point_indexes.push_back(figure.points.size()-2);
                figure.faces.push_back(face);
                currentPosition = point;
            }
            // Si e A0
            if (l_system.draw(replacementRule[i]) == false) {
                currentPosition += currentAngle[0];
            }
        }
    }
    return figure;
}