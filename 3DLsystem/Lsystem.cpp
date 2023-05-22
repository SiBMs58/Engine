//
// Created by Siebe Mees on 17/03/2023.
//

#include <fstream>
#include "Lsystem.h"
#include "Wireframe3D.h"

void Lsystem::applyDrawFunctions(Figure &figure, const vector<double> vectorColorFigure, const double scale, const double rotateX, const double rotateY, const double rotateZ, const vector<double> vectorCenter, const vector<double> figureColor) {
    // Color
    Color colorFigure;

    colorFigure.red = vectorColorFigure[0] * 255;
    colorFigure.green = vectorColorFigure[1] * 255;
    colorFigure.blue = vectorColorFigure[2] * 255;

    figure.color = colorFigure;

    // Schalen
    Matrix scaleMatrix = figure.scaleFigure(scale);
    // RotateX
    Matrix rotateXMatrix = figure.rotateX(rotateX);
    // RotateT
    Matrix rotateYMatrix = figure.rotateY(rotateY);
    // RotateZ
    Matrix rotateZMatrix = figure.rotateZ(rotateZ);
    // Translatie over vector
    Vector3D translatieVector = Vector3D::point(vectorCenter[0], vectorCenter[1], vectorCenter[2]);
    Matrix translationMatrix = figure.translate(translatieVector);

    Matrix transformationMatrix =
            scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * translationMatrix;
    applyTransformation(figure, transformationMatrix);

    // Color
    figure.color.red = figureColor[0] * 255;
    figure.color.green = figureColor[1] * 255;
    figure.color.blue = figureColor[2] * 255;
}

vector<Figure> Lsystem::generateFigures(const ini::Configuration &configuration) {
    // Maak een wirefame aan
    Wireframe3D wireframe;
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    // Maak figures aan
    vector<Figure> figures;
    for (int i = 0; i < nrFigures; ++i) {
        // Maak een figuur aan
        Figure figure;
        string figureType = configuration["Figure"+to_string(i)]["type"].as_string_or_die();
        bool fractal = false;
        if (figureType == "Cube"){
            figure = wireframe.createCube();
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
            figure = wireframe.createTetrahedron();
        } else if (figureType == "Icosahedron"){
            figure = wireframe.createIcosahedron();
        } else if (figureType == "Octahedron"){
            figure = wireframe.createOctahedron();
        } else if (figureType == "Dodecahedron"){
            figure = wireframe.createDodecahedron();
        } else if (figureType == "Cone") {
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure"+to_string(i)]["height"].as_double_or_die();
            figure = wireframe.createCone(n, height);
        } else if (figureType == "Cylinder"){
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure"+to_string(i)]["height"].as_double_or_die();
            figure = wireframe.createCylinder(n, height);
        } else if (figureType == "Sphere"){
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            figure = wireframe.createSphere(1 ,n);
        } else if (figureType == "Torus"){
            double r = configuration["Figure"+to_string(i)]["r"].as_double_or_die();
            double R = configuration["Figure"+to_string(i)]["R"].as_double_or_die();
            double m = configuration["Figure"+to_string(i)]["m"].as_double_or_die();
            double n = configuration["Figure"+to_string(i)]["n"].as_double_or_die();
            figure = wireframe.createTorus(r, R, n ,m);
        } else if (figureType == "BuckyBall") {
            figure = wireframe.createBuckeyBall();
        } else if (figureType == "LineDrawing") {
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
        } else if (figureType == "FractalTetrahedron") {
            Figure tetrahedron = wireframe.createTetrahedron();
            Figures3D fractalTetrahedron;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(tetrahedron, fractalTetrahedron, nrIterations, fractalScale);
            figures = fractalTetrahedron;
            fractal = true;
        } else if (figureType == "FractalCube") {
            Figure cube = wireframe.createCube();
            Figures3D fractalCube;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(cube, fractalCube, nrIterations, fractalScale);
            figures = fractalCube;
            fractal = true;
        } else if (figureType == "FractalIcosahedron") {
            Figure icosahedron = wireframe.createIcosahedron();
            Figures3D fractalIcosahedron;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(icosahedron, fractalIcosahedron, nrIterations, fractalScale);
            figures = fractalIcosahedron;
            fractal = true;
        } else if (figureType == "FractalOctahedron") {
            Figure octahedron = wireframe.createOctahedron();
            Figures3D fractalOctahedron;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(octahedron, fractalOctahedron, nrIterations, fractalScale);
            figures = fractalOctahedron;
            fractal = true;
        } else if (figureType == "FractalDodecahedron") {
            Figure dodecahedron = wireframe.createDodecahedron();
            Figures3D fractalDodecahedron;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(dodecahedron, fractalDodecahedron, nrIterations, fractalScale);
            figures = fractalDodecahedron;
            fractal = true;
        } else if (figureType == "FractalBuckyBall") {
            Figure buckyBall = wireframe.createBuckeyBall();
            Figures3D fractalBuckyBall;
            int nrIterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die();
            int fractalScale = configuration["Figure"+to_string(i)]["fractalScale"].as_int_or_die();
            generateFractal(buckyBall, fractalBuckyBall, nrIterations, fractalScale);
            figures = fractalBuckyBall;
            fractal = true;
        }

        vector<double> vectorColorFigure = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_die();
        double scale = configuration["Figure"+to_string(i)]["scale"].as_double_or_die();
        double rotateX = configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die();
        double rotateY = configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die();
        double rotateZ = configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die();
        vector<double> vectorCenter = configuration["Figure"+to_string(i)]["center"].as_double_tuple_or_die();
        vector<double> figureColor = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_die();

        if (fractal) {
            for (Figure &figure : figures) {
                applyDrawFunctions(figure, vectorColorFigure, scale, rotateX, rotateY, rotateZ, vectorCenter, figureColor);
            }
        } else {
            applyDrawFunctions(figure, vectorColorFigure, scale, rotateX, rotateY, rotateZ, vectorCenter, figureColor);

            figures.push_back(figure);
        }
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
                Line2D line = Line2D(point1, point2, figures[i].color, vectorPoint1.z, vectorPoint2.z);
                lines.push_back(line);
            }
            // Add line from last point_index to first one
            int point_index1 = figures[i].faces[j].point_indexes[figures[i].faces[j].point_indexes.size()-1];
            int point_index2 = figures[i].faces[j].point_indexes[0];
            Vector3D vectorPoint1 = Vector3D::point(figures[i].points[point_index1].x, figures[i].points[point_index1].y, figures[i].points[point_index1].z);
            Vector3D vectorPoint2 = Vector3D::point(figures[i].points[point_index2].x, figures[i].points[point_index2].y, figures[i].points[point_index2].z);
            Point2D point1 = generatePoint2D(vectorPoint1, 1);
            Point2D point2 = generatePoint2D(vectorPoint2, 1);

            Line2D line = Line2D(point1, point2, figures[i].color, vectorPoint1.z, vectorPoint1.z);
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

vector<Face> Lsystem::triangulate(const Face &face) {
    vector<Face> triangles;
    for (int i = 1; i < face.point_indexes.size()-1; ++i) {
        Face newFace;
        newFace.point_indexes.push_back(face.point_indexes[0]);
        newFace.point_indexes.push_back(face.point_indexes[i]);
        newFace.point_indexes.push_back(face.point_indexes[i+1]);
        triangles.push_back(newFace);
    }

    return triangles;
}

img::EasyImage Lsystem::drawZbufTriangles(vector<Figure> &figures, const int size, const vector<double> &backgroundColor) {
    Lsystem lsystem3D;
    // 1. Triangulatie
    for (Figure &fig : figures) {
        vector<Face> triangleFaces;
        for (Face &face : fig.faces) {
            vector<Face> triangles = lsystem3D.triangulate(face);
            triangleFaces.insert(triangleFaces.end(), triangles.begin(), triangles.end());
        }
        fig.faces = triangleFaces;
    }

    Lines2D lines = doProjection(figures);
    Lsystem2D lsystem2D;
    map<string, double> resultCalculations = lsystem2D.calculateImageSizeAndScales(lines, size);
    // 2. Bereken de image grootte
    double imageX = resultCalculations["imageX"];
    double imageY = resultCalculations["imageY"];

    ZBuffer zbuffer = ZBuffer(ceil(imageX), ceil(imageY));
    img::EasyImage image(imageX, imageY);
    // backgroundColor van de image
    for(unsigned int i = 0; i < image.get_width(); i++) {
        for(unsigned int j = 0; j < image.get_height(); j++) {
            image(i,j).red = backgroundColor[0]*255;
            image(i,j).green = backgroundColor[1]*255;
            image(i,j).blue = backgroundColor[2]*255;
        }
    }

    // 3. Bereken de afstand van het projectievlak tot de oorsprong
    double d = resultCalculations["d"];
    // 4. Bereken de offsets dx en dy
    double dx = resultCalculations["dx"];
    double dy = resultCalculations["dy"];
    // Het Z-Buffer algoritme zélf wordt driehoek per driehoek toegepast
    for (Figure &fig : figures) {
        for (int i = 0; i < fig.faces.size(); ++i) {
            int point_index1 = fig.faces[i].point_indexes[0];
            int point_index2 = fig.faces[i].point_indexes[1];
            int point_index3 = fig.faces[i].point_indexes[2];
            Vector3D vectorPoint1 = Vector3D::point(fig.points[point_index1].x, fig.points[point_index1].y, fig.points[point_index1].z);
            Vector3D vectorPoint2 = Vector3D::point(fig.points[point_index2].x, fig.points[point_index2].y, fig.points[point_index2].z);
            Vector3D vectorPoint3 = Vector3D::point(fig.points[point_index3].x, fig.points[point_index3].y, fig.points[point_index3].z);

            zbuffer.draw_zbuf_triag(zbuffer, image, vectorPoint1, vectorPoint2, vectorPoint3, d, dx, dy, fig.color);

        }
    }

    return image;
}

void Lsystem::generateFractal(Figure &fig, Figures3D &fractal, const int nr_iterations, double scale) {

    /*
     *     fractal.emplace_back(figr);

    while (nr_iterations > 0) {
        Figures3D newFractals;
        for (Figure fig: fractal) {
            for (int point = 0; point < fig.points.size(); ++point) {

                Figure fig2 = fig;
                fig2.applyTransformation(TransMatrix::scaleFigure(1 / scale));

                Vector3D pi = fig.points[point];
                Vector3D pib = fig2.points[point];

                fig2.applyTransformation(TransMatrix::translate(pi - pib));

                newFractals.emplace_back(fig2);
            }
        }
        fractal = newFractals;
        nr_iterations--;
    }
     */

    fractal.push_back(fig);

    for (int i = 0; i < nr_iterations; ++i) {
        Figures3D newFractals;
        for (Figure fig: fractal) {
            for (int j = 0; j < fig.points.size(); ++j) {

                Figure fig2 = fig;
                Matrix scaleMatrix = fig2.scaleFigure(1/scale);
                applyTransformation(fig2, scaleMatrix);

                Vector3D pi = fig.points[j];
                Vector3D pib = fig2.points[j];

                Matrix transMatrix = fig2.translate(pi - pib);
                applyTransformation(fig2, transMatrix);

                newFractals.emplace_back(fig2);
            }
        }
        fractal = newFractals;
        /*// 1. Kopieer de oorspronkelijke figuur
        Figure copyFig = fig;
        // 2. Schaal de punten van de nieuwe figuur
        copyFig.scaleFigure(scale);

        fractal.push_back(copyFig);*/
    }

}