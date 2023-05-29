//
// Created by Siebe Mees on 02/04/2023.
//

#include "Wireframe3D.h"

Figure Wireframe3D::createCube() {
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

Figure Wireframe3D::createTetrahedron() {
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

Figure Wireframe3D::createIcosahedron() {
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

Figure Wireframe3D::createOctahedron() {
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

Figure Wireframe3D::createDodecahedron() {
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

Figure Wireframe3D::createCone(const int n, const double h) {
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

Figure Wireframe3D::createCylinder(const int n, const double h) {
    // Maak een figuur aan
    Figure figure;
    // Punten
    Face topFace;
    for (int i = 0; i < n+1; ++i) {
        Vector3D pointBottom = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), 0);
        topFace.point_indexes.push_back(i);
        figure.points.push_back(pointBottom);
    }
    figure.faces.push_back(topFace);
    Face bottomFace;
    for (int i = 0; i < n+1; ++i) {
        Vector3D pointTop = Vector3D::point(cos((2*i*M_PI)/n), sin((2*i*M_PI)/n), h);
        bottomFace.point_indexes.push_back(i+n+1);
        figure.points.push_back(pointTop);
    }
    figure.faces.push_back(bottomFace);
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

Figure Wireframe3D::createSphere(const double radius, const int n) {
    // Maak een figuur aan
    Figure figure;
    // 1. Genereer een icosahedron
    figure = createIcosahedron();
    // 2. Deel elke driehoek op in kleinere driehoeken
    int counter = 0;
    // 3. Herhaal stap 2 in totaal ‘n’ keer
    while (counter < n) {
        Figure bol;
        for (unsigned int j = 0; j < figure.faces.size(); ++j) {
            vector<int> pointIndexes = figure.faces[j].point_indexes;
            // Punten
            Vector3D pointA;
            pointA = figure.points[pointIndexes[0]];
            bol.points.push_back(pointA);
            Vector3D pointB;
            pointB = figure.points[pointIndexes[1]];
            bol.points.push_back(pointB);
            Vector3D pointC;
            pointC = figure.points[pointIndexes[2]];
            bol.points.push_back(pointC);
            Vector3D pointD;
            pointD = (pointA+pointB)/2;
            bol.points.push_back(pointD);
            Vector3D pointE;
            pointE = (pointA+pointC)/2;
            bol.points.push_back(pointE);
            Vector3D pointF;
            pointF = (pointB+pointC)/2;
            bol.points.push_back(pointF);
            // Faces
            Face face1; // ADE
            face1.point_indexes.push_back(bol.points.size()-6);
            face1.point_indexes.push_back(bol.points.size()-3);
            face1.point_indexes.push_back(bol.points.size()-2);
            bol.faces.push_back(face1);
            Face face2; // BFD
            face2.point_indexes.push_back(bol.points.size()-5);
            face2.point_indexes.push_back(bol.points.size()-1);
            face2.point_indexes.push_back(bol.points.size()-3);
            bol.faces.push_back(face2);
            Face face3; // CEF
            face3.point_indexes.push_back(bol.points.size()-4);
            face3.point_indexes.push_back(bol.points.size()-2);
            face3.point_indexes.push_back(bol.points.size()-1);
            bol.faces.push_back(face3);
            Face face4; // DFE
            face4.point_indexes.push_back(bol.points.size()-3);
            face4.point_indexes.push_back(bol.points.size()-1);
            face4.point_indexes.push_back(bol.points.size()-2);
            bol.faces.push_back(face4);
        }
        figure = bol;
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

Figure Wireframe3D::createTorus(const double r, const double R, const int n, const int m) {
    // Maak een figuur aan
    Figure figure;
    for (unsigned int i = 0; i <= n; ++i) {
        for (unsigned int j = 0; j <= m; ++j) {
            double u = (2*i*M_PI)/n;
            double v = (2*j*M_PI)/m;
            double Xuv = (R+r*cos(v))*cos(u);
            double Yuv = (R+r*cos(v))*sin(u);
            double Zuv = r*sin(v);
            // Punten
            Vector3D point = Vector3D::point(Xuv,Yuv,Zuv);
            figure.points.push_back(point);
            if (i <= n-1){
                // Faces
                Face face;
                face.point_indexes.push_back(figure.points.size()-1);
                face.point_indexes.push_back(figure.points.size()+m);
                face.point_indexes.push_back(figure.points.size()+m+1);
                face.point_indexes.push_back(figure.points.size());
                figure.faces.push_back(face);
            }
        }
    }
    // Fix bug last line
    figure.faces.pop_back();
    return figure;
}

Figure Wireframe3D::createBuckyBall() {
    // 1. Maak een icosahedron
    Figure buckyBall = createIcosahedron();

    // 2. Je kan de buckyball construeren door elk van de 20 driehoeken van de icosahedron op te delen in een gelijkzijdige zeshoek en drie driehoeken.
    vector<Vector3D> newPoints;
    vector<Face> newFaces;
    for (Face& face : buckyBall.faces) {
        Face currentFace = face;
        int newPointsSize = newPoints.size();
        for (int i = 0; i < 3; ++i) {
            newPoints.push_back(buckyBall.points[currentFace.point_indexes[i % 3]] + (1.0 / 3) * (buckyBall.points[currentFace.point_indexes[(i + 1) % 3]] - (buckyBall.points[currentFace.point_indexes[i % 3]])));
            newPoints.push_back(buckyBall.points[currentFace.point_indexes[i % 3]] + (2.0 / 3) * (buckyBall.points[currentFace.point_indexes[(i + 1) % 3]] - (buckyBall.points[currentFace.point_indexes[i % 3]])));
        }
        Face newFace; // Creëert point indexes e.g. 1, 2, 3, 4, 5, 6
        for (int i = newPointsSize; i < newPoints.size(); ++i) {
            newFace.point_indexes.push_back(i+1);
        }
        newFaces.push_back(newFace);
    }
    buckyBall.points = newPoints;
    buckyBall.faces = newFaces;

    // 3. Er on staan hierdoor piramides met vijf zijdes
    Face vijfhoek1;
    vijfhoek1.point_indexes = {1, 6, 12, 18, 24};
    Face vijfhoek2;
    vijfhoek2.point_indexes = {28, 82, 31, 3, 2};
    Face vijfhoek3;
    vijfhoek3.point_indexes = {9, 5, 4, 34, 42};
    Face vijfhoek4;
    vijfhoek4.point_indexes = {10, 46, 54, 15, 11};
    Face vijfhoek5;
    vijfhoek5.point_indexes = {17, 16, 58, 66, 21};
    Face vijfhoek6;
    vijfhoek6.point_indexes = {78, 27, 23, 22, 70};
    Face vijfhoek7;
    vijfhoek7.point_indexes = {33, 32, 117, 95, 39};
    Face vijfhoek8;
    vijfhoek8.point_indexes = {41, 40, 92, 51, 45};
    Face vijfhoek9;
    vijfhoek9.point_indexes = {53, 52, 107, 63, 57};
    Face vijfhoek10;
    vijfhoek10.point_indexes = {69, 65, 64, 104, 75};
    Face vijfhoek11;
    vijfhoek11.point_indexes = {110, 87, 81, 77, 76};
    Face vijfhoek12;
    vijfhoek12.point_indexes = {96, 109, 103, 108, 91};
    buckyBall.faces.push_back(vijfhoek1);
    buckyBall.faces.push_back(vijfhoek2);
    buckyBall.faces.push_back(vijfhoek3);
    buckyBall.faces.push_back(vijfhoek4);
    buckyBall.faces.push_back(vijfhoek5);
    buckyBall.faces.push_back(vijfhoek6);
    buckyBall.faces.push_back(vijfhoek7);
    buckyBall.faces.push_back(vijfhoek8);
    buckyBall.faces.push_back(vijfhoek9);
    buckyBall.faces.push_back(vijfhoek10);
    buckyBall.faces.push_back(vijfhoek11);
    buckyBall.faces.push_back(vijfhoek12);

    // 4. Fix bug indexering starts with 0
    for (Face& faceBuckyBall : buckyBall.faces){
        for (int& pointIndexes : faceBuckyBall.point_indexes) {
            pointIndexes--;
        }
    }

    return buckyBall;
}