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
        double rotateX = configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die();
        Matrix rotateXMatrix = figure.rotateX(rotateX);
        double rotateY = configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die();
        Matrix rotateYMatrix = figure.rotateY(rotateY);
        double rotateZ = configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die();
        Matrix rotatezMatrix = figure.rotateZ(rotateZ);

        figures.push_back(figure);
    }
    return figures;
}

