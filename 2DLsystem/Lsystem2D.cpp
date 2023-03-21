//
// Created by Siebe Mees on 10/03/2023.
//

#include "Lsystem2D.h"

string getReplacementRule(const LParser::LSystem2D &Lsystem) {
    string initiator = Lsystem.get_initiator();
    unsigned int nrIterations = Lsystem.get_nr_iterations();
    string replacementRule="";
    if (nrIterations > 0){
        // Loop nr of iterations
        for (unsigned int i = 0; i < nrIterations; ++i) {
            replacementRule ="";
            for (unsigned int i = 0; i < initiator.size(); ++i) {
                if (initiator[i] == '-' || initiator[i] == '+' || initiator[i] == ')' || initiator[i] == '('){
                    replacementRule += initiator[i];
                } else {
                    replacementRule += Lsystem.get_replacement(initiator[i]);
                }
            }
            initiator = replacementRule;
        }
    } else {
        for (unsigned int i = 0; i < initiator.size(); ++i) {
            if (initiator[i] == '-' || initiator[i] == '+' || initiator[i] == ')' || initiator[i] == '('){
                replacementRule += initiator[i];
            } else {
                replacementRule += Lsystem.get_replacement(initiator[i]);
            }
        }
    }
    return replacementRule;
}

img::EasyImage Lsystem2D::draw2DLines(Lines2D &lines, const int size, const vector<double> &backgroundColor) {
    // 1. Bepaal Xmin, Xmax, Ymin en Ymax
    double Xmax = lines[0].p1.x;
    double Xmin = lines[0].p1.x;
    double Ymax = lines[0].p1.y;
    double Ymin = lines[0].p1.y;
    for (int i = 0; i < lines.size(); ++i) {
        if (Xmax < lines[i].p1.x){
            Xmax = lines[i].p1.x;
        }
        if (Xmin > lines[i].p1.x){
            Xmin = lines[i].p1.x;
        }
        if (Xmax < lines[i].p2.x){
            Xmax = lines[i].p2.x;
        }
        if (Xmin > lines[i].p2.x){
            Xmin = lines[i].p2.x;
        }
        if (Ymax < lines[i].p1.y){
            Ymax = lines[i].p1.y;
        }
        if (Ymin > lines[i].p1.y){
            Ymin = lines[i].p1.y;
        }
        if (Ymax < lines[i].p2.y){
            Ymax = lines[i].p2.y;
        }
        if (Ymin > lines[i].p2.y){
            Ymin = lines[i].p2.y;
        }
    }
    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    // 2. Bereken de grootte van de image
    double imageX = size*((Xrange)/max(Xrange, Yrange));
    double imageY = size*((Yrange)/max(Xrange, Yrange));
    // 2.1 Maak image
    img::EasyImage image(imageX, imageY);
    // 2.2 Set backgroundColor van de image
    for(unsigned int i = 0; i < image.get_width(); i++) {
        for(unsigned int j = 0; j < image.get_height(); j++) {
            image(i,j).red = backgroundColor[0]*255;
            image(i,j).green = backgroundColor[1]*255;
            image(i,j).blue = backgroundColor[2]*255;
        }
    }
    // 3. Schaal de lijntekening
    // 3.1 Bereken de schaalfactor
    double d = 0.95*(imageX/Xrange);
    // 3.2 Vermenigvuldig de coördinaten van alle punten met d
    for (int i = 0; i < lines.size(); ++i) {
        lines[i].p1.x = lines[i].p1.x*d;
        lines[i].p1.y = lines[i].p1.y*d;
        lines[i].p2.x = lines[i].p2.x*d;
        lines[i].p2.y = lines[i].p2.y*d;
    }
    // 4. Verschuif de lijntekening
    // 4.1 Bereken
    double DCx = d*((Xmin+Xmax)/2);
    double DCy = d*((Ymin+Ymax)/2);
    double dx = (imageX/2)-DCx;
    double dy = (imageY/2)-DCy;
    // 4.2 Tel bij alle punten van de lijntekening (dx, dy) op.
    for (int i = 0; i < lines.size(); ++i) {
        lines[i].p1.x = lines[i].p1.x+dx;
        lines[i].p1.y = lines[i].p1.y+dy;
        lines[i].p2.x = lines[i].p2.x+dx;
        lines[i].p2.y = lines[i].p2.y+dy;
    }
    // 5. Rond de coördinaten van de punten af, en teken de lijnen op de image
    for (int i = 0; i < lines.size(); ++i) {
        img::Color lineColor = img::Color(lines[i].color.red, lines[i].color.green, lines[i].color.blue);
        image.draw_line(lround(lines[i].p1.x), lround(lines[i].p1.y), lround(lines[i].p2.x), lround(lines[i].p2.y), lineColor);
    }
    return image;
}

Lines2D Lsystem2D::drawLSystem(const LParser::LSystem2D &l_system, const vector<double> &lineColor) {
    // Maak lines2D aan
    Lines2D lines;
    double startingAngle = l_system.get_starting_angle();
    // Convert to radians
    startingAngle = startingAngle*M_PI/180;
    set<char> alphabet = l_system.get_alphabet();
    double angle = l_system.get_angle();
    // Convert to radians
    angle = angle*M_PI/180;
    // Get total sting based on initiator
    string replacementRule = getReplacementRule(l_system);
    // 1. We starten in een willekeurige positie, typisch positie (x, y) = (0, 0), en in de richting ↵ = ↵0 radialen.
    Point2D startingPoint;
    startingPoint.x = 0;
    startingPoint.y = 0;
    double currentAngle = startingAngle;
    // Make for a lineColor
    Color color;
    color.red = lineColor[0]*255;
    color.green = lineColor[1]*255;
    color.blue = lineColor[2]*255;
    // 2. We overlopen de string S = s1s2 . . . sk van links naar rechts en voeren afhankelijk van het huidige symbool si volgende actie uit:
    for (int i = 0; i < replacementRule.size(); ++i) {
        // Si = +
        if (replacementRule[i] == '+') {
            currentAngle += angle;
        }
        // Si = -
        else if (replacementRule[i] == '-') {
            currentAngle -= angle;
        }
        else if (replacementRule[i] == '('){
            positionStack.push(startingPoint);
            angleStack.push(currentAngle);
        }
        else if (replacementRule[i] == ')'){
            startingPoint = positionStack.top();
            positionStack.pop();
            currentAngle = angleStack.top();
            angleStack.pop();
        }
        else{
            // Si e A1
            if (l_system.draw(replacementRule[i])) {
                Point2D point;
                point.x = startingPoint.x+cos(currentAngle);
                point.y = startingPoint.y+sin(currentAngle);
                Line2D line = Line2D(startingPoint, point, color);
                lines.push_back(line);
                startingPoint = point;
            }
            // Si e A0
            if (l_system.draw(replacementRule[i]) == false) {
                startingPoint.x = startingPoint.x+cos(currentAngle);
                startingPoint.y = startingPoint.y+sin(currentAngle);
            }
        }

    }
    return lines;
}