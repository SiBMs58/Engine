//
// Created by Siebe Mees on 10/03/2023.
//

#ifndef ENGINE_LSYSTEM2D_H
#define ENGINE_LSYSTEM2D_H

#include <string>
#include <stack>
#include <cmath>
#include <map>
#include <vector>
#include <limits>
#include "../l_parser.h"
#include "2DLine/Line2D.h"
#include "../easy_image.h"
#include "./ZBuffer.h"

using namespace std;

using Lines2D = vector<Line2D>;

class Lsystem2D {
private:
    stack<Point2D> positionStack;
    stack<double> angleStack;
public:
    Lines2D drawLSystem(const LParser::LSystem2D &l_system, const vector<double> &lineColor);
    img::EasyImage draw2DLines(Lines2D &lines, const int size, const vector<double> &backgroundColor);
    img::EasyImage drawZbufLines(Lines2D &lines, const int size, const vector<double> &backgroundColor);
    map<string, double> calculateImageSizeAndScales(Lines2D &lines, const int size);
};


#endif //ENGINE_LSYSTEM2D_H
