//
// Created by Siebe Mees on 10/03/2023.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H

#include "Point2D.h"
#include "Color.h"

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;

    Line2D(const Point2D &p1, const Point2D &p2, const Color &color) : p1(p1), p2(p2), color(color) {}
};


#endif //ENGINE_LINE2D_H
