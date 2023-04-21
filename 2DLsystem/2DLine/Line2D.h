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

    Line2D(const Point2D &p1, const Point2D &p2, const Color &color, double z1, double z2) : p1(p1), p2(p2),
                                                                                             color(color), z1(z1),
                                                                                             z2(z2) {}

    double z1;
    double z2;
};




#endif //ENGINE_LINE2D_H
