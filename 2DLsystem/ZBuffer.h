//
// Created by Siebe Mees on 21/04/2023.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include <vector>
#include <sstream>
#include <utility>
#include "../easy_image.h"
#include "./2DLine/Line2D.h"
#include "../vector3d.h"
#include "../3DLsystem/3DFigure/Light.h"

class ZBuffer: public std::vector<std::vector<double>> {
public:
    //Constructor: maakt een Z-Buffer van de correcte
    //grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height);

    void draw_zbuf_line(ZBuffer &, img::EasyImage &, const unsigned int x0, const unsigned int y0, const double z0,const unsigned int x1, const unsigned int y1, const double z1, Color &color);
    void draw_zbuf_triag(ZBuffer &, img::EasyImage &, const Vector3D &A, const Vector3D &B, const Vector3D &C,
                         double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoefficient, Lights3D& lights, Vector3D& eye);
};


#endif //ENGINE_ZBUFFER_H
