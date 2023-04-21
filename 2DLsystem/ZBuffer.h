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

class ZBuffer: public std::vector<std::vector<double>> {
public:
//Constructor: maakt een Z-Buffer van de correcte
//grootte aan en initialiseert alle velden op +inf
ZBuffer(const int width, const int height);

void draw_zbuf_line(ZBuffer &, img::EasyImage &, const unsigned int x0, const unsigned int y0, const double z0,const unsigned int x1, const unsigned int y1, const double z1, Color &color);

};


#endif //ENGINE_ZBUFFER_H
