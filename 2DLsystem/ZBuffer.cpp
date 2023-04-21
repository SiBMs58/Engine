//
// Created by Siebe Mees on 21/04/2023.
//

#include "ZBuffer.h"


ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = std::numeric_limits<double>::infinity();
    //double negInf = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            this->push_back({posInf, posInf});
        }
    }
}

void ZBuffer::draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, const double z0, unsigned int x1, unsigned int y1, const double z1, Color &color) {
    img::Color c = img::Color(color.red, color.green, color.blue);

    double a = 

    if (x0 >= image.get_width() || y0 >= image.get_height() || x1 >= image.get_width() || y1 > image.get_height()) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << image.get_width() << " and height " << image.get_height();
        throw std::runtime_error(ss.str());
    }
    if (x0 == x1)
    {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            image(x0, i) = c;
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            image(i, y0) = c;
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                image(x0 + i, (unsigned int) round(y0 + m * i)) = c;
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                image((unsigned int) round(x0 + (i / m)), y0 + i) = c;
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                image((unsigned int) round(x0 - (i / m)), y0 - i) = c;
            }
        }
    }
}