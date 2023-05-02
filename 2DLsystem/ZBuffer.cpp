//
// Created by Siebe Mees on 21/04/2023.
//

#include "ZBuffer.h"
#include <iostream>


ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = std::numeric_limits<double>::infinity();
    for (int i = 0; i < width; ++i) {
        // Columns
        vector<double> col;
        for (int j = 0; j < height; ++j) {
            // Row
            col.push_back(posInf);
        }
        this->push_back(col);
    }
}

void ZBuffer::draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, Color &color) {

    img::Color c = img::Color(color.red, color.green, color.blue);

    /**
     * a = length of the line
     * s = i, according to cursus
     * i = iterator of drawing the line (based on draw_line function)
     */

    if (x0 >= image.get_width() || y0 >= image.get_height() || x1 >= image.get_width() || y1 > image.get_height()) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << image.get_width() << " and height " << image.get_height();
        throw std::runtime_error(ss.str());
    }
    if (x0 == x1)
    {
        //special case for x0 == x1
        unsigned int a = std::max(y0, y1) - std::min(y0, y1)+1;
        unsigned int s = a;
        unsigned int i = std::min(y0, y1);
        while (i <= std::max(y0, y1)) {
            double p = s / a;
            if (y0 < y1){
                std::swap(z0, z1);
                std::swap(x0, x1);
            }
            double Zi = p / z0 + (1 - p) / z1;
            if (Zi < zbuffer[x0][i]) {
                image(x0, i) = c;
                zbuffer[x0][i] = Zi;
            }
            i++;
            s--;
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        unsigned int a = std::max(x0, x1) - std::min(x0, x1)+1;
        unsigned int s = a;
        unsigned int i = std::min(x0, x1);
        while (i <= std::max(x0, x1)) {
            double p = s/a;
            if (x0 < x1) {
                std::swap(z0, z1);
                std::swap(y0, y1);
            }
            double Zi = p/z0 + (1-p)/z1;
            if (Zi < zbuffer[i][y0]) {
                image(i, y0) = c;
                zbuffer[i][y0] = Zi;
            }
            i++;
            s--;
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            unsigned int a = (x1 - x0)+1;
            unsigned int s = a;
            unsigned int i = 0;
            while (i <= (x1 - x0)) {
                double p = s/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[x0 + i][(unsigned int) round(y0 + m * i)]) {
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = c;
                    zbuffer[x0 + i][(unsigned int) round(y0 + m * i)] = Zi;
                }
                s--;
                i++;
            }
        }
        else if (m > 1.0)
        {
            unsigned int a = (y1 - y0)+1;
            unsigned int s = a;
            unsigned int i = 0;
            while (i <= (y1 - y0)) {
                double p = i/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i]) {
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = c;
                    zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i] = Zi;
                }
                i++;
                s--;
            }
        }
        else if (m < -1.0)
        {
            unsigned int a = (y0 - y1)+1;
            unsigned int s = a;
            unsigned int i = 0;
            while (i <= (y0 - y1)) {
                double p = i/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[(unsigned int) round(x0 - (i / m))][y0 - i]) {
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = c;
                    zbuffer[(unsigned int) round(x0 - (i / m))][y0 - i] = Zi;
                }
                i++;
                s--;
            }
        }
    }

}